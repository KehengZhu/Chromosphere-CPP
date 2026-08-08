#include "eos.hpp"
#include "profiling.hpp"
#include "parallel.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace chromosphere {
namespace {

constexpr double LN_10     = 2.30258509299404568402;
constexpr double PI        = 3.14159265358979323846;

const std::map<std::string, std::string>& expected_metadata() {
    static const std::map<std::string, std::string> values = {
        {"format_version", "1"},
        {"table_id", "chromosphere2026_gamma1_hydrogen_v1"},
        {"material", "pure_H"},
        {"excitation", "0"},
        {"fermi_gas", "0"},
        {"coulomb_correction", "0"},
        {"ground_stat_weight", "1"},
        {"saha_prefactor", "coefficient_one"},
        {"axis_1", "log10_T_K"},
        {"axis_2", "log10_nH_m-3"},
        {"value", "Gamma1"},
        {"k_b_J_K", "1.380649e-23"},
        {"m_e_kg", "9.1093837015e-31"},
        {"m_H_kg", "1.6726219e-27"},
        {"h_J_s", "6.62607015e-34"},
        {"chi_H_J", "2.179872361e-18"},
        {"generator", "util/eos/tabulate_gamma.f90"},
        {"generation_date", "2026-07-19"},
        {"source_revision", "crash_gamma_eos_stage0_2_v1"},
    };
    return values;
}

std::size_t parse_dimension(const std::map<std::string, std::string>& metadata,
                            const std::string& key, const std::string& path) {
    const auto found = metadata.find(key);
    if (found == metadata.end())
        throw std::runtime_error("Gamma1 table: missing metadata key '" + key
                                 + "' in '" + path + "'");
    std::size_t consumed = 0;
    unsigned long value = 0;
    try {
        value = std::stoul(found->second, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("Gamma1 table: invalid dimension " + key + "='"
                                 + found->second + "' in '" + path + "'");
    }
    if (consumed != found->second.size() || value < 2)
        throw std::runtime_error("Gamma1 table: invalid dimension " + key + "='"
                                 + found->second + "' in '" + path + "'");
    return static_cast<std::size_t>(value);
}

struct Row {
    double log_t_10;
    double log_n_10;
    double gamma1;
};

bool eos_operation_counting_enabled = false;
struct alignas(64) EosCountSlot { EosOperationCounts counts; };
std::array<EosCountSlot, kMaximumParallelThreads> eos_operation_slots;

inline EosOperationCounts& eos_operation_counter() noexcept {
    return eos_operation_slots[static_cast<std::size_t>(parallel_thread_index())].counts;
}

inline double require_finite_positive(double value, const char* name) {
    if (!(value > 0.0) || !std::isfinite(value))
        throw std::domain_error(std::string(name) + " must be positive and finite");
    return value;
}

inline double temperature_log(double temperature) {
    if (eos_operation_counting_enabled) ++eos_operation_counter().temperature_logs;
    return std::log(temperature);
}

inline double n_h_log(double n_h) {
    if (eos_operation_counting_enabled) ++eos_operation_counter().n_h_logs;
    return std::log(n_h);
}

struct EosCoordinates {
    double temperature;
    double n_h;
    double log_temperature;
    double log_n_h;
};

EosCoordinates eos_coordinates(double n_h, double temperature) {
    require_finite_positive(n_h, "n_H");
    require_finite_positive(temperature, "temperature");
    return EosCoordinates{temperature, n_h, temperature_log(temperature),
                          n_h_log(n_h)};
}

EosCoordinates eos_coordinates_with_log_n(double n_h, double log_n_h,
                                           double temperature) {
    require_finite_positive(n_h, "n_H");
    require_finite_positive(temperature, "temperature");
    return EosCoordinates{temperature, n_h, temperature_log(temperature), log_n_h};
}

EosCoordinates eos_coordinates_with_logs(double n_h, double temperature,
                                          double log_temperature,
                                          double log_n_h) {
    require_finite_positive(n_h, "n_H");
    require_finite_positive(temperature, "temperature");
    if (!std::isfinite(log_temperature) || !std::isfinite(log_n_h))
        throw std::domain_error("EOS logarithmic coordinates must be finite");
    return EosCoordinates{temperature, n_h, log_temperature, log_n_h};
}

inline double saha_log_c() {
    static const double value = 1.5 * std::log(2.0 * PI * eos_constants::m_e
                                               * eos_constants::k_b /
                                               (eos_constants::h * eos_constants::h));
    return value;
}

double saha_from_coordinates(const EosCoordinates& c) {
    const double log_a = saha_log_c() + 1.5 * c.log_temperature
                       - eos_constants::chi_h /
                         (eos_constants::k_b * c.temperature)
                       - c.log_n_h;
    if (log_a >= 0.0)
        return 2.0 / (1.0 + std::sqrt(1.0 + 4.0 * std::exp(-log_a)));
    const double u = std::exp(0.5 * log_a);
    return 2.0 * u / (u + std::sqrt(u * u + 4.0));
}

double saha_dx_d_temperature(double x, double temperature) {
    const double theta = eos_constants::chi_h /
                         (eos_constants::k_b * temperature);
    return x * (1.0 - x) * (1.5 + theta) /
           ((2.0 - x) * temperature);
}

// ---------------------------------------------------------------------------
// One authoritative caloric evaluation. EosCoordinates carries the validated
// log-space coordinates, so Saha, density guards and Gamma1 reuse the same two
// logarithms instead of recomputing them independently.
// ---------------------------------------------------------------------------
struct CaloricEval {
    double n_h;
    double log_n_h;
    double log_temperature;
    double x;
    double pressure;
    double internal_energy;
    double heat_capacity;   // NaN unless produced by caloric_eval_full
    double n_e() const { return x * n_h; }
    double n_hi() const { return (1.0 - x) * n_h; }
};

CaloricEval caloric_eval_from_coordinates(const EosCoordinates& c) {
    CaloricEval ev;
    ev.n_h = c.n_h;
    ev.log_n_h = c.log_n_h;
    ev.log_temperature = c.log_temperature;
    ev.x = saha_from_coordinates(c);
    ev.pressure = (1.0 + ev.x) * ev.n_h * eos_constants::k_b * c.temperature;
    ev.internal_energy = 1.5 * ev.pressure + ev.x * ev.n_h * eos_constants::chi_h;
    ev.heat_capacity = std::numeric_limits<double>::quiet_NaN();
    return ev;
}

CaloricEval caloric_eval(double rho_total, double temperature) {
    const double n_h = rho_total / eos_constants::m_h;
    return caloric_eval_from_coordinates(eos_coordinates(n_h, temperature));
}

CaloricEval caloric_eval_with_log_n(double rho_total, double temperature,
                                    double n_h, double log_n_h) {
    (void)rho_total;
    return caloric_eval_from_coordinates(
        eos_coordinates_with_log_n(n_h, log_n_h, temperature));
}

CaloricEval caloric_eval_with_logs(double rho_total, double temperature,
                                   double log_temperature) {
    const double n_h = rho_total / eos_constants::m_h;
    return caloric_eval_from_coordinates(
        eos_coordinates_with_logs(n_h, temperature, log_temperature, n_h_log(n_h)));
}

// The heat capacity depends on the Saha fraction only through x, so it can be
// completed from an energy-only evaluation without a second Saha solve.
double caloric_capacity(double n_h, double x, double temperature) {
    const double dx_d_t = saha_dx_d_temperature(x, temperature);
    return n_h * (1.5 * eos_constants::k_b * (1.0 + x)
                  + (1.5 * eos_constants::k_b * temperature
                     + eos_constants::chi_h) * dx_d_t);
}

CaloricEval caloric_eval_full(double rho_total, double temperature) {
    CaloricEval ev = caloric_eval(rho_total, temperature);
    ev.heat_capacity = caloric_capacity(ev.n_h, ev.x, temperature);
    return ev;
}

CaloricEval caloric_eval_full_with_log_n(double rho_total, double temperature,
                                         double n_h, double log_n_h) {
    CaloricEval ev = caloric_eval_with_log_n(rho_total, temperature, n_h, log_n_h);
    ev.heat_capacity = caloric_capacity(ev.n_h, ev.x, temperature);
    return ev;
}

CaloricState caloric_state_of(const CaloricEval& ev, double temperature) {
    const double capacity = std::isnan(ev.heat_capacity)
        ? caloric_capacity(ev.n_h, ev.x, temperature) : ev.heat_capacity;
    return CaloricState{ev.n_h, ev.x, ev.n_e(), ev.n_hi(), ev.pressure,
                        ev.internal_energy, capacity};
}

double equilibrium_heat_capacity_impl(double rho_total, double temperature) {
    return caloric_eval_full(rho_total, temperature).heat_capacity;
}

std::string trim(const std::string& input) {
    const std::string whitespace = " \t\r\n";
    const std::size_t first = input.find_first_not_of(whitespace);
    if (first == std::string::npos) return std::string();
    const std::size_t last = input.find_last_not_of(whitespace);
    return input.substr(first, last - first + 1);
}

std::size_t lower_cell(const std::vector<double>& axis, double value) {
    if (value <= axis.front()) return 0;
    if (value >= axis.back()) return axis.size() - 2;
    return static_cast<std::size_t>(
        std::upper_bound(axis.begin(), axis.end(), value) - axis.begin() - 1);
}

void push_unique(std::vector<double>& values, double value) {
    if (values.empty() || value != values.back()) values.push_back(value);
}

} // namespace

void set_eos_operation_counting(bool enabled) noexcept {
    eos_operation_counting_enabled = enabled;
}

void reset_eos_operation_counts() noexcept {
    for (auto& slot : eos_operation_slots)
        slot.counts = EosOperationCounts{};
}

EosOperationCounts eos_operation_counts() noexcept {
    EosOperationCounts merged;
    for (const auto& slot : eos_operation_slots) {
        merged.gamma1_queries += slot.counts.gamma1_queries;
        merged.temperature_logs += slot.counts.temperature_logs;
        merged.n_h_logs += slot.counts.n_h_logs;
    }
    return merged;
}

double saha_ionization_fraction_n_h(double n_h, double temperature) {
    return saha_from_coordinates(eos_coordinates(n_h, temperature));
}

double equilibrium_density_from_pressure(double pressure, double temperature) {
    require_finite_positive(pressure, "pressure");
    require_finite_positive(temperature, "temperature");

    // log S(T), written with exactly the same constant and term order as the
    // log_a expression in saha_ionization_fraction_n_h above (log_a = log S - log n_H).
    const double log_s = saha_log_c() + 1.5 * temperature_log(temperature)
                       - eos_constants::chi_h / (eos_constants::k_b * temperature);
    const double y = pressure / (eos_constants::k_b * temperature);
    const double d = std::log(y) - log_s;

    // x = sqrt(S/(S+y)) = sqrt(1/(1+e^d)). Two algebraically identical branches,
    // each chosen so its exponential argument is non-positive: no overflow for
    // either y >> S (deep neutral) or y << S (fully ionized).
    double x;
    if (d <= 0.0) {
        x = 1.0 / std::sqrt(1.0 + std::exp(d));
    } else {
        const double e = std::exp(-d);                  // = S/y
        x = std::sqrt(e / (1.0 + e));
    }
    const double n_h = y / (1.0 + x);
    return n_h * eos_constants::m_h;
}

double saha_ionization_fraction(double rho_total, double temperature) {
    require_finite_positive(rho_total, "rho_total");
    return saha_ionization_fraction_n_h(rho_total / eos_constants::m_h, temperature);
}

double equilibrium_internal_energy(double rho_total, double temperature) {
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    return caloric_eval(rho_total, temperature).internal_energy;
}

double equilibrium_heat_capacity(double rho_total, double temperature) {
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    return equilibrium_heat_capacity_impl(rho_total, temperature);
}

CaloricState equilibrium_caloric_state(double rho_total, double temperature) {
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    return caloric_state_of(caloric_eval_full(rho_total, temperature),
                            temperature);
}

CaloricState equilibrium_caloric_state_from_logs(
    double rho_total, double temperature,
    double log_temperature, double log_n_h) {
    require_finite_positive(rho_total, "rho_total");
    const double n_h = rho_total/eos_constants::m_h;
    const EosCoordinates coordinates = eos_coordinates_with_logs(
        n_h, temperature, log_temperature, log_n_h);
    CaloricEval ev = caloric_eval_from_coordinates(coordinates);
    ev.heat_capacity = caloric_capacity(ev.n_h, ev.x, temperature);
    return caloric_state_of(ev, temperature);
}

namespace {

GammaState gamma_state_from_eval(const EosGammaTable& table,
                                 const CaloricEval& ev, double temperature,
                                 bool debug_clamp) {
    return GammaState{
        ev.x, 1.0 + ev.pressure / ev.internal_energy,
        table.gamma1_from_logs(temperature, ev.n_h, ev.log_temperature,
                               ev.log_n_h, debug_clamp)};
}

/// Outcome of one caloric inversion. When `has_state` is set, `eval` is a valid
/// caloric evaluation AT `temperature`, so the caller can build the GammaState
/// (and the partial pressures) without paying for another Saha solve.
struct InversionResult {
    double temperature;
    bool has_state;
    CaloricEval eval;
};

InversionResult invert_caloric(const EosGammaTable& table, double rho_total,
                               double internal_energy,
                               double temperature_guess, bool debug_clamp);

} // namespace

GammaState gamma_state(const EosGammaTable& table, double rho_total,
                       double temperature, bool debug_clamp) {
    if (table.empty()) throw std::logic_error("Gamma1 table is empty");
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    const double n_h = rho_total / eos_constants::m_h;
    const EosCoordinates coordinates = eos_coordinates(n_h, temperature);
    table.require_n_h_in_bounds_from_log(n_h, coordinates.log_n_h, debug_clamp);
    return gamma_state_from_eval(
        table, caloric_eval_from_coordinates(coordinates), temperature, debug_clamp);
}

double temperature_from_rho_eint(const EosGammaTable& table, double rho_total,
                                 double internal_energy,
                                 double temperature_guess, bool debug_clamp) {
    return invert_caloric(table, rho_total, internal_energy, temperature_guess,
                          debug_clamp).temperature;
}

namespace {

InversionResult invert_caloric(const EosGammaTable& table, double rho_total,
                               double internal_energy,
                               double temperature_guess, bool debug_clamp) {
    if (table.empty()) throw std::logic_error("Gamma1 table is empty");
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(internal_energy, "internal_energy");
    const double n_h = rho_total / eos_constants::m_h;
    require_finite_positive(n_h, "n_H");
    const double log_n_h = n_h_log(n_h);
    table.require_n_h_in_bounds_from_log(n_h, log_n_h, debug_clamp);

    profile_note_inversion_call();
    const double scale = std::max(std::numeric_limits<double>::min(), internal_energy);

    double lower = table.min_temperature();
    double upper = table.max_temperature();

    // ---- zero/one-update fast path -------------------------------------------
    // Measured on the production column, ~97% of all inversions converge with at
    // most ONE Newton update from the supplied guess, yet each still paid for the
    // two table-wide bracket-endpoint energies (a further ~470M Saha solves per
    // 20 s run) and evaluated energy and heat capacity as two separate Saha
    // solves. Here the guess is evaluated ONCE (energy + capacity fused), the
    // step is taken with exactly the safeguards the full solver applies, and the
    // candidate is evaluated once. Anything else falls through to the untouched
    // safeguarded Newton/bisection solver below.
    //
    // The fast path is gated on the guess AND the accepted candidate lying inside
    // a 1e-9-relative inset of the table temperature range. That margin is many
    // orders above the 64-eps window in which the bracket block below returns the
    // endpoint temperature itself or throws out_of_range, so no call that the full
    // path would have terminated at an endpoint (or rejected) can be answered here
    // instead — the fast path only short-circuits interior roots.
    const double inset_lower = lower * (1.0 + 1.0e-9);
    const double inset_upper = upper * (1.0 - 1.0e-9);
    if (std::isfinite(temperature_guess)
        && temperature_guess > inset_lower && temperature_guess < inset_upper) {
        const CaloricEval ev = caloric_eval_full_with_log_n(
            rho_total, temperature_guess, n_h, log_n_h);
        const double residual = ev.internal_energy - internal_energy;
        if (std::abs(residual) <= 2.0e-13*scale) {
            profile_note_inversion_initial_guess();
            profile_note_inversion_iterations(0);
            return InversionResult{temperature_guess, true, ev};
        }
        // The full solver's first iteration reproduces exactly this state: it
        // narrows one bracket end to the guess, then takes T - residual/C_V.
        double fast_lower = lower, fast_upper = upper;
        if (residual > 0.0) fast_upper = temperature_guess;
        else                fast_lower = temperature_guess;
        const double candidate = temperature_guess - residual/ev.heat_capacity;
        const double edge = 4.0*std::numeric_limits<double>::epsilon()
                          * std::max(1.0, std::abs(temperature_guess));
        const bool accepted =
               std::isfinite(candidate)
            && candidate > fast_lower + edge
            && candidate < fast_upper - edge
            && (fast_upper - fast_lower
                > 2.0e-13*std::max(1.0, std::abs(temperature_guess)))
            && std::abs(2.0*residual) < std::abs((upper - lower)*ev.heat_capacity)
            && candidate > inset_lower && candidate < inset_upper;
        if (accepted) {
            const CaloricEval next = caloric_eval_with_log_n(
                rho_total, candidate, n_h, log_n_h);
            if (std::abs(next.internal_energy - internal_energy) <= 2.0e-13*scale) {
                profile_note_inversion_iterations(1);
                return InversionResult{candidate, true, next};
            }
        }
    }

    // ---- safeguarded Newton / bisection (unchanged) ---------------------------
    // A supplied temperature is used two ways: as an exact-hit fast path here, and
    // (below) as the Newton seed over the table-wide bracket. Seeding the BRACKET
    // from the guess as well was tried and is a pessimisation — a ~6%-wide bracket
    // leaves Newton no room, so any slight overshoot is forced to bisect: 9.8
    // iterations and 7.5 bisections per call, versus 1.13 iterations and ZERO
    // bisections for the same guess used only as a seed. The wide bracket is a
    // safety net; it should stay wide.
    if (std::isfinite(temperature_guess)
        && temperature_guess > table.min_temperature()
        && temperature_guess < table.max_temperature()) {
        // NB: the guess evaluation is not a bracket-endpoint evaluation and is
        // deliberately not counted as one — the accepted fast path must still
        // report bracket_evaluations == 0.
        const CaloricEval guessed = caloric_eval_with_log_n(
            rho_total, temperature_guess, n_h, log_n_h);
        if (std::abs(guessed.internal_energy-internal_energy) <= 2.0e-13*scale) {
            profile_note_inversion_initial_guess();
            profile_note_inversion_iterations(0);
            return InversionResult{temperature_guess, true, guessed};
        }
    }

    {
        profile_note_inversion_bracket_evaluations(2);
        const CaloricEval ev_lower = caloric_eval_with_log_n(
            rho_total, lower, n_h, log_n_h);
        const CaloricEval ev_upper = caloric_eval_with_log_n(
            rho_total, upper, n_h, log_n_h);
        const double e_lower = ev_lower.internal_energy;
        const double e_upper = ev_upper.internal_energy;
        const double lower_tolerance = 64.0 * std::numeric_limits<double>::epsilon()
                                     * std::max(internal_energy, e_lower);
        const double upper_tolerance = 64.0 * std::numeric_limits<double>::epsilon()
                                     * std::max(internal_energy, e_upper);
        if (internal_energy < e_lower - lower_tolerance
            || internal_energy > e_upper + upper_tolerance) {
            std::ostringstream message;
            message << "EOS temperature inversion out of bounds: rho=" << rho_total
                    << " kg/m^3, e_int=" << internal_energy << " J/m^3; valid e_int=["
                    << e_lower << ',' << e_upper << "] over T=[" << lower << ','
                    << upper << "] K";
            throw std::out_of_range(message.str());
        }
        if (std::abs(internal_energy - e_lower) <= lower_tolerance) {
            profile_note_inversion_iterations(0);
            return InversionResult{lower, true, ev_lower};
        }
        if (std::abs(internal_energy - e_upper) <= upper_tolerance) {
            profile_note_inversion_iterations(0);
            return InversionResult{upper, true, ev_upper};
        }
    }

    double temperature = temperature_guess;
    if (!std::isfinite(temperature) || temperature <= lower || temperature >= upper) {
        temperature = internal_energy / (1.5 * n_h * eos_constants::k_b);
        temperature = std::max(lower, std::min(upper, temperature));
        if (temperature <= lower || temperature >= upper)
            temperature = 0.5 * (lower + upper);
    }

    double last_heat_capacity = std::numeric_limits<double>::quiet_NaN();
    double last_candidate = std::numeric_limits<double>::quiet_NaN();
    // Safeguarded Newton (Numerical Recipes "rtsafe"). Track the bracket width from
    // the PREVIOUS iteration so the slow-progress test can require each Newton step
    // to beat a bisection; that is what makes the bounce-forever case impossible
    // without also rejecting good steps.
    double previous_bracket_width = upper - lower;
    for (int iteration = 0; iteration < 128; ++iteration) {
        // Energy and heat capacity from ONE Saha evaluation. The residual test
        // below can still short-circuit before the capacity is needed, but the
        // extra work in that case is a handful of flops — no second Saha solve.
        const CaloricEval ev = caloric_eval_full_with_log_n(
            rho_total, temperature, n_h, log_n_h);
        const double value = ev.internal_energy;
        const double residual = value - internal_energy;
        if (std::abs(residual) <= 2.0e-13 * scale) {
            profile_note_inversion_iterations(static_cast<std::uint64_t>(iteration));
            return InversionResult{temperature, true, ev};
        }
        if (residual > 0.0) upper = temperature;
        else lower = temperature;

        if (upper - lower <= 2.0e-13 * std::max(1.0, std::abs(temperature))) {
            profile_note_inversion_iterations(static_cast<std::uint64_t>(iteration));
            return InversionResult{0.5 * (lower + upper), false, CaloricEval()};
        }

        const double heat_capacity = ev.heat_capacity;
        double candidate = temperature - residual / heat_capacity;
        const double bracket_width = upper - lower;
        // Reject a Newton step only if it leaves the bracket (up to a few ulp of the
        // endpoints — the roundoff bounce this guard was added for) or if it fails to
        // beat halving the previous bracket.
        //
        // The previous guard rejected any candidate inside the outer 10% of the
        // bracket on EITHER side. e_int(T) is solved over the whole table range, and
        // the root sits in the lowest ~0.04% of it (T ~ 7 kK in [3.2e3, 1e7] K), so
        // that test threw away essentially every early Newton step and forced ~10
        // bisections just to shrink the bracket enough for Newton to be allowed to
        // run. Measured on the h1600 ns=1000 column: 29.0 iterations and 25.8
        // bisections per call, 21.4M calls in 30 s of simulated time — the inversion
        // was ~100% of runtime across decode/rhs/projection/conduction.
        const double edge = 4.0 * std::numeric_limits<double>::epsilon()
                          * std::max(1.0, std::abs(temperature));
        const bool out_of_bracket = !std::isfinite(candidate)
                                 || candidate <= lower + edge
                                 || candidate >= upper - edge;
        // Slow-progress test: require the Newton step to beat halving the previous
        // bracket.
        const bool too_slow = !(std::abs(2.0 * residual)
                                < std::abs(previous_bracket_width * heat_capacity));
        if (out_of_bracket || too_slow) {
            profile_note_inversion_bisection();
            candidate = 0.5 * (lower + upper);
        }
        previous_bracket_width = bracket_width;
        last_heat_capacity = heat_capacity;
        last_candidate = candidate;
        temperature = candidate;
    }
    std::ostringstream message;
    message << "EOS temperature inversion did not converge: rho=" << rho_total
            << " kg/m^3, e_int=" << internal_energy << " J/m^3, T="
            << temperature << " K, bracket=[" << lower << ',' << upper
            << "] K, residual="
            << (equilibrium_internal_energy(rho_total, temperature)-internal_energy)
            << " J/m^3, Cv=" << last_heat_capacity << " J/m^3/K, candidate="
            << last_candidate << " K";
    throw std::runtime_error(message.str());
}

} // namespace

namespace {

struct DecodedCaloricResult {
    CaloricMixtureThermo thermo;
    CaloricEval eval;
};

DecodedCaloricResult decode_caloric_impl(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double temperature_guess, bool debug_clamp,
    CaloricState* caloric_at_temperature) {
    require_finite_positive(rho_i, "rho_i");
    require_finite_positive(rho_n, "rho_n");
    if (!std::isfinite(momentum_i) || !std::isfinite(momentum_n)
        || !std::isfinite(energy_i) || !std::isfinite(energy_n)
        || !std::isfinite(phi_of_this_state))
        throw std::domain_error("mixture state values must be finite");

    const double rho = require_finite_positive(rho_i + rho_n, "rho_total");
    const double kinetic_i = 0.5 * momentum_i * momentum_i / rho_i;
    const double kinetic_n = 0.5 * momentum_n * momentum_n / rho_n;
    const double internal_energy = energy_i + energy_n - kinetic_i - kinetic_n
                                 - rho * phi_of_this_state;
    require_finite_positive(internal_energy, "decoded internal_energy");
    const InversionResult inverted = invert_caloric(
        table, rho, internal_energy, temperature_guess, debug_clamp);
    const double temperature = inverted.temperature;
    const CaloricEval ev = inverted.has_state ? inverted.eval
                                              : caloric_eval(rho, temperature);
    if (caloric_at_temperature)
        *caloric_at_temperature = caloric_state_of(ev, temperature);
    const double n_h = ev.n_h;
    const double p_e = ev.x * n_h * eos_constants::k_b * temperature;
    const double p_i = 2.0 * p_e;
    const double p_n = (1.0 - ev.x) * n_h
                     * eos_constants::k_b * temperature;
    const CaloricMixtureThermo thermo{
        rho, temperature, ev.x, rho_i / rho, n_h, ev.x * n_h,
        (1.0 - ev.x) * n_h, p_i, p_n, p_e, internal_energy};
    return DecodedCaloricResult{thermo, ev};
}

} // namespace

CaloricMixtureThermo decode_equilibrium_caloric_mixture(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double temperature_guess, bool debug_clamp,
    CaloricState* caloric_at_temperature) {
    return decode_caloric_impl(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n,
        phi_of_this_state, temperature_guess, debug_clamp,
        caloric_at_temperature).thermo;
}

MixtureThermo decode_equilibrium_mixture(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double temperature_guess, bool debug_clamp,
    CaloricState* caloric_at_temperature) {
    const DecodedCaloricResult decoded = decode_caloric_impl(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n,
        phi_of_this_state, temperature_guess, debug_clamp,
        caloric_at_temperature);
    const CaloricMixtureThermo& th = decoded.thermo;
    const GammaState gamma = gamma_state_from_eval(
        table, decoded.eval, th.T, debug_clamp);
    return MixtureThermo{
        th.rho, th.T, th.x_eq, th.x_row, th.n_H, th.n_e, th.n_HI,
        th.p_i, th.p_n, th.p_e, gamma.gamma_sound, th.internal_energy};
}

namespace {

struct ProjectedRowsResult {
    ProjectedMixtureRows rows;
    CaloricMixtureThermo thermo;
    CaloricEval eval;
};

ProjectedRowsResult project_rows_impl(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double trace_fraction_floor,
    double temperature_guess, bool debug_clamp) {
    // In equilibrium single-fluid mode the two density rows are only a
    // predictor representation: this projection immediately reconstructs them
    // from the conserved total density and the Saha equilibrium fraction. A
    // tiny trace-row undershoot from the explicit hydro predictor therefore
    // must not abort an otherwise valid conservative state. Keep the strict
    // per-row positivity requirement in the ordinary EOS decoders; here the
    // physically relevant requirements are finite rows and positive total rho.
    if (!std::isfinite(rho_i) || !std::isfinite(rho_n)
        || !std::isfinite(momentum_i) || !std::isfinite(momentum_n)
        || !std::isfinite(energy_i) || !std::isfinite(energy_n)
        || !std::isfinite(phi_of_this_state))
        throw std::domain_error("mixture state values must be finite");
    if (!(trace_fraction_floor > 0.0) || !(trace_fraction_floor < 0.5)
        || !std::isfinite(trace_fraction_floor))
        throw std::domain_error("trace_fraction_floor must be finite and in (0,0.5)");

    const double rho = require_finite_positive(rho_i + rho_n, "rho_total");
    const double momentum = momentum_i + momentum_n;
    if (!std::isfinite(momentum))
        throw std::domain_error("total momentum must be finite");
    const double velocity = momentum / rho;
    const double total_energy = energy_i + energy_n;
    if (!std::isfinite(total_energy))
        throw std::domain_error("total energy must be finite");
    const double internal_energy = total_energy - 0.5 * rho * velocity * velocity
                                 - rho * phi_of_this_state;
    require_finite_positive(internal_energy, "projected internal_energy");

    const InversionResult inverted = invert_caloric(
        table, rho, internal_energy, temperature_guess, debug_clamp);
    const double temperature = inverted.temperature;
    const CaloricEval ev = inverted.has_state ? inverted.eval
                                              : caloric_eval(rho, temperature);
    const double n_h = ev.n_h;
    const double x_row = std::max(trace_fraction_floor,
                                 std::min(1.0 - trace_fraction_floor, ev.x));
    const double projected_rho_i = x_row * rho;
    const double projected_rho_n = (1.0 - x_row) * rho;
    const double p_e = ev.x * n_h * eos_constants::k_b * temperature;
    const double p_i = 2.0 * p_e;
    const double p_n = (1.0 - ev.x) * n_h
                     * eos_constants::k_b * temperature;
    const double common_specific = 0.5 * velocity * velocity + phi_of_this_state;
    const double projected_energy_i = 1.5 * p_i
                                    + ev.x * n_h * eos_constants::chi_h
                                    + projected_rho_i * common_specific;
    const double projected_energy_n = total_energy - projected_energy_i;
    const double projected_energy_e = 1.5 * p_e;
    const CaloricMixtureThermo thermo{
        rho, temperature, ev.x, x_row, n_h, ev.x*n_h,
        (1.0-ev.x)*n_h, p_i, p_n, p_e, internal_energy};
    const ProjectedMixtureRows rows{
        projected_rho_i, projected_rho_n, projected_rho_i*velocity,
        projected_rho_n*velocity, projected_energy_i, projected_energy_n,
        projected_energy_e, temperature};
    return ProjectedRowsResult{rows, thermo, ev};
}

} // namespace

ProjectedMixtureRows project_equilibrium_rows(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double trace_fraction_floor,
    double temperature_guess, bool debug_clamp) {
    return project_rows_impl(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n,
        phi_of_this_state, trace_fraction_floor, temperature_guess,
        debug_clamp).rows;
}

ProjectedMixture project_equilibrium_single_fluid(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double trace_fraction_floor,
    double temperature_guess, bool debug_clamp) {
    const ProjectedRowsResult projected = project_rows_impl(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n,
        phi_of_this_state, trace_fraction_floor, temperature_guess, debug_clamp);
    const CaloricMixtureThermo& th = projected.thermo;
    const GammaState gamma = gamma_state_from_eval(
        table, projected.eval, th.T, debug_clamp);
    const MixtureThermo thermo{
        th.rho, th.T, th.x_eq, th.x_row, th.n_H, th.n_e, th.n_HI,
        th.p_i, th.p_n, th.p_e, gamma.gamma_sound, th.internal_energy};
    const ProjectedMixtureRows& rows = projected.rows;
    return ProjectedMixture{
        rows.rho_i, rows.rho_n, rows.momentum_i, rows.momentum_n,
        rows.energy_i, rows.energy_n, rows.energy_e, thermo};
}

namespace {

ProjectedRowsResult pack_rows_impl(
    const EosGammaTable& table, double rho, double momentum,
    double total_energy, double temperature, double phi,
    double trace_fraction_floor, double residual_relative_tolerance,
    bool debug_clamp) {
    require_finite_positive(rho, "rho_total");
    require_finite_positive(temperature, "temperature");
    if (!std::isfinite(momentum) || !std::isfinite(total_energy)
        || !std::isfinite(phi))
        throw std::domain_error("known-temperature equilibrium inputs must be finite");
    if (!(trace_fraction_floor > 0.0) || !(trace_fraction_floor < 0.5)
        || !std::isfinite(trace_fraction_floor))
        throw std::domain_error("trace_fraction_floor must be finite and in (0,0.5)");
    if (!(residual_relative_tolerance >= 0.0)
        || !std::isfinite(residual_relative_tolerance))
        throw std::domain_error(
            "known-temperature residual tolerance must be finite and non-negative");

    const double velocity = momentum/rho;
    const double internal_authoritative = total_energy
        - 0.5*rho*velocity*velocity-rho*phi;
    require_finite_positive(internal_authoritative,
                            "known-temperature internal_energy");
    const CaloricEval ev = caloric_eval(rho, temperature);
    const double relative_residual = std::abs(
        ev.internal_energy-internal_authoritative)
        / std::max(internal_authoritative, std::numeric_limits<double>::min());
    if (relative_residual > residual_relative_tolerance) {
        std::ostringstream message;
        message << "known-temperature equilibrium residual exceeds tolerance: rho="
                << rho << ", T=" << temperature << ", relative_residual="
                << relative_residual << ", tolerance="
                << residual_relative_tolerance;
        throw std::runtime_error(message.str());
    }

    const double n_h = ev.n_h;
    table.require_n_h_in_bounds_from_log(n_h, ev.log_n_h, debug_clamp);
    const double x_row = std::max(trace_fraction_floor,
        std::min(1.0-trace_fraction_floor, ev.x));
    const double rho_i = x_row*rho;
    const double rho_n = rho-rho_i;
    const double p_e = ev.x*n_h*eos_constants::k_b*temperature;
    const double p_i = 2.0*p_e;
    const double p_n = (1.0-ev.x)*n_h*eos_constants::k_b*temperature;
    const double specific_carried = 0.5*velocity*velocity+phi;
    const double energy_i = 1.5*p_i+ev.x*n_h*eos_constants::chi_h
                          + rho_i*specific_carried;
    const double energy_n = total_energy-energy_i;
    const CaloricMixtureThermo thermo{
        rho, temperature, ev.x, x_row, n_h, ev.x*n_h,
        (1.0-ev.x)*n_h, p_i, p_n, p_e, internal_authoritative};
    const ProjectedMixtureRows rows{
        rho_i, rho_n, rho_i*velocity, rho_n*velocity,
        energy_i, energy_n, 1.5*p_e, temperature};
    return ProjectedRowsResult{rows, thermo, ev};
}

} // namespace

ProjectedMixtureRows pack_equilibrium_rows_from_known_temperature(
    const EosGammaTable& table, double rho, double momentum,
    double total_energy, double temperature, double phi,
    double trace_fraction_floor, double residual_relative_tolerance,
    bool debug_clamp) {
    return pack_rows_impl(
        table, rho, momentum, total_energy, temperature, phi,
        trace_fraction_floor, residual_relative_tolerance, debug_clamp).rows;
}

ProjectedMixture pack_equilibrium_from_known_temperature(
    const EosGammaTable& table, double rho, double momentum,
    double total_energy, double temperature, double phi,
    double trace_fraction_floor, double residual_relative_tolerance,
    bool debug_clamp) {
    const ProjectedRowsResult packed = pack_rows_impl(
        table, rho, momentum, total_energy, temperature, phi,
        trace_fraction_floor, residual_relative_tolerance, debug_clamp);
    const CaloricMixtureThermo& th = packed.thermo;
    const GammaState gamma = gamma_state_from_eval(
        table, packed.eval, temperature, debug_clamp);
    const MixtureThermo thermo{
        th.rho, th.T, th.x_eq, th.x_row, th.n_H, th.n_e, th.n_HI,
        th.p_i, th.p_n, th.p_e, gamma.gamma_sound, th.internal_energy};
    const ProjectedMixtureRows& rows = packed.rows;
    return ProjectedMixture{
        rows.rho_i, rows.rho_n, rows.momentum_i, rows.momentum_n,
        rows.energy_i, rows.energy_n, rows.energy_e, thermo};
}

namespace {

MixtureFaceState face_state_from_eval(
    const EosGammaTable& table, double rho_total, double velocity,
    double temperature, double phi_of_face, double trace_fraction_floor,
    bool debug_clamp, const CaloricEval& ev) {
    const GammaState gamma = gamma_state_from_eval(
        table, ev, temperature, debug_clamp);
    const double n_h = ev.n_h;
    const double x_row = std::max(trace_fraction_floor,
                                 std::min(1.0 - trace_fraction_floor, ev.x));
    const double rho_i = x_row * rho_total;
    const double rho_n = (1.0 - x_row) * rho_total;
    const double p_e = ev.x * n_h * eos_constants::k_b * temperature;
    const double p_i = 2.0 * p_e;
    const double p_n = (1.0 - ev.x) * n_h
                     * eos_constants::k_b * temperature;
    const double e_int = ev.internal_energy;
    const double specific = 0.5 * velocity * velocity + phi_of_face;
    const double energy_i = 1.5 * p_i + ev.x * n_h * eos_constants::chi_h
                          + rho_i * specific;
    const double energy_n = e_int + rho_total * specific - energy_i;
    const MixtureThermo thermo{
        rho_total, temperature, ev.x, x_row, n_h, ev.x*n_h,
        (1.0-ev.x)*n_h, p_i, p_n, p_e, gamma.gamma_sound, e_int};
    const ProjectedMixture conserved{
        rho_i, rho_n, rho_i*velocity, rho_n*velocity, energy_i, energy_n,
        1.5*p_e, thermo};
    const double p_total = p_i + p_n;
    const double sound_speed = std::sqrt(gamma.gamma_sound * p_total / rho_total);
    const double dx_d_t = saha_dx_d_temperature(ev.x, temperature);
    const double dp_d_t_rho = n_h * eos_constants::k_b
        * (1.0 + ev.x + temperature * dx_d_t);
    const double de_d_t_rho = caloric_capacity(n_h, ev.x, temperature);
    const double dp_deint_rho = dp_d_t_rho / de_d_t_rho;
    return MixtureFaceState{{rho_total, velocity, temperature}, conserved,
                            p_total, sound_speed, dp_deint_rho};
}

void validate_face_inputs(const EosGammaTable& table, double rho_total,
                          double velocity, double temperature,
                          double phi_of_face, double trace_fraction_floor) {
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    if (!std::isfinite(velocity) || !std::isfinite(phi_of_face))
        throw std::domain_error("mixture face velocity and potential must be finite");
    if (!(trace_fraction_floor > 0.0) || !(trace_fraction_floor < 0.5)
        || !std::isfinite(trace_fraction_floor))
        throw std::domain_error("trace_fraction_floor must be finite and in (0,0.5)");
    if (table.empty()) throw std::logic_error("Gamma1 table is empty");
}

} // namespace

MixtureFaceState equilibrium_mixture_face_state(
    const EosGammaTable& table, double rho_total, double velocity,
    double temperature, double phi_of_face, double trace_fraction_floor,
    bool debug_clamp) {
    validate_face_inputs(table, rho_total, velocity, temperature,
                         phi_of_face, trace_fraction_floor);
    const double n_h = rho_total / eos_constants::m_h;
    const EosCoordinates coordinates = eos_coordinates(n_h, temperature);
    table.require_n_h_in_bounds_from_log(n_h, coordinates.log_n_h, debug_clamp);
    return face_state_from_eval(
        table, rho_total, velocity, temperature, phi_of_face,
        trace_fraction_floor, debug_clamp,
        caloric_eval_from_coordinates(coordinates));
}

MixtureFaceState equilibrium_mixture_face_state_from_logs(
    const EosGammaTable& table, double log_rho_total, double velocity,
    double log_temperature, double phi_of_face, double trace_fraction_floor,
    bool debug_clamp) {
    if (!std::isfinite(log_rho_total) || !std::isfinite(log_temperature))
        throw std::domain_error("mixture face logarithms must be finite");
    const double rho_total = std::exp(log_rho_total);
    const double temperature = std::exp(log_temperature);
    validate_face_inputs(table, rho_total, velocity, temperature,
                         phi_of_face, trace_fraction_floor);
    const double n_h = rho_total / eos_constants::m_h;
    const double log_n_h = n_h_log(n_h);
    const EosCoordinates coordinates = eos_coordinates_with_logs(
        n_h, temperature, log_temperature, log_n_h);
    table.require_n_h_in_bounds_from_log(n_h, log_n_h, debug_clamp);
    return face_state_from_eval(
        table, rho_total, velocity, temperature, phi_of_face,
        trace_fraction_floor, debug_clamp,
        caloric_eval_from_coordinates(coordinates));
}

std::array<double, 7> equilibrium_mixture_flux(const MixtureFaceState& face) {
    const ProjectedMixture& u = face.conserved;
    const double v = face.primitive.velocity;
    return {{u.momentum_i,
             u.momentum_n,
             u.rho_i*v*v + u.thermo.p_i,
             u.rho_n*v*v + u.thermo.p_n,
             (u.energy_i + u.thermo.p_i)*v,
             (u.energy_n + u.thermo.p_n)*v,
             u.energy_e*v}};
}

EosGammaTable EosGammaTable::load(const std::string& path) {
    std::ifstream input(path.c_str());
    if (!input) throw std::runtime_error("Gamma1 table: cannot open '" + path + "'");

    std::map<std::string, std::string> metadata;
    std::vector<Row> rows;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::string text = trim(line);
        if (text.empty()) continue;
        if (text[0] == '#') {
            const std::string field = trim(text.substr(1));
            const std::size_t equals = field.find('=');
            if (equals == std::string::npos || equals == 0 || equals + 1 >= field.size())
                throw std::runtime_error("Gamma1 table: malformed metadata at " + path
                                         + ":" + std::to_string(line_number));
            const std::string key = trim(field.substr(0, equals));
            const std::string value = trim(field.substr(equals + 1));
            if (key.empty() || value.empty() || !metadata.emplace(key, value).second)
                throw std::runtime_error("Gamma1 table: duplicate or empty metadata key '"
                                         + key + "' at " + path + ":"
                                         + std::to_string(line_number));
            continue;
        }
        std::istringstream parser(text);
        Row row;
        std::string extra;
        if (!(parser >> row.log_t_10 >> row.log_n_10 >> row.gamma1) || (parser >> extra)) {
            throw std::runtime_error("Gamma1 table: expected exactly 3 numeric columns at "
                                     + path + ":" + std::to_string(line_number));
        }
        if (!std::isfinite(row.log_t_10) || !std::isfinite(row.log_n_10)
            || !(row.gamma1 > 0.0) || !std::isfinite(row.gamma1)) {
            throw std::runtime_error("Gamma1 table: non-finite axis or non-positive Gamma1 at "
                                     + path + ":" + std::to_string(line_number));
        }
        rows.push_back(row);
    }
    for (const auto& item : expected_metadata()) {
        const auto found = metadata.find(item.first);
        if (found == metadata.end())
            throw std::runtime_error("Gamma1 table: missing metadata key '" + item.first
                                     + "' in '" + path + "'");
        if (found->second != item.second)
            throw std::runtime_error("Gamma1 table: unsupported metadata " + item.first
                                     + "='" + found->second + "' in '" + path + "'");
    }
    const std::size_t declared_n_t = parse_dimension(metadata, "nT", path);
    const std::size_t declared_n_n = parse_dimension(metadata, "nN", path);
    if (metadata.size() != expected_metadata().size() + 2)
        throw std::runtime_error("Gamma1 table: unexpected metadata key in '" + path + "'");
    if (rows.empty()) throw std::runtime_error("Gamma1 table: no data rows in '" + path + "'");

    // The generator writes density-major blocks with temperature ascending.
    std::vector<double> log_t_10;
    for (const Row& row : rows) {
        if (!log_t_10.empty() && row.log_n_10 != rows.front().log_n_10) break;
        push_unique(log_t_10, row.log_t_10);
    }
    if (log_t_10.size() < 2)
        throw std::runtime_error("Gamma1 table: temperature axis needs at least 2 points");
    std::vector<double> log_n_10;
    for (std::size_t offset = 0; offset < rows.size(); offset += log_t_10.size())
        push_unique(log_n_10, rows[offset].log_n_10);
    if (log_n_10.size() < 2 || rows.size() != log_t_10.size() * log_n_10.size())
        throw std::runtime_error("Gamma1 table: incomplete rectangular grid in '" + path + "'");
    if (log_t_10.size() != declared_n_t || log_n_10.size() != declared_n_n)
        throw std::runtime_error("Gamma1 table: metadata dimensions disagree with data in '"
                                 + path + "'");
    if (!std::is_sorted(log_t_10.begin(), log_t_10.end())
        || std::adjacent_find(log_t_10.begin(), log_t_10.end(),
                              std::greater_equal<double>()) != log_t_10.end()
        || !std::is_sorted(log_n_10.begin(), log_n_10.end())
        || std::adjacent_find(log_n_10.begin(), log_n_10.end(),
                              std::greater_equal<double>()) != log_n_10.end()) {
        throw std::runtime_error("Gamma1 table: axes must be strictly increasing");
    }
    for (std::size_t i_n = 0; i_n < log_n_10.size(); ++i_n) {
        for (std::size_t i_t = 0; i_t < log_t_10.size(); ++i_t) {
            const Row& row = rows[i_n * log_t_10.size() + i_t];
            if (row.log_n_10 != log_n_10[i_n] || row.log_t_10 != log_t_10[i_t])
                throw std::runtime_error("Gamma1 table: inconsistent rectangular row order");
        }
    }

    EosGammaTable table;
    table.log_temperature_.reserve(log_t_10.size());
    table.log_n_h_.reserve(log_n_10.size());
    for (double value : log_t_10) table.log_temperature_.push_back(value * LN_10);
    for (double value : log_n_10) table.log_n_h_.push_back(value * LN_10);
    table.gamma1_.reserve(rows.size());
    for (const Row& row : rows) table.gamma1_.push_back(row.gamma1);
    return table;
}

double EosGammaTable::min_temperature() const {
    if (empty()) throw std::logic_error("Gamma1 table is empty");
    return std::exp(log_temperature_.front());
}

double EosGammaTable::temperature_at(std::size_t index) const {
    if (index >= log_temperature_.size())
        throw std::out_of_range("Gamma1 temperature-axis index out of range");
    return std::exp(log_temperature_[index]);
}

double EosGammaTable::n_h_at(std::size_t index) const {
    if (index >= log_n_h_.size())
        throw std::out_of_range("Gamma1 density-axis index out of range");
    return std::exp(log_n_h_[index]);
}
double EosGammaTable::max_temperature() const {
    if (empty()) throw std::logic_error("Gamma1 table is empty");
    return std::exp(log_temperature_.back());
}
double EosGammaTable::min_n_h() const {
    if (empty()) throw std::logic_error("Gamma1 table is empty");
    return std::exp(log_n_h_.front());
}
double EosGammaTable::max_n_h() const {
    if (empty()) throw std::logic_error("Gamma1 table is empty");
    return std::exp(log_n_h_.back());
}

bool EosGammaTable::contains(double temperature, double n_h) const noexcept {
    if (empty() || !(temperature > 0.0) || !(n_h > 0.0)
        || !std::isfinite(temperature) || !std::isfinite(n_h)) return false;
    const double log_t = std::log(temperature);
    const double log_n = std::log(n_h);
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
                           * std::max({1.0, std::abs(log_t), std::abs(log_n)});
    return log_t >= log_temperature_.front() - tolerance
        && log_t <= log_temperature_.back() + tolerance
        && log_n >= log_n_h_.front() - tolerance
        && log_n <= log_n_h_.back() + tolerance;
}

void EosGammaTable::require_n_h_in_bounds(double n_h, bool debug_clamp) const {
    if (empty()) throw std::logic_error("Gamma1 table is empty");
    require_finite_positive(n_h, "n_H");
    require_n_h_in_bounds_from_log(n_h, n_h_log(n_h), debug_clamp);
}

void EosGammaTable::require_n_h_in_bounds_from_log(
    double n_h, double log_n_h, bool debug_clamp) const {
    if (empty()) throw std::logic_error("Gamma1 table is empty");
    require_finite_positive(n_h, "n_H");
    if (!std::isfinite(log_n_h))
        throw std::domain_error("log(n_H) must be finite");
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
                           * std::max(1.0, std::abs(log_n_h));
    const bool outside = log_n_h < log_n_h_.front() - tolerance
                      || log_n_h > log_n_h_.back() + tolerance;
    if (outside && !debug_clamp) {
        std::ostringstream message;
        message << "EOS density out of bounds: rho=" << n_h * eos_constants::m_h
                << " kg/m^3, n_H=" << n_h << " m^-3; valid n_H=["
                << min_n_h() << ',' << max_n_h() << "] m^-3";
        throw std::out_of_range(message.str());
    }
}

double EosGammaTable::gamma1(double temperature, double n_h, bool debug_clamp) const {
    if (empty()) throw std::logic_error("Gamma1 table is empty");
    require_finite_positive(temperature, "temperature");
    require_finite_positive(n_h, "n_H");
    return gamma1_from_logs(temperature, n_h, temperature_log(temperature),
                            n_h_log(n_h), debug_clamp);
}

double EosGammaTable::gamma1_from_logs(
    double temperature, double n_h, double log_temperature, double log_n_h,
    bool debug_clamp) const {
    if (empty()) throw std::logic_error("Gamma1 table is empty");
    require_finite_positive(temperature, "temperature");
    require_finite_positive(n_h, "n_H");
    if (!std::isfinite(log_temperature) || !std::isfinite(log_n_h))
        throw std::domain_error("Gamma1 logarithmic coordinates must be finite");
    if (eos_operation_counting_enabled) ++eos_operation_counter().gamma1_queries;
    double log_t = log_temperature;
    double log_n = log_n_h;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
                           * std::max({1.0, std::abs(log_t), std::abs(log_n)});
    const bool outside = log_t < log_temperature_.front() - tolerance
                      || log_t > log_temperature_.back() + tolerance
                      || log_n < log_n_h_.front() - tolerance
                      || log_n > log_n_h_.back() + tolerance;
    if (outside && !debug_clamp) {
        std::ostringstream message;
        message << "Gamma1 table query out of bounds: T=" << temperature
                << " K, n_H=" << n_h << " m^-3; valid T=[" << min_temperature()
                << ',' << max_temperature() << "], n_H=[" << min_n_h() << ','
                << max_n_h() << ']';
        throw std::out_of_range(message.str());
    }
    log_t = std::max(log_temperature_.front(), std::min(log_temperature_.back(), log_t));
    log_n = std::max(log_n_h_.front(), std::min(log_n_h_.back(), log_n));
    const std::size_t i_t = lower_cell(log_temperature_, log_t);
    const std::size_t i_n = lower_cell(log_n_h_, log_n);
    const double f_t = (log_t - log_temperature_[i_t]) /
                       (log_temperature_[i_t + 1] - log_temperature_[i_t]);
    const double f_n = (log_n - log_n_h_[i_n]) /
                       (log_n_h_[i_n + 1] - log_n_h_[i_n]);
    const std::size_t n_t = log_temperature_.size();
    const double g00 = gamma1_[i_n * n_t + i_t];
    const double g01 = gamma1_[i_n * n_t + i_t + 1];
    const double g10 = gamma1_[(i_n + 1) * n_t + i_t];
    const double g11 = gamma1_[(i_n + 1) * n_t + i_t + 1];
    return (1.0 - f_n) * ((1.0 - f_t) * g00 + f_t * g01)
         + f_n * ((1.0 - f_t) * g10 + f_t * g11);
}

} // namespace chromosphere

