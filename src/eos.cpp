#include "eos.hpp"
#include "profiling.hpp"

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

double require_finite_positive(double value, const char* name) {
    if (!(value > 0.0) || !std::isfinite(value))
        throw std::domain_error(std::string(name) + " must be positive and finite");
    return value;
}

double saha_dx_d_temperature(double x, double temperature) {
    const double theta = eos_constants::chi_h /
                         (eos_constants::k_b * temperature);
    return x * (1.0 - x) * (1.5 + theta) /
           ((2.0 - x) * temperature);
}

double equilibrium_heat_capacity_impl(double rho_total, double temperature) {
    const double n_h = rho_total / eos_constants::m_h;
    const double x = saha_ionization_fraction_n_h(n_h, temperature);
    const double dx_d_t = saha_dx_d_temperature(x, temperature);
    return n_h * (1.5 * eos_constants::k_b * (1.0 + x)
                  + (1.5 * eos_constants::k_b * temperature
                     + eos_constants::chi_h) * dx_d_t);
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

double saha_ionization_fraction_n_h(double n_h, double temperature) {
    require_finite_positive(n_h, "n_H");
    require_finite_positive(temperature, "temperature");

    // x^2/(1-x)=a. Work entirely through log(a), selecting algebraically
    // equivalent roots that do not overflow for either a >> 1 or a << 1.
    const double log_c = 1.5 * std::log(2.0 * PI * eos_constants::m_e
                                        * eos_constants::k_b /
                                        (eos_constants::h * eos_constants::h));
    const double log_a = log_c + 1.5 * std::log(temperature)
                       - eos_constants::chi_h / (eos_constants::k_b * temperature)
                       - std::log(n_h);
    if (log_a >= 0.0) {
        return 2.0 / (1.0 + std::sqrt(1.0 + 4.0 * std::exp(-log_a)));
    }
    const double u = std::exp(0.5 * log_a);
    return 2.0 * u / (u + std::sqrt(u * u + 4.0));
}

double saha_ionization_fraction(double rho_total, double temperature) {
    require_finite_positive(rho_total, "rho_total");
    return saha_ionization_fraction_n_h(rho_total / eos_constants::m_h, temperature);
}

double equilibrium_internal_energy(double rho_total, double temperature) {
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    const double n_h = rho_total / eos_constants::m_h;
    const double x = saha_ionization_fraction_n_h(n_h, temperature);
    const double pressure = (1.0 + x) * n_h * eos_constants::k_b * temperature;
    return 1.5 * pressure + x * n_h * eos_constants::chi_h;
}

double equilibrium_heat_capacity(double rho_total, double temperature) {
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    return equilibrium_heat_capacity_impl(rho_total, temperature);
}

GammaState gamma_state(const EosGammaTable& table, double rho_total,
                       double temperature, bool debug_clamp) {
    if (table.empty()) throw std::logic_error("Gamma1 table is empty");
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    const double n_h = rho_total / eos_constants::m_h;
    table.require_n_h_in_bounds(n_h, debug_clamp);
    const double x = saha_ionization_fraction_n_h(n_h, temperature);
    const double pressure = (1.0 + x) * n_h * eos_constants::k_b * temperature;
    const double internal_energy = 1.5 * pressure + x * n_h * eos_constants::chi_h;
    return GammaState{x, 1.0 + pressure / internal_energy,
                      table.gamma1(temperature, n_h, debug_clamp)};
}

double temperature_from_rho_eint(const EosGammaTable& table, double rho_total,
                                 double internal_energy,
                                 double temperature_guess, bool debug_clamp) {
    if (table.empty()) throw std::logic_error("Gamma1 table is empty");
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(internal_energy, "internal_energy");
    table.require_n_h_in_bounds(rho_total / eos_constants::m_h, debug_clamp);

    profile_note_inversion_call();
    const double scale = std::max(std::numeric_limits<double>::min(), internal_energy);

    double lower = table.min_temperature();
    double upper = table.max_temperature();

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
        const double guessed_energy = equilibrium_internal_energy(
            rho_total, temperature_guess);
        if (std::abs(guessed_energy-internal_energy) <= 2.0e-13*scale) {
            profile_note_inversion_initial_guess();
            profile_note_inversion_iterations(0);
            return temperature_guess;
        }
    }

    {
        profile_note_inversion_bracket_evaluations(2);
        const double e_lower = equilibrium_internal_energy(rho_total, lower);
        const double e_upper = equilibrium_internal_energy(rho_total, upper);
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
            return lower;
        }
        if (std::abs(internal_energy - e_upper) <= upper_tolerance) {
            profile_note_inversion_iterations(0);
            return upper;
        }
    }

    double temperature = temperature_guess;
    if (!std::isfinite(temperature) || temperature <= lower || temperature >= upper) {
        const double n_h = rho_total / eos_constants::m_h;
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
        const double value = equilibrium_internal_energy(rho_total, temperature);
        const double residual = value - internal_energy;
        if (std::abs(residual) <= 2.0e-13 * scale) {
            profile_note_inversion_iterations(static_cast<std::uint64_t>(iteration));
            return temperature;
        }
        if (residual > 0.0) upper = temperature;
        else lower = temperature;

        if (upper - lower <= 2.0e-13 * std::max(1.0, std::abs(temperature))) {
            profile_note_inversion_iterations(static_cast<std::uint64_t>(iteration));
            return 0.5 * (lower + upper);
        }

        const double heat_capacity = equilibrium_heat_capacity_impl(rho_total, temperature);
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

MixtureThermo decode_equilibrium_mixture(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double temperature_guess, bool debug_clamp) {
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
    const double temperature = temperature_from_rho_eint(
        table, rho, internal_energy, temperature_guess, debug_clamp);
    const GammaState gamma = gamma_state(table, rho, temperature, debug_clamp);
    const double n_h = rho / eos_constants::m_h;
    const double p_e = gamma.x_eq * n_h * eos_constants::k_b * temperature;
    const double p_i = 2.0 * p_e;
    const double p_n = (1.0 - gamma.x_eq) * n_h
                     * eos_constants::k_b * temperature;
    return MixtureThermo{rho, temperature, gamma.x_eq, rho_i / rho, n_h,
                         gamma.x_eq * n_h, (1.0 - gamma.x_eq) * n_h,
                         p_i, p_n, p_e, gamma.gamma_sound, internal_energy};
}

ProjectedMixture project_equilibrium_single_fluid(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double trace_fraction_floor,
    double temperature_guess, bool debug_clamp) {
    require_finite_positive(rho_i, "rho_i");
    require_finite_positive(rho_n, "rho_n");
    if (!std::isfinite(momentum_i) || !std::isfinite(momentum_n)
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

    const double temperature = temperature_from_rho_eint(
        table, rho, internal_energy, temperature_guess, debug_clamp);
    const GammaState gamma = gamma_state(table, rho, temperature, debug_clamp);
    const double n_h = rho / eos_constants::m_h;
    const double x_row = std::max(trace_fraction_floor,
                                 std::min(1.0 - trace_fraction_floor, gamma.x_eq));
    const double projected_rho_i = x_row * rho;
    const double projected_rho_n = (1.0 - x_row) * rho;
    const double p_e = gamma.x_eq * n_h * eos_constants::k_b * temperature;
    const double p_i = 2.0 * p_e;
    const double p_n = (1.0 - gamma.x_eq) * n_h
                     * eos_constants::k_b * temperature;
    const double common_specific = 0.5 * velocity * velocity + phi_of_this_state;
    const double projected_energy_i = 1.5 * p_i
                                    + gamma.x_eq * n_h * eos_constants::chi_h
                                    + projected_rho_i * common_specific;
    // Use the conserved total as the final row remainder. The temperature solve
    // makes this equal to the analytic neutral mapping to root tolerance, while
    // the remainder form preserves E_I+E_N in finite precision.
    const double projected_energy_n = total_energy - projected_energy_i;
    const double projected_energy_e = 1.5 * p_e;
    const MixtureThermo thermo{
        rho, temperature, gamma.x_eq, x_row, n_h, gamma.x_eq*n_h,
        (1.0-gamma.x_eq)*n_h, p_i, p_n, p_e, gamma.gamma_sound,
        internal_energy};
    return ProjectedMixture{
        projected_rho_i, projected_rho_n, projected_rho_i*velocity,
        projected_rho_n*velocity, projected_energy_i, projected_energy_n,
        projected_energy_e, thermo};
}

ProjectedMixture pack_equilibrium_from_known_temperature(
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
    const double internal_eos = equilibrium_internal_energy(rho, temperature);
    const double relative_residual = std::abs(internal_eos-internal_authoritative)
        / std::max(internal_authoritative, std::numeric_limits<double>::min());
    if (relative_residual > residual_relative_tolerance) {
        std::ostringstream message;
        message << "known-temperature equilibrium residual exceeds tolerance: rho="
                << rho << ", T=" << temperature << ", relative_residual="
                << relative_residual << ", tolerance="
                << residual_relative_tolerance;
        throw std::runtime_error(message.str());
    }

    const GammaState gamma = gamma_state(table, rho, temperature, debug_clamp);
    const double n_h = rho/eos_constants::m_h;
    const double x_row = std::max(trace_fraction_floor,
        std::min(1.0-trace_fraction_floor, gamma.x_eq));
    const double rho_i = x_row*rho;
    const double rho_n = rho-rho_i;
    const double p_e = gamma.x_eq*n_h*eos_constants::k_b*temperature;
    const double p_i = 2.0*p_e;
    const double p_n = (1.0-gamma.x_eq)*n_h*eos_constants::k_b*temperature;
    const double specific_carried = 0.5*velocity*velocity+phi;
    const double energy_i = 1.5*p_i+gamma.x_eq*n_h*eos_constants::chi_h
                          + rho_i*specific_carried;
    const double energy_n = total_energy-energy_i;
    const MixtureThermo thermo{
        rho, temperature, gamma.x_eq, x_row, n_h, gamma.x_eq*n_h,
        (1.0-gamma.x_eq)*n_h, p_i, p_n, p_e, gamma.gamma_sound,
        internal_authoritative};
    return ProjectedMixture{
        rho_i, rho_n, rho_i*velocity, rho_n*velocity,
        energy_i, energy_n, 1.5*p_e, thermo};
}

MixtureFaceState equilibrium_mixture_face_state(
    const EosGammaTable& table, double rho_total, double velocity,
    double temperature, double phi_of_face, double trace_fraction_floor,
    bool debug_clamp) {
    require_finite_positive(rho_total, "rho_total");
    require_finite_positive(temperature, "temperature");
    if (!std::isfinite(velocity) || !std::isfinite(phi_of_face))
        throw std::domain_error("mixture face velocity and potential must be finite");
    if (!(trace_fraction_floor > 0.0) || !(trace_fraction_floor < 0.5)
        || !std::isfinite(trace_fraction_floor))
        throw std::domain_error("trace_fraction_floor must be finite and in (0,0.5)");

    const GammaState gamma = gamma_state(table, rho_total, temperature, debug_clamp);
    const double n_h = rho_total / eos_constants::m_h;
    const double x_row = std::max(trace_fraction_floor,
                                 std::min(1.0 - trace_fraction_floor, gamma.x_eq));
    const double rho_i = x_row * rho_total;
    const double rho_n = (1.0 - x_row) * rho_total;
    const double p_e = gamma.x_eq * n_h * eos_constants::k_b * temperature;
    const double p_i = 2.0 * p_e;
    const double p_n = (1.0 - gamma.x_eq) * n_h
                     * eos_constants::k_b * temperature;
    const double e_int = equilibrium_internal_energy(rho_total, temperature);
    const double specific = 0.5 * velocity * velocity + phi_of_face;
    const double energy_i = 1.5 * p_i + gamma.x_eq * n_h * eos_constants::chi_h
                          + rho_i * specific;
    const double energy_n = e_int + rho_total * specific - energy_i;
    const MixtureThermo thermo{
        rho_total, temperature, gamma.x_eq, x_row, n_h, gamma.x_eq*n_h,
        (1.0-gamma.x_eq)*n_h, p_i, p_n, p_e, gamma.gamma_sound, e_int};
    const ProjectedMixture conserved{
        rho_i, rho_n, rho_i*velocity, rho_n*velocity, energy_i, energy_n,
        1.5*p_e, thermo};
    const double p_total = p_i + p_n;
    const double sound_speed = std::sqrt(gamma.gamma_sound * p_total / rho_total);
    return MixtureFaceState{{rho_total, velocity, temperature}, conserved,
                            p_total, sound_speed};
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
    const double log_n = std::log(n_h);
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon()
                           * std::max(1.0, std::abs(log_n));
    const bool outside = log_n < log_n_h_.front() - tolerance
                      || log_n > log_n_h_.back() + tolerance;
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
    double log_t = std::log(temperature);
    double log_n = std::log(n_h);
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

