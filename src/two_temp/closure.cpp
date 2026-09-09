/*!
 * @file two_temp/closure.cpp
 * @brief EXPERIMENTAL two-temperature solver: thermodynamic closure, state
 *        decode/pack, and the algebraic face-state builder.
 * @ingroup two_temp_solver
 *
 * See two_temp/two_temp.hpp for the governing equations. Nothing here reaches
 * into the release solver or the two-fluid research solver; the only shared code
 * is the Grid (chromosphere.hpp) and the equilibrium Saha closure (eos.hpp).
 *
 * ## The identity this file is built on
 *
 * The release caloric closure is
 *
 *     e_int(rho, T) = (3/2)(1 + x) n_H k_B T + x n_H chi_H,
 *
 * whose "1" is precisely the HEAVY translational term (3/2) n_H k_B T — protons
 * and neutrals together number n_H regardless of x. Splitting the mixture into a
 * heavy pool and an electron pool that owns the ionization reservoir therefore
 * needs no new thermodynamics at all:
 *
 *     E_e(rho, T_e)   = e_int(rho, T_e) - (3/2) n_H k_B T_e
 *     dE_e/dT_e       = C_V^eff(rho, T_e) - (3/2) n_H k_B
 *     e_h(rho, T_i)   = (3/2) n_H k_B T_i
 *
 * Both electron quantities are taken as exactly those differences of the
 * PUBLIC release EOS functions, so the two solvers can never drift apart in
 * their thermodynamics, and at T_e = T_i the two pools re-sum to e_int bit for
 * bit up to one floating-point subtraction and addition.
 */

#include "two_temp/two_temp.hpp"
#include "parallel.hpp"
#include "profiling.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>

namespace chromosphere {

namespace {

/// Heavy translational internal energy (3/2) n_H k_B T [J m^-3]. n_H is the total
/// heavy-particle density rho/m_H: protons and neutral hydrogen together always
/// number n_H, independent of the ionization fraction.
inline double heavy_internal_energy(double rho, double temperature) {
    return 1.5*(rho/eos_constants::m_h)*eos_constants::k_b*temperature;
}

/// Heavy heat capacity at constant density, d e_h/dT_i = (3/2) n_H k_B
/// [J m^-3 K^-1]. Constant in T_i — the heavy pool carries no ionization energy.
inline double heavy_capacity(double rho) {
    return 1.5*(rho/eos_constants::m_h)*eos_constants::k_b;
}

inline arma::SizeMat state_size(const Grid& grid) {
    return arma::size(grid.ns, num_of_two_temp_eq);
}

} // namespace

// ============================================================================
// Electron caloric closure
// ============================================================================

double two_temp_electron_energy(double rho, double T_e) {
    return equilibrium_internal_energy(rho, T_e) - heavy_internal_energy(rho, T_e);
}

double two_temp_electron_capacity(double rho, double T_e) {
    return equilibrium_heat_capacity(rho, T_e) - heavy_capacity(rho);
}

double two_temp_density_from_pressure(const EosGammaTable& table, double pressure,
                                      double T_e, double T_i) {
    if (!(pressure > 0.0) || !(T_e > 0.0) || !(T_i > 0.0))
        throw std::domain_error(
            "two-temperature ghost pressure and temperatures must be positive");
    auto pressure_at = [&](double n_h) {
        const double x = saha_ionization_fraction_n_h(n_h, T_e);
        return n_h*eos_constants::k_b*(T_i + x*T_e);
    };
    double lower = table.min_n_h();
    double upper = table.max_n_h();
    if (pressure < pressure_at(lower) || pressure > pressure_at(upper))
        throw std::out_of_range(
            "two-temperature ghost pressure implies a density outside the EOS table");
    // Bisection in the log of n_H: the bracket is exact and the closure is
    // monotone, so 100 halvings take the relative width below 1e-30. This runs
    // once per boundary rung per step, so a closed form buys nothing.
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double middle = std::sqrt(lower*upper);
        if (pressure_at(middle) < pressure) lower = middle; else upper = middle;
        if (upper - lower <= 1.0e-14*upper) break;
    }
    return 0.5*(lower + upper)*eos_constants::m_h;
}

double two_temp_temperature_from_electron_energy(
    const EosGammaTable& table, double rho, double e_elec,
    double temperature_guess, bool* clamped) {
    if (clamped) *clamped = false;
    if (!(rho > 0.0) || !std::isfinite(rho))
        throw std::domain_error("two-temperature decode needs a positive density");
    if (!std::isfinite(e_elec))
        throw std::domain_error("two-temperature decode got a non-finite E_e");

    double lower = table.min_temperature();
    double upper = table.max_temperature();
    const double e_lo = two_temp_electron_energy(rho, lower);
    const double e_hi = two_temp_electron_energy(rho, upper);
    // E_e is strictly increasing in T_e (both (3/2) p_e and x n_H chi_H are), so
    // the table bracket is the exact bracket and anything outside it is a clamp.
    if (e_elec <= e_lo) { if (clamped) *clamped = true; return lower; }
    if (e_elec >= e_hi) { if (clamped) *clamped = true; return upper; }

    // Safeguarded Newton with GUARANTEED bisection progress (Numerical Recipes
    // `rtsafe`). E_e is strictly increasing, so `lower` always carries a negative
    // residual and `upper` a positive one, and the bracket is an invariant.
    //
    // The second safeguard test — "Newton is not at least halving the bracket" —
    // is the one that matters here and is not optional. Near full ionization the
    // caloric curve has a very sharp knee: at x -> 1 the capacity collapses to the
    // translational value (order 1e-6 J m^-3 K^-1) while just below the knee, in
    // the partial-ionization zone, dx/dT_e makes it larger by orders of magnitude.
    // A Newton step taken from the flat side therefore overshoots the knee by
    // thousands of kelvin, and an inversion safeguarded ONLY against leaving the
    // bracket oscillates between the two bracket ends forever: each iterate lands
    // strictly inside by ~1e-10 K, so the out-of-bracket test never fires and the
    // bracket shrinks by ~1e-10 K per pass. That is a real, observed failure — it
    // aborted a 4000 s run at rho = 8.26e-11, x = 0.999965 — not a hypothetical.
    double temperature =
        (std::isfinite(temperature_guess) && temperature_guess > lower
         && temperature_guess < upper)
        ? temperature_guess : std::sqrt(lower*upper);
    double residual = two_temp_electron_energy(rho, temperature) - e_elec;
    double capacity = two_temp_electron_capacity(rho, temperature);
    if (residual < 0.0) lower = temperature; else upper = temperature;
    double step_previous = upper - lower;
    double step = step_previous;
    for (int iteration = 0; iteration < 100; ++iteration) {
        const bool newton_leaves_bracket =
            ((temperature - upper)*capacity - residual)
          * ((temperature - lower)*capacity - residual) > 0.0;
        const double previous = temperature;
        if (!(capacity > 0.0) || newton_leaves_bracket
            || std::abs(2.0*residual) > std::abs(step_previous*capacity)) {
            step_previous = step;
            step = 0.5*(upper - lower);
            temperature = lower + step;
        } else {
            step_previous = step;
            step = residual/capacity;
            temperature -= step;
        }
        // Either the bracket collapsed to adjacent doubles or the step is below
        // tolerance: both are convergence, and the first cannot be improved on.
        if (previous == temperature
            || std::abs(step) <= 1.0e-12*temperature) return temperature;
        residual = two_temp_electron_energy(rho, temperature) - e_elec;
        capacity = two_temp_electron_capacity(rho, temperature);
        if (residual < 0.0) lower = temperature; else upper = temperature;
    }
    {
        char buffer[512];
        std::snprintf(buffer, sizeof buffer,
                      "two-temperature electron-energy inversion did not "
                      "converge: rho=%.17g E_e=%.17g T=%.17g lower=%.17g "
                      "upper=%.17g x=%.6g E_e(T)=%.17g cap=%.6g",
                      rho, e_elec, temperature, lower, upper,
                      saha_ionization_fraction_n_h(rho/eos_constants::m_h,
                                                   temperature),
                      two_temp_electron_energy(rho, temperature),
                      two_temp_electron_capacity(rho, temperature));
        throw std::runtime_error(buffer);
    }
}

// ============================================================================
// Decode and pack
// ============================================================================

TwoTempThermo two_temp_decode(const EosGammaTable& table, double rho,
                              double momentum, double energy, double e_elec,
                              double phi_of_this_state,
                              double temperature_guess,
                              bool* clamped_Te, bool* clamped_Ti) {
    if (clamped_Ti) *clamped_Ti = false;
    if (!(rho > 0.0) || !std::isfinite(rho))
        throw std::domain_error("two-temperature decode needs a positive density");
    const double n_h = rho/eos_constants::m_h;
    table.require_n_h_in_bounds(n_h);

    TwoTempThermo out{};
    out.rho = rho;
    out.n_H = n_h;
    out.T_e = two_temp_temperature_from_electron_energy(
        table, rho, e_elec, temperature_guess, clamped_Te);
    out.x = saha_ionization_fraction_n_h(n_h, out.T_e);
    out.n_e = out.x*n_h;
    out.n_HI = (1.0 - out.x)*n_h;
    out.p_e = out.n_e*eos_constants::k_b*out.T_e;
    // E_e is re-evaluated AT the decoded T_e rather than reusing the stored row,
    // so the reported electron pool is exactly on the Saha manifold and the two
    // pools sum to the stored E to round-off.
    out.e_elec = two_temp_electron_energy(rho, out.T_e);

    // The heavy internal energy is what the TOTAL energy row has left over.
    const double e_heavy = energy - out.e_elec
                         - 0.5*momentum*momentum/rho - rho*phi_of_this_state;
    const double capacity = heavy_capacity(rho);
    double T_i = e_heavy/capacity;
    const double t_min = table.min_temperature();
    const double t_max = table.max_temperature();
    if (!std::isfinite(T_i) || T_i < t_min) {
        if (clamped_Ti) *clamped_Ti = true;
        T_i = t_min;
    } else if (T_i > t_max) {
        if (clamped_Ti) *clamped_Ti = true;
        T_i = t_max;
    }
    out.T_i = T_i;
    out.e_heavy = capacity*T_i;
    out.p_i = n_h*eos_constants::k_b*T_i;
    out.p = out.p_e + out.p_i;
    out.sound_speed = std::sqrt((5.0/3.0)*out.p/rho);
    return out;
}

void two_temp_pack_ghost(const Grid& grid, Vec& ghost, double rho,
                         double velocity, double T_e, double T_i, double phi) {
    if (ghost.n_elem != num_of_two_temp_eq) ghost.zeros(num_of_two_temp_eq);
    grid.eos_gamma_table.require_n_h_in_bounds(rho/eos_constants::m_h);
    const double e_elec = two_temp_electron_energy(rho, T_e);
    const double e_heavy = heavy_internal_energy(rho, T_i);
    ghost(tt::RHO)    = static_cast<Real>(rho);
    ghost(tt::MOM)    = static_cast<Real>(rho*velocity);
    ghost(tt::ENERGY) = static_cast<Real>(
        two_temp_total_energy(rho, rho*velocity, e_heavy, e_elec, phi));
    ghost(tt::E_ELEC) = static_cast<Real>(e_elec);
}

void two_temp_pack_cell(const Grid& grid, Vec& state, arma::uword cell,
                        double rho, double velocity, double T_e, double T_i,
                        double phi) {
    Vec packed(num_of_two_temp_eq, arma::fill::zeros);
    two_temp_pack_ghost(grid, packed, rho, velocity, T_e, T_i, phi);
    const auto sz = state_size(grid);
    for (arma::uword row = 0; row < num_of_two_temp_eq; ++row)
        state(arma::sub2ind(sz, cell, row)) = packed(row);
}

// ============================================================================
// Face state
// ============================================================================

TwoTempThermo two_temp_face_state(const EosGammaTable& table, double rho,
                                  double velocity, double pressure,
                                  double T_e, bool* clamped) {
    if (clamped) *clamped = false;
    if (!(rho > 0.0) || !(pressure > 0.0) || !(T_e > 0.0)
        || !std::isfinite(velocity))
        throw std::domain_error("two-temperature face state is not admissible");
    const double n_h = rho/eos_constants::m_h;
    table.require_n_h_in_bounds(n_h);

    TwoTempThermo out{};
    out.rho = rho;
    out.n_H = n_h;
    out.T_e = T_e;
    out.x = saha_ionization_fraction_n_h(n_h, T_e);
    out.n_e = out.x*n_h;
    out.n_HI = (1.0 - out.x)*n_h;
    out.p = pressure;
    out.p_e = out.n_e*eos_constants::k_b*T_e;
    // The TOTAL pressure is the authoritative mechanical variable of the
    // reconstruction and is preserved exactly; if the independently limited
    // electron temperature would claim all of it, cap the ELECTRON share instead
    // and count the event. This is the same "guard, count, report" discipline the
    // release Godunov flux applies to its Rusanov fallbacks. x and n_e are left at
    // the pre-cap T_e, so a clamped face is admissible but not on-manifold and a
    // nonzero count invalidates the run.
    const double p_e_max = kFacePeFraction*pressure;
    if (out.p_e > p_e_max) {
        if (clamped) *clamped = true;
        out.p_e = p_e_max;
        out.T_e = (out.n_e > 0.0) ? out.p_e/(out.n_e*eos_constants::k_b) : T_e;
    }
    out.p_i = out.p - out.p_e;
    out.T_i = out.p_i/(n_h*eos_constants::k_b);
    out.e_elec = 1.5*out.p_e + out.x*n_h*eos_constants::chi_h;
    out.e_heavy = 1.5*out.p_i;
    out.sound_speed = std::sqrt((5.0/3.0)*pressure/rho);
    return out;
}

// ============================================================================
// Collision and conduction closures
// ============================================================================
//
// The two conductivities are the SAME expressions the release conduction operator
// uses (physics.hpp::physical_kappa_e / physical_kappa_n, transcribed here so
// src/two_temp/ depends only on chromosphere.hpp and eos.hpp, exactly as
// src/single_fluid/integrator.cpp does). The experiment changes only which
// temperature each channel is evaluated at, so kappa_e + kappa_n reproduces the
// release total conductivity exactly whenever T_e = T_i.

double two_temp_kappa_e(double n_e, double n_HI, double T_e) {
    return 9.2048e-12*n_e*std::pow(T_e, 2.5)
         / (n_e + 2.836e-11*n_HI*T_e*T_e);
}

double two_temp_kappa_i(double n_e, double n_HI, double T_i) {
    return 0.0342006*n_HI*T_i
         / (1.20613*n_e*std::sqrt(2.0*T_i)
            + 1.70573*n_HI*std::sqrt(T_i));
}

namespace {

/// Coulomb logarithm lnLambda_ei, NRL Plasma Formulary, Z = 1, clamped to
/// [5, 30]. Scalar transcription of physics.hpp::coulomb_log_ei.
double coulomb_log_ei_scalar(double n_e, double T_e) {
    constexpr double k_b_eV = 8.617333262e-5;   // eV/K
    const double ne_cm3 = std::max(n_e, 1.0)*1.0e-6;
    const double T_eV = k_b_eV*std::max(T_e, 1.0);
    const double sqrt_ne = std::sqrt(ne_cm3);
    const double value = (T_eV < 10.0)
        ? 23.0 - std::log(sqrt_ne*std::pow(T_eV, -1.5))
        : 24.0 - std::log(sqrt_ne/T_eV);
    return std::min(30.0, std::max(5.0, value));
}

} // namespace

double two_temp_exchange_conductance(double n_e, double n_HI, double T_e) {
    const double T = std::max(T_e, 1.0);
    // Electron-ion (Spitzer 1962 / NRL) energy-equilibration frequency
    // nu_ei = 2.030e-43 n_e lnLambda / (k_B T_e)^{3/2}, and the electron-neutral
    // one nu_en = (2 m_e/m_H) * 1.55e-15 n_HI sqrt(T_e), both scalar
    // transcriptions of physics.hpp. Both channels are kept because the lower
    // chromosphere is only weakly ionized (n_HI >> n_e), where the neutral channel
    // is the faster route to a common temperature.
    const double kT = eos_constants::k_b*T;
    const double nu_ei = 2.030e-43*n_e*coulomb_log_ei_scalar(n_e, T)
                       / std::pow(kT, 1.5);
    const double nu_en = (2.0*eos_constants::m_e/eos_constants::m_h)
                       * 1.55e-15*n_HI*std::sqrt(T);
    // The frequencies are defined against the TRANSLATIONAL electron heat
    // capacity (they come from dT_e/dt = nu (T_i - T_e)), so the conductance uses
    // (3/2) n_e k_B and NOT the ionization-inflated dE_e/dT_e.
    return 1.5*n_e*eos_constants::k_b*(nu_ei + nu_en);
}

// ============================================================================
// Field decode
// ============================================================================

void TwoTempField::resize(arma::uword ns_in) {
    ns = ns_in;
    cells.resize(ns_in);
    extended_primitive.set_size((ns_in + 4)*num_of_two_temp_eq);
}

void two_temp_decode_into(const Grid& grid, const Vec& state,
                          TwoTempField& out, TwoTempStats& stats,
                          const TwoTempField* previous) {
    ProfileScope timer(ProfileRegion::Decode);
    if (grid.eos_gamma_table.empty())
        throw std::logic_error("the two-temperature solver requires a Gamma1 table");
    if (state.n_elem != grid.ns*num_of_two_temp_eq)
        throw std::invalid_argument("two-temperature state size mismatch");

    out.resize(grid.ns);
    const arma::uword ext_n = grid.ns + 4;
    const auto ext_size = arma::size(ext_n, num_of_two_temp_eq);
    const auto sz = state_size(grid);
    auto put = [&](arma::uword ext_i, const TwoTempThermo& th, double momentum) {
        auto slot = [&](arma::uword k) { return arma::sub2ind(ext_size, ext_i, k); };
        out.extended_primitive(slot(TT_LOG_RHO)) = std::log(th.rho);
        out.extended_primitive(slot(TT_U))       = momentum/th.rho;
        out.extended_primitive(slot(TT_LOG_P))   = std::log(th.p);
        out.extended_primitive(slot(TT_LOG_TE))  = std::log(th.T_e);
    };

    const bool has_guesses = previous && previous->ns == grid.ns
        && previous->cells.size() == grid.ns;
    // Guard counters are accumulated per thread and merged once, so the shared
    // stats object is never written concurrently.
    std::array<std::uint64_t, kMaximumParallelThreads> te_clamps{};
    std::array<std::uint64_t, kMaximumParallelThreads> ti_clamps{};
    std::array<double, kMaximumParallelThreads> decoupling{};
    ParallelFailure failure;
    parallel_for_cells(grid.ns, [&](std::size_t raw_i) {
        const arma::uword i = static_cast<arma::uword>(raw_i);
        const std::size_t tid = static_cast<std::size_t>(parallel_thread_index());
        try {
            auto at = [&](arma::uword row) {
                return static_cast<double>(state(arma::sub2ind(sz, i, row)));
            };
            const double guess = has_guesses
                ? previous->cells[i].T_e
                : std::numeric_limits<double>::quiet_NaN();
            bool clamp_e = false, clamp_i = false;
            out.cells[i] = two_temp_decode(
                grid.eos_gamma_table, at(tt::RHO), at(tt::MOM), at(tt::ENERGY),
                at(tt::E_ELEC), two_temp_cell_phi(grid, i), guess,
                &clamp_e, &clamp_i);
            if (clamp_e) ++te_clamps[tid];
            if (clamp_i) ++ti_clamps[tid];
            const TwoTempThermo& th = out.cells[i];
            const double gap = std::abs(th.T_e - th.T_i)/th.T_i;
            decoupling[tid] = std::max(decoupling[tid], gap);
            put(i + 2, th, at(tt::MOM));
        } catch (...) {
            failure.capture(i, std::current_exception());
        }
    });
    failure.rethrow_lowest();
    for (std::size_t slot = 0; slot < kMaximumParallelThreads; ++slot) {
        stats.decode_Te_clamps += te_clamps[slot];
        stats.decode_Ti_clamps += ti_clamps[slot];
        stats.max_relative_decoupling =
            std::max(stats.max_relative_decoupling, decoupling[slot]);
    }

    auto decode_ghost = [&](const Vec& ghost, double phi, arma::uword ext_i) {
        bool clamp_e = false, clamp_i = false;
        const TwoTempThermo th = two_temp_decode(
            grid.eos_gamma_table, ghost(tt::RHO), ghost(tt::MOM),
            ghost(tt::ENERGY), ghost(tt::E_ELEC), phi,
            std::numeric_limits<double>::quiet_NaN(), &clamp_e, &clamp_i);
        if (clamp_e) ++stats.decode_Te_clamps;
        if (clamp_i) ++stats.decode_Ti_clamps;
        put(ext_i, th, static_cast<double>(ghost(tt::MOM)));
    };
    const double phi_inner = grid.phi_g_imh(0);
    const double phi_outer = grid.phi_g_iph(grid.ns - 1);
    decode_ghost(grid.tt_inner_boundary1, phi_inner, 0);
    decode_ghost(grid.tt_inner_boundary0, phi_inner, 1);
    decode_ghost(grid.tt_outer_boundary0, phi_outer, grid.ns + 2);
    decode_ghost(grid.tt_outer_boundary1, phi_outer, grid.ns + 3);
}

} // namespace chromosphere
