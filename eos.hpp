/*! \file eos.hpp
 *  Classical pure-H Saha utility and the Stage-2 CRASH Gamma1 table loader.
 *
 *  Stages 0--4 establish the thermodynamic inputs, read-only mixture decode,
 *  and explicit conservative projection. Integrator/flux activation begins at
 *  Stage 5 and requires equilibrium-manifold reconstruction.
 */
#pragma once

#include <cstddef>
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace chromosphere {

/// Single source of truth for the classical pure-hydrogen EOS constants (SI).
namespace eos_constants {
constexpr double k_b   = 1.380649e-23;
constexpr double m_e   = 9.1093837015e-31;
constexpr double m_h   = 1.6726219e-27;
constexpr double h     = 6.62607015e-34;
constexpr double chi_h = 2.179872361e-18;
} // namespace eos_constants

/// Textbook pure-hydrogen Saha fraction x(rho,T), evaluated in double precision.
/// Throws std::domain_error unless rho and T are positive and finite.
double saha_ionization_fraction(double rho_total, double temperature);

/// Same closure parameterized by hydrogen-nuclei density n_H [m^-3].
double saha_ionization_fraction_n_h(double n_h, double temperature);

/// Equilibrium indices at one total-density/temperature state. gamma_energy is
/// analytic from the same Saha caloric closure; gamma_sound comes from CRASH.
struct GammaState {
    double x_eq;
    double gamma_energy;
    double gamma_sound;
};

/// Read-only thermodynamic interpretation of the two stored carrier rows.
/// x_eq drives physics; x_row records the stored mass split only.
struct MixtureThermo {
    double rho;
    double T;
    double x_eq;
    double x_row;
    double n_H;
    double n_e;
    double n_HI;
    double p_i;
    double p_n;
    double p_e;
    double gamma1;
    double internal_energy;
};

/// Seven conserved carrier rows after a Stage-4 equilibrium projection.
/// E_E is diagnostic and is not included in the conserved E_I+E_N total.
struct ProjectedMixture {
    double rho_i;
    double rho_n;
    double momentum_i;
    double momentum_n;
    double energy_i;
    double energy_n;
    double energy_e;
    MixtureThermo thermo;
};

/// Primitive variables limited by the gamma-table MUSCL path. Density and
/// temperature are stored linearly here; reconstruction uses their logarithms.
struct MixturePrimitive {
    double rho;
    double velocity;
    double temperature;
};

/// Fully decoded, on-manifold face state. Flux and wave-speed evaluation use
/// these fields directly, avoiding a temperature inversion after reconstruction.
struct MixtureFaceState {
    MixturePrimitive primitive;
    ProjectedMixture conserved;
    double p_total;
    double sound_speed;
};

/// Rectangular log(T)-log(n_H) table of the equilibrium sound-speed index Gamma1.
/// An empty value is the default/off state.  Only the stripped three-column,
/// classical pure-H runtime format emitted by util/eos/tabulate_gamma.f90 is read.
class EosGammaTable {
public:
    EosGammaTable() = default;

    static EosGammaTable load(const std::string& path);

    bool empty() const noexcept { return gamma1_.empty(); }
    std::size_t temperature_size() const noexcept { return log_temperature_.size(); }
    std::size_t density_size() const noexcept { return log_n_h_.size(); }

    /// Physical axis coordinates, primarily for exhaustive validation tools.
    double temperature_at(std::size_t index) const;
    double n_h_at(std::size_t index) const;

    double min_temperature() const;
    double max_temperature() const;
    double min_n_h() const;
    double max_n_h() const;

    /// True only for positive finite values inside both table axes.
    bool contains(double temperature, double n_h) const noexcept;

    /// Enforce the production density domain. The explicit debug bypass is
    /// paired with gamma1(..., debug_clamp=true); production callers must not use it.
    void require_n_h_in_bounds(double n_h, bool debug_clamp = false) const;

    /// Bilinear interpolation in log(T),log(n_H).  Out-of-bounds is a hard
    /// error unless debug_clamp is explicitly true.
    double gamma1(double temperature, double n_h, bool debug_clamp = false) const;

private:
    std::vector<double> log_temperature_; // natural log
    std::vector<double> log_n_h_;         // natural log
    std::vector<double> gamma1_;          // density-major [i_n * nT + i_T]
};

/// Analytic equilibrium total internal-energy density, including ionization.
double equilibrium_internal_energy(double rho_total, double temperature);

/// Analytic effective heat capacity d e_int(rho,T)/dT at fixed density.
double equilibrium_heat_capacity(double rho_total, double temperature);

/// Analytic equilibrium pressure and the two thermodynamic indices.
GammaState gamma_state(const EosGammaTable& table, double rho_total,
                       double temperature, bool debug_clamp = false);

/// Invert the analytic caloric EOS on the table's temperature bracket. The
/// optional guess accelerates convergence but is never required for correctness.
double temperature_from_rho_eint(
    const EosGammaTable& table, double rho_total, double internal_energy,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false);

/// Read-only decode. It subtracts the two original kinetic energies and the
/// state-carried gravitational potential; it never rewrites carrier rows.
MixtureThermo decode_equilibrium_mixture(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false);

/// Conservative equilibrium projection. It preserves total mass, momentum,
/// and E_I+E_N, thermalizes relative drift through the center-of-mass kinetic
/// energy, and uses trace_fraction_floor only for the carrier-row mass split.
ProjectedMixture project_equilibrium_single_fluid(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double trace_fraction_floor = 1.0e-8,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false);

/// Pack an equilibrium state when a nonlinear source solve already knows T.
/// The supplied total energy remains authoritative and the neutral row is its
/// conservative remainder after constructing the charged row.
ProjectedMixture pack_equilibrium_from_known_temperature(
    const EosGammaTable& table, double rho_total, double total_momentum,
    double total_energy_authoritative, double temperature,
    double phi_of_this_state, double trace_fraction_floor,
    double residual_relative_tolerance = 2.0e-11,
    bool debug_clamp = false);

/// Construct one equilibrium face state algebraically from reconstructed
/// mixture variables. Both sides of an interface must pass the same face phi.
MixtureFaceState equilibrium_mixture_face_state(
    const EosGammaTable& table, double rho_total, double velocity,
    double temperature, double phi_of_face,
    double trace_fraction_floor = 1.0e-8, bool debug_clamp = false);

/// Physical seven-row flux of a predecoded equilibrium face state.
std::array<double, 7> equilibrium_mixture_flux(const MixtureFaceState& face);

} // namespace chromosphere
