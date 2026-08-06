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
#include <cstdint>
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

/// Focused EOS operation counters used by regression tests and one-off profiling.
/// Counting is disabled by default so production hot paths pay only a predictable
/// false branch. These counters are not part of the runtime physics.
struct EosOperationCounts {
    std::uint64_t gamma1_queries = 0;
    std::uint64_t temperature_logs = 0;
    std::uint64_t n_h_logs = 0;
};

void set_eos_operation_counting(bool enabled) noexcept;
void reset_eos_operation_counts() noexcept;
EosOperationCounts eos_operation_counts() noexcept;

/// Textbook pure-hydrogen Saha fraction x(rho,T), evaluated in double precision.
/// Throws std::domain_error unless rho and T are positive and finite.
double saha_ionization_fraction(double rho_total, double temperature);

/// Same closure parameterized by hydrogen-nuclei density n_H [m^-3].
double saha_ionization_fraction_n_h(double n_h, double temperature);

/// Algebraic inverse of the pure-H Saha pressure closure. With
///   S(T) = (2 pi m_e k_B T / h^2)^{3/2} exp(-chi_H / k_B T),
///   x^2/(1-x) = S/n_H,  p = (1+x) n_H k_B T,
/// eliminating n_H gives x = sqrt(S/(S+y)) and n_H = y/(1+x) for y = p/(k_B T);
/// the returned total mass density is m_H n_H. Evaluated in the log domain so it
/// is exact in both the fully-neutral (y >> S) and fully-ionized (y << S) limits.
/// Throws std::domain_error unless pressure and temperature are positive/finite.
double equilibrium_density_from_pressure(double pressure, double temperature);

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

    /// Internal log-aware counterpart. Physical values are still validated;
    /// the supplied logarithm must already correspond to the same validated n_H.
    void require_n_h_in_bounds_from_log(double n_h, double log_n_h,
                                        bool debug_clamp = false) const;

    /// Bilinear interpolation in log(T),log(n_H).  Out-of-bounds is a hard
    /// error unless debug_clamp is explicitly true.
    double gamma1(double temperature, double n_h, bool debug_clamp = false) const;

    /// Internal log-aware counterpart used by hot EOS paths after validation.
    double gamma1_from_logs(double temperature, double n_h,
                            double log_temperature, double log_n_h,
                            bool debug_clamp = false) const;

private:
    std::vector<double> log_temperature_; // natural log
    std::vector<double> log_n_h_;         // natural log
    std::vector<double> gamma1_;          // density-major [i_n * nT + i_T]
};

/// Analytic equilibrium total internal-energy density, including ionization.
double equilibrium_internal_energy(double rho_total, double temperature);

/// Analytic effective heat capacity d e_int(rho,T)/dT at fixed density.
double equilibrium_heat_capacity(double rho_total, double temperature);

/// Every caloric quantity of the pure-H equilibrium closure at one (rho, T),
/// produced by a SINGLE Saha evaluation. Each field is computed with exactly the
/// expression of the standalone function it replaces, so a caller that used to
/// invoke saha_ionization_fraction_n_h, equilibrium_internal_energy and
/// equilibrium_heat_capacity separately gets bit-for-bit the same values.
struct CaloricState {
    double n_h;
    double x;
    double n_e;
    double n_hi;
    double pressure;
    double internal_energy;
    double heat_capacity;
};

/// Equilibrium thermodynamics that are valid without a Gamma1 sound-speed query.
/// This is deliberately separate from MixtureThermo, whose gamma1 field is always
/// valid by contract.
struct CaloricMixtureThermo {
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
    double internal_energy;
};

/// Conserved rows plus the solved temperature, for projection/packing callers
/// that do not consume Gamma1 or the full MixtureThermo object.
struct ProjectedMixtureRows {
    double rho_i;
    double rho_n;
    double momentum_i;
    double momentum_n;
    double energy_i;
    double energy_n;
    double energy_e;
    double temperature;
};

/// Fused caloric evaluation. Throws std::domain_error unless rho and T are
/// positive and finite (same guards as the functions it fuses).
CaloricState equilibrium_caloric_state(double rho_total, double temperature);

/// Log-aware counterpart for hot internal callers and exact equivalence tests.
/// The supplied logs must describe the same validated physical state.
CaloricState equilibrium_caloric_state_from_logs(
    double rho_total, double temperature,
    double log_temperature, double log_n_h);

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
/// When `caloric_at_temperature` is non-null it receives the caloric evaluation
/// AT the decoded temperature — the one the inversion already performed. A
/// nonlinear source stage whose first residual pass sits at exactly that
/// (rho, T) can therefore start without repeating the Saha solve.
MixtureThermo decode_equilibrium_mixture(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false,
    CaloricState* caloric_at_temperature = nullptr);

/// Caloric-only decode for internal source and predictor stages that do not use
/// acoustic Gamma1. The returned object contains no partially-valid field.
CaloricMixtureThermo decode_equilibrium_caloric_mixture(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false,
    CaloricState* caloric_at_temperature = nullptr);

/// Conservative equilibrium projection. It preserves total mass, momentum,
/// and E_I+E_N, thermalizes relative drift through the center-of-mass kinetic
/// energy, and uses trace_fraction_floor only for the carrier-row mass split.
ProjectedMixture project_equilibrium_single_fluid(
    const EosGammaTable& table, double rho_i, double rho_n,
    double momentum_i, double momentum_n, double energy_i, double energy_n,
    double phi_of_this_state, double trace_fraction_floor = 1.0e-8,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false);

/// Projection rows for the packed solver path, without a Gamma1 query.
ProjectedMixtureRows project_equilibrium_rows(
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

/// Known-temperature packing rows for the conduction solver. This performs one
/// authoritative fused caloric evaluation and does not query Gamma1.
ProjectedMixtureRows pack_equilibrium_rows_from_known_temperature(
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

/// Log-aware face construction for MUSCL, which already carries log(rho) and
/// log(T). The physical values are reconstructed once and the supplied log(T)
/// is reused by Saha and Gamma1 interpolation.
MixtureFaceState equilibrium_mixture_face_state_from_logs(
    const EosGammaTable& table, double log_rho_total, double velocity,
    double log_temperature, double phi_of_face,
    double trace_fraction_floor = 1.0e-8, bool debug_clamp = false);

/// Physical seven-row flux of a predecoded equilibrium face state.
std::array<double, 7> equilibrium_mixture_flux(const MixtureFaceState& face);

} // namespace chromosphere
