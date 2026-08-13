/*! \file eos.hpp
 *  Classical pure-H Saha equilibrium closure and the CRASH Gamma1 table loader.
 *
 *  This is the thermodynamic closure of the release single-fluid field-aligned
 *  model. The authoritative state of one cell or face is the three-component
 *  mixture state (rho, rho u, E); every carrier quantity — the ionization
 *  fraction x, the electron density n_e, the neutral density n_HI, the electron
 *  partial pressure p_e — is DERIVED from it through Saha equilibrium and is a
 *  diagnostic, never an independently advanced variable.
 */
#pragma once

#include <cstddef>
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace chromosphere {

/** @addtogroup eos
 *  @{
 */

/// Single source of truth for the classical pure-hydrogen EOS constants (SI).
/// EosGammaTable::load() rejects any Gamma1 table whose metadata header does not
/// declare exactly these five values, so a loaded table and this closure are
/// guaranteed to have been built from the same physical constants.
namespace eos_constants {
constexpr double k_b   = 1.380649e-23;    ///< Boltzmann constant [J/K]
constexpr double m_e   = 9.1093837015e-31;///< electron rest mass [kg]
/// Mass per hydrogen nucleus [kg]. The pure-H mixture carries all of its mass in
/// the nuclei, so the hydrogen-nuclei density is n_H = rho_total / m_h.
constexpr double m_h   = 1.6726219e-27;
constexpr double h     = 6.62607015e-34;  ///< Planck constant [J s]
/// Hydrogen ground-state ionization potential chi_H [J] (13.5984 eV). It sets the
/// Saha exponential and is the per-ionization energy stored in e_int.
constexpr double chi_h = 2.179872361e-18;
} // namespace eos_constants

/// Focused EOS operation counters used by regression tests and one-off profiling.
/// Counting is disabled by default so production hot paths pay only a predictable
/// false branch. These counters are not part of the runtime physics.
struct EosOperationCounts {
    std::uint64_t gamma1_queries = 0;   ///< Bilinear Gamma1 table interpolations performed.
    std::uint64_t temperature_logs = 0; ///< Evaluations of log(T) inside the EOS.
    std::uint64_t n_h_logs = 0;         ///< Evaluations of log(n_H) inside the EOS.
};

/// Enable or disable EOS operation counting process-wide. Disabled by default;
/// the hot paths then only test one boolean. Not thread-scoped: the flag is
/// global while the counters themselves are per-thread.
void set_eos_operation_counting(bool enabled) noexcept;

/// Zero the per-thread EOS operation counters of every thread slot.
void reset_eos_operation_counts() noexcept;

/// Sum of the per-thread EOS operation counters over all thread slots.
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
    double x_eq;          ///< Saha equilibrium ionization fraction n_e/n_H [-]
    /// Caloric index gamma_E = 1 + p/e_int [-], with e_int the total internal
    /// energy density INCLUDING the ionization energy. This is the analytic index
    /// that relates pressure to internal energy; it is NOT the acoustic index.
    double gamma_energy;
    /// CRASH-tabulated adiabatic (acoustic) index Gamma1 [-], interpolated from
    /// EosGammaTable. It sets the sound speed c = sqrt(Gamma1 p / rho) and is a
    /// distinct quantity from gamma_energy.
    double gamma_sound;
};

/// Equilibrium thermodynamics of one mixture state, decoded from (rho, rho u, E).
/// Every carrier quantity below (x, n_e, n_HI, p_e) is derived from Saha
/// equilibrium at (rho, T); none of them is an independent conserved variable.
struct MixtureThermo {
    double rho;               ///< total mass density [kg/m^3] (the conserved row)
    double T;                 ///< common mixture temperature [K], from the caloric inversion
    double x;                 ///< Saha ionization fraction n_e/n_H [-]
    double n_H;               ///< hydrogen-nuclei density rho/m_H [m^-3]
    double n_e;               ///< electron number density = x n_H [m^-3]
    double n_HI;              ///< neutral-hydrogen number density = (1-x) n_H [m^-3]
    double p;                 ///< total pressure (1+x) n_H k_B T [Pa]
    double p_e;               ///< electron partial pressure x n_H k_B T [Pa]
    double gamma1;            ///< CRASH adiabatic (acoustic) index [-]
    double internal_energy;   ///< e_int(rho,T) including the ionization energy [J/m^3]
};

/// Caloric-only counterpart for stages that do not need the acoustic Gamma1.
/// The returned object contains no partially-valid field.
struct CaloricMixtureThermo {
    double rho;               ///< total mass density [kg/m^3] (the conserved row)
    double T;                 ///< common mixture temperature [K], from the caloric inversion
    double x;                 ///< Saha ionization fraction n_e/n_H [-]
    double n_H;               ///< hydrogen-nuclei density rho/m_H [m^-3]
    double n_e;               ///< electron number density = x n_H [m^-3]
    double n_HI;              ///< neutral-hydrogen number density = (1-x) n_H [m^-3]
    double p;                 ///< total pressure (1+x) n_H k_B T [Pa]
    double p_e;               ///< electron partial pressure x n_H k_B T [Pa]
    double internal_energy;   ///< e_int(rho,T) including the ionization energy [J/m^3]
};

/// Fully decoded, on-manifold face state of the release solver. Flux and
/// wave-speed evaluation use these fields directly, so no temperature inversion
/// is needed after the numerical flux.
struct MixtureFaceState {
    double rho;               ///< total mass density at the face [kg/m^3]
    double velocity;          ///< field-aligned velocity at the face [m/s]
    double temperature;       ///< common mixture temperature at the face [K]
    double pressure;          ///< total p [Pa]
    double internal_energy;   ///< e_int(rho,T) including ionization energy [J/m^3]
    double energy;            ///< E = e_int + 1/2 rho v^2 + rho phi [J/m^3]
    double sound_speed;       ///< sqrt(Gamma1 p / rho) [m/s]
    double gamma1;            ///< CRASH adiabatic (acoustic) index Gamma1 [-]
    /// General-EOS pressure response at fixed total density: the partial
    /// derivative of p with respect to e_int at fixed rho, b = (dp/de_int)_rho
    /// [dimensionless]. It is evaluated from the same Saha closure as the ratio
    /// (dp/dT)_rho / (de_int/dT)_rho. This is distinct from Gamma1: Gamma1 fixes
    /// the isentropic acoustic speed, while b fixes the energy component of the
    /// Euler characteristic vectors.
    double dp_deint_rho;
};

/// Rectangular log(T)-log(n_H) table of the equilibrium sound-speed index Gamma1.
/// An empty value is the default/off state.  Only the stripped three-column,
/// classical pure-H runtime format emitted by util/eos/tabulate_gamma.f90 is read.
class EosGammaTable {
public:
    EosGammaTable() = default;

    /// Read the runtime Gamma1 table from `path`. The file is a comment-header of
    /// `# key = value` metadata followed by density-major rows of
    /// (log10 T [K], log10 n_H [m^-3], Gamma1). Every metadata field — including
    /// the five physical constants and the two axis lengths — must match the
    /// values this build expects, and the data must form a complete rectangular
    /// grid with strictly increasing axes; anything else throws std::runtime_error.
    static EosGammaTable load(const std::string& path);

    /// True when no table has been loaded. An empty table is the default/off
    /// state, and every query method throws std::logic_error on it.
    bool empty() const noexcept { return gamma1_.empty(); }
    /// Number of points on the temperature axis.
    std::size_t temperature_size() const noexcept { return log_temperature_.size(); }
    /// Number of points on the hydrogen-nuclei-density axis.
    std::size_t density_size() const noexcept { return log_n_h_.size(); }

    /// Physical axis coordinates, primarily for exhaustive validation tools.
    /// Temperature [K] at temperature-axis index; throws std::out_of_range.
    double temperature_at(std::size_t index) const;
    /// Hydrogen-nuclei density [m^-3] at density-axis index; throws std::out_of_range.
    double n_h_at(std::size_t index) const;

    /// Lowest tabulated temperature [K]. Also the lower bracket end of the
    /// caloric temperature inversion.
    double min_temperature() const;
    /// Highest tabulated temperature [K]. Also the upper bracket end of the
    /// caloric temperature inversion.
    double max_temperature() const;
    /// Lowest tabulated hydrogen-nuclei density [m^-3].
    double min_n_h() const;
    /// Highest tabulated hydrogen-nuclei density [m^-3].
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
    /// Temperature axis as NATURAL log(T [K]), strictly increasing. The file
    /// stores log10; load() converts once so every hot query stays in natural logs.
    std::vector<double> log_temperature_;
    /// Density axis as NATURAL log(n_H [m^-3]), strictly increasing.
    std::vector<double> log_n_h_;
    /// Gamma1 values [-] on the rectangular axis product, density-major:
    /// element (i_n, i_T) lives at index i_n * temperature_size() + i_T.
    std::vector<double> gamma1_;
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
    double n_h;               ///< hydrogen-nuclei density rho/m_H [m^-3]
    double x;                 ///< Saha ionization fraction n_e/n_H [-]
    double n_e;               ///< electron number density = x n_H [m^-3]
    double n_hi;              ///< neutral-hydrogen number density = (1-x) n_H [m^-3]
    double pressure;          ///< total pressure (1+x) n_H k_B T [Pa]
    /// Total internal energy density e_int = 3p/2 + x n_H chi_H [J/m^3]; the
    /// second term is the ionization energy stored in the mixture.
    double internal_energy;
    /// Effective heat capacity at constant density C_V^eff = d e_int/dT [J/m^3/K].
    /// It includes the d x/dT ionization contribution, which dominates across the
    /// ionization zone, and is the Newton derivative of the caloric inversion.
    double heat_capacity;
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

/// Invert the equilibrium pressure closure p(rho,T) = (1+x(rho,T)) n_H k_B T for
/// the temperature at a KNOWN total density. This is the counterpart of
/// equilibrium_density_from_pressure and uses exactly the same authoritative
/// pressure as MixtureFaceState::pressure.
///
/// The closure supplies its own exact bracket: with T1 = p m_H/(rho k_B) the
/// fully-neutral temperature, x in (0,1) gives T in (T1/2, T1] for every physical
/// state, and dp/dT|_rho = n_H k_B (1 + x + T dx/dT) > 0 makes p(T) strictly
/// monotone there. The solver is therefore a safeguarded Newton iteration that can
/// never leave a valid factor-of-two bracket and falls back to bisection on any
/// rejected step. Throws std::domain_error unless rho and p are positive/finite,
/// and std::runtime_error if the iteration fails to converge.
double equilibrium_temperature_from_density_pressure(
    double rho_total, double pressure,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN());

/// Invert the analytic caloric EOS on the table's temperature bracket. The
/// optional guess accelerates convergence but is never required for correctness.
double temperature_from_rho_eint(
    const EosGammaTable& table, double rho_total, double internal_energy,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false);

/// Total energy of a mixture state, E = e_int + 1/2 rho v^2 + rho phi.
inline double mixture_total_energy(double rho, double momentum,
                                   double internal_energy, double phi) {
    return internal_energy + 0.5 * momentum * momentum / rho + rho * phi;
}

/// Read-only decode of one conserved mixture state (rho, rho u, E). It subtracts
/// the kinetic energy and the state-carried gravitational potential, inverts the
/// caloric EOS for T, and reports the derived Saha carrier quantities.
/// When `caloric_at_temperature` is non-null it receives the caloric evaluation
/// AT the decoded temperature — the one the inversion already performed. A
/// nonlinear source stage whose first residual pass sits at exactly that
/// (rho, T) can therefore start without repeating the Saha solve.
MixtureThermo decode_equilibrium_mixture(
    const EosGammaTable& table, double rho, double momentum, double energy,
    double phi_of_this_state,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false,
    CaloricState* caloric_at_temperature = nullptr);

/// Caloric-only decode for internal source and predictor stages that do not use
/// acoustic Gamma1. The returned object contains no partially-valid field.
CaloricMixtureThermo decode_equilibrium_caloric_mixture(
    const EosGammaTable& table, double rho, double momentum, double energy,
    double phi_of_this_state,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN(),
    bool debug_clamp = false,
    CaloricState* caloric_at_temperature = nullptr);

/// Construct one equilibrium face state algebraically from reconstructed
/// mixture variables. Both sides of an interface must pass the same face phi.
MixtureFaceState equilibrium_mixture_face_state(
    const EosGammaTable& table, double rho_total, double velocity,
    double temperature, double phi_of_face, bool debug_clamp = false);

/// Log-aware face construction for MUSCL, which already carries log(rho) and
/// log(T). The physical values are reconstructed once and the supplied log(T)
/// is reused by Saha and Gamma1 interpolation.
MixtureFaceState equilibrium_mixture_face_state_from_logs(
    const EosGammaTable& table, double log_rho_total, double velocity,
    double log_temperature, double phi_of_face, bool debug_clamp = false);

/// Pressure-based counterpart of equilibrium_mixture_face_state_from_logs, used
/// by the (log rho, V, log p) MUSCL reconstruction. The face temperature is
/// recovered by equilibrium_temperature_from_density_pressure and every remaining
/// face quantity is then produced by the SAME authoritative closure the
/// temperature-based builder uses, so the two differ only in which pair of
/// reconstructed variables is authoritative.
///
/// `temperature_guess` is an OPTIONAL initial guess for that inversion and
/// nothing else. It is accepted only if it lies strictly inside the exact
/// physical bracket (T_neutral/2, T_neutral]; otherwise the generic
/// 0.75*T_neutral start is used. The safeguarded Newton/bisection still
/// converges to the same root to the same tolerance, so the returned face state
/// is independent of the hint.
MixtureFaceState equilibrium_mixture_face_state_from_log_pressure(
    const EosGammaTable& table, double log_rho_total, double velocity,
    double log_pressure, double phi_of_face, bool debug_clamp = false,
    double temperature_guess = std::numeric_limits<double>::quiet_NaN());

/// Physical field-aligned flux (rho u, rho u^2 + p, (E + p) u) of a predecoded
/// equilibrium face state, in the (rho, rho u, E) row order of the release state.
std::array<double, 3> equilibrium_mixture_flux(const MixtureFaceState& face);

/** @} */

} // namespace chromosphere
