/*!
 * @file two_temp/two_temp.hpp
 * @brief EXPERIMENTAL two-temperature (T_e != T_i) field-aligned solver.
 *        Conserved state U = (rho, rho u, E, E_e).
 * @ingroup two_temp_solver
 *
 * two_temp/two_temp.hpp — the EXPERIMENTAL two-temperature solver.
 *
 * **Status: EXPERIMENTAL. This is not the release path and never runs from a
 * `model_column` command line.** It is a third, separate solver, selected only by
 * the `model_column_2t` scenario, and it exists to answer one question: how far
 * apart do the electron and heavy-particle temperatures actually drift in the
 * modelled 1600–2153 km column, where the release assumes a single temperature?
 * The release solver in `single_fluid/` is untouched by this file: no code in
 * `src/single_fluid/` reads any two-temperature `Grid` member, so the release
 * arithmetic is unaffected by this solver's presence. (The reverse edge exists
 * and is deliberate: `two_temp/hydro.cpp` reuses `single_fluid/exact_rs.hpp`,
 * the standalone port of the SWMF ideal-gas Riemann solver.)
 *
 * ## Physical model
 *
 * One fluid, one velocity, TWO temperatures, on the same prescribed flux tube of
 * cross-section A(s) with A(s) B(s) = const. Electrons are massless and comoving
 * (quasi-neutrality plus zero field-aligned current), so there is still exactly
 * one momentum equation. Protons and neutral hydrogen are lumped into a single
 * HEAVY-particle population at temperature T_h (written `T_i` in the outputs and
 * the prose, following the usual convention) — they exchange energy through
 * charge exchange and elastic collisions on a timescale far shorter than
 * anything resolved here, and the release model has no ion/neutral drift either.
 *
 *   d_t(rho)     + (1/A) d_s(A rho u)            = 0
 *   d_t(rho u)   + (1/A) d_s(A (rho u^2 + p))    = p d_s(ln A) - rho d_s(phi_g)
 *   d_t(E)       + (1/A) d_s(A (E + p) u)        = -(1/A) d_s(A (q_e + q_h))
 *   d_t(E_e)     + (1/A) d_s(A E_e u)            = -p_e (1/A) d_s(A u)
 *                                                  -(1/A) d_s(A q_e) + Q_ei
 *
 * ## Closure
 *
 *   n_H  = rho / m_H                 (heavy-particle density: protons + neutrals)
 *   x    = x_Saha(n_H, T_e)          (ionization equilibrium is an ELECTRON process)
 *   n_e  = x n_H,   n_HI = (1 - x) n_H
 *   p_e  = n_e k_B T_e               (electron partial pressure)
 *   p_h  = n_H k_B T_h               (heavy partial pressure; n_p + n_HI = n_H)
 *   p    = p_e + p_h = n_H k_B (T_h + x T_e)
 *   E_e  = (3/2) p_e + x n_H chi_H   (electron pool OWNS the ionization reservoir)
 *   E    = (3/2) p_h + E_e + 1/2 rho u^2 + rho phi_g
 *   q_e  = -kappa_e(n_e, n_HI, T_e) d_s T_e     (Spitzer electron channel)
 *   q_h  = -kappa_n(n_e, n_HI, T_h) d_s T_h     (neutral-hydrogen channel)
 *   Q_ei = (3/2) n_e k_B (nu_ei + nu_en) (T_h - T_e)
 *
 * At T_e = T_h = T this degenerates EXACTLY to the release closure:
 * p = (1+x) n_H k_B T and E_e + (3/2) p_h = (3/2) p + x n_H chi_H = e_int.
 *
 * ## Why this particular split, and why the total-energy row is kept
 *
 * `E` retains its release meaning — the TOTAL energy — so mass, momentum and
 * total energy stay in exact conservation form and the release numerical flux
 * applies to those three rows unchanged. The heavy internal energy is then a
 * DERIVED quantity, e_h = E - E_e - 1/2 rho u^2 - rho phi_g, which makes
 * T_h = (2/3) e_h / (n_H k_B) purely algebraic — no inversion at all. This
 * "total energy + electron energy" split is the standard two-temperature MHD
 * formulation (AWSoM/BATS-R-US; van der Holst et al. 2010, Sokolov et al. 2021)
 * and is the same one this project's historical two-fluid `cons::E_E` row already
 * uses (see chromosphere.hpp).
 *
 * The ionization energy is placed in the ELECTRON pool because hydrogen
 * ionization and recombination in this regime are electron-impact processes: the
 * energy chi_H per ionization is taken from, and returned to, the electron
 * thermal reservoir. That is also why Saha equilibrium is evaluated at T_e.
 *
 * ## Acoustics: unchanged from the release, and that is a result, not an accident
 *
 * Both species are monatomic, so at frozen composition and frozen energy
 * partition the TOTAL internal energy still splits exactly as
 *
 *     e_int = e_h + e_e = (3/2) p_h + (3/2) p_e + e_ion = p/(5/3 - 1) + e_ion,
 *
 * independent of HOW the thermal energy is shared between electrons and heavies.
 * The release SWMF-style exact-Riemann Godunov flux at gamma = 5/3 with the
 * ionization energy in the passive offset E0 is therefore exactly as valid here
 * as it is for the release, the frozen sound speed is still sqrt(5/3 p/rho), and
 * the acoustic CFL condition is unchanged. The extra physical content of the 2T
 * system rides on a SECOND characteristic at lambda = u: the eigenvalues of the
 * 4x4 system are {u - c, u, u, u + c}, the doubled `u` carrying the total entropy
 * and the electron/heavy energy PARTITION. That is what fixes the boundary
 * conditions — see two_temp/conduction.cpp and the study recap.
 *
 * ## What is NOT here
 *
 * No finite-rate ionization, no radiative cooling, no TRAC, no beam, no
 * volumetric heating, no ion/neutral drift, no anisotropic ion temperature. The
 * timestep is exactly
 *
 *   U^n -> explicit MUSCL-Hancock/Godunov hydro -> U* -> implicit coupled
 *          (T_e, T_h) conduction + collisional exchange -> U^{n+1}
 *
 * i.e. the release's two-stage Lie split, with the single scalar conduction solve
 * replaced by a 2x2-block one.
 */
#pragma once

#include "chromosphere.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace chromosphere {

/** @defgroup two_temp_solver Experimental two-temperature solver
 *  @brief Separate electron and heavy-particle temperatures on the release column.
 *  @{
 */

// ============================================================================
// State indices
// ============================================================================

/// Conserved rows of the experimental two-temperature state
/// U = (rho, rho u, E, E_e).
///
///   RHO      total mass density                                   [kg m^-3]
///   MOM      field-aligned momentum density rho u                 [kg m^-2 s^-1]
///   ENERGY   TOTAL energy density
///            E = (3/2) p_h + E_e + 1/2 rho u^2 + rho phi_g        [J m^-3]
///   E_ELEC   electron internal energy INCLUDING the ionization
///            reservoir, E_e = (3/2) p_e + x n_H chi_H             [J m^-3]
///
/// The first three rows carry exactly the release meanings, which is what lets
/// the release numerical flux be reused on them verbatim.
namespace tt {
    const arma::uword RHO    = 0;  ///< total mass density rho                 [kg m^-3]
    const arma::uword MOM    = 1;  ///< field-aligned momentum density rho u   [kg m^-2 s^-1]
    const arma::uword ENERGY = 2;  ///< TOTAL energy density E                 [J m^-3]
    const arma::uword E_ELEC = 3;  ///< electron internal energy E_e           [J m^-3]
}
/// Number of conserved rows in the two-temperature state; the packed state length
/// is `ns * num_of_two_temp_eq`. Deliberately different from both
/// `num_of_mixture_eq` (3) and `num_of_eq` (7): the three solvers do not share a
/// state width.
const arma::uword num_of_two_temp_eq = 4;

// ============================================================================
// Thermodynamics
// ============================================================================

/// Full two-temperature thermodynamic decode of one cell, ghost or face state.
/// Every carrier quantity is DERIVED from (rho, T_e); none is independently
/// advanced. `T_i` is the heavy-particle temperature shared by protons and
/// neutral hydrogen.
struct TwoTempThermo {
    double rho;      ///< total mass density [kg m^-3]
    double T_e;      ///< electron temperature [K]
    double T_i;      ///< heavy-particle (proton + neutral) temperature [K]
    double x;        ///< Saha ionization fraction n_e/n_H at T_e [-]
    double n_H;      ///< heavy-particle density rho/m_H [m^-3]
    double n_e;      ///< electron density x n_H [m^-3]
    double n_HI;     ///< neutral-hydrogen density (1-x) n_H [m^-3]
    double p_e;      ///< electron partial pressure n_e k_B T_e [Pa]
    double p_i;      ///< heavy partial pressure n_H k_B T_i [Pa]
    double p;        ///< total pressure p_e + p_i [Pa]
    double e_elec;   ///< electron internal energy (3/2)p_e + x n_H chi_H [J m^-3]
    double e_heavy;  ///< heavy internal energy (3/2) p_i [J m^-3]
    double sound_speed;  ///< frozen-composition speed sqrt(5/3 p/rho) [m s^-1]
};

/// Electron internal energy density at one (rho, T_e), including the ionization
/// reservoir: E_e = (3/2) x n_H k_B T_e + x n_H chi_H [J m^-3]. Strictly
/// increasing in T_e at fixed rho, which is what makes the decode inversion below
/// unconditionally bracketable.
double two_temp_electron_energy(double rho, double T_e);

/// Analytic derivative dE_e/dT_e at fixed rho [J m^-3 K^-1]. It contains the
/// dx/dT_e ionization contribution, which dominates across the partial-ionization
/// zone and is the Newton derivative of both the decode inversion and the
/// implicit electron conduction row.
double two_temp_electron_capacity(double rho, double T_e);

/// Invert the two-temperature pressure closure
/// p = n_H k_B (T_i + x(n_H, T_e) T_e) for the total mass density m_H n_H at a
/// known pair of temperatures. p is strictly increasing in n_H (the heavy term
/// linearly, and n_e = x n_H monotonically through Saha), so the table's density
/// axis is an exact bracket and a safeguarded bisection is unconditionally
/// convergent. Used by the boundary closure to build a ghost density from an
/// imposed pressure and two extrapolated temperatures; a request outside the
/// table's density range throws.
double two_temp_density_from_pressure(const EosGammaTable& table, double pressure,
                                      double T_e, double T_i);

/// Invert `two_temp_electron_energy` for T_e on the Gamma1 table's temperature
/// bracket. Safeguarded Newton with bisection fallback; the optional guess only
/// accelerates it. Values outside the bracket are CLAMPED to it and reported
/// through `clamped` so the event is counted rather than silent.
double two_temp_temperature_from_electron_energy(
    const EosGammaTable& table, double rho, double e_elec,
    double temperature_guess, bool* clamped = nullptr);

/// Decode one conserved two-temperature state. `phi_of_this_state` is the same
/// gravitational-potential gauge the release decode uses. Both temperature floors
/// are reported through `clamped_*` rather than silently applied.
TwoTempThermo two_temp_decode(const EosGammaTable& table, double rho,
                              double momentum, double energy, double e_elec,
                              double phi_of_this_state,
                              double temperature_guess,
                              bool* clamped_Te = nullptr,
                              bool* clamped_Ti = nullptr);

/// Build one on-manifold FACE state from the reconstructed primitive set
/// (rho, u, p, T_e). Fully algebraic — no temperature inversion anywhere:
///   x = x_Saha(n_H, T_e), p_e = x n_H k_B T_e, p_i = p - p_e,
///   T_i = p_i/(n_H k_B).
/// If the reconstruction leaves p_e >= p (which would imply T_i <= 0) the
/// ELECTRON temperature is reduced until p_e = `kFacePeFraction` * p, because the
/// total pressure is the authoritative mechanical variable and must be preserved
/// exactly. `clamped` reports that event.
TwoTempThermo two_temp_face_state(const EosGammaTable& table, double rho,
                                  double velocity, double pressure,
                                  double T_e, bool* clamped = nullptr);

/// Largest electron fraction of the total pressure a reconstructed face state is
/// allowed to carry [-]. One minus this is the fraction guaranteed to the heavy
/// species, so T_i stays strictly positive.
constexpr double kFacePeFraction = 0.999999;

/// Total energy of a two-temperature state,
/// E = (3/2) p_i + E_e + 1/2 rho u^2 + rho phi [J m^-3].
inline double two_temp_total_energy(double rho, double momentum,
                                    double e_heavy, double e_elec, double phi) {
    return e_heavy + e_elec + 0.5*momentum*momentum/rho + rho*phi;
}

/// Cell-centre gravitational potential [J kg^-1], the gauge every decode uses.
/// Deliberately the same expression as the release `mixture_cell_phi`, restated
/// here so `src/two_temp/` needs no header from `src/single_fluid/`.
inline double two_temp_cell_phi(const Grid& grid, arma::uword i) {
    return 0.5*(static_cast<double>(grid.phi_g_imh(i))
              + static_cast<double>(grid.phi_g_iph(i)));
}

// ============================================================================
// Collision and conduction closures
// ============================================================================

/// Spitzer ELECTRON field-aligned heat conductivity [W m^-1 K^-1], evaluated at
/// the ELECTRON temperature. Identical expression to the release
/// `physical_kappa_e`; the only change is which temperature it is asked for.
double two_temp_kappa_e(double n_e, double n_HI, double T_e);

/// Neutral-hydrogen field-aligned heat conductivity [W m^-1 K^-1], evaluated at
/// the HEAVY-particle temperature. Identical expression to the release
/// `physical_kappa_n`. Proton (Spitzer ion) conduction is ~2.3 % of kappa_e for
/// hydrogen and is deliberately NOT added, so that kappa_e + kappa_n reproduces
/// the release total conductivity exactly at T_e = T_i and the experiment changes
/// only WHICH temperature each existing channel acts on.
double two_temp_kappa_i(double n_e, double n_HI, double T_i);

/// Electron-heavy temperature-equilibration conductance [W m^-3 K^-1] such that
/// Q_ei = g_ei (T_i - T_e) heats the electrons. It is the translational electron
/// heat capacity (3/2) n_e k_B times the total equilibration frequency
/// nu_ei + nu_en, using the same Spitzer/NRL electron-ion rate and the same
/// Vranjes & Krstic (2013) electron-neutral cross-section as the two-fluid
/// solver's `physics.hpp::nu_ei` / `nu_en`. Both channels matter: the lower
/// chromosphere is weakly ionized, so n_HI >> n_e there and electron-neutral
/// collisions carry a large share of the electron thermalization.
double two_temp_exchange_conductance(double n_e, double n_HI, double T_e);

// ============================================================================
// Diagnostics
// ============================================================================

/// Running tally of every guard the two-temperature solver can trip, so a
/// clamped state is diagnosable rather than silent. Reported once at the end of
/// the run by the driver.
struct TwoTempStats {
    std::uint64_t face_pe_clamps = 0;    ///< Reconstructed faces whose p_e was capped below p.
    std::uint64_t decode_Te_clamps = 0;  ///< Decodes whose E_e fell outside the caloric bracket.
    std::uint64_t decode_Ti_clamps = 0;  ///< Decodes whose derived e_heavy hit the T_i floor.
    std::uint64_t godunov_faces = 0;     ///< Faces the exact Riemann solve was asked for.
    std::uint64_t godunov_fallbacks = 0; ///< Faces that fell back to face-local Rusanov.
    /// Largest |T_e - T_i| / T_i seen at any cell over the whole run [-]. This is
    /// the headline measurement of the study.
    double max_relative_decoupling = 0.0;
    /// Largest number of Newton passes the coupled conduction solve needed [-].
    std::uint64_t max_conduction_iterations = 0;

    /// Total number of clamp events of any kind.
    std::uint64_t clamps() const {
        return face_pe_clamps + decode_Te_clamps + decode_Ti_clamps;
    }
};

// ============================================================================
// Decoded field
// ============================================================================

/// Decode of one immutable two-temperature state plus its two ghost layers at
/// each end, in the same extended (ns+4) indexing the release uses.
struct TwoTempField {
    std::vector<TwoTempThermo> cells;  ///< physical cells, index 0..ns-1
    /// (ns+4) x 4 reconstruction primitives (log rho, u, log p, log T_e).
    arma::Col<double> extended_primitive;
    arma::uword ns = 0;  ///< cell count this field was decoded for

    /// Size the buffers for an `ns`-cell mesh.
    void resize(arma::uword ns_in);
};

/// Reconstruction slots of the two-temperature MUSCL step. Slot 3 is the
/// ELECTRON thermal variable; slot 2 remains the total pressure, so the
/// authoritative mechanical quantity is still reconstructed directly and the
/// electron/heavy split is the only thing the new slot changes.
enum TwoTempSlot : arma::uword {
    TT_LOG_RHO = 0, TT_U = 1, TT_LOG_P = 2, TT_LOG_TE = 3
};

/// One side of every interface, as the two-temperature flux needs it.
struct TwoTempFaceArrays {
    arma::Col<double> conserved;  ///< ns*4, (rho, rho u, E, E_e)
    arma::Col<double> flux;       ///< ns*4, the physical field-aligned flux
    arma::Col<double> velocity;   ///< ns [m s^-1]
    arma::Col<double> pressure;   ///< ns, TOTAL pressure [Pa]
    arma::Col<double> pressure_e; ///< ns, electron partial pressure [Pa]
    arma::Col<double> sound_speed;///< ns, frozen sqrt(5/3 p/rho) [m s^-1]

    /// Size every array above for an `ns`-cell mesh.
    void resize(arma::uword ns);
};

/// Caller-owned reusable workspace of the two-temperature solver. The driver owns
/// exactly one and passes it to every stage, which keeps this experimental solver
/// entirely out of the shared `Grid` and off the release path.
struct TwoTempWorkspace {
    TwoTempField decoded;    ///< decode of the state at the start of the step
    TwoTempField predicted;  ///< decode of the MUSCL-Hancock predictor state
    std::array<arma::Col<double>, 24> work;    ///< untyped reconstruction scratch
    std::array<TwoTempFaceArrays, 4> faces;    ///< the four one-sided face states
    arma::Col<double> flux_iph;   ///< numerical flux at the upper face i+1/2, ns*4
    arma::Col<double> flux_imh;   ///< numerical flux at the lower face i-1/2, ns*4
    arma::Col<double> source;     ///< momentum + electron-work source, ns*4
    arma::Col<double> u_face_iph; ///< numerical face velocity at i+1/2, ns [m s^-1]
    arma::Col<double> u_face_imh; ///< numerical face velocity at i-1/2, ns [m s^-1]
    Vec predicted_state;                       ///< packed predictor state
    Vec rhs;                                   ///< assembled explicit RHS

    /// Non-uniform MUSCL reconstruction weights and slope-ratio metrics. Functions
    /// of the STATIC mesh alone, so they are rebuilt only when
    /// `Grid::metrics_generation()` changes.
    ///@{
    arma::Col<double> W1, W3, W4;
    arma::Col<double> metric_r, metric_r_ip1, metric_r_im1;
    std::uint64_t weights_generation = 0;
    bool weights_valid = false;
    ///@}

    /// Per-cell scratch of the coupled implicit conduction solve.
    ///@{
    std::vector<double> c_rho, c_Te, c_Ti, c_ne, c_nHI, c_nH;
    std::vector<double> c_Ee_old, c_eh_old, c_Ee_at_T, c_cap_e, c_cap_i;
    std::vector<double> c_ke, c_ki, c_gE_L, c_gE_R, c_gI_L, c_gI_R, c_gex;
    std::vector<double> c_Ee_target, c_eh_target;
    /// 2x2-block tridiagonal Newton system: `blk_a`/`blk_c` are the diagonal
    /// sub/super blocks (electron and heavy channels never couple across cells),
    /// `blk_b` is the full 2x2 diagonal block carrying the collisional exchange.
    std::vector<double> blk_a, blk_c;                 ///< 2 per cell
    std::vector<double> blk_b;                        ///< 4 per cell, row-major
    std::vector<double> blk_rhs, blk_delta;           ///< 2 per cell
    ///@}

    /// Running tally of the guards; accumulated over the whole run.
    TwoTempStats stats;

    /// Size every buffer above for an `ns`-cell mesh.
    void resize(arma::uword ns);
};

// ============================================================================
// Solver stages
// ============================================================================

/// Decode `state` and its four ghost layers into `output`.
void two_temp_decode_into(const Grid& grid, const Vec& state,
                          TwoTempField& output, TwoTempStats& stats,
                          const TwoTempField* previous = nullptr);

/// Pack one two-temperature cell/ghost from primitive data (rho, u, T_e, T_i).
/// The packed state is on the Saha manifold at T_e by construction.
void two_temp_pack_ghost(const Grid& grid, Vec& ghost, double rho,
                         double velocity, double T_e, double T_i, double phi);
/// Same, into cell `cell` of a packed two-temperature state vector.
void two_temp_pack_cell(const Grid& grid, Vec& state, arma::uword cell,
                        double rho, double velocity, double T_e, double T_i,
                        double phi);

/// Explicit field-aligned RHS: MUSCL reconstruction of
/// (log rho, u, log p, log T_e) with the MC3/Koren limiter, a MUSCL-Hancock
/// predictor half step applying the same sources as the corrector, the release
/// SWMF exact-Riemann Godunov flux on the three release rows (face-local Rusanov
/// fallback), the electron energy upwinded on the contact with the exact mass
/// flux, and the geometric pressure/gravity plus electron-compression sources.
Vec two_temp_rhs_explicit(const Grid& grid, const Vec& state,
                          const TwoTempField& decoded, double dt_predictor,
                          TwoTempWorkspace& work);

/// CFL-limited timestep from |u| + sqrt(5/3 p/rho). Identical rule to the
/// release Godunov path: the frozen acoustic speed depends only on the TOTAL
/// pressure, so splitting the temperature does not change it. The stiff
/// collisional exchange imposes no step limit because it is solved implicitly.
Vec two_temp_timestep(const Grid& grid, const Vec& state,
                      const TwoTempField& decoded);

/// Coupled implicit backward-Euler conduction and collisional exchange:
/// Spitzer kappa_e on T_e, neutral-hydrogen kappa_n on T_i, and the pointwise
/// Q_ei = g_ei (T_i - T_e) exchange, solved together by Newton on
/// (T_e, T_i) with a 2x2-block tridiagonal Jacobian. Mass and momentum are
/// untouched, the exchange is antisymmetric, and both conductive fluxes are
/// re-summed into the total energy row, so the stage conserves total energy by
/// construction.
///
/// Boundary conditions, from the characteristic budget of the hyperbolic system
/// and the physics of the two conduction channels:
///   * outer face — DIRICHLET on T_e at `Grid::outer_conduction_temperature`
///     (the TR/coronal conductive reservoir is electron-conducted), NEUMANN
///     (zero flux) on T_i, because the fully ionized plasma above the domain has
///     no heavy-particle conductive channel to supply;
///   * inner face — Dirichlet on BOTH temperatures at the truncation reservoir
///     temperature, which is collisionally equilibrated there, unless
///     `Grid::inner_conduction_neumann` makes the base insulating for both.
Vec two_temp_apply_conduction(const Grid& grid, const Vec& state, double dt,
                              TwoTempWorkspace& work);

/// The complete experimental two-temperature timestep:
///   U^n -> MUSCL-Hancock/Godunov hydro -> U* -> implicit coupled
///          (T_e, T_i) conduction + exchange -> U^{n+1}
Vec two_temp_advance(Grid& grid, const Vec& state, const Vec& dt_i,
                     const TwoTempField& decoded, TwoTempWorkspace& work);

/** @} */

} // namespace chromosphere
