/*!
 * @file scenarios/model_column.hpp
 * @brief model_column (RELEASE) and model_gentle (historical two-fluid) — the
 *        field-aligned chromosphere-to-corona column.
 * @ingroup scenarios
 *
 * model_column — the field-aligned chromosphere→corona column. One IC/BC
 * implementation, TWO scenarios that never cross solvers:
 *
 *   * `model_column`  — the DEFAULT and the RELEASE scenario. It always loads
 *     the production Gamma1 table, so it always runs the single-fluid
 *     equilibrium-mixture solver (src/single_fluid/) with the state
 *     U = (rho, rho u, E) and the timestep
 *         U^n -> MUSCL-Hancock/Godunov hydro -> U* -> implicit physical conduction -> U^{n+1}.
 *     Every historical two-fluid knob listed below is REJECTED outright in this
 *     mode, so `model_column` names exactly one production model.
 *
 *   * `model_gentle`  — the HISTORICAL two-fluid research preset. It never loads
 *     a Gamma1 table, so it always runs the seven-row two-fluid solver
 *     (src/two_fluid/) with a fixed adiabatic index, finite-rate ionization, a
 *     radiative sink, TRAC and the mesh-scaled artificial conduction term.
 *
 * A single straight field line (B = 1, gravity on) initialized from the real
 * Model C7 atmosphere (Avrett & Loeser 2008, Table 26; via c7_full_profile in
 * model_c7.cpp), driven to conduction-driven ("gentle", Antiochos & Sturrock
 * 1978) evaporation.
 *
 * The numerical scheme is the well-balanced, high-order configuration validated in
 * docs/gentle_evaporation_downflow.md (the "bestwb" set) — ALWAYS ON here, no toggles:
 *   - well-balanced explicit reconstruction (φ_g-consistent shifted pressures),
 *   - log-space MUSCL reconstruction of the exponentially-stratified slots,
 *   - the asymmetric MC3/Koren limiter (β = 2),
 *   - equilibrium-reference ("δ-form") well-balancing (V = 0 exact for any strat.),
 *   - inner (photospheric) discrete-HSE reservoir well-balanced through BOTH ghosts,
 *     density included (ρ ∝ p at T₀) so V = 0 is an exact fixed point for mass too.
 * The release additionally uses the SWMF exact-Riemann Godunov flux (frozen
 * composition, gamma = 5/3) and the
 * (ln rho, u, ln p) primitive set. Together these hold a hydrostatic column at
 * V ≈ 0 to round-off, so any residual flow is physical (conduction-driven
 * evaporation), not truncation drainage.
 *
 * IC: the C7 temperature and ionization FRACTION with the total density re-integrated
 * hydrostatically from the base (clean V ≈ 0 start).
 *
 * Env knobs valid in BOTH modes (numerics are fixed; only physics is tunable):
 *   ISO_NS            grid cells (release preset: 500 coarse-equivalent)
 *   ISO_H_BASE        domain base height [km]   (release preset: 1600)
 *   ISO_DH            domain height span [km]   (release preset: 553)
 *   ISO_HEAT_FLUX     0 = Stage-1 relaxation (conduction off) / 1 = conduction on
 *   ISO_RELAX_TIME    Stage-1 adiabatic relaxation phase length [s] (default 0)
 *   ISO_TJUMP_A/_B    top ghost-T jump factors a, b (default 1)
 *   ISO_T_TOP         absolute imposed top-boundary temperature [K] (release: 22000)
 *   ISO_HYDRO_T_DECOUPLE  hydro ghost T free-floats; conduction keeps the fixed
 *                     PHYSICAL-FACE wall (release preset: 1)
 *   ISO_VCAP          outer outflow Mach cap (default 0.1)
 *   ISO_INNER_T_NEUMANN  conduction lower BC: 0 = Dirichlet (default) / 1 = Neumann
 *   ISO_REFINE_*      static local mesh refinement (scenarios/mesh.hpp)
 *
 * HISTORICAL two-fluid-only knobs — rejected by the release scenario:
 *   ISO_GAMMA         fixed adiabatic index of the two-fluid solver
 *   ISO_TWO_FLUID     evolve ion & neutral as separate fluids
 *   ISO_IONIZATION    Stage-E ionization/recombination network
 *   ISO_COOLING       optically-thick + thin radiative sink (Stage R)
 *   ISO_CORONA / ISO_CORONA_TOP_KM   extend into a resolved ~1 MK corona
 *   ISO_TRAC / _TCHROM / _TCMAXFRAC  TRAC TR-broadening
 *   ISO_QFLUX / _ENHANCE / _T_ON / _RAMP / _0   imposed coronal q(T) flux + ramp
 *   ISO_CHEAT / _E0 / _TMAX_MK / _SH_FRAC / _S0_MM / _ENHANCE / _T_ON / _RAMP
 *                     ambient volumetric coronal heating H(s) + ramp
 *   ISO_TBOOST / _TC / _W            one-time IC coronal superheat
 *   ISO_NUMERICAL_DIFFUSIVITY_MULT   mesh-scaled artificial conduction
 *
 * Reference-only numerical overrides (regression/comparison, never production):
 *   ISO_RIEMANN=roe-local (equilibrium-Gamma1 Roe linearization, the previous
 *   release flux), ISO_RIEMANN=rusanov, ISO_RECONSTRUCTION=lnrho-v-lnt,
 *   ISO_LIMITER
 */

#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Initial condition shared by `model_column` (release) and `model_gentle`
/// (historical two-fluid).
///
/// Loads the `GAMMA_TABLE` when one is given — a loaded Gamma1 table is what
/// selects the single-fluid release solver, an empty one leaves the two-fluid
/// solver — and in release mode REJECTS every historical two-fluid knob listed
/// in the file header, so a release run cannot be silently reconfigured. Then it
/// builds the mesh (optionally statically refined, which may resize the Grid),
/// sets the field-line geometry and gravity, samples the Model C7 temperature
/// and ionization fraction, re-integrates the density hydrostatically from the
/// base for a clean V ≈ 0 start, fills both ghost layers via
/// model_column_update_bc(), and returns the initial conserved state — three
/// `mix::` rows in release mode, seven `cons::` rows in two-fluid mode.
///
/// The release leg is REFERENCE-FREE: it leaves `Grid::eq_wb` false and stores
/// no frozen equilibrium state. Only the two-fluid leg sets eq_wb and freezes
/// this IC as `Grid::eq_state`.
Vec  model_column_ic(Grid& grid);

/// Per-step boundary refresh for both legs of the scenario.
///
/// Decodes the top and base interior cells with whichever state convention is
/// active (mixture decode for the release, cons2prim algebra for the two-fluid
/// path) and rewrites the two outer and two inner ghost layers:
///
///   * outer — subsonic outflow: a hydrostatic ghost ladder anchored on the fixed
///     reservoir back-pressure `kOuterPRef` (captured once, the single externally
///     imposed condition of a subsonic outflow face), the top-cell velocity
///     extrapolated but Mach-capped at `ISO_VCAP` × c_s, and — in the release
///     preset, where the hydro ghost temperature free-floats with the top cell —
///     the 22,000 K thermal datum handed to the conduction stage at the PHYSICAL
///     outer face instead;
///   * inner — the quasi-static stratified reservoir at the 1600 km truncation
///     temperature `kInnerTRef` with V = 0, laid out on a second-order
///     (trapezoidal, EOS-closed) ghost ladder.
///
/// It also advances the time-dependent two-fluid drivers (the imposed q(T) flux
/// ramp and the coronal-heating ramp), which are inactive in release mode.
void model_column_update_bc(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
