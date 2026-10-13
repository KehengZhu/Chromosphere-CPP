/*!
 * @file scenarios/model_column.hpp
 * @brief model_column (RELEASE) and model_column_2t (EXPERIMENTAL) — the
 *        field-aligned chromosphere-to-corona column.
 * @ingroup scenarios
 *
 * model_column — the field-aligned chromosphere→corona column. One IC/BC
 * implementation for the two Saha/Gamma1 scenarios:
 *
 *   * `model_column`  — the DEFAULT and the RELEASE scenario. It always loads
 *     the production Gamma1 table, so it always runs the single-fluid
 *     equilibrium-mixture solver (src/single_fluid/) with the state
 *     U = (rho, rho u, E) and the timestep
 *         U^n -> MUSCL-Hancock/Godunov hydro -> U* -> implicit physical conduction -> U^{n+1}.
 *     Every historical two-fluid knob listed below is REJECTED outright, so
 *     `model_column` names exactly one production model.
 *
 *   * `model_column_2t` — the EXPERIMENTAL two-temperature twin (src/two_temp/).
 *
 * The HISTORICAL two-fluid research preset `model_gentle`, which used to share
 * this IC/BC, is a separate legacy scenario (legacy/scenarios/model_gentle.hpp,
 * built only with -DCHROMO_ENABLE_LEGACY=ON).
 *
 * A single straight field line (B = 1, gravity on) initialized from the real
 * Model C7 atmosphere (Avrett & Loeser 2008, Table 26; via c7_full_profile and
 * c7_full_temperature_pchip in c7_profile.cpp), driven to conduction-driven
 * ("gentle", Antiochos & Sturrock 1978) evaporation.
 *
 * The numerical scheme is fixed here, no toggles (full statement in
 * docs/standards/release_facts.md):
 *   - MUSCL-Hancock reconstruction of (ln rho, u, ln p) with the asymmetric
 *     MC3/Koren limiter (β = 2); the ln p slot is hydrostatic: the limiter acts
 *     only on the departure from discrete hydrostatic balance, and the gravity
 *     source uses the same discrete relation. The face temperature is recovered
 *     from the reconstructed (rho, p) through the same Saha closure,
 *   - the SWMF exact-Riemann Godunov flux at frozen composition (gamma = 5/3),
 *   - the momentum source (gravity) applied in both the Hancock predictor and
 *     the corrector,
 *   - an EOS-closed hydrostatic ladder, in that same discrete relation, for the
 *     two ghosts below the 1600 km mid-chromospheric truncation.
 * The scheme is reference-free (no frozen equilibrium is stored or subtracted)
 * and exactly balanced at rest: a column in the discrete hydrostatic relation is
 * a fixed point to round-off (N=500, t=0: max|F^rho| = 1.5e-19 kg m^-2 s^-1;
 * the plain ln p scheme it replaced left 6.0e-14, second order). Exact balance
 * holds for B = 1, the current model. Evidence:
 * docs/studies/numerics/hydrostatic_reconstruction_and_mc3_fix.md.
 *
 * IC: the C7 temperature and ionization FRACTION with the total density re-integrated
 * hydrostatically from the base (clean V ≈ 0 start).
 *
 * Env knobs (numerics are fixed; only physics is tunable):
 *   ISO_NS            grid cells (release preset: 500 coarse-equivalent)
 *   ISO_H_BASE        domain base height [km]   (release preset: 1600)
 *   ISO_DH            domain height span [km]   (release preset: 553)
 *   ISO_HEAT_FLUX     0 = Stage-1 relaxation (conduction off) / 1 = conduction on
 *   ISO_RELAX_TIME    Stage-1 adiabatic relaxation phase length [s] (default 0)
 *   ISO_TJUMP_A/_B    top ghost-T jump factors a, b (default 1)
 *   ISO_T_TOP         absolute imposed top-boundary temperature [K] (release: 22000)
 *   ISO_HYDRO_T_DECOUPLE  1 = hydro ghost T free-floats and conduction keeps the
 *                     fixed wall at the PHYSICAL FACE; 0 = one fixed ghost T
 *                     serves hydro and conduction (release preset: 0)
 *   ISO_VCAP          outer outflow Mach cap (default 0.1)
 *   ISO_INNER_T_NEUMANN  conduction lower BC: 0 = Dirichlet (default) / 1 = Neumann
 *   ISO_REFINE_*      static local mesh refinement (scenarios/mesh.hpp)
 *
 * HISTORICAL two-fluid-only knobs — rejected here; they configure the legacy
 * model_gentle scenario:
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

#include "core/chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Initial condition of `model_column` (release), also the base of
/// `model_column_2t`.
///
/// Loads the `GAMMA_TABLE` when one is given and throws if no Gamma1 table is
/// loaded (the fixed-gamma two-fluid column is the legacy `model_gentle`
/// scenario). It REJECTS every historical two-fluid knob listed in the file
/// header, so a release run cannot be silently reconfigured. Then it
/// builds the mesh (optionally statically refined, which may resize the Grid),
/// sets the field-line geometry and gravity, samples the Model C7 temperature
/// and ionization fraction, re-integrates the density hydrostatically from the
/// base for a clean V ≈ 0 start, fills both ghost layers via
/// model_column_update_bc(), and returns the initial conserved state — three
/// `mix::` rows (four `tt::` rows on the two-temperature leg).
///
/// The release is REFERENCE-FREE: it leaves `Grid::eq_wb` false and stores
/// no frozen equilibrium state.
Vec  model_column_ic(Grid& grid);

/// Per-step boundary refresh of `model_column` (and, through
/// model_column_2t_update_bc(), of `model_column_2t`).
///
/// Decodes the top and base interior cells with whichever state convention is
/// active (the mixture decode for the release, the two-temperature decode on
/// the experimental leg) and rewrites the two outer and two inner ghost layers:
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
void model_column_update_bc(Grid& grid, const Vec& xn);

/// EXPERIMENTAL. Initial condition of `model_column_2t`, the two-temperature
/// (T_e != T_i) test-study twin of `model_column`.
///
/// It is the SAME atmosphere, mesh, gravity, numerics and boundary geometry as
/// the release column — that is the point of the study, so nothing about the
/// configuration is allowed to differ — with exactly three changes:
///
///   * the conserved state is the four `tt::` rows (rho, rho u, E, E_e) of
///     `src/two_temp/`, and `Grid::two_temperature` is set so the driver selects
///     that solver;
///   * every cell and ghost is initialized with T_e = T_i equal to the release
///     temperature. That is not a simplification: the measured electron-heavy
///     equilibration time in this column is 0.09-3.9 ms
///     (`util/release/tau_ei_equilibration.py`), so a common initial temperature is the
///     physically correct starting state and any later split is generated by the
///     solver rather than imposed by the IC;
///   * `TT_OUTER_TI` selects the outer-face heavy-particle conduction condition,
///     `neumann` (default, zero heavy conductive flux from the fully ionized
///     plasma above the domain) or `dirichlet` (the controlled comparison, the
///     electron wall temperature imposed on T_i as well).
///
/// Both outer-face hydrodynamic temperatures free-float with the top cell in
/// either case, because a subsonic outflow face admits exactly ONE externally
/// imposed condition and the reservoir back-pressure already spends it.
Vec  model_column_2t_ic(Grid& grid);

/// EXPERIMENTAL. Per-step boundary refresh of `model_column_2t`: the same
/// closure as `model_column_update_bc` — fixed reservoir back-pressure at the
/// top, quasi-static stratified reservoir at the 1600 km base — rewritten for the
/// four-row two-temperature ghosts, with BOTH temperatures zero-gradient
/// extrapolated at the outer face and the ghost density inverted from the
/// two-temperature pressure closure at that pair of temperatures.
void model_column_2t_update_bc(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
