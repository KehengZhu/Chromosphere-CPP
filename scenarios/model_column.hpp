/*!
 * @file scenarios/model_column.hpp
 * @brief model_column (RELEASE) — the field-aligned chromosphere-to-corona
 *        column.
 * @ingroup scenarios
 *
 * model_column — the field-aligned chromosphere→corona column, the DEFAULT and
 * the RELEASE scenario. It always loads the production Gamma1 table, so it
 * always runs the single-fluid equilibrium-mixture solver (src/single_fluid/)
 * with the state U = (rho, rho u, E) and the timestep
 *     U^n -> MUSCL-Hancock/Godunov hydro -> U* -> implicit physical conduction -> U^{n+1}.
 * Every historical two-fluid knob listed below is REJECTED outright, so
 * `model_column` names exactly one production model.
 *
 * Two scenarios that used to share this IC/BC are legacy now (built only with
 * -DCHROMO_ENABLE_LEGACY=ON): the HISTORICAL two-fluid research preset
 * `model_gentle` (legacy/scenarios/model_gentle.hpp) and the two-temperature
 * test-study twin `model_column_2t` (legacy/scenarios/model_column_2t.hpp).
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

/// Initial condition of `model_column` (release).
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
/// `mix::` rows per cell.
///
/// The release is REFERENCE-FREE: it leaves `Grid::eq_wb` false and stores
/// no frozen equilibrium state.
Vec  model_column_ic(Grid& grid);

/// Per-step boundary refresh of `model_column`.
///
/// Decodes the top and base interior cells and rewrites the two outer and two
/// inner ghost layers:
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

/** @} */

} // namespace chromosphere
