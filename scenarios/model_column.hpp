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
 * Environment (release names; the full table is the configuration reference):
 *   CHROMO_MESH_NS              coarse-equivalent cell count (preset 500)
 *   CHROMO_MESH_REFINE_*        static local refinement (scenarios/mesh.hpp;
 *                               preset outer, factor 4, start 500 km, grade 20 km)
 *   CHROMO_DOMAIN_BASE_KM       lower truncation height [km] (preset 1600)
 *   CHROMO_DOMAIN_SPAN_KM       domain span [km] (preset 553)
 *   CHROMO_CONDUCTION           1 = conduction on (preset) / 0 = off
 *   CHROMO_BC_T_TOP_K           fixed outer ghost temperature [K] (preset 22000;
 *                               0 = the IC top-cell temperature)
 *   CHROMO_BC_VCAP_MACH         outer outflow Mach cap (default 0.1)
 *   CHROMO_IC_P_TOP_PA          rebuild the IC so the top cell sits at this p [Pa]
 *   CHROMO_EOS_TABLE            Gamma1 table path (preset: the production table)
 * Diagnostic / reference only, never production:
 *   CHROMO_BC_OUTER_P           reservoir (release) | hydrostatic
 *   CHROMO_COND_FACE            harmonic (release) | kirchhoff | arithmetic
 *   CHROMO_NUM_RIEMANN          godunov (release) | roe (the previous release
 *                               flux, the equilibrium-Gamma1 Roe linearization)
 *
 * Every ISO_* or GRID_REFINE_* variable, TT_OUTER_TI and GAMMA_TABLE stop a
 * release run:
 * the renamed ones name their CHROMO_* replacement, and the historical
 * two-fluid knobs (ISO_GAMMA, ISO_TWO_FLUID, ISO_COOLING, ISO_TRAC, ...) point
 * to the legacy model_gentle scenario.
 */

#pragma once

#include "core/chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Initial condition of `model_column` (release).
///
/// Loads the `CHROMO_EOS_TABLE` when one is given and throws if no Gamma1 table is
/// loaded (the fixed-gamma two-fluid column is the legacy `model_gentle`
/// scenario). It REJECTS every retired ISO_* name listed in the file header,
/// so a release run cannot be silently reconfigured. Then it
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
///     extrapolated but Mach-capped at `CHROMO_BC_VCAP_MACH` × c_s, and the
///     fixed 22,000 K ghost temperature, which conduction reads as its Dirichlet
///     datum at the ghost centre;
///   * inner — the quasi-static stratified reservoir at the 1600 km truncation
///     temperature `kInnerTRef` with V = 0, laid out on an EOS-closed ghost
///     ladder in the discrete hydrostatic relation of the reconstruction.
void model_column_update_bc(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
