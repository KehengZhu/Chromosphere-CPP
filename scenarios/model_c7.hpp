/*!
 * @file scenarios/model_c7.hpp
 * @brief Model C7 quiet-Sun chromosphere scenario and the shared C7 atmosphere library.
 * @ingroup scenarios
 *
 * Model C7 chromosphere scenario (writeup Table 1).
 *
 * `model_c7_ic` configures the Grid's cell-centered fields (ds_i, B, boundary
 * buffers) and returns the initial conserved-state vector. `model_c7_update_bc`
 * applies the model_c7-specific boundary conditions:
 *
 *   * Outer face (lower TR, h ≈ 2153 km) — two modes. In BOTH of them the coronal
 *     heat enters as the imposed Stage-D Neumann flux q(T), never as a Dirichlet
 *     hot wall: on this coarse grid a fixed coronal ghost temperature would
 *     over-conduct (κ_e ∝ T^{5/2}) and over-pressurize the top cell.
 *       default (tr_jump_bc=false): fixed lower-TR reservoir (RTV 1978) — ρ,T
 *         Neumann (continuous pressure across the constant-pressure TR base),
 *         V,U Mach-capped outflow (U = V).
 *       "New explanation" (tr_jump_bc=true): hydrostatic ghost pressure (centered
 *         HSE ⇒ V = 0 is a boundary fixed point, no boundary downflow), ghost
 *         temperature CONTINUOUS from the top cell (T_g0 = a·T_top,
 *         T_g1 = b·T_g0, with a = b = 1 by default so there is no jump —
 *         C7_TJUMP_A / C7_TJUMP_B are an experimentation lever only), and an EOS
 *         ghost density at that (p, T); V,U Mach-capped outflow. The ghost
 *         temperature sets the EOS density and the hydro flux; it does NOT drive
 *         conduction.
 *
 *   * Inner face — photospheric discrete-HSE V=U=0 reservoir (Pandey 2024):
 *       ρ_i, ρ_n        → Dirichlet (pinned at the C7 base captured in IC)
 *       p_i, p_n        → hydrostatic ghost p = p_0 + ρ_0 g Δs
 *       V, U            → 0 (exact momentum fixed point)
 */

#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Build the Model C7 initial conserved-state vector and populate the Grid's
/// cell sizes, B-field, gravity potential, and boundary buffers.
///
/// `extended` selects the domain top: false (default) ⇒ the standard column to
/// h ≈ 2.153 Mm (T ≈ 23 kK, lower TR) with the RTV-estimated coronal flux q(T) —
/// the baseline model_c7 / model_flare use this. true ⇒ extend through the real
/// Avrett & Loeser (2008) Table 26 rows up to h ≈ 2.628 Mm (T ≈ 0.29 MK, upper
/// TR), with q(T) computed self-consistently from C7's own conductive gradient
/// at the new top — the corona-as-boundary configuration for gentle
/// conduction-driven evaporation (docs/studies/evaporation/gentle_evaporation_plan.md v2; model_gentle).
///
/// `tr_jump_bc` selects the upper boundary closure. false (default) ⇒ the
/// Mach-capped RTV reservoir (model_flare / analytic_canopy / model_gentle keep
/// this). true ⇒ the "New explanation" upper BC
/// (docs/studies/evaporation/gentle_evaporation_downflow.md): well-balanced reconstruction +
/// hydrostatic ghost pressure + a ghost temperature extrapolated continuously
/// from the top cell + EOS ghost density — set by the bare `model_c7` scenario.
/// Both modes keep the imposed Neumann flux q(T) as the conductive driver;
/// `tr_jump_bc` does NOT switch it off. The ghost-temperature factors are
/// env-tunable but default to C7_TJUMP_A = C7_TJUMP_B = 1, i.e. no jump.
Vec model_c7_ic(Grid& grid, bool extended = false, bool tr_jump_bc = false);

/// Sample the complete Model C7 atmosphere (Avrett & Loeser 2008 Table 26,
/// photosphere h ≈ −10 km → corona h ≈ 68 Mm) at height `h_km`: temperature `T`
/// linearly, electron density `n_e` and neutral-hydrogen density `n_HI`
/// log-linearly. Concatenates the photosphere rows (h < 1003 km, absent from the
/// model_c7 chromosphere/corona table) with MODEL_C7; clamps to the table ends.
/// Used by model_column's corona testbed to set a realistic frozen ionization
/// profile (so the corona is ionized and conducts via Spitzer κ_e).
void c7_full_profile(float h_km, float& T, float& n_e, float& n_HI);

/// Shape-preserving C1 interpolation of the complete Model C7 temperature
/// table. Unlike c7_full_profile's historical piecewise-linear temperature,
/// this has a continuous dT/dh and therefore does not seed knot-scale jumps in
/// conductive heat flux. Used only by gamma-table model_column initialization;
/// fixed-gamma scenarios retain their existing interpolation exactly.
double c7_full_temperature_pchip(double h_km);

/// Route-B photoionization closure (docs/studies/eos-ionization/photoionization_c7_inversion_plan.md,
/// writeup §3.1): given cell-centered (T, n_i, n_n) and heights [km], return the
/// per-cell photoionization rate P_phot that makes the profile a fixed point of
/// the active Stage-E ionization network, blended above ~1500 km to the frozen
/// FAL-C rate (Chae 2021). Shared by model_c7_ic and model_column_ic.
Vec c7_route_b_photoionization(const Grid& grid, const Vec& T_c,
                               const Vec& ni_c, const Vec& nn_c,
                               const Vec& h_cell);

/// Refresh ghost cells. Applies the coronal q(T) time ramp when one is
/// configured, lays down the default ghost layer with `apply_open_bcs`, then
/// OVERWRITES both faces: the outer face with whichever of the two lower-TR
/// closures `tr_jump_bc` selected, and the inner face with the photospheric
/// discrete-HSE reservoir (Dirichlet ρ_i and ρ_n pinned at the C7 base captured
/// in the IC, hydrostatic ghost pressures, V = U = 0).
/// See the file-level comment above for the full BC table.
void model_c7_update_bc(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
