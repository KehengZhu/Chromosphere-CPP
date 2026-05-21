/*!
 * Model C7 chromosphere scenario (writeup Table 1).
 *
 * `model_c7_ic` configures the Grid's cell-centered fields (ds_i, B, boundary
 * buffers) and returns the initial conserved-state vector. `model_c7_update_bc`
 * applies the model_c7-specific boundary conditions:
 *
 *   * Outer face — Neumann ρ (via `apply_open_bcs`) with two cascades
 *     overlaid afterwards:
 *       Velocity damping: V_g0 = V/2,  V_g1 = V/4   (same for U)
 *       Temperature      : T_g0 = T,    T_g1 = 2·T  (sharp jump at face ns+½
 *                          to mimic the chromosphere→corona transition region)
 *
 *   * Inner face — photospheric reservoir BC:
 *       ρ_n, T_n        → Dirichlet (pinned at the C7 base values captured
 *                         during `model_c7_ic`)
 *       n_i, T_i, V, U  → Neumann (zero-gradient copy of cell 0)
 */

#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

/// Build the Model C7 initial conserved-state vector and populate the Grid's
/// cell sizes, B-field, gravity potential, and boundary buffers.
Vec model_c7_ic(Grid& grid);

/// Refresh ghost cells: `apply_open_bcs` for the outer face followed by the
/// model_c7 velocity (V/2, V/4) and temperature (T, 2T) cascades, plus the
/// photospheric inner BC (Dirichlet ρ_n, T_n; Neumann n_i, T_i, V, U).
/// See the file-level comment above for the full BC table.
void model_c7_update_bc(Grid& grid, const Vec& xn);

} // namespace chromosphere
