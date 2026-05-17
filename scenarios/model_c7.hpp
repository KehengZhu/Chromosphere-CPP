/*!
 * Model C7 chromosphere scenario (writeup Table 1).
 *
 * `model_c7_ic` configures the Grid's cell-centered fields (ds_i, B, boundary
 * buffers) and returns the initial conserved-state vector. `model_c7_update_bc`
 * refreshes the outer ghost cells each step, halving the velocity to damp
 * outflow.
 */

#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

/// Build the Model C7 initial conserved-state vector and populate the Grid's
/// cell sizes, B-field, gravity potential, and boundary buffers.
Vec model_c7_ic(Grid& grid);

/// Recompute outer ghost cells: densities and temperatures are pinned to the
/// initial Model C7 values; velocity is halved each call.
void model_c7_update_bc(Grid& grid, const Vec& xn);

} // namespace chromosphere
