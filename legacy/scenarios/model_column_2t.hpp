/*!
 * @file legacy/scenarios/model_column_2t.hpp
 * @brief LEGACY: the two-temperature (T_e != T_i) twin of the release column.
 * @ingroup scenarios
 *
 * Built only with -DCHROMO_ENABLE_LEGACY=ON. The solver lives in
 * legacy/two_temp/; the study it served is
 * docs/studies/conduction/two_temperature_test_study.md. The IC/BC is a frozen
 * copy of the column implementation as it stood at commit 7e06466, so the study's
 * runs reproduce bit for bit with the historical `ISO_*` environment names.
 */
#pragma once

#include "core/chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// LEGACY. Initial condition of `model_column_2t`: the same atmosphere, mesh,
/// gravity, numerics and boundary geometry as the release column of commit
/// 7e06466, with the conserved state in the four `tt::` rows
/// (rho, rho u, E, E_e), T_e = T_i at t = 0, and `Grid::two_temperature` set
/// so the driver selects the legacy/two_temp/ solver. `TT_OUTER_TI` selects
/// the outer heavy-particle conduction condition, `neumann` (default) or
/// `dirichlet`. Both outer-face hydrodynamic temperatures free-float with the
/// top cell (preset `ISO_HYDRO_T_DECOUPLE=1`).
Vec  model_column_2t_ic(Grid& grid);

/// LEGACY. Per-step boundary refresh of `model_column_2t`: fixed reservoir
/// back-pressure at the top, quasi-static stratified reservoir at the 1600 km
/// base, written into the four-row two-temperature ghosts.
void model_column_2t_update_bc(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
