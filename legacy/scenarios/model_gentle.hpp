/*!
 * @file legacy/scenarios/model_gentle.hpp
 * @brief LEGACY: `model_gentle`, the historical two-fluid field-aligned
 *        chromosphere-to-corona column.
 * @ingroup scenarios
 *
 * The historical two-fluid research preset of the C7 column: a fixed adiabatic
 * index, a neutral fluid, finite-rate ionization, a radiative sink, TRAC and
 * the mesh-scaled artificial conduction term, on the seven-row two-fluid solver
 * (legacy/two_fluid/). It never loads a Gamma1 table. It used to share one
 * IC/BC implementation with the release `model_column`; it now lives here,
 * frozen, and is built only with -DCHROMO_ENABLE_LEGACY=ON.
 *
 * Its environment knobs are the column knobs listed in
 * scenarios/model_column.hpp plus the two-fluid-only ones the release rejects:
 * ISO_GAMMA, ISO_TWO_FLUID, ISO_IONIZATION, ISO_COOLING, ISO_CORONA /
 * ISO_CORONA_TOP_KM, ISO_TRAC*, ISO_QFLUX*, ISO_CHEAT*, ISO_TBOOST*,
 * ISO_NUMERICAL_DIFFUSIVITY_MULT.
 */

#pragma once

#include "core/chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Initial condition of the legacy two-fluid `model_gentle` column.
///
/// Throws if the Grid carries a Gamma1 table. Builds the (optionally statically
/// refined) mesh, samples the Model C7 temperature and ionization fraction,
/// re-integrates the density hydrostatically at the fixed adiabatic index,
/// packs the seven `cons::` rows, captures the first-order inner ghost ladder,
/// configures the two-fluid research stages from the environment, fills both
/// ghost layers via model_gentle_update_bc(), and freezes the result as
/// `Grid::eq_state` with `Grid::eq_wb` set.
Vec  model_gentle_ic(Grid& grid);

/// Per-step boundary refresh of `model_gentle`: the imposed q(T) flux ramp, the
/// relaxation → Stage-2 switch, the fixed-back-pressure hydrostatic outer ghost
/// ladder with Mach-capped outflow, and the first-order inner reservoir ghosts
/// (V = U = 0), all in the two-fluid `cons::` convention.
void model_gentle_update_bc(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
