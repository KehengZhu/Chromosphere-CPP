/*!
 * @file scenarios/model_flare.hpp
 * @brief Explosive chromospheric evaporation scenario (nonthermal electron beam).
 * @ingroup scenarios
 *
 * Explosive chromospheric evaporation scenario (Fisher, Canfield & McClymont
 * 1985, ApJ 289, 425).
 *
 * Reuses the Model C7 chromosphere column (`model_c7_ic`) and its boundary
 * closures (`model_c7_update_bc`), then adds a transient nonthermal-electron-
 * beam volumetric heating source in the upper chromosphere and opens the outer
 * boundary so the evaporated plasma can flow out supersonically. With a beam
 * flux above the Fisher threshold (F_crit ≈ 7×10⁶ W m⁻²) the upper chromosphere
 * cannot radiate the deposited energy, is heated to coronal temperatures, and
 * drives an explosive upflow into the corona with a condensation downflow below.
 *
 * See docs/explosive_evaporation_plan.md.
 */

#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Build the flare IC: identical to model_c7_ic, then enable beam heating, force
/// ionization + radiative cooling on (the radiative sink defines the explosive
/// threshold), and open the outer boundary (free supersonic outflow). The beam
/// flux is read from the FLARE_BEAM_FLUX environment variable when set
/// (W/m²; default 5×10⁷ = explosive), so the explosive and gentle-control runs
/// share one build.
Vec model_flare_ic(Grid& grid);

/// Boundary refresh — delegates to model_c7_update_bc (which honors
/// grid.outer_free_outflow set by model_flare_ic).
void model_flare_update_bc(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
