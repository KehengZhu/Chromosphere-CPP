/*!
 * PFSS field-line scenario.
 *
 * Initializes the grid from a tabulated data file produced by
 * util/extract_field_line.py (PFSS extrapolation of a synoptic magnetogram,
 * traced from a chosen footpoint to the top of the chromospheric domain).
 *
 * File format: see scenarios/data_file_parser.hpp.
 *
 * Boundary update delegates to `apply_open_bcs` (scenario.hpp): reflecting
 * wall at the inner face, zero-gradient outflow at the outer face. The
 * [GHOSTS] section of the data file is therefore only used to seed the
 * pre-step ghost state in `pfss_ic`; the first `update_bc` call replaces it.
 */

#pragma once

#include <string>

#include "chromosphere.hpp"

namespace chromosphere {

/// Read ns from the file's [META] section. Called before Grid::init.
arma::uword pfss_peek_ns(const std::string& data_path);

/// Populate the Grid from the file and return the initial conserved state.
///
/// Opt-in flare mode (env PFSS_FLARE=1): overlay the model_flare beam-heating
/// physics (Fisher et al. 1985 explosive evaporation) on this field line —
/// forces ionization + radiative cooling, adds numerical diffusivity / TRAC /
/// RTV q(T), and switches on the nonthermal beam. Beam knobs share the FLARE_*
/// env names with model_flare_ic. Default (unset) leaves the scenario unchanged.
Vec pfss_ic(Grid& grid, const std::string& data_path);

/// Refresh ghost cells via `apply_open_bcs`: inner reflecting wall, outer
/// zero-gradient outflow.
void pfss_update_bc(Grid& grid, const Vec& xn);

} // namespace chromosphere
