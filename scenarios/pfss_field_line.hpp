/*!
 * PFSS field-line scenario.
 *
 * Initializes the grid from a tabulated data file produced by
 * util/extract_field_line.py (PFSS extrapolation of a synoptic magnetogram,
 * traced from a chosen footpoint to the top of the chromospheric domain).
 *
 * File format: see scenarios/data_file_parser.hpp.
 *
 * Boundary update follows the model_c7 convention: outer ghost densities and
 * temperatures stay pinned to the values written in the file's [GHOSTS]
 * section; outer velocities are halved each call to damp outflow.
 */

#pragma once

#include <string>

#include "chromosphere.hpp"

namespace chromosphere {

/// Read ns from the file's [META] section. Called before Grid::init.
arma::uword pfss_peek_ns(const std::string& data_path);

/// Populate the Grid from the file and return the initial conserved state.
Vec pfss_ic(Grid& grid, const std::string& data_path);

/// Refresh outer ghost cells (halve velocity, keep n/T pinned).
void pfss_update_bc(Grid& grid, const Vec& xn);

} // namespace chromosphere
