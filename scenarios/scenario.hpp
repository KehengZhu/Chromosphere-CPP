/*!
 * @file scenarios/scenario.hpp
 * @brief Scenario dispatcher.
 * @ingroup scenarios
 *
 * Scenario dispatcher.
 *
 * A Scenario packages the three things chromo_main needs from any
 * initial-condition + boundary-condition source: how big the grid should be
 * (peek_ns), how to populate it (ic), and how to refresh the boundary each
 * step (update_bc). make_scenario() returns the requested Scenario by name.
 */

#pragma once

#include <functional>
#include <string>

#include "core/chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// One complete initial-condition + boundary-condition source, in the form
/// chromo_main consumes it.
///
/// A Scenario packages the three things the driver needs from ANY IC/BC source,
/// tabulated or analytic, plus a name for logging and sidecar metadata:
///
///   * peek_ns    — how many cells the grid should have. Queried BEFORE
///                  Grid::init(), because the Grid must be allocated first.
///   * ic         — fill the freshly allocated Grid (geometry, gravity, ghost
///                  buffers) and return the initial conserved state vector.
///   * update_bc  — refresh the Grid's boundary buffers from the current state,
///                  once per timestep.
///
/// Holding them in std::function keeps the driver independent of which scenario
/// is running; make_scenario() binds the free functions of one scenario
/// translation unit into this struct.
struct Scenario {
    /// Returns the grid resolution this scenario expects (a tabulated
    /// scenario reads it from its data file; an analytic one returns a
    /// sensible default). Called before Grid::init().
    std::function<arma::uword()> peek_ns;

    /// Populate the Grid's cell-centered + boundary fields and return the
    /// initial conserved-state vector.
    std::function<Vec(Grid&)> ic;

    /// Refresh the boundary buffers from the current state. Called once
    /// per timestep before advance_Euler_state.
    std::function<void(Grid&, const Vec&)> update_bc;

    /// Scenario identifier as spelled on the command line (e.g. "model_column").
    /// Used for log lines and run-metadata sidecars.
    std::string name;
};

/// Build a Scenario by name: `model_column` (release) or `model_column_2t`
/// (experimental), or — only in a -DCHROMO_ENABLE_LEGACY=ON build — one of the
/// legacy two-fluid scenarios `model_c7`, `model_flare`, `model_gentle`,
/// `analytic_canopy`, `pfss_field_line`. `data_path` is required when the
/// scenario reads a tabulated input file (e.g. "pfss_field_line"); ignored
/// otherwise. Throws std::runtime_error on an unknown name, a missing required
/// data_path, or a legacy name in a build without the legacy code (the message
/// names -DCHROMO_ENABLE_LEGACY=ON).
Scenario make_scenario(const std::string& name, const std::string& data_path = "");

/** @} */

} // namespace chromosphere
