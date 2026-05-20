/*!
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

#include "chromosphere.hpp"

namespace chromosphere {

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

    std::string name;
};

/// Build a Scenario by name. `data_path` is required when the scenario reads
/// a tabulated input file (e.g. "pfss_field_line"); ignored otherwise. Throws
/// std::runtime_error on unknown name or missing required data_path.
Scenario make_scenario(const std::string& name, const std::string& data_path = "");

} // namespace chromosphere
