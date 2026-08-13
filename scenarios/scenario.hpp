/*!
 * @file scenarios/scenario.hpp
 * @brief Scenario dispatcher and the shared open-hyperbolic boundary update.
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

#include "chromosphere.hpp"

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

/// Build a Scenario by name. `data_path` is required when the scenario reads
/// a tabulated input file (e.g. "pfss_field_line"); ignored otherwise. Throws
/// std::runtime_error on unknown name or missing required data_path.
Scenario make_scenario(const std::string& name, const std::string& data_path = "");

/// Open-hyperbolic boundary update shared by all built-in scenarios.
///
///   * Inner face (s = s_{-1/2}) — closed reflecting wall:
///       ρ(-1) = ρ(0),  T(-1) = T(0),  V(-1) = -V(0)
///       ρ(-2) = ρ(1),  T(-2) = T(1),  V(-2) = -V(1)
///     Anti-symmetric velocity guarantees zero mass flux across s = 0; the
///     ρ and T mirrors keep the pressure gradient there well posed.
///
///   * Outer face (s = s_{ns-1/2}) — pure Neumann outflow on every variable:
///       ρ(ns)   = ρ(ns+1) = ρ(ns-1)
///       T(ns)   = T(ns+1) = T(ns-1)
///       V(ns)   = V(ns+1) = V(ns-1)
///       U(ns)   = U(ns+1) = U(ns-1)
///     No impedance / dissipation at the face — well-posed only for IC near
///     HSE or with an added sponge layer; transient runs (e.g. C7 IC under
///     gravity) may go unstable.
///
/// Ghost energies are rebuilt with the same gauge as the IC: φ_g(-1) =
/// grid.phi_g_imh[0] and φ_g(ns) = grid.phi_g_iph[ns-1].
void apply_open_bcs(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
