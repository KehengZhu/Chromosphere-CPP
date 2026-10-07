/*!
 * @file legacy/scenarios/legacy_scenarios.hpp
 * @brief LEGACY: name dispatch for the frozen two-fluid scenarios, and the open
 *        boundary update they share.
 * @ingroup scenarios
 *
 * Built only with -DCHROMO_ENABLE_LEGACY=ON. make_scenario() forwards every
 * legacy scenario name here; without the option it rejects those names with a
 * message naming the option.
 */

#pragma once

#include <string>

#include "scenarios/scenario.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Build a legacy two-fluid Scenario by name: `model_c7`, `model_flare`,
/// `model_gentle`, `analytic_canopy` or `pfss_field_line` (which requires
/// `data_path`). Applies each scenario's override-preserving environment
/// preset. Throws std::runtime_error for any other name or a missing
/// `data_path`.
Scenario make_legacy_scenario(const std::string& name, const std::string& data_path);

/// Open-hyperbolic boundary update shared by the legacy two-fluid scenarios
/// (model_c7, model_flare, analytic_canopy, pfss_field_line).
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
