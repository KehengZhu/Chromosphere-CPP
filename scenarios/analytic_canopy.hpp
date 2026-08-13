/*!
 * @file scenarios/analytic_canopy.hpp
 * @brief Exponential magnetic-canopy B(z) overlay on the C7 thermodynamic profile.
 * @ingroup scenarios
 *
 * Analytic canopy scenario.
 *
 * C7 thermodynamic profile (T, n_e, n_n) overlaid with the consensus
 * exponential-canopy magnetic-field recipe
 *
 *     B(z) = B_inf + (B_0 - B_inf) * exp(-(z - z_base) / H_B)
 *
 * Defaults (see Bellot Rubio & Orozco Suárez 2019; Kontar et al. 2008;
 * Martínez-Sykora et al. 2019):
 *   B_0   = 100 G  = 1.0e-2 T   (footpoint network concentration)
 *   B_inf =  15 G  = 1.5e-3 T   (canopy-merged chromospheric value)
 *   H_B   = 300 km                (chromospheric magnetic scale height)
 *
 * Domain and IC otherwise identical to model_c7, so this scenario differs
 * from model_c7 only through the ∂(1/B)/∂s source in flux.cpp:74.
 */

#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Grid resolution this scenario expects, queried before Grid::init(). The
/// canopy shares the model_c7 domain, so it uses the same fixed 100 cells.
arma::uword analytic_canopy_peek_ns();

/// Initial condition: run model_c7_ic (thermodynamics, gravity, ghost cells and
/// the conserved state), then overwrite the magnetic geometry with the
/// exponential canopy B(s) — face values B_imh / B_iph, the cell average B_i and
/// the flux-tube source coefficient d(1/B)/ds — and re-broadcast the packed Grid
/// caches. The field is vertical here, so arc length s equals height above the
/// C7 base. Returns the (unmodified) initial conserved state.
Vec  analytic_canopy_ic(Grid& grid);

/// Boundary refresh, delegated verbatim to model_c7_update_bc: the outer BC is a
/// function of (n, T, V) only and no term in it involves B.
void analytic_canopy_update_bc(Grid& grid, const Vec& xn);

/** @} */

} // namespace chromosphere
