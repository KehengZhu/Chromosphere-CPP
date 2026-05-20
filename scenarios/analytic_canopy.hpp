/*!
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

arma::uword analytic_canopy_peek_ns();

Vec  analytic_canopy_ic(Grid& grid);
void analytic_canopy_update_bc(Grid& grid, const Vec& xn);

} // namespace chromosphere
