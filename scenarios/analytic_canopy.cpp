/*!
 * @file scenarios/analytic_canopy.cpp
 * @brief Exponential magnetic-canopy B(z) overlay on the C7 thermodynamic profile.
 * @ingroup scenarios
 */
#include "analytic_canopy.hpp"
#include "model_c7.hpp"

#include <cmath>

namespace chromosphere {

namespace {
// Consensus chromospheric canopy parameters. See header for citations.
constexpr Real kB0_T   = 1.0e-2f;   // 100 G footpoint
constexpr Real kBinf_T = 1.5e-3f;   // 15 G canopy-merged
constexpr Real kHB_m   = 3.0e5f;    // 300 km scale height
} // namespace

arma::uword analytic_canopy_peek_ns() { return 100; }

Vec analytic_canopy_ic(Grid& grid) {
    // Reuse the full C7 setup (thermodynamics, gravity, ghost cells, conserved
    // state). model_c7_ic leaves B=1, dinvB_ds=0; we now override.
    Vec xn = model_c7_ic(grid);

    // Walk up the field line using grid.ds_i to get the height at each face.
    // Field is vertical here, so arc length s == height above the C7 base.
    Real s_imh = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const Real s_iph = s_imh + grid.ds_i(i);
        const Real Bimh  = kBinf_T + (kB0_T - kBinf_T) * std::exp(-s_imh / kHB_m);
        const Real Biph  = kBinf_T + (kB0_T - kBinf_T) * std::exp(-s_iph / kHB_m);
        grid.B_imh(i)      = Bimh;
        grid.B_iph(i)      = Biph;
        grid.B_i (i)       = 0.5f * (Bimh + Biph);
        grid.dinvB_ds_i(i) = (1.0f / Biph - 1.0f / Bimh) / grid.ds_i(i);
        s_imh = s_iph;
    }

    grid.broadcast();  // re-sync the packed _state caches after B/dinvB changes
    return xn;
}

void analytic_canopy_update_bc(Grid& grid, const Vec& xn) {
    // The outer BC depends only on (n, T, V) — none of which involve B —
    // so the C7 refresh logic transfers verbatim.
    model_c7_update_bc(grid, xn);
}

} // namespace chromosphere
