#include "mesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>

namespace chromosphere {

namespace {

double env_double(const char* key, double fallback) {
    if (const char* e = std::getenv(key)) {
        try { return std::stod(std::string(e)); } catch (...) {}
    }
    return fallback;
}

// Append a UNIFORM segment (seg_start, seg_end] discretized into `n` cells of
// exactly (seg_end-seg_start)/n so both endpoints land on faces. The starting
// face (seg_start) is assumed already present as faces.back(); only the n
// interior/end faces are appended, with the last snapped to seg_end exactly.
void append_uniform(std::vector<double>& faces, double seg_start, double seg_end,
                    arma::uword n) {
    const double w = (seg_end - seg_start) / static_cast<double>(n);
    for (arma::uword j = 1; j <= n; ++j)
        faces.push_back(seg_start + w * static_cast<double>(j));
    faces.back() = seg_end;
}

// Append the fine→coarse TRANSITION segment (seg_start, seg_end] via an
// exponential cell-density mapping. With stretch k = ln(ds_end/ds_start) the
// node j sits at seg_start + tau·(e^{k j/N} − 1)/(e^k − 1), so cell widths grow
// geometrically with adjacent ratio (ds_end/ds_start)^{1/N}. N is chosen large
// enough to (a) keep that ratio ≤ MESH_MAX_RATIO and (b) make the first/last
// widths continuous with the fine/coarse neighbours (enough cells for the ratio
// cap). Endpoints seg_start and seg_end are exact faces by construction.
void append_transition(std::vector<double>& faces, double seg_start, double seg_end,
                       double ds_start, double ds_end) {
    const double tau = seg_end - seg_start;
    // Degenerate / non-growing transition ⇒ a uniform bridge at the mean width.
    if (!(ds_end > ds_start * (1.0 + 1.0e-9)) || tau <= 0.0) {
        arma::uword n = static_cast<arma::uword>(
            std::max(1.0, std::round(tau / std::max(ds_start, 1.0e-30))));
        append_uniform(faces, seg_start, seg_end, n);
        return;
    }
    const double k     = std::log(ds_end / ds_start);
    // Continuity-preserving count (first width ≈ ds_start): N = k / ln(1+(b-a)/tau).
    const double n_cont = k / std::log1p((ds_end - ds_start) / tau);
    // Ratio-cap count: (ds_end/ds_start)^{1/N} ≤ MESH_MAX_RATIO ⇒ N ≥ k/ln(cap).
    const double n_cap  = k / std::log(MESH_MAX_RATIO);
    arma::uword N = static_cast<arma::uword>(
        std::max(1.0, std::ceil(std::max(n_cont, n_cap))));
    const double denom = std::expm1(k);   // e^k − 1 = ds_end/ds_start − 1
    for (arma::uword j = 1; j <= N; ++j) {
        const double frac = static_cast<double>(j) / static_cast<double>(N);
        faces.push_back(seg_start + tau * std::expm1(k * frac) / denom);
    }
    faces.back() = seg_end;
}

} // namespace

RefineParams refine_params_from_env(bool iso_alias) {
    RefineParams rp;
    rp.factor        = env_double("GRID_REFINE_FACTOR", rp.factor);
    rp.s_lo_km       = env_double("GRID_REFINE_S_LO_KM", rp.s_lo_km);
    rp.s_hi_km       = env_double("GRID_REFINE_S_HI_KM", rp.s_hi_km);
    rp.transition_km = env_double("GRID_REFINE_TRANSITION_KM", rp.transition_km);
    if (iso_alias) {
        rp.factor        = env_double("ISO_REFINE_FACTOR", rp.factor);
        rp.s_lo_km       = env_double("ISO_REFINE_S_LO_KM", rp.s_lo_km);
        rp.s_hi_km       = env_double("ISO_REFINE_S_HI_KM", rp.s_hi_km);
        rp.transition_km = env_double("ISO_REFINE_TRANSITION_KM", rp.transition_km);
    }
    return rp;
}

std::vector<double> build_static_mesh_faces(double length_m, arma::uword ns_coarse,
                                            const RefineParams& rp) {
    const double coarse_ds = length_m / static_cast<double>(ns_coarse);

    const double s_lo  = rp.s_lo_km * 1000.0;
    const double s_hi  = rp.s_hi_km * 1000.0;
    const double s_end = s_hi + rp.transition_km * 1000.0;   // top of transition
    const double fine_ds = coarse_ds / rp.factor;

    // Uniform fallback: refinement off, or the refined+transition window does not
    // fit strictly inside the domain, or a degenerate window. Byte-for-byte the
    // k·coarse_ds grid a scenario builds when refinement is disabled.
    const bool fits = rp.enabled() && fine_ds > 0.0 &&
                      s_lo >= 0.0 && s_hi > s_lo && s_end > s_hi && s_end < length_m;
    if (!fits) {
        std::vector<double> faces(ns_coarse + 1);
        for (arma::uword k = 0; k <= ns_coarse; ++k)
            faces[k] = coarse_ds * static_cast<double>(k);
        faces.back() = length_m;
        return faces;
    }

    std::vector<double> faces;
    faces.push_back(0.0);

    // [0, s_lo] — coarse (only present when the refined region does not start at
    // the inner face). Retains the outer coarse spacing below s_lo.
    if (s_lo > 0.0) {
        arma::uword n = static_cast<arma::uword>(std::max(1.0, std::round(s_lo / coarse_ds)));
        append_uniform(faces, 0.0, s_lo, n);
    }

    // [s_lo, s_hi] — fine, ≈ coarse_ds/factor, exact endpoints.
    {
        arma::uword n = static_cast<arma::uword>(
            std::max(1.0, std::round((s_hi - s_lo) / fine_ds)));
        append_uniform(faces, s_lo, s_hi, n);
    }

    // [s_hi, s_end] — exponential transition fine→coarse, ratio ≤ MESH_MAX_RATIO.
    append_transition(faces, s_hi, s_end, fine_ds, coarse_ds);

    // [s_end, length] — coarse, retains the original outer spacing.
    {
        arma::uword n = static_cast<arma::uword>(
            std::max(1.0, std::round((length_m - s_end) / coarse_ds)));
        append_uniform(faces, s_end, length_m, n);
    }

    return faces;
}

} // namespace chromosphere
