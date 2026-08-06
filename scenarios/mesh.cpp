#include "mesh.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
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

// `lower` | `outer`, case-insensitive; anything else keeps `fallback`.
RefineProfile env_profile(const char* key, RefineProfile fallback) {
    const char* e = std::getenv(key);
    if (!e) return fallback;
    std::string v(e);
    for (char& c : v) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (v == "outer") return RefineProfile::Outer;
    if (v == "lower") return RefineProfile::Lower;
    return fallback;
}

// Append a UNIFORM segment (seg_start, seg_end] whose cells are AT MOST
// `target_ds` wide (ceil, so a segment is never coarser than its target — that is
// what keeps the ratio at a junction with a graded segment bounded). Both
// endpoints land on faces; the starting face is assumed already present as
// faces.back() and the last appended face is snapped to seg_end exactly.
void append_uniform(std::vector<double>& faces, double seg_start, double seg_end,
                    double target_ds) {
    const double span = seg_end - seg_start;
    const arma::uword n = static_cast<arma::uword>(
        std::max(1.0, std::ceil(span / std::max(target_ds, 1.0e-30) - 1.0e-9)));
    const double w = span / static_cast<double>(n);
    for (arma::uword j = 1; j <= n; ++j)
        faces.push_back(seg_start + w * static_cast<double>(j));
    faces.back() = seg_end;
}

// An exponential (geometric) grade bridging cell width `fine_ds` to `coarse_ds`.
//
// Cell widths are fine_ds·r^m, m = 0…n-1, with r = (coarse_ds/fine_ds)^{1/n}, so
// the narrow end matches `fine_ds` EXACTLY and the wide end is coarse_ds/r. Every
// adjacent ratio — inside the grade and at BOTH junctions — is then r, which is
// ≤ MESH_MAX_RATIO by construction.
//
// That fixes both the width and the cell count: `width` is an OUTPUT, and the
// requested transition width acts as a MINIMUM. `n` is the smallest count whose
// grade is at least `tau_min_m` wide (the width fine_ds·(f−1)/(r−1) grows
// monotonically with n), so a requested width the ratio cap cannot support is
// widened to the nearest one it can, instead of silently breaking the cap.
struct Grade {
    arma::uword n     = 0;
    double      ratio = 1.0;
    double      width = 0.0;
};

Grade make_grade(double fine_ds, double coarse_ds, double tau_min_m) {
    Grade g;
    const double f = coarse_ds / fine_ds;
    if (!(f > 1.0 + 1.0e-9) || !(fine_ds > 0.0)) return g;       // no grade needed
    auto width_of = [&](arma::uword n) {
        const double r = std::pow(f, 1.0 / static_cast<double>(n));
        return fine_ds * (f - 1.0) / (r - 1.0);                  // Σ fine·r^m, m<n
    };
    // Ratio cap: r = f^{1/n} ≤ MESH_MAX_RATIO ⇒ n ≥ ln f / ln cap.
    arma::uword n = static_cast<arma::uword>(
        std::max(1.0, std::ceil(std::log(f) / std::log(MESH_MAX_RATIO))));
    while (width_of(n) < tau_min_m && n < 1000000u) ++n;
    g.n     = n;
    g.ratio = std::pow(f, 1.0 / static_cast<double>(n));
    g.width = width_of(n);
    return g;
}

// Append the graded transition of width `g.width` starting at `seg_start`, with
// the cells running fine→coarse when `outward_coarsening` and coarse→fine
// otherwise. Node offsets use the closed-form partial sums (no accumulation
// drift); the final face is snapped to seg_start + g.width exactly.
void append_grade(std::vector<double>& faces, double seg_start, double fine_ds,
                  const Grade& g, bool outward_coarsening) {
    const double c = fine_ds / (g.ratio - 1.0);      // Σ_{m<j} fine·r^m = c·(r^j − 1)
    for (arma::uword j = 1; j <= g.n; ++j) {
        const double from_fine_end =
            c * (std::pow(g.ratio, static_cast<double>(j)) - 1.0);
        faces.push_back(seg_start + (outward_coarsening ? from_fine_end
                                                        : g.width - c *
              (std::pow(g.ratio, static_cast<double>(g.n - j)) - 1.0)));
    }
    faces.back() = seg_start + g.width;
}

} // namespace

RefineParams refine_params_from_env(bool iso_alias) {
    RefineParams rp;
    rp.factor        = env_double("GRID_REFINE_FACTOR", rp.factor);
    rp.s_lo_km       = env_double("GRID_REFINE_S_LO_KM", rp.s_lo_km);
    rp.s_hi_km       = env_double("GRID_REFINE_S_HI_KM", rp.s_hi_km);
    rp.transition_km = env_double("GRID_REFINE_TRANSITION_KM", rp.transition_km);
    rp.profile       = env_profile("GRID_REFINE_PROFILE", rp.profile);
    if (iso_alias) {
        rp.factor        = env_double("ISO_REFINE_FACTOR", rp.factor);
        rp.s_lo_km       = env_double("ISO_REFINE_S_LO_KM", rp.s_lo_km);
        rp.s_hi_km       = env_double("ISO_REFINE_S_HI_KM", rp.s_hi_km);
        rp.transition_km = env_double("ISO_REFINE_TRANSITION_KM", rp.transition_km);
        rp.profile       = env_profile("ISO_REFINE_PROFILE", rp.profile);
    }
    return rp;
}

std::vector<double> build_static_mesh_faces(double length_m, arma::uword ns_coarse,
                                            const RefineParams& rp) {
    const double coarse_ds = length_m / static_cast<double>(ns_coarse);

    const double s_lo    = rp.s_lo_km * 1000.0;
    const double s_hi    = rp.s_hi_km * 1000.0;
    const double fine_ds = coarse_ds / rp.factor;
    const bool   outer   = (rp.profile == RefineProfile::Outer);

    // The graded transition's width is set by the ratio cap; the requested
    // ISO_REFINE_TRANSITION_KM is its MINIMUM (see make_grade).
    const Grade grade = (rp.enabled() && fine_ds > 0.0)
                      ? make_grade(fine_ds, coarse_ds, rp.transition_km * 1000.0)
                      : Grade{};
    const double tau   = grade.width;
    const double s_end = s_hi + tau;                          // top of transition (Lower)
    const double s_mid = s_lo - tau;                          // foot of transition (Outer)

    // Uniform fallback: refinement off, or the refined+transition window does not
    // fit strictly inside the domain, or a degenerate window. Byte-for-byte the
    // k·coarse_ds grid a scenario builds when refinement is disabled.
    const bool fits = rp.enabled() && fine_ds > 0.0 && grade.n > 0 &&
                      (outer ? (s_mid > 0.0 && s_lo < length_m)
                             : (s_lo >= 0.0 && s_hi > s_lo && s_end < length_m));
    if (!fits) {
        std::vector<double> faces(ns_coarse + 1);
        for (arma::uword k = 0; k <= ns_coarse; ++k)
            faces[k] = coarse_ds * static_cast<double>(k);
        faces.back() = length_m;
        return faces;
    }

    std::vector<double> faces;
    faces.push_back(0.0);

    if (outer) {
        // OUTER profile: coarse → grade coarse→fine → fine up to the domain top.
        append_uniform(faces, 0.0, s_mid, coarse_ds);                 // [0, s_lo − τ]
        append_grade(faces, s_mid, fine_ds, grade, /*outward_coarsening=*/false);
        append_uniform(faces, s_lo, length_m, fine_ds);               // [s_lo, length]
        return faces;
    }

    // [0, s_lo] — coarse (only present when the refined region does not start at
    // the inner face). Retains the outer coarse spacing below s_lo.
    if (s_lo > 0.0) append_uniform(faces, 0.0, s_lo, coarse_ds);

    // [s_lo, s_hi] — fine, ≈ coarse_ds/factor, exact endpoints.
    append_uniform(faces, s_lo, s_hi, fine_ds);

    // [s_hi, s_hi + τ] — exponential transition fine→coarse.
    append_grade(faces, s_hi, fine_ds, grade, /*outward_coarsening=*/true);

    // [s_hi + τ, length] — coarse, retains the original outer spacing.
    append_uniform(faces, s_end, length_m, coarse_ds);

    return faces;
}

} // namespace chromosphere
