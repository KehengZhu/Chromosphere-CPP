/*!
 * @file scenarios/mesh.hpp
 * @brief Shared STATIC local-refinement mesh construction (not AMR).
 * @ingroup scenarios
 *
 * Shared static local-refinement mesh construction.
 *
 * Builds a metric-aware, STATIC non-uniform finite-volume face grid used by the
 * built-in static scenarios. This is static local refinement, NOT AMR: the mesh
 * is generated once at IC time and never regridded, remapped, or adapted at
 * runtime. PFSS field lines carry their own externally-specified mesh and never
 * go through this builder.
 *
 * The refinement is opt-in and reduces to the exact uniform grid when disabled
 * (factor 1), so uniform runs stay byte-for-byte identical.
 */

#pragma once

#include <armadillo>
#include <vector>

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Which end of the domain the refined band is anchored to.
enum class RefineProfile {
    /// coarse → fine [s_lo,s_hi] → graded fine→coarse transition → coarse.
    /// The historical (lower-domain) profile; `s_hi` and the transition ABOVE it
    /// are both meaningful.
    Lower,
    /// coarse → graded coarse→fine transition [s_lo−τ,s_lo] → fine [s_lo,L].
    /// A compact refined band that reaches the OUTER boundary (the upper TR).
    /// `s_hi_km` is unused: the fine region always ends at the domain top, and
    /// the transition sits immediately BELOW `s_lo`.
    Outer,
};

/// Static local-refinement controls. Read from the environment; the documented
/// lower-domain diagnostic profile refines 0–700 km by 4× and transitions back to
/// the outer coarse spacing smoothly over 700–800 km, while the documented
/// outer/TR profile refines from s_lo up to the domain top with the transition in
/// [s_lo−τ, s_lo]. Coordinates are field-line distance [km] from the domain INNER
/// face (s = 0 at the inner boundary).
struct RefineParams {
    double factor        = 1.0;  ///< fine/coarse resolution ratio; 1 ⇒ disabled
    double s_lo_km       = 0.0;  ///< inner edge of the refined region [km]
    double s_hi_km       = 700.0;  ///< outer edge of the refined region [km] (Lower only)
    double transition_km = 100.0;  ///< width of the graded coarse↔fine transition [km]
    /// Which end of the domain the refined band is anchored to, i.e. how
    /// s_lo_km / s_hi_km / transition_km are laid out (see RefineProfile).
    RefineProfile profile = RefineProfile::Lower;

    /// True when refinement is actually requested. `factor` is compared against
    /// 1 with a small tolerance, so a nominally-1 factor read from the
    /// environment still yields the exact uniform mesh.
    bool enabled() const { return factor > 1.0 + 1.0e-9; }
};

/// Read `GRID_REFINE_FACTOR` / `_S_LO_KM` / `_S_HI_KM` / `_TRANSITION_KM` /
/// `_PROFILE` (`lower` | `outer`, case-insensitive). When `iso_alias` is true the
/// ISO_REFINE_* variables OVERRIDE their GRID_REFINE_* counterparts
/// (model_column). Missing variables keep the struct defaults.
RefineParams refine_params_from_env(bool iso_alias);

/// Build canonical cell FACES [m], measured from the inner face: the returned
/// vector has (ncells+1) entries, faces.front() == 0 and faces.back() == length_m
/// exactly. `ns_coarse` is the coarse-grid-equivalent cell count; the coarse
/// spacing is length_m/ns_coarse and refinement only ADDS cells (it never coarsens
/// the rest of the column). When rp is disabled — or the refined window does not
/// fit inside the domain — a uniform ns_coarse grid is returned.
///
/// RefineProfile::Lower: [0,s_lo] coarse, [s_lo,s_hi] fine (≈coarse/factor),
/// [s_hi,s_hi+transition] an exponential (geometric) grade fine→coarse, then
/// [s_hi+transition,length] coarse. Exact faces at 0, s_lo, s_hi, s_hi+transition.
///
/// RefineProfile::Outer: [0,s_lo−transition] coarse, [s_lo−transition,s_lo] an
/// exponential grade coarse→fine, then [s_lo,length] fine (≈coarse/factor) all the
/// way to the outer boundary. Exact faces at 0, s_lo−transition, s_lo, length.
///
/// Both grades keep the adjacent cell-width ratio ≤ MESH_MAX_RATIO.
std::vector<double> build_static_mesh_faces(double length_m,
                                            arma::uword ns_coarse,
                                            const RefineParams& rp);

/// Adjacent cell-width ratio cap enforced across the transition segment.
constexpr double MESH_MAX_RATIO = 1.1;

/** @} */

} // namespace chromosphere
