/*!
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

// Static local-refinement controls. Read from the environment; the documented
// diagnostic profile refines the lower domain 0–700 km by 4× and transitions
// back to the outer coarse spacing smoothly over 700–800 km. Coordinates are
// field-line distance [km] from the domain INNER face (s = 0 at the inner
// boundary), so the 0–700 km region matches the photosphere-anchored scenarios.
struct RefineParams {
    double factor        = 1.0;    // fine/coarse resolution ratio; 1 ⇒ disabled
    double s_lo_km       = 0.0;    // inner edge of the refined region [km]
    double s_hi_km       = 700.0;  // outer edge of the refined region [km]
    double transition_km = 100.0;  // width of the fine→coarse transition [km]

    bool enabled() const { return factor > 1.0 + 1.0e-9; }
};

// Read GRID_REFINE_FACTOR / _S_LO_KM / _S_HI_KM / _TRANSITION_KM. When
// `iso_alias` is true the ISO_REFINE_* variables OVERRIDE their GRID_REFINE_*
// counterparts (model_column). Missing variables keep the struct defaults.
RefineParams refine_params_from_env(bool iso_alias);

// Build canonical cell FACES [m], measured from the inner face: the returned
// vector has (ncells+1) entries, faces.front() == 0 and faces.back() == length_m
// exactly. `ns_coarse` is the coarse-grid-equivalent cell count; the coarse
// spacing is length_m/ns_coarse and refinement only ADDS lower-domain cells (it
// never coarsens the rest of the column). When rp is disabled — or the refined
// window does not fit inside the domain — a uniform ns_coarse grid is returned.
//
// The refined column is [0,s_lo] coarse, [s_lo,s_hi] fine (≈coarse/factor),
// [s_hi,s_hi+transition] an exponential (geometric) grade with adjacent
// cell-width ratio ≤ MESH_MAX_RATIO, then [s_hi+transition,length] coarse. The
// segment boundaries 0, s_lo, s_hi and s_hi+transition are exact faces.
std::vector<double> build_static_mesh_faces(double length_m,
                                            arma::uword ns_coarse,
                                            const RefineParams& rp);

// Adjacent cell-width ratio cap enforced across the transition segment.
constexpr double MESH_MAX_RATIO = 1.1;

} // namespace chromosphere
