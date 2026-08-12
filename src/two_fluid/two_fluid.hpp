/*!
 * two_fluid.hpp — the HISTORICAL two-fluid / multi-temperature /
 * finite-rate-ionization research solver. NOT the release path.
 *
 * The conserved state is the seven-row ion/neutral/electron carrier vector
 * `cons::` declared in chromosphere.hpp. Every entry point below rejects a Grid
 * that carries a Gamma1 table, i.e. a Grid configured for the single-fluid
 * release solver in single_fluid/mixture.hpp.
 *
 * Dependency direction: this header depends only on the shared infrastructure
 * in chromosphere.hpp. It never includes or calls the single-fluid solver, and
 * the single-fluid solver never includes or calls anything declared here.
 *
 * This header contains code derived from CoMFi (MIT License).
 */

#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

// Used by ip1/im1/...: SLICE means "treat input as one scalar field of length
// ns", CUBE means "treat input as the packed seven-row state".
const arma::uword SLICE = 1;
const arma::uword CUBE  = num_of_eq;

// Outer/inner ghost-cell handling for ip2/im2: copy interior (Neumann) when
// true, otherwise pull from the two-fluid ghost buffers in Grid.
const bool USE_NEUMANN_BC = false;


// ============================================================================
// LEGACY TWO-FLUID SOLVER (non-release)
//
// Everything below advances the seven-row ion/neutral/electron carrier state.
// It is the historical research path (two-fluid drift, three-temperature
// electrons, finite-rate ionization) and is NOT used by the canonical
// `model_column` release, which runs the single-fluid solver in `mixture.hpp`.
// These entry points reject a Grid that carries a Gamma1 table.
// ============================================================================

// ----------------------------------------------------------------------------
// Packed-state helpers
// ----------------------------------------------------------------------------

/// Extract one scalar field (length ns) from a packed state.
Vec get_scalar(const Grid& grid, const Vec& xn_state, arma::uword index);

/// Insert a scalar field (length ns) into the given slot of a packed state
/// (length n_state, other slots zero).
Vec scalar_to(const Grid& grid, const Vec& xn_i, arma::uword index);

/// Index shifts. `nk = CUBE` treats the input as packed state; `nk = SLICE`
/// treats it as one scalar field of length ns.
Vec ip1(const Grid& grid, const Vec& xn, arma::uword nk = CUBE);
Vec im1(const Grid& grid, const Vec& xn, arma::uword nk = CUBE);
Vec ip2(const Grid& grid, const Vec& xn, arma::uword nk = CUBE);
Vec im2(const Grid& grid, const Vec& xn, arma::uword nk = CUBE);

// ============================================================================
// Equation of state
// ============================================================================

Vec cons2prim(const Grid& grid, const Vec& cons_state);
Vec prim2cons(const Grid& grid, const Vec& prim_state);

// ============================================================================
// Numerics: flux, source, spectral radius, MUSCL limiter
// ============================================================================

/// Minmod limiter (writeup eq 15): φ(r) = max(0, min(1, r)).
Vec flux_lim(const Vec& r);

/// MC3 / Koren limiter (BATSRUS ModFaceValue 'mc3') — the ASYMMETRIC third-order
/// (κ=1/3) monotonized-central limiter. In the code's ratio form r = Δ₋/Δ₊ the two
/// faces of a cell take DIFFERENT limited slopes (unlike the symmetric minmod):
///   '+' faces (Lxn = u_i + ½φΔ₊):  φ₊(r) = max(0, min(β r, β, (2r+1)/3))
///   '−' faces (Rxn = u_i − ½φΔ₊):  φ₋(r) = max(0, min(β r, β, (r+2)/3))
/// β = grid.limiter_beta (2 ⇒ classic Koren). Non-finite r ⇒ 0 (first order at
/// flats/extrema, matching flux_lim). Selected by grid.mc3_limiter in rhs.
Vec flux_lim_mc3_plus (const Vec& r, float beta);
Vec flux_lim_mc3_minus(const Vec& r, float beta);

/// Cell-centered flux F(U).
Vec cal_flux_state(const Grid& grid, const Vec& xn_state);

/// Per-cell spectral radius ρ(∂F/∂U) broadcast across all num_of_eq equations.
Vec cal_spectral_radius_state(const Grid& grid, const Vec& xn_state);

/// Explicit source term (pressure-area + gravity).
Vec cal_source_state(const Grid& grid, const Vec& xn_state);

// ============================================================================
// RHS of the conservation laws (writeup §4)
// ============================================================================

/// Explicit RHS: TVD-MUSCL + Rusanov flux differencing + cal_source_state.
Vec rhs_explicit_state(const Grid& grid, const Vec& xn_state);

/// Implicit RHS: ion-neutral drag, collisional + frictional heating, and
/// field-aligned heat conduction (writeup §4.3–4.4).
Vec rhs_implicit_state(const Grid& grid, const Vec& xn_state);

// ============================================================================
// Time stepping
// ============================================================================

/// Per-cell maximum signal speed (extracted from the spectral radius broadcast).
Vec cal_max_v_i(const Grid& grid, const Vec& xn_state);

/// Uniform CFL-limited timestep.
Vec cal_dt_i(const Grid& grid, const Vec& xn_state);

/// Semi-implicit (backward-Euler) integrator. Explicit MUSCL+Rusanov step
/// for R_E, then point-implicit drag + frictional heating, point-implicit
/// ion–neutral temperature equilibration, a tridiagonal heat-conduction
/// solve per species, and (if grid.enable_ionization) the point-implicit
/// ionization Stage E (writeup §3.7, §5.3). Mutates grid.dt_state.
Vec advance_Euler_state(Grid& grid, const Vec& xn_state, const Vec& dt_i);

/// Stage E: backward-Euler ionization / recombination on a single cell, scalar
/// quadratic solve in ionization fraction f = ρ_i / (ρ_i + ρ_n). Operates on
/// primitive state in place; intended for the operator-split integrator but
/// exposed for direct testing. (writeup §5.3.)
void apply_ionization_stage(const Grid& grid, Vec& prim_state, float dt);

/// Stage R (writeup §2.4): backward-Euler optically-thick chromospheric
/// radiative cooling on the folded ion+electron thermal pressure. Carlsson
/// & Leenaarts 2012 recipe summed over H I + Ca II + Mg II; tables in
/// physics.hpp::cl2012. Activated by Grid::enable_radiative_cooling.
void apply_radiative_cooling_stage(const Grid& grid, Vec& prim_state, float dt);

/// Flare beam-heating stage: deposits grid.beam_flux into the upper-chromosphere
/// electron/ion thermal pool (p_i) over a finite layer, gated by the temporal
/// window around grid.sim_time, to drive chromospheric evaporation (Fisher et
/// al. 1985). No-op unless grid.enable_beam_heating. Mutates p_i in place.
void apply_beam_heating_stage(const Grid& grid, Vec& prim_state, float dt);

/// Ambient coronal heating stage: deposits the steady footpoint-anchored
/// volumetric heating H(s) (physics.hpp::coronal_heating_rate) into the charged
/// thermal pool (shared with neutrals by heat capacity; electrons when
/// enable_Te), driving the conductive flux that sustains the corona and, when
/// ramped up (coronal_heat_enhance > 1), gentle chromospheric evaporation
/// (Antiochos & Sturrock 1978). No-op unless grid.enable_coronal_heating.
void apply_coronal_heating_stage(const Grid& grid, Vec& prim_state, float dt);

/// TRAC adaptive cutoff temperature T_c (Johnston et al. 2020 Eq. 8): the
/// maximum temperature among grid cells where the charged-fluid TR is
/// under-resolved (L_R/L_T > 1/2, with L_T = T/|dT/ds|, L_R = Δs), clamped to
/// [grid.trac_T_chrom, 0.2 T_peak]. Used by the conduction and cooling stages to
/// broaden the unresolved TR. Exposed for testing.
float compute_trac_cutoff_T(const Grid& grid, const Vec& prim_state);

/// Pure-explicit forward-Euler step — only the MUSCL+Rusanov R_E predictor of
/// advance_Euler_state, with the implicit drag / temperature / conduction
/// stages skipped (equivalent to R_I ≡ 0). Mutates grid.dt_state.
Vec advance_Euler_explicit_state(Grid& grid, const Vec& xn_state, const Vec& dt_i);

/// Explicit RK4. Mutates grid.dt_state.
Vec advance_RK4(Grid& grid, const Vec& xn_state, const Vec& dt_i);

// ============================================================================
// Debug
// ============================================================================

void print_xn(const Grid& grid, const Vec& xn);

} // namespace chromosphere
