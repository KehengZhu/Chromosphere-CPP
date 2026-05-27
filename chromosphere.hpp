/*!
 * Chromosphere model — public API.
 *
 * All solver state lives in a Grid struct that the caller owns; functions
 * take `const Grid&` (or `Grid&` for those that mutate scratch buffers).
 * No global mutable state — multiple Grids can coexist.
 *
 * This header contains code derived from CoMFi (MIT License).
 */

#pragma once

#include <armadillo>

namespace chromosphere {

/// Armadillo float column vector — the workhorse type for cell-wise data.
typedef arma::Col<float> Vec;

// ============================================================================
// State indices
// ============================================================================

const arma::uword num_of_eq = 6;

// Conserved-variable indices (writeup eq 61).
namespace cons {
    const arma::uword RHO_I = 0; // ion mass density       ρ_i
    const arma::uword RHO_N = 1; // neutral mass density   ρ_n
    const arma::uword MOM_I = 2; // ion momentum density   ρ_i V
    const arma::uword MOM_N = 3; // neutral momentum dens. ρ_n U
    const arma::uword E_I   = 4; // ion total energy       e_i
    const arma::uword E_N   = 5; // neutral total energy   e_n
}

// Primitive-variable indices.
namespace prim {
    const arma::uword RHO_I = 0;
    const arma::uword RHO_N = 1;
    const arma::uword V     = 2;
    const arma::uword U     = 3;
    const arma::uword P_I   = 4;
    const arma::uword P_N   = 5;
}

// Used by ip1/im1/...: SLICE means "treat input as one scalar field of length ns",
// CUBE means "treat input as the packed state of length ns*num_of_eq".
const arma::uword SLICE = 1;
const arma::uword CUBE  = num_of_eq;

// Outer/inner ghost-cell handling for ip2/im2: copy interior (Neumann) when true,
// otherwise pull from the ghost buffers in Grid.
const bool USE_NEUMANN_BC = false;

// ============================================================================
// Grid: owns all solver state.
// ============================================================================

struct Grid {
    // --- sizes -------------------------------------------------------------
    arma::uword ns      = 0;
    arma::uword n_state = 0;            // ns * num_of_eq
    float       CFL     = 0.25f;

    // --- physical constants (SI) ------------------------------------------
    float gamma_mono = 5.0f / 3.0f;
    float m_i        = 1.6726219e-27f;
    float m_n        = 1.6726219e-27f;
    float m_e        = 9.10938356e-31f;
    float g          = 0.27395e3f;
    float mu_0       = 4.0f * static_cast<float>(arma::datum::pi) * 1.0e-7f;
    float k_b        = 1.380649e-23f;
    float q_e        = 1.602176634e-19f;
    float chi_H_J    = 2.179872361e-18f;     // hydrogen ionization potential, 13.6 eV in J

    // --- runtime toggles ---------------------------------------------------
    // Enable Stage E (hydrogen ionization / recombination) in advance_Euler_state.
    // Default off so the existing test suite remains a clean regression baseline;
    // flip to true to activate the writeup §5.3 ionization stage.
    bool enable_ionization = false;

    // Enable optically-thick chromospheric radiative cooling (CL2012 recipe,
    // writeup §2.4) as Stage R in the operator-split integrator. Default off
    // so the existing test suite remains a clean regression baseline; flip
    // to true to activate.
    bool enable_radiative_cooling = false;

    // Effective hydrogen photoionization rate from the ground state [s^-1],
    // modeling the Lyα-excitation + Balmer-continuum two-step channel
    // (Carlsson & Stein 2002 establish this as the dominant chromospheric
    // ionization path; the Lyman continuum itself contributes negligibly).
    // C&S 2002 report a chromospheric ionization/recombination timescale of
    // 10^3–10^5 s; we adopt the geometric mean ~10^4 s as a single-rate
    // closure suitable for the 2-state (ground+continuum) reduction used
    // here. See Leenaarts et al. 2007 for the "fixed radiative rates"
    // implementation paradigm in non-equilibrium chromospheric MHD codes.
    //
    // The scalar `photoionization_rate` is the uniform fallback used when
    // the cell-centered `photoionization_rate_i` array is empty (i.e. no
    // height-dependent calibration has been performed by the scenario IC).
    // Scenarios that want a height-dependent profile (e.g. model_c7_ic,
    // which calibrates against C7's tabulated f to make C7 a fixed point
    // of the Stage E quadratic) populate `photoionization_rate_i` of
    // length `ns`; the Stage E driver picks that up automatically via
    // `physics.hpp::photoionization_rate_P(grid)`.
    float photoionization_rate = 1.0e-4f;
    Vec   photoionization_rate_i;

    // --- cell-centered & face arrays (length ns) --------------------------
    Vec ds_i;
    Vec B_i, B_imh, B_iph;
    Vec dinvB_ds_i;
    Vec phi_g_imh, phi_g_iph;

    // --- ghost buffers (length num_of_eq) ---------------------------------
    Vec outer_boundary0_i, outer_boundary1_i;
    Vec inner_boundary0_i, inner_boundary1_i;

    // --- packed broadcasts of the above (length n_state) ------------------
    // Cached versions of B_i, B_imh, B_iph, dinvB_ds_i, ds_i, and a scratch
    // dt buffer broadcast across all num_of_eq equations. Recomputed by
    // broadcast() whenever the cell-centered fields change.
    Vec B_state, B_state_imh, B_state_iph;
    Vec dinvB_ds_state, ds_state, dt_state;

    /// Allocate arrays for ns cells, fill physical constants, zero the fields.
    void init(arma::uword ns_in, float CFL_in);

    /// Recompute the packed `_state` broadcasts from ds_i, B_i, B_imh, B_iph,
    /// dinvB_ds_i. Call after editing any of those fields.
    void broadcast();
};

// ============================================================================
// Packed-state helpers
// ============================================================================

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
