#include "chromosphere.hpp"
#include "physics.hpp"

#include <iostream>

namespace chromosphere {

// ============================================================================
// Time-step calculation
// ============================================================================

Vec cal_max_v_i(const Grid& grid, const Vec& xn_state) {
    Vec spc_r_state = cal_spectral_radius_state(grid, xn_state);
    return get_scalar(grid, spc_r_state, 0);
}

Vec cal_dt_i(const Grid& grid, const Vec& xn_state) {
    Vec dt_i = grid.CFL * grid.ds_i / cal_max_v_i(grid, xn_state);
    float dt_min_i = arma::min(dt_i);
    return dt_min_i * arma::ones<Vec>(grid.ns);
}

// ============================================================================
// Time integrators
// ============================================================================
//
// Both integrators mutate grid.dt_state because rhs_explicit_state reads it
// during the MUSCL prediction step. The pattern is: broadcast dt_i across all
// num_of_eq slots, then call the RHS.
//

static void broadcast_dt(Grid& grid, const Vec& dt_i) {
    grid.dt_state.zeros(grid.n_state);
    for (arma::uword i = 0; i < num_of_eq; ++i) {
        grid.dt_state += scalar_to(grid, dt_i, i);
    }
}

// ----------------------------------------------------------------------------
// Stage B (writeup §3.7 Eqs. (drag-w)--(drag-partition)): point-implicit drag
// + frictional heating. Exact for the linear drift ODE with α frozen at U*.
// Mutates V, U, p_i, p_n in place; ρ_i, ρ_n unchanged.
// ----------------------------------------------------------------------------
static void apply_drag_stage(const Grid& grid, Vec& prim_state, float dt) {
    Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    Vec V     = get_scalar(grid, prim_state, prim::V);
    Vec U     = get_scalar(grid, prim_state, prim::U);
    Vec p_i   = get_scalar(grid, prim_state, prim::P_I);
    Vec p_n   = get_scalar(grid, prim_state, prim::P_N);

    const Vec n_i   = rho_i / grid.m_i;
    const Vec n_n   = rho_n / grid.m_n;
    const Vec T_i   = p_i / (2.0 * n_i * grid.k_b);
    const Vec T_n   = p_n / (n_n * grid.k_b);

    const Vec alpha = rho_i % nu_in(grid, n_n, T_i, T_n);
    const Vec w_old = V - U;
    const Vec rho_tot = rho_i + rho_n;
    const Vec V_cm  = (rho_i % V + rho_n % U) / rho_tot;
    const Vec mu    = rho_i % rho_n / rho_tot;
    const Vec lambda_drag = alpha % rho_tot / (rho_i % rho_n);

    const Vec w_new = w_old / (1.0 + dt * lambda_drag);
    const Vec V_new = V_cm + (rho_n / rho_tot) % w_new;
    const Vec U_new = V_cm - (rho_i / rho_tot) % w_new;

    const Vec dE_drag = 0.5 * mu % (w_old % w_old - w_new % w_new);
    const float f_i = grid.m_n / (grid.m_i + grid.m_n);
    const float f_n = grid.m_i / (grid.m_i + grid.m_n);
    const Vec p_i_new = p_i + (2.0f / 3.0f) * f_i * dE_drag;
    const Vec p_n_new = p_n + (2.0f / 3.0f) * f_n * dE_drag;

    prim_state.zeros();
    prim_state += scalar_to(grid, rho_i,   prim::RHO_I);
    prim_state += scalar_to(grid, rho_n,   prim::RHO_N);
    prim_state += scalar_to(grid, V_new,   prim::V);
    prim_state += scalar_to(grid, U_new,   prim::U);
    prim_state += scalar_to(grid, p_i_new, prim::P_I);
    prim_state += scalar_to(grid, p_n_new, prim::P_N);
}

// ----------------------------------------------------------------------------
// Stage C (writeup §3.7 Eqs. (temp-relax)--(temp-reconstruct)): point-implicit
// temperature equilibration between ions/electrons and neutrals. Exact for
// the linear drift ODE with α frozen at the post-drag state.
// ----------------------------------------------------------------------------
static void apply_temperature_stage(const Grid& grid, Vec& prim_state, float dt) {
    Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    Vec V     = get_scalar(grid, prim_state, prim::V);
    Vec U     = get_scalar(grid, prim_state, prim::U);
    Vec p_i   = get_scalar(grid, prim_state, prim::P_I);
    Vec p_n   = get_scalar(grid, prim_state, prim::P_N);

    const Vec n_i   = rho_i / grid.m_i;
    const Vec n_n   = rho_n / grid.m_n;
    const Vec T_i   = p_i / (2.0 * n_i * grid.k_b);
    const Vec T_n   = p_n / (n_n * grid.k_b);

    const Vec alpha = rho_i % nu_in(grid, n_n, T_i, T_n);
    const Vec C_i = 3.0 * n_i * grid.k_b;
    const Vec C_n = 1.5 * n_n * grid.k_b;
    const Vec C_sum = C_i + C_n;

    const Vec beta = 3.0 * grid.k_b * alpha / (grid.m_i + grid.m_n);
    const Vec lambda_T = beta % C_sum / (C_i % C_n);

    const Vec T_cm  = (C_i % T_i + C_n % T_n) / C_sum;
    const Vec dT    = (T_i - T_n) / (1.0 + dt * lambda_T);
    const Vec T_i_new = T_cm + (C_n / C_sum) % dT;
    const Vec T_n_new = T_cm - (C_i / C_sum) % dT;

    const Vec p_i_new = 2.0 * n_i * grid.k_b % T_i_new;
    const Vec p_n_new =       n_n * grid.k_b % T_n_new;

    prim_state.zeros();
    prim_state += scalar_to(grid, rho_i,   prim::RHO_I);
    prim_state += scalar_to(grid, rho_n,   prim::RHO_N);
    prim_state += scalar_to(grid, V,       prim::V);
    prim_state += scalar_to(grid, U,       prim::U);
    prim_state += scalar_to(grid, p_i_new, prim::P_I);
    prim_state += scalar_to(grid, p_n_new, prim::P_N);
}

// ----------------------------------------------------------------------------
// Thomas algorithm (in-place destruction of b, d): solves
//   a[i] x[i-1] + b[i] x[i] + c[i] x[i+1] = d[i],  i = 0..N-1,
// where a[0] and c[N-1] are unused. No pivoting; the matrix is strictly
// diagonally dominant for the backward-Euler diffusion operator below.
// ----------------------------------------------------------------------------
static Vec thomas_solve(const Vec& a, Vec b, const Vec& c, Vec d) {
    const arma::uword N = b.n_elem;
    for (arma::uword i = 1; i < N; ++i) {
        const float m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }
    Vec x(N);
    x[N - 1] = d[N - 1] / b[N - 1];
    for (arma::uword i = N - 1; i > 0; --i) {
        x[i - 1] = (d[i - 1] - c[i - 1] * x[i]) / b[i - 1];
    }
    return x;
}

// ----------------------------------------------------------------------------
// Stage D (writeup §3.7 Eqs. (cond-tridiag)--(cond-coeffs)): backward-Euler
// heat conduction with conductivities lagged at the post-equilibration state.
// One tridiagonal solve per species (charged fluid uses κ_e + κ_i, neutrals
// use κ_n). Ghost temperatures are folded into the RHS as Dirichlet data,
// mirroring the φ_g extrapolation used by rhs_implicit_state.
// ----------------------------------------------------------------------------
static Vec ghost_extended_T(const Grid& grid,
                            const Vec& cons_state,
                            arma::uword cons_RHO,
                            arma::uword cons_MOM,
                            arma::uword cons_E,
                            float k_per_n,   // 2*k_B for ions, k_B for neutrals
                            float m_species,
                            Vec& T_ghost_inner_out,
                            Vec& T_ghost_outer_out)
{
    // T at all ns cells, computed cell-locally from cons + cell-centered phi_g.
    const Vec rho   = get_scalar(grid, cons_state, cons_RHO);
    const Vec mom   = get_scalar(grid, cons_state, cons_MOM);
    const Vec etot  = get_scalar(grid, cons_state, cons_E);
    const Vec vel   = mom / rho;
    const Vec phi_g = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    const Vec p     = 2.0/3.0 * etot - 1.0/3.0 * rho % vel % vel - 2.0/3.0 * rho % phi_g;
    const Vec n     = rho / m_species;
    const Vec T     = p / (k_per_n * n);

    // T at ghost positions (i = -1 and i = ns). Same construction as the
    // ip1/im1 paths in rhs_implicit_state: cons from the ghost cube, phi_g
    // taken at the face (Neumann mirror beyond the boundary face).
    const Vec cons_ip1 = ip1(grid, cons_state);
    const Vec cons_im1 = im1(grid, cons_state);

    auto T_at_shifted = [&](const Vec& cons_shift, bool outer) {
        const Vec rho_s = get_scalar(grid, cons_shift, cons_RHO);
        const Vec mom_s = get_scalar(grid, cons_shift, cons_MOM);
        const Vec e_s   = get_scalar(grid, cons_shift, cons_E);
        const Vec v_s   = mom_s / rho_s;
        Vec phi_g_s     = outer
            ? ip1(grid, 0.5 * (grid.phi_g_iph + grid.phi_g_imh), SLICE)
            : im1(grid, 0.5 * (grid.phi_g_iph + grid.phi_g_imh), SLICE);
        if (outer) phi_g_s[grid.ns - 1] = grid.phi_g_iph[grid.ns - 1];
        else       phi_g_s[0]           = grid.phi_g_imh[0];
        const Vec p_s = 2.0/3.0 * e_s - 1.0/3.0 * rho_s % v_s % v_s - 2.0/3.0 * rho_s % phi_g_s;
        const Vec n_s = rho_s / m_species;
        return Vec(p_s / (k_per_n * n_s));
    };

    Vec T_ip1 = T_at_shifted(cons_ip1, true);
    Vec T_im1 = T_at_shifted(cons_im1, false);

    T_ghost_outer_out = T_ip1.tail(1);  // value at i = ns (only slot used)
    T_ghost_inner_out = T_im1.head(1);  // value at i = -1 (only slot used)
    return T;
}

static void apply_conduction_stage(const Grid& grid, Vec& prim_state, float dt) {
    // Conduction operates on temperatures; rebuild cons once to reuse the
    // existing ip1/im1 ghost machinery, then convert back.
    Vec cons_state = prim2cons(grid, prim_state);

    // Cell-centered densities/temperatures.
    const Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    const Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    const Vec n_i = rho_i / grid.m_i;
    const Vec n_n = rho_n / grid.m_n;

    Vec T_ghost_in_i, T_ghost_out_i, T_ghost_in_n, T_ghost_out_n;
    Vec T_i = ghost_extended_T(grid, cons_state,
                               cons::RHO_I, cons::MOM_I, cons::E_I,
                               2.0f * grid.k_b, grid.m_i,
                               T_ghost_in_i, T_ghost_out_i);
    Vec T_n = ghost_extended_T(grid, cons_state,
                               cons::RHO_N, cons::MOM_N, cons::E_N,
                               grid.k_b, grid.m_n,
                               T_ghost_in_n, T_ghost_out_n);

    // Face conductivities (lagged at current T). Reuse the same averaging as
    // rhs_implicit_state, with ghost values from the shifted cubes.
    const Vec K_e  = kappa_e(n_i, n_n, T_i);
    const Vec K_n  = kappa_n(n_i, n_n, T_i, T_n);

    // Build per-cell K_iph, K_imh by averaging cell-centered K with shifted
    // K (the shifts pull in the ghost K through the ip1/im1 path).
    const Vec cons_ip1 = ip1(grid, cons_state);
    const Vec cons_im1 = im1(grid, cons_state);
    const Vec rho_i_ip1 = get_scalar(grid, cons_ip1, cons::RHO_I);
    const Vec rho_n_ip1 = get_scalar(grid, cons_ip1, cons::RHO_N);
    const Vec rho_i_im1 = get_scalar(grid, cons_im1, cons::RHO_I);
    const Vec rho_n_im1 = get_scalar(grid, cons_im1, cons::RHO_N);
    const Vec n_i_ip1 = rho_i_ip1 / grid.m_i;
    const Vec n_n_ip1 = rho_n_ip1 / grid.m_n;
    const Vec n_i_im1 = rho_i_im1 / grid.m_i;
    const Vec n_n_im1 = rho_n_im1 / grid.m_n;

    // Shifted T at neighbour cells. ghost_extended_T already gave T at ns
    // cells from the unshifted cube; use the shifted cubes here for the K_e
    // and K_n inputs.
    Vec T_ghost_unused_a, T_ghost_unused_b;
    Vec T_i_ip1_full = ghost_extended_T(grid, cons_ip1,
                                        cons::RHO_I, cons::MOM_I, cons::E_I,
                                        2.0f * grid.k_b, grid.m_i,
                                        T_ghost_unused_a, T_ghost_unused_b);
    Vec T_n_ip1_full = ghost_extended_T(grid, cons_ip1,
                                        cons::RHO_N, cons::MOM_N, cons::E_N,
                                        grid.k_b, grid.m_n,
                                        T_ghost_unused_a, T_ghost_unused_b);
    Vec T_i_im1_full = ghost_extended_T(grid, cons_im1,
                                        cons::RHO_I, cons::MOM_I, cons::E_I,
                                        2.0f * grid.k_b, grid.m_i,
                                        T_ghost_unused_a, T_ghost_unused_b);
    Vec T_n_im1_full = ghost_extended_T(grid, cons_im1,
                                        cons::RHO_N, cons::MOM_N, cons::E_N,
                                        grid.k_b, grid.m_n,
                                        T_ghost_unused_a, T_ghost_unused_b);

    const Vec K_e_ip1 = kappa_e(n_i_ip1, n_n_ip1, T_i_ip1_full);
    const Vec K_e_im1 = kappa_e(n_i_im1, n_n_im1, T_i_im1_full);
    const Vec K_n_ip1 = kappa_n(n_i_ip1, n_n_ip1, T_i_ip1_full, T_n_ip1_full);
    const Vec K_n_im1 = kappa_n(n_i_im1, n_n_im1, T_i_im1_full, T_n_im1_full);

    const Vec ds_i_ip1 = ip1(grid, grid.ds_i, SLICE);
    const Vec ds_i_im1 = im1(grid, grid.ds_i, SLICE);
    const Vec ds_iph   = 0.5 * (grid.ds_i + ds_i_ip1);
    const Vec ds_imh   = 0.5 * (ds_i_im1  + grid.ds_i);

    // ----- ion (charged) row ----- K_ion_total = κ_e (κ_i is zero in this code)
    const Vec K_face_iph_i = 0.5 * (K_e + K_e_ip1);
    const Vec K_face_imh_i = 0.5 * (K_e + K_e_im1);
    const Vec C_i = 3.0 * n_i * grid.k_b;
    const Vec g_iph_i = grid.B_i / grid.ds_i % (K_face_iph_i / grid.B_iph) / ds_iph;
    const Vec g_imh_i = grid.B_i / grid.ds_i % (K_face_imh_i / grid.B_imh) / ds_imh;
    Vec R_coef_i = dt * g_iph_i / C_i;
    Vec L_coef_i = dt * g_imh_i / C_i;

    Vec a_i(grid.ns, arma::fill::zeros);
    Vec b_i(grid.ns);
    Vec c_i(grid.ns, arma::fill::zeros);
    Vec d_i = T_i;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        a_i[i] = -L_coef_i[i];
        c_i[i] = -R_coef_i[i];
        b_i[i] = 1.0f + L_coef_i[i] + R_coef_i[i];
    }
    // Inner BC (Dirichlet at i = -1): fold L_0 * T_ghost into RHS, zero a[0].
    d_i[0] += L_coef_i[0] * T_ghost_in_i[0];
    a_i[0] = 0.0f;
    // Outer BC (Dirichlet at i = ns): fold R_{ns-1} * T_ghost into RHS, zero c[ns-1].
    d_i[grid.ns - 1] += R_coef_i[grid.ns - 1] * T_ghost_out_i[0];
    c_i[grid.ns - 1] = 0.0f;
    const Vec T_i_new = thomas_solve(a_i, b_i, c_i, d_i);

    // ----- neutral row -----
    const Vec K_face_iph_n = 0.5 * (K_n + K_n_ip1);
    const Vec K_face_imh_n = 0.5 * (K_n + K_n_im1);
    const Vec C_n_vec = 1.5 * n_n * grid.k_b;
    const Vec g_iph_n = grid.B_i / grid.ds_i % (K_face_iph_n / grid.B_iph) / ds_iph;
    const Vec g_imh_n = grid.B_i / grid.ds_i % (K_face_imh_n / grid.B_imh) / ds_imh;
    Vec R_coef_n = dt * g_iph_n / C_n_vec;
    Vec L_coef_n = dt * g_imh_n / C_n_vec;

    Vec a_n(grid.ns, arma::fill::zeros);
    Vec b_n(grid.ns);
    Vec c_n(grid.ns, arma::fill::zeros);
    Vec d_n = T_n;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        a_n[i] = -L_coef_n[i];
        c_n[i] = -R_coef_n[i];
        b_n[i] = 1.0f + L_coef_n[i] + R_coef_n[i];
    }
    d_n[0] += L_coef_n[0] * T_ghost_in_n[0];
    a_n[0] = 0.0f;
    d_n[grid.ns - 1] += R_coef_n[grid.ns - 1] * T_ghost_out_n[0];
    c_n[grid.ns - 1] = 0.0f;
    const Vec T_n_new = thomas_solve(a_n, b_n, c_n, d_n);

    // Rebuild pressures and store back into prim_state.
    const Vec p_i_new = 2.0 * n_i * grid.k_b % T_i_new;
    const Vec p_n_new =       n_n * grid.k_b % T_n_new;
    const Vec V = get_scalar(grid, prim_state, prim::V);
    const Vec U = get_scalar(grid, prim_state, prim::U);

    prim_state.zeros();
    prim_state += scalar_to(grid, rho_i,   prim::RHO_I);
    prim_state += scalar_to(grid, rho_n,   prim::RHO_N);
    prim_state += scalar_to(grid, V,       prim::V);
    prim_state += scalar_to(grid, U,       prim::U);
    prim_state += scalar_to(grid, p_i_new, prim::P_I);
    prim_state += scalar_to(grid, p_n_new, prim::P_N);
}

// ----------------------------------------------------------------------------
// Stage E (writeup §5.3): point-implicit hydrogen ionization / recombination.
// The local ODE for f ≡ ρ_i / (ρ_i + ρ_n) under quasi-neutrality is
//   df/dt = n_tot [ f(1-f) S_i(T_e) − f² α_r(T_e) ],
// which backward Euler at lagged T_e reduces to the cell-local quadratic
//   A (f^{n+1})² + B f^{n+1} − f^n = 0,
//   A ≡ Δt n_tot (S_i + α_r),   B ≡ 1 − Δt n_tot S_i.
// The positive root gives f^{n+1}; integrated ionizations / recombinations
// then drive the conservative mass / momentum / energy updates of §3.4.
// Mutates V, U, p_i, p_n, ρ_i, ρ_n in place.
// ----------------------------------------------------------------------------
void apply_ionization_stage(const Grid& grid, Vec& prim_state, float dt) {
    const Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    const Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    const Vec V     = get_scalar(grid, prim_state, prim::V);
    const Vec U     = get_scalar(grid, prim_state, prim::U);
    const Vec p_i   = get_scalar(grid, prim_state, prim::P_I);
    const Vec p_n   = get_scalar(grid, prim_state, prim::P_N);

    const float m = grid.m_i;  // m_i = m_n for hydrogen
    const Vec n_i = rho_i / m;
    const Vec n_n = rho_n / m;
    const Vec T_i = p_i / (2.0 * n_i * grid.k_b);   // T_e = T_i under quasi-neutrality
    const Vec T_n = p_n / (n_n * grid.k_b);

    // Rate coefficients at lagged T_e.
    const Vec S = ionization_rate_S(grid, T_i);
    const Vec a = recombination_rate_alpha(grid, T_i);

    // Backward-Euler quadratic A f² + B f - f^old = 0.
    const Vec rho_tot = rho_i + rho_n;
    const Vec n_tot   = rho_tot / m;
    const Vec f_old   = rho_i / rho_tot;
    const Vec Acoef = dt * n_tot % (S + a);
    const Vec Bcoef = 1.0 - dt * n_tot % S;

    Vec f_new(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        // Degenerate A → 0 (no rate activity): linear fallback f = f_old / B.
        if (Acoef(i) < 1e-30f) {
            const float B_i = Bcoef(i);
            f_new(i) = (B_i > 0.0f) ? (f_old(i) / B_i) : f_old(i);
        } else {
            const float B_i = Bcoef(i);
            const float disc = B_i * B_i + 4.0f * Acoef(i) * f_old(i);
            // disc ≥ 0 always when f_old, A ≥ 0; guard against rounding.
            const float sqrt_disc = std::sqrt(std::max(disc, 0.0f));
            f_new(i) = (-B_i + sqrt_disc) / (2.0f * Acoef(i));
        }
        if (f_new(i) < 0.0f) f_new(i) = 0.0f;
        if (f_new(i) > 1.0f) f_new(i) = 1.0f;
    }

    // Integrated ionization / recombination counts (per m³) over Δt. Evaluated
    // at post-step densities and lagged T_e (writeup §5.3, plan §3.4 issue d).
    const Vec rho_i_new = f_new % rho_tot;
    const Vec rho_n_new = (1.0 - f_new) % rho_tot;
    const Vec n_i_new   = rho_i_new / m;
    const Vec n_n_new   = rho_n_new / m;
    const Vec Gamma_ion = dt * n_i_new % n_n_new % S;
    const Vec Gamma_rec = dt * n_i_new % n_i_new % a;

    // Momentum bookkeeping (lagged velocities).
    const Vec rhoV_i_old = rho_i % V;
    const Vec rhoU_n_old = rho_n % U;
    const Vec dmom       = m * (Gamma_ion % U - Gamma_rec % V);
    const Vec rhoV_i_new = rhoV_i_old + dmom;
    const Vec rhoU_n_new = rhoU_n_old - dmom;
    const Vec V_new = rhoV_i_new / rho_i_new;
    const Vec U_new = rhoU_n_new / rho_n_new;

    // Energy bookkeeping (writeup §5.3, eqs. stageE-ei / stageE-en). We update
    // in *thermal* form: convert pressures to thermal energy ε = (3/2) p, apply
    // the source-driven Δε, then rebuild pressures from the new ε and KE.
    const Vec eps_i_old = 1.5 * p_i;
    const Vec eps_n_old = 1.5 * p_n;
    const Vec KE_i_old  = 0.5 * rho_i % V % V;
    const Vec KE_n_old  = 0.5 * rho_n % U % U;
    const Vec KE_i_new  = 0.5 * rho_i_new % V_new % V_new;
    const Vec KE_n_new  = 0.5 * rho_n_new % U_new % U_new;

    const Vec th_inherit = 1.5 * grid.k_b * T_n;   // per-particle thermal carried by new ion
    const Vec th_leave   = 1.5 * grid.k_b * T_i;   // per-particle thermal leaving with recomb. neutral
    const Vec KE_inherit = 0.5 * m * U % U;        // per-particle KE inherited by new ion (from neutral)
    const Vec KE_leave   = 0.5 * m * V % V;        // per-particle KE leaving ion fluid

    // Δ(total energy_i) = Q^{e_i}_ion · Δt  (writeup §2.3 boxed)
    //                   = Γ_ion^Δ (3/2 k_B T_n + ½ m U²)
    //                   − Γ_rec^Δ (3/2 k_B T_i + ½ m V²)
    //                   − Γ_ion^Δ χ_H
    const Vec dE_i = Gamma_ion % (th_inherit + KE_inherit)
                   - Gamma_rec % (th_leave   + KE_leave)
                   - Gamma_ion * grid.chi_H_J;
    const Vec dE_n = -Gamma_ion % (th_inherit + KE_inherit)
                   +  Gamma_rec % (th_leave   + KE_leave);

    // Thermal-energy update: Δε = (Q_ion · Δt) − ΔKE. The boxed Q^{e_i}_ion
    // accounts only for thermal+KE flux between species (Meier & Shumlak 2012);
    // gravity contributes Δρ_s · φ_g to each species' conserved energy through
    // the mass conservation, which is "tracked silently" via ρ_new φ_g in
    // prim2cons() and does *not* enter the thermal update.
    const Vec eps_i_new = eps_i_old + dE_i - (KE_i_new - KE_i_old);
    const Vec eps_n_new = eps_n_old + dE_n - (KE_n_new - KE_n_old);

    // Floor pressures at a small positive value if χ_H drain exceeded the
    // available ion thermal energy. This is a non-conservative regularization
    // (writeup §5.3); the lost energy is the radiative-heating term that a
    // future PR will supply.
    const float p_floor = 1.0e-15f;
    Vec p_i_new = (2.0f / 3.0f) * eps_i_new;
    Vec p_n_new = (2.0f / 3.0f) * eps_n_new;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        if (p_i_new(i) < p_floor) p_i_new(i) = p_floor;
        if (p_n_new(i) < p_floor) p_n_new(i) = p_floor;
    }

    prim_state.zeros();
    prim_state += scalar_to(grid, rho_i_new, prim::RHO_I);
    prim_state += scalar_to(grid, rho_n_new, prim::RHO_N);
    prim_state += scalar_to(grid, V_new,     prim::V);
    prim_state += scalar_to(grid, U_new,     prim::U);
    prim_state += scalar_to(grid, p_i_new,   prim::P_I);
    prim_state += scalar_to(grid, p_n_new,   prim::P_N);
}

// ----------------------------------------------------------------------------
// Backward-Euler integrator (writeup §3.7, §5.3): explicit MUSCL step for R_E,
// then four stiff sub-steps (B→C→D→E) that solve R_I exactly inside one Δt.
// ----------------------------------------------------------------------------
Vec advance_Euler_state(Grid& grid, const Vec& xn_state, const Vec& dt_i) {
    broadcast_dt(grid, dt_i);

    Vec U_star = xn_state + grid.dt_state % rhs_explicit_state(grid, xn_state);
    Vec prim   = cons2prim(grid, U_star);

    const float dt = dt_i(0);  // uniform by cal_dt_i construction
    apply_drag_stage(grid, prim, dt);
    apply_temperature_stage(grid, prim, dt);
    apply_conduction_stage(grid, prim, dt);
    if (grid.enable_ionization) {
        apply_ionization_stage(grid, prim, dt);
    }

    return prim2cons(grid, prim);
}

// Pure-explicit step: identical to the first half of advance_Euler_state with
// the implicit drag / temperature / conduction stages omitted. Useful for
// comparison runs that need R_I ≡ 0.
Vec advance_Euler_explicit_state(Grid& grid, const Vec& xn_state, const Vec& dt_i) {
    broadcast_dt(grid, dt_i);
    return xn_state + grid.dt_state % rhs_explicit_state(grid, xn_state);
}

Vec advance_RK4(Grid& grid, const Vec& xn_state, const Vec& dt_i) {
    broadcast_dt(grid, dt_i);
    const Vec k1 = grid.dt_state % rhs_explicit_state(grid, xn_state);
    const Vec k2 = grid.dt_state % rhs_explicit_state(grid, xn_state + 0.5 * k1);
    const Vec k3 = grid.dt_state % rhs_explicit_state(grid, xn_state + 0.5 * k2);
    const Vec k4 = grid.dt_state % rhs_explicit_state(grid, xn_state + k3);
    return xn_state + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

// ============================================================================
// Debug
// ============================================================================

void print_xn(const Grid& grid, const Vec& xn) {
    for (arma::uword i = 0; i < grid.ns; ++i) {
        for (arma::uword j = 0; j < num_of_eq; ++j) {
            std::cout << xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, j)) << " ";
        }
        std::cout << std::endl;
    }
}

} // namespace chromosphere
