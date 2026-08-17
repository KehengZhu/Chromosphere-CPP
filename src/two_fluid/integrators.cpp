/*!
 * @file two_fluid/integrators.cpp
 * @brief Two-fluid time integrators and the operator-split source stages.
 * @ingroup two_fluid_solver
 *
 * advance_Euler_state() is the semi-implicit driver. Its stages, in order:
 *
 *  - **Stage A** — explicit MUSCL + Rusanov step (rhs_explicit_state()).
 *  - **Stage B** — point-implicit ion-neutral drag and frictional heating.
 *  - **Stage C** — point-implicit temperature equilibration (three-way
 *    nu_ei / nu_en / nu_in relaxation when Grid::enable_Te).
 *  - **Stage D** — tridiagonal field-aligned heat conduction per species.
 *  - **Stage E** — point-implicit hydrogen ionization / recombination.
 *  - **Stage R** — backward-Euler radiative cooling.
 *  - beam, coronal-heating and positivity-floor stages, each gated by its
 *    own Grid flag.
 *
 * NONE of this is release physics: the release timestep is the two stages in
 * single_fluid/integrator.cpp.
 */
#include "two_fluid/two_fluid.hpp"
#include "physics.hpp"
#include "profiling.hpp"
#include "parallel.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace chromosphere {

// ============================================================================
// Time-step calculation
// ============================================================================

Vec cal_max_v_i(const Grid& grid, const Vec& xn_state) {
    Vec spc_r_state = cal_spectral_radius_state(grid, xn_state);
    return get_scalar(grid, spc_r_state, 0);
}

Vec cal_dt_i(const Grid& grid, const Vec& xn_state) {
    if (!grid.eos_gamma_table.empty())
        throw std::logic_error(
            "the legacy two-fluid timestep cannot size a Gamma1-table Grid; "
            "use mixture_timestep");
    ProfileScope timer(ProfileRegion::Cfl);
    Vec dt_i = grid.CFL * grid.ds_i / cal_max_v_i(grid, xn_state);
    float dt_min_i = arma::min(dt_i);
    TimestepLimiter active_limiter = TimestepLimiter::Acoustic;

    // Beam-heating timescale limit. In the low-density chromosphere (n ~ 10^16
    // m^-3) the heat capacity is tiny, so the flare beam's heating time
    // τ_heat = ε / Q can be SHORTER than the hydro CFL step. Operator splitting
    // then over-heats a cell in one step (ΔT ~ 10^5 K) and breaks. Cap dt so the
    // beam raises the total thermal energy by at most BEAM_HEAT_CFL per step.
    if (grid.enable_beam_heating) {
        const Vec prim  = cons2prim(grid, xn_state);
        const Vec rho_i = get_scalar(grid, prim, prim::RHO_I);
        const Vec rho_n = get_scalar(grid, prim, prim::RHO_N);
        const Vec p_i   = get_scalar(grid, prim, prim::P_I);
        const Vec p_n   = get_scalar(grid, prim, prim::P_N);
        const Vec p_e   = get_scalar(grid, prim, prim::P_E);
        const Vec n_i   = rho_i / grid.m_i;
        const Vec n_n   = rho_n / grid.m_n;
        const Vec Q     = beam_heating_rate(grid, n_i, n_n);   // W/m^3
        const float Qmax = Q.max();
        if (Qmax > 0.0f) {
            const float BEAM_HEAT_CFL = 0.1f;
            // Three-T: the beam deposits into electrons alone, so the heating
            // timescale that matters is ε_e/Q (a much smaller pool). Single-T:
            // the deposited heat is shared, so use the total thermal energy.
            const Vec eps = grid.enable_Te ? Vec(grid.inv_gm1() * p_e)
                                           : Vec(grid.inv_gm1() * (p_i + p_n));
            // dt ≤ BEAM_HEAT_CFL · min(ε/Q) over cells receiving the beam.
            Vec tau = eps / arma::clamp(Q, 1.0e-30f, arma::datum::inf);
            // Only constrain where the beam actually deposits (Q>0).
            for (arma::uword i = 0; i < grid.ns; ++i)
                if (Q(i) <= 0.0f) tau(i) = arma::datum::inf;
            const float dt_beam = BEAM_HEAT_CFL * arma::min(tau);
            if (dt_beam < dt_min_i) {
                dt_min_i = dt_beam;
                active_limiter = TimestepLimiter::BeamHeating;
            }
        }
    }

    // Coronal-heating timescale limit (same rationale as the beam guard). The
    // steady footpoint heating is gentle (τ_heat = ε/H ~ tens of s in the tenuous
    // corona, far longer in the dense chromosphere), so this rarely binds — but a
    // ramped/localized enhancement (Phase 3) can shorten it, so cap dt to a small
    // fractional energy gain per step where H > 0, exactly as for the beam.
    if (grid.enable_coronal_heating) {
        const Vec prim = cons2prim(grid, xn_state);
        const Vec p_i  = get_scalar(grid, prim, prim::P_I);
        const Vec p_n  = get_scalar(grid, prim, prim::P_N);
        const Vec p_e  = get_scalar(grid, prim, prim::P_E);
        const Vec H    = coronal_heating_rate(grid);   // W/m^3
        const float Hmax = H.max();
        if (Hmax > 0.0f) {
            const float HEAT_CFL = 0.1f;
            const Vec eps = grid.enable_Te ? Vec(grid.inv_gm1() * p_e)
                                           : Vec(grid.inv_gm1() * (p_i + p_n));
            Vec tau = eps / arma::clamp(H, 1.0e-30f, arma::datum::inf);
            for (arma::uword i = 0; i < grid.ns; ++i)
                if (H(i) <= 0.0f) tau(i) = arma::datum::inf;
            const float dt_heat = HEAT_CFL * arma::min(tau);
            if (dt_heat < dt_min_i) {
                dt_min_i = dt_heat;
                active_limiter = TimestepLimiter::CoronalHeating;
            }
        }
    }
    profile_note_timestep_limiter(active_limiter);
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
// TRAC adaptive cutoff temperature (Johnston et al. 2020 Eq. 8). T_c is the
// highest charged-fluid temperature among cells where the TR is under-resolved
// (L_R/L_T > δ = 1/2, with L_T = T/|dT/ds| the temperature length scale and
// L_R = Δs the grid spacing), clamped to [trac_T_chrom, 0.2 T_peak]. When no
// cell is under-resolved it defaults to the floor trac_T_chrom.
// ----------------------------------------------------------------------------
float compute_trac_cutoff_T(const Grid& grid, const Vec& prim_state) {
    const Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    const Vec p_i   = get_scalar(grid, prim_state, prim::P_I);
    const Vec p_e   = get_scalar(grid, prim_state, prim::P_E);
    const Vec n_i   = rho_i / grid.m_i;
    // TRAC broadens the (electron) heat-conduction operator, so it keys off the
    // temperature that drives it: T_e in the three-T model, T_charged otherwise.
    const Vec T     = grid.enable_Te ? Vec(p_e / (n_i * grid.k_b))
                                     : Vec(p_i / (2.0f * n_i * grid.k_b));

    const float Tpeak    = arma::max(T);
    const float Tc_upper = std::max(grid.trac_Tc_max_frac * Tpeak, grid.trac_T_chrom);
    float Tc = grid.trac_T_chrom;
    const arma::uword ns = grid.ns;
    for (arma::uword i = 0; i < ns; ++i) {
        float dTds;
        // Center-to-center separations (grid.ds_iph_i/ds_imh_i mirror the boundary
        // widths, so the one-sided ends match the legacy Δs on a uniform mesh).
        if (i == 0)            dTds = (T(1) - T(0)) / grid.ds_iph_i(0);
        else if (i == ns - 1)  dTds = (T(i) - T(i - 1)) / grid.ds_imh_i(i);
        else                   dTds = (T(i + 1) - T(i - 1)) / (grid.ds_iph_i(i) + grid.ds_imh_i(i));
        // Under-resolved ⟺ L_R/L_T = Δs|dT/ds|/T > 1/2.
        if (std::fabs(dTds) * grid.ds_i(i) > 0.5f * T(i)) {
            if (T(i) > Tc) Tc = T(i);
        }
    }
    Tc = std::min(Tc, Tc_upper);
    Tc = std::max(Tc, grid.trac_T_chrom);

    // Cutoff-temperature limiter (Johnston 2020 Appendix A.2): a sudden jump in
    // T_c changes κ/Λ over many cells at once, producing a conduction/radiation
    // shock. Limit the per-step change to a small ratio (both directions) so the
    // broadening region grows/shrinks gradually; the integrated quantities are
    // unaffected. grid.trac_cutoff_T holds the previous step's value.
    const float prev      = grid.trac_cutoff_T;
    const float max_ratio = 1.03f;
    if (prev >= grid.trac_T_chrom) {
        Tc = std::min(Tc, prev * max_ratio);
        Tc = std::max(Tc, prev / max_ratio);
        Tc = std::min(Tc, Tc_upper);
        Tc = std::max(Tc, grid.trac_T_chrom);
    }
    return Tc;
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
    const Vec p_e = get_scalar(grid, prim_state, prim::P_E);   // carried untouched

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

    // Single-fluid limit: rigidly lock the drift to zero (V→U→V_cm), the exact
    // dt·λ_drag→∞ limit. The full relative kinetic energy 0.5·μ·w_old² is then
    // thermalized below via dE_drag, conserving total momentum and energy.
    const Vec w_new = grid.single_fluid
                          ? Vec(arma::zeros<Vec>(w_old.n_elem))
                          : Vec(w_old / (1.0 + dt * lambda_drag));
    const Vec V_new = V_cm + (rho_n / rho_tot) % w_new;
    const Vec U_new = V_cm - (rho_i / rho_tot) % w_new;

    const Vec dE_drag = 0.5 * mu % (w_old % w_old - w_new % w_new);
    const float f_i = grid.m_n / (grid.m_i + grid.m_n);
    const float f_n = grid.m_i / (grid.m_i + grid.m_n);
    const Vec p_i_new = p_i + grid.gm1() * f_i * dE_drag;
    const Vec p_n_new = p_n + grid.gm1() * f_n * dE_drag;

    // Frictional heating goes to the combined charged pool (P_I = total charged
    // pressure); with the electron pressure P_E carried unchanged through this
    // stage, the heat lands on the protons (p_proton = P_I − P_E) — ion–neutral
    // friction does not heat electrons. P_E is preserved verbatim.
    prim_state.zeros();
    prim_state += scalar_to(grid, rho_i,   prim::RHO_I);
    prim_state += scalar_to(grid, rho_n,   prim::RHO_N);
    prim_state += scalar_to(grid, V_new,   prim::V);
    prim_state += scalar_to(grid, U_new,   prim::U);
    prim_state += scalar_to(grid, p_i_new, prim::P_I);
    prim_state += scalar_to(grid, p_n_new, prim::P_N);
    prim_state += scalar_to(grid, p_e,     prim::P_E);
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
    Vec p_e   = get_scalar(grid, prim_state, prim::P_E);

    const Vec n_i   = rho_i / grid.m_i;
    const Vec n_n   = rho_n / grid.m_n;

    if (!grid.enable_Te) {
        // ---- Single-temperature baseline (unchanged 2-way T_charged ↔ T_n) ----
        const Vec T_i   = p_i / (2.0 * n_i * grid.k_b);   // T_charged
        const Vec T_n   = p_n / (n_n * grid.k_b);

        const Vec alpha = rho_i % nu_in(grid, n_n, T_i, T_n);
        // Heat capacities per volume C = nk/(γ−1); the charged row carries ions
        // + electrons (n_e = n_i ⇒ 2 n_i), the neutral row n_n. (β below is a
        // kinetic collisional coefficient — see Grid::gm1 note — left at 3 k_B.)
        const Vec C_i = 2.0f * grid.inv_gm1() * n_i * grid.k_b;
        const Vec C_n = grid.inv_gm1() * n_n * grid.k_b;
        const Vec C_sum = C_i + C_n;

        const Vec beta = 3.0 * grid.k_b * alpha / (grid.m_i + grid.m_n);
        const Vec lambda_T = beta % C_sum / (C_i % C_n);

        const Vec T_cm  = (C_i % T_i + C_n % T_n) / C_sum;
        // Single-fluid limit: lock both temperatures to the common T_cm (dt·λ_T→∞).
        const Vec dT    = grid.single_fluid
                              ? Vec(arma::zeros<Vec>(T_i.n_elem))
                              : Vec((T_i - T_n) / (1.0 + dt * lambda_T));
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
        prim_state += scalar_to(grid, p_e,     prim::P_E);
        return;
    }

    // ---- Three-temperature point-implicit relaxation (T_e, T_i, T_n) --------
    // Generalizes the 2-way solve to a symmetric 3×3 backward-Euler relaxation
    // with three pairwise conductances g_ab [W m⁻³ K⁻¹] (g_ab = g_ba ⇒ energy
    // conserved exactly): electron–ion g_ei = C_e ν_ei, electron–neutral
    // g_en = C_e ν_en, ion–neutral g_in = 3 k_B α/(m_i+m_n) (the same β as the
    // 2-way solve, now connecting protons↔neutrals only — electrons reach the
    // neutrals through the separate ν_en channel). Heat capacities split the
    // old combined charged C = 3 n_i k_B into equal electron and proton halves.
    const Vec p_proton = p_i - p_e;                       // derived proton pressure
    const Vec T_e = p_e      / (n_i * grid.k_b);
    const Vec T_i = p_proton / (n_i * grid.k_b);
    const Vec T_n = p_n      / (n_n * grid.k_b);

    const Vec C_e = grid.inv_gm1() * n_i * grid.k_b;
    const Vec C_i = grid.inv_gm1() * n_i * grid.k_b;
    const Vec C_n = grid.inv_gm1() * n_n * grid.k_b;

    const Vec alpha = rho_i % nu_in(grid, n_n, T_i, T_n);
    const Vec g_in  = 3.0 * grid.k_b * alpha / (grid.m_i + grid.m_n);
    const Vec g_ei  = C_e % nu_ei(grid, n_i, T_e);
    const Vec g_en  = C_e % nu_en(grid, n_n, T_e);

    Vec T_e_new(grid.ns), T_i_new(grid.ns), T_n_new(grid.ns);
    for (arma::uword k = 0; k < grid.ns; ++k) {
        // Symmetric 3×3 system M·[T_e,T_i,T_n]ᵀ = [C_e T_e, C_i T_i, C_n T_n]ᵀ.
        //   M_ee = C_e + dt(g_ei+g_en),  M_ei = -dt g_ei,  M_en = -dt g_en
        //   M_ii = C_i + dt(g_ei+g_in),  M_in = -dt g_in
        //   M_nn = C_n + dt(g_en+g_in)
        const float ge = g_ei(k), gn = g_en(k), gi = g_in(k);
        const float a11 = C_e(k) + dt*(ge+gn), a12 = -dt*ge, a13 = -dt*gn;
        const float a22 = C_i(k) + dt*(ge+gi), a23 = -dt*gi;
        const float a33 = C_n(k) + dt*(gn+gi);
        const float b1 = C_e(k)*T_e(k), b2 = C_i(k)*T_i(k), b3 = C_n(k)*T_n(k);
        // Solve via cofactor expansion (symmetric, diagonally dominant ⇒ regular).
        const float c11 =  a22*a33 - a23*a23;
        const float c12 = -(a12*a33 - a23*a13);
        const float c13 =  a12*a23 - a22*a13;
        const float c22 =  a11*a33 - a13*a13;
        const float c23 = -(a11*a23 - a12*a13);
        const float c33 =  a11*a22 - a12*a12;
        const float det = a11*c11 + a12*c12 + a13*c13;
        const float inv = (std::fabs(det) > 0.0f) ? 1.0f/det : 0.0f;
        T_e_new(k) = inv*(c11*b1 + c12*b2 + c13*b3);
        T_i_new(k) = inv*(c12*b1 + c22*b2 + c23*b3);
        T_n_new(k) = inv*(c13*b1 + c23*b2 + c33*b3);
    }

    const Vec p_e_new      = n_i * grid.k_b % T_e_new;
    const Vec p_proton_new = n_i * grid.k_b % T_i_new;
    const Vec p_n_new      = n_n * grid.k_b % T_n_new;
    const Vec p_i_new      = p_proton_new + p_e_new;      // total charged pressure

    prim_state.zeros();
    prim_state += scalar_to(grid, rho_i,    prim::RHO_I);
    prim_state += scalar_to(grid, rho_n,    prim::RHO_N);
    prim_state += scalar_to(grid, V,        prim::V);
    prim_state += scalar_to(grid, U,        prim::U);
    prim_state += scalar_to(grid, p_i_new,  prim::P_I);
    prim_state += scalar_to(grid, p_n_new,  prim::P_N);
    prim_state += scalar_to(grid, p_e_new,  prim::P_E);
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
    const Vec p     = grid.gm1() * etot - grid.half_gm1() * rho % vel % vel - grid.gm1() * rho % phi_g;
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
        const Vec p_s = grid.gm1() * e_s - grid.half_gm1() * rho_s % v_s % v_s - grid.gm1() * rho_s % phi_g_s;
        const Vec n_s = rho_s / m_species;
        return Vec(p_s / (k_per_n * n_s));
    };

    Vec T_ip1 = T_at_shifted(cons_ip1, true);
    Vec T_im1 = T_at_shifted(cons_im1, false);

    T_ghost_outer_out = T_ip1.tail(1);  // value at i = ns (only slot used)
    T_ghost_inner_out = T_im1.head(1);  // value at i = -1 (only slot used)
    return T;
}

// Electron-temperature variant of ghost_extended_T for the three-temperature
// conduction path. The electron internal energy E_E carries NO kinetic or
// gravitational part, so T_e = (2/3 E_E)/(n_i k_B) = p_e/(n_i k_B) directly —
// no KE/φ subtraction. Ghost values come from the shifted cubes' E_E and ρ_i.
static Vec ghost_extended_Te(const Grid& grid, const Vec& cons_state,
                             Vec& T_ghost_inner_out, Vec& T_ghost_outer_out) {
    auto Te_of = [&](const Vec& cs) {
        const Vec rho_i = get_scalar(grid, cs, cons::RHO_I);
        const Vec e_e   = get_scalar(grid, cs, cons::E_E);
        const Vec n_i   = rho_i / grid.m_i;
        return Vec(grid.gm1() * e_e / (n_i * grid.k_b));
    };
    const Vec T     = Te_of(cons_state);
    const Vec T_ip1 = Te_of(ip1(grid, cons_state));
    const Vec T_im1 = Te_of(im1(grid, cons_state));
    T_ghost_outer_out = T_ip1.tail(1);
    T_ghost_inner_out = T_im1.head(1);
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

    // The "charged row" conducts either the combined charged temperature
    // T_charged = p_total/(2 n_i k_B) with heat capacity 3 n_i k_B (single-T
    // baseline) or the electron temperature T_e = p_e/(n_i k_B) with heat
    // capacity 1.5 n_i k_B (three-temperature model). Spitzer κ_e ∝ T^{5/2} is
    // evaluated at whichever temperature drives the row; the neutral row is
    // identical in both. The proton conduction κ_i ≈ κ_e/43 is dropped (plan §3).
    const bool  Te_on = grid.enable_Te;
    // C_charged = cfac · n_i k_B. Single-T conducts the combined charged pool
    // (ions+electrons ⇒ 2 n_i): cfac = 2/(γ−1). Three-T conducts electrons only
    // (n_i): cfac = 1/(γ−1). (γ=5/3 ⇒ 3.0 and 1.5 respectively.)
    const float cfac  = Te_on ? grid.inv_gm1() : 2.0f * grid.inv_gm1();
    const Vec   p_i_orig = get_scalar(grid, prim_state, prim::P_I);
    const Vec   p_e_orig = get_scalar(grid, prim_state, prim::P_E);

    Vec T_ghost_in_i, T_ghost_out_i, T_ghost_in_n, T_ghost_out_n;
    Vec T_i = Te_on
        ? ghost_extended_Te(grid, cons_state, T_ghost_in_i, T_ghost_out_i)
        : ghost_extended_T(grid, cons_state, cons::RHO_I, cons::MOM_I, cons::E_I,
                           2.0f * grid.k_b, grid.m_i, T_ghost_in_i, T_ghost_out_i);
    Vec T_n = ghost_extended_T(grid, cons_state,
                               cons::RHO_N, cons::MOM_N, cons::E_N,
                               grid.k_b, grid.m_n,
                               T_ghost_in_n, T_ghost_out_n);

    // Conduction-only outer thermal wall (chromosphere.hpp
    // ::outer_conduction_temperature_override). Replace the hydro outer ghost's
    // temperature on BOTH conduction rows so the Dirichlet datum folded into the
    // RHS below — and the outer ghost conductivity assembled from the shifted
    // cubes further down — stay pinned to the wall while the hydro ghost's own
    // temperature is free to float. Ghost DENSITIES are untouched.
    if (grid.outer_conduction_temperature_override) {
        T_ghost_out_i[0] = grid.outer_conduction_temperature;
        T_ghost_out_n[0] = grid.outer_conduction_temperature;
    }

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
    Vec T_i_ip1_full = Te_on
        ? ghost_extended_Te(grid, cons_ip1, T_ghost_unused_a, T_ghost_unused_b)
        : ghost_extended_T(grid, cons_ip1, cons::RHO_I, cons::MOM_I, cons::E_I,
                           2.0f * grid.k_b, grid.m_i, T_ghost_unused_a, T_ghost_unused_b);
    Vec T_n_ip1_full = ghost_extended_T(grid, cons_ip1,
                                        cons::RHO_N, cons::MOM_N, cons::E_N,
                                        grid.k_b, grid.m_n,
                                        T_ghost_unused_a, T_ghost_unused_b);
    Vec T_i_im1_full = Te_on
        ? ghost_extended_Te(grid, cons_im1, T_ghost_unused_a, T_ghost_unused_b)
        : ghost_extended_T(grid, cons_im1, cons::RHO_I, cons::MOM_I, cons::E_I,
                           2.0f * grid.k_b, grid.m_i, T_ghost_unused_a, T_ghost_unused_b);
    Vec T_n_im1_full = ghost_extended_T(grid, cons_im1,
                                        cons::RHO_N, cons::MOM_N, cons::E_N,
                                        grid.k_b, grid.m_n,
                                        T_ghost_unused_a, T_ghost_unused_b);

    // Last entry of the ip1-shifted temperatures IS the outer ghost; pin it to the
    // conduction wall too so κ at the outer face matches the Dirichlet datum.
    if (grid.outer_conduction_temperature_override) {
        T_i_ip1_full[grid.ns - 1] = grid.outer_conduction_temperature;
        T_n_ip1_full[grid.ns - 1] = grid.outer_conduction_temperature;
    }

    const Vec K_e_ip1 = kappa_e(n_i_ip1, n_n_ip1, T_i_ip1_full);
    const Vec K_e_im1 = kappa_e(n_i_im1, n_n_im1, T_i_im1_full);
    const Vec K_n_ip1 = kappa_n(n_i_ip1, n_n_ip1, T_i_ip1_full, T_n_ip1_full);
    const Vec K_n_im1 = kappa_n(n_i_im1, n_n_im1, T_i_im1_full, T_n_im1_full);

    const Vec ds_i_ip1 = ip1(grid, grid.ds_i, SLICE);
    const Vec ds_i_im1 = im1(grid, grid.ds_i, SLICE);
    const Vec ds_iph   = 0.5 * (grid.ds_i + ds_i_ip1);
    const Vec ds_imh   = 0.5 * (ds_i_im1  + grid.ds_i);

    // Isotropic numerical diffusion (Pandey 2024 §3.3): add K_num = χ·C to each
    // face conductivity, where C is the face-averaged heat capacity per volume
    // (ion C = 3 n_i k_B, neutral C = 1.5 n_n k_B). Folding it into K turns the
    // existing backward-Euler tridiagonal into an unconditionally-stable solve
    // for the combined Spitzer + isotropic-diffusion operator at no extra cost.
    Vec chi_iph(grid.ns), chi_imh(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        chi_iph(i) = numerical_diffusivity_at_face(grid, grid.ds_iph_i(i));
        chi_imh(i) = numerical_diffusivity_at_face(grid, grid.ds_imh_i(i));
    }
    // K_num = χ · C_face. Charged-row capacity is cfac·n_i k_B (cfac = 3 single-T,
    // 1.5 three-T), so the face value is 0.5·cfac·k_B·(n_i + n_i_shift).
    const Vec Knum_iph_i = (0.5f * cfac * grid.k_b) * (n_i + n_i_ip1) % chi_iph;
    const Vec Knum_imh_i = (0.5f * cfac * grid.k_b) * (n_i + n_i_im1) % chi_imh;
    const Vec Knum_iph_n = (0.5f * grid.inv_gm1() * grid.k_b) * (n_n + n_n_ip1) % chi_iph;
    const Vec Knum_imh_n = (0.5f * grid.inv_gm1() * grid.k_b) * (n_n + n_n_im1) % chi_imh;

    // TRAC broadening (Johnston 2020): below the adaptive cutoff T_c the charged-
    // fluid conductivity is enhanced by ε=(T_c/T)^{5/2} (→ κ' = κ(T_c), constant),
    // spreading the unresolved TR onto the coarse grid. Applied per cell at the
    // cell temperature (and at the shifted-neighbour temperatures) before face
    // averaging. Only κ_e is broadened; the neutral row is chromospheric.
    Vec Kb = K_e, Kb_ip1 = K_e_ip1, Kb_im1 = K_e_im1;
    if (grid.enable_trac) {
        Kb     %= trac_broadening_factor(grid, T_i);
        Kb_ip1 %= trac_broadening_factor(grid, T_i_ip1_full);
        Kb_im1 %= trac_broadening_factor(grid, T_i_im1_full);
    }

    // ----- ion (charged) row ----- K_ion_total = κ_e (κ_i is zero in this code)
    // Spitzer κ_e is combined across a face by the legacy arithmetic average on a
    // uniform mesh, or the width-weighted series-resistance combination on a
    // refined mesh (physics.hpp::face_conductivity_series). The isotropic
    // numerical diffusivity K_num keeps its face average either way.
    const Vec K_face_iph_i = (grid.uniform_mesh ? Vec(0.5 * (Kb + Kb_ip1))
                                : face_conductivity_series(Kb, Kb_ip1, grid.ds_i, ds_i_ip1))
                             + Knum_iph_i;
    const Vec K_face_imh_i = (grid.uniform_mesh ? Vec(0.5 * (Kb + Kb_im1))
                                : face_conductivity_series(Kb, Kb_im1, grid.ds_i, ds_i_im1))
                             + Knum_imh_i;
    const Vec C_i = cfac * n_i * grid.k_b;
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
    // Inner BC. Default Dirichlet at i = -1: fold L_0 * T_ghost into RHS, zero
    // a[0]. Neumann (grid.inner_conduction_neumann): zero conductive flux through
    // the base — drop the inner-face coupling entirely (no L_0 ghost term, no L_0
    // in the diagonal), so dT/ds = 0 at the base and the photosphere is insulating.
    a_i[0] = 0.0f;
    if (grid.inner_conduction_neumann) {
        b_i[0] = 1.0f + R_coef_i[0];
    } else {
        d_i[0] += L_coef_i[0] * T_ghost_in_i[0];
    }
    // Outer BC. Either:
    //   * Neumann fixed-flux (impose_outer_heat_flux): the top face carries the
    //     imposed coronal conductive flux q = outer_heat_flux [W/m²] into the
    //     cell. Drop the R_coef coupling (no Dirichlet ghost) and add the flux
    //     as a source ΔT = dt·q/(C·ds). Robust on the coarse grid — does not
    //     over-conduct through a hot T-wall. (RTV q(T); plan upper-BC step 2.)
    //   * Dirichlet (default): fold R_{ns-1}·T_ghost into the RHS. T_ghost_out is
    //     the hydro outer ghost's temperature, or the pinned conduction wall when
    //     outer_conduction_temperature_override is set (substituted above).
    const arma::uword nlast = grid.ns - 1;
    if (grid.impose_outer_heat_flux) {
        b_i[nlast] = 1.0f + L_coef_i[nlast];
        d_i[nlast] = T_i[nlast]
                   + dt * grid.outer_heat_flux / (C_i[nlast] * grid.ds_i[nlast]);
    } else {
        d_i[nlast] += R_coef_i[nlast] * T_ghost_out_i[0];
    }
    c_i[nlast] = 0.0f;
    const Vec T_i_new = thomas_solve(a_i, b_i, c_i, d_i);

    // ----- neutral row ----- (same series-resistance treatment on a refined mesh)
    const Vec K_face_iph_n = (grid.uniform_mesh ? Vec(0.5 * (K_n + K_n_ip1))
                                : face_conductivity_series(K_n, K_n_ip1, grid.ds_i, ds_i_ip1))
                             + Knum_iph_n;
    const Vec K_face_imh_n = (grid.uniform_mesh ? Vec(0.5 * (K_n + K_n_im1))
                                : face_conductivity_series(K_n, K_n_im1, grid.ds_i, ds_i_im1))
                             + Knum_imh_n;
    const Vec C_n_vec = grid.inv_gm1() * n_n * grid.k_b;
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
    a_n[0] = 0.0f;
    if (grid.inner_conduction_neumann) {
        b_n[0] = 1.0f + R_coef_n[0];
    } else {
        d_n[0] += L_coef_n[0] * T_ghost_in_n[0];
    }
    d_n[grid.ns - 1] += R_coef_n[grid.ns - 1] * T_ghost_out_n[0];
    c_n[grid.ns - 1] = 0.0f;
    Vec T_n_new = thomas_solve(a_n, b_n, c_n, d_n);
    Vec T_i_new2 = T_i_new;

    // Flare robustness: at the explosive evaporation front the conduction
    // tridiagonal can occasionally yield a non-finite temperature (extreme
    // κ_e∝T^{5/2} contrasts across a near-discontinuity on the coarse grid). Fall
    // back to the pre-conduction temperature for any such cell rather than
    // letting the NaN propagate. Gated on beam heating → steady runs unchanged.
    if (grid.enable_beam_heating) {
        for (arma::uword i = 0; i < grid.ns; ++i) {
            if (!std::isfinite(T_i_new2(i))) T_i_new2(i) = T_i(i);
            if (!std::isfinite(T_n_new(i)))  T_n_new(i)  = T_n(i);
        }
    }

    // Rebuild pressures and store back into prim_state.
    const Vec p_n_new = n_n * grid.k_b % T_n_new;
    Vec p_i_new, p_e_new;
    if (!Te_on) {
        // Single-T baseline: T_i_new2 is the combined charged temperature.
        p_i_new = 2.0 * n_i * grid.k_b % T_i_new2;        // total charged pressure
        p_e_new = p_e_orig;                               // carried (re-slaved later)
    } else {
        // Three-T: T_i_new2 is the electron temperature. Only electrons conducted,
        // so the total charged pressure shifts by exactly the electron pressure
        // change (protons unchanged); E_I stays the energy-conserving total.
        p_e_new = n_i * grid.k_b % T_i_new2;
        p_i_new = p_i_orig + (p_e_new - p_e_orig);
    }
    const Vec V = get_scalar(grid, prim_state, prim::V);
    const Vec U = get_scalar(grid, prim_state, prim::U);

    prim_state.zeros();
    prim_state += scalar_to(grid, rho_i,   prim::RHO_I);
    prim_state += scalar_to(grid, rho_n,   prim::RHO_N);
    prim_state += scalar_to(grid, V,       prim::V);
    prim_state += scalar_to(grid, U,       prim::U);
    prim_state += scalar_to(grid, p_i_new, prim::P_I);
    prim_state += scalar_to(grid, p_n_new, prim::P_N);
    prim_state += scalar_to(grid, p_e_new, prim::P_E);
}

// ----------------------------------------------------------------------------
// Stage R (writeup §2.4, §5.6): optically-thick chromospheric radiative
// cooling. Carlsson & Leenaarts 2012 recipe summing H I + Ca II + Mg II
// losses; tables in physics.hpp::cl2012. Operates as a backward-Euler sink
// on the folded ion+electron thermal pressure: p_i ← p_i / (1 + Δt Q/ε_i)
// at lagged Q evaluated from the post-conduction T_e. This preserves
// positivity unconditionally — when Δt Q ≫ ε_i, the cell relaxes asymptot-
// ically to a small positive pressure rather than overshooting to a floor.
// Mass / velocities are unaffected; only p_i is modified.
// ----------------------------------------------------------------------------
void apply_radiative_cooling_stage(const Grid& grid, Vec& prim_state, float dt) {
    const Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    const Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    const Vec p_i   = get_scalar(grid, prim_state, prim::P_I);
    const Vec p_e   = get_scalar(grid, prim_state, prim::P_E);

    const Vec n_i = rho_i / grid.m_i;
    const Vec n_n = rho_n / grid.m_n;
    // Radiation is an electron-pool loss (Λ ∝ n_e and the bound-level excitation
    // is electron-collisional). Single-T baseline: cool the combined charged pool
    // at T_charged = p_total/(2 n_i k_B). Three-T: cool the electron pool at the
    // true T_e = p_e/(n_i k_B), and mirror the loss into the total energy P_I.
    const Vec T_e = grid.enable_Te ? Vec(p_e / (n_i * grid.k_b))
                                   : Vec(p_i / (2.0f * n_i * grid.k_b));

    // Total radiative sink: optically-thick chromospheric (CL2012, H/Ca/Mg) plus
    // optically-thin TR/coronal losses n_e n_H Λ(T). The two are stitched by a
    // ~2e4 K switch inside radiative_loss_thin so neither double-counts; together
    // they span the chromosphere through the lower TR / corona, giving the TR the
    // radiative sink to re-radiate the imposed coronal conductive flux q(T).
    // Optically-thin loss is TRAC-broadened (divided by ε): Λ' = Λ(T/T_c)^{5/2}
    // below T_c, so κ'Λ' = κΛ and the TR-integrated radiation is conserved on the
    // coarse grid (Johnston 2020). The optically-thick CL2012 chromospheric term
    // is left unbroadened (it lives below the TRAC region).
    Vec Q_thin = radiative_loss_thin(grid, n_i, n_n, T_e);
    if (grid.enable_trac) Q_thin /= trac_broadening_factor(grid, T_e);
    const Vec Q   = radiative_loss_thick(grid, n_i, n_n, T_e)    // W/m^3, +ve = cooling
                  + Q_thin;
    // Pool that radiates: combined charged (single-T) or electron-only (three-T).
    const Vec p_pool = grid.enable_Te ? p_e : p_i;
    const Vec eps = grid.inv_gm1() * p_pool;                    // thermal energy density
    // Backward-Euler relaxation: ε^{n+1} = ε^n / (1 + Δt Q/ε^n) with Q frozen.
    const Vec denom = 1.0f + dt * Q / arma::clamp(eps, 1.0e-30f, arma::datum::inf);
    const Vec p_pool_new = p_pool / arma::clamp(denom, 1.0e-6f, arma::datum::inf);
    const Vec dp = p_pool_new - p_pool;
    if (grid.enable_Te) {
        prim_state += scalar_to(grid, dp, prim::P_E);   // electrons radiate
        prim_state += scalar_to(grid, dp, prim::P_I);   // mirror into total energy
    } else {
        prim_state += scalar_to(grid, dp, prim::P_I);
    }
}

/// Vacuum-floor / single-fluid-collapse stage for the flare scenario. Explosive
/// evaporation fully ionizes the heated gas, so the neutral density n_n collapses
/// toward zero in the evaporated column. There the neutral velocity decode
/// U = ρ_n U / ρ_n becomes 0/0 noise and the spectral radius |U| + c_s,n blows up
/// → NaN (the ρ_n→0 ill-conditioning of two-fluid codes). We therefore floor both
/// densities, floor the pressures positive, cap the per-fluid speeds at a generous
/// physical ceiling, and — where the gas is ionized — COLLAPSE the neutral onto
/// the charged fluid (U = V, T_n = T_e), the exact single-fluid limit (Gómez
/// Míguez et al. 2024; the same slaving the C7 outer BC uses).
///
/// The collapse triggers on the ionization fraction f = n_i/(n_i+n_n) exceeding
/// F_SLAVE, NOT on an absolute n_n floor. This is the key to a smooth neutral
/// profile: f varies smoothly through the hot column, so the whole ionized region
/// is slaved consistently. An absolute-floor trigger instead let neighbouring
/// cells flip in and out of the slaved state every step — the explicit hydro
/// re-noises the near-zero-mass neutral momentum each step, pushing ρ_n a hair
/// above/below the floor — producing the cell-to-cell checkerboard (sawtooth) seen
/// in u and p_n. Gated on enable_beam_heating so steady scenarios and the existing
/// test suite are untouched. Mutates prim in place.
///
/// Internal to the historical two-fluid integrator: it is called from the
/// operator-split step (guarded by Grid::floors_active()) and has no declaration
/// in two_fluid.hpp.
void apply_flare_floor_stage(const Grid& grid, Vec& prim_state) {
    const float m_i = grid.m_i, m_n = grid.m_n, k_b = grid.k_b;
    const float RHO_I_FLOOR = m_i * 1.0e10f;   // n_i ≥ 10^10 m^-3
    const float RHO_N_FLOOR = m_n * 1.0e10f;   // n_n ≥ 10^10 m^-3
    const float P_FLOOR     = 1.0e-8f;         // Pa (chromospheric p ~ 10^-2 Pa)
    const float V_MAX       = 1.5e6f;          // m/s, above the 2.35 c_s evaporation ceiling
    const float T_MAX       = 5.0e7f;          // K, ceiling — beyond coronal validity; bounds
                                               // κ_e∝T^{5/2} & TRAC ε so conduction stays finite
    const float F_SLAVE     = 0.99f;           // ionization fraction above which the neutral is a
                                               // collisionally-locked trace species → single fluid

    Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    Vec V     = get_scalar(grid, prim_state, prim::V);
    Vec U     = get_scalar(grid, prim_state, prim::U);
    Vec p_i   = get_scalar(grid, prim_state, prim::P_I);
    Vec p_n   = get_scalar(grid, prim_state, prim::P_N);

    rho_i = arma::clamp(rho_i, RHO_I_FLOOR, arma::datum::inf);
    p_i   = arma::clamp(p_i,   P_FLOOR,     arma::datum::inf);
    V     = arma::clamp(V,    -V_MAX, V_MAX);

    const Vec n_i = rho_i / m_i;
    // Temperature ceiling on the charged fluid: p_i ≤ 2 n_i k_B T_MAX.
    p_i = arma::min(p_i, 2.0f * n_i * (k_b * T_MAX));
    const Vec T_e = p_i / (2.0f * n_i * k_b);     // electron/ion temperature (≤ T_MAX)
    for (arma::uword i = 0; i < grid.ns; ++i) {
        // Floor both densities positive so the V = ρV_i/ρ_i and U = ρU_n/ρ_n decodes
        // and the spectral radii stay well-conditioned, then pick the regime.
        if (rho_n(i) < RHO_N_FLOOR) rho_n(i) = RHO_N_FLOOR;
        const float n_n_i = rho_n(i) / m_n;
        const float f_ion = rho_i(i) / (rho_i(i) + rho_n(i));   // ionization fraction
        if (f_ion >= F_SLAVE) {                  // ion-dominated → neutral is trace: slave neutral
            U(i)   = V(i);                        // U = V
            p_n(i) = n_n_i * k_b * T_e(i);        // T_n = T_e
        } else {                                  // partially/fully neutral → two-fluid, bound both
            if (U(i) >  V_MAX) U(i) =  V_MAX;
            if (U(i) < -V_MAX) U(i) = -V_MAX;
            if (p_n(i) < P_FLOOR) p_n(i) = P_FLOOR;
            float T_n_i = p_n(i) / (n_n_i * k_b);
            if (T_n_i > T_MAX) { T_n_i = T_MAX; p_n(i) = n_n_i * k_b * T_MAX; }   // T_n ≤ T_MAX
            // Trace-ion overheating cap. When the ion is a minority species (f small)
            // its T_i = p_i/(2 n_i k_B) blows up if n_i craters faster than p_i can
            // relax — recombination at the condensation front. Cap the trace ion's
            // temperature at the dominant neutral's T_n (the collisional single-fluid
            // limit). DOWNWARD-only, with a 2× deadband, so it fires ONLY on genuine
            // overheating (T_i > 2 T_n): the cool quiescent lower chromosphere
            // (T_i ≈ T_n, also low f) is left untouched and V is never overwritten —
            // unlike a full slave, which would rewrite the whole cool column and
            // violate energy conservation (it crashed the run at beam onset).
            if (f_ion <= 1.0f - F_SLAVE) {
                const float p_i_cap = 2.0f * (rho_i(i) / m_i) * k_b * T_n_i;
                if (p_i(i) > 2.0f * p_i_cap) p_i(i) = p_i_cap;
            }
        }
    }

    // Electron pressure (three-T only): keep T_e in (0, T_MAX] and p_e < p_i so the
    // derived proton pressure p_proton = p_i − p_e stays positive. Single-T leaves
    // P_E untouched (re-slaved at end of step).
    Vec p_e = get_scalar(grid, prim_state, prim::P_E);
    if (grid.enable_Te) {
        for (arma::uword i = 0; i < grid.ns; ++i) {
            const float n_i_i = rho_i(i) / m_i;
            const float p_e_cap = n_i_i * (k_b * T_MAX);     // T_e ≤ T_MAX
            if (p_e(i) < P_FLOOR)            p_e(i) = P_FLOOR;
            if (p_e(i) > p_e_cap)            p_e(i) = p_e_cap;
            if (p_e(i) > p_i(i) - P_FLOOR)   p_e(i) = p_i(i) - P_FLOOR;
            if (p_e(i) < P_FLOOR)            p_e(i) = P_FLOOR;
        }
    }

    prim_state += scalar_to(grid, rho_i - get_scalar(grid, prim_state, prim::RHO_I), prim::RHO_I);
    prim_state += scalar_to(grid, rho_n - get_scalar(grid, prim_state, prim::RHO_N), prim::RHO_N);
    prim_state += scalar_to(grid, V     - get_scalar(grid, prim_state, prim::V),     prim::V);
    prim_state += scalar_to(grid, U     - get_scalar(grid, prim_state, prim::U),     prim::U);
    prim_state += scalar_to(grid, p_i   - get_scalar(grid, prim_state, prim::P_I),   prim::P_I);
    prim_state += scalar_to(grid, p_n   - get_scalar(grid, prim_state, prim::P_N),   prim::P_N);
    prim_state += scalar_to(grid, p_e   - get_scalar(grid, prim_state, prim::P_E),   prim::P_E);
}

// ----------------------------------------------------------------------------
// Flare beam-heating stage. Adds the nonthermal-electron-beam volumetric
// heating Q_beam [W/m³] (physics.hpp::beam_heating_rate) to the gas over one
// step. The deposited energy is shared between the charged and neutral fluids in
// proportion to their heat capacities per volume, C_i = 3 n_i k_b (ions +
// electrons, since ε_i = 3/2 p_i = 3 n_i k_b T) and C_n = 3/2 n_n k_b, so both
// fluids gain the SAME ΔT = Δt Q / (C_i + C_n). This is essential in the weakly
// ionized chromosphere: there n_e ≪ n_n, so dumping the whole flux into the tiny
// electron pool alone would spike T_e by many× in one step (numerically fatal,
// and unphysical — a precipitating beam heats the bulk gas, which the fast
// electron–neutral collisions then share; the temperature-equilibration stage
// enforces the same on longer scales). Resulting pressure increments are
// Δp_i = 2 n_i k_b ΔT and Δp_n = n_n k_b ΔT. Q_beam is T-independent so this
// explicit additive update is exact and positivity-preserving. When TRAC is
// active the heating is broadened by the same factor ε = (T_c/T)^{5/2} as the
// optically-thin loss, conserving κQ across the unresolved TR (Johnston 2020).
// Mutates p_i and p_n in place. (Nonthermal collisional ionization by the beam
// is omitted; the thermal heating drives Stage E ionization as T rises.)
// ----------------------------------------------------------------------------
void apply_beam_heating_stage(const Grid& grid, Vec& prim_state, float dt) {
    if (!grid.enable_beam_heating) return;
    const Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    const Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    const Vec n_i   = rho_i / grid.m_i;
    const Vec n_n   = rho_n / grid.m_n;
    const float k_b = grid.k_b;

    Vec Q_beam = beam_heating_rate(grid, n_i, n_n);   // W/m^3, +ve = heating
    if (grid.enable_trac) {
        const Vec p_i = get_scalar(grid, prim_state, prim::P_I);
        const Vec p_e = get_scalar(grid, prim_state, prim::P_E);
        const Vec T_e = grid.enable_Te ? Vec(p_e / (n_i * k_b))
                                       : Vec(p_i / (2.0f * n_i * k_b));
        Q_beam /= trac_broadening_factor(grid, T_e);
    }
    if (grid.enable_Te) {
        // Three-temperature: the thick-target beam deposits its energy into the
        // ambient ELECTRONS (Q_beam → ε_e). The subsequent ν_ei equilibration
        // (Stage C, next step) shares it with protons/neutrals; until then T_e can
        // run well above T_i at the impulsive onset — the decoupling the plan
        // targets. The total energy P_I gains the same Δε so E_I stays the
        // conserved total. (cal_dt_i caps Δt on the electron pool so the spike is
        // bounded.) ΔT_e = Δt Q /(1.5 n_i k_B); Δp_e = (2/3)·1.5 n_i k_B ΔT_e = Δt Q.
        const Vec dp_e = grid.gm1() * dt * Q_beam;        // = (γ−1) Δε_e
        prim_state += scalar_to(grid, dp_e, prim::P_E);
        prim_state += scalar_to(grid, dp_e, prim::P_I);
        return;
    }
    // Single-temperature baseline: share Q_beam between the charged and neutral
    // pools by heat capacity so both gain the SAME ΔT (essential in the weakly
    // ionized chromosphere; see header note). ΔT = Δt Q / (C_i + C_n).
    const Vec C_i  = 2.0f * grid.inv_gm1() * n_i * k_b;
    const Vec C_n  = grid.inv_gm1() * n_n * k_b;
    const Vec dT   = (dt * Q_beam) / (C_i + C_n);
    const Vec dp_i = 2.0f * n_i % (k_b * dT);          // = (2/3) C_i ΔT
    const Vec dp_n =        n_n % (k_b * dT);          // = (2/3) C_n ΔT
    prim_state += scalar_to(grid, dp_i, prim::P_I);
    prim_state += scalar_to(grid, dp_n, prim::P_N);
}

// ----------------------------------------------------------------------------
// Ambient coronal heating stage. Adds the steady footpoint-anchored volumetric
// heating H(s) [W/m³] (physics.hpp::coronal_heating_rate) over one step, with
// the SAME pool partitioning and TRAC broadening as the flare beam: in the
// single-temperature baseline the heat is shared between the charged and neutral
// fluids by heat capacity so both gain the same ΔT = Δt H / (C_i + C_n); with
// enable_Te it is deposited into the electron pool (electrons carry the Spitzer
// conduction that delivers the flux into the TR), mirrored into the total energy
// P_I so E_I stays the conserved total. H is T-independent, so this explicit
// additive update is exact and positivity-preserving. No-op unless
// grid.enable_coronal_heating. Mutates the thermal pressures in place.
// ----------------------------------------------------------------------------
void apply_coronal_heating_stage(const Grid& grid, Vec& prim_state, float dt) {
    if (!grid.enable_coronal_heating) return;
    const Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    const Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    const Vec n_i   = rho_i / grid.m_i;
    const Vec n_n   = rho_n / grid.m_n;
    const float k_b = grid.k_b;

    Vec Q = coronal_heating_rate(grid);                // W/m³, +ve = heating
    if (arma::all(Q <= 0.0f)) return;
    if (grid.enable_trac) {
        const Vec p_i = get_scalar(grid, prim_state, prim::P_I);
        const Vec p_e = get_scalar(grid, prim_state, prim::P_E);
        const Vec T_e = grid.enable_Te ? Vec(p_e / (n_i * k_b))
                                       : Vec(p_i / (2.0f * n_i * k_b));
        Q /= trac_broadening_factor(grid, T_e);
    }
    if (grid.enable_Te) {
        // Deposit into the electrons (Δp_e = Δt Q); mirror into the total P_I.
        const Vec dp_e = grid.gm1() * dt * Q;          // = (γ−1) Δε_e
        prim_state += scalar_to(grid, dp_e, prim::P_E);
        prim_state += scalar_to(grid, dp_e, prim::P_I);
        return;
    }
    // Single-temperature baseline: share by heat capacity → equal ΔT.
    const Vec C_i  = 2.0f * grid.inv_gm1() * n_i * k_b;
    const Vec C_n  = grid.inv_gm1() * n_n * k_b;
    const Vec dT   = (dt * Q) / (C_i + C_n);
    const Vec dp_i = 2.0f * n_i % (k_b * dT);          // = (2/3) C_i ΔT
    const Vec dp_n =        n_n % (k_b * dT);          // = (2/3) C_n ΔT
    prim_state += scalar_to(grid, dp_i, prim::P_I);
    prim_state += scalar_to(grid, dp_n, prim::P_N);
}

// ----------------------------------------------------------------------------
// Stage E (writeup §5.3): point-implicit hydrogen ionization / recombination,
// Route B network (docs/studies/eos-ionization/photoionization_c7_inversion_plan.md). The local ODE
// for f ≡ ρ_i / (ρ_i + ρ_n) under quasi-neutrality (n_e = n_i = f n_tot) is
//   df/dt = (1-f) P_phot                          [photoionization,    always on]
//         + f(1-f) n_tot S_CR                     [multilevel collisional, always on]
//         + f(1-f) n_tot S_i                      [Voronov direct,  flag-gated, off by default]
//         - f²      n_tot  α_r                     [radiative recomb, always on]
//         - f³      n_tot² κ_c                     [three-body recomb, flag-gated, off by default]
// Default network (grid.enable_direct_collisional_ionization = false,
// grid.enable_threebody_recombination = false): P_phot + S_CR + α_r only.
// S_i is subdominant by 2–3 orders vs. S_CR everywhere in the chromosphere.
// α_c is within 1–2 orders of α_r below ~700 km but negligible above.
// Backward Euler at lagged T_e reduces to the cell-local polynomial
//   a3 f³ + a2 f² + a1 f + a0 = 0,
//   a3 = Δt n_tot² κ_c,                   (= 0 when three-body is off → quadratic)
//   a2 = Δt n_tot (S_i + S_CR + α_r),
//   a1 = 1 + Δt P_phot - Δt n_tot (S_i + S_CR),
//   a0 = -(f^n + Δt P_phot).
// g(0)=a0≤0 and g(1)≥0 always, so the physical root lies in [0,1]; it is
// found by a safeguarded Newton iteration seeded from the κ_c→0 quadratic
// root (which recovers the pre-Route-B solve exactly when κ_c=0).
// χ_H bookkeeping: collisional ionization (S_i and S_CR) drains χ_H from the
// electron thermal pool; three-body recombination is super-elastic and returns
// χ_H to it (the S_CR drain and α_c return cancel at LTE by detailed balance);
// photoionization adds ions but drains no χ_H (energy comes from the absorbed
// UV photon), and radiative recombination radiates χ_H away (case B, optically
// thin). Mutates V, U, p_i, p_n, ρ_i, ρ_n in place.
// ----------------------------------------------------------------------------

// Solve a3 f³ + a2 f² + a1 f + a0 = 0 for the physical root in [0,1].
// The Stage-E coefficients guarantee g(0)=a0≤0 and g(1)=a3+a2+a1+a0≥0, so a
// root is bracketed by [0,1]. Safeguarded Newton (Newton step when it stays in
// the bracket, bisection otherwise) — robust through the degenerate quadratic
// and linear limits (a3,a2 → 0). `f_seed` is the κ_c→0 quadratic root.
static float solve_ion_cubic(float a3, float a2, float a1, float a0, float f_seed) {
    auto g  = [&](float f) { return ((a3 * f + a2) * f + a1) * f + a0; };
    auto dg = [&](float f) { return (3.0f * a3 * f + 2.0f * a2) * f + a1; };
    if (g(0.0f) >= 0.0f) return 0.0f;     // root at/below 0 (f^n = 0, no photoion)
    if (g(1.0f) <= 0.0f) return 1.0f;     // root at/above 1
    float lo = 0.0f, hi = 1.0f;           // g(lo) < 0 < g(hi)
    float f = (f_seed > 0.0f && f_seed < 1.0f) ? f_seed : 0.5f;
    for (int it = 0; it < 60; ++it) {
        const float gf = g(f);
        if (gf > 0.0f) hi = f; else lo = f;
        const float d = dg(f);
        float f_next = (d != 0.0f) ? f - gf / d : 0.5f * (lo + hi);
        if (!(f_next > lo && f_next < hi)) f_next = 0.5f * (lo + hi);
        if (std::fabs(f_next - f) < 1.0e-7f) return f_next;
        f = f_next;
    }
    return f;
}

void apply_ionization_stage(const Grid& grid, Vec& prim_state, float dt) {
    const Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    const Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    const Vec V     = get_scalar(grid, prim_state, prim::V);
    const Vec U     = get_scalar(grid, prim_state, prim::U);
    const Vec p_i   = get_scalar(grid, prim_state, prim::P_I);
    const Vec p_n   = get_scalar(grid, prim_state, prim::P_N);
    const Vec p_e_old = get_scalar(grid, prim_state, prim::P_E);

    const float m = grid.m_i;  // m_i = m_n for hydrogen
    const Vec n_i = rho_i / m;
    const Vec n_n = rho_n / m;
    const Vec T_i = p_i / (2.0 * n_i * grid.k_b);   // rate temperature (T_charged proxy)
    const Vec T_n = p_n / (n_n * grid.k_b);

    // Rate coefficients at lagged T_e (Route B network).
    // Default network: P_phot + S_CR (ionization), α_r (recombination).
    // S_i and κ_c are zeroed unless their Grid flags are set (see chromosphere.hpp).
    const Vec S       = grid.enable_direct_collisional_ionization
                      ? ionization_rate_S(grid, T_i)
                      : Vec(grid.ns, arma::fill::zeros);       // S_i, direct Voronov (off by default)
    const Vec a       = recombination_rate_alpha(grid, T_i);   // α_r, case-B radiative
    const Vec P       = photoionization_rate_P(grid);          // P_phot
    const Vec S_cr    = ionization_rate_S_CR(grid, T_i);       // multilevel (= κ_c Φ, n_e-free)
    const Vec kappa_c = grid.enable_threebody_recombination
                      ? recombination_rate_kappa_c(grid, T_i)
                      : Vec(grid.ns, arma::fill::zeros);       // κ_c for α_c=κ_c n_e (off by default)

    // Backward-Euler cubic a3 f³ + a2 f² + a1 f + a0 = 0 (see header).
    const Vec rho_tot = rho_i + rho_n;
    const Vec n_tot   = rho_tot / m;
    const Vec f_old   = rho_i / rho_tot;
    const Vec Stot    = S + S_cr;                              // total collisional ionization
    // a3 = dt n_tot² κ_c — group κ_c in first: n_tot² alone overflows float32
    // (n_tot ~ 1e19 → n_tot² ~ 1e38 ≈ FLT_MAX), but n_tot·κ_c is ~1e-18.
    const Vec a3v = dt * n_tot % (n_tot % kappa_c);
    const Vec a2v = dt * n_tot % (Stot + a);
    const Vec a1v = 1.0 + dt * P - dt * n_tot % Stot;
    const Vec a0v = -(f_old + dt * P);

    Vec f_new(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        // Seed: κ_c→0 quadratic root a2 f² + a1 f - (-a0) = 0 (the pre-Route-B
        // solve). Linear fallback when a2→0 (photoionization-only limit).
        const float A = a2v(i), B = a1v(i), C = -a0v(i);
        float seed;
        if (A < 1e-30f) {
            seed = (B > 0.0f) ? (C / B) : f_old(i);
        } else {
            const float disc = B * B + 4.0f * A * C;
            seed = (-B + std::sqrt(std::max(disc, 0.0f))) / (2.0f * A);
        }
        seed = std::min(std::max(seed, 0.0f), 1.0f);
        f_new(i) = solve_ion_cubic(a3v(i), a2v(i), a1v(i), a0v(i), seed);
    }

    // Keep f strictly interior so neither ρ_i nor ρ_n underflows to zero. At the
    // hot top (corona-ghost T) Route B fully ionizes the gas (f→1); a vanishing
    // neutral density would make U_new = ρU_n/ρ_n (and cons2prim downstream)
    // divide by ~0. F_FLOOR bounds the minority species to a dynamically
    // negligible but numerically safe fraction.
    constexpr float F_FLOOR = 1.0e-8f;
    f_new = arma::clamp(f_new, F_FLOOR, 1.0f - F_FLOOR);

    // Integrated channel counts (per m³) over Δt, evaluated at post-step
    // densities and lagged T_e (writeup §5.3, plan §3.4 issue d). Collisional
    // ionization (S_i, S_CR) uses the n_e·n_n product; photoionization is the
    // (1-f)·P branch; recombination (α_r radiative, α_c three-body) uses n_e·n_i.
    const Vec rho_i_new = f_new % rho_tot;
    const Vec rho_n_new = (1.0 - f_new) % rho_tot;
    const Vec n_i_new   = rho_i_new / m;
    const Vec n_n_new   = rho_n_new / m;
    const Vec a_c        = recombination_rate_alpha_c(grid, T_i, n_i_new);  // = κ_c n_e
    const Vec Gamma_coll = dt * n_i_new % n_n_new % S;       // ground collisional — drains χ_H
    const Vec Gamma_mlvl = dt * n_i_new % n_n_new % S_cr;    // multilevel collisional — drains χ_H
    const Vec Gamma_phot = dt * n_n_new % P;                 // photoionization — no χ_H drain
    const Vec Gamma_ion  = Gamma_coll + Gamma_mlvl + Gamma_phot;   // total ionizations
    const Vec Gamma_rec_rad = dt * n_i_new % n_i_new % a;    // radiative — χ_H radiated away
    const Vec Gamma_rec_3b  = dt * n_i_new % n_i_new % a_c;  // three-body — super-elastic, returns χ_H
    const Vec Gamma_rec  = Gamma_rec_rad + Gamma_rec_3b;     // total recombinations

    // Momentum bookkeeping (lagged velocities).
    const Vec rhoV_i_old = rho_i % V;
    const Vec rhoU_n_old = rho_n % U;
    const Vec dmom       = m * (Gamma_ion % U - Gamma_rec % V);
    const Vec rhoV_i_new = rhoV_i_old + dmom;
    const Vec rhoU_n_new = rhoU_n_old - dmom;
    Vec V_new = rhoV_i_new / rho_i_new;
    Vec U_new = rhoU_n_new / rho_n_new;
    // When a species is at the F_FLOOR limit its momentum/density ratio is
    // ill-conditioned (ρU_n/ρ_n with ρ_n→0). Physically the trace fluid is
    // collisionally locked to the dominant one, so set its velocity to the
    // dominant fluid's rather than the (numerically noisy) division result.
    constexpr float F_GUARD = 1.0e-7f;   // a few × F_FLOOR
    for (arma::uword i = 0; i < grid.ns; ++i) {
        if (rho_n_new(i) < F_GUARD * rho_tot(i)) U_new(i) = V_new(i);
        if (rho_i_new(i) < F_GUARD * rho_tot(i)) V_new(i) = U_new(i);
    }

    // Energy bookkeeping (writeup §5.3, eqs. stageE-ei / stageE-en). We update
    // in *thermal* form: convert pressures to thermal energy ε = (3/2) p, apply
    // the source-driven Δε, then rebuild pressures from the new ε and KE.
    const Vec eps_i_old = grid.inv_gm1() * p_i;
    const Vec eps_n_old = grid.inv_gm1() * p_n;
    const Vec KE_i_old  = 0.5 * rho_i % V % V;
    const Vec KE_n_old  = 0.5 * rho_n % U % U;
    const Vec KE_i_new  = 0.5 * rho_i_new % V_new % V_new;
    const Vec KE_n_new  = 0.5 * rho_n_new % U_new % U_new;

    const Vec th_inherit = grid.inv_gm1() * grid.k_b * T_n;   // per-particle thermal carried by new ion
    const Vec th_leave   = grid.inv_gm1() * grid.k_b * T_i;   // per-particle thermal leaving with recomb. neutral
    const Vec KE_inherit = 0.5 * m * U % U;        // per-particle KE inherited by new ion (from neutral)
    const Vec KE_leave   = 0.5 * m * V % V;        // per-particle KE leaving ion fluid

    // Δ(total energy_i) = Q^{e_i}_ion · Δt  (writeup §2.3 boxed; Route B)
    //                   = Γ_ion^Δ  (3/2 k_B T_n + ½ m U²)
    //                   − Γ_rec^Δ  (3/2 k_B T_i + ½ m V²)
    //                   − (Γ_coll + Γ_mlvl)^Δ χ_H        [collisional ioniz drain]
    //                   + Γ_rec,3b^Δ χ_H                 [three-body super-elastic]
    // χ_H exchange with the electron thermal pool: both collisional ionization
    // channels (direct S_i and multilevel S_CR) remove χ_H per event; three-body
    // recombination is super-elastic and returns χ_H to the electrons (the
    // S_CR drain and α_c return cancel at LTE by detailed balance). Photoioniz-
    // ations are powered by absorbed UV photons (χ_H from the radiation field,
    // not the plasma); radiative recombination radiates χ_H away (case B,
    // optically thin). Thermal/KE inheritance applies to all new ions.
    const Vec dE_i = Gamma_ion  % (th_inherit + KE_inherit)
                   - Gamma_rec  % (th_leave   + KE_leave)
                   - (Gamma_coll + Gamma_mlvl - Gamma_rec_3b) * grid.chi_H_J;
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
    Vec p_i_new = grid.gm1() * eps_i_new;
    Vec p_n_new = grid.gm1() * eps_n_new;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        if (p_i_new(i) < p_floor) p_i_new(i) = p_floor;
        if (p_n_new(i) < p_floor) p_n_new(i) = p_floor;
    }

    // Electron pressure. Single-T baseline: carry it (re-slaved each step). Three-T:
    // the χ_H ionization cost is an ELECTRON-pool term (collisional ionization drains
    // χ_H from the free electrons; three-body recombination returns it) — already in
    // the total-energy dE_i above, so attributing it to P_E leaves the protons
    // (p_proton = P_I − P_E) carrying only the thermal/KE exchange with the neutrals.
    Vec p_e_new = p_e_old;
    if (grid.enable_Te) {
        const Vec chi_eps = -(Gamma_coll + Gamma_mlvl - Gamma_rec_3b) * grid.chi_H_J;
        p_e_new = p_e_old + grid.gm1() * chi_eps;
        for (arma::uword i = 0; i < grid.ns; ++i) {
            if (p_e_new(i) < p_floor)                 p_e_new(i) = p_floor;
            if (p_e_new(i) > p_i_new(i) - p_floor)    p_e_new(i) = p_i_new(i) - p_floor;
        }
    }

    prim_state.zeros();
    prim_state += scalar_to(grid, rho_i_new, prim::RHO_I);
    prim_state += scalar_to(grid, rho_n_new, prim::RHO_N);
    prim_state += scalar_to(grid, V_new,     prim::V);
    prim_state += scalar_to(grid, U_new,     prim::U);
    prim_state += scalar_to(grid, p_i_new,   prim::P_I);
    prim_state += scalar_to(grid, p_n_new,   prim::P_N);
    prim_state += scalar_to(grid, p_e_new,   prim::P_E);
}

// Equilibrium-reference well-balancing (chromosphere.hpp::eq_wb): on the first
// step, cache R_eq = rhs_explicit_state(eq_state) so it can be subtracted every
// step, making eq_state an exact discrete fixed point. Computed here (not in the
// const rhs_explicit_state) because it mutates grid; the empty() guard inside
// rhs_explicit_state ensures this capture is itself uncorrected. dt_state is set
// by broadcast_dt before this runs; R_eq's residual dt-dependence is 2nd-order
// small at equilibrium (the predictor half-step is a near-no-op at balance).
static void ensure_eq_residual(Grid& grid) {
    if (grid.eq_wb && grid.eq_residual.is_empty() && !grid.eq_state.is_empty()) {
        grid.eq_residual = rhs_explicit_state(grid, grid.eq_state);
    }
}

// ----------------------------------------------------------------------------
// Backward-Euler split integrator (writeup §3.7, §5.3). The operational order is
// A: explicit MUSCL/Rusanov hydro RHS; B: drag; C: ion-neutral temperature
// equilibration; D: implicit conduction; E: ionization (when enabled); then
// optional beam/coronal heating, radiative cooling, and floors/electron slaving.
// In particular, Stage D changes pressure at fixed velocity and the next Stage-A
// RHS differentiates that updated pressure field.
// ----------------------------------------------------------------------------
Vec advance_Euler_state(Grid& grid, const Vec& xn_state, const Vec& dt_i) {
    if (!grid.eos_gamma_table.empty())
        throw std::logic_error(
            "the legacy two-fluid integrator cannot advance a Gamma1-table Grid; "
            "use mixture_advance");
    broadcast_dt(grid, dt_i);
    ensure_eq_residual(grid);

    Vec U_star = xn_state + grid.dt_state % rhs_explicit_state(grid, xn_state);
    Vec prim   = cons2prim(grid, U_star);

    const float dt = dt_i(0);  // uniform by cal_dt_i construction
    // Flare: sanitize the explicit predictor before any source stage reads it.
    // The MUSCL update can drive p_i/p_n (hence T) negative at the steep
    // evaporation front on the coarse grid, which would feed a negative T into
    // the conduction tridiagonal and produce NaN. The floor clips densities and
    // pressures positive and slaves the trace neutral; no-op for steady runs.
    if (grid.floors_active()) apply_flare_floor_stage(grid, prim);
    // TRAC: recompute the adaptive cutoff T_c from the current profile so the
    // conduction and cooling stages broaden the (currently) unresolved TR.
    if (grid.enable_trac) grid.trac_cutoff_T = compute_trac_cutoff_T(grid, prim);
    apply_drag_stage(grid, prim, dt);
    apply_temperature_stage(grid, prim, dt);
    // Stage D heat conduction. Default-on; the isentropic relaxation experiment
    // (model_column Stage 1) disables it so the adiabatic atmosphere stays a
    // conduction-free Euler steady state (docs/studies/evaporation/gentle_evaporation_downflow.md).
    if (grid.enable_conduction) {
        apply_conduction_stage(grid, prim, dt);
    }
    if (grid.enable_ionization) {
        apply_ionization_stage(grid, prim, dt);
    }
    // Beam heating BEFORE cooling so the heating-vs-radiation competition that
    // defines the explosive/gentle threshold (Fisher et al. 1985) is captured:
    // the cooling stage then radiates whatever fraction of the deposited heat
    // the gas can.
    if (grid.enable_beam_heating) {
        apply_beam_heating_stage(grid, prim, dt);
    }
    // Ambient coronal heating BEFORE cooling, like the beam: the RTV balance is
    // the competition between this heating, conduction, and the radiative sink.
    // Always-on for a steady run (no temporal gate beyond the optional Phase-3
    // ramp inside coronal_heating_rate); no-op when enable_coronal_heating is off.
    if (grid.enable_coronal_heating) {
        apply_coronal_heating_stage(grid, prim, dt);
    }
    if (grid.enable_radiative_cooling) {
        apply_radiative_cooling_stage(grid, prim, dt);
    }
    // Vacuum floor / trace-neutral slaving — keeps the fully-ionized evaporated
    // column well-conditioned (ρ_n→0). No-op for steady scenarios.
    if (grid.floors_active()) {
        apply_flare_floor_stage(grid, prim);
    }

    // Single-temperature limit: slave the electron pool to half the charged
    // thermal energy (p_e = ½ p_total ⇒ T_e ≡ T_i = T_charged). E_E is the 7th
    // variable that is always carried but, with ENABLE_TE=0, never feeds back —
    // this overwrite makes the reported T_e equal T_i to round-off and recovers
    // the pre-T_e single-temperature baseline exactly.
    if (!grid.enable_Te) {
        const Vec p_i = get_scalar(grid, prim, prim::P_I);
        const Vec p_e = get_scalar(grid, prim, prim::P_E);
        prim += scalar_to(grid, 0.5f * p_i - p_e, prim::P_E);
    }

    return prim2cons(grid, prim);
}

// Pure-explicit step: identical to the first half of advance_Euler_state with
// the implicit drag / temperature / conduction stages omitted. Useful for
// comparison runs that need R_I ≡ 0.
Vec advance_Euler_explicit_state(Grid& grid, const Vec& xn_state, const Vec& dt_i) {
    if (!grid.eos_gamma_table.empty())
        throw std::logic_error(
            "advance_Euler_explicit_state is a legacy two-fluid entry point and "
            "cannot advance a Gamma1-table Grid");
    broadcast_dt(grid, dt_i);
    ensure_eq_residual(grid);
    return xn_state + grid.dt_state % rhs_explicit_state(grid, xn_state);
}

Vec advance_RK4(Grid& grid, const Vec& xn_state, const Vec& dt_i) {
    if (!grid.eos_gamma_table.empty())
        throw std::logic_error(
            "advance_RK4 is a legacy two-fluid entry point and cannot advance a "
            "Gamma1-table Grid");
    broadcast_dt(grid, dt_i);
    ensure_eq_residual(grid);
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
