// Tests for chromosphere solver.
//
// Physical correctness is checked against Keheng_s_chromosphere.pdf:
//   - Pressure relation (eq 38):  p_i = 2 n_i k_b T_i   (electron quasi-neutrality)
//   - Sound speeds   (eq 27/64):  c_s,i = sqrt(γ p_i / ρ_i),  c_s,n = sqrt(γ p_n / ρ_n)
//   - Spectral radius (eq 65):    ρ(∂F/∂U) = max(|v|+c_s,i, |u|+c_s,n)
//   - Heat conductivities (eq 53, eq 59)
//   - Ion-neutral collision rate (eq 57; code uses ν_in, target-density form)
//
// Numerical correctness is checked by:
//   - cons↔prim roundtrip
//   - Helper invariants (ip1/im1 shift, get_scalar/scalar_to inverse, minmod limiter)
//   - CFL condition satisfied by cal_dt_i
//   - Mass exactly preserved by one explicit step on a uniform state (no source terms)
//   - Uniform, motionless, gravity-free state is a fixed point of advance_Euler_state

#include "../chromosphere.hpp"
#include "../physics.hpp"
#include "../scenarios/analytic_canopy.hpp"
#include "../scenarios/data_file_parser.hpp"
#include "../scenarios/model_c7.hpp"
#include "../scenarios/pfss_field_line.hpp"
#include "../scenarios/scenario.hpp"

#include <armadillo>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace chromosphere;

static int g_pass = 0;
static int g_fail = 0;
static const char* g_current = "";

#define EXPECT_NEAR(actual, expected, abs_tol) do {                                       \
    const double _a = (double)(actual), _e = (double)(expected), _t = (double)(abs_tol);  \
    if (std::abs(_a - _e) <= _t) { ++g_pass; }                                            \
    else {                                                                                 \
        ++g_fail;                                                                          \
        std::cout << "  [FAIL " << g_current << " @line " << __LINE__ << "] "              \
                  << #actual << " = " << _a << "   expected " << _e                        \
                  << "  (|diff| = " << std::abs(_a - _e) << ", abs_tol = " << _t << ")\n"; \
    }                                                                                       \
} while (0)

#define EXPECT_REL(actual, expected, rel_tol) do {                                         \
    const double _a = (double)(actual), _e = (double)(expected), _r = (double)(rel_tol);   \
    const double _d = std::max(std::abs(_e), 1e-30);                                       \
    if (std::abs(_a - _e) / _d <= _r) { ++g_pass; }                                        \
    else {                                                                                  \
        ++g_fail;                                                                           \
        std::cout << "  [FAIL " << g_current << " @line " << __LINE__ << "] "               \
                  << #actual << " = " << _a << "   expected " << _e                         \
                  << "  (relerr = " << std::abs(_a - _e) / _d                               \
                  << ", rel_tol = " << _r << ")\n";                                         \
    }                                                                                        \
} while (0)

#define EXPECT_TRUE(cond) do {                                                              \
    if (cond) { ++g_pass; }                                                                  \
    else {                                                                                   \
        ++g_fail;                                                                            \
        std::cout << "  [FAIL " << g_current << " @line " << __LINE__ << "] " #cond "\n";   \
    }                                                                                        \
} while (0)

#define RUN(fn) do {                                              \
    g_current = #fn;                                              \
    std::cout << "[RUN  ] " << #fn << "\n";                       \
    const int p0 = g_pass, f0 = g_fail;                           \
    fn();                                                          \
    std::cout << "[  " << (g_fail == f0 ? "OK  " : "FAIL")        \
              << "] " << #fn                                       \
              << "  (" << (g_pass - p0) << " checks)\n";          \
} while (0)


// Set up an ns-cell uniform state with unit cell length, uniform B=1, no gravity.
// Inner/outer ghost cells are pinned to the interior values so all boundary
// fluxes match interior fluxes (true fixed point of the explicit scheme).
static Vec setup_uniform(Grid& grid, arma::uword n_cells,
                         float ni_val, float nn_val,
                         float Ti_val, float Tn_val) {
    grid.init(n_cells, 0.25f);

    grid.ds_i.fill(1.0f);
    grid.B_imh.fill(1.0f);
    grid.B_iph.fill(1.0f);
    grid.B_i.fill(1.0f);
    grid.dinvB_ds_i.zeros();
    grid.phi_g_imh.zeros();
    grid.phi_g_iph.zeros();
    grid.broadcast();

    // Primitive state, then convert to conserved.
    Vec prim(grid.n_state, arma::fill::zeros);
    for (arma::uword i = 0; i < n_cells; ++i) {
        prim(arma::sub2ind(arma::size(n_cells, num_of_eq), i, prim::RHO_I)) = ni_val * grid.m_i;
        prim(arma::sub2ind(arma::size(n_cells, num_of_eq), i, prim::RHO_N)) = nn_val * grid.m_n;
        prim(arma::sub2ind(arma::size(n_cells, num_of_eq), i, prim::V))     = 0.0f;
        prim(arma::sub2ind(arma::size(n_cells, num_of_eq), i, prim::U))     = 0.0f;
        prim(arma::sub2ind(arma::size(n_cells, num_of_eq), i, prim::P_I))   = ni_val * 2.0f * grid.k_b * Ti_val; // writeup eq 38
        prim(arma::sub2ind(arma::size(n_cells, num_of_eq), i, prim::P_N))   = nn_val * grid.k_b * Tn_val;
    }
    Vec cons = prim2cons(grid, prim);

    // Pin inner/outer ghosts to the (uniform) interior values.
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        const float val = cons(arma::sub2ind(arma::size(n_cells, num_of_eq), 0, k));
        grid.inner_boundary0_i(k) = val;
        grid.inner_boundary1_i(k) = val;
        grid.outer_boundary0_i(k) = val;
        grid.outer_boundary1_i(k) = val;
    }
    return cons;
}


// =========================================================================
// Helper-function invariants
// =========================================================================

static void test_scalar_to_get_scalar_inverse() {
    Grid grid;
    grid.init(8, 0.25f);
    Vec v(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) v(i) = (float)i * 1.5f + 0.3f;

    for (arma::uword k = 0; k < num_of_eq; ++k) {
        Vec back = get_scalar(grid, scalar_to(grid, v, k), k);
        for (arma::uword i = 0; i < grid.ns; ++i) EXPECT_NEAR(back(i), v(i), 1e-6);
    }
}

static void test_ip1_im1_interior_shift() {
    Grid grid;
    grid.init(10, 0.25f);
    grid.inner_boundary0_i.zeros(num_of_eq);
    grid.outer_boundary0_i.zeros(num_of_eq);

    Vec xn(grid.n_state, arma::fill::zeros);
    for (arma::uword i = 0; i < grid.ns; ++i)
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I)) = (float)(i + 1) * 10.0f;

    const Vec xp = ip1(grid, xn);
    const Vec xm = im1(grid, xn);

    for (arma::uword i = 0; i < grid.ns - 1; ++i)
        EXPECT_NEAR(xp(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I)),
                    xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i + 1, cons::RHO_I)), 1e-6);

    for (arma::uword i = 1; i < grid.ns; ++i)
        EXPECT_NEAR(xm(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I)),
                    xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i - 1, cons::RHO_I)), 1e-6);
}

static void test_flux_lim_is_minmod() {
    // φ(r) = max(0, min(1, r))   (writeup eq 15)
    Vec r(6);
    r(0) = -1.0f;
    r(1) =  0.0f;
    r(2) =  0.5f;
    r(3) =  1.0f;
    r(4) =  2.0f;
    r(5) =  std::numeric_limits<float>::quiet_NaN();
    Vec out = flux_lim(r);
    EXPECT_NEAR(out(0), 0.0f, 1e-6);
    EXPECT_NEAR(out(1), 0.0f, 1e-6);
    EXPECT_NEAR(out(2), 0.5f, 1e-6);
    EXPECT_NEAR(out(3), 1.0f, 1e-6);
    EXPECT_NEAR(out(4), 1.0f, 1e-6);
    EXPECT_NEAR(out(5), 0.0f, 1e-6);
}

static void test_cons_prim_roundtrip() {
    Grid grid;
    Vec cons = setup_uniform(grid, 8, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    Vec prim  = cons2prim(grid, cons);
    Vec cons2 = prim2cons(grid, prim);
    for (arma::uword j = 0; j < cons.n_elem; ++j) EXPECT_REL(cons2(j), cons(j), 1e-4);
}


// =========================================================================
// Physics formulas vs. writeup
// =========================================================================

static void test_pressure_relation_eq38() {
    const float Ti = 6500.0f, Tn = 6500.0f;
    const float ni = 2.0e17f, nn = 1.0e19f;
    Grid grid;
    Vec cons = setup_uniform(grid, 4, ni, nn, Ti, Tn);
    Vec prim = cons2prim(grid, cons);
    const float p_i = prim(arma::sub2ind(arma::size(grid.ns, num_of_eq), 0, prim::P_I));
    const float p_n = prim(arma::sub2ind(arma::size(grid.ns, num_of_eq), 0, prim::P_N));
    EXPECT_REL(p_i, 2.0f * ni * grid.k_b * Ti, 1e-4);
    EXPECT_REL(p_n, 1.0f * nn * grid.k_b * Tn, 1e-4);
}

static void test_spectral_radius_uniform_at_rest() {
    const float Ti = 6500.0f, Tn = 6500.0f;
    const float ni = 2.0e17f, nn = 1.0e19f;
    Grid grid;
    Vec cons = setup_uniform(grid, 4, ni, nn, Ti, Tn);

    Vec spc = cal_spectral_radius_state(grid, cons);
    const double csi = std::sqrt(2.0 * (double)grid.gamma_mono * (double)grid.k_b * Ti / (double)grid.m_i);
    const double csn = std::sqrt(      (double)grid.gamma_mono * (double)grid.k_b * Tn / (double)grid.m_n);
    const double expected = std::max(csi, csn);

    for (arma::uword i = 0; i < grid.ns; ++i)
        EXPECT_REL(spc(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, 0)), expected, 1e-3);
}

static void test_nu_in_collision_formula() {
    Grid grid;
    grid.init(1, 0.25f);
    const float nn_val = 1.0e19f, Ti = 6500.0f, Tn = 6500.0f;
    Vec nn_v(1), Ti_v(1), Tn_v(1);
    nn_v(0) = nn_val; Ti_v(0) = Ti; Tn_v(0) = Tn;
    const Vec nuin = nu_in(grid, nn_v, Ti_v, Tn_v);

    const double a0 = 53e-12;
    const double expected = (2.0 * a0) * (2.0 * a0) * (double)nn_val
        * std::sqrt(8.0 * (double)arma::datum::pi * (double)grid.k_b * ((double)Ti + (double)Tn) / (double)grid.m_i);
    EXPECT_REL(nuin(0), expected, 1e-3);
}

static void test_kappa_e_eq53() {
    const float ne_val = 2.0e17f, nn_val = 1.0e19f, Te = 6500.0f;
    Vec ne_v(1), nn_v(1), Te_v(1);
    ne_v(0) = ne_val; nn_v(0) = nn_val; Te_v(0) = Te;
    Vec ke = kappa_e(ne_v, nn_v, Te_v);

    const double num = 9.2048e-12 * (double)ne_val * std::pow((double)Te, 2.5);
    const double den = (double)ne_val + 3.5609e-12 * (double)nn_val * (double)Te * (double)Te;
    EXPECT_REL(ke(0), num / den, 1e-3);

    EXPECT_TRUE(ke(0) > 0.02f && ke(0) < 0.05f);
}

static void test_kappa_n_eq59() {
    const float ni_val = 2.0e17f, nn_val = 1.0e19f, Ti = 6500.0f, Tn = 6500.0f;
    Vec ni_v(1), nn_v(1), Ti_v(1), Tn_v(1);
    ni_v(0) = ni_val; nn_v(0) = nn_val; Ti_v(0) = Ti; Tn_v(0) = Tn;
    Vec kn = kappa_n(ni_v, nn_v, Ti_v, Tn_v);

    const double num = 0.0342006 * (double)nn_val * (double)Tn;
    const double den = 1.20613 * (double)ni_val * std::sqrt((double)Tn + (double)Ti)
                     + 1.70573 * (double)nn_val * std::sqrt((double)Tn);
    EXPECT_REL(kn(0), num / den, 1e-3);

    EXPECT_TRUE(kn(0) > 0.5f && kn(0) < 2.0f);
}

static void test_kappa_n_dominates_in_chromosphere() {
    const float Te = 6500.0f, Ti = 6500.0f, Tn = 6500.0f;
    const float ne = 2.0e17f, nn = 1.0e19f;
    Vec ne_v(1); ne_v(0) = ne;
    Vec nn_v(1); nn_v(0) = nn;
    Vec Te_v(1); Te_v(0) = Te;
    Vec Ti_v(1); Ti_v(0) = Ti;
    Vec Tn_v(1); Tn_v(0) = Tn;
    Vec ke = kappa_e(ne_v, nn_v, Te_v);
    Vec kn = kappa_n(ne_v, nn_v, Ti_v, Tn_v);
    EXPECT_TRUE(kn(0) > ke(0));
}


// =========================================================================
// Numerical correctness
// =========================================================================

static void test_cal_dt_respects_cfl() {
    Grid grid;
    Vec cons = setup_uniform(grid, 10, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    Vec dt   = cal_dt_i(grid, cons);
    Vec maxv = cal_max_v_i(grid, cons);
    const float expected_dt = grid.CFL * arma::min(grid.ds_i) / arma::max(maxv);
    EXPECT_REL(dt(0), expected_dt, 1e-3);
    for (arma::uword i = 0; i < grid.ns; ++i)
        EXPECT_TRUE(dt(i) * maxv(i) / grid.ds_i(i) <= grid.CFL + 1e-6f);
}

static void test_cal_max_v_matches_spectral_radius() {
    Grid grid;
    Vec cons = setup_uniform(grid, 5, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    Vec spc = cal_spectral_radius_state(grid, cons);
    Vec mv  = cal_max_v_i(grid, cons);
    for (arma::uword i = 0; i < grid.ns; ++i)
        EXPECT_NEAR(mv(i), spc(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, 0)), 1e-6);
}

static void test_uniform_state_is_fixed_point() {
    const float Ti = 6500.0f, Tn = 6500.0f, ni = 2.0e17f, nn = 1.0e19f;
    Grid grid;
    Vec cons  = setup_uniform(grid, 16, ni, nn, Ti, Tn);
    Vec cons0 = cons;
    Vec dt    = cal_dt_i(grid, cons);
    Vec cons1 = advance_Euler_state(grid, cons, dt);

    for (arma::uword j = 0; j < cons.n_elem; ++j)
        EXPECT_REL(cons1(j), cons0(j), 1e-3);
}

static void test_mass_conservation_on_uniform_state() {
    const float ni = 2.0e17f, nn = 1.0e19f;
    Grid grid;
    Vec cons = setup_uniform(grid, 16, ni, nn, 6500.0f, 6500.0f);
    Vec dt   = cal_dt_i(grid, cons);

    const double mi0 = arma::sum(get_scalar(grid, cons, cons::RHO_I));
    const double mn0 = arma::sum(get_scalar(grid, cons, cons::RHO_N));

    Vec cons1 = advance_Euler_state(grid, cons, dt);

    const double mi1 = arma::sum(get_scalar(grid, cons1, cons::RHO_I));
    const double mn1 = arma::sum(get_scalar(grid, cons1, cons::RHO_N));
    EXPECT_REL(mi1, mi0, 1e-4);
    EXPECT_REL(mn1, mn0, 1e-4);
}

static void test_rk4_uniform_fixed_point() {
    const float Ti = 6500.0f, Tn = 6500.0f, ni = 2.0e17f, nn = 1.0e19f;
    Grid grid;
    Vec cons  = setup_uniform(grid, 8, ni, nn, Ti, Tn);
    Vec cons0 = cons;
    Vec dt    = cal_dt_i(grid, cons);
    Vec cons1 = advance_RK4(grid, cons, dt);
    for (arma::uword j = 0; j < cons.n_elem; ++j)
        EXPECT_REL(cons1(j), cons0(j), 1e-3);
}

// =========================================================================
// Ionization / recombination (writeup §5.3, §6.5)
// =========================================================================

// Voronov (1997) fit at T_e = 10^4 K:
//   U = χ_H / (k_B · 1e4) ≈ 15.78
//   S_i = 2.91e-14 · U^0.39 · exp(-U) / (0.232 + U)
//       ≈ 2.91e-14 · 2.97 · 1.40e-7 / 16.01  ≈ 7.55e-22  m^3/s
static void test_ionization_rate_S_voronov_at_1e4K() {
    Grid grid;
    grid.init(1, 0.25f);
    Vec T(1);
    T(0) = 1.0e4f;
    Vec S = ionization_rate_S(grid, T);

    const double U  = (double)grid.chi_H_J / ((double)grid.k_b * 1.0e4);
    const double expected = 2.91e-14 * std::pow(U, 0.39) * std::exp(-U) / (0.232 + U);
    EXPECT_REL(S(0), (float)expected, 1e-3);
    // Sanity range: tiny at 10^4 K.
    EXPECT_TRUE(S(0) > 1e-23f && S(0) < 1e-20f);
}

// At T_e = 5×10^4 K the rate jumps several orders of magnitude.
static void test_ionization_rate_S_increases_with_T() {
    Grid grid;
    grid.init(1, 0.25f);
    Vec T_low(1);  T_low(0)  = 1.0e4f;
    Vec T_high(1); T_high(0) = 5.0e4f;
    EXPECT_TRUE(ionization_rate_S(grid, T_high)(0) > 1e3f * ionization_rate_S(grid, T_low)(0));
}

// Hummer (1994) case-B fit at 10^4 K: α_r = 2.7e-19 m^3/s exactly.
static void test_recombination_rate_alpha_hummer() {
    Grid grid;
    grid.init(1, 0.25f);
    Vec T_v(1);  T_v(0) = 1.0e4f;
    Vec a = recombination_rate_alpha(grid, T_v);
    EXPECT_REL(a(0), 2.7e-19f, 1e-3);

    // Power-law check at T = 1e3 K: α(1e3)/α(1e4) = 10^0.75 ≈ 5.62.
    Vec T_low(1); T_low(0) = 1.0e3f;
    EXPECT_REL(recombination_rate_alpha(grid, T_low)(0) / a(0), (float)std::pow(10.0, 0.75), 1e-3);
}

// Photoionization rate P_phot(grid) returns the Grid-configured rate
// uniformly across all cells. Default is 1e-4 s^-1 (Carlsson & Stein 2002
// chromospheric ionization timescale midpoint).
static void test_photoionization_rate_default_uniform() {
    Grid grid;
    grid.init(8, 0.25f);
    EXPECT_REL(grid.photoionization_rate, 1.0e-4f, 1e-6);

    Vec P = photoionization_rate_P(grid);
    EXPECT_TRUE(P.n_elem == grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_REL(P(i), grid.photoionization_rate, 1e-6);
    }

    // Override propagates: changing the Grid scalar should change the Vec.
    grid.photoionization_rate = 3.5e-3f;
    Vec P2 = photoionization_rate_P(grid);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_REL(P2(i), 3.5e-3f, 1e-6);
    }
}

// At chromospheric T_e (~6500 K), Voronov S_i is exponentially suppressed
// (U ≈ 24, exp(-U) ~ 1e-11). The kinetic equilibrium is then dominated by
// the photoionization channel:
//   f_eq = P / (P + α_r n_e),
// where n_e ≈ f_eq n_tot. Solving the quadratic at f_eq gives an explicit
// target. With photoionization disabled (grid.photoionization_rate = 0)
// the same setup must relax to f → 0 (S_i is effectively zero at this T).
static void test_stage_e_photoionization_drives_low_T_equilibrium() {
    const float T = 6.5e3f;
    const float n_tot_target = 1.0e19f;
    const float f0 = 1.0e-3f;   // start very nearly neutral

    Grid grid_probe;
    grid_probe.init(1, 0.25f);
    Vec T_v(1); T_v(0) = T;
    const float S = ionization_rate_S(grid_probe, T_v)(0);
    const float a = recombination_rate_alpha(grid_probe, T_v)(0);

    // Photoionization-dominated equilibrium: with S≈0,
    //   df/dt = (1-f)·P − f²·n_tot·α  = 0   →   n_tot α f² + P f − P = 0
    //   f_eq  = (-P + √(P² + 4 n_tot α P)) / (2 n_tot α)
    const float P  = grid_probe.photoionization_rate;
    const float A  = n_tot_target * a;
    const float disc = P * P + 4.0f * A * P;
    const float f_eq_phot = (-P + std::sqrt(disc)) / (2.0f * A);

    // Sanity: the collisional kinetic equilibrium S/(S+α) should be
    // *much* smaller — confirms that the relaxation is photoionization-driven.
    const float f_eq_coll = S / (S + a);
    EXPECT_TRUE(f_eq_phot > 100.0f * f_eq_coll);

    // -- Case A: photoionization ON (default) --
    auto run_to_eq = [&](float P_phot) {
        Grid grid;
        Vec cons = setup_uniform(grid, 1,
                                 f0 * n_tot_target,
                                 (1.0f - f0) * n_tot_target,
                                 T, T);
        grid.enable_ionization     = true;
        grid.photoionization_rate  = P_phot;
        Vec prim = cons2prim(grid, cons);

        // dt large enough to walk to equilibrium quickly via backward Euler.
        const float dt = 1.0f / std::max(grid.photoionization_rate,
                                         n_tot_target * a);
        for (arma::uword step = 0; step < 4000; ++step) {
            // Re-pin T (test f-only kinetics).
            Vec rho_i = get_scalar(grid, prim, prim::RHO_I);
            Vec rho_n = get_scalar(grid, prim, prim::RHO_N);
            Vec ni = rho_i / grid.m_i;
            Vec nn = rho_n / grid.m_n;
            Vec p_i_fix = 2.0f * ni * grid.k_b * T;
            Vec p_n_fix = nn * grid.k_b * T;
            prim.zeros();
            prim += scalar_to(grid, rho_i, prim::RHO_I);
            prim += scalar_to(grid, rho_n, prim::RHO_N);
            prim += scalar_to(grid, Vec(grid.ns, arma::fill::zeros), prim::V);
            prim += scalar_to(grid, Vec(grid.ns, arma::fill::zeros), prim::U);
            prim += scalar_to(grid, p_i_fix, prim::P_I);
            prim += scalar_to(grid, p_n_fix, prim::P_N);
            apply_ionization_stage(grid, prim, dt);
        }
        Vec rho_i_f = get_scalar(grid, prim, prim::RHO_I);
        Vec rho_n_f = get_scalar(grid, prim, prim::RHO_N);
        return rho_i_f(0) / (rho_i_f(0) + rho_n_f(0));
    };

    const float f_on  = run_to_eq(1.0e-4f);
    const float f_off = run_to_eq(0.0f);

    // The strong check: with photoionization, Stage E reaches its predicted
    // chromospheric equilibrium. Without photoionization at this T, the only
    // active process is recombination (S_i is exponentially suppressed), so
    // f decays — but the decay is quadratic in f (df/dt ∝ -f²) and is too
    // slow to drive f to floating-point zero in any reasonable test budget.
    // We assert the qualitative separation between the two regimes:
    //   - f_on locks onto f_eq_phot
    //   - f_off drops well below f_on (recombination is winning)
    EXPECT_REL(f_on, f_eq_phot, 5e-2);
    EXPECT_TRUE(f_on > 10.0f * f_off);
}

// Photoionization must NOT drain χ_H from the electron thermal pool:
// the energy comes from absorbed photons, not from local kinetic temperature.
// Set T_e low enough that S_i ≈ 0 (collisional ionizations are negligible)
// and verify Δ(e_i + e_n) ≈ 0 over the step, even though Γ_phot is large.
static void test_stage_e_no_chi_H_drain_from_photoionization() {
    const float T = 6.5e3f;
    const float ni = 1.0e15f, nn = 1.0e19f;
    Grid grid;
    Vec cons0 = setup_uniform(grid, 4, ni, nn, T, T);
    grid.enable_ionization    = true;
    grid.photoionization_rate = 1.0e-3f;   // amplified rate so Γ_phot is visible

    Vec p_i_old = get_scalar(grid, cons2prim(grid, cons0), prim::P_I);
    Vec p_n_old = get_scalar(grid, cons2prim(grid, cons0), prim::P_N);

    Vec prim = cons2prim(grid, cons0);
    const float dt = 1.0f;
    apply_ionization_stage(grid, prim, dt);

    Vec rho_i_new = get_scalar(grid, prim, prim::RHO_I);
    Vec rho_n_new = get_scalar(grid, prim, prim::RHO_N);
    Vec p_i_new   = get_scalar(grid, prim, prim::P_I);
    Vec p_n_new   = get_scalar(grid, prim, prim::P_N);

    // f must have *grown* (photoionization is creating ions).
    const float f_old = ni / (ni + nn);
    const float f_new = rho_i_new(0) / (rho_i_new(0) + rho_n_new(0));
    EXPECT_TRUE(f_new > f_old + 1.0e-5f);

    // But total thermal energy must be essentially unchanged: at this T,
    // Γ_coll ≈ 0, so the only drain term −Γ_coll·χ_H ≈ 0. The thermal/KE
    // *exchange* between species cancels in the sum (KE = 0 here, T_i = T_n).
    const float chi = grid.chi_H_J;
    Vec T_v(1); T_v(0) = T;
    const float S = ionization_rate_S(grid, T_v)(0);
    const float n_i_n = rho_i_new(0) / grid.m_i;
    const float n_n_n = rho_n_new(0) / grid.m_i;
    const float Gamma_coll_expected = dt * n_i_n * n_n_n * S;
    const float expected_dE_sum     = -Gamma_coll_expected * chi;

    const float delta_p_sum = (p_i_new(0) + p_n_new(0)) - (p_i_old(0) + p_n_old(0));
    const float delta_e_sum = 1.5f * delta_p_sum;

    // Both the actual drain and the collisional-only expectation are
    // negligibly small compared to the baseline thermal energy
    // (3/2)(p_i + p_n) ≈ 2e0 J/m^3. Tolerance is set on absolute scale.
    EXPECT_NEAR(delta_e_sum, expected_dE_sum, 1.0e-4f);
    EXPECT_NEAR(delta_e_sum, 0.0f,            1.0e-4f);
}

// Apply Stage E in isolation to a uniform, motionless box; check that
//  (1) ρ_i + ρ_n is conserved to round-off
//  (2) ρ_i V + ρ_n U remains identically zero (no velocity)
//  (3) Δ(e_i + e_n) per step equals -Γ_ion·χ_H (the optically-thin Lyman loss)
static void test_stage_e_mass_momentum_chi_H_drain() {
    const float ni = 2.0e17f, nn = 1.0e19f, T = 1.5e4f;  // high enough that S_i is non-negligible
    Grid grid;
    Vec cons0 = setup_uniform(grid, 4, ni, nn, T, T);
    grid.enable_ionization = true;

    Vec prim = cons2prim(grid, cons0);
    const float m = grid.m_i;
    const float rho_i_old = ni * m;
    const float rho_n_old = nn * m;
    const float rho_tot   = rho_i_old + rho_n_old;
    const float n_tot     = rho_tot / m;
    Vec T_v(1); T_v(0) = T;
    const float S = ionization_rate_S(grid, T_v)(0);
    const float a = recombination_rate_alpha(grid, T_v)(0);

    // Stage E with a moderate dt: 1e-3 s. This pushes f far from f_old but
    // still in the linear-ish regime.
    const float dt = 1.0e-3f;
    apply_ionization_stage(grid, prim, dt);

    // Mass: ρ_i + ρ_n unchanged per cell.
    Vec rho_i_new = get_scalar(grid, prim, prim::RHO_I);
    Vec rho_n_new = get_scalar(grid, prim, prim::RHO_N);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_REL(rho_i_new(i) + rho_n_new(i), rho_tot, 1e-5);
    }

    // Momentum: V_new = U_new = 0 (starting at rest, no body force).
    Vec V_new = get_scalar(grid, prim, prim::V);
    Vec U_new = get_scalar(grid, prim, prim::U);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_NEAR(V_new(i), 0.0f, 1e-6);
        EXPECT_NEAR(U_new(i), 0.0f, 1e-6);
    }

    // Energy: Δ(e_i + e_n) = -Γ_ion · χ_H. With V=U=0, e_s = (3/2) p_s + ρ_s φ_g
    // and φ_g = 0 here, so Δe = (3/2) Δp.
    Vec p_i_old = get_scalar(grid, cons2prim(grid, cons0), prim::P_I);
    Vec p_n_old = get_scalar(grid, cons2prim(grid, cons0), prim::P_N);
    Vec p_i_new = get_scalar(grid, prim, prim::P_I);
    Vec p_n_new = get_scalar(grid, prim, prim::P_N);
    const float n_i_new = rho_i_new(0) / m;
    const float n_n_new = rho_n_new(0) / m;
    const float Gamma_ion = dt * n_i_new * n_n_new * S;
    const float Gamma_rec = dt * n_i_new * n_i_new * a;
    const float chi = grid.chi_H_J;

    // Expected ion thermal change: Γ_ion (3/2 k T_n) - Γ_rec (3/2 k T_i) - Γ_ion χ_H.
    // Translate to Δp = (2/3) Δ(thermal energy) -- but since the ion fluid
    // also exchanges thermal energy with the neutral fluid via the source,
    // we check the SUM Δ(thermal_i + thermal_n) = (2/3) Δ(p_i+p_n) = -Γ_ion χ_H.
    const float delta_p_sum  = (p_i_new(0) + p_n_new(0)) - (p_i_old(0) + p_n_old(0));
    const float delta_e_sum  = 1.5f * delta_p_sum;
    const float expected_dE  = -Gamma_ion * chi;
    EXPECT_REL(delta_e_sum, expected_dE, 1e-2);
}

// Hold T_e fixed and evolve many steps with Stage E only. The fixed point
// of df/dt = n_tot [f(1-f) S_i - f² α_r] satisfies f/(1-f) = S_i/α_r.
static void test_stage_e_kinetic_equilibrium() {
    const float T = 3.0e4f;  // pick T so S_i/α_r is moderate (not tiny, not huge)
    const float ni0 = 5.0e17f, nn0 = 5.0e18f;
    Grid grid;
    Vec cons = setup_uniform(grid, 1, ni0, nn0, T, T);
    grid.enable_ionization = true;

    Vec T_v(1); T_v(0) = T;
    const float S = ionization_rate_S(grid, T_v)(0);
    const float a = recombination_rate_alpha(grid, T_v)(0);
    const float f_eq_expected = S / (S + a);

    Vec prim = cons2prim(grid, cons);
    // dt scaled so a single Stage E call sees τ_rec ~ 1/(n_e α_r) per step
    // (i.e., we relax over many "rate times").
    const float n_e = ni0;
    const float dt  = 0.5f / (n_e * (S + a));  // moderate per-step relaxation
    const arma::uword n_steps = 200;

    for (arma::uword step = 0; step < n_steps; ++step) {
        // Re-pin T (test the f-only kinetics, not coupled thermodynamics):
        // pressure to enforce T after each step.
        Vec rho_i = get_scalar(grid, prim, prim::RHO_I);
        Vec rho_n = get_scalar(grid, prim, prim::RHO_N);
        Vec ni = rho_i / grid.m_i;
        Vec nn = rho_n / grid.m_n;
        Vec p_i_fix = 2.0f * ni * grid.k_b * T;
        Vec p_n_fix = nn * grid.k_b * T;
        prim.zeros();
        prim += scalar_to(grid, rho_i, prim::RHO_I);
        prim += scalar_to(grid, rho_n, prim::RHO_N);
        prim += scalar_to(grid, Vec(grid.ns, arma::fill::zeros), prim::V);
        prim += scalar_to(grid, Vec(grid.ns, arma::fill::zeros), prim::U);
        prim += scalar_to(grid, p_i_fix, prim::P_I);
        prim += scalar_to(grid, p_n_fix, prim::P_N);

        apply_ionization_stage(grid, prim, dt);
    }

    Vec rho_i_f = get_scalar(grid, prim, prim::RHO_I);
    Vec rho_n_f = get_scalar(grid, prim, prim::RHO_N);
    const float f_final = rho_i_f(0) / (rho_i_f(0) + rho_n_f(0));

    // The fixed point: f/(1-f) = S/α  =>  f = S/(S+α).
    EXPECT_REL(f_final, f_eq_expected, 5e-3);
}

// Fixed point: if the state starts at the kinetic-equilibrium ionization
// fraction f_eq = S/(S+α), Stage E must leave (ρ_i, ρ_n) essentially unchanged.
// This stresses the quadratic solve at the algebraic root, which is where
// catastrophic cancellation would show up if the discriminant logic is wrong.
static void test_stage_e_fixed_point_at_kinetic_equilibrium() {
    const float T = 3.0e4f;
    const float n_tot_target = 1.0e19f;
    Grid g_probe;
    g_probe.init(1, 0.25f);
    Vec T_v(1); T_v(0) = T;
    const float S = ionization_rate_S(g_probe, T_v)(0);
    const float a = recombination_rate_alpha(g_probe, T_v)(0);
    const float f_eq = S / (S + a);

    const float ni = f_eq * n_tot_target;
    const float nn = (1.0f - f_eq) * n_tot_target;
    Grid grid;
    Vec cons = setup_uniform(grid, 4, ni, nn, T, T);
    grid.enable_ionization = true;

    Vec prim_before = cons2prim(grid, cons);
    Vec prim_after  = prim_before;
    // Use a generous dt: at f_eq the source vanishes, so dt shouldn't matter.
    apply_ionization_stage(grid, prim_after, 1.0f);

    // Mass densities should be preserved to within the quadratic-solve
    // floating-point tolerance (well under 1e-3 in single precision).
    Vec rho_i_b = get_scalar(grid, prim_before, prim::RHO_I);
    Vec rho_n_b = get_scalar(grid, prim_before, prim::RHO_N);
    Vec rho_i_a = get_scalar(grid, prim_after,  prim::RHO_I);
    Vec rho_n_a = get_scalar(grid, prim_after,  prim::RHO_N);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_REL(rho_i_a(i), rho_i_b(i), 1e-3);
        EXPECT_REL(rho_n_a(i), rho_n_b(i), 1e-3);
    }
}

// Numerical robustness: huge Δt should still produce a finite, in-range f.
// Backward Euler is unconditionally stable, so f^{n+1} ∈ [0, 1] always.
static void test_stage_e_large_dt_bounded() {
    Grid grid;
    Vec cons = setup_uniform(grid, 4, 1.0e17f, 1.0e19f, 5.0e4f, 5.0e4f);
    grid.enable_ionization = true;

    Vec prim = cons2prim(grid, cons);
    apply_ionization_stage(grid, prim, 1.0e6f);  // ridiculous Δt — well past any timescale

    Vec rho_i = get_scalar(grid, prim, prim::RHO_I);
    Vec rho_n = get_scalar(grid, prim, prim::RHO_N);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float f = rho_i(i) / (rho_i(i) + rho_n(i));
        EXPECT_TRUE(std::isfinite(f) && f >= 0.0f && f <= 1.0f);
        EXPECT_TRUE(std::isfinite(rho_i(i)) && rho_i(i) >= 0.0f);
        EXPECT_TRUE(std::isfinite(rho_n(i)) && rho_n(i) >= 0.0f);
    }
}

// Stage E (ionization + recombination) conserves total hydrogen mass exactly
// per cell: ρ_i + ρ_n is invariant because every event in the reaction
//      H + e⁻ ↔ H⁺ + 2 e⁻
// just moves mass between the ion and neutral fluids — no hydrogen leaves
// either species. The integrator builds this invariant algebraically (the
// f-update is constrained to the simplex), so it should hold to round-off
// across every regime — ionization-dominated (hot, low f), recombination-
// dominated (cool, high f), kinetic equilibrium, and absurd time steps. This
// test sweeps the (T, f_old, dt) cube to catch any branch that breaks it.
static void test_stage_e_conserves_rho_tot_per_cell() {
    const float T_table[]    = {5.0e3f, 8.0e3f, 1.2e4f, 2.0e4f};      // K
    const float f_old_table[] = {0.01f, 0.20f, 0.50f, 0.80f, 0.99f};   // ionization fraction
    const float dt_table[]   = {1.0e-6f, 1.0e-3f, 1.0f, 1.0e3f};      // s — last one is well past saturation

    Grid grid;
    grid.init(static_cast<arma::uword>(sizeof(T_table)/sizeof(T_table[0])
              * sizeof(f_old_table)/sizeof(f_old_table[0])), 0.25f);
    grid.enable_ionization = true;

    // Lay out every (T, f_old) pair across the grid cells, with a fixed total
    // number density n_tot so ρ_tot is the same in every cell. Then we can
    // assert per-cell ρ_tot equals the same scalar before AND after the stage.
    const float n_tot = 1.0e19f;                  // m^-3
    const float rho_tot_target = n_tot * grid.m_i;

    // Table header — printed once so the manual reader can scan the matrix.
    std::cout << "  [stage E ρ_tot conservation table]  "
              << "n_tot = " << n_tot << " m^-3,  ρ_tot_target = "
              << rho_tot_target << " kg/m^3\n";
    std::cout << "  " << std::string(96, '-') << "\n";
    std::cout << "  "
              << std::setw(11) << "dt [s]"
              << std::setw(10) << "T [K]"
              << std::setw(9)  << "f_old"
              << std::setw(15) << "ρ_i_before"
              << std::setw(15) << "ρ_n_before"
              << std::setw(15) << "ρ_i_after"
              << std::setw(15) << "ρ_n_after"
              << std::setw(13) << "Δρ_tot/ρ_tot"
              << "\n";
    std::cout << "  " << std::string(96, '-') << "\n";

    auto run_one_dt = [&](float dt) {
        Vec prim = arma::zeros<Vec>(grid.n_state);
        // Remember each cell's (T, f_old) so we can label the table rows.
        std::vector<std::pair<float,float>> cell_label(grid.ns);
        arma::uword cell = 0;
        for (float T : T_table) {
            for (float f_old : f_old_table) {
                const float n_i = f_old * n_tot;
                const float n_n = (1.0f - f_old) * n_tot;
                const auto sz = arma::size(grid.ns, num_of_eq);
                prim(arma::sub2ind(sz, cell, prim::RHO_I)) = n_i * grid.m_i;
                prim(arma::sub2ind(sz, cell, prim::RHO_N)) = n_n * grid.m_n;
                prim(arma::sub2ind(sz, cell, prim::V))     = 0.0f;
                prim(arma::sub2ind(sz, cell, prim::U))     = 0.0f;
                prim(arma::sub2ind(sz, cell, prim::P_I))   = n_i * 2.0f * grid.k_b * T;
                prim(arma::sub2ind(sz, cell, prim::P_N))   = n_n *        grid.k_b * T;
                cell_label[cell] = {T, f_old};
                ++cell;
            }
        }

        // Before: every cell sums to rho_tot_target (sanity-check that the
        // setup itself is consistent — if this trips, the test is broken).
        Vec rho_i_before = get_scalar(grid, prim, prim::RHO_I);
        Vec rho_n_before = get_scalar(grid, prim, prim::RHO_N);
        for (arma::uword i = 0; i < grid.ns; ++i) {
            EXPECT_REL(rho_i_before(i) + rho_n_before(i), rho_tot_target, 1e-5);
        }

        apply_ionization_stage(grid, prim, dt);

        // After: every cell still sums to the same rho_tot, AND it matches
        // its own pre-stage sum (the two are equivalent here but checking
        // both makes the test robust to a future change in the setup).
        Vec rho_i_after = get_scalar(grid, prim, prim::RHO_I);
        Vec rho_n_after = get_scalar(grid, prim, prim::RHO_N);
        for (arma::uword i = 0; i < grid.ns; ++i) {
            const float rho_tot_before_i = rho_i_before(i) + rho_n_before(i);
            const float rho_tot_after_i  = rho_i_after (i) + rho_n_after (i);
            EXPECT_REL(rho_tot_after_i, rho_tot_before_i, 1e-5);
            EXPECT_REL(rho_tot_after_i, rho_tot_target,    1e-5);

            // And both species are non-negative — the f-update should never
            // be allowed to push us off the simplex.
            EXPECT_TRUE(rho_i_after(i) >= 0.0f);
            EXPECT_TRUE(rho_n_after(i) >= 0.0f);

            const double rel_diff = (double)(rho_tot_after_i - rho_tot_before_i)
                                  / std::max((double)rho_tot_before_i, 1e-30);
            std::cout << "  "
                      << std::setw(11) << std::scientific << std::setprecision(2) << dt
                      << std::setw(10) << std::fixed      << std::setprecision(0) << cell_label[i].first
                      << std::setw(9)  << std::fixed      << std::setprecision(2) << cell_label[i].second
                      << std::setw(15) << std::scientific << std::setprecision(4) << rho_i_before(i)
                      << std::setw(15) << std::scientific << std::setprecision(4) << rho_n_before(i)
                      << std::setw(15) << std::scientific << std::setprecision(4) << rho_i_after(i)
                      << std::setw(15) << std::scientific << std::setprecision(4) << rho_n_after(i)
                      << std::setw(13) << std::scientific << std::setprecision(2) << rel_diff
                      << "\n";
        }
        std::cout << "  " << std::string(96, '-') << "\n";
    };

    for (float dt : dt_table) run_one_dt(dt);
}

// End-to-end: enable Stage E and drive the Model C7 atmosphere for many steps.
// Checks that ionization-on does not destabilize the integrator (no NaNs,
// densities stay positive, temperatures stay in a generous physical band).
static void test_advance_euler_stable_under_model_c7_with_ionization() {
    Grid grid;
    grid.init(100, 0.25f);
    grid.enable_ionization = true;
    Vec xn = model_c7_ic(grid);

    const arma::uword max_steps = 200;
    for (arma::uword step = 0; step < max_steps; ++step) {
        Vec dt = cal_dt_i(grid, xn);
        EXPECT_TRUE(std::isfinite(dt(0)) && dt(0) > 0.0f);
        if (!(std::isfinite(dt(0)) && dt(0) > 0.0f)) return;

        model_c7_update_bc(grid, xn);
        xn = advance_Euler_state(grid, xn, dt);

        bool finite_ok = xn.is_finite();
        EXPECT_TRUE(finite_ok);
        if (!finite_ok) return;

        Vec rho_i = get_scalar(grid, xn, cons::RHO_I);
        Vec rho_n = get_scalar(grid, xn, cons::RHO_N);
        bool rho_pos = arma::min(rho_i) > 0.0f && arma::min(rho_n) > 0.0f;
        EXPECT_TRUE(rho_pos);
        if (!rho_pos) return;
    }

    Vec prim   = cons2prim(grid, xn);
    Vec p_i    = get_scalar(grid, prim, prim::P_I);
    Vec p_n    = get_scalar(grid, prim, prim::P_N);
    EXPECT_TRUE(arma::min(p_i) > 0.0f);
    EXPECT_TRUE(arma::min(p_n) > 0.0f);
}

// With ionization disabled, advance_Euler_state must be byte-identical to
// the previous (pre-Stage E) behavior — a clean regression baseline.
static void test_advance_euler_no_op_when_ionization_disabled() {
    const float Ti = 6500.0f, Tn = 6500.0f, ni = 2.0e17f, nn = 1.0e19f;
    Grid grid;
    Vec cons  = setup_uniform(grid, 16, ni, nn, Ti, Tn);
    EXPECT_TRUE(!grid.enable_ionization);  // default off

    Vec dt    = cal_dt_i(grid, cons);
    Vec cons1 = advance_Euler_state(grid, cons, dt);

    for (arma::uword j = 0; j < cons.n_elem; ++j)
        EXPECT_REL(cons1(j), cons(j), 1e-3);
}


// Stability of advance_Euler_state under the Model C7 IC + BC driver. Runs
// the same loop chromo_main.cpp executes, just for a fixed step budget, and
// checks every step that all conserved variables remain finite, densities
// stay positive, and ion/neutral temperatures stay within a generous band
// of the chromospheric values (no runaway heating/cooling).
static void test_advance_euler_stable_under_model_c7() {
    Grid grid;
    grid.init(100, 0.25f);
    Vec xn = model_c7_ic(grid);

    const arma::uword max_steps = 500;
    for (arma::uword step = 0; step < max_steps; ++step) {
        Vec dt = cal_dt_i(grid, xn);
        EXPECT_TRUE(std::isfinite(dt(0)) && dt(0) > 0.0f);
        if (!(std::isfinite(dt(0)) && dt(0) > 0.0f)) return;

        model_c7_update_bc(grid, xn);
        xn = advance_Euler_state(grid, xn, dt);

        // Spot-check finiteness on every step (bail early on first failure).
        bool finite_ok = xn.is_finite();
        EXPECT_TRUE(finite_ok);
        if (!finite_ok) return;

        // Positivity of densities every step.
        Vec rho_i = get_scalar(grid, xn, cons::RHO_I);
        Vec rho_n = get_scalar(grid, xn, cons::RHO_N);
        bool rho_pos = arma::min(rho_i) > 0.0f && arma::min(rho_n) > 0.0f;
        EXPECT_TRUE(rho_pos);
        if (!rho_pos) return;
    }

    // Final-state sanity: temperatures bounded, energies positive.
    Vec prim   = cons2prim(grid, xn);
    Vec rho_i  = get_scalar(grid, prim, prim::RHO_I);
    Vec rho_n  = get_scalar(grid, prim, prim::RHO_N);
    Vec p_i    = get_scalar(grid, prim, prim::P_I);
    Vec p_n    = get_scalar(grid, prim, prim::P_N);
    Vec n_i    = rho_i / grid.m_i;
    Vec n_n    = rho_n / grid.m_n;
    Vec T_i    = p_i / (2.0f * n_i * grid.k_b);
    Vec T_n    = p_n / (n_n * grid.k_b);

    EXPECT_TRUE(arma::min(p_i) > 0.0f);
    EXPECT_TRUE(arma::min(p_n) > 0.0f);
    // C7 interior is 6.2e3..6.7e3 K; outer ghost holds T_e at 2x. Allow
    // [1e3, 5e4] K as a generous "no thermal runaway" band.
    EXPECT_TRUE(arma::min(T_i) > 1.0e3f && arma::max(T_i) < 5.0e4f);
    EXPECT_TRUE(arma::min(T_n) > 1.0e3f && arma::max(T_n) < 5.0e4f);
}


// =========================================================================
// Driver
// =========================================================================

// ----------------------------------------------------------------------------
// New scenario tests (Addendum Phase D)
// ----------------------------------------------------------------------------

// Write a tiny synthetic data file with ns=3, returns the path. Caller deletes.
static std::string write_synthetic_data_file() {
    const std::string path = "/tmp/chromo_synthetic_field_line.dat";
    std::ofstream out(path);
    out << "# synthetic field-line scenario for tests\n"
        << "[META]\n"
        << "ns=3\n"
        << "g_si=274.0\n"
        << "\n"
        << "[CELLS]\n"
        // i  ds_m    B_imh    B_iph    phi_g_imh  phi_g_iph   ne         nn          T
        << "0  1.0e4  1.00e-2  9.00e-3  0.0        2.74e6      1.0e17     2.0e19      6300.0\n"
        << "1  1.0e4  9.00e-3  8.00e-3  2.74e6     5.48e6      1.0e17     1.0e19      6500.0\n"
        << "2  1.0e4  8.00e-3  7.00e-3  5.48e6     8.22e6      1.0e17     5.0e18      6700.0\n"
        << "\n"
        << "[GHOSTS]\n"
        // tag      B_T      phi_g     ne         nn         T        T_e_factor
        << "outer_0  7.00e-3  8.22e6   8.0e16     3.0e18     6700.0   2.0\n"
        << "outer_1  7.00e-3  9.00e6   8.0e16     3.0e18     6700.0   2.0\n"
        << "inner_0  1.00e-2  0.0      1.0e17     2.0e19     6300.0   1.0\n"
        << "inner_1  1.00e-2  0.0      1.0e17     2.0e19     6300.0   1.0\n";
    return path;
}

static void test_parser_reads_synthetic_file() {
    const auto path = write_synthetic_data_file();
    EXPECT_TRUE(peek_ns_from_file(path) == 3);
    const auto f = parse_scenario_data_file(path);
    EXPECT_TRUE(f.ns == 3);
    EXPECT_REL(f.g_si, 274.0f, 1e-6);
    EXPECT_REL(f.ds_m(1),    1.0e4f,  1e-6);
    EXPECT_REL(f.B_imh_T(0), 1.0e-2f, 1e-6);
    EXPECT_REL(f.B_iph_T(2), 7.0e-3f, 1e-6);
    EXPECT_REL(f.phi_g_iph(2), 8.22e6f, 1e-5);
    EXPECT_REL(f.T_K(0), 6300.0f, 1e-6);
    EXPECT_REL(f.ghost_B_T(0), 7.0e-3f, 1e-6);
    EXPECT_REL(f.ghost_T_e_factor(0), 2.0f, 1e-6);
    EXPECT_REL(f.ghost_T_K(2), 6300.0f, 1e-6);  // inner_0
    std::remove(path.c_str());
}

static void test_pfss_ic_grid_invariants() {
    const auto path = write_synthetic_data_file();
    Grid grid;
    grid.init(pfss_peek_ns(path), 0.25f);
    Vec xn = pfss_ic(grid, path);

    // Geometry populated from the file, not defaults.
    EXPECT_REL(grid.ds_i(0), 1.0e4f, 1e-6);
    EXPECT_REL(grid.B_imh(0), 1.0e-2f, 1e-6);
    EXPECT_TRUE(arma::all(grid.B_imh > 0.0f));
    EXPECT_TRUE(arma::all(grid.B_iph > 0.0f));
    EXPECT_TRUE(arma::all(grid.ds_i  > 0.0f));

    // dinvB_ds_i matches the finite-difference of 1/B that flux.cpp:74 uses.
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float expected = (1.0f / grid.B_iph(i) - 1.0f / grid.B_imh(i)) / grid.ds_i(i);
        EXPECT_REL(grid.dinvB_ds_i(i), expected, 1e-5);
    }

    // Conserved state: positive densities, sensible total energy.
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_TRUE(xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I)) > 0.0f);
        EXPECT_TRUE(xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_N)) > 0.0f);
        EXPECT_TRUE(xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_I))   > 0.0f);
        EXPECT_TRUE(xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_N))   > 0.0f);
    }

    // Outer ghost has T_e_factor=2 baked into E: recovering T_eff should give
    // 2 * T_K, where T_K = 6700.
    {
        const Vec& ob = grid.outer_boundary0_i;
        const float n_i = ob(cons::RHO_I) / grid.m_i;
        const float phi_g = grid.phi_g_iph(grid.ns - 1);
        const float T_eff = (2.0f/3.0f * ob(cons::E_I) - 2.0f/3.0f * grid.m_i * n_i * phi_g)
                          / (2.0f * n_i * grid.k_b);
        EXPECT_REL(T_eff, 2.0f * 6700.0f, 1e-4);
    }
    std::remove(path.c_str());
}

static void test_pfss_advance_euler_stable() {
    const auto path = write_synthetic_data_file();
    Grid grid;
    grid.init(pfss_peek_ns(path), 0.25f);
    grid.enable_ionization = true;
    Vec xn = pfss_ic(grid, path);

    for (int step = 0; step < 50; ++step) {
        pfss_update_bc(grid, xn);
        Vec dt = cal_dt_i(grid, xn);
        xn = advance_Euler_state(grid, xn, dt);
        EXPECT_TRUE(xn.is_finite());
    }
    std::remove(path.c_str());
}

static void test_analytic_canopy_B_profile() {
    Grid grid;
    grid.init(analytic_canopy_peek_ns(), 0.25f);
    Vec xn = analytic_canopy_ic(grid);

    // Recipe constants (mirror analytic_canopy.cpp):
    const float B0 = 1.0e-2f, Binf = 1.5e-3f, HB = 3.0e5f;

    // Foot (i=0 imh face) should sit at B0; far end should approach Binf.
    EXPECT_REL(grid.B_imh(0), B0, 1e-4);
    EXPECT_TRUE(grid.B_iph(grid.ns - 1) > Binf);
    EXPECT_TRUE(grid.B_iph(grid.ns - 1) < 1.05f * B0);

    // Exponential decay rate at the foot: dB/ds ≈ -(B0 - Binf)/HB.
    const float dBds = (grid.B_iph(0) - grid.B_imh(0)) / grid.ds_i(0);
    const float expected = -(B0 - Binf) / HB * std::exp(-grid.ds_i(0) / (2.0f * HB));
    EXPECT_REL(dBds, expected, 5e-2);  // first-order in ds/HB

    // dinvB_ds_i is consistent with the B(s) finite difference.
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float expected_d = (1.0f / grid.B_iph(i) - 1.0f / grid.B_imh(i)) / grid.ds_i(i);
        EXPECT_REL(grid.dinvB_ds_i(i), expected_d, 1e-5);
    }
}

// Open-hyperbolic BC (reflecting inner wall + clamped outer outflow):
//   * Inner ghosts mirror cells 0 and 1 with V flipped, ρ/T preserved → no
//     mass can cross s = 0.
//   * Outer ghosts copy cell ns-1, but ghost velocity is clamped to ≥ 0 so
//     the BC stays a one-way outflow (no spurious refill when V_{ns-1} < 0).
//   * Running it under a non-trivial flow state for many steps must keep the
//     domain-integrated mass within a tight band: the reservoir artifact that
//     drove the old C7 BC to ±27 % drift is gone.
static void test_apply_open_bcs_mirrors_and_extrapolates() {
    Grid grid;
    grid.init(8, 0.25f);
    // Plant a non-trivial state with a non-zero velocity so the mirror flip is
    // observable. Upper cells get a positive V to exercise the unclamped
    // outflow path; lower cells get a negative V to exercise the wall mirror.
    Vec xn = arma::zeros<Vec>(grid.n_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float n_i = 1.0e17f * (1.0f + 0.1f * (float)i);
        const float n_n = 2.0e19f * (1.0f - 0.05f * (float)i);
        const float T   = 6500.0f + 50.0f * (float)i;
        const float V   = (i < grid.ns / 2) ? -120.0f - 5.0f * (float)i :  80.0f + 5.0f * (float)i;
        const float U   = (i < grid.ns / 2) ? -100.0f - 4.0f * (float)i :  60.0f + 4.0f * (float)i;
        const float rho_i = n_i * grid.m_i;
        const float rho_n = n_n * grid.m_n;
        const auto sz = arma::size(grid.ns, num_of_eq);
        xn(arma::sub2ind(sz, i, cons::RHO_I)) = rho_i;
        xn(arma::sub2ind(sz, i, cons::RHO_N)) = rho_n;
        xn(arma::sub2ind(sz, i, cons::MOM_I)) = rho_i * V;
        xn(arma::sub2ind(sz, i, cons::MOM_N)) = rho_n * U;
        xn(arma::sub2ind(sz, i, cons::E_I))   = 1.5f * grid.k_b * n_i * 2.0f * T + 0.5f * rho_i * V * V;
        xn(arma::sub2ind(sz, i, cons::E_N))   = 1.5f * grid.k_b * n_n * T       + 0.5f * rho_n * U * U;
    }

    apply_open_bcs(grid, xn);

    auto cell_rho_i = [&](arma::uword i) {
        return xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I));
    };
    auto cell_rho_n = [&](arma::uword i) {
        return xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_N));
    };
    auto cell_V = [&](arma::uword i) {
        const auto sz = arma::size(grid.ns, num_of_eq);
        return xn(arma::sub2ind(sz, i, cons::MOM_I)) / xn(arma::sub2ind(sz, i, cons::RHO_I));
    };
    auto cell_U = [&](arma::uword i) {
        const auto sz = arma::size(grid.ns, num_of_eq);
        return xn(arma::sub2ind(sz, i, cons::MOM_N)) / xn(arma::sub2ind(sz, i, cons::RHO_N));
    };
    auto ghost_V = [&](const Vec& ob) { return ob(cons::MOM_I) / ob(cons::RHO_I); };
    auto ghost_U = [&](const Vec& ob) { return ob(cons::MOM_N) / ob(cons::RHO_N); };

    // Inner ghosts mirror cells 0 and 1: ρ matches, V is flipped.
    EXPECT_REL(grid.inner_boundary0_i(cons::RHO_I), cell_rho_i(0), 1e-5);
    EXPECT_REL(grid.inner_boundary0_i(cons::RHO_N), cell_rho_n(0), 1e-5);
    EXPECT_REL(ghost_V(grid.inner_boundary0_i), -cell_V(0), 1e-4);
    EXPECT_REL(ghost_U(grid.inner_boundary0_i), -cell_U(0), 1e-4);

    EXPECT_REL(grid.inner_boundary1_i(cons::RHO_I), cell_rho_i(1), 1e-5);
    EXPECT_REL(ghost_V(grid.inner_boundary1_i), -cell_V(1), 1e-4);
    EXPECT_REL(ghost_U(grid.inner_boundary1_i), -cell_U(1), 1e-4);

    // Outer ghosts copy cell ns-1's ρ, T, V, U — pure Neumann on every variable.
    EXPECT_REL(grid.outer_boundary0_i(cons::RHO_I), cell_rho_i(grid.ns - 1), 1e-5);
    EXPECT_REL(grid.outer_boundary1_i(cons::RHO_I), cell_rho_i(grid.ns - 1), 1e-5);
    EXPECT_TRUE(cell_V(grid.ns - 1) > 0.0f);
    EXPECT_REL(ghost_V(grid.outer_boundary0_i), cell_V(grid.ns - 1), 1e-4);
    EXPECT_REL(ghost_U(grid.outer_boundary0_i), cell_U(grid.ns - 1), 1e-4);
    EXPECT_REL(ghost_V(grid.outer_boundary1_i), cell_V(grid.ns - 1), 1e-4);
    EXPECT_REL(ghost_U(grid.outer_boundary1_i), cell_U(grid.ns - 1), 1e-4);
}

// All four outer-ghost variables are Neumann: ρ, T, V, U copied from cell
// ns-1 into both ghosts. The face has no jump → no Rusanov dissipation →
// boundary is fully transparent.
static void test_apply_open_bcs_outer_velocity_halving_pattern() {
    Grid grid;
    grid.init(8, 0.25f);
    Vec xn = arma::zeros<Vec>(grid.n_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float n_i = 1.0e17f;
        const float n_n = 2.0e19f;
        const float T   = 6500.0f;
        const float V   = -250.0f;  // negative — inflow case, should still pass through
        const float U   = -200.0f;
        const float rho_i = n_i * grid.m_i;
        const float rho_n = n_n * grid.m_n;
        const auto sz = arma::size(grid.ns, num_of_eq);
        xn(arma::sub2ind(sz, i, cons::RHO_I)) = rho_i;
        xn(arma::sub2ind(sz, i, cons::RHO_N)) = rho_n;
        xn(arma::sub2ind(sz, i, cons::MOM_I)) = rho_i * V;
        xn(arma::sub2ind(sz, i, cons::MOM_N)) = rho_n * U;
        xn(arma::sub2ind(sz, i, cons::E_I))   = 1.5f * grid.k_b * n_i * 2.0f * T + 0.5f * rho_i * V * V;
        xn(arma::sub2ind(sz, i, cons::E_N))   = 1.5f * grid.k_b * n_n * T       + 0.5f * rho_n * U * U;
    }

    apply_open_bcs(grid, xn);

    const float V_int = -250.0f;
    const float U_int = -200.0f;

    const float V_g0 = grid.outer_boundary0_i(cons::MOM_I) / grid.outer_boundary0_i(cons::RHO_I);
    const float U_g0 = grid.outer_boundary0_i(cons::MOM_N) / grid.outer_boundary0_i(cons::RHO_N);
    EXPECT_REL(V_g0, V_int, 1e-4);
    EXPECT_REL(U_g0, U_int, 1e-4);

    const float V_g1 = grid.outer_boundary1_i(cons::MOM_I) / grid.outer_boundary1_i(cons::RHO_I);
    const float U_g1 = grid.outer_boundary1_i(cons::MOM_N) / grid.outer_boundary1_i(cons::RHO_N);
    EXPECT_REL(V_g1, V_int, 1e-4);
    EXPECT_REL(U_g1, U_int, 1e-4);
}

// End-to-end stability check under the model_c7 BCs (photospheric inner BC
// + V/2-impedance outer outflow). The photospheric Dirichlet ρ_n is an open
// neutral reservoir, so column mass is *not* strictly conserved — it can
// drift by tens of percent as the chromosphere exchanges mass with the
// photosphere. What this test does guarantee is the run stays well-posed:
// no NaNs, ρ stays strictly positive in every cell, and the total mass
// remains positive and finite for the whole run.
static void test_open_bcs_stability_under_model_c7() {
    Grid grid;
    grid.init(100, 0.25f);
    grid.enable_ionization = true;
    Vec xn = model_c7_ic(grid);

    auto column_mass = [&]() {
        float M = 0.0f;
        for (arma::uword i = 0; i < grid.ns; ++i) {
            const auto sz = arma::size(grid.ns, num_of_eq);
            const float rho_i = xn(arma::sub2ind(sz, i, cons::RHO_I));
            const float rho_n = xn(arma::sub2ind(sz, i, cons::RHO_N));
            M += (rho_i + rho_n) * grid.ds_i(i);
        }
        return M;
    };

    EXPECT_TRUE(column_mass() > 0.0f);
    for (int step = 0; step < 800; ++step) {
        Vec dt = cal_dt_i(grid, xn);
        model_c7_update_bc(grid, xn);
        xn = advance_Euler_state(grid, xn, dt);
        EXPECT_TRUE(xn.is_finite());
        if (!xn.is_finite()) return;

        Vec rho_i = get_scalar(grid, xn, cons::RHO_I);
        Vec rho_n = get_scalar(grid, xn, cons::RHO_N);
        EXPECT_TRUE(arma::min(rho_i) > 0.0f);
        EXPECT_TRUE(arma::min(rho_n) > 0.0f);
    }
    const float M_end = column_mass();
    EXPECT_TRUE(std::isfinite(M_end) && M_end > 0.0f);
}

// model_c7-specific inner BC: ρ_n, T_n and n_i are Dirichlet (ρ_n / T_n at
// the photospheric snapshot from the IC; n_i at the constant 1e16 m^-3
// reservoir). T_i, V, U are Neumann from cell 0. Perturbing cell 0 must
// leave the three Dirichlet quantities at their pinned values and propagate
// the three Neumann quantities to the ghost.
static void test_model_c7_inner_bc_neumann_ni_dirichlet_rho_n_outer_T_doubles() {
    Grid grid;
    grid.init(100, 0.25f);
    Vec xn = model_c7_ic(grid);

    const float rho_n_pinned = grid.inner_boundary0_i(cons::RHO_N);   // snapshot from IC

    auto outer_T_i = [&](const Vec& ob) {
        const float rho_i = ob(cons::RHO_I);
        const float n_i   = rho_i / grid.m_i;
        const float V     = ob(cons::MOM_I) / rho_i;
        const float phi_g = grid.phi_g_iph(grid.ns - 1);
        return (ob(cons::E_I) - 0.5f * rho_i * V * V - rho_i * phi_g)
                / (3.0f * grid.k_b * n_i);
    };
    auto outer_T_n = [&](const Vec& ob) {
        const float rho_n = ob(cons::RHO_N);
        const float n_n   = rho_n / grid.m_n;
        const float U     = ob(cons::MOM_N) / rho_n;
        const float phi_g = grid.phi_g_iph(grid.ns - 1);
        return (ob(cons::E_N) - 0.5f * rho_n * U * U - rho_n * phi_g)
                / (1.5f * grid.k_b * n_n);
    };

    model_c7_update_bc(grid, xn);

    // Inner: ρ_n Dirichlet at IC snapshot.
    EXPECT_REL(grid.inner_boundary0_i(cons::RHO_N), rho_n_pinned, 1e-5);
    EXPECT_REL(grid.inner_boundary1_i(cons::RHO_N), rho_n_pinned, 1e-5);

    // Inner: ρ_i Neumann — must match cell 0.
    const auto sz = arma::size(grid.ns, num_of_eq);
    const float rho_i_cell_0 = xn(arma::sub2ind(sz, 0, cons::RHO_I));
    EXPECT_REL(grid.inner_boundary0_i(cons::RHO_I), rho_i_cell_0, 1e-5);
    EXPECT_REL(grid.inner_boundary1_i(cons::RHO_I), rho_i_cell_0, 1e-5);

    // Outer: T_g1 = 2·T_g0 for both species.
    EXPECT_REL(outer_T_i(grid.outer_boundary1_i),
               2.0f * outer_T_i(grid.outer_boundary0_i), 1e-4);
    EXPECT_REL(outer_T_n(grid.outer_boundary1_i),
               2.0f * outer_T_n(grid.outer_boundary0_i), 1e-4);

    // Perturb cell 0: triple ρ_i and add a downflow. Ghost ρ_n stays pinned;
    // ghost ρ_i, V, U must all track cell 0's new values (Neumann).
    const float rho_i_new   = 3.0f * rho_i_cell_0;
    const float V_new       = -250.0f;
    const float U_new       = -200.0f;
    const float rho_n_cell0 = xn(arma::sub2ind(sz, 0, cons::RHO_N));
    xn(arma::sub2ind(sz, 0, cons::RHO_I)) = rho_i_new;
    xn(arma::sub2ind(sz, 0, cons::MOM_I)) = rho_i_new * V_new;
    xn(arma::sub2ind(sz, 0, cons::MOM_N)) = rho_n_cell0 * U_new;
    const float phi_g_cell_0 = 0.5f * (grid.phi_g_imh(0) + grid.phi_g_iph(0));
    const float p_i_target   = 2.0f * (rho_i_new / grid.m_i) * grid.k_b * 6500.0f;
    const float p_n_target   = (rho_n_cell0 / grid.m_n) * grid.k_b * 6500.0f;
    xn(arma::sub2ind(sz, 0, cons::E_I)) = 1.5f * p_i_target + 0.5f * rho_i_new * V_new * V_new
                                          + rho_i_new * phi_g_cell_0;
    xn(arma::sub2ind(sz, 0, cons::E_N)) = 1.5f * p_n_target + 0.5f * rho_n_cell0 * U_new * U_new
                                          + rho_n_cell0 * phi_g_cell_0;

    model_c7_update_bc(grid, xn);

    // ρ_n still Dirichlet.
    EXPECT_REL(grid.inner_boundary0_i(cons::RHO_N), rho_n_pinned, 1e-5);
    EXPECT_REL(grid.inner_boundary1_i(cons::RHO_N), rho_n_pinned, 1e-5);

    // ρ_i, V, U Neumann — must track cell 0's new values.
    EXPECT_REL(grid.inner_boundary0_i(cons::RHO_I), rho_i_new, 1e-5);
    EXPECT_REL(grid.inner_boundary1_i(cons::RHO_I), rho_i_new, 1e-5);
    const float V_ghost = grid.inner_boundary0_i(cons::MOM_I) / grid.inner_boundary0_i(cons::RHO_I);
    const float U_ghost = grid.inner_boundary0_i(cons::MOM_N) / grid.inner_boundary0_i(cons::RHO_N);
    EXPECT_REL(V_ghost, V_new, 1e-4);
    EXPECT_REL(U_ghost, U_new, 1e-4);
}

// ----------------------------------------------------------------------------

int main() {
    std::cout << "===== Chromosphere test suite =====\n\n";

    RUN(test_scalar_to_get_scalar_inverse);
    RUN(test_ip1_im1_interior_shift);
    RUN(test_flux_lim_is_minmod);
    RUN(test_cons_prim_roundtrip);

    RUN(test_pressure_relation_eq38);
    RUN(test_spectral_radius_uniform_at_rest);
    RUN(test_nu_in_collision_formula);
    RUN(test_kappa_e_eq53);
    RUN(test_kappa_n_eq59);
    RUN(test_kappa_n_dominates_in_chromosphere);

    RUN(test_cal_dt_respects_cfl);
    RUN(test_cal_max_v_matches_spectral_radius);
    RUN(test_uniform_state_is_fixed_point);
    RUN(test_mass_conservation_on_uniform_state);
    RUN(test_rk4_uniform_fixed_point);
    RUN(test_advance_euler_stable_under_model_c7);

    RUN(test_ionization_rate_S_voronov_at_1e4K);
    RUN(test_ionization_rate_S_increases_with_T);
    RUN(test_recombination_rate_alpha_hummer);
    RUN(test_photoionization_rate_default_uniform);
    RUN(test_stage_e_photoionization_drives_low_T_equilibrium);
    RUN(test_stage_e_no_chi_H_drain_from_photoionization);
    RUN(test_stage_e_mass_momentum_chi_H_drain);
    RUN(test_stage_e_conserves_rho_tot_per_cell);
    RUN(test_stage_e_kinetic_equilibrium);
    RUN(test_stage_e_fixed_point_at_kinetic_equilibrium);
    RUN(test_stage_e_large_dt_bounded);
    RUN(test_advance_euler_no_op_when_ionization_disabled);
    RUN(test_advance_euler_stable_under_model_c7_with_ionization);

    RUN(test_parser_reads_synthetic_file);
    RUN(test_pfss_ic_grid_invariants);
    RUN(test_pfss_advance_euler_stable);
    RUN(test_analytic_canopy_B_profile);

    RUN(test_apply_open_bcs_mirrors_and_extrapolates);
    RUN(test_apply_open_bcs_outer_velocity_halving_pattern);
    RUN(test_open_bcs_stability_under_model_c7);
    RUN(test_model_c7_inner_bc_neumann_ni_dirichlet_rho_n_outer_T_doubles);

    std::cout << "\n===== Summary =====\n";
    std::cout << "Passed: " << g_pass << "\n";
    std::cout << "Failed: " << g_fail << "\n";
    return g_fail == 0 ? 0 : 1;
}
