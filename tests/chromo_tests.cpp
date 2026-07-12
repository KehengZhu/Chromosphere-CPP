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
#include "../scenarios/mesh.hpp"
#include "../scenarios/model_c7.hpp"
#include "../scenarios/model_column.hpp"
#include "../scenarios/pfss_field_line.hpp"
#include "../scenarios/scenario.hpp"

#include <armadillo>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
        // Electron partial pressure p_e = n_e k_B T_e with T_e = T_i (= ½ P_I), so
        // the state is self-consistent for the 7-variable (T_e ≠ T_i) model. With
        // ENABLE_TE off this matches the value advance_Euler_state slaves each step.
        prim(arma::sub2ind(arma::size(n_cells, num_of_eq), i, prim::P_E))   = ni_val * grid.k_b * Ti_val;
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

static void test_flux_lim_mc3() {
    // MC3 / Koren limiter (BATSRUS 'mc3'), asymmetric third-order (β=2):
    //   φ₊(r) = max(0, min(2r, 2, (2r+1)/3))    ['+'/right-face terms]
    //   φ₋(r) = max(0, min(2r, 2, (r+2)/3))     ['−'/left-face  terms]
    Vec r(7);
    r(0) = -1.0f; r(1) = 0.0f; r(2) = 0.5f; r(3) = 1.0f;
    r(4) = 2.0f;  r(5) = 10.0f;
    r(6) = std::numeric_limits<float>::quiet_NaN();

    Vec p = flux_lim_mc3_plus(r, 2.0f);
    EXPECT_NEAR(p(0), 0.0f,        1e-6);   // r<0 (sign disagreement) ⇒ 0
    EXPECT_NEAR(p(1), 0.0f,        1e-6);
    EXPECT_NEAR(p(2), 2.0f/3.0f,   1e-6);   // (2·0.5+1)/3
    EXPECT_NEAR(p(3), 1.0f,        1e-6);   // full slope at r=1
    EXPECT_NEAR(p(4), 5.0f/3.0f,   1e-6);   // (2·2+1)/3
    EXPECT_NEAR(p(5), 2.0f,        1e-6);   // β cap
    EXPECT_NEAR(p(6), 0.0f,        1e-6);   // NaN ⇒ 0

    Vec m = flux_lim_mc3_minus(r, 2.0f);
    EXPECT_NEAR(m(0), 0.0f,        1e-6);
    EXPECT_NEAR(m(1), 0.0f,        1e-6);
    EXPECT_NEAR(m(2), 5.0f/6.0f,   1e-6);   // (0.5+2)/3
    EXPECT_NEAR(m(3), 1.0f,        1e-6);   // matches φ₊ at r=1 (symmetric there)
    EXPECT_NEAR(m(4), 4.0f/3.0f,   1e-6);   // (2+2)/3
    EXPECT_NEAR(m(5), 2.0f,        1e-6);   // β cap
    EXPECT_NEAR(m(6), 0.0f,        1e-6);
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
    // e–n term coefficient uses the accurate Krstic & Schultz e–H cross section
    // (Vranjes & Krstic 2013): σ_en ≈ 2.5e-19 m² → 2.836e-11 (was 3.5609e-12).
    const double den = (double)ne_val + 2.836e-11 * (double)nn_val * (double)Te * (double)Te;
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

// ----------------------------------------------------------------------------
// Route B network helpers (writeup §3.1). The Stage-E ionization balance is
//   df/dt = (1-f) P + f(1-f) n_tot (S_i+S_CR) - f² n_tot α_r - f³ n_tot² κ_c.
// network_dfdt evaluates it; network_equilibrium_f finds the root in (0,1) by
// bisection. These replace the obsolete two-coefficient f_eq = S/(S+α): the
// multilevel collisional channel S_CR (Rydberg ladder) dominates direct
// Voronov S_i in the chromosphere, so the fixed point now depends on the full
// network (and on n_e, via the three-body f³ term).
// ----------------------------------------------------------------------------
static float network_dfdt(const Grid& grid, float T, float n_tot, float f, float P) {
    Vec Tv(1); Tv(0) = T;
    // Respect the same grid flags as apply_ionization_stage so the helper
    // stays consistent with the integrator's active network.
    const float Si   = grid.enable_direct_collisional_ionization
                     ? ionization_rate_S(grid, Tv)(0) : 0.0f;
    const float ar   = recombination_rate_alpha(grid, Tv)(0);
    const float Scr  = ionization_rate_S_CR(grid, Tv)(0);
    const float kc   = grid.enable_threebody_recombination
                     ? recombination_rate_kappa_c(grid, Tv)(0) : 0.0f;
    const float Stot = Si + Scr;
    // Group κ_c into the three-body product before squaring n_tot (n_tot²
    // alone overflows float32 for n_tot ≳ 1e19).
    return (1.0f - f) * P + f * (1.0f - f) * n_tot * Stot
         - f * f * n_tot * ar - (f * f * f) * n_tot * (n_tot * kc);
}
static float network_equilibrium_f(const Grid& grid, float T, float n_tot, float P) {
    // df/dt > 0 at f→0+ (when ionization can win) and < 0 at f=1, so the
    // physical equilibrium is bracketed by (ε, 1].
    float lo = 1.0e-12f, hi = 1.0f;
    for (int it = 0; it < 200; ++it) {
        const float mid = 0.5f * (lo + hi);
        if (network_dfdt(grid, T, n_tot, mid, P) > 0.0f) lo = mid; else hi = mid;
    }
    return 0.5f * (lo + hi);
}

// At chromospheric T_e (~6500 K) the direct Voronov S_i is exponentially
// suppressed, but the multilevel collisional channel S_CR is NOT — so the gas
// has a real collisional ionization route even without photoionization. Adding
// photoionization raises the equilibrium ionization fraction. Verify that
// Stage E relaxes to the full-network equilibrium in both cases, and that the
// photoionization-on equilibrium exceeds the photoionization-off one.
static void test_stage_e_photoionization_drives_low_T_equilibrium() {
    const float T = 6.5e3f;
    const float n_tot_target = 1.0e19f;
    const float f0 = 1.0e-3f;   // start very nearly neutral

    Grid grid_probe;
    grid_probe.init(1, 0.25f);

    // Full-network equilibria with photoionization on (P=1e-4) and off (P=0).
    const float eq_on  = network_equilibrium_f(grid_probe, T, n_tot_target, 1.0e-4f);
    const float eq_off = network_equilibrium_f(grid_probe, T, n_tot_target, 0.0f);
    // Photoionization is an extra ionization source, so it raises f_eq.
    EXPECT_TRUE(eq_on > eq_off);

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
        Vec Tv(1); Tv(0) = T;
        const float a_probe = recombination_rate_alpha(grid, Tv)(0);
        const float dt = 1.0f / std::max(grid.photoionization_rate,
                                         n_tot_target * a_probe);
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

    // Each run must lock onto its own full-network equilibrium, and the
    // photoionization-on equilibrium must exceed the photoionization-off one
    // (photoionization is an additional ionization source). Note f_off is NOT
    // ~0: the multilevel collisional channel S_CR keeps the gas partially
    // ionized even with the radiation field switched off.
    EXPECT_REL(f_on,  eq_on,  5e-2);
    EXPECT_REL(f_off, eq_off, 5e-2);
    EXPECT_TRUE(f_on > f_off);
}

// `model_c7_ic` populates `grid.photoionization_rate_i` with the Route B
// closure (writeup §3.1; docs/photoionization_c7_inversion_plan.md): below
// ~1500 km, P_phot is the C7-context equilibrium inversion of the *full Route B
// network*, so C7 is a self-consistent Stage-E fixed point there; above ~1500
// km it is blended (tanh window) to the frozen Chae (2021) FAL-C rate, the NEQ
// region where over-ionization is real physics and C7 is deliberately not
// pinned. The model column spans 1003–1989 km, so the blend bisects it.
// Verify: (a) the per-cell array is populated; (b) every rate is finite,
// positive, and bounded, and the deep-NEQ top cell matches the frozen FAL-C
// rate; (c) Stage E on the C7 IC is a NEAR-FIXED-POINT in the deep equilibrium
// region (h<1200 km, drift→0 by construction) and only a small bounded
// relaxation in the NEQ region (h>1700 km) — no runaway anywhere.
static void test_model_c7_photoionization_uses_route_b_closure() {
    Grid grid;
    grid.init(100, 0.25f);
    Vec xn = model_c7_ic(grid);

    // Reconstruct cell-center heights (model_c7_ic builds faces uniformly over
    // [1003, 2153] km — chromosphere + lower TR — with nF = ns + 5).
    const float h0 = 1003.0f, h1 = 2153.0f;
    const arma::uword nF = grid.ns + 5;
    Vec h_cell(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float hFi  = h0 + (h1 - h0) * (float)i       / (float)(nF - 1);
        const float hFi1 = h0 + (h1 - h0) * (float)(i + 1) / (float)(nF - 1);
        h_cell(i) = 0.5f * (hFi + hFi1);
    }

    // (a) per-cell rate populated with ns entries; helper returns it.
    EXPECT_TRUE(grid.photoionization_rate_i.n_elem == grid.ns);
    Vec P = photoionization_rate_P(grid);
    EXPECT_TRUE(P.n_elem == grid.ns);

    // (b) finite, positive, bounded; the deep-NEQ top cell ≈ frozen FAL-C rate.
    Vec P_falc = photoionization_rate_chae(h_cell);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_TRUE(std::isfinite(P(i)) && P(i) > 0.0f && P(i) <= 1.0e-2f);
    }
    EXPECT_REL(P(grid.ns - 1), P_falc(grid.ns - 1), 0.1);   // top ≈ pure Chae

    // (c) Stage E relaxation on the C7 IC, split by region.
    Vec prim_before = cons2prim(grid, xn);
    Vec rho_i_b = get_scalar(grid, prim_before, prim::RHO_I);
    Vec rho_n_b = get_scalar(grid, prim_before, prim::RHO_N);

    grid.enable_ionization = true;
    Vec prim_after = prim_before;
    apply_ionization_stage(grid, prim_after, 1.0f);

    Vec rho_i_a = get_scalar(grid, prim_after, prim::RHO_I);
    Vec rho_n_a = get_scalar(grid, prim_after, prim::RHO_N);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float f_b = rho_i_b(i) / (rho_i_b(i) + rho_n_b(i));
        const float f_a = rho_i_a(i) / (rho_i_a(i) + rho_n_a(i));
        EXPECT_TRUE(f_a > 0.0f && f_a < 1.0f);             // no runaway
        if (h_cell(i) < 1200.0f) {
            EXPECT_TRUE(std::fabs(f_a - f_b) < 5.0e-3f);   // near fixed point (inversion)
        } else if (h_cell(i) > 1700.0f) {
            EXPECT_TRUE(std::fabs(f_a - f_b) < 0.1f);      // bounded NEQ relaxation
        }
    }
}

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

    // Photoionization contributes NO χ_H drain (its energy is supplied by the
    // absorbed UV photon). The only χ_H exchange with the thermal pool is the
    // collisional balance: the two collisional ionization channels (direct S_i
    // and multilevel S_CR) drain χ_H, while three-body recombination returns it
    // (super-elastic). So Δ(e_i+e_n) must equal −(Γ_coll+Γ_mlvl−Γ_3b)·χ_H even
    // though a large photoionization rate is creating ions. The thermal/KE
    // exchange between species cancels in the sum (KE = 0 here, T_i = T_n).
    const float chi = grid.chi_H_J;
    Vec T_v(1); T_v(0) = T;
    const float S    = ionization_rate_S(grid, T_v)(0);
    const float Scr  = ionization_rate_S_CR(grid, T_v)(0);
    const float kc   = recombination_rate_kappa_c(grid, T_v)(0);
    const float n_i_n = rho_i_new(0) / grid.m_i;
    const float n_n_n = rho_n_new(0) / grid.m_i;
    const float Gamma_coll = dt * n_i_n * n_n_n * S;
    const float Gamma_mlvl = dt * n_i_n * n_n_n * Scr;
    const float Gamma_3b   = dt * n_i_n * n_i_n * (kc * n_i_n);
    const float expected_dE_sum = -(Gamma_coll + Gamma_mlvl - Gamma_3b) * chi;

    const float delta_p_sum = (p_i_new(0) + p_n_new(0)) - (p_i_old(0) + p_n_old(0));
    const float delta_e_sum = 1.5f * delta_p_sum;

    // The drain matches the collisional balance (independent of the large
    // photoionization rate). Absolute tolerance set well below |expected_dE_sum|.
    EXPECT_NEAR(delta_e_sum, expected_dE_sum, 1.0e-5f);
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

    // Stage E with a small dt: the multilevel S_CR drain is large at 1.5e4 K,
    // so dt is kept small enough that the χ_H drain stays well below the
    // available thermal energy (no pressure-floor clamp), keeping the clean
    // Δ(e_i+e_n) = −(drain) balance testable.
    const float dt = 1.0e-5f;
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
    const float Scr = ionization_rate_S_CR(grid, T_v)(0);
    const float kc  = recombination_rate_kappa_c(grid, T_v)(0);
    const float chi = grid.chi_H_J;

    // The thermal/KE exchange between species cancels in the SUM Δ(e_i+e_n),
    // leaving only the χ_H balance with the electron thermal pool: collisional
    // ionization (S_i and multilevel S_CR) drains χ_H, three-body recombination
    // returns it (super-elastic). Radiative recombination radiates χ_H away and
    // photoionization is photon-powered — neither touches the thermal pool. So
    //   Δ(e_i+e_n) = -(Γ_coll + Γ_mlvl - Γ_3b)·χ_H.
    const float Gamma_coll = dt * n_i_new * n_n_new * S;
    const float Gamma_mlvl = dt * n_i_new * n_n_new * Scr;
    const float Gamma_3b   = dt * n_i_new * n_i_new * (kc * n_i_new);
    const float delta_p_sum  = (p_i_new(0) + p_n_new(0)) - (p_i_old(0) + p_n_old(0));
    const float delta_e_sum  = 1.5f * delta_p_sum;
    const float expected_dE  = -(Gamma_coll + Gamma_mlvl - Gamma_3b) * chi;
    EXPECT_REL(delta_e_sum, expected_dE, 1e-2);
}

// Hold T_e fixed and evolve many steps with Stage E only (photoionization off).
// The full Route B network relaxes to its own ionization equilibrium, where
// ionization = recombination including the multilevel S_CR and three-body α_c
// channels (no longer the two-coefficient f = S_i/(S_i+α_r)).
static void test_stage_e_kinetic_equilibrium() {
    const float T = 3.0e4f;  // pick T so S_i/α_r is moderate (not tiny, not huge)
    const float ni0 = 5.0e17f, nn0 = 5.0e18f;
    const float n_tot = ni0 + nn0;
    Grid grid;
    Vec cons = setup_uniform(grid, 1, ni0, nn0, T, T);
    grid.enable_ionization    = true;
    grid.photoionization_rate = 0.0f;   // isolate the collisional/recombination network

    Vec T_v(1); T_v(0) = T;
    const float S = ionization_rate_S(grid, T_v)(0);
    const float a = recombination_rate_alpha(grid, T_v)(0);

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

    // Compare against the full-network equilibrium (ionization = recombination).
    const float f_eq_expected = network_equilibrium_f(grid, T, n_tot, 0.0f);
    EXPECT_REL(f_final, f_eq_expected, 1e-2);
}

// Fixed point: if the state starts at the full-network ionization equilibrium,
// Stage E must leave (ρ_i, ρ_n) essentially unchanged. This stresses the cubic
// solve at its algebraic root, where catastrophic cancellation would show up
// if the safeguarded-Newton logic were wrong.
static void test_stage_e_fixed_point_at_kinetic_equilibrium() {
    const float T = 8.0e3f;   // moderate f_eq (~0.3), so both species are well
                              // resolved in float32 (avoids the 1-f≈0 minority
                              // precision loss that T=3e4, f_eq→1 would cause)
    const float n_tot_target = 1.0e19f;
    // Enable the full network (including three-body κ_c) so the Stage-E cubic
    // is exercised at its algebraic root — the cubic solver is the target here.
    Grid g_probe;
    g_probe.init(1, 0.25f);
    g_probe.photoionization_rate                  = 0.0f;
    g_probe.enable_direct_collisional_ionization  = true;
    g_probe.enable_threebody_recombination        = true;
    const float f_eq = network_equilibrium_f(g_probe, T, n_tot_target, 0.0f);

    const float ni = f_eq * n_tot_target;
    const float nn = (1.0f - f_eq) * n_tot_target;
    Grid grid;
    Vec cons = setup_uniform(grid, 4, ni, nn, T, T);
    grid.enable_ionization                       = true;
    grid.photoionization_rate                    = 0.0f;   // match the P=0 equilibrium above
    grid.enable_direct_collisional_ionization    = true;
    grid.enable_threebody_recombination          = true;

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
                prim(arma::sub2ind(sz, cell, prim::P_E))   = n_i *        grid.k_b * T;  // p_e = ½ P_I
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
static void test_model_c7_bc_discrete_hse_inner_mach_capped_outer() {
    Grid grid;
    grid.init(100, 0.25f);
    Vec xn = model_c7_ic(grid);
    const auto sz = arma::size(grid.ns, num_of_eq);

    // IC-pinned inner reservoir (Dirichlet ρ) and decoded cell-0 pressure.
    const float rho_i_pinned = grid.inner_boundary0_i(cons::RHO_I);
    const float rho_n_pinned = grid.inner_boundary0_i(cons::RHO_N);
    const float phi_g_inner  = grid.phi_g_imh(0);
    auto ghost_p = [&](const Vec& ob, arma::uword E, arma::uword RHO, arma::uword MOM) {
        const float rho = ob(RHO);
        const float V   = ob(MOM) / rho;
        return 2.0f/3.0f * (ob(E) - 0.5f * rho * V * V - rho * phi_g_inner);
    };
    // Cell-0 pressures for the discrete-HSE comparison.
    const float rho_i_0 = xn(arma::sub2ind(sz, 0, cons::RHO_I));
    const float n_i_0   = rho_i_0 / grid.m_i;
    const float phi_g_c0 = 0.5f * (grid.phi_g_imh(0) + grid.phi_g_iph(0));
    const float V0       = xn(arma::sub2ind(sz, 0, cons::MOM_I)) / rho_i_0;
    const float p_i_cell0 = 2.0f/3.0f * (xn(arma::sub2ind(sz, 0, cons::E_I))
                            - 0.5f * rho_i_0 * V0 * V0 - rho_i_0 * phi_g_c0);

    model_c7_update_bc(grid, xn);

    // --- Inner: discrete-HSE V=U=0 reservoir, Dirichlet ρ ---
    EXPECT_REL(grid.inner_boundary0_i(cons::RHO_I), rho_i_pinned, 1e-5);
    EXPECT_REL(grid.inner_boundary0_i(cons::RHO_N), rho_n_pinned, 1e-5);
    EXPECT_REL(grid.inner_boundary1_i(cons::RHO_I), rho_i_pinned, 1e-5);
    EXPECT_TRUE(grid.inner_boundary0_i(cons::MOM_I) == 0.0f);   // V = 0
    EXPECT_TRUE(grid.inner_boundary0_i(cons::MOM_N) == 0.0f);   // U = 0
    // Hydrostatic ghost pressure: p_ghost = p_0 + ρ_0 g ds0 > p_0.
    const float p_i_ghost = ghost_p(grid.inner_boundary0_i, cons::E_I, cons::RHO_I, cons::MOM_I);
    EXPECT_TRUE(p_i_ghost > p_i_cell0);
    const float ds0 = grid.ds_i(0);
    EXPECT_REL(p_i_ghost, p_i_cell0 + rho_i_0 * grid.g * ds0, 1e-3);

    // --- Outer: Mach-capped outflow, U locked to V ---
    // Force a large outflow in cell ns-1; the ghost V must be capped at
    // 0.05 c_s(T_TR) and the neutral ghost U must equal V.
    const arma::uword L = grid.ns - 1;
    const float rho_i_L = xn(arma::sub2ind(sz, L, cons::RHO_I));
    const float rho_n_L = xn(arma::sub2ind(sz, L, cons::RHO_N));
    const float V_big   = 5.0e4f;   // 50 km/s, far above the cap
    xn(arma::sub2ind(sz, L, cons::MOM_I)) = rho_i_L * V_big;
    xn(arma::sub2ind(sz, L, cons::MOM_N)) = rho_n_L * V_big;

    model_c7_update_bc(grid, xn);

    const float T_TR  = 2.310e4f;   // table-top TR-base temperature
    const float c_s   = std::sqrt(2.0f * grid.gamma_mono * grid.k_b * T_TR / grid.m_i);
    const float V_cap = 0.05f * c_s;
    const float Vg = grid.outer_boundary0_i(cons::MOM_I) / grid.outer_boundary0_i(cons::RHO_I);
    const float Ug = grid.outer_boundary0_i(cons::MOM_N) / grid.outer_boundary0_i(cons::RHO_N);
    EXPECT_REL(Vg, V_cap, 1e-3);    // capped (positive outflow)
    EXPECT_REL(Ug, Vg, 1e-5);       // neutral locked to charge fluid
    // Outer ρ is Neumann (tracks cell ns-1).
    EXPECT_REL(grid.outer_boundary0_i(cons::RHO_I), rho_i_L, 1e-5);
}

// "New explanation" model_c7 upper BC, HYBRID form (tr_jump_bc=true,
// docs/gentle_evaporation_downflow.md): turns on well_balanced and imposes at the
// top face a hydrostatic ghost pressure + EOS ghost density with a CONTINUOUS
// ghost temperature (T_ghost0 = a·T_top, a=b=1 defaults — no jump). The coronal
// heat keeps entering via the imposed Neumann flux q(T) (NOT a Dirichlet T-jump,
// which over-conducts on the coarse grid), so impose_outer_heat_flux is left
// tracking enable_radiative_cooling rather than forced off.
static void test_model_c7_tr_jump_bc() {
    Grid grid;
    grid.init(100, 0.25f);
    Vec xn = model_c7_ic(grid, /*extended=*/false, /*tr_jump_bc=*/true);
    const auto sz = arma::size(grid.ns, num_of_eq);

    // IC turns on the hybrid BC + well_balanced; the imposed flux is untouched
    // (it tracks cooling — off here in the default-constructed grid).
    EXPECT_TRUE(grid.c7_tr_jump_bc);
    EXPECT_TRUE(grid.well_balanced);
    EXPECT_TRUE(grid.impose_outer_heat_flux == grid.enable_radiative_cooling);

    const arma::uword nl = grid.ns - 1;
    const float ds        = grid.ds_i(nl);
    const float phi_g_out = grid.phi_g_iph(nl);

    // Decode an interior cell's ion pressure / density / temperature (cons2prim).
    auto cell_piT = [&](arma::uword i, float& p_i, float& rho_i, float& T_i) {
        rho_i = xn(arma::sub2ind(sz, i, cons::RHO_I));
        const float V     = xn(arma::sub2ind(sz, i, cons::MOM_I)) / rho_i;
        const float E_i   = xn(arma::sub2ind(sz, i, cons::E_I));
        const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        p_i = 2.0f/3.0f * E_i - 1.0f/3.0f * rho_i * V * V - 2.0f/3.0f * rho_i * phi_g;
        T_i = p_i / (2.0f * (rho_i / grid.m_i) * grid.k_b);
    };
    float p_i_top, rho_i_top, T_top, p_i_2, rho_i_2, T_2;
    cell_piT(nl, p_i_top, rho_i_top, T_top);
    cell_piT(nl - 1, p_i_2, rho_i_2, T_2);

    model_c7_update_bc(grid, xn);

    // Decode the ghost ion pressure / temperature (V=0 expected at relaxation IC).
    auto ghost_pi_T = [&](const Vec& ob, float& p_i, float& T_i) {
        const float rho_i = ob(cons::RHO_I);
        const float V     = ob(cons::MOM_I) / rho_i;
        p_i = 2.0f/3.0f * (ob(cons::E_I) - 0.5f * rho_i * V * V - rho_i * phi_g_out);
        T_i = p_i / (2.0f * (rho_i / grid.m_i) * grid.k_b);   // p_i = 2 n_i k T
    };
    float p_i_g0, T_i_g0, p_i_g1, T_i_g1;
    ghost_pi_T(grid.outer_boundary0_i, p_i_g0, T_i_g0);
    ghost_pi_T(grid.outer_boundary1_i, p_i_g1, T_i_g1);

    // (2) CONTINUOUS ghost temperature (a=b=1 ⇒ both ghosts at the top-cell T).
    EXPECT_REL(T_i_g0, T_top, 1e-3);
    EXPECT_REL(T_i_g1, T_top, 1e-3);

    // (1) hydrostatic ghost pressure: (p_2 − p_ghost0)/(2Δs) = ρ_top g.
    EXPECT_REL(p_i_g0, p_i_2 - 2.0f * ds * rho_i_top * grid.g, 1e-3);

    // (3) EOS density at the imposed (p, T): n_i = p_i / (2 k T).
    const float n_i_g0 = grid.outer_boundary0_i(cons::RHO_I) / grid.m_i;
    EXPECT_REL(n_i_g0, p_i_g0 / (2.0f * grid.k_b * T_i_g0), 1e-3);

    // Electron energy packs T_e = T_i (ε_e = ¾ p_i).
    EXPECT_REL(grid.outer_boundary0_i(cons::E_E), 0.75f * p_i_g0, 1e-3);
}

// TRAC broadening factor ε(T): 1 outside [T_b, T_c), (T_c/T)^{5/2} inside, and
// the conservation property κ'=κε is constant (=κ(T_c)) so κ'Λ' = κΛ.
static void test_trac_broadening_conserves_kappa_lambda() {
    Grid grid;
    grid.init(6, 0.25f);
    grid.enable_trac   = true;
    grid.trac_T_chrom  = 2.0e4f;
    grid.trac_cutoff_T = 1.0e5f;
    Vec T(6);
    T(0) = 1.0e4f; T(1) = 2.0e4f; T(2) = 5.0e4f;
    T(3) = 1.0e5f; T(4) = 2.0e5f; T(5) = 8.0e3f;
    Vec eps = trac_broadening_factor(grid, T);

    EXPECT_REL(eps(0), 1.0f, 1e-4);   // below T_b
    EXPECT_REL(eps(5), 1.0f, 1e-4);   // below T_b
    EXPECT_REL(eps(3), 1.0f, 1e-4);   // T == T_c → not broadened
    EXPECT_REL(eps(4), 1.0f, 1e-4);   // above T_c
    EXPECT_REL(eps(1), std::pow(1.0e5f / 2.0e4f, 2.5f), 1e-3);  // inside
    EXPECT_REL(eps(2), std::pow(1.0e5f / 5.0e4f, 2.5f), 1e-3);

    // Spitzer-like κ ∝ T^{5/2}: broadened κ' = κ·ε must equal κ(T_c), constant.
    const float kappa_Tc = std::pow(1.0e5f, 2.5f);
    EXPECT_REL(std::pow(T(1), 2.5f) * eps(1), kappa_Tc, 1e-3);
    EXPECT_REL(std::pow(T(2), 2.5f) * eps(2), kappa_Tc, 1e-3);

    // Disabled → all ones (clean baseline).
    grid.enable_trac = false;
    Vec eps_off = trac_broadening_factor(grid, T);
    for (arma::uword i = 0; i < 6; ++i) EXPECT_REL(eps_off(i), 1.0f, 1e-5);
}

// TRAC adaptive cutoff: returns the floor for a resolved profile, rises above it
// for an under-resolved (steep) TR, bounded by 0.2 T_peak, and the per-step
// limiter caps the rate of change.
static void test_trac_cutoff_detection_and_limiter() {
    Grid grid;
    grid.init(10, 0.25f);
    grid.ds_i.fill(1.0e4f);
    grid.B_i.fill(1.0f); grid.B_imh.fill(1.0f); grid.B_iph.fill(1.0f);
    grid.dinvB_ds_i.zeros(); grid.phi_g_imh.zeros(); grid.phi_g_iph.zeros();
    grid.broadcast();
    grid.enable_trac  = true;
    grid.trac_T_chrom = 2.0e4f;

    auto make_prim = [&](const float* Tp) {
        Vec prim(grid.n_state, arma::fill::zeros);
        const float ni = 1.0e16f, nn = 1.0e10f;
        const auto sz = arma::size(10, num_of_eq);
        for (arma::uword i = 0; i < 10; ++i) {
            prim(arma::sub2ind(sz, i, prim::RHO_I)) = ni * grid.m_i;
            prim(arma::sub2ind(sz, i, prim::RHO_N)) = nn * grid.m_n;
            prim(arma::sub2ind(sz, i, prim::P_I))   = ni * 2.0f * grid.k_b * Tp[i];
            prim(arma::sub2ind(sz, i, prim::P_N))   = nn * grid.k_b * Tp[i];
        }
        return prim;
    };

    // Smooth, well-resolved profile → nothing under-resolved → T_c at floor.
    float smooth[10];
    for (int i = 0; i < 10; ++i) smooth[i] = 6.0e3f + 100.0f * i;
    grid.trac_cutoff_T = 2.0e4f;
    EXPECT_REL(compute_trac_cutoff_T(grid, make_prim(smooth)), grid.trac_T_chrom, 1e-4);

    // Steep TR jump (under-resolved). T_peak = 1.2e5 → cutoff bound 0.2·T_peak = 2.4e4.
    float steep[10] = {6.0e3f,6.0e3f,6.0e3f,6.0e3f,6.0e3f,6.0e3f,2.0e4f,6.0e4f,1.0e5f,1.2e5f};
    grid.trac_cutoff_T = 2.0e4f;
    const float Tc1 = compute_trac_cutoff_T(grid, make_prim(steep));
    EXPECT_TRUE(Tc1 >= 2.0e4f && Tc1 <= 2.0e4f * 1.0301f);   // limiter caps the rise
    for (int k = 0; k < 200; ++k)
        grid.trac_cutoff_T = compute_trac_cutoff_T(grid, make_prim(steep));
    EXPECT_TRUE(grid.trac_cutoff_T > 2.05e4f);                // TRAC engaged
    EXPECT_TRUE(grid.trac_cutoff_T <= 0.2f * 1.2e5f + 1.0f);  // bounded by 0.2 T_peak
}

// Flare beam heating (physics.hpp::beam_heating_rate): the volumetric heating
// is zero before onset, deposits only inside the height window during the
// flat-top, and its column integral equals the beam energy flux (∫φ ds = 1).
static void test_beam_heating_rate_profile() {
    Grid grid;
    grid.init(10, 0.25f);
    grid.ds_i.fill(1.0e4f);                 // 10 km cells; centers 1008,1018,..,1098 km
    grid.B_i.ones(); grid.B_imh.ones(); grid.B_iph.ones();
    grid.dinvB_ds_i.zeros(); grid.phi_g_imh.zeros(); grid.phi_g_iph.zeros();
    grid.broadcast();
    grid.enable_beam_heating = true;
    grid.beam_flux     = 5.0e7f;            // W/m^2
    grid.beam_t_on     = 2.0f;
    grid.beam_duration = 10.0f;
    grid.beam_ramp     = 1.0f;
    grid.beam_h_lo_km  = 1020.0f;
    grid.beam_h_hi_km  = 1070.0f;           // selects cell centers 1028..1068 (i=2..6)

    Vec n_i(10), n_n(10);
    n_i.fill(1.0e16f); n_n.fill(2.0e16f);

    // Before onset → no heating.
    grid.sim_time = 0.0f;
    EXPECT_TRUE(arma::max(beam_heating_rate(grid, n_i, n_n)) == 0.0f);

    // Flat-top (t well inside [t_on+ramp, t_off-ramp]) → g = 1.
    grid.sim_time = 6.0f;
    Vec Q = beam_heating_rate(grid, n_i, n_n);
    EXPECT_TRUE(Q(0) == 0.0f && Q(1) == 0.0f);          // below window
    EXPECT_TRUE(Q(7) == 0.0f && Q(9) == 0.0f);          // above window
    for (int i = 2; i <= 6; ++i) EXPECT_TRUE(Q(i) > 0.0f);   // inside window
    // Column-integrated heating equals the imposed flux at flat-top.
    EXPECT_REL(arma::dot(Q, grid.ds_i), grid.beam_flux, 1e-3);

    // Disabled → identically zero (clean baseline for other scenarios).
    grid.enable_beam_heating = false;
    EXPECT_TRUE(arma::max(beam_heating_rate(grid, n_i, n_n)) == 0.0f);
}

// Flare beam heating stage shares the deposited energy by heat capacity so both
// fluids gain the SAME ΔT (the weakly-ionized chromosphere would otherwise spike
// T_e). Verifies ΔT_i == ΔT_n and that ΔT matches Δt·Q/(C_i + C_n).
static void test_beam_heating_partitions_by_heat_capacity() {
    Grid grid;
    grid.init(4, 0.25f);
    grid.ds_i.fill(1.0e4f);
    grid.B_i.ones(); grid.B_imh.ones(); grid.B_iph.ones();
    grid.dinvB_ds_i.zeros(); grid.phi_g_imh.zeros(); grid.phi_g_iph.zeros();
    grid.broadcast();
    grid.enable_beam_heating = true;
    grid.enable_trac   = false;
    grid.beam_flux     = 1.0e7f;
    grid.beam_t_on     = 0.0f;
    grid.beam_duration = 10.0f;
    grid.beam_ramp     = 1.0e-6f;           // ~instant flat top
    grid.beam_h_lo_km  = 1000.0f;
    grid.beam_h_hi_km   = 1100.0f;          // whole little grid in window
    grid.sim_time      = 5.0f;

    const float ni = 1.0e16f, nn = 5.0e16f, T0 = 6.0e3f;
    Vec prim(grid.n_state, arma::fill::zeros);
    const auto sz = arma::size(4, num_of_eq);
    for (arma::uword i = 0; i < 4; ++i) {
        prim(arma::sub2ind(sz, i, prim::RHO_I)) = ni * grid.m_i;
        prim(arma::sub2ind(sz, i, prim::RHO_N)) = nn * grid.m_n;
        prim(arma::sub2ind(sz, i, prim::P_I))   = 2.0f * ni * grid.k_b * T0;
        prim(arma::sub2ind(sz, i, prim::P_N))   =        nn * grid.k_b * T0;
    }

    const float dt = 0.01f;
    Vec Q = beam_heating_rate(grid, Vec(4, arma::fill::value(ni)),
                                    Vec(4, arma::fill::value(nn)));
    apply_beam_heating_stage(grid, prim, dt);

    for (arma::uword i = 0; i < 4; ++i) {
        const float pi = prim(arma::sub2ind(sz, i, prim::P_I));
        const float pn = prim(arma::sub2ind(sz, i, prim::P_N));
        const float Ti = pi / (2.0f * ni * grid.k_b);
        const float Tn = pn /        (nn * grid.k_b);
        const float dTi = Ti - T0, dTn = Tn - T0;
        EXPECT_REL(dTi, dTn, 1e-3);                              // equal ΔT both fluids
        const float C = (3.0f * ni + 1.5f * nn) * grid.k_b;     // total heat capacity
        EXPECT_REL(dTi, dt * Q(i) / C, 1e-3);                   // = Δt Q / (C_i+C_n)
    }
}

// Ambient coronal heating: footpoint-anchored exponential H(s) = E0·exp(−d/s_H),
// one-sided from the inner footpoint (open / half loop) or two-sided (full loop),
// plus the optional Phase-3 time ramp. Verifies the spatial shape, the full-loop
// symmetry, the ramp, and the disabled = identically-zero baseline.
static void test_coronal_heating_rate_profile() {
    Grid grid;
    grid.init(10, 0.25f);
    grid.ds_i.fill(1.0e4f);                 // 10 km cells; arc centers 5,15,..,95 km
    grid.B_i.ones(); grid.B_imh.ones(); grid.B_iph.ones();
    grid.dinvB_ds_i.zeros(); grid.phi_g_imh.zeros(); grid.phi_g_iph.zeros();
    grid.broadcast();

    // Disabled → identically zero (clean baseline for steady scenarios / tests).
    EXPECT_TRUE(arma::max(coronal_heating_rate(grid)) == 0.0f);

    grid.enable_coronal_heating = true;
    grid.coronal_heat_E0 = 1.0e-3f;         // W/m^3
    grid.coronal_heat_sH = 2.0e4f;          // 20 km scale length
    grid.coronal_heat_s0 = 0.0f;

    // One-sided: max at the inner footpoint, monotonically decaying upward.
    Vec H = coronal_heating_rate(grid);
    EXPECT_REL(H(0), 1.0e-3f * std::exp(-5.0e3f / 2.0e4f), 1e-4);
    for (arma::uword i = 1; i < 10; ++i) EXPECT_TRUE(H(i) < H(i - 1));

    // Two-sided (full loop): symmetric about the apex, minimum at the center.
    grid.coronal_heat_two_sided = true;
    Vec H2 = coronal_heating_rate(grid);
    EXPECT_REL(H2(0), H2(9), 1e-4);         // both footpoints equal
    EXPECT_REL(H2(4), H2(5), 1e-4);         // cells straddling the apex (the minimum)
    EXPECT_TRUE(H2(4) < H2(0) && H2(5) < H2(9));

    // Phase-3 ramp: amp = 1 before t_on, → enhance after t_on + ramp.
    grid.coronal_heat_two_sided = false;
    grid.coronal_heat_enhance = 3.0f;
    grid.coronal_heat_t_on    = 2.0f;
    grid.coronal_heat_ramp    = 1.0f;
    grid.sim_time = 0.0f;
    EXPECT_REL(coronal_heating_rate(grid)(0), H(0), 1e-4);          // pre-ramp = steady
    grid.sim_time = 5.0f;
    EXPECT_REL(coronal_heating_rate(grid)(0), 3.0f * H(0), 1e-4);   // fully ramped ×3
}

// Ambient coronal heating stage shares the deposited energy by heat capacity in
// the single-T baseline so both fluids gain the SAME ΔT = Δt H / (C_i + C_n),
// mirroring the beam-heating partition (apply_coronal_heating_stage).
static void test_coronal_heating_partitions_by_heat_capacity() {
    Grid grid;
    grid.init(4, 0.25f);
    grid.ds_i.fill(1.0e4f);
    grid.B_i.ones(); grid.B_imh.ones(); grid.B_iph.ones();
    grid.dinvB_ds_i.zeros(); grid.phi_g_imh.zeros(); grid.phi_g_iph.zeros();
    grid.broadcast();
    grid.enable_coronal_heating = true;
    grid.enable_trac     = false;
    grid.coronal_heat_E0 = 1.0e-3f;
    grid.coronal_heat_sH = 1.0e9f;          // ≫ grid ⇒ ~uniform heating

    const float ni = 1.0e16f, nn = 5.0e16f, T0 = 6.0e3f;
    Vec prim(grid.n_state, arma::fill::zeros);
    const auto sz = arma::size(4, num_of_eq);
    for (arma::uword i = 0; i < 4; ++i) {
        prim(arma::sub2ind(sz, i, prim::RHO_I)) = ni * grid.m_i;
        prim(arma::sub2ind(sz, i, prim::RHO_N)) = nn * grid.m_n;
        prim(arma::sub2ind(sz, i, prim::P_I))   = 2.0f * ni * grid.k_b * T0;
        prim(arma::sub2ind(sz, i, prim::P_N))   =        nn * grid.k_b * T0;
    }

    const float dt = 1.0f;
    Vec Q = coronal_heating_rate(grid);
    apply_coronal_heating_stage(grid, prim, dt);

    for (arma::uword i = 0; i < 4; ++i) {
        const float pi = prim(arma::sub2ind(sz, i, prim::P_I));
        const float pn = prim(arma::sub2ind(sz, i, prim::P_N));
        const float Ti = pi / (2.0f * ni * grid.k_b);
        const float Tn = pn /        (nn * grid.k_b);
        const float dTi = Ti - T0, dTn = Tn - T0;
        EXPECT_REL(dTi, dTn, 1e-3);                             // equal ΔT both fluids
        const float C = (3.0f * ni + 1.5f * nn) * grid.k_b;     // total heat capacity
        EXPECT_REL(dTi, dt * Q(i) / C, 1e-3);                   // = Δt H / (C_i+C_n)
    }
}

// =========================================================================
// Static local-refinement mesh builder
// =========================================================================

static void test_mesh_disabled_is_uniform() {
    // factor 1 ⇒ refinement off ⇒ exact k·coarse_ds grid.
    RefineParams rp;  // defaults: factor 1
    const double L = 1000.0e3;  // 1000 km
    const arma::uword nc = 250;
    const std::vector<double> f = build_static_mesh_faces(L, nc, rp);
    EXPECT_TRUE(f.size() == nc + 1);
    const double ds = L / static_cast<double>(nc);
    for (arma::uword k = 0; k <= nc; ++k)
        EXPECT_NEAR(f[k], ds * static_cast<double>(k), 1.0e-3);
    EXPECT_NEAR(f.back(), L, 1.0e-6);
}

static void test_mesh_refined_diagnostic_profile() {
    // Documented diagnostic profile: 4×, 0–700 km, 100 km transition.
    RefineParams rp; rp.factor = 4.0; rp.s_lo_km = 0.0; rp.s_hi_km = 700.0; rp.transition_km = 100.0;
    const double L = 1303.0e3;
    const arma::uword nc = 2000;
    const std::vector<double> f = build_static_mesh_faces(L, nc, rp);
    const double coarse_ds = L / static_cast<double>(nc);

    // Exact domain endpoints.
    EXPECT_NEAR(f.front(), 0.0, 1.0e-9);
    EXPECT_NEAR(f.back(),  L,   1.0e-6);

    // Exact 700 km and 800 km faces exist.
    double d700 = 1.0e30, d800 = 1.0e30;
    for (double x : f) { d700 = std::min(d700, std::fabs(x - 700.0e3));
                         d800 = std::min(d800, std::fabs(x - 800.0e3)); }
    EXPECT_NEAR(d700, 0.0, 1.0e-3);
    EXPECT_NEAR(d800, 0.0, 1.0e-3);

    // Strictly positive widths; adjacent ratio ≤ 1.1 everywhere.
    double max_ratio = 1.0;
    bool all_pos = true;
    for (arma::uword i = 0; i + 1 < f.size(); ++i) {
        const double w = f[i + 1] - f[i];
        if (!(w > 0.0)) all_pos = false;
        if (i > 0) {
            const double wp = f[i] - f[i - 1];
            const double r  = std::max(w / wp, wp / w);
            max_ratio = std::max(max_ratio, r);
        }
    }
    EXPECT_TRUE(all_pos);
    EXPECT_TRUE(max_ratio <= 1.1 + 1.0e-6);

    // ×4 fine spacing in the refined region (first cell) and retained outer
    // coarse spacing (last cell).
    const double fine_w = f[1] - f[0];
    const double outer_w = f.back() - f[f.size() - 2];
    EXPECT_REL(fine_w, coarse_ds / 4.0, 0.02);
    EXPECT_REL(outer_w, coarse_ds, 0.02);
}

static void test_refine_params_env_aliases() {
    // ISO_REFINE_* override GRID_REFINE_*; both parsed. Restore env after.
    setenv("GRID_REFINE_FACTOR", "2", 1);
    setenv("GRID_REFINE_S_HI_KM", "500", 1);
    RefineParams g = refine_params_from_env(/*iso_alias=*/false);
    EXPECT_NEAR(g.factor, 2.0, 1e-9);
    EXPECT_NEAR(g.s_hi_km, 500.0, 1e-9);
    setenv("ISO_REFINE_FACTOR", "4", 1);
    RefineParams a = refine_params_from_env(/*iso_alias=*/true);
    EXPECT_NEAR(a.factor, 4.0, 1e-9);       // ISO_ override
    EXPECT_NEAR(a.s_hi_km, 500.0, 1e-9);    // falls through to GRID_
    unsetenv("GRID_REFINE_FACTOR"); unsetenv("GRID_REFINE_S_HI_KM"); unsetenv("ISO_REFINE_FACTOR");
}

// =========================================================================
// Irregular-grid metric caches + manufactured operator tests
// =========================================================================

// Build an ns-cell state on a non-uniform cell-width array. B = 1; gravity via a
// constant g populating φ_g from arc length. Interior ρ,T uniform; ghosts pinned
// to the matching end cell so a constant state is a discrete fixed point.
static Vec setup_irregular(Grid& grid, const std::vector<float>& ds,
                           float ni_val, float nn_val, float Ti_val, float Tn_val) {
    grid.init(static_cast<arma::uword>(ds.size()), 0.25f);
    grid.B_imh.fill(1.0f); grid.B_iph.fill(1.0f); grid.B_i.fill(1.0f);
    grid.dinvB_ds_i.zeros();
    float s = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        grid.ds_i(i) = ds[i];
        grid.phi_g_imh(i) = 0.0f;
        grid.phi_g_iph(i) = 0.0f;
        s += ds[i];
    }
    grid.broadcast();

    Vec prim(grid.n_state, arma::fill::zeros);
    const auto sz = arma::size(grid.ns, num_of_eq);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        prim(arma::sub2ind(sz, i, prim::RHO_I)) = ni_val * grid.m_i;
        prim(arma::sub2ind(sz, i, prim::RHO_N)) = nn_val * grid.m_n;
        prim(arma::sub2ind(sz, i, prim::P_I))   = ni_val * 2.0f * grid.k_b * Ti_val;
        prim(arma::sub2ind(sz, i, prim::P_N))   = nn_val * grid.k_b * Tn_val;
        prim(arma::sub2ind(sz, i, prim::P_E))   = ni_val * grid.k_b * Ti_val;
    }
    Vec cons = prim2cons(grid, prim);
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        grid.inner_boundary0_i(k) = cons(arma::sub2ind(sz, 0, k));
        grid.inner_boundary1_i(k) = cons(arma::sub2ind(sz, 0, k));
        grid.outer_boundary0_i(k) = cons(arma::sub2ind(sz, grid.ns - 1, k));
        grid.outer_boundary1_i(k) = cons(arma::sub2ind(sz, grid.ns - 1, k));
    }
    return cons;
}

// A graded, strictly non-uniform width list (fine base → coarse top).
static std::vector<float> graded_widths() {
    std::vector<float> ds;
    float w = 100.0f;
    for (int i = 0; i < 20; ++i) { ds.push_back(w); w *= 1.08f; }
    return ds;
}

static void test_metric_caches_center_to_center() {
    Grid grid;
    const std::vector<float> ds = graded_widths();
    setup_irregular(grid, ds, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    EXPECT_TRUE(!grid.uniform_mesh);
    // ds_iph_i / ds_imh_i = mean of adjacent widths, mirrored at the ends.
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float dsp1 = (i + 1 < grid.ns) ? ds[i + 1] : ds[i];
        const float dsm1 = (i > 0) ? ds[i - 1] : ds[i];
        EXPECT_REL(grid.ds_iph_i(i), 0.5f * (ds[i] + dsp1), 1e-5);
        EXPECT_REL(grid.ds_imh_i(i), 0.5f * (dsm1 + ds[i]), 1e-5);
    }
    // Cell centers accumulate half-widths.
    float s = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_REL(grid.s_i(i), s + 0.5f * ds[i], 1e-5);
        s += ds[i];
    }
    EXPECT_REL(grid.s_face(grid.ns), s, 1e-5);
}

static void test_irregular_uniform_state_preserved() {
    // A constant (V=0, uniform ρ,p, gravity-free) state has zero explicit RHS on
    // a non-uniform mesh — the metric MUSCL must reproduce constants exactly.
    Grid grid;
    Vec cons = setup_irregular(grid, graded_widths(), 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    grid.dt_state.fill(1.0e-6f);   // small predictor step
    Vec R = rhs_explicit_state(grid, cons);
    EXPECT_TRUE(!R.has_nan());
    EXPECT_NEAR(arma::max(arma::abs(R)), 0.0, 1e-3);
}

static void test_irregular_linear_pressure_gradient() {
    // Linear charged pressure p_i(s) = p0 + b·s, V=0, gravity-free, constant ρ.
    // The metric MUSCL reconstructs each face exactly, so the momentum RHS in the
    // interior equals the analytic −dp/ds = −b. (A uniform-spacing reconstruction
    // would recover the WRONG gradient on this irregular grid.)
    Grid grid;
    const std::vector<float> ds = graded_widths();
    setup_irregular(grid, ds, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    const auto sz = arma::size(grid.ns, num_of_eq);
    const float p0 = 1.0f, b = 3.0e-6f;    // Pa, Pa/m
    const float rho_i = 2.0e17f * grid.m_i, rho_n = 1.0e19f * grid.m_n;
    Vec prim(grid.n_state, arma::fill::zeros);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float pI = p0 + b * grid.s_i(i);
        prim(arma::sub2ind(sz, i, prim::RHO_I)) = rho_i;
        prim(arma::sub2ind(sz, i, prim::RHO_N)) = rho_n;
        prim(arma::sub2ind(sz, i, prim::P_I))   = pI;
        prim(arma::sub2ind(sz, i, prim::P_N))   = 0.5f * pI;
        prim(arma::sub2ind(sz, i, prim::P_E))   = 0.5f * pI;
    }
    Vec cons = prim2cons(grid, prim);
    // Linearly-continued ghosts so the boundary stencil is consistent too.
    auto ghost = [&](float s_ghost) {
        Vec pr(grid.n_state, arma::fill::zeros);
        const float pI = p0 + b * s_ghost;
        pr(arma::sub2ind(sz, 0, prim::RHO_I)) = rho_i; pr(arma::sub2ind(sz, 0, prim::RHO_N)) = rho_n;
        pr(arma::sub2ind(sz, 0, prim::P_I)) = pI; pr(arma::sub2ind(sz, 0, prim::P_N)) = 0.5f * pI;
        pr(arma::sub2ind(sz, 0, prim::P_E)) = 0.5f * pI;
        Vec c = prim2cons(grid, pr);
        Vec g(num_of_eq); for (arma::uword k = 0; k < num_of_eq; ++k) g(k) = c(arma::sub2ind(sz, 0, k));
        return g;
    };
    grid.inner_boundary0_i = ghost(grid.s_i(0) - grid.ds_imh_i(0));
    grid.inner_boundary1_i = ghost(grid.s_i(0) - 2.0f * grid.ds_imh_i(0));
    grid.outer_boundary0_i = ghost(grid.s_i(grid.ns - 1) + grid.ds_iph_i(grid.ns - 1));
    grid.outer_boundary1_i = ghost(grid.s_i(grid.ns - 1) + 2.0f * grid.ds_iph_i(grid.ns - 1));
    grid.dt_state.fill(1.0e-9f);   // suppress predictor O(dt²)
    Vec R = rhs_explicit_state(grid, cons);
    // Interior cells only (fully-interior stencil).
    for (arma::uword i = 3; i + 3 < grid.ns; ++i)
        EXPECT_REL(R(arma::sub2ind(sz, i, cons::MOM_I)), -b, 5e-3);
}

static void test_irregular_cfl_selection() {
    Grid grid;
    const std::vector<float> ds = graded_widths();
    Vec cons = setup_irregular(grid, ds, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    Vec dt   = cal_dt_i(grid, cons);
    Vec maxv = cal_max_v_i(grid, cons);
    Vec ratio(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) ratio(i) = grid.ds_i(i) / maxv(i);
    const float expected = grid.CFL * arma::min(ratio);
    EXPECT_REL(dt(0), expected, 1e-3);
    for (arma::uword i = 0; i < grid.ns; ++i)
        EXPECT_TRUE(dt(i) * maxv(i) / grid.ds_i(i) <= grid.CFL + 1e-6f);
}

static void test_irregular_conduction_conserves_energy() {
    // Non-uniform mesh, V=U=0, T_i=T_n a linear ramp, matched (flat) ghosts so
    // there is no boundary flux. The conservative conduction term then telescopes
    // to zero: Σ (C_i + C_n)·ds_i/B_i ≈ 0. This exercises the width-weighted
    // series-resistance face conductivity and the center-to-center distances.
    Grid grid;
    const std::vector<float> ds = graded_widths();
    setup_irregular(grid, ds, 2.0e17f, 1.0e19f, 1.0e4f, 1.0e4f);
    const auto sz = arma::size(grid.ns, num_of_eq);
    const float ni = 2.0e17f, nn = 1.0e19f;
    Vec prim(grid.n_state, arma::fill::zeros);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float T = 1.0e4f + 200.0f * static_cast<float>(i);   // linear ramp
        prim(arma::sub2ind(sz, i, prim::RHO_I)) = ni * grid.m_i;
        prim(arma::sub2ind(sz, i, prim::RHO_N)) = nn * grid.m_n;
        prim(arma::sub2ind(sz, i, prim::P_I))   = ni * 2.0f * grid.k_b * T;
        prim(arma::sub2ind(sz, i, prim::P_N))   = nn * grid.k_b * T;
        prim(arma::sub2ind(sz, i, prim::P_E))   = ni * grid.k_b * T;
    }
    Vec cons = prim2cons(grid, prim);
    for (arma::uword k = 0; k < num_of_eq; ++k) {   // flat (matched) ghosts ⇒ zero boundary flux
        grid.inner_boundary0_i(k) = cons(arma::sub2ind(sz, 0, k));
        grid.inner_boundary1_i(k) = cons(arma::sub2ind(sz, 0, k));
        grid.outer_boundary0_i(k) = cons(arma::sub2ind(sz, grid.ns - 1, k));
        grid.outer_boundary1_i(k) = cons(arma::sub2ind(sz, grid.ns - 1, k));
    }
    Vec RI = rhs_implicit_state(grid, cons);
    float net = 0.0f, scale = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float cE = RI(arma::sub2ind(sz, i, cons::E_I)) + RI(arma::sub2ind(sz, i, cons::E_N));
        net   += cE * grid.ds_i(i) / grid.B_i(i);
        scale += std::fabs(cE) * grid.ds_i(i) / grid.B_i(i);
    }
    EXPECT_TRUE(scale > 0.0f);                        // conduction actually active
    EXPECT_TRUE(std::fabs(net) <= 1e-4f * scale);     // conserved to round-off
}

static void test_irregular_boundary_hse() {
    // Boundary/interior HSE consistency on a REFINED (non-uniform) mesh: model_column
    // always uses equilibrium-reference well-balancing, so the hydrostatic C7 IC is an
    // EXACT discrete steady state and one full step leaves V = 0 to round-off — even at
    // the fine base cells whose HSE ghost is built with the actual center-to-ghost
    // distance ds_i(0). A mis-set ghost weight or wrong metric reconstruction would
    // break the fixed point.
    setenv("ISO_REFINE_FACTOR", "4", 1);
    setenv("ISO_REFINE_S_HI_KM", "700", 1);
    setenv("ISO_REFINE_TRANSITION_KM", "100", 1);
    Grid grid; grid.init(200, 0.25f);
    Vec xn = model_column_ic(grid);
    EXPECT_TRUE(!grid.uniform_mesh);
    Vec dt = cal_dt_i(grid, xn);
    model_column_update_bc(grid, xn);
    Vec out = advance_Euler_state(grid, xn, dt);
    EXPECT_TRUE(!out.has_nan());
    const auto sz = arma::size(grid.ns, num_of_eq);
    float vmax = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float v = out(arma::sub2ind(sz, i, cons::MOM_I)) /
                        out(arma::sub2ind(sz, i, cons::RHO_I));
        vmax = std::max(vmax, std::fabs(v));
    }
    EXPECT_TRUE(vmax < 1.0e-2f);   // exact V=0 fixed point on the irregular mesh
    unsetenv("ISO_REFINE_FACTOR"); unsetenv("ISO_REFINE_S_HI_KM"); unsetenv("ISO_REFINE_TRANSITION_KM");
}

static void test_uniform_regression_refinement_unset() {
    // With refinement UNSET, model_column builds a uniform grid and the step is
    // bitwise-reproducible (deterministic) — the regression guard that the static-mesh
    // machinery is a no-op on the uniform path (the best-WB reconstruction is always on).
    unsetenv("ISO_REFINE_FACTOR"); unsetenv("GRID_REFINE_FACTOR");
    auto run_once = [&]() {
        Grid grid; grid.init(200, 0.25f);
        Vec xn = model_column_ic(grid);
        Vec dt = cal_dt_i(grid, xn);
        model_column_update_bc(grid, xn);
        Vec out = advance_Euler_state(grid, xn, dt);
        return std::make_pair(grid.uniform_mesh, out);
    };
    auto a = run_once();
    auto b = run_once();
    EXPECT_TRUE(a.first);                                  // uniform_mesh detected
    EXPECT_TRUE(!a.second.has_nan());
    EXPECT_TRUE(a.second.n_elem == b.second.n_elem);
    EXPECT_TRUE(arma::approx_equal(a.second, b.second, "absdiff", 0.0));  // bitwise identical
}

static void test_refined_mesh_ic_runs() {
    // End-to-end: model_column with the refinement profile builds a larger,
    // non-uniform grid (finer at the base), and a full semi-implicit step stays
    // finite with velocities small (HSE-consistent refined boundary).
    setenv("ISO_REFINE_FACTOR", "4", 1);
    setenv("ISO_REFINE_S_LO_KM", "0", 1);
    setenv("ISO_REFINE_S_HI_KM", "700", 1);
    setenv("ISO_REFINE_TRANSITION_KM", "100", 1);
    Grid grid; grid.init(200, 0.25f);
    Vec xn = model_column_ic(grid);
    EXPECT_TRUE(grid.ns > 200);            // refinement added lower-domain cells
    EXPECT_TRUE(!grid.uniform_mesh);
    EXPECT_TRUE(grid.ds_i.min() > 0.0f);
    EXPECT_REL(arma::min(grid.ds_i), arma::max(grid.ds_i) / 4.0f, 0.05);  // fine ≈ coarse/4
    Vec dt = cal_dt_i(grid, xn);
    model_column_update_bc(grid, xn);
    Vec out = advance_Euler_state(grid, xn, dt);
    EXPECT_TRUE(!out.has_nan());
    const auto sz = arma::size(grid.ns, num_of_eq);
    float vmax = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float rho_i = out(arma::sub2ind(sz, i, cons::RHO_I));
        const float mom_i = out(arma::sub2ind(sz, i, cons::MOM_I));
        vmax = std::max(vmax, std::fabs(mom_i / rho_i));
    }
    EXPECT_TRUE(vmax < 1.0e3f);             // no boundary blowup on the first step
    unsetenv("ISO_REFINE_FACTOR");
    unsetenv("ISO_REFINE_S_LO_KM"); unsetenv("ISO_REFINE_S_HI_KM"); unsetenv("ISO_REFINE_TRANSITION_KM");
}

// ----------------------------------------------------------------------------

int main() {
    std::cout << "===== Chromosphere test suite =====\n\n";

    RUN(test_scalar_to_get_scalar_inverse);
    RUN(test_ip1_im1_interior_shift);
    RUN(test_flux_lim_is_minmod);
    RUN(test_flux_lim_mc3);
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
    RUN(test_model_c7_photoionization_uses_route_b_closure);
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
    RUN(test_model_c7_bc_discrete_hse_inner_mach_capped_outer);
    RUN(test_model_c7_tr_jump_bc);
    RUN(test_trac_broadening_conserves_kappa_lambda);
    RUN(test_trac_cutoff_detection_and_limiter);
    RUN(test_beam_heating_rate_profile);
    RUN(test_beam_heating_partitions_by_heat_capacity);
    RUN(test_coronal_heating_rate_profile);
    RUN(test_coronal_heating_partitions_by_heat_capacity);

    // Static local refinement: mesh builder, metric caches, irregular-grid
    // manufactured operators, and the refinement-unset regression.
    RUN(test_mesh_disabled_is_uniform);
    RUN(test_mesh_refined_diagnostic_profile);
    RUN(test_refine_params_env_aliases);
    RUN(test_metric_caches_center_to_center);
    RUN(test_irregular_uniform_state_preserved);
    RUN(test_irregular_linear_pressure_gradient);
    RUN(test_irregular_cfl_selection);
    RUN(test_irregular_conduction_conserves_energy);
    RUN(test_irregular_boundary_hse);
    RUN(test_uniform_regression_refinement_unset);
    RUN(test_refined_mesh_ic_runs);

    std::cout << "\n===== Summary =====\n";
    std::cout << "Passed: " << g_pass << "\n";
    std::cout << "Failed: " << g_fail << "\n";
    return g_fail == 0 ? 0 : 1;
}
