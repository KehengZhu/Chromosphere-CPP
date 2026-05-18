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
#include "../scenarios/model_c7.hpp"

#include <armadillo>
#include <cmath>
#include <iostream>

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

    std::cout << "\n===== Summary =====\n";
    std::cout << "Passed: " << g_pass << "\n";
    std::cout << "Failed: " << g_fail << "\n";
    return g_fail == 0 ? 0 : 1;
}
