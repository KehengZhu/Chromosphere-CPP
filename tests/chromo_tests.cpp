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
//   - Uniform, motionless, gravity-free state is a fixed point of advance_Euler_ii

#include "../chromosphere.hpp"
#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>

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
static Vec setup_uniform(uword n_cells, float ni_val, float nn_val,
                         float Ti_val, float Tn_val) {
    chromo_init(n_cells, 0.25f);

    ds_i.fill(1.0f);
    B_imh.fill(1.0f);
    B_iph.fill(1.0f);
    B_i.fill(1.0f);
    dinvB_ds_i.zeros();
    gPotential_imh.zeros();
    gPotential_iph.zeros();

    // Broadcast 1D quantities into the (ns*num_of_eq) layout
    B_iimh.zeros(num_of_elem);
    B_iiph.zeros(num_of_elem);
    B_ii.zeros(num_of_elem);
    ds_ii.zeros(num_of_elem);
    dinvB_ds_ii.zeros(num_of_elem);
    for (uword k = 0; k < num_of_eq; ++k) {
        B_iimh      += scalar_to(B_imh,      k);
        B_iiph      += scalar_to(B_iph,      k);
        B_ii        += scalar_to(B_i,        k);
        ds_ii       += scalar_to(ds_i,       k);
        dinvB_ds_ii += scalar_to(dinvB_ds_i, k);
    }

    // Primitive state, then convert to conserved.
    Vec prim(num_of_elem, fill::zeros);
    for (uword i = 0; i < n_cells; ++i) {
        prim(sub2ind(size(n_cells, num_of_eq), i, PNI)) = ni_val * m_i;
        prim(sub2ind(size(n_cells, num_of_eq), i, PNN)) = nn_val * m_n;
        prim(sub2ind(size(n_cells, num_of_eq), i, PV))  = 0.0f;
        prim(sub2ind(size(n_cells, num_of_eq), i, PU))  = 0.0f;
        prim(sub2ind(size(n_cells, num_of_eq), i, PPI)) = ni_val * 2.0f * k_b * Ti_val; // writeup eq 38
        prim(sub2ind(size(n_cells, num_of_eq), i, PPN)) = nn_val * k_b * Tn_val;
    }
    Vec cons = prim2cons(prim);

    // Pin inner/outer ghosts to the (uniform) interior values.
    for (uword k = 0; k < num_of_eq; ++k) {
        const float val = cons(sub2ind(size(n_cells, num_of_eq), 0, k));
        inner_boundary0_i(k) = val;
        inner_boundary1_i(k) = val;
        outer_boundary0_i(k) = val;
        outer_boundary1_i(k) = val;
    }
    return cons;
}


// =========================================================================
// Helper-function invariants
// =========================================================================

static void test_scalar_to_get_scalar_inverse() {
    chromo_init(8, 0.25f);
    Vec v(ns);
    for (uword i = 0; i < ns; ++i) v(i) = (float)i * 1.5f + 0.3f;

    for (uword k = 0; k < num_of_eq; ++k) {
        Vec back = get_scalar(scalar_to(v, k), k);
        for (uword i = 0; i < ns; ++i) EXPECT_NEAR(back(i), v(i), 1e-6);
    }
}

static void test_ip1_im1_interior_shift() {
    chromo_init(10, 0.25f);
    inner_boundary0_i.zeros(num_of_eq);
    outer_boundary0_i.zeros(num_of_eq);

    Vec xn(num_of_elem, fill::zeros);
    for (uword i = 0; i < ns; ++i)
        xn(sub2ind(size(ns, num_of_eq), i, CNI)) = (float)(i + 1) * 10.0f;

    const Vec xp = ip1(xn);
    const Vec xm = im1(xn);

    // ip1 shifts values from i+1 down into slot i (interior only)
    for (uword i = 0; i < ns - 1; ++i)
        EXPECT_NEAR(xp(sub2ind(size(ns, num_of_eq), i, CNI)),
                    xn(sub2ind(size(ns, num_of_eq), i + 1, CNI)), 1e-6);

    // im1 shifts values from i-1 up into slot i (interior only)
    for (uword i = 1; i < ns; ++i)
        EXPECT_NEAR(xm(sub2ind(size(ns, num_of_eq), i, CNI)),
                    xn(sub2ind(size(ns, num_of_eq), i - 1, CNI)), 1e-6);
}

static void test_flux_lim_is_minmod() {
    // Code labels it ospre in docstring, but the implementation is minmod:
    //   φ(r) = max(0, min(1, r))   (writeup eq 15)
    Vec r(6);
    r(0) = -1.0f;  // negative -> 0
    r(1) =  0.0f;  // zero     -> 0
    r(2) =  0.5f;  // 0 < r < 1 -> r
    r(3) =  1.0f;  // r == 1    -> 1
    r(4) =  2.0f;  // r  > 1    -> 1
    r(5) =  std::numeric_limits<float>::quiet_NaN(); // NaN -> 0 (sentinel)
    Vec out = flux_lim(r);
    EXPECT_NEAR(out(0), 0.0f, 1e-6);
    EXPECT_NEAR(out(1), 0.0f, 1e-6);
    EXPECT_NEAR(out(2), 0.5f, 1e-6);
    EXPECT_NEAR(out(3), 1.0f, 1e-6);
    EXPECT_NEAR(out(4), 1.0f, 1e-6);
    EXPECT_NEAR(out(5), 0.0f, 1e-6);
}

static void test_cons_prim_roundtrip() {
    Vec cons = setup_uniform(8, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    Vec prim  = cons2prim(cons);
    Vec cons2 = prim2cons(prim);
    for (uword j = 0; j < cons.n_elem; ++j) EXPECT_REL(cons2(j), cons(j), 1e-4);
}


// =========================================================================
// Physics formulas vs. writeup
// =========================================================================

static void test_pressure_relation_eq38() {
    // p_i = 2 n_i k_b T_i (electron quasi-neutrality folds Te=Ti, ne=ni).
    // p_n = n_n k_b T_n. Recover from conserved state via cons2prim
    // and check directly.
    const float Ti = 6500.0f, Tn = 6500.0f;
    const float ni = 2.0e17f,  nn = 1.0e19f;
    Vec cons = setup_uniform(4, ni, nn, Ti, Tn);
    Vec prim = cons2prim(cons);
    const float p_i = prim(sub2ind(size(ns, num_of_eq), 0, PPI));
    const float p_n = prim(sub2ind(size(ns, num_of_eq), 0, PPN));
    EXPECT_REL(p_i, 2.0f * ni * k_b * Ti, 1e-4);
    EXPECT_REL(p_n, 1.0f * nn * k_b * Tn, 1e-4);
}

static void test_spectral_radius_uniform_at_rest() {
    // For motionless state, spectral radius = max(c_s,i, c_s,n)
    // c_s,i = sqrt(γ p_i / ρ_i) = sqrt(2 γ k_b T_i / m_i)    (writeup eq 64/65)
    // c_s,n = sqrt(γ p_n / ρ_n) = sqrt(  γ k_b T_n / m_n)
    const float Ti = 6500.0f, Tn = 6500.0f;
    const float ni = 2.0e17f,  nn = 1.0e19f;
    Vec cons = setup_uniform(4, ni, nn, Ti, Tn);

    Vec spc = find_spectral_radius_ii(cons);
    const double csi = std::sqrt(2.0 * (double)gammamono * (double)k_b * Ti / (double)m_i);
    const double csn = std::sqrt(      (double)gammamono * (double)k_b * Tn / (double)m_n);
    const double expected = std::max(csi, csn);

    for (uword i = 0; i < ns; ++i)
        EXPECT_REL(spc(sub2ind(size(ns, num_of_eq), i, 0)), expected, 1e-3);
}

static void test_nu_in_collision_formula() {
    // Code: ν_in = (2 a₀)² n_n √(8π k_b (T_i+T_n) / m_i)        (writeup eq 57, target-density form)
    chromo_init(1, 0.25f);
    const float nn_val = 1.0e19f, Ti = 6500.0f, Tn = 6500.0f;
    Vec nn_v(1), Ti_v(1), Tn_v(1);
    nn_v(0) = nn_val; Ti_v(0) = Ti; Tn_v(0) = Tn;
    const Vec nuin = nu_in(nn_v, Ti_v, Tn_v);

    const double a0 = 53e-12;
    const double expected = (2.0 * a0) * (2.0 * a0) * (double)nn_val
        * std::sqrt(8.0 * (double)pi * (double)k_b * ((double)Ti + (double)Tn) / (double)m_i);
    EXPECT_REL(nuin(0), expected, 1e-3);
}

static void test_kappa_e_eq53() {
    // κ_e = 9.2048e-12 · n_e · T_e^(5/2) / (n_e + 3.5609e-12 · n_n · T_e²)
    chromo_init(1, 0.25f);
    const float ne_val = 2.0e17f, nn_val = 1.0e19f, Te = 6500.0f;
    Vec ne_v(1), nn_v(1), Te_v(1);
    ne_v(0) = ne_val; nn_v(0) = nn_val; Te_v(0) = Te;
    Vec ke = kappa_e(ne_v, nn_v, Te_v);

    const double num = 9.2048e-12 * (double)ne_val * std::pow((double)Te, 2.5);
    const double den = (double)ne_val + 3.5609e-12 * (double)nn_val * (double)Te * (double)Te;
    EXPECT_REL(ke(0), num / den, 1e-3);

    // Sanity: writeup eq 54 gives ~0.031 W/(m·K) at chromosphere conditions
    EXPECT_TRUE(ke(0) > 0.02f && ke(0) < 0.05f);
}

static void test_kappa_n_eq59() {
    // κ_n = 0.0342006 · n_n · T_n / (1.20613 · n_i · √(T_n+T_i) + 1.70573 · n_n · √T_n)
    chromo_init(1, 0.25f);
    const float ni_val = 2.0e17f, nn_val = 1.0e19f, Ti = 6500.0f, Tn = 6500.0f;
    Vec ni_v(1), nn_v(1), Ti_v(1), Tn_v(1);
    ni_v(0) = ni_val; nn_v(0) = nn_val; Ti_v(0) = Ti; Tn_v(0) = Tn;
    Vec kn = kappa_n(ni_v, nn_v, Ti_v, Tn_v);

    const double num = 0.0342006 * (double)nn_val * (double)Tn;
    const double den = 1.20613 * (double)ni_val * std::sqrt((double)Tn + (double)Ti)
                     + 1.70573 * (double)nn_val * std::sqrt((double)Tn);
    EXPECT_REL(kn(0), num / den, 1e-3);

    // Neutral conductivity dominates in chromosphere (writeup §3.2): κ_n ~ O(1) W/(m·K)
    EXPECT_TRUE(kn(0) > 0.5f && kn(0) < 2.0f);
}

static void test_kappa_n_dominates_in_chromosphere() {
    // Writeup §3.2: at T ~ 6500K, n_n >> n_e, neutral heat conduction dominates electron.
    chromo_init(1, 0.25f);
    const float Te = 6500.0f, Ti = 6500.0f, Tn = 6500.0f;
    const float ne = 2.0e17f, nn = 1.0e19f; // Model-C7-ish bottom-of-chromosphere
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
    Vec cons = setup_uniform(10, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    Vec dt = cal_dt_i(cons);
    Vec maxv = get_max_v_i(cons);
    const float expected_dt = CFL * arma::min(ds_i) / arma::max(maxv);
    // cal_dt_i returns a uniform vector = min over all cells
    EXPECT_REL(dt(0), expected_dt, 1e-3);
    // CFL constraint must hold pointwise: dt * max_v / ds <= CFL
    for (uword i = 0; i < ns; ++i)
        EXPECT_TRUE(dt(i) * maxv(i) / ds_i(i) <= CFL + 1e-6f);
}

static void test_get_max_v_matches_spectral_radius() {
    Vec cons = setup_uniform(5, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    Vec spc = find_spectral_radius_ii(cons);
    Vec mv  = get_max_v_i(cons);
    for (uword i = 0; i < ns; ++i)
        EXPECT_NEAR(mv(i), spc(sub2ind(size(ns, num_of_eq), i, 0)), 1e-6);
}

static void test_uniform_state_is_fixed_point() {
    // Uniform, motionless, no gravity, B uniform → exact fixed point of explicit step.
    // With dinvB_ds=0 the source term cal_S_ii vanishes; with uniform L=R reconstructions
    // the Rusanov flux is identical across all faces.
    const float Ti = 6500.0f, Tn = 6500.0f, ni = 2.0e17f, nn = 1.0e19f;
    Vec cons  = setup_uniform(16, ni, nn, Ti, Tn);
    Vec cons0 = cons;
    Vec dt    = cal_dt_i(cons);
    Vec cons1 = advance_Euler_ii(cons, dt);

    for (uword j = 0; j < cons.n_elem; ++j)
        EXPECT_REL(cons1(j), cons0(j), 1e-3);
}

static void test_mass_conservation_on_uniform_state() {
    // Total ion and neutral mass preserved exactly across one explicit step
    // on the uniform fixed-point setup.
    const float ni = 2.0e17f, nn = 1.0e19f;
    Vec cons = setup_uniform(16, ni, nn, 6500.0f, 6500.0f);
    Vec dt   = cal_dt_i(cons);

    const double mi0 = arma::sum(get_scalar(cons, CNI));
    const double mn0 = arma::sum(get_scalar(cons, CNN));

    Vec cons1 = advance_Euler_ii(cons, dt);

    const double mi1 = arma::sum(get_scalar(cons1, CNI));
    const double mn1 = arma::sum(get_scalar(cons1, CNN));
    EXPECT_REL(mi1, mi0, 1e-4);
    EXPECT_REL(mn1, mn0, 1e-4);
}

static void test_rk4_uniform_fixed_point() {
    // RK4 must also leave the uniform state untouched.
    const float Ti = 6500.0f, Tn = 6500.0f, ni = 2.0e17f, nn = 1.0e19f;
    Vec cons  = setup_uniform(8, ni, nn, Ti, Tn);
    Vec cons0 = cons;
    Vec dt    = cal_dt_i(cons);
    Vec cons1 = advance_RK4(cons, dt);
    for (uword j = 0; j < cons.n_elem; ++j)
        EXPECT_REL(cons1(j), cons0(j), 1e-3);
}


// =========================================================================
// Driver
// =========================================================================

int main() {
    std::cout << "===== Chromosphere test suite =====\n\n";

    // Helper invariants
    RUN(test_scalar_to_get_scalar_inverse);
    RUN(test_ip1_im1_interior_shift);
    RUN(test_flux_lim_is_minmod);
    RUN(test_cons_prim_roundtrip);

    // Physics vs. writeup
    RUN(test_pressure_relation_eq38);
    RUN(test_spectral_radius_uniform_at_rest);
    RUN(test_nu_in_collision_formula);
    RUN(test_kappa_e_eq53);
    RUN(test_kappa_n_eq59);
    RUN(test_kappa_n_dominates_in_chromosphere);

    // Numerics
    RUN(test_cal_dt_respects_cfl);
    RUN(test_get_max_v_matches_spectral_radius);
    RUN(test_uniform_state_is_fixed_point);
    RUN(test_mass_conservation_on_uniform_state);
    RUN(test_rk4_uniform_fixed_point);

    std::cout << "\n===== Summary =====\n";
    std::cout << "Passed: " << g_pass << "\n";
    std::cout << "Failed: " << g_fail << "\n";
    return g_fail == 0 ? 0 : 1;
}
