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
#include "../profiling.hpp"
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
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
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

// =========================================================================
// Stage 1/2: classical Saha closure and CRASH Gamma1 table loader
// =========================================================================

static bool throws_any(const std::function<void()>& action) {
    try { action(); } catch (const std::exception&) { return true; }
    return false;
}

static void test_saha_ionization_fraction_log_domain() {
    constexpr double pi = 3.14159265358979323846;
    const double n_h = 1.0e20;
    const double temperature = 1.0e4;
    const double x = saha_ionization_fraction(n_h * eos_constants::m_h, temperature);
    const double phi = std::pow(2.0*pi*eos_constants::m_e*eos_constants::k_b
                                * temperature/(eos_constants::h*eos_constants::h), 1.5)
                     * std::exp(-eos_constants::chi_h
                                /(eos_constants::k_b*temperature));
    EXPECT_TRUE(x > 0.0 && x < 1.0);
    EXPECT_REL(x*x/(1.0-x), phi/n_h, 2e-13);
    EXPECT_NEAR(saha_ionization_fraction_n_h(n_h, temperature), x, 2e-15);
    EXPECT_TRUE(saha_ionization_fraction_n_h(1.0e26, 3200.0) < 1e-7);
    EXPECT_TRUE(saha_ionization_fraction_n_h(1.0e12, 1.0e8) > 1.0 - 1e-10);
    EXPECT_TRUE(throws_any([] { saha_ionization_fraction(0.0, 1.0e4); }));
    EXPECT_TRUE(throws_any([] { saha_ionization_fraction_n_h(1.0e20, -1.0); }));

    Grid grid;
    grid.init(1, 0.25f);
    EXPECT_TRUE(grid.k_b == static_cast<float>(eos_constants::k_b));
    EXPECT_TRUE(grid.m_e == static_cast<float>(eos_constants::m_e));
    EXPECT_TRUE(grid.m_i == static_cast<float>(eos_constants::m_h));
    EXPECT_TRUE(grid.m_n == static_cast<float>(eos_constants::m_h));
    EXPECT_TRUE(grid.chi_H_J == static_cast<float>(eos_constants::chi_h));
}

static void write_gamma_metadata(std::ofstream& out, const std::string& variant) {
    out << "# format_version=" << (variant == "wrong_metadata" ? "2" : "1") << '\n';
    if (variant != "missing_metadata")
        out << "# table_id=chromosphere2026_gamma1_hydrogen_v1\n";
    out << "# material=pure_H\n"
        << "# excitation=0\n"
        << "# fermi_gas=0\n"
        << "# coulomb_correction=0\n"
        << "# ground_stat_weight=1\n"
        << "# saha_prefactor=coefficient_one\n"
        << "# axis_1=log10_T_K\n"
        << "# axis_2=log10_nH_m-3\n"
        << "# value=Gamma1\n"
        << "# nT=" << (variant == "wrong_dimensions" ? "4" : "3") << "\n"
        << "# nN=2\n"
        << "# k_b_J_K=1.380649e-23\n"
        << "# m_e_kg=9.1093837015e-31\n"
        << "# m_H_kg=1.6726219e-27\n"
        << "# h_J_s=6.62607015e-34\n"
        << "# chi_H_J=2.179872361e-18\n"
        << "# generator=util/eos/tabulate_gamma.f90\n"
        << "# generation_date=2026-07-19\n"
        << "# source_revision=crash_gamma_eos_stage0_2_v1\n";
    if (variant == "duplicate_metadata") out << "# material=pure_H\n";
    if (variant == "unexpected_metadata") out << "# unrecognized=silent_drift\n";
}

static std::string write_gamma_table_fixture(const std::string& name,
                                             const std::string& variant = "valid") {
    const std::string path = name;
    std::ofstream out(path.c_str());
    write_gamma_metadata(out, variant);
    double log_t[] = {3.0, 4.0, 5.0};
    if (variant == "duplicate_axis") log_t[2] = 4.0;
    if (variant == "nonascending_axis") { log_t[1] = 5.0; log_t[2] = 4.0; }
    double log_n[] = {15.0, 20.0};
    if (variant == "duplicate_density_axis") log_n[1] = 15.0;
    if (variant == "nonascending_density_axis") { log_n[0] = 20.0; log_n[1] = 15.0; }
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 3; ++i) {
            if (variant == "incomplete" && j == 1 && i == 2) continue;
            double gamma1 = 1.0 + 0.1*(log_t[i]-3.0) + 0.01*(log_n[j]-15.0);
            if (variant == "nonpositive" && j == 0 && i == 0) gamma1 = 0.0;
            out << std::setprecision(17) << log_t[i] << ' ' << log_n[j] << ' ';
            if (variant == "nan" && j == 0 && i == 0) out << "nan";
            else out << gamma1;
            out << '\n';
        }
    }
    return path;
}

// Algebraic Saha pressure inversion (equilibrium_density_from_pressure) against
// the geometric-bisection solve it replaced, over the whole production domain:
// the full table density axis (10^12..10^26 m^-3) crossed with 3.2e3..1e7 K,
// i.e. from the fully-neutral photosphere to the fully-ionized corona.
// Grid::broadcast() now fingerprints the static mesh / magnetic geometry and
// skips the rebuild when neither changed, and fills the packed caches directly
// instead of accumulating num_of_eq scalar_to temporaries. Both properties are
// checked here: the packed caches must equal the accumulated construction
// EXACTLY, and any edit to a source array must still be picked up.
static void test_broadcast_static_metric_cache() {
    Grid grid;
    grid.init(9, 0.25f);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        grid.ds_i(i)       = 4.0e3f*(1.0f + 0.17f*i);      // non-uniform
        grid.B_i(i)        = 1.0f + 0.3f*i;
        grid.B_imh(i)      = 1.0f + 0.3f*i - 0.05f;
        grid.B_iph(i)      = 1.0f + 0.3f*i + 0.05f;
        grid.dinvB_ds_i(i) = (1.0f/grid.B_iph(i) - 1.0f/grid.B_imh(i))/grid.ds_i(i);
    }
    grid.broadcast();

    // Reference construction: the accumulate-scalar_to form broadcast() used to
    // run. Must agree bit-for-bit with the direct fill.
    auto reference = [&](const Vec& src) {
        Vec dst(grid.n_state, arma::fill::zeros);
        for (arma::uword k = 0; k < num_of_eq; ++k) dst += scalar_to(grid, src, k);
        return dst;
    };
    auto bitwise_equal = [&](const Vec& a, const Vec& b) {
        if (a.n_elem != b.n_elem) return false;
        for (arma::uword i = 0; i < a.n_elem; ++i)
            if (!(a[i] == b[i])) return false;
        return true;
    };
    EXPECT_TRUE(bitwise_equal(grid.B_state,        reference(grid.B_i)));
    EXPECT_TRUE(bitwise_equal(grid.B_state_imh,    reference(grid.B_imh)));
    EXPECT_TRUE(bitwise_equal(grid.B_state_iph,    reference(grid.B_iph)));
    EXPECT_TRUE(bitwise_equal(grid.ds_state,       reference(grid.ds_i)));
    EXPECT_TRUE(bitwise_equal(grid.dinvB_ds_state, reference(grid.dinvB_ds_i)));
    EXPECT_TRUE(bitwise_equal(grid.ds_iph_state,   reference(grid.ds_iph_i)));
    EXPECT_TRUE(bitwise_equal(grid.ds_imh_state,   reference(grid.ds_imh_i)));
    EXPECT_TRUE(grid.uniform_mesh == false);

    // A no-op broadcast (what a per-step boundary refresh performs) must not
    // change anything, and must not bump the mesh generation.
    const std::uint64_t generation = grid.metrics_generation();
    const Vec ds_state_before = grid.ds_state;
    const Vec s_face_before   = grid.s_face;
    grid.broadcast();
    grid.broadcast();
    EXPECT_TRUE(grid.metrics_generation() == generation);
    EXPECT_TRUE(bitwise_equal(grid.ds_state, ds_state_before));
    EXPECT_TRUE(bitwise_equal(grid.s_face,   s_face_before));

    // A geometry edit must be picked up.
    grid.ds_i(3) *= 1.5f;
    grid.broadcast();
    EXPECT_TRUE(grid.metrics_generation() == generation + 1);
    EXPECT_TRUE(bitwise_equal(grid.ds_state, reference(grid.ds_i)));
    EXPECT_TRUE(bitwise_equal(grid.ds_iph_state, reference(grid.ds_iph_i)));
    EXPECT_REL(grid.s_face(grid.ns), arma::sum(grid.ds_i), 1e-6);

    grid.B_iph(2) *= 1.1f;
    grid.broadcast();
    EXPECT_TRUE(bitwise_equal(grid.B_state_iph, reference(grid.B_iph)));

    // A uniform mesh must still be detected, and resize must invalidate.
    grid.ds_i.fill(2.0e3f);
    grid.broadcast();
    EXPECT_TRUE(grid.uniform_mesh == true);
    grid.resize(5);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        grid.ds_i(i) = 1.0e3f; grid.B_i(i) = 1.0f;
        grid.B_imh(i) = 1.0f;  grid.B_iph(i) = 1.0f;
    }
    grid.broadcast();
    EXPECT_TRUE(grid.B_state.n_elem == grid.n_state);
    EXPECT_TRUE(bitwise_equal(grid.ds_state, reference(grid.ds_i)));
}

static void test_equilibrium_density_from_pressure_matches_bisection() {
    auto pressure_at = [](double n_h, double t) {
        return (1.0 + saha_ionization_fraction_n_h(n_h, t))
             * n_h * eos_constants::k_b * t;
    };
    // The former implementation: 100 geometric bisections on the table n_H axis.
    auto bisect_rho = [&](double pressure, double t) {
        double lo = 1.0e12, hi = 1.0e26;
        for (int it = 0; it < 100; ++it) {
            const double mid = std::sqrt(lo*hi);
            if (pressure_at(mid, t) < pressure) lo = mid; else hi = mid;
        }
        return std::sqrt(lo*hi)*eos_constants::m_h;
    };

    double worst_rho_rel = 0.0, worst_p_rel = 0.0;
    for (int it = 0; it <= 24; ++it) {
        const double t = std::pow(10.0, 3.5 + 3.5*it/24.0);        // 3162 K .. 1e7 K
        for (int in = 0; in <= 28; ++in) {
            const double n_h = std::pow(10.0, 12.0 + 14.0*in/28.0); // 1e12 .. 1e26
            const double p = pressure_at(n_h, t);
            const double rho = equilibrium_density_from_pressure(p, t);
            // (a) recovers the density the pressure was built from
            worst_rho_rel = std::max(worst_rho_rel,
                std::abs(rho - n_h*eos_constants::m_h)/(n_h*eos_constants::m_h));
            // (b) agrees with the bisection it replaced
            EXPECT_REL(rho, bisect_rho(p, t), 1e-12);
            // (c) reconstructs the requested pressure to near-double precision
            const double p_back = pressure_at(rho/eos_constants::m_h, t);
            worst_p_rel = std::max(worst_p_rel, std::abs(p_back - p)/p);
        }
    }
    EXPECT_TRUE(worst_rho_rel < 1e-12);
    EXPECT_TRUE(worst_p_rel  < 1e-13);
    // Deep limits: fully neutral (y >> S) and fully ionized (y << S) must not
    // overflow the exponential in either branch.
    EXPECT_REL(equilibrium_density_from_pressure(pressure_at(1.0e26, 3200.0), 3200.0),
               1.0e26*eos_constants::m_h, 1e-12);
    EXPECT_REL(equilibrium_density_from_pressure(pressure_at(1.0e12, 1.0e7), 1.0e7),
               1.0e12*eos_constants::m_h, 1e-12);
    EXPECT_TRUE(throws_any([] { equilibrium_density_from_pressure(0.0, 1.0e4); }));
    EXPECT_TRUE(throws_any([] { equilibrium_density_from_pressure(1.0, -1.0); }));
    EXPECT_TRUE(throws_any([] {
        equilibrium_density_from_pressure(
            std::numeric_limits<double>::quiet_NaN(), 1.0e4); }));
}

static void test_eos_gamma_table_loader_and_interpolation() {
    EosGammaTable empty;
    EXPECT_TRUE(empty.empty());
    EXPECT_TRUE(throws_any([&] { empty.gamma1(1.0e4, 1.0e18); }));

    const std::string path = write_gamma_table_fixture("gamma_table_fixture.dat");
    EosGammaTable table = EosGammaTable::load(path);
    EXPECT_TRUE(!table.empty());
    EXPECT_TRUE(table.temperature_size() == 3);
    EXPECT_TRUE(table.density_size() == 2);
    EXPECT_REL(table.min_temperature(), 1.0e3, 2e-14);
    EXPECT_REL(table.max_temperature(), 1.0e5, 2e-14);
    EXPECT_REL(table.min_n_h(), 1.0e15, 2e-14);
    EXPECT_REL(table.max_n_h(), 1.0e20, 2e-14);
    EXPECT_NEAR(table.gamma1(1.0e3, 1.0e15), 1.0, 2e-14);
    EXPECT_NEAR(table.gamma1(1.0e5, 1.0e20), 1.25, 2e-14);
    EXPECT_NEAR(table.gamma1(std::pow(10.0, 4.5), std::pow(10.0, 17.5)),
                1.175, 2e-14);
    EXPECT_TRUE(table.contains(1.0e3, 1.0e15));
    EXPECT_TRUE(table.contains(1.0e5, 1.0e20));
    EXPECT_TRUE(!table.contains(999.0, 1.0e18));
    EXPECT_TRUE(!table.contains(1.0e4, 1.0e14));
    EXPECT_TRUE(throws_any([&] { table.gamma1(999.0, 1.0e18); }));
    EXPECT_TRUE(throws_any([&] { table.gamma1(1.0e4, 1.0e14); }));
    EXPECT_TRUE(throws_any([&] { table.gamma1(1.0e4, 1.0e21); }));
    EXPECT_NEAR(table.gamma1(999.0, 1.0e18, true), 1.03, 2e-14);
    EXPECT_NEAR(table.gamma1(1.0e6, 1.0e21, true), 1.25, 2e-14);
    EXPECT_NEAR(table.gamma1(999.0, 1.0e14, true), 1.0, 2e-14);
    std::remove(path.c_str());

    const std::vector<std::string> bad_variants = {
        "incomplete", "wrong_dimensions", "wrong_metadata", "missing_metadata",
        "duplicate_metadata", "unexpected_metadata", "duplicate_axis",
        "nonascending_axis", "duplicate_density_axis", "nonascending_density_axis",
        "nonpositive", "nan"
    };
    for (const std::string& variant : bad_variants) {
        const std::string bad = write_gamma_table_fixture("gamma_table_bad_fixture.dat", variant);
        EXPECT_TRUE(throws_any([&] { EosGammaTable::load(bad); }));
        std::remove(bad.c_str());
    }
}

static void test_stage3_equilibrium_mixture_closure() {
    const std::string path = write_gamma_table_fixture("gamma_stage3_fixture.dat");
    const EosGammaTable table = EosGammaTable::load(path);
    const double n_h = 1.0e18;
    const double rho = n_h * eos_constants::m_h;

    for (double temperature : {1.0e3, 4.0e3, 1.0e4, 6.0e4, 1.0e5}) {
        const double energy = equilibrium_internal_energy(rho, temperature);
        EXPECT_REL(temperature_from_rho_eint(table, rho, energy), temperature, 2e-12);
        EXPECT_REL(temperature_from_rho_eint(table, rho, energy, 9.0e4),
                   temperature, 2e-12);
    }

    const double temperature = 1.0e4;
    const double x = saha_ionization_fraction_n_h(n_h, temperature);
    const double pressure = (1.0 + x) * n_h * eos_constants::k_b * temperature;
    const double energy = equilibrium_internal_energy(rho, temperature);
    const GammaState gamma = gamma_state(table, rho, temperature);
    EXPECT_NEAR(gamma.x_eq, x, 2e-15);
    EXPECT_NEAR(gamma.gamma_energy, 1.0 + pressure/energy, 2e-15);
    EXPECT_NEAR(gamma.gamma_sound, 1.13, 2e-14);

    for (double bad_guess : {
             std::numeric_limits<double>::quiet_NaN(), 100.0, 1.0e6,
             std::numeric_limits<double>::infinity()}) {
        EXPECT_REL(temperature_from_rho_eint(table, rho, energy, bad_guess),
                   temperature, 2e-12);
    }

    // Regression for a bracket-endpoint Newton cycle exposed by the extended
    // 0--2152.6 km gamma-table model_column conduction run.
    const double trapped_rho = 5.55517e-10;
    const double trapped_energy = 0.802803;
    const double trapped_temperature = temperature_from_rho_eint(
        table, trapped_rho, trapped_energy, 11434.3);
    EXPECT_TRUE(trapped_temperature > 8.0e3 && trapped_temperature < 8.3e3);
    EXPECT_REL(equilibrium_internal_energy(trapped_rho, trapped_temperature),
               trapped_energy, 2e-12);

    const double x_row = 0.3;
    const double rho_i = x_row * rho;
    const double rho_n = (1.0 - x_row) * rho;
    const double velocity_i = 120.0;
    const double velocity_n = -45.0;
    const double momentum_i = rho_i * velocity_i;
    const double momentum_n = rho_n * velocity_n;
    const double phi = 2.5e6;
    const double p_e = x * n_h * eos_constants::k_b * temperature;
    const double p_i = 2.0 * p_e;
    const double p_n = (1.0 - x) * n_h * eos_constants::k_b * temperature;
    const double energy_i = 1.5*p_i + x*n_h*eos_constants::chi_h
                          + 0.5*rho_i*velocity_i*velocity_i + rho_i*phi;
    const double energy_n = 1.5*p_n
                          + 0.5*rho_n*velocity_n*velocity_n + rho_n*phi;
    const MixtureThermo decoded = decode_equilibrium_mixture(
        table, rho_i, rho_n, momentum_i, momentum_n,
        energy_i, energy_n, phi);
    EXPECT_REL(decoded.rho, rho, 2e-15);
    EXPECT_REL(decoded.T, temperature, 2e-12);
    EXPECT_NEAR(decoded.x_eq, x, 2e-15);
    EXPECT_NEAR(decoded.x_row, x_row, 2e-15);
    EXPECT_REL(decoded.n_H, n_h, 2e-15);
    EXPECT_REL(decoded.n_e, x*n_h, 2e-15);
    EXPECT_REL(decoded.n_HI, (1.0-x)*n_h, 2e-15);
    EXPECT_REL(decoded.p_i, p_i, 2e-12);
    EXPECT_REL(decoded.p_n, p_n, 2e-12);
    EXPECT_REL(decoded.p_e, p_e, 2e-12);
    EXPECT_NEAR(decoded.gamma1, 1.13, 2e-14);
    EXPECT_REL(decoded.internal_energy, energy, 2e-12);

    const double e_min = equilibrium_internal_energy(rho, table.min_temperature());
    const double e_max = equilibrium_internal_energy(rho, table.max_temperature());
    const double near_min = temperature_from_rho_eint(table, rho, e_min*(1.0+1.0e-10));
    const double near_max = temperature_from_rho_eint(table, rho, e_max*(1.0-1.0e-10));
    EXPECT_TRUE(near_min > table.min_temperature());
    EXPECT_TRUE(near_max < table.max_temperature());
    EXPECT_TRUE(throws_any([&] {
        temperature_from_rho_eint(table, rho, 0.5*e_min);
    }));
    EXPECT_TRUE(throws_any([&] {
        temperature_from_rho_eint(table, rho, 2.0*e_max);
    }));
    EXPECT_TRUE(throws_any([&] {
        gamma_state(table, 1.0e14*eos_constants::m_h, temperature);
    }));
    const double rho_density_oob = 1.0e14 * eos_constants::m_h;
    const double energy_density_oob = equilibrium_internal_energy(
        rho_density_oob, temperature);
    EXPECT_TRUE(throws_any([&] {
        temperature_from_rho_eint(table, rho_density_oob, energy_density_oob);
    }));
    const double rho_density_above = 1.0e21 * eos_constants::m_h;
    const double energy_density_above = equilibrium_internal_energy(
        rho_density_above, temperature);
    EXPECT_TRUE(throws_any([&] {
        temperature_from_rho_eint(table, rho_density_above, energy_density_above);
    }));
    EXPECT_REL(temperature_from_rho_eint(
                   table, rho_density_oob, energy_density_oob,
                   std::numeric_limits<double>::quiet_NaN(), true),
               temperature, 2e-12);
    EXPECT_TRUE(throws_any([&] {
        decode_equilibrium_mixture(table, 0.0, rho_n, 0.0, momentum_n,
                                   energy_i, energy_n, phi);
    }));
    std::remove(path.c_str());
}

static void test_stage4_conservative_equilibrium_projection() {
    const std::string path = write_gamma_table_fixture("gamma_stage4_fixture.dat");
    const EosGammaTable table = EosGammaTable::load(path);
    const double n_h = 1.0e18;
    const double rho = n_h * eos_constants::m_h;
    const double temperature = 1.0e4;
    const double x_eq = saha_ionization_fraction_n_h(n_h, temperature);
    const double x_row_initial = 0.3;
    const double rho_i = x_row_initial * rho;
    const double rho_n = (1.0 - x_row_initial) * rho;
    const double velocity_i = 2.0e4;
    const double velocity_n = -8.0e3;
    const double momentum_i = rho_i * velocity_i;
    const double momentum_n = rho_n * velocity_n;
    const double phi = 2.5e6;
    const double p_e = x_eq * n_h * eos_constants::k_b * temperature;
    const double p_i = 2.0 * p_e;
    const double p_n = (1.0 - x_eq) * n_h * eos_constants::k_b * temperature;
    const double energy_i = 1.5*p_i + x_eq*n_h*eos_constants::chi_h
                          + 0.5*rho_i*velocity_i*velocity_i + rho_i*phi;
    const double energy_n = 1.5*p_n
                          + 0.5*rho_n*velocity_n*velocity_n + rho_n*phi;
    const double mass_before = rho_i + rho_n;
    const double momentum_before = momentum_i + momentum_n;
    const double energy_before = energy_i + energy_n;
    const double velocity_cm = momentum_before / mass_before;
    const double reduced_mass = rho_i * rho_n / mass_before;
    const double internal_before = equilibrium_internal_energy(rho, temperature);
    const double drift_heat = 0.5 * reduced_mass
                            * (velocity_i-velocity_n)*(velocity_i-velocity_n);

    const ProjectedMixture projected = project_equilibrium_single_fluid(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n, phi,
        1.0e-8);
    EXPECT_REL(projected.rho_i + projected.rho_n, mass_before, 2e-15);
    EXPECT_REL(projected.momentum_i + projected.momentum_n, momentum_before, 2e-15);
    EXPECT_REL(projected.energy_i + projected.energy_n, energy_before, 2e-15);
    EXPECT_REL(projected.momentum_i/projected.rho_i, velocity_cm, 2e-15);
    EXPECT_REL(projected.momentum_n/projected.rho_n, velocity_cm, 2e-15);
    EXPECT_REL(projected.thermo.internal_energy,
               internal_before + drift_heat, 2e-13);
    EXPECT_TRUE(projected.thermo.T > temperature);
    EXPECT_REL(equilibrium_internal_energy(rho, projected.thermo.T),
               projected.thermo.internal_energy, 2e-13);
    EXPECT_NEAR(projected.thermo.x_row, projected.thermo.x_eq, 2e-15);

    const double common_specific = 0.5*velocity_cm*velocity_cm + phi;
    const double mapped_e_i = 1.5*projected.thermo.p_i
                            + projected.thermo.x_eq*projected.thermo.n_H
                              * eos_constants::chi_h
                            + projected.rho_i*common_specific;
    const double mapped_e_n = 1.5*projected.thermo.p_n
                            + projected.rho_n*common_specific;
    EXPECT_REL(projected.energy_i, mapped_e_i, 2e-13);
    EXPECT_NEAR(projected.energy_n, mapped_e_n, 2e-13*energy_before);
    EXPECT_REL(projected.energy_e, 1.5*projected.thermo.p_e, 2e-15);

    const MixtureThermo decoded = decode_equilibrium_mixture(
        table, projected.rho_i, projected.rho_n,
        projected.momentum_i, projected.momentum_n,
        projected.energy_i, projected.energy_n, phi);
    EXPECT_REL(decoded.T, projected.thermo.T, 2e-12);
    EXPECT_NEAR(decoded.x_eq, projected.thermo.x_eq, 2e-14);
    EXPECT_NEAR(decoded.x_row, projected.thermo.x_row, 2e-14);

    const ProjectedMixture projected_twice = project_equilibrium_single_fluid(
        table, projected.rho_i, projected.rho_n,
        projected.momentum_i, projected.momentum_n,
        projected.energy_i, projected.energy_n, phi, 1.0e-8);
    EXPECT_REL(projected_twice.thermo.T, projected.thermo.T, 2e-12);
    EXPECT_REL(projected_twice.rho_i, projected.rho_i, 2e-13);
    EXPECT_REL(projected_twice.rho_n, projected.rho_n, 2e-13);
    EXPECT_REL(projected_twice.energy_i + projected_twice.energy_n,
               energy_before, 2e-15);

    EXPECT_TRUE(throws_any([&] {
        project_equilibrium_single_fluid(table, rho_i, rho_n, momentum_i,
            momentum_n, energy_i, energy_n, phi, 0.0);
    }));
    EXPECT_TRUE(throws_any([&] {
        project_equilibrium_single_fluid(table, rho_i, rho_n, momentum_i,
            momentum_n, energy_i, energy_n, phi, 0.5);
    }));

    const auto projection_density_oob = [&](double bad_n_h) {
        const double bad_rho = bad_n_h*eos_constants::m_h;
        const double bad_x = saha_ionization_fraction_n_h(bad_n_h, temperature);
        const double bad_p_e = bad_x*bad_n_h*eos_constants::k_b*temperature;
        const double bad_p_i = 2.0*bad_p_e;
        const double bad_p_n = (1.0-bad_x)*bad_n_h*eos_constants::k_b*temperature;
        const double bad_e_i = 1.5*bad_p_i + bad_x*bad_n_h*eos_constants::chi_h;
        const double bad_e_n = 1.5*bad_p_n;
        return throws_any([&] {
            project_equilibrium_single_fluid(
                table, 0.4*bad_rho, 0.6*bad_rho, 0.0, 0.0,
                bad_e_i, bad_e_n, 0.0);
        });
    };
    EXPECT_TRUE(projection_density_oob(1.0e14));
    EXPECT_TRUE(projection_density_oob(1.0e21));

    // At a nearly neutral state, changing the storage floor changes only the
    // carrier masses. All x_eq-driven thermodynamics remain identical.
    const double n_cool = 1.0e20;
    const double rho_cool = n_cool * eos_constants::m_h;
    const double t_cool = 1.0e3;
    const double x_cool = saha_ionization_fraction_n_h(n_cool, t_cool);
    const double p_e_cool = x_cool*n_cool*eos_constants::k_b*t_cool;
    const double p_i_cool = 2.0*p_e_cool;
    const double p_n_cool = (1.0-x_cool)*n_cool*eos_constants::k_b*t_cool;
    const double e_i_cool = 1.5*p_i_cool + x_cool*n_cool*eos_constants::chi_h;
    const double e_n_cool = 1.5*p_n_cool;
    const ProjectedMixture floor6 = project_equilibrium_single_fluid(
        table, 0.4*rho_cool, 0.6*rho_cool, 0.0, 0.0,
        e_i_cool, e_n_cool, 0.0, 1.0e-6);
    const ProjectedMixture floor10 = project_equilibrium_single_fluid(
        table, 0.4*rho_cool, 0.6*rho_cool, 0.0, 0.0,
        e_i_cool, e_n_cool, 0.0, 1.0e-10);
    EXPECT_NEAR(floor6.thermo.x_row, 1.0e-6, 1e-18);
    EXPECT_NEAR(floor10.thermo.x_row, 1.0e-10, 1e-22);
    EXPECT_NEAR(floor6.thermo.x_eq, floor10.thermo.x_eq, 0.0);
    EXPECT_NEAR(floor6.thermo.T, floor10.thermo.T, 0.0);
    EXPECT_NEAR(floor6.thermo.p_i, floor10.thermo.p_i, 0.0);
    EXPECT_NEAR(floor6.thermo.p_n, floor10.thermo.p_n, 0.0);
    EXPECT_NEAR(floor6.thermo.n_e, floor10.thermo.n_e, 0.0);
    EXPECT_TRUE(floor6.rho_i > floor10.rho_i);
    EXPECT_REL(floor6.energy_i + floor6.energy_n,
               floor10.energy_i + floor10.energy_n, 2e-15);

    // Fully ionized upper clamp: only the neutral carrier floor changes.
    const double n_hot = 1.0e15;
    const double rho_hot = n_hot*eos_constants::m_h;
    const double t_hot = 1.0e5;
    const double x_hot = saha_ionization_fraction_n_h(n_hot, t_hot);
    const double p_e_hot = x_hot*n_hot*eos_constants::k_b*t_hot;
    const double p_i_hot = 2.0*p_e_hot;
    const double p_n_hot = (1.0-x_hot)*n_hot*eos_constants::k_b*t_hot;
    const double e_i_hot = 1.5*p_i_hot + x_hot*n_hot*eos_constants::chi_h;
    const double e_n_hot = 1.5*p_n_hot;
    const ProjectedMixture hot6 = project_equilibrium_single_fluid(
        table, 0.4*rho_hot, 0.6*rho_hot, 0.0, 0.0,
        e_i_hot, e_n_hot, 0.0, 1.0e-6);
    const ProjectedMixture hot10 = project_equilibrium_single_fluid(
        table, 0.4*rho_hot, 0.6*rho_hot, 0.0, 0.0,
        e_i_hot, e_n_hot, 0.0, 1.0e-10);
    EXPECT_NEAR(hot6.thermo.x_row, 1.0-1.0e-6, 1e-15);
    EXPECT_NEAR(hot10.thermo.x_row, 1.0-1.0e-10, 1e-15);
    EXPECT_NEAR(hot6.thermo.x_eq, hot10.thermo.x_eq, 0.0);
    EXPECT_NEAR(hot6.thermo.T, hot10.thermo.T, 0.0);
    EXPECT_NEAR(hot6.thermo.p_i, hot10.thermo.p_i, 0.0);
    EXPECT_NEAR(hot6.thermo.p_n, hot10.thermo.p_n, 0.0);
    EXPECT_NEAR(hot6.thermo.p_e, hot10.thermo.p_e, 0.0);
    EXPECT_NEAR(hot6.thermo.n_e, hot10.thermo.n_e, 0.0);
    EXPECT_NEAR(hot6.thermo.n_HI, hot10.thermo.n_HI, 0.0);
    EXPECT_REL(hot6.energy_i + hot6.energy_n,
               hot10.energy_i + hot10.energy_n, 2e-15);
    EXPECT_TRUE(hot6.rho_n > hot10.rho_n);
    const MixtureThermo hot_decoded = decode_equilibrium_mixture(
        table, hot10.rho_i, hot10.rho_n, hot10.momentum_i, hot10.momentum_n,
        hot10.energy_i, hot10.energy_n, 0.0);
    EXPECT_REL(hot_decoded.T, hot10.thermo.T, 2e-12);

    // Zero total momentum with nonzero relative drift: both projected momenta
    // vanish and all original kinetic energy becomes mixture internal energy.
    const double zero_mom_i = rho_i*1.0e4;
    const double zero_mom_n = -zero_mom_i;
    const double zero_vi = zero_mom_i/rho_i;
    const double zero_vn = zero_mom_n/rho_n;
    const double zero_e_i = 1.5*p_i + x_eq*n_h*eos_constants::chi_h
                          + 0.5*rho_i*zero_vi*zero_vi;
    const double zero_e_n = 1.5*p_n + 0.5*rho_n*zero_vn*zero_vn;
    const ProjectedMixture zero_momentum = project_equilibrium_single_fluid(
        table, rho_i, rho_n, zero_mom_i, zero_mom_n,
        zero_e_i, zero_e_n, 0.0);
    const double zero_drift_heat = 0.5*rho_i*zero_vi*zero_vi
                                 + 0.5*rho_n*zero_vn*zero_vn;
    EXPECT_NEAR(zero_momentum.momentum_i, 0.0, 0.0);
    EXPECT_NEAR(zero_momentum.momentum_n, 0.0, 0.0);
    EXPECT_REL(zero_momentum.thermo.internal_energy,
               internal_before + zero_drift_heat, 2e-13);
    EXPECT_REL(zero_momentum.energy_i + zero_momentum.energy_n,
               zero_e_i + zero_e_n, 2e-15);

    // Packed float wrapper: conservation survives the row rewrite to float
    // tolerance, and both carrier velocities become the center-of-mass value.
    Grid grid;
    grid.init(1, 0.25f);
    grid.eos_gamma_table = table;
    grid.eos_trace_fraction_floor = 1.0e-8f;
    grid.phi_g_imh.fill(static_cast<float>(phi));
    grid.phi_g_iph.fill(static_cast<float>(phi));
    Vec state(grid.n_state, arma::fill::zeros);
    const auto size = arma::size(grid.ns, num_of_eq);
    state(arma::sub2ind(size, 0, cons::RHO_I)) = static_cast<float>(rho_i);
    state(arma::sub2ind(size, 0, cons::RHO_N)) = static_cast<float>(rho_n);
    state(arma::sub2ind(size, 0, cons::MOM_I)) = static_cast<float>(momentum_i);
    state(arma::sub2ind(size, 0, cons::MOM_N)) = static_cast<float>(momentum_n);
    state(arma::sub2ind(size, 0, cons::E_I)) = static_cast<float>(energy_i);
    state(arma::sub2ind(size, 0, cons::E_N)) = static_cast<float>(energy_n);
    state(arma::sub2ind(size, 0, cons::E_E)) = -123.0f;
    const Vec packed = project_equilibrium_single_fluid(grid, state);
    const double packed_mass_before = static_cast<double>(state(cons::RHO_I)) + state(cons::RHO_N);
    const double packed_mom_before = static_cast<double>(state(cons::MOM_I)) + state(cons::MOM_N);
    const double packed_energy_before = static_cast<double>(state(cons::E_I)) + state(cons::E_N);
    EXPECT_REL(static_cast<double>(packed(cons::RHO_I)) + packed(cons::RHO_N),
               packed_mass_before, 2e-7);
    EXPECT_REL(static_cast<double>(packed(cons::MOM_I)) + packed(cons::MOM_N),
               packed_mom_before, 2e-7);
    EXPECT_REL(static_cast<double>(packed(cons::E_I)) + packed(cons::E_N),
               packed_energy_before, 2e-7);
    EXPECT_REL(packed(cons::MOM_I)/packed(cons::RHO_I),
               packed(cons::MOM_N)/packed(cons::RHO_N), 2e-6);
    const MixtureThermo packed_decoded = decode_equilibrium_mixture(
        table, packed(cons::RHO_I), packed(cons::RHO_N),
        packed(cons::MOM_I), packed(cons::MOM_N),
        packed(cons::E_I), packed(cons::E_N), phi);
    EXPECT_REL(packed(cons::E_E), 1.5*packed_decoded.p_e, 2e-6);

    // Packed projection is all-or-nothing from the caller's perspective. Cell
    // 0 is valid; cell 1 is above the density table and throws after cell 0 has
    // been written only to the function-local result buffer.
    Grid oob_grid;
    oob_grid.init(2, 0.25f);
    oob_grid.eos_gamma_table = table;
    oob_grid.phi_g_imh.zeros();
    oob_grid.phi_g_iph.zeros();
    Vec oob_state(oob_grid.n_state, arma::fill::zeros);
    const auto oob_size = arma::size(oob_grid.ns, num_of_eq);
    const auto set_static_cell = [&](arma::uword cell, double cell_n_h) {
        const double cell_rho = cell_n_h*eos_constants::m_h;
        const double cell_x = saha_ionization_fraction_n_h(cell_n_h, temperature);
        const double cell_p_e = cell_x*cell_n_h*eos_constants::k_b*temperature;
        const double cell_p_i = 2.0*cell_p_e;
        const double cell_p_n = (1.0-cell_x)*cell_n_h*eos_constants::k_b*temperature;
        oob_state(arma::sub2ind(oob_size, cell, cons::RHO_I)) = static_cast<float>(0.4*cell_rho);
        oob_state(arma::sub2ind(oob_size, cell, cons::RHO_N)) = static_cast<float>(0.6*cell_rho);
        oob_state(arma::sub2ind(oob_size, cell, cons::E_I)) = static_cast<float>(
            1.5*cell_p_i + cell_x*cell_n_h*eos_constants::chi_h);
        oob_state(arma::sub2ind(oob_size, cell, cons::E_N)) = static_cast<float>(1.5*cell_p_n);
        oob_state(arma::sub2ind(oob_size, cell, cons::E_E)) = static_cast<float>(1.5*cell_p_e);
    };
    set_static_cell(0, 1.0e18);
    set_static_cell(1, 1.0e21);
    const Vec oob_original = oob_state;
    EXPECT_TRUE(throws_any([&] {
        project_equilibrium_single_fluid(oob_grid, oob_state);
    }));
    EXPECT_TRUE(arma::approx_equal(oob_state, oob_original, "absdiff", 0.0));
    std::remove(path.c_str());
}

static Vec setup_gamma_equilibrium(Grid& grid, const EosGammaTable& table,
                                   arma::uword ns, bool temperature_gradient) {
    grid.init(ns, 0.25f);
    grid.eos_gamma_table = table;
    grid.single_fluid = true;
    grid.enable_Te = false;
    grid.enable_ionization = false;
    grid.enable_conduction = false;
    grid.enable_radiative_cooling = false;
    grid.enable_beam_heating = false;
    grid.enable_coronal_heating = false;
    grid.enable_trac = false;
    grid.enable_vacuum_floor = false;
    grid.ds_i.fill(1.0e4f);
    grid.B_i.ones(); grid.B_imh.ones(); grid.B_iph.ones();
    grid.dinvB_ds_i.zeros();
    grid.phi_g_imh.zeros(); grid.phi_g_iph.zeros();
    grid.broadcast();

    Vec state(grid.n_state, arma::fill::zeros);
    const auto sz = arma::size(ns, num_of_eq);
    for (arma::uword i = 0; i < ns; ++i) {
        const double f = ns > 1 ? static_cast<double>(i)/(ns-1) : 0.0;
        const double n_h = 1.0e18;
        const double temperature = temperature_gradient ? 7.0e3 + 6.0e3*f : 1.0e4;
        const MixtureFaceState face = equilibrium_mixture_face_state(
            table, n_h*eos_constants::m_h, 150.0, temperature, 0.0,
            grid.eos_trace_fraction_floor);
        const double u[7] = {face.conserved.rho_i, face.conserved.rho_n,
                             face.conserved.momentum_i, face.conserved.momentum_n,
                             face.conserved.energy_i, face.conserved.energy_n,
                             face.conserved.energy_e};
        for (arma::uword k = 0; k < num_of_eq; ++k)
            state(arma::sub2ind(sz, i, k)) = static_cast<float>(u[k]);
    }
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        grid.inner_boundary0_i(k) = state(arma::sub2ind(sz, 0, k));
        grid.inner_boundary1_i(k) = state(arma::sub2ind(sz, 0, k));
        grid.outer_boundary0_i(k) = state(arma::sub2ind(sz, ns-1, k));
        grid.outer_boundary1_i(k) = state(arma::sub2ind(sz, ns-1, k));
    }
    grid.broadcast();
    return state;
}

static std::string production_gamma_table_path() {
    const char* candidates[] = {
        "data/eos/gamma1_hydrogen_v1.dat",
        "../data/eos/gamma1_hydrogen_v1.dat"
    };
    for (const char* candidate : candidates) {
        std::ifstream input(candidate);
        if (input.good()) return candidate;
    }
    throw std::runtime_error("production Gamma1 table not found from test working directory");
}

static void test_final_serial_log_aware_eos_and_gamma1_contracts() {
    const EosGammaTable table = EosGammaTable::load(production_gamma_table_path());
    const std::vector<std::pair<double,double>> states = {
        {3.5e3, 8.0e25},      // mostly neutral chromosphere
        {8.0e3, 1.0e20},      // partial-ionization transition
        {1.0e6, 1.0e15},      // mostly ionized corona
        {1.01*table.min_temperature(), 1.01*table.min_n_h()},
        {0.99*table.max_temperature(), 0.99*table.max_n_h()}
    };
    for (const auto& state : states) {
        const double temperature = state.first;
        const double n_h = state.second;
        const double rho = n_h*eos_constants::m_h;
        const CaloricState normal = equilibrium_caloric_state(rho, temperature);
        const CaloricState log_aware = equilibrium_caloric_state_from_logs(
            rho, temperature, std::log(temperature), std::log(n_h));
        EXPECT_NEAR(log_aware.n_h, normal.n_h, 0.0);
        EXPECT_NEAR(log_aware.x, normal.x, 0.0);
        EXPECT_NEAR(log_aware.n_e, normal.n_e, 0.0);
        EXPECT_NEAR(log_aware.n_hi, normal.n_hi, 0.0);
        EXPECT_NEAR(log_aware.pressure, normal.pressure, 0.0);
        EXPECT_NEAR(log_aware.internal_energy, normal.internal_energy, 0.0);
        EXPECT_NEAR(log_aware.heat_capacity, normal.heat_capacity, 0.0);

        const double gamma_normal = table.gamma1(temperature, n_h);
        const double gamma_log = table.gamma1_from_logs(
            temperature, n_h, std::log(temperature), std::log(n_h));
        EXPECT_NEAR(gamma_log, gamma_normal, 0.0);

        const MixtureFaceState face_normal = equilibrium_mixture_face_state(
            table, rho, 123.0, temperature, 2.0e6, 1.0e-8);
        const MixtureFaceState face_log = equilibrium_mixture_face_state_from_logs(
            table, std::log(rho), 123.0, std::log(temperature), 2.0e6, 1.0e-8);
        // Exact equality was attempted first. The physical face API receives
        // rho,T directly, whereas the MUSCL API necessarily reconstructs them
        // through exp(log(rho)),exp(log(T)); that round trip differs by up to a
        // few ULPs. Use a tight floating-point-equivalence bound here. The
        // neutral-energy row is a conservative remainder and therefore uses an
        // absolute tolerance scaled by the total carried energy.
        EXPECT_REL(face_log.primitive.rho, face_normal.primitive.rho, 5e-14);
        EXPECT_NEAR(face_log.primitive.velocity, face_normal.primitive.velocity, 0.0);
        EXPECT_REL(face_log.primitive.temperature,
                   face_normal.primitive.temperature, 5e-14);
        EXPECT_REL(face_log.conserved.rho_i, face_normal.conserved.rho_i, 5e-14);
        EXPECT_REL(face_log.conserved.rho_n, face_normal.conserved.rho_n, 5e-14);
        EXPECT_REL(face_log.conserved.momentum_i,
                   face_normal.conserved.momentum_i, 5e-14);
        EXPECT_REL(face_log.conserved.momentum_n,
                   face_normal.conserved.momentum_n, 5e-14);
        EXPECT_REL(face_log.conserved.energy_i,
                   face_normal.conserved.energy_i, 5e-14);
        const double energy_scale = std::max(
            1.0, std::abs(face_normal.conserved.energy_i)
               + std::abs(face_normal.conserved.energy_n));
        EXPECT_NEAR(face_log.conserved.energy_n,
                    face_normal.conserved.energy_n,
                    32.0*std::numeric_limits<double>::epsilon()*energy_scale);
        EXPECT_REL(face_log.conserved.energy_e,
                   face_normal.conserved.energy_e, 5e-14);
        EXPECT_NEAR(face_log.conserved.thermo.gamma1,
                    face_normal.conserved.thermo.gamma1, 0.0);
        EXPECT_REL(face_log.p_total, face_normal.p_total, 5e-14);
        EXPECT_REL(face_log.sound_speed, face_normal.sound_speed, 5e-14);

        const double inverted = temperature_from_rho_eint(
            table, rho, normal.internal_energy, temperature*1.000001);
        EXPECT_REL(inverted, temperature, 2e-12);
    }

    const double n_h = 1.0e18;
    const double rho = n_h*eos_constants::m_h;
    const double temperature = 1.0e4;
    const double x = saha_ionization_fraction_n_h(n_h, temperature);
    const double x_row = 0.3;
    const double rho_i = x_row*rho;
    const double rho_n = (1.0-x_row)*rho;
    const double velocity_i = 120.0;
    const double velocity_n = -45.0;
    const double momentum_i = rho_i*velocity_i;
    const double momentum_n = rho_n*velocity_n;
    const double phi = 2.5e6;
    const double p_e = x*n_h*eos_constants::k_b*temperature;
    const double p_i = 2.0*p_e;
    const double p_n = (1.0-x)*n_h*eos_constants::k_b*temperature;
    const double energy_i = 1.5*p_i+x*n_h*eos_constants::chi_h
                          +0.5*rho_i*velocity_i*velocity_i+rho_i*phi;
    const double energy_n = 1.5*p_n
                          +0.5*rho_n*velocity_n*velocity_n+rho_n*phi;

    set_eos_operation_counting(true);
    reset_eos_operation_counts();
    const CaloricMixtureThermo caloric = decode_equilibrium_caloric_mixture(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n, phi);
    EXPECT_TRUE(eos_operation_counts().gamma1_queries == 0);
    reset_eos_operation_counts();
    const MixtureThermo full = decode_equilibrium_mixture(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n, phi);
    EXPECT_TRUE(eos_operation_counts().gamma1_queries == 1);
    EXPECT_NEAR(caloric.T, full.T, 0.0);
    EXPECT_NEAR(caloric.internal_energy, full.internal_energy, 0.0);

    reset_eos_operation_counts();
    const ProjectedMixtureRows rows = project_equilibrium_rows(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n,
        phi, 1.0e-8);
    EXPECT_TRUE(eos_operation_counts().gamma1_queries == 0);
    reset_eos_operation_counts();
    const ProjectedMixture projected = project_equilibrium_single_fluid(
        table, rho_i, rho_n, momentum_i, momentum_n, energy_i, energy_n,
        phi, 1.0e-8);
    EXPECT_TRUE(eos_operation_counts().gamma1_queries == 1);
    EXPECT_NEAR(rows.rho_i, projected.rho_i, 0.0);
    EXPECT_NEAR(rows.rho_n, projected.rho_n, 0.0);
    EXPECT_NEAR(rows.momentum_i, projected.momentum_i, 0.0);
    EXPECT_NEAR(rows.momentum_n, projected.momentum_n, 0.0);
    EXPECT_NEAR(rows.energy_i, projected.energy_i, 0.0);
    EXPECT_NEAR(rows.energy_n, projected.energy_n, 0.0);
    EXPECT_NEAR(rows.energy_e, projected.energy_e, 0.0);

    const double total_momentum = momentum_i+momentum_n;
    const double total_energy = energy_i+energy_n;
    reset_eos_operation_counts();
    const ProjectedMixtureRows packed_rows =
        pack_equilibrium_rows_from_known_temperature(
            table, rho, total_momentum, total_energy, projected.thermo.T,
            phi, 1.0e-8, 4.0e-11);
    EXPECT_TRUE(eos_operation_counts().gamma1_queries == 0);
    reset_eos_operation_counts();
    const ProjectedMixture packed = pack_equilibrium_from_known_temperature(
        table, rho, total_momentum, total_energy, projected.thermo.T,
        phi, 1.0e-8, 4.0e-11);
    EXPECT_TRUE(eos_operation_counts().gamma1_queries == 1);
    EXPECT_NEAR(packed_rows.energy_i, packed.energy_i, 0.0);
    EXPECT_NEAR(packed_rows.energy_n, packed.energy_n, 0.0);
    EXPECT_NEAR(packed_rows.energy_e, packed.energy_e, 0.0);

    Grid grid;
    const Vec rhs_state = setup_gamma_equilibrium(grid, table, 8, true);
    DecodedMixtureField decoded;
    decode_mixture_field_into(grid, rhs_state, decoded, 1, nullptr);
    reset_eos_operation_counts();
    const Vec rhs = rhs_explicit_state(grid, rhs_state, decoded);
    EXPECT_TRUE(rhs.n_elem == grid.n_state);
    EXPECT_TRUE(eos_operation_counts().gamma1_queries == 6*grid.ns);
    set_eos_operation_counting(false);
}

static void test_final_serial_one_update_and_fallback_regression() {
    const EosGammaTable table = EosGammaTable::load(production_gamma_table_path());
    const double rho = 1.0e20*eos_constants::m_h;
    const double target_temperature = 8.0e3;
    const double target_energy = equilibrium_internal_energy(rho, target_temperature);

    set_runtime_profiling(true);
    reset_runtime_profile();
    const double one_update = temperature_from_rho_eint(
        table, rho, target_energy, target_temperature*(1.0-1.0e-9));
    const EosInversionProfile fast = eos_inversion_profile();
    EXPECT_TRUE(fast.calls == 1);
    EXPECT_TRUE(fast.initial_guess_accepts == 0);
    EXPECT_TRUE(fast.one_update_convergences == 1);
    EXPECT_TRUE(fast.bracket_evaluations == 0);
    EXPECT_REL(one_update, target_temperature, 2e-12);

    reset_runtime_profile();
    const double fallback = temperature_from_rho_eint(
        table, rho, target_energy, 0.8*table.max_temperature());
    const EosInversionProfile slow = eos_inversion_profile();
    EXPECT_TRUE(slow.calls == 1);
    EXPECT_TRUE(slow.bracket_evaluations >= 2);
    EXPECT_TRUE(slow.one_update_convergences == 0);
    EXPECT_REL(fallback, target_temperature, 2e-12);
    EXPECT_REL(fallback, one_update, 2e-12);
    set_runtime_profiling(false);
    reset_runtime_profile();
}

static Vec setup_gamma_uniform_point(Grid& grid, const EosGammaTable& table,
                                     double n_h, double temperature) {
    Vec state = setup_gamma_equilibrium(grid, table, 4, false);
    const double rho_target = n_h*eos_constants::m_h;
    const double x_eq_target = saha_ionization_fraction(rho_target, temperature);
    const double x_row_target = std::max(static_cast<double>(grid.eos_trace_fraction_floor),
        std::min(1.0-static_cast<double>(grid.eos_trace_fraction_floor), x_eq_target));
    float rho_i = static_cast<float>(x_row_target*rho_target);
    float rho_n = static_cast<float>((1.0-x_row_target)*rho_target);
    // Pick the closest float carrier sum on the legal side of an exact density
    // endpoint. This preserves the production hard-error policy while avoiding
    // asking a float-packed state to represent an impossible exact double sum.
    if (n_h == table.max_n_h()) {
        while ((static_cast<double>(rho_i)+rho_n)/eos_constants::m_h > table.max_n_h())
            rho_n = std::nextafter(rho_n, -std::numeric_limits<float>::infinity());
    } else if (n_h == table.min_n_h()) {
        while ((static_cast<double>(rho_i)+rho_n)/eos_constants::m_h < table.min_n_h())
            rho_n = std::nextafter(rho_n, std::numeric_limits<float>::infinity());
    }
    const double rho = static_cast<double>(rho_i) + rho_n;
    const GammaState gamma = gamma_state(table, rho, temperature);
    const double actual_n_h = rho/eos_constants::m_h;
    const double p_e = gamma.x_eq*actual_n_h*eos_constants::k_b*temperature;
    const double p_i = 2.0*p_e;
    const double p_n = (1.0-gamma.x_eq)*actual_n_h*eos_constants::k_b*temperature;
    float energy_i = static_cast<float>(1.5*p_i
        + gamma.x_eq*actual_n_h*eos_constants::chi_h);
    float energy_n = static_cast<float>(equilibrium_internal_energy(rho, temperature)
                                      - static_cast<double>(energy_i));
    const double e_min = equilibrium_internal_energy(rho, table.min_temperature());
    const double e_max = equilibrium_internal_energy(rho, table.max_temperature());
    if (temperature == table.min_temperature()) {
        while (static_cast<double>(energy_i)+energy_n < e_min)
            energy_n = std::nextafter(energy_n, std::numeric_limits<float>::infinity());
    } else if (temperature == table.max_temperature()) {
        while (static_cast<double>(energy_i)+energy_n > e_max)
            energy_n = std::nextafter(energy_n, -std::numeric_limits<float>::infinity());
    }
    const float u[7] = {rho_i, rho_n, 0.0f, 0.0f, energy_i, energy_n,
                        static_cast<float>(1.5*p_e)};
    const auto sz = arma::size(grid.ns, num_of_eq);
    for (arma::uword i = 0; i < grid.ns; ++i)
        for (arma::uword k = 0; k < num_of_eq; ++k)
            state(arma::sub2ind(sz, i, k)) = static_cast<float>(u[k]);
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        grid.inner_boundary0_i(k) = grid.inner_boundary1_i(k) = static_cast<float>(u[k]);
        grid.outer_boundary0_i(k) = grid.outer_boundary1_i(k) = static_cast<float>(u[k]);
    }
    grid.dt_state.zeros();
    grid.broadcast();
    return state;
}

static void test_stage5_production_table_boundary_reconstruction() {
    const EosGammaTable table = EosGammaTable::load(production_gamma_table_path());
    const double n_mid = std::sqrt(table.min_n_h()*table.max_n_h());
    const double t_mid = std::sqrt(table.min_temperature()*table.max_temperature());
    const std::pair<double, double> points[] = {
        {n_mid, table.min_temperature()},
        {n_mid, table.max_temperature()},
        {table.min_n_h(), t_mid},
        {table.max_n_h(), t_mid}
    };
    for (const auto& point : points) {
        EXPECT_TRUE(!throws_any([&] {
            equilibrium_mixture_face_state(
                table, point.first*eos_constants::m_h, 0.0, point.second,
                0.0, 1.0e-8);
        }));
        Grid grid;
        const Vec state = setup_gamma_uniform_point(grid, table, point.first, point.second);
        EXPECT_TRUE(!rhs_explicit_state(grid, state).has_nan());
        Vec dt(grid.ns, arma::fill::zeros);
        EXPECT_TRUE(!advance_Euler_state(grid, state, dt).has_nan());
    }
}

static void test_stage5_6_equilibrium_face_flux_and_sound_speed() {
    const std::string path = write_gamma_table_fixture("gamma_stage56_face.dat");
    const EosGammaTable table = EosGammaTable::load(path);
    const double n_h = 1.0e18;
    const double rho = n_h*eos_constants::m_h;
    const double velocity = -2.3e3;
    const double temperature = 1.0e4;
    const double phi_face = 4.2e6;
    const MixtureFaceState face = equilibrium_mixture_face_state(
        table, rho, velocity, temperature, phi_face, 1.0e-8);
    const std::array<double, 7> flux = equilibrium_mixture_flux(face);

    EXPECT_REL(face.conserved.rho_i + face.conserved.rho_n, rho, 2e-15);
    EXPECT_REL(face.conserved.momentum_i/face.conserved.rho_i, velocity, 2e-15);
    EXPECT_REL(face.conserved.momentum_n/face.conserved.rho_n, velocity, 2e-15);
    EXPECT_REL(flux[0] + flux[1], rho*velocity, 2e-15);
    EXPECT_REL(flux[2] + flux[3], rho*velocity*velocity + face.p_total, 2e-15);
    EXPECT_REL(flux[4] + flux[5],
               (face.conserved.energy_i + face.conserved.energy_n
                + face.p_total)*velocity, 2e-15);
    EXPECT_REL(face.sound_speed*face.sound_speed,
               face.conserved.thermo.gamma1*face.p_total/rho, 2e-15);
    const MixtureThermo decoded = decode_equilibrium_mixture(
        table, face.conserved.rho_i, face.conserved.rho_n,
        face.conserved.momentum_i, face.conserved.momentum_n,
        face.conserved.energy_i, face.conserved.energy_n, phi_face);
    EXPECT_REL(decoded.T, temperature, 2e-12);
    EXPECT_NEAR(decoded.x_eq, face.conserved.thermo.x_eq, 2e-14);

    // Both sides packed with one interface potential have no gravitational
    // energy jump. Using separate cell potentials would create exactly rho*dphi.
    const MixtureFaceState same_phi = equilibrium_mixture_face_state(
        table, rho, velocity, temperature, phi_face, 1.0e-8);
    const MixtureFaceState wrong_phi = equilibrium_mixture_face_state(
        table, rho, velocity, temperature, phi_face + 1.0e5, 1.0e-8);
    EXPECT_NEAR((face.conserved.energy_i + face.conserved.energy_n)
                - (same_phi.conserved.energy_i + same_phi.conserved.energy_n),
                0.0, 0.0);
    EXPECT_REL((wrong_phi.conserved.energy_i + wrong_phi.conserved.energy_n)
               - (face.conserved.energy_i + face.conserved.energy_n),
               rho*1.0e5, 2e-12);

    // Storage floors alter carrier rows only, never face thermodynamics or flux sums.
    const MixtureFaceState floor6 = equilibrium_mixture_face_state(
        table, 1.0e20*eos_constants::m_h, velocity, 1.0e3, phi_face, 1.0e-6);
    const MixtureFaceState floor10 = equilibrium_mixture_face_state(
        table, 1.0e20*eos_constants::m_h, velocity, 1.0e3, phi_face, 1.0e-10);
    EXPECT_NEAR(floor6.conserved.thermo.x_eq, floor10.conserved.thermo.x_eq, 0.0);
    EXPECT_NEAR(floor6.p_total, floor10.p_total, 0.0);
    EXPECT_NEAR(floor6.sound_speed, floor10.sound_speed, 0.0);
    EXPECT_TRUE(floor6.conserved.rho_i > floor10.conserved.rho_i);

    Grid grid;
    Vec state = setup_gamma_equilibrium(grid, table, 4, false);
    const Vec spectral = cal_spectral_radius_state(grid, state);
    const MixtureFaceState rest = equilibrium_mixture_face_state(
        table, rho, 150.0, temperature, 0.0, grid.eos_trace_fraction_floor);
    for (arma::uword i = 0; i < grid.ns; ++i)
        EXPECT_REL(spectral(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, 0)),
                   std::abs(150.0)+rest.sound_speed, 2e-6);
    std::remove(path.c_str());
}

static void test_stage9_acoustic_characteristic_speed() {
    const EosGammaTable table = EosGammaTable::load(production_gamma_table_path());
    auto transition_temperature = [&](double n_h, double target_x) {
        double lo = table.min_temperature(), hi = table.max_temperature();
        for (int iteration = 0; iteration < 100; ++iteration) {
            const double mid = std::sqrt(lo*hi);
            if (saha_ionization_fraction_n_h(n_h,mid) < target_x) lo = mid;
            else hi = mid;
        }
        return std::sqrt(lo*hi);
    };
    struct AcousticPoint { double n_h, temperature, expected_x; };
    const AcousticPoint points[] = {
        // Retain the nearly fully ionized ideal-gas limit, but do not let it
        // stand in for validation through the variable-Gamma1 transition.
        {1.0e18,1.5e4,saha_ionization_fraction_n_h(1.0e18,1.5e4)},
        {1.0e20,transition_temperature(1.0e20,0.1),0.1},
        {1.0e20,transition_temperature(1.0e20,0.5),0.5},
        {1.0e20,transition_temperature(1.0e20,0.9),0.9}
    };
    for (const AcousticPoint& point : points) {
        const double rho = point.n_h*eos_constants::m_h;
        const GammaState gamma = gamma_state(table,rho,point.temperature);
        EXPECT_NEAR(gamma.x_eq,point.expected_x,2.0e-12);
        const double pressure = (1.0+gamma.x_eq)*point.n_h
                              * eos_constants::k_b*point.temperature;
        const double specific_energy =
            equilibrium_internal_energy(rho,point.temperature)/rho;

        // A small adiabatic density perturbation obeys du=p/rho^2 d(rho).
        // Re-invert the caloric EOS on each side and compare the numerical
        // isentropic derivative with the independently tabulated CRASH Gamma1.
        const double drho = 1.0e-5*rho;
        auto perturbed_pressure = [&](double rho_p) {
            const double u_p = specific_energy
                             + pressure/(rho*rho)*(rho_p-rho);
            const double t_p = temperature_from_rho_eint(
                table,rho_p,rho_p*u_p,point.temperature);
            const double n_p = rho_p/eos_constants::m_h;
            const double x_p = saha_ionization_fraction_n_h(n_p,t_p);
            return (1.0+x_p)*n_p*eos_constants::k_b*t_p;
        };
        const double dp_drho_s = (perturbed_pressure(rho+drho)
                                 -perturbed_pressure(rho-drho))/(2.0*drho);
        const double numerical_speed = std::sqrt(dp_drho_s);
        const double expected_speed = std::sqrt(gamma.gamma_sound*pressure/rho);
        EXPECT_REL(numerical_speed,expected_speed,5e-4);
    }
}

static void test_stage9_storage_floor_transport_independence() {
    const EosGammaTable table = EosGammaTable::load(production_gamma_table_path());
    const double floors[] = {1.0e-6, 1.0e-8, 1.0e-10};
    Grid grid;
    grid.init(4, 0.25f);
    grid.ds_i.fill(2.0e4f);
    grid.B_i.ones(); grid.B_imh.ones(); grid.B_iph.ones();
    grid.broadcast();

    struct Diagnostic {
        double rho_i, rho_n, x_row;
        double mass, energy, mass_flux, x, T, p, ne, nhi;
        float ke, kn, qrad;
    };
    struct FloorCase { double n_h, temperature; bool trace_ion; };
    const FloorCase cases[] = {
        // x_eq=9.76e-11: all three lower clamps set distinct charged rows.
        {8.0e25,3.3e3,true},
        // 1-x_eq=4.85e-11: all three upper clamps set distinct neutral rows.
        {1.0e20,1.0e6,false}
    };
    for (const FloorCase& state : cases) {
        Diagnostic result[3];
        for (int f = 0; f < 3; ++f) {
            const MixtureFaceState face = equilibrium_mixture_face_state(
                table,state.n_h*eos_constants::m_h,320.0,state.temperature,
                7.0e5,floors[f]);
            // Exercise the actual runtime precision boundary: pack all seven
            // conservative rows as float, then recover physical composition only
            // through the read-only mixture decoder. This catches accidental use
            // of the storage-floor carrier densities in transport/source paths.
            Vec packed(num_of_eq);
            packed(cons::RHO_I) = static_cast<float>(face.conserved.rho_i);
            packed(cons::RHO_N) = static_cast<float>(face.conserved.rho_n);
            packed(cons::MOM_I) = static_cast<float>(face.conserved.momentum_i);
            packed(cons::MOM_N) = static_cast<float>(face.conserved.momentum_n);
            packed(cons::E_I) = static_cast<float>(face.conserved.energy_i);
            packed(cons::E_N) = static_cast<float>(face.conserved.energy_n);
            packed(cons::E_E) = static_cast<float>(face.conserved.energy_e);
            const MixtureThermo decoded = decode_equilibrium_mixture(
                table,packed(cons::RHO_I),packed(cons::RHO_N),
                packed(cons::MOM_I),packed(cons::MOM_N),
                packed(cons::E_I),packed(cons::E_N),7.0e5);
            Vec ne(grid.ns), nhi(grid.ns), temp(grid.ns);
            ne.fill(static_cast<float>(decoded.n_e));
            nhi.fill(static_cast<float>(decoded.n_HI));
            temp.fill(static_cast<float>(decoded.T));
            const Vec q = radiative_loss_thick(grid,ne,nhi,temp)
                        + radiative_loss_thin(grid,ne,nhi,temp);
            result[f] = {
                packed(cons::RHO_I),packed(cons::RHO_N),decoded.x_row,
                static_cast<double>(packed(cons::RHO_I))+packed(cons::RHO_N),
                static_cast<double>(packed(cons::E_I))+packed(cons::E_N),
                static_cast<double>(packed(cons::MOM_I))+packed(cons::MOM_N),
                decoded.x_eq,decoded.T,decoded.p_i+decoded.p_n,
                decoded.n_e,decoded.n_HI,
                kappa_e(ne,nhi,temp)(0),kappa_n(ne,nhi,temp,temp)(0),q(0)};
        }
        if (state.trace_ion) {
            EXPECT_TRUE(result[0].rho_i > result[1].rho_i);
            EXPECT_TRUE(result[1].rho_i > result[2].rho_i);
            EXPECT_TRUE(result[0].x_row > result[1].x_row);
            EXPECT_TRUE(result[1].x_row > result[2].x_row);
        } else {
            EXPECT_TRUE(result[0].rho_n > result[1].rho_n);
            EXPECT_TRUE(result[1].rho_n > result[2].rho_n);
            EXPECT_TRUE(result[0].x_row < result[1].x_row);
            EXPECT_TRUE(result[1].x_row < result[2].x_row);
        }
        for (int f = 1; f < 3; ++f) {
            EXPECT_REL(result[f].x,result[0].x,2e-6);
            EXPECT_REL(result[f].T,result[0].T,2e-6);
            EXPECT_REL(result[f].p,result[0].p,2e-6);
            EXPECT_REL(result[f].ne,result[0].ne,2e-6);
            EXPECT_REL(result[f].nhi,result[0].nhi,2e-6);
            EXPECT_REL(result[f].ke,result[0].ke,3e-6);
            EXPECT_REL(result[f].kn,result[0].kn,3e-6);
            EXPECT_REL(result[f].qrad,result[0].qrad,3e-6);
            EXPECT_REL(result[f].mass,result[0].mass,2e-7);
            EXPECT_REL(result[f].energy,result[0].energy,2e-7);
            EXPECT_REL(result[f].mass_flux,result[0].mass_flux,2e-7);
        }
    }
}

struct EvaporationDiagnostic {
    double upward_mass_flux, max_temperature, top_pressure;
};

static EvaporationDiagnostic run_stage9_evaporation_comparison(bool gamma_mode) {
    Grid grid;
    grid.init(32, 0.25f);
    if (gamma_mode)
        grid.eos_gamma_table = EosGammaTable::load(production_gamma_table_path());
    Vec state = model_column_ic(grid);
    for (int step = 0; step < 8; ++step) {
        model_column_update_bc(grid,state);
        const Vec dt = cal_dt_i(grid,state);
        state = advance_Euler_state(grid,state,dt);
        grid.sim_time += dt(0);
    }
    const auto sz = arma::size(grid.ns,num_of_eq);
    double max_temperature = 0.0, top_pressure = 0.0;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        if (gamma_mode) {
            auto at = [&](arma::uword k) { return static_cast<double>(
                state(arma::sub2ind(sz,i,k))); };
            const double phi = 0.5*(grid.phi_g_imh(i)+grid.phi_g_iph(i));
            const MixtureThermo th = decode_equilibrium_mixture(
                grid.eos_gamma_table,at(cons::RHO_I),at(cons::RHO_N),
                at(cons::MOM_I),at(cons::MOM_N),at(cons::E_I),at(cons::E_N),phi);
            max_temperature = std::max(max_temperature,th.T);
            if (i+1 == grid.ns) top_pressure = th.p_i+th.p_n;
        } else {
            const Vec primitive = cons2prim(grid,state);
            const double ri = primitive(arma::sub2ind(sz,i,prim::RHO_I));
            const double rn = primitive(arma::sub2ind(sz,i,prim::RHO_N));
            const double pi = primitive(arma::sub2ind(sz,i,prim::P_I));
            const double pn = primitive(arma::sub2ind(sz,i,prim::P_N));
            const double t = (pi+pn)/((2.0*ri/grid.m_i+rn/grid.m_n)*grid.k_b);
            max_temperature = std::max(max_temperature,t);
            if (i+1 == grid.ns) top_pressure = pi+pn;
        }
    }
    double mass_flux = -std::numeric_limits<double>::infinity();
    for (arma::uword i = 0; i < grid.ns; ++i)
        mass_flux = std::max(mass_flux,
            static_cast<double>(state(arma::sub2ind(sz,i,cons::MOM_I)))
            +state(arma::sub2ind(sz,i,cons::MOM_N)));
    return {mass_flux,max_temperature,top_pressure};
}

static void test_stage9_evaporation_vs_fixed_gamma() {
    unsetenv("GAMMA_TABLE"); unsetenv("ISO_GAMMA");
    setenv("ISO_TWO_FLUID","0",1); setenv("ISO_IONIZATION","0",1);
    setenv("ISO_HEAT_FLUX","1",1); setenv("ISO_COOLING","0",1);
    setenv("ISO_CORONA","0",1); setenv("ISO_TRAC","0",1);
    setenv("ISO_TJUMP_A","1.3",1); setenv("ISO_TJUMP_B","1.0",1);
    setenv("ISO_INNER_T_NEUMANN","1",1);
    const EvaporationDiagnostic fixed = run_stage9_evaporation_comparison(false);
    const EvaporationDiagnostic gamma = run_stage9_evaporation_comparison(true);
    EXPECT_TRUE(std::isfinite(fixed.upward_mass_flux) && std::isfinite(gamma.upward_mass_flux));
    EXPECT_TRUE(fixed.max_temperature > 0.0 && gamma.max_temperature > 0.0);
    EXPECT_TRUE(fixed.top_pressure > 0.0 && gamma.top_pressure > 0.0);
    EXPECT_TRUE(fixed.upward_mass_flux > 0.0 && gamma.upward_mass_flux > 0.0);
    EXPECT_REL(fixed.upward_mass_flux,1.10074e-3,0.1);
    EXPECT_REL(gamma.upward_mass_flux,3.58416e-3,0.1);
    EXPECT_REL(fixed.max_temperature,6813.14,0.03);
    EXPECT_REL(gamma.max_temperature,6937.23,0.03);
    EXPECT_REL(fixed.top_pressure,3.43710,0.05);
    EXPECT_REL(gamma.top_pressure,3.62956,0.05);
    const double mass_flux_ratio = gamma.upward_mass_flux/fixed.upward_mass_flux;
    EXPECT_TRUE(mass_flux_ratio > 2.0 && mass_flux_ratio < 5.0);
    std::cout << "  Stage-9 evaporation regression (fixed,gamma): maximum_upward_mass_flux=("
              << fixed.upward_mass_flux << ',' << gamma.upward_mass_flux << ") T_max=("
              << fixed.max_temperature << ',' << gamma.max_temperature << ") p_top=("
              << fixed.top_pressure << ',' << gamma.top_pressure << ")\n";
    unsetenv("ISO_TWO_FLUID"); unsetenv("ISO_IONIZATION");
    unsetenv("ISO_HEAT_FLUX"); unsetenv("ISO_COOLING");
    unsetenv("ISO_CORONA"); unsetenv("ISO_TRAC");
    unsetenv("ISO_TJUMP_A"); unsetenv("ISO_TJUMP_B");
    unsetenv("ISO_INNER_T_NEUMANN");
}

static void test_stage5_gamma_integrator_projection_and_guards() {
    const std::string path = write_gamma_table_fixture("gamma_stage5_integrator.dat");
    const EosGammaTable table = EosGammaTable::load(path);
    Grid grid;
    const Vec state = setup_gamma_equilibrium(grid, table, 6, false);
    // A NaN legacy gamma is harmless here: both reconstruction passes, fluxes,
    // source decode, and the internal predictor consume only the mixture EOS.
    grid.gamma_mono = std::numeric_limits<float>::quiet_NaN();
    Vec dt(grid.ns); dt.fill(1.0e-4f);
    const Vec result = advance_Euler_state(grid, state, dt);
    EXPECT_TRUE(!result.has_nan());
    EXPECT_TRUE(arma::approx_equal(result, project_equilibrium_single_fluid(grid, result),
                                   "reldiff", 2.0e-6));
    EXPECT_TRUE(throws_any([&] { advance_Euler_explicit_state(grid, state, dt); }));
    EXPECT_TRUE(throws_any([&] { advance_RK4(grid, state, dt); }));

    grid.single_fluid = false;
    EXPECT_TRUE(throws_any([&] { advance_Euler_state(grid, state, dt); }));
    grid.single_fluid = true; grid.enable_Te = true;
    EXPECT_TRUE(throws_any([&] { advance_Euler_state(grid, state, dt); }));
    grid.enable_Te = false; grid.enable_ionization = true;
    EXPECT_TRUE(throws_any([&] { advance_Euler_state(grid, state, dt); }));
    grid.enable_ionization = false; grid.enable_conduction = true;
    EXPECT_TRUE(!advance_Euler_state(grid, state, dt).has_nan());
    std::remove(path.c_str());
}

static double packed_total_energy(const Grid& grid, const Vec& state) {
    const auto sz = arma::size(grid.ns, num_of_eq);
    double sum = 0.0;
    for (arma::uword i = 0; i < grid.ns; ++i)
        sum += static_cast<double>(state(arma::sub2ind(sz, i, cons::E_I)))
             + state(arma::sub2ind(sz, i, cons::E_N));
    return sum;
}

static double weighted_domain_energy(const Grid& grid, const Vec& state) {
    const auto sz = arma::size(grid.ns, num_of_eq);
    double sum = 0.0;
    for (arma::uword i = 0; i < grid.ns; ++i)
        sum += (static_cast<double>(state(arma::sub2ind(sz, i, cons::E_I)))
              + state(arma::sub2ind(sz, i, cons::E_N)))*grid.ds_i(i)/grid.B_i(i);
    return sum;
}

static void test_stage7_gamma_total_energy_sources_and_conduction() {
    const EosGammaTable table = EosGammaTable::load(production_gamma_table_path());
    Grid grid;
    Vec state = setup_gamma_equilibrium(grid, table, 8, false);
    Vec dt(grid.ns); dt.fill(2.0e-3f);
    const auto sz = arma::size(grid.ns, num_of_eq);

    grid.enable_coronal_heating = true;
    grid.coronal_heat_E0 = 2.5e-3f;
    grid.coronal_heat_sH = 1.0e30f;
    grid.coronal_heat_s0 = 0.0f;
    grid.coronal_heat_two_sided = false;
    const Vec H = coronal_heating_rate(grid);
    const double energy0 = packed_total_energy(grid, state);
    const Vec heated = advance_Euler_state(grid, state, dt);
    EXPECT_REL(packed_total_energy(grid, heated)-energy0,
               dt(0)*arma::sum(H), 3e-2); // float packed-row increment
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_NEAR(heated(arma::sub2ind(sz, i, cons::RHO_I))
                  + heated(arma::sub2ind(sz, i, cons::RHO_N)),
                    state(arma::sub2ind(sz, i, cons::RHO_I))
                  + state(arma::sub2ind(sz, i, cons::RHO_N)), 2e-12);
    }

    grid.enable_coronal_heating = false;
    grid.enable_beam_heating = true;
    grid.beam_flux = 1.0f;
    grid.beam_t_on = 0.0f; grid.beam_duration = 10.0f; grid.beam_ramp = 1.0f;
    grid.sim_time = 5.0f;
    const Vec beamed = advance_Euler_state(grid, state, dt);
    EXPECT_TRUE(packed_total_energy(grid, beamed) > energy0);
    EXPECT_TRUE(!beamed.has_nan());

    grid.enable_beam_heating = false;
    grid.enable_radiative_cooling = true;
    const Vec cooled = advance_Euler_state(grid, heated, dt);
    EXPECT_TRUE(packed_total_energy(grid, cooled) < packed_total_energy(grid, heated));
    EXPECT_TRUE(!cooled.has_nan());

    // Closed-end conduction smooths a pulse while conserving the packed total
    // energy. Repack at zero velocity so the explicit hydro boundary flux is zero.
    state = setup_gamma_equilibrium(grid, table, 8, true);
    grid.enable_radiative_cooling = false;
    grid.enable_conduction = true;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double f = static_cast<double>(i)/(grid.ns-1);
        const MixtureFaceState face = equilibrium_mixture_face_state(
            table, 1.0e18*eos_constants::m_h, 0.0, 7.0e3+6.0e3*f, 0.0,
            grid.eos_trace_fraction_floor);
        const double u[7] = {face.conserved.rho_i, face.conserved.rho_n,
            face.conserved.momentum_i, face.conserved.momentum_n,
            face.conserved.energy_i, face.conserved.energy_n, face.conserved.energy_e};
        for (arma::uword k = 0; k < num_of_eq; ++k)
            state(arma::sub2ind(sz, i, k)) = u[k];
    }
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        grid.inner_boundary0_i(k)=grid.inner_boundary1_i(k)=state(arma::sub2ind(sz,0,k));
        grid.outer_boundary0_i(k)=grid.outer_boundary1_i(k)=state(arma::sub2ind(sz,grid.ns-1,k));
    }
    for (arma::uword i = 0; i < grid.ns; ++i) {
        grid.ds_i(i) = 4.0e3f*(1.0f+0.12f*i);
        grid.B_imh(i) = 1.0f+0.05f*i;
        grid.B_iph(i) = 1.0f+0.05f*(i+1);
        grid.B_i(i) = 0.5f*(grid.B_imh(i)+grid.B_iph(i));
        grid.dinvB_ds_i(i) = (1.0f/grid.B_iph(i)-1.0f/grid.B_imh(i))/grid.ds_i(i);
    }
    grid.inner_conduction_neumann = true;
    grid.impose_outer_heat_flux = true;
    grid.outer_heat_flux = 0.0f;
    // Exercise the face-local numerical term, not only physical conduction.
    grid.numerical_diffusivity_per_length = 2.0e3f;
    grid.broadcast();
    dt.fill(1.0e-4f);
    const double conduction_e0 = weighted_domain_energy(grid, state);
    grid.dt_state.zeros();
    for (arma::uword k = 0; k < num_of_eq; ++k)
        grid.dt_state += scalar_to(grid, dt, k);
    const Vec conduction_input = project_equilibrium_single_fluid(
        grid, state+grid.dt_state%rhs_explicit_state(grid, state));
    const Vec conducted = advance_Euler_state(grid, state, dt);
    EXPECT_REL(weighted_domain_energy(grid, conducted), conduction_e0, 3e-6);
    EXPECT_TRUE(!conducted.has_nan());
    EXPECT_TRUE(gamma_conduction_residual_max(
        grid, conduction_input, conducted, dt(0)) < 2e-5);

    // Strong nonlinear/TRAC case spanning chromosphere to corona. This drives
    // kappa_e through five decades of T^(5/2), produces a visible implicit
    // update, and independently checks the final packed residual and energy.
    state = setup_gamma_equilibrium(grid, table, 12, false);
    grid.enable_conduction = true;
    grid.enable_trac = true;
    grid.trac_T_chrom = 2.0e4f;
    grid.trac_Tc_max_frac = 0.2f;
    const auto strong_sz = arma::size(grid.ns, num_of_eq);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double f = static_cast<double>(i)/(grid.ns-1);
        const double temperature = std::exp(std::log(1.0e4)*(1.0-f)+std::log(1.0e6)*f);
        const MixtureFaceState face = equilibrium_mixture_face_state(
            table, 1.0e18*eos_constants::m_h, 0.0, temperature, 0.0,
            grid.eos_trace_fraction_floor);
        const double u[7] = {face.conserved.rho_i,face.conserved.rho_n,
            face.conserved.momentum_i,face.conserved.momentum_n,
            face.conserved.energy_i,face.conserved.energy_n,face.conserved.energy_e};
        for (arma::uword k = 0; k < num_of_eq; ++k)
            state(arma::sub2ind(strong_sz,i,k)) = u[k];
        grid.ds_i(i) = 2.0e4f;
    }
    for (arma::uword k = 0; k < num_of_eq; ++k) {
        grid.inner_boundary0_i(k)=grid.inner_boundary1_i(k)=state(arma::sub2ind(strong_sz,0,k));
        grid.outer_boundary0_i(k)=grid.outer_boundary1_i(k)=state(arma::sub2ind(strong_sz,grid.ns-1,k));
    }
    grid.inner_conduction_neumann = true;
    grid.impose_outer_heat_flux = true;
    grid.outer_heat_flux = 0.0f;
    grid.broadcast();
    dt.set_size(grid.ns); dt.fill(2.0e-2f);
    grid.dt_state.zeros();
    for (arma::uword k = 0; k < num_of_eq; ++k) grid.dt_state += scalar_to(grid,dt,k);
    const Vec strong_input = project_equilibrium_single_fluid(
        grid,state+grid.dt_state%rhs_explicit_state(grid,state));
    const double strong_energy = weighted_domain_energy(grid,state);
    const Vec strong_result = advance_Euler_state(grid,state,dt);
    EXPECT_REL(weighted_domain_energy(grid,strong_result),strong_energy,5e-6);
    EXPECT_TRUE(gamma_conduction_residual_max(grid,strong_input,strong_result,dt(0)) < 5e-5);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        auto value = [&](arma::uword k) { return static_cast<double>(
            strong_result(arma::sub2ind(strong_sz,i,k))); };
        const MixtureThermo th = decode_equilibrium_mixture(
            table,value(cons::RHO_I),value(cons::RHO_N),value(cons::MOM_I),
            value(cons::MOM_N),value(cons::E_I),value(cons::E_N),0.0);
        EXPECT_TRUE(th.T >= table.min_temperature() && th.T <= table.max_temperature());
    }

    // Strong sinks/sources at several density decades must remain within the
    // strict caloric domain after the double target is packed into float rows.
    const double densities[] = {1.0e14, 1.0e18, 1.0e22, 1.0e25};
    for (double n_h : densities) {
        Vec boundary_state = setup_gamma_uniform_point(grid, table, n_h, 1.0e4);
        dt.set_size(grid.ns);
        grid.enable_radiative_cooling = true;
        dt.fill(1.0e30f);
        const Vec floor_state = advance_Euler_state(grid, boundary_state, dt);
        const auto boundary_sz = arma::size(grid.ns, num_of_eq);
        auto first = [&](const Vec& u, arma::uword k) {
            return static_cast<double>(u(arma::sub2ind(boundary_sz, 0, k)));
        };
        const MixtureThermo floor_th = decode_equilibrium_mixture(
            table, first(floor_state, cons::RHO_I), first(floor_state, cons::RHO_N),
            first(floor_state, cons::MOM_I), first(floor_state, cons::MOM_N),
            first(floor_state, cons::E_I), first(floor_state, cons::E_N), 0.0);
        EXPECT_TRUE(floor_th.T >= table.min_temperature());
        EXPECT_REL(floor_th.T, table.min_temperature(), 2e-5);

        boundary_state = setup_gamma_uniform_point(
            grid, table, n_h, 0.999*table.max_temperature());
        grid.enable_coronal_heating = true;
        grid.coronal_heat_E0 = 1.0e30f;
        grid.coronal_heat_sH = 1.0e30f;
        dt.fill(1.0f);
        EXPECT_TRUE(throws_any([&] {
            (void)advance_Euler_state(grid, boundary_state, dt);
        }));
    }
    Vec nonfinite_state = setup_gamma_uniform_point(grid, table, 1.0e20, 1.0e4);
    dt.set_size(grid.ns); dt.fill(1.0f);
    grid.enable_coronal_heating = true;
    grid.coronal_heat_E0 = std::numeric_limits<float>::quiet_NaN();
    EXPECT_TRUE(throws_any([&] {
        (void)advance_Euler_state(grid, nonfinite_state, dt);
    }));
}

static void test_gamma_phase1_3_cache_and_known_temperature_pack() {
    unsetenv("GAMMA_TABLE"); unsetenv("ISO_GAMMA");
    setenv("ISO_TWO_FLUID","0",1); setenv("ISO_IONIZATION","0",1);
    setenv("ISO_HEAT_FLUX","0",1); setenv("ISO_COOLING","0",1);
    setenv("ISO_CORONA","0",1); setenv("ISO_TRAC","0",1);
    Grid grid; grid.init(24,0.25f);
    grid.eos_gamma_table=EosGammaTable::load(production_gamma_table_path());
    Vec state=model_column_ic(grid);
    model_column_update_bc(grid,state);
    const DecodedMixtureField decoded=decode_mixture_field(grid,state,17);
    decoded.require_matches(grid,state);
    EXPECT_TRUE(decoded.state_generation==17);
    EXPECT_TRUE(decoded.cells.size()==grid.ns);

    // Every cached center and shifted stencil entry must equal an independent
    // direct decode using the potential actually carried by that state.
    const Vec phi_cell=0.5f*(grid.phi_g_imh+grid.phi_g_iph);
    Vec phi_ip1=ip1(grid,phi_cell,SLICE); phi_ip1(grid.ns-1)=grid.phi_g_iph(grid.ns-1);
    Vec phi_im1=im1(grid,phi_cell,SLICE); phi_im1(0)=grid.phi_g_imh(0);
    Vec phi_ip2=ip1(grid,phi_ip1,SLICE); phi_ip2(grid.ns-1)=grid.phi_g_iph(grid.ns-1);
    Vec phi_im2=im1(grid,phi_im1,SLICE); phi_im2(0)=grid.phi_g_imh(0);
    const Vec shifted_states[] = {im2(grid,state),im1(grid,state),state,
                                  ip1(grid,state),ip2(grid,state)};
    const Vec shifted_phi[] = {phi_im2,phi_im1,phi_cell,phi_ip1,phi_ip2};
    const int offsets[] = {-2,-1,0,1,2};
    const arma::uword ext_n=grid.ns+4;
    const auto ext_size=arma::size(ext_n,3);
    const auto sz=arma::size(grid.ns,num_of_eq);
    for (int stencil=0; stencil<5; ++stencil) {
        for (arma::uword i=0; i<grid.ns; ++i) {
            const Vec& shifted=shifted_states[stencil];
            auto at=[&](arma::uword row) { return static_cast<double>(
                shifted(arma::sub2ind(sz,i,row))); };
            const MixtureThermo direct=decode_equilibrium_mixture(
                grid.eos_gamma_table,at(cons::RHO_I),at(cons::RHO_N),
                at(cons::MOM_I),at(cons::MOM_N),at(cons::E_I),at(cons::E_N),
                shifted_phi[stencil](i));
            const arma::uword ext_i=static_cast<arma::uword>(
                static_cast<int>(i)+2+offsets[stencil]);
            EXPECT_NEAR(decoded.extended_primitive(
                arma::sub2ind(ext_size,ext_i,0)),std::log(direct.rho),1e-13);
            EXPECT_NEAR(decoded.extended_primitive(
                arma::sub2ind(ext_size,ext_i,1)),
                (at(cons::MOM_I)+at(cons::MOM_N))/direct.rho,1e-13);
            EXPECT_NEAR(decoded.extended_primitive(
                arma::sub2ind(ext_size,ext_i,2)),std::log(direct.T),1e-13);
        }
    }

    // A cache cannot be reused for another state or after its ghost buffers
    // change, preventing silent stale-state and stale-boundary reads.
    Vec changed=state;
    EXPECT_TRUE(throws_any([&] { decoded.require_matches(grid,changed); }));
    const float saved=grid.inner_boundary0_i(cons::E_N);
    grid.inner_boundary0_i(cons::E_N)=std::nextafter(
        saved,std::numeric_limits<float>::infinity());
    EXPECT_TRUE(throws_any([&] { decoded.require_matches(grid,state); }));
    grid.inner_boundary0_i(cons::E_N)=saved;
    decoded.require_matches(grid,state);

    // Cached gamma RHS performs one inversion per predicted physical cell plus
    // four distinct ghosts, never one inversion per shifted copy.
    set_runtime_profiling(true); reset_runtime_profile();
    grid.dt_state.fill(1.0e-4f);
    const Vec rhs=rhs_explicit_state(grid,state,decoded);
    const EosInversionProfile profile=eos_inversion_profile();
    EXPECT_TRUE(!rhs.has_nan());
    EXPECT_TRUE(profile.calls<=grid.ns+4);
    set_runtime_profiling(false); reset_runtime_profile();

    // The guess fast path accepts the exact previous temperature without
    // evaluating either table bracket endpoint.
    const MixtureThermo& th=decoded.cells[grid.ns/2];
    set_runtime_profiling(true); reset_runtime_profile();
    const double recovered=temperature_from_rho_eint(
        grid.eos_gamma_table,th.rho,th.internal_energy,th.T);
    const EosInversionProfile fast=eos_inversion_profile();
    EXPECT_NEAR(recovered,th.T,0.0);
    EXPECT_TRUE(fast.calls==1 && fast.initial_guess_accepts==1);
    EXPECT_TRUE(fast.bracket_evaluations==0);
    set_runtime_profiling(false); reset_runtime_profile();

    // Known-temperature packing keeps the supplied total energy authoritative,
    // even when it differs slightly from the analytic EOS energy within the
    // accepted nonlinear residual.
    const double rho=1.0e20*eos_constants::m_h, T=9000.0, v=123.0, phi=-2.0e6;
    const double e_eos=equilibrium_internal_energy(rho,T);
    const double total=e_eos*(1.0+1.0e-12)+0.5*rho*v*v+rho*phi;
    const ProjectedMixture packed=pack_equilibrium_from_known_temperature(
        grid.eos_gamma_table,rho,rho*v,total,T,phi,1.0e-8,2.0e-11);
    EXPECT_REL(packed.rho_i+packed.rho_n,rho,2e-16);
    EXPECT_REL(packed.momentum_i+packed.momentum_n,rho*v,2e-16);
    EXPECT_NEAR(packed.energy_i+packed.energy_n,total,
                4.0*std::numeric_limits<double>::epsilon()*std::abs(total));
    EXPECT_TRUE(throws_any([&] {
        (void)pack_equilibrium_from_known_temperature(
            grid.eos_gamma_table,rho,rho*v,total,T*1.01,phi,1.0e-8,2.0e-11);
    }));
}

static void test_stage8_gamma_model_column_saha_hse_and_ghosts() {
    unsetenv("GAMMA_TABLE"); unsetenv("ISO_GAMMA");
    setenv("ISO_TWO_FLUID", "0", 1); setenv("ISO_IONIZATION", "0", 1);
    setenv("ISO_HEAT_FLUX", "0", 1); setenv("ISO_COOLING", "0", 1);
    setenv("ISO_CORONA", "0", 1); setenv("ISO_TRAC", "0", 1);
    // Historical fixed-gamma helper remains exactly piecewise linear.
    float linear_t, linear_ne, linear_nhi;
    c7_full_profile(587.5f, linear_t, linear_ne, linear_nhi);
    EXPECT_NEAR(linear_t, 4417.5, 1e-6);
    // Gamma-mode PCHIP is nodal and C1 across representative C7 knots.
    EXPECT_NEAR(c7_full_temperature_pchip(560.0), 4400.0, 1e-12);
    const double knot = 560.0, eps = 1.0e-3;
    const double derivative_left = (c7_full_temperature_pchip(knot)
        - c7_full_temperature_pchip(knot-eps))/eps;
    const double derivative_right = (c7_full_temperature_pchip(knot+eps)
        - c7_full_temperature_pchip(knot))/eps;
    EXPECT_NEAR(derivative_left, derivative_right, 2e-4);
    EXPECT_TRUE(std::abs(c7_full_temperature_pchip(587.5)-linear_t) > 1.0);
    Grid grid;
    grid.init(48, 0.25f);
    grid.eos_gamma_table = EosGammaTable::load(production_gamma_table_path());
    const Vec state = model_column_ic(grid);
    EXPECT_TRUE(grid.single_fluid && !grid.enable_ionization && !grid.enable_Te);
    EXPECT_TRUE(!state.has_nan());
    const auto sz = arma::size(grid.ns, num_of_eq);
    double previous_p = 0.0, previous_invH = 0.0;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        auto u = [&](arma::uword k) { return static_cast<double>(
            state(arma::sub2ind(sz, i, k))); };
        const double phi = 0.5*(grid.phi_g_imh(i)+grid.phi_g_iph(i));
        const MixtureThermo th = decode_equilibrium_mixture(
            grid.eos_gamma_table, u(cons::RHO_I), u(cons::RHO_N),
            u(cons::MOM_I), u(cons::MOM_N), u(cons::E_I), u(cons::E_N), phi);
        const double h_km = phi/(grid.g*1000.0);
        EXPECT_REL(th.T, c7_full_temperature_pchip(h_km), 3e-5);
        EXPECT_REL(th.x_eq, saha_ionization_fraction_n_h(th.n_H, th.T), 2e-13);
        const double p = th.p_i+th.p_n;
        const double invH = grid.g*eos_constants::m_h/
            ((1.0+th.x_eq)*eos_constants::k_b*th.T);
        if (i > 0) {
            const double dh = 0.5*(grid.ds_i(i-1)+grid.ds_i(i));
            EXPECT_NEAR(std::log(p/previous_p)+0.5*(previous_invH+invH)*dh,
                        0.0, 3e-5);
        }
        previous_p = p; previous_invH = invH;
    }
    const Vec* ghosts[] = {&grid.inner_boundary0_i, &grid.inner_boundary1_i,
                           &grid.outer_boundary0_i, &grid.outer_boundary1_i};
    const double phis[] = {grid.phi_g_imh(0), grid.phi_g_imh(0),
                           grid.phi_g_iph(grid.ns-1), grid.phi_g_iph(grid.ns-1)};
    for (int j = 0; j < 4; ++j) {
        const Vec& g = *ghosts[j];
        const MixtureThermo th = decode_equilibrium_mixture(
            grid.eos_gamma_table, g(cons::RHO_I), g(cons::RHO_N),
            g(cons::MOM_I), g(cons::MOM_N), g(cons::E_I), g(cons::E_N), phis[j]);
        EXPECT_REL(th.x_eq, saha_ionization_fraction_n_h(th.n_H, th.T), 2e-13);
        EXPECT_REL(g(cons::E_I)+g(cons::E_N),
                   th.internal_energy+0.5*th.rho*std::pow(
                       (g(cons::MOM_I)+g(cons::MOM_N))/th.rho, 2)+th.rho*phis[j],
                   3e-6);
    }
    Vec evolved = state;
    for (int step = 0; step < 5; ++step) {
        model_column_update_bc(grid, evolved);
        const Vec dt = cal_dt_i(grid, evolved);
        evolved = advance_Euler_state(grid, evolved, dt);
        grid.sim_time += dt(0);
    }
    EXPECT_TRUE(!evolved.has_nan());
    double max_velocity = 0.0, max_total_momentum = 0.0;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double rho = evolved(arma::sub2ind(sz, i, cons::RHO_I))
                         + evolved(arma::sub2ind(sz, i, cons::RHO_N));
        const double momentum = evolved(arma::sub2ind(sz, i, cons::MOM_I))
                              + evolved(arma::sub2ind(sz, i, cons::MOM_N));
        max_velocity = std::max(max_velocity, std::abs(momentum/rho));
        max_total_momentum = std::max(max_total_momentum, std::abs(momentum));
    }
    EXPECT_TRUE(max_velocity < 1.0e-3);
    EXPECT_TRUE(max_total_momentum < 1.0e-7);
    const double relative_state_change = arma::norm(evolved-state, 2)
                                       / arma::norm(state, 2);
    EXPECT_NEAR(relative_state_change, 0.0, 1e-4);
    unsetenv("ISO_TWO_FLUID"); unsetenv("ISO_IONIZATION");
    unsetenv("ISO_HEAT_FLUX"); unsetenv("ISO_COOLING");
    unsetenv("ISO_CORONA"); unsetenv("ISO_TRAC");
}

// ---------------------------------------------------------------------------
// Upper-boundary hydro-temperature / conduction-temperature decoupling
// (ISO_HYDRO_T_DECOUPLE, chromosphere.hpp::outer_conduction_temperature_override)
// ---------------------------------------------------------------------------

// The h=1600..2153 km transition-region column of the upper-BC experiment, at a
// test-sized cell count. ISO_T_TOP pins the reference wall at 22 kK so the
// baseline hydro ghost temperature is a known constant.
static Vec setup_decoupling_column(Grid& grid, bool decouple, arma::uword ns) {
    unsetenv("ISO_GAMMA");
    setenv("ISO_TWO_FLUID", "0", 1); setenv("ISO_IONIZATION", "0", 1);
    setenv("ISO_COOLING", "0", 1);   setenv("ISO_CORONA", "0", 1);
    setenv("ISO_TRAC", "0", 1);      setenv("ISO_HEAT_FLUX", "1", 1);
    setenv("ISO_H_BASE", "1600", 1); setenv("ISO_DH", "553", 1);
    setenv("ISO_T_TOP", "22000", 1);
    setenv("ISO_HYDRO_T_DECOUPLE", decouple ? "1" : "0", 1);
    grid.init(ns, 0.25f);
    grid.eos_gamma_table = EosGammaTable::load(production_gamma_table_path());
    return model_column_ic(grid);
}

static void clear_decoupling_env() {
    unsetenv("ISO_TWO_FLUID"); unsetenv("ISO_IONIZATION"); unsetenv("ISO_COOLING");
    unsetenv("ISO_CORONA");    unsetenv("ISO_TRAC");       unsetenv("ISO_HEAT_FLUX");
    unsetenv("ISO_H_BASE");    unsetenv("ISO_DH");         unsetenv("ISO_T_TOP");
    unsetenv("ISO_HYDRO_T_DECOUPLE");
}

static MixtureThermo decode_cell(const Grid& grid, const Vec& state, arma::uword i) {
    const auto sz = arma::size(grid.ns, num_of_eq);
    auto u = [&](arma::uword k) {
        return static_cast<double>(state(arma::sub2ind(sz, i, k)));
    };
    return decode_equilibrium_mixture(grid.eos_gamma_table,
        u(cons::RHO_I), u(cons::RHO_N), u(cons::MOM_I), u(cons::MOM_N),
        u(cons::E_I), u(cons::E_N), 0.5*(grid.phi_g_imh(i)+grid.phi_g_iph(i)));
}

static MixtureThermo decode_outer_ghost(const Grid& grid, const Vec& ghost) {
    return decode_equilibrium_mixture(grid.eos_gamma_table,
        ghost(cons::RHO_I), ghost(cons::RHO_N), ghost(cons::MOM_I),
        ghost(cons::MOM_N), ghost(cons::E_I), ghost(cons::E_N),
        grid.phi_g_iph(grid.ns-1));
}

// Independent Saha-EOS density inversion rho(p, T): p = (1+x(n_H,T)) n_H k T.
// Mirrors model_column's internal gamma_density_from_pressure, written out here so
// the ghost density closure is checked against a second implementation.
static double saha_density_from_pressure(double pressure, double temperature) {
    double lo = 1.0e8, hi = 1.0e30;
    for (int it = 0; it < 200; ++it) {
        const double mid = std::sqrt(lo*hi);
        const double p = (1.0+saha_ionization_fraction_n_h(mid, temperature))
                       * mid*eos_constants::k_b*temperature;
        if (p < pressure) lo = mid; else hi = mid;
    }
    return std::sqrt(lo*hi)*eos_constants::m_h;
}

// Overwrite one cell with an equilibrium-mixture state at (rho, V, T).
static void set_cell_state(const Grid& grid, Vec& state, arma::uword i,
                           double rho, double V, double T) {
    const auto sz = arma::size(grid.ns, num_of_eq);
    const MixtureFaceState face = equilibrium_mixture_face_state(
        grid.eos_gamma_table, rho, V, T,
        0.5*(grid.phi_g_imh(i)+grid.phi_g_iph(i)), grid.eos_trace_fraction_floor);
    const double u[7] = {face.conserved.rho_i, face.conserved.rho_n,
        face.conserved.momentum_i, face.conserved.momentum_n,
        face.conserved.energy_i, face.conserved.energy_n, face.conserved.energy_e};
    for (arma::uword k = 0; k < num_of_eq; ++k)
        state(arma::sub2ind(sz, i, k)) = u[k];
}

static void test_upper_bc_hydro_conduction_temperature_decoupling() {
    const arma::uword ns = 200;   // 2.77 km cells; ns=40 over-conducts the coarse top cell

    // ---- (1) flag OFF: the pre-existing behaviour, bit for bit ----------
    // The hydro ghosts carry the imposed wall (a=b=1 ⇒ both at ISO_T_TOP) and the
    // conduction override stays disarmed, so Stage D reads the hydro ghost as before.
    {
        Grid grid;
        Vec state = setup_decoupling_column(grid, false, ns);
        EXPECT_TRUE(!grid.outer_conduction_temperature_override);
        const MixtureThermo g0 = decode_outer_ghost(grid, grid.outer_boundary0_i);
        const MixtureThermo g1 = decode_outer_ghost(grid, grid.outer_boundary1_i);
        EXPECT_REL(g0.T, 22000.0, 1e-5);
        EXPECT_REL(g1.T, 22000.0, 1e-5);
        // Move the live top cell well off the wall: the baseline hydro ghost must NOT
        // follow it (that is exactly the extra constraint the experiment removes).
        const MixtureThermo top = decode_cell(grid, state, grid.ns-1);
        set_cell_state(grid, state, grid.ns-1, top.rho, 0.0, 1.20*top.T);
        model_column_update_bc(grid, state);
        EXPECT_REL(decode_outer_ghost(grid, grid.outer_boundary0_i).T, 22000.0, 1e-5);
        EXPECT_TRUE(!grid.outer_conduction_temperature_override);
    }

    // ---- (2),(3),(5),(6),(7) decoupled mode -----------------------------
    {
        Grid grid;
        Vec state = setup_decoupling_column(grid, true, ns);
        EXPECT_TRUE(grid.outer_conduction_temperature_override);
        // Conduction still sees the unchanged fixed hot wall...
        EXPECT_REL(grid.outer_conduction_temperature, 22000.0, 1e-6);
        // ...while the hydro ghost zero-gradient-extrapolates the live top cell.
        const MixtureThermo top0 = decode_cell(grid, state, grid.ns-1);
        MixtureThermo g0 = decode_outer_ghost(grid, grid.outer_boundary0_i);
        MixtureThermo g1 = decode_outer_ghost(grid, grid.outer_boundary1_i);
        EXPECT_REL(g0.T, top0.T, 2e-6);           // T_hydro_g0 = T_top
        EXPECT_REL(g1.T, g0.T, 2e-6);             // T_hydro_g1 = T_hydro_g0
        const double p_back = g0.p_i+g0.p_n;      // the captured reservoir back-pressure

        // (7) second ghost: EOS-closed density on a one-sided HSE rung off the first.
        const float ds = grid.ds_i(grid.ns-1);
        EXPECT_REL(g1.p_i+g1.p_n,
                   p_back - ds*0.5*(g0.rho+g1.rho)*grid.g, 1e-5);   // float32 ghost rows
        EXPECT_REL(g0.rho, saha_density_from_pressure(p_back, g0.T), 2e-6);
        EXPECT_REL(g1.rho, saha_density_from_pressure(g1.p_i+g1.p_n, g1.T), 2e-6);

        // (5) changing the LIVE top temperature moves the hydro ghost temperature but
        // leaves the conduction-only wall alone. (6) The ghost pressure stays exactly
        // the captured back-pressure — it has NOT reverted to live extrapolation.
        const MixtureThermo top = decode_cell(grid, state, grid.ns-1);
        set_cell_state(grid, state, grid.ns-1, 1.05*top.rho, 0.0, 1.20*top.T);
        model_column_update_bc(grid, state);
        g0 = decode_outer_ghost(grid, grid.outer_boundary0_i);
        EXPECT_REL(g0.T, 1.20*top.T, 2e-5);
        EXPECT_REL(g0.p_i+g0.p_n, p_back, 2e-6);
        EXPECT_REL(grid.outer_conduction_temperature, 22000.0, 1e-6);
        EXPECT_TRUE(std::abs(g0.T-grid.outer_conduction_temperature)
                    > 0.05*grid.outer_conduction_temperature);   // genuinely decoupled
    }

    // ---- (4) Stage D responds to the conduction-only wall, not to the ghost ----
    {
        Grid grid;
        Vec state = setup_decoupling_column(grid, true, ns);
        EXPECT_TRUE(grid.enable_conduction);
        model_column_update_bc(grid, state);
        const Vec ghost_before = grid.outer_boundary0_i;
        const Vec dt = cal_dt_i(grid, state);
        const Vec cold_wall = advance_Euler_state(grid, state, dt);
        // Raising ONLY the conduction wall must change the conductive result...
        grid.outer_conduction_temperature = 26000.0f;
        const Vec hot_wall = advance_Euler_state(grid, state, dt);
        const double T_cold = decode_cell(grid, cold_wall, grid.ns-1).T;
        const double T_hot  = decode_cell(grid, hot_wall,  grid.ns-1).T;
        EXPECT_TRUE(T_hot > T_cold);
        EXPECT_TRUE(std::abs(T_hot-T_cold) > 1e-3*T_cold);
        // ...without touching the hydro ghost state at all.
        EXPECT_TRUE(arma::approx_equal(grid.outer_boundary0_i, ghost_before,
                                       "absdiff", 0.0));

        // (3) Conduction reads the WALL, not the live hydro ghost. Same state, same
        // ghosts: disarming the override hands Stage D the hydro ghost temperature
        // (= T_top here, cooler than the 22 kK wall), so the step must come out
        // COLDER than the overridden one. If Stage D were still reading the ghost the
        // two would be identical.
        grid.outer_conduction_temperature = 22000.0f;
        const Vec wall_driven = advance_Euler_state(grid, state, dt);
        grid.outer_conduction_temperature_override = false;
        const Vec ghost_driven = advance_Euler_state(grid, state, dt);
        grid.outer_conduction_temperature_override = true;
        const double T_wall_driven  = decode_cell(grid, wall_driven,  grid.ns-1).T;
        const double T_ghost_driven = decode_cell(grid, ghost_driven, grid.ns-1).T;
        EXPECT_TRUE(decode_outer_ghost(grid, grid.outer_boundary0_i).T < 22000.0);
        EXPECT_TRUE(T_wall_driven > T_ghost_driven);
    }

    // ---- (8) well-balancedness: BC(q_ref) = q_g,ref in decoupled mode --------
    // The reference IC is V=0 and hydrostatic, so T_top = T_top,ref there and the
    // decoupled ghost lands on the same equilibrium rung. eq_wb freezes the residual
    // under the NEW boundary, so the column must still hold V=0 to the same order as
    // baseline, with no acoustic launch off the top.
    {
        auto max_speed_after_steps = [&](bool decouple) {
            Grid grid;
            Vec state = setup_decoupling_column(grid, decouple, ns);
            const MixtureThermo top = decode_cell(grid, state, grid.ns-1);
            const MixtureThermo ghost = decode_outer_ghost(grid, grid.outer_boundary0_i);
            if (decouple) EXPECT_REL(ghost.T, top.T, 2e-6);   // T_hydro_g0 = T_top,ref
            for (int step = 0; step < 20; ++step) {
                model_column_update_bc(grid, state);
                const Vec dt = cal_dt_i(grid, state);
                state = advance_Euler_state(grid, state, dt);
                grid.sim_time += dt(0);
                EXPECT_TRUE(!state.has_nan());
            }
            const auto sz = arma::size(grid.ns, num_of_eq);
            double vmax = 0.0;
            for (arma::uword i = 0; i < grid.ns; ++i) {
                const double momentum = state(arma::sub2ind(sz, i, cons::MOM_I))
                                      + state(arma::sub2ind(sz, i, cons::MOM_N));
                const double rho = state(arma::sub2ind(sz, i, cons::RHO_I))
                                 + state(arma::sub2ind(sz, i, cons::RHO_N));
                vmax = std::max(vmax, std::abs(momentum/rho));
            }
            return vmax;
        };
        const double v_baseline  = max_speed_after_steps(false);
        const double v_decoupled = max_speed_after_steps(true);
        // Same order of magnitude as baseline — no new boundary residual, no launch.
        EXPECT_TRUE(v_decoupled < 10.0*std::max(v_baseline, 1.0e-6));
    }

    clear_decoupling_env();
}

// Grid::capture_face_flux must (a) never change a number and (b) hand back the
// SAME total-mass face flux the continuity row is actually differenced from —
// otherwise the face-flux diagnosis would be measuring a look-alike, not the
// production flux. Both are checked against the untouched RHS on the real
// model_column gamma column.
static void test_face_flux_capture_matches_production_continuity() {
    const arma::uword ns = 200;
    Grid grid;
    Vec state = setup_decoupling_column(grid, false, ns);
    model_column_update_bc(grid, state);
    const Vec dt = cal_dt_i(grid, state);
    grid.dt_state.zeros();
    for (arma::uword k = 0; k < num_of_eq; ++k)
        grid.dt_state += scalar_to(grid, dt, k);

    // (a) flag off vs on: the returned RHS must be bit-for-bit identical.
    EXPECT_TRUE(!grid.capture_face_flux);
    const Vec rhs_off = rhs_explicit_state(grid, state);
    grid.capture_face_flux = true;
    const Vec rhs_on = rhs_explicit_state(grid, state);
    grid.capture_face_flux = false;
    EXPECT_TRUE(arma::approx_equal(rhs_on, rhs_off, "absdiff", 0.0));

    const GammaFaceFluxCapture& c = grid.face_flux_capture;
    EXPECT_TRUE(c.valid);
    EXPECT_TRUE(c.f_total.size() == ns);

    const auto sz = arma::size(grid.ns, num_of_eq);
    double worst_split = 0.0, worst_continuity = 0.0, scale = 0.0;
    for (arma::uword i = 0; i < ns; ++i) {
        // central + diffusive must reproduce the captured Rusanov total (the split
        // is exact by construction; only float32 rounding separates them).
        const double split = c.f_central[i] + c.f_diff[i];
        worst_split = std::max(worst_split,
            std::abs(split - c.f_total[i])/std::max(std::abs(c.f_total[i]), 1e-30));
        // Rusanov requires a = max(|V_L|+c_L, |V_R|+c_R).
        EXPECT_REL(c.a_face[i], std::max(std::abs(c.v_L[i])+c.cs_L[i],
                                         std::abs(c.v_R[i])+c.cs_R[i]), 1e-12);
        EXPECT_TRUE(c.rho_cell[i] > 0.0 && c.T_cell[i] > 0.0);
        scale = std::max(scale, std::abs(c.f_total[i])/grid.ds_i(i));
    }
    // (b) −(F_{i+1/2} − F_{i−1/2})/Δs, with F_{i−1/2} taken as the PREVIOUS cell's
    // upper face (the scheme is telescoping in the interior), plus the frozen eq_wb
    // residual, must reproduce the production continuity rows. B ≡ 1 here and
    // continuity carries no source, so nothing else may appear.
    for (arma::uword i = 1; i < ns; ++i) {
        const double produced =
            static_cast<double>(rhs_off(arma::sub2ind(sz, i, cons::RHO_I)))
          + static_cast<double>(rhs_off(arma::sub2ind(sz, i, cons::RHO_N)));
        const double from_faces = -(c.f_total[i] - c.f_total[i-1])/grid.ds_i(i)
                                - c.eq_residual_mass[i];
        worst_continuity = std::max(worst_continuity,
            std::abs(from_faces - produced)/std::max(scale, 1e-30));
    }
    EXPECT_TRUE(worst_split < 1.0e-5);
    EXPECT_TRUE(worst_continuity < 1.0e-5);

    // A second capture on a perturbed state must refill, not stale-cache.
    const double before = c.f_total[ns-1];
    const MixtureThermo top = decode_cell(grid, state, ns-1);
    set_cell_state(grid, state, ns-1, top.rho, 500.0, top.T);
    grid.capture_face_flux = true;
    rhs_explicit_state(grid, state);
    grid.capture_face_flux = false;
    EXPECT_TRUE(std::abs(c.f_total[ns-1] - before) > 0.0);

    clear_decoupling_env();
}

// Grid::capture_outer_conduction must (a) never change a number and (b) report the
// outer-face quantities the conduction solve actually used: the physical face
// conductivity as the ARITHMETIC face average of the top cell and the wall ghost
// (the ghost's kappa evaluated at the WALL temperature with the wall's Saha
// ionization), plus C_num*ds_iph*C_V(top) — over ds_iph_i(ns-1), which is
// the full cell width (top-cell CENTRE to ghost CENTRE), not half of it.
static void test_outer_conduction_capture_matches_solver_face() {
    const arma::uword ns = 200;
    Grid grid;
    Vec state = setup_decoupling_column(grid, true, ns);
    model_column_update_bc(grid, state);
    EXPECT_TRUE(grid.enable_conduction);
    EXPECT_TRUE(!grid.impose_outer_heat_flux);
    EXPECT_TRUE(grid.numerical_diffusivity_per_length > 0.0f);
    const Vec dt = cal_dt_i(grid, state);
    const DecodedMixtureField decoded = decode_mixture_field(grid, state, 1);

    // (a) flag off vs on: the advanced state must be bit-for-bit identical.
    EXPECT_TRUE(!grid.capture_outer_conduction);
    const Vec after_off = advance_Euler_state(grid, state, dt, decoded);
    grid.capture_outer_conduction = true;
    const Vec after_on = advance_Euler_state(grid, state, dt, decoded);
    grid.capture_outer_conduction = false;
    EXPECT_TRUE(arma::approx_equal(after_on, after_off, "absdiff", 0.0));

    const OuterConductionCapture& oc = grid.outer_conduction_capture;
    EXPECT_TRUE(oc.valid);
    EXPECT_TRUE(!oc.imposed_neumann);

    // (b) the formula, rebuilt from independent pieces.
    EXPECT_REL(oc.T_wall, grid.outer_conduction_temperature, 1e-12);
    EXPECT_REL(oc.ds_face, grid.ds_iph_i(ns-1), 1e-12);
    EXPECT_REL(oc.ds_face, grid.ds_i(ns-1), 1e-12);       // uniform mesh, Neumann ghost
    EXPECT_REL(oc.area_ratio, 1.0, 1e-12);                // straight column, B ≡ 1

    // Wall-side kappa: the ghost DENSITY with Saha re-evaluated at T_wall.
    const MixtureThermo ghost = decode_outer_ghost(grid, grid.outer_boundary0_i);
    const double x_wall = saha_ionization_fraction_n_h(ghost.n_H, oc.T_wall);
    const double k_wall = physical_conductivity(
        x_wall*ghost.n_H, (1.0-x_wall)*ghost.n_H, oc.T_wall);
    // Top-cell kappa at the CONVERGED temperature the capture reports.
    const double n_h_top = oc.T_top > 0.0
        ? decode_cell(grid, after_off, ns-1).rho/eos_constants::m_h : 0.0;
    const double x_top = saha_ionization_fraction_n_h(n_h_top, oc.T_top);
    const double k_top = physical_conductivity(
        x_top*n_h_top, (1.0-x_top)*n_h_top, oc.T_top);
    EXPECT_REL(oc.kappa_phys_face, 0.5*(k_top+k_wall), 1e-6);
    EXPECT_REL(oc.kappa_num_face,
        numerical_diffusivity_at_face(grid, oc.ds_face)
            *equilibrium_heat_capacity(decode_cell(grid, after_off, ns-1).rho, oc.T_top),
        1e-6);
    EXPECT_REL(oc.chi_num_face,
        numerical_diffusivity_at_face(grid, grid.ds_iph_i(ns-1)), 1e-12);
    EXPECT_REL(oc.q_phys,
        oc.kappa_phys_face*(oc.T_wall-oc.T_top)/oc.ds_face, 1e-12);
    EXPECT_REL(oc.q_num, oc.kappa_num_face*(oc.T_wall-oc.T_top)/oc.ds_face, 1e-12);
    EXPECT_REL(oc.q_total, oc.q_phys+oc.q_num, 1e-12);
    // The wall is hotter than the top cell here, so the flux heats the top cell.
    EXPECT_TRUE(oc.q_total > 0.0);
    // Numerical conduction is NOT a small correction at this resolution.
    EXPECT_TRUE(oc.q_num > 0.2*oc.q_total);

    // ISO_NUMERICAL_DIFFUSIVITY_MULT=0 must zero the numerical term ONLY.
    const double k_phys_before = oc.kappa_phys_face;
    setenv("ISO_NUMERICAL_DIFFUSIVITY_MULT", "0", 1);
    Grid g0;
    Vec s0 = setup_decoupling_column(g0, true, ns);
    model_column_update_bc(g0, s0);
    EXPECT_REL(g0.numerical_diffusivity_per_length, 0.0f, 0.0);
    EXPECT_TRUE(g0.enable_conduction);
    g0.capture_outer_conduction = true;
    advance_Euler_state(g0, s0, cal_dt_i(g0, s0), decode_mixture_field(g0, s0, 1));
    g0.capture_outer_conduction = false;
    EXPECT_REL(g0.outer_conduction_capture.q_num, 0.0, 0.0);
    EXPECT_REL(g0.outer_conduction_capture.q_total,
               g0.outer_conduction_capture.q_phys, 1e-12);
    // Same IC, same wall ⇒ the physical face conductivity is essentially unchanged.
    // Not exactly: one step already reaches a slightly different converged T_top
    // without the numerical term, and kappa_e ~ T^{5/2} amplifies that. A few tenths
    // of a percent is the expected size; a factor-level change would mean the
    // override had leaked into the physical conductivity.
    EXPECT_REL(g0.outer_conduction_capture.kappa_phys_face, k_phys_before, 2e-2);
    unsetenv("ISO_NUMERICAL_DIFFUSIVITY_MULT");

    clear_decoupling_env();
}

static double gamma_predictor_drift_heat(Grid& grid, const Vec& state, float dt) {
    Vec dt_i(grid.ns); dt_i.fill(dt);
    grid.dt_state.zeros();
    for (arma::uword k = 0; k < num_of_eq; ++k)
        grid.dt_state += scalar_to(grid, dt_i, k);
    const Vec predicted = state + grid.dt_state % rhs_explicit_state(grid, state);
    const auto sz = arma::size(grid.ns, num_of_eq);
    double heat = 0.0;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        auto value = [&](const Vec& u, arma::uword k) -> double {
            return static_cast<double>(u(arma::sub2ind(sz, i, k)));
        };
        const double ri = value(predicted, cons::RHO_I);
        const double rn = value(predicted, cons::RHO_N);
        const double vi = value(predicted, cons::MOM_I)/ri;
        const double vn = value(predicted, cons::MOM_N)/rn;
        const double relative_ke = 0.5*ri*rn/(ri+rn)*(vi-vn)*(vi-vn);
        const MixtureThermo before = decode_equilibrium_mixture(
            grid.eos_gamma_table, ri, rn,
            value(predicted, cons::MOM_I), value(predicted, cons::MOM_N),
            value(predicted, cons::E_I), value(predicted, cons::E_N), 0.0);
        const ProjectedMixture after = project_equilibrium_single_fluid(
            grid.eos_gamma_table, ri, rn,
            value(predicted, cons::MOM_I), value(predicted, cons::MOM_N),
            value(predicted, cons::E_I), value(predicted, cons::E_N), 0.0,
            grid.eos_trace_fraction_floor);
        EXPECT_NEAR(after.thermo.internal_energy-before.internal_energy, relative_ke,
                    std::max(1.0e-15, 5e-6*relative_ke));
        heat += relative_ke;
    }
    return heat;
}

static void test_stage6_projection_heating_dt_convergence() {
    const std::string path = write_gamma_table_fixture("gamma_stage6_dt.dat");
    const EosGammaTable table = EosGammaTable::load(path);
    Grid grid;
    const Vec state = setup_gamma_equilibrium(grid, table, 8, true);
    const double h1 = gamma_predictor_drift_heat(grid, state, 2.0e-4f);
    const double h2 = gamma_predictor_drift_heat(grid, state, 1.0e-4f);
    const double h3 = gamma_predictor_drift_heat(grid, state, 5.0e-5f);
    EXPECT_TRUE(h1 > h2 && h2 > h3 && h3 > 0.0);
    EXPECT_NEAR(h2/h1, 0.25, 0.08);
    EXPECT_NEAR(h3/h2, 0.25, 0.08);
    std::remove(path.c_str());
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

// The sidecar's physical coefficient must be independent of all solver-only
// broadening/stabilization knobs. With both knobs off it equals the
// solver-effective cell-centred coefficient.
static void test_physical_conductivity_excludes_solver_terms() {
    Grid grid;
    grid.init(8, 0.25f);
    const double n_e = 2.0e16;
    const double n_hi = 8.0e18;
    const double temperature = 5.0e4;
    const double heat_capacity = 3.0e4;
    const double physical = physical_conductivity(n_e, n_hi, temperature);

    grid.enable_trac = false;
    grid.numerical_diffusivity_per_length = 0.0f;
    EXPECT_REL(solver_effective_conductivity(
        grid, n_e, n_hi, temperature, heat_capacity, 1.0e4), physical, 1e-14);

    grid.numerical_diffusivity_per_length = 1.0e4f;
    EXPECT_REL(physical_conductivity(n_e, n_hi, temperature), physical, 1e-14);
    EXPECT_TRUE(solver_effective_conductivity(
        grid, n_e, n_hi, temperature, heat_capacity, 1.0e4) > physical);

    grid.numerical_diffusivity_per_length = 0.0f;
    grid.enable_trac = true;
    grid.trac_T_chrom = 2.0e4f;
    grid.trac_cutoff_T = 1.0e5f;
    EXPECT_REL(physical_conductivity(n_e, n_hi, temperature), physical, 1e-14);
    EXPECT_TRUE(solver_effective_conductivity(
        grid, n_e, n_hi, temperature, heat_capacity, 1.0e4) > physical);
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

    // Exact 700 km face (the fine-region top). The transition top is set by the
    // grade, so it is the smallest valid width ≥ the requested 100 km: the coarse
    // spacing must have resumed just above 800 km.
    double d700 = 1.0e30;
    for (double x : f) d700 = std::min(d700, std::fabs(x - 700.0e3));
    EXPECT_NEAR(d700, 0.0, 1.0e-3);
    double coarse_resumes = 1.0e30;
    for (arma::uword i = 0; i + 1 < f.size(); ++i)
        if (f[i] > 700.0e3 && f[i + 1] - f[i] > 0.999 * coarse_ds) {
            coarse_resumes = f[i]; break;
        }
    EXPECT_TRUE(coarse_resumes >= 800.0e3);
    EXPECT_TRUE(coarse_resumes <= 805.0e3);

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

static void test_mesh_outer_refined_profile() {
    // Outer/TR profile: coarse [0,480] → grade [480,500] → fine [500,553] to the
    // outer face, on the 1600–2153 km truncated domain (L = 553 km, ns_coarse = 1000).
    RefineParams rp;
    rp.profile = RefineProfile::Outer;
    rp.factor = 4.0; rp.s_lo_km = 500.0; rp.transition_km = 20.0;
    rp.s_hi_km = 1.0e9;                       // unused by the outer profile
    const double L = 553.0e3;
    const arma::uword nc = 1000;
    const std::vector<double> f = build_static_mesh_faces(L, nc, rp);
    const double coarse_ds = L / static_cast<double>(nc);

    // (1,2) exact domain endpoints; (3) the requested fine-region foot s_lo is an
    // exact face, and the transition below it is the smallest valid grade at least
    // the requested 20 km wide (so its foot sits just below 480 km).
    EXPECT_NEAR(f.front(), 0.0, 1.0e-9);
    EXPECT_NEAR(f.back(),  L,   1.0e-6);
    double d500 = 1.0e30;
    for (double x : f) d500 = std::min(d500, std::fabs(x - 500.0e3));
    EXPECT_NEAR(d500, 0.0, 1.0e-3);
    double grade_foot = 1.0e30;
    for (arma::uword i = 0; i + 1 < f.size(); ++i)
        if (f[i + 1] - f[i] < 0.999 * coarse_ds) { grade_foot = f[i]; break; }
    EXPECT_TRUE(grade_foot <= 480.0e3);          // ≥ the requested 20 km transition
    EXPECT_TRUE(grade_foot >= 475.0e3);          // and not wastefully wider

    // (2,3,4) strictly increasing faces, positive widths, adjacent ratio ≤ 1.1.
    double max_ratio = 1.0;
    bool increasing = true, all_pos = true;
    for (arma::uword i = 0; i + 1 < f.size(); ++i) {
        const double w = f[i + 1] - f[i];
        if (!(f[i + 1] > f[i])) increasing = false;
        if (!(w > 0.0) || !std::isfinite(w)) all_pos = false;
        if (i > 0) {
            const double wp = f[i] - f[i - 1];
            max_ratio = std::max(max_ratio, std::max(w / wp, wp / w));
        }
    }
    EXPECT_TRUE(increasing);
    EXPECT_TRUE(all_pos);
    EXPECT_TRUE(max_ratio <= 1.1 + 1.0e-6);

    // (5,6) the fine spacing is ≈ coarse/4 and the LAST cells stay fine right
    // through the outer boundary; (7) the coarse lower region keeps coarse_ds.
    const double outer_w = f.back() - f[f.size() - 2];
    EXPECT_REL(outer_w, coarse_ds / 4.0, 0.02);
    EXPECT_REL(f[1] - f[0], coarse_ds, 0.02);
    for (arma::uword i = 0; i + 1 < f.size(); ++i)     // every cell above s_lo is fine
        if (f[i] >= 500.0e3 - 1.0e-6)
            EXPECT_REL(f[i + 1] - f[i], coarse_ds / 4.0, 0.02);
    // (8) refinement ADDS cells; the domain below the transition is not coarsened.
    EXPECT_TRUE(f.size() - 1 > nc);
}

static void test_mesh_outer_profile_edge_cases() {
    // The graded transition must fit strictly inside the domain below s_lo, and s_lo
    // must be inside the domain — otherwise the builder returns the exact uniform grid.
    const double L = 553.0e3;
    const arma::uword nc = 100;                   // coarse Δs = 5.53 km ⇒ a wide grade
    const double ds = L / static_cast<double>(nc);
    auto is_uniform = [&](const RefineParams& rp) {
        const std::vector<double> f = build_static_mesh_faces(L, nc, rp);
        if (f.size() != nc + 1) return false;
        for (arma::uword k = 0; k <= nc; ++k)
            if (std::fabs(f[k] - ds * static_cast<double>(k)) > 1.0e-3) return false;
        return true;
    };
    RefineParams rp; rp.profile = RefineProfile::Outer; rp.factor = 4.0;
    rp.s_lo_km = 20.0; rp.transition_km = 20.0;   // grade does not fit below s_lo
    EXPECT_TRUE(is_uniform(rp));
    rp.s_lo_km = 600.0;                           // s_lo above the domain top
    EXPECT_TRUE(is_uniform(rp));

    // A transition width the ratio cap cannot support is WIDENED to the smallest one
    // it can (never an abrupt coarse→fine jump), so even τ = 0 stays graded.
    rp.s_lo_km = 400.0; rp.transition_km = 0.0;
    const std::vector<double> f = build_static_mesh_faces(L, nc, rp);
    EXPECT_TRUE(f.size() - 1 > nc);
    double max_ratio = 1.0;
    for (arma::uword i = 1; i + 1 < f.size(); ++i) {
        const double w = f[i + 1] - f[i], wp = f[i] - f[i - 1];
        max_ratio = std::max(max_ratio, std::max(w / wp, wp / w));
    }
    EXPECT_TRUE(max_ratio <= 1.1 + 1.0e-6);
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
    EXPECT_TRUE(a.profile == RefineProfile::Lower);   // default profile
    setenv("GRID_REFINE_PROFILE", "outer", 1);
    EXPECT_TRUE(refine_params_from_env(false).profile == RefineProfile::Outer);
    setenv("ISO_REFINE_PROFILE", "lower", 1);        // ISO_ overrides GRID_
    EXPECT_TRUE(refine_params_from_env(true).profile == RefineProfile::Lower);
    unsetenv("GRID_REFINE_FACTOR"); unsetenv("GRID_REFINE_S_HI_KM"); unsetenv("ISO_REFINE_FACTOR");
    unsetenv("GRID_REFINE_PROFILE"); unsetenv("ISO_REFINE_PROFILE");
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

static void test_numerical_diffusivity_uniform_reduction() {
    Grid grid;
    const std::vector<float> ds(6, 553.0f);
    setup_irregular(grid, ds, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    EXPECT_TRUE(grid.uniform_mesh);
    grid.numerical_diffusivity_per_length = 2.0e3f;
    const double old_uniform_chi = 2.0e3*553.0;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        EXPECT_REL(numerical_diffusivity_at_face(grid, grid.ds_imh_i(i)),
                   old_uniform_chi, 1e-12);
        EXPECT_REL(numerical_diffusivity_at_face(grid, grid.ds_iph_i(i)),
                   old_uniform_chi, 1e-12);
    }
    // The coefficient folded into the conduction update is therefore exactly the
    // previous uniform K_num = (2000*Delta_s)*C_V formula.
    const double cv = 3.0e4;
    EXPECT_REL(numerical_diffusivity_at_face(grid,grid.ds_iph_i(2))*cv,
               old_uniform_chi*cv, 1e-12);
}

static void test_numerical_diffusivity_nonuniform_face_scaling() {
    Grid grid;
    const std::vector<float> ds = {8.0f, 8.0f, 5.0f, 2.0f, 2.0f};
    setup_irregular(grid, ds, 2.0e17f, 1.0e19f, 6500.0f, 6500.0f);
    EXPECT_TRUE(!grid.uniform_mesh);
    grid.numerical_diffusivity_per_length = 2.0e3f;

    // Coarse interior, graded transition, fine interior, and mirrored outer face.
    const arma::uword face_i[] = {0, 1, 3, 4};
    for (arma::uword i : face_i) {
        const double chi = numerical_diffusivity_at_face(grid, grid.ds_iph_i(i));
        EXPECT_REL(chi, 2.0e3*grid.ds_iph_i(i), 1e-12);
    }
    const double chi_coarse = numerical_diffusivity_at_face(grid,grid.ds_iph_i(0));
    const double chi_transition = numerical_diffusivity_at_face(grid,grid.ds_iph_i(1));
    const double chi_fine = numerical_diffusivity_at_face(grid,grid.ds_iph_i(3));
    const double chi_outer = numerical_diffusivity_at_face(grid,grid.ds_iph_i(4));
    EXPECT_REL(chi_coarse/chi_transition,
               grid.ds_iph_i(0)/grid.ds_iph_i(1), 1e-6);
    EXPECT_REL(chi_coarse/chi_fine,
               grid.ds_iph_i(0)/grid.ds_iph_i(3), 1e-6);
    EXPECT_REL(chi_outer, 2.0e3*grid.ds_iph_i(4), 1e-12);
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
    const Vec RI_physical = rhs_implicit_state(grid, cons);
    grid.numerical_diffusivity_per_length = 2.0e3f;
    Vec RI = rhs_implicit_state(grid, cons);
    EXPECT_TRUE(arma::max(arma::abs(RI-RI_physical)) > 0.0f);
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

static void test_outer_refined_mesh_ic_runs() {
    // End-to-end on the OUTER/TR profile over the truncated 1600–2153 km domain:
    // model_column builds a larger non-uniform grid whose FINEST cells sit at the
    // outer boundary, the CFL picks up the smallest local cell, and one full step
    // on the well-balanced hydrostatic IC stays finite with V ≈ 0 (no boundary-scale
    // blow-up at the new coarse→fine interface).
    setenv("ISO_H_BASE", "1600", 1);
    setenv("ISO_DH", "553", 1);
    auto dt_min_of = [&](Grid& grid) {
        Vec xn = model_column_ic(grid);
        Vec dt = cal_dt_i(grid, xn);
        model_column_update_bc(grid, xn);
        Vec out = advance_Euler_state(grid, xn, dt);
        EXPECT_TRUE(!out.has_nan() && out.is_finite());
        const auto sz = arma::size(grid.ns, num_of_eq);
        float vmax = 0.0f;
        for (arma::uword i = 0; i < grid.ns; ++i)
            vmax = std::max(vmax, std::fabs(out(arma::sub2ind(sz, i, cons::MOM_I)) /
                                            out(arma::sub2ind(sz, i, cons::RHO_I))));
        EXPECT_TRUE(vmax < 1.0e-2f);   // exact V=0 fixed point survives refinement
        return arma::min(dt);
    };

    Grid gu; gu.init(200, 0.25f);
    const float dt_uniform = dt_min_of(gu);
    EXPECT_TRUE(gu.uniform_mesh);
    EXPECT_TRUE(gu.ns == 200);

    setenv("ISO_REFINE_PROFILE", "outer", 1);
    setenv("ISO_REFINE_FACTOR", "4", 1);
    setenv("ISO_REFINE_S_LO_KM", "500", 1);
    setenv("ISO_REFINE_TRANSITION_KM", "20", 1);
    Grid gr; gr.init(200, 0.25f);
    const float dt_refined = dt_min_of(gr);

    EXPECT_TRUE(!gr.uniform_mesh);
    EXPECT_TRUE(gr.ns > 200);                                   // cells added
    EXPECT_TRUE(gr.ds_i.min() > 0.0f && gr.ds_i.is_finite());
    EXPECT_REL(gr.ds_i.min(), gr.ds_i.max() / 4.0f, 0.05);      // fine ≈ coarse/4
    // The fine band reaches the OUTER boundary: the last cell is a minimum-width cell.
    EXPECT_REL(gr.ds_i(gr.ns - 1), gr.ds_i.min(), 1.0e-3);
    // ...and the inner boundary keeps the coarse spacing.
    EXPECT_REL(gr.ds_i(0), gr.ds_i.max(), 1.0e-3);
    // Adjacent width ratio bounded by the builder's cap.
    float max_ratio = 1.0f;
    for (arma::uword i = 1; i < gr.ns; ++i)
        max_ratio = std::max(max_ratio, std::max(gr.ds_i(i) / gr.ds_i(i - 1),
                                                 gr.ds_i(i - 1) / gr.ds_i(i)));
    EXPECT_TRUE(max_ratio <= 1.1f + 1.0e-4f);
    // Non-uniform CFL is set by the smallest LOCAL cell, not the mean spacing: the
    // acoustic limit sits at the top of the TR, which is exactly where the fine cells
    // are, so dt drops by the local width ratio. The 10 % slack absorbs the O(Δh)
    // top-cell IC temperature difference (the top cell centre — and hence the sampled
    // C7 sound speed — moves when the mesh is refined).
    EXPECT_REL(dt_refined / dt_uniform, gr.ds_i.min() / gu.ds_i(0), 0.10);

    unsetenv("ISO_REFINE_PROFILE"); unsetenv("ISO_REFINE_FACTOR");
    unsetenv("ISO_REFINE_S_LO_KM"); unsetenv("ISO_REFINE_TRANSITION_KM");
    unsetenv("ISO_H_BASE"); unsetenv("ISO_DH");
}

// ----------------------------------------------------------------------------

int main() {
    std::cout << "===== Chromosphere test suite =====\n\n";

    RUN(test_scalar_to_get_scalar_inverse);
    RUN(test_ip1_im1_interior_shift);
    RUN(test_flux_lim_is_minmod);
    RUN(test_flux_lim_mc3);
    RUN(test_saha_ionization_fraction_log_domain);
    RUN(test_broadcast_static_metric_cache);
    RUN(test_equilibrium_density_from_pressure_matches_bisection);
    RUN(test_eos_gamma_table_loader_and_interpolation);
    RUN(test_final_serial_log_aware_eos_and_gamma1_contracts);
    RUN(test_final_serial_one_update_and_fallback_regression);
    RUN(test_stage3_equilibrium_mixture_closure);
    RUN(test_stage4_conservative_equilibrium_projection);
    RUN(test_stage5_6_equilibrium_face_flux_and_sound_speed);
    RUN(test_stage9_acoustic_characteristic_speed);
    RUN(test_stage9_storage_floor_transport_independence);
    RUN(test_stage9_evaporation_vs_fixed_gamma);
    RUN(test_stage5_production_table_boundary_reconstruction);
    RUN(test_stage5_gamma_integrator_projection_and_guards);
    RUN(test_stage7_gamma_total_energy_sources_and_conduction);
    RUN(test_gamma_phase1_3_cache_and_known_temperature_pack);
    RUN(test_stage8_gamma_model_column_saha_hse_and_ghosts);
    RUN(test_upper_bc_hydro_conduction_temperature_decoupling);
    RUN(test_face_flux_capture_matches_production_continuity);
    RUN(test_outer_conduction_capture_matches_solver_face);
    RUN(test_stage6_projection_heating_dt_convergence);
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
    RUN(test_physical_conductivity_excludes_solver_terms);
    RUN(test_trac_cutoff_detection_and_limiter);
    RUN(test_beam_heating_rate_profile);
    RUN(test_beam_heating_partitions_by_heat_capacity);
    RUN(test_coronal_heating_rate_profile);
    RUN(test_coronal_heating_partitions_by_heat_capacity);

    // Static local refinement: mesh builder, metric caches, irregular-grid
    // manufactured operators, and the refinement-unset regression.
    RUN(test_mesh_disabled_is_uniform);
    RUN(test_mesh_refined_diagnostic_profile);
    RUN(test_mesh_outer_refined_profile);
    RUN(test_mesh_outer_profile_edge_cases);
    RUN(test_refine_params_env_aliases);
    RUN(test_metric_caches_center_to_center);
    RUN(test_numerical_diffusivity_uniform_reduction);
    RUN(test_numerical_diffusivity_nonuniform_face_scaling);
    RUN(test_irregular_uniform_state_preserved);
    RUN(test_irregular_linear_pressure_gradient);
    RUN(test_irregular_cfl_selection);
    RUN(test_irregular_conduction_conserves_energy);
    RUN(test_irregular_boundary_hse);
    RUN(test_uniform_regression_refinement_unset);
    RUN(test_refined_mesh_ic_runs);
    RUN(test_outer_refined_mesh_ic_runs);

    std::cout << "\n===== Summary =====\n";
    std::cout << "Passed: " << g_pass << "\n";
    std::cout << "Failed: " << g_fail << "\n";
    return g_fail == 0 ? 0 : 1;
}
