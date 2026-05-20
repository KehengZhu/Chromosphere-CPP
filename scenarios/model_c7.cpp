#include "model_c7.hpp"

namespace chromosphere {

// Natural cubic spline interpolation, port of interp1 in ModChromosphereTest.f90
static Vec interp1_spline(const Vec& xData, const Vec& yData, const Vec& xVal) {
    const arma::uword n = xData.n_elem;
    Vec h(n - 1);
    Vec alpha(n, arma::fill::zeros);
    Vec a(n), b(n - 1), c(n, arma::fill::zeros), d(n - 1);
    Vec l(n, arma::fill::zeros), mu(n, arma::fill::zeros), z(n, arma::fill::zeros);

    for (arma::uword i = 0; i < n - 1; ++i) h(i) = xData(i + 1) - xData(i);
    for (arma::uword i = 1; i < n - 1; ++i)
        alpha(i) = 3.0f / h(i) * (yData(i + 1) - yData(i))
                 - 3.0f / h(i - 1) * (yData(i) - yData(i - 1));

    l(0) = 1.0f;
    for (arma::uword i = 1; i < n - 1; ++i) {
        l(i)  = 2.0f * (xData(i + 1) - xData(i - 1)) - h(i - 1) * mu(i - 1);
        mu(i) = h(i) / l(i);
        z(i)  = (alpha(i) - h(i - 1) * z(i - 1)) / l(i);
    }
    l(n - 1) = 1.0f;

    for (long j = (long)n - 2; j >= 0; --j) {
        c(j) = z(j) - mu(j) * c(j + 1);
        b(j) = (yData(j + 1) - yData(j)) / h(j) - h(j) * (c(j + 1) + 2.0f * c(j)) / 3.0f;
        d(j) = (c(j + 1) - c(j)) / (3.0f * h(j));
        a(j) = yData(j);
    }

    Vec yVal(xVal.n_elem);
    for (arma::uword i = 0; i < xVal.n_elem; ++i) {
        arma::uword j = 0;
        while (j < n - 1 && xVal(i) > xData(j + 1)) ++j;
        if (j >= n - 1) j = n - 2;
        const float dx = xVal(i) - xData(j);
        yVal(i) = a(j) + b(j) * dx + c(j) * dx * dx + d(j) * dx * dx * dx;
    }
    return yVal;
}

// Model C7 chromosphere profile (writeup Table 1):
// height [km], n_e [m^-3], n_n [m^-3], T [K]
static const float MODEL_C7[15][4] = {
    {1.003e+03f, 1.903e+17f, 2.693e+19f, 6.225e+03f},
    {1.032e+03f, 2.021e+17f, 2.179e+19f, 6.315e+03f},
    {1.065e+03f, 2.091e+17f, 1.722e+19f, 6.400e+03f},
    {1.101e+03f, 2.104e+17f, 1.340e+19f, 6.474e+03f},
    {1.143e+03f, 1.995e+17f, 1.009e+19f, 6.531e+03f},
    {1.214e+03f, 1.760e+17f, 6.322e+18f, 6.576e+03f},
    {1.299e+03f, 1.489e+17f, 3.641e+18f, 6.598e+03f},
    {1.398e+03f, 1.438e+17f, 1.880e+18f, 6.610e+03f},
    {1.520e+03f, 1.423e+17f, 7.956e+17f, 6.623e+03f},
    {1.617e+03f, 1.267e+17f, 4.012e+17f, 6.633e+03f},
    {1.722e+03f, 1.027e+17f, 1.916e+17f, 6.643e+03f},
    {1.820e+03f, 7.969e+16f, 9.843e+16f, 6.652e+03f},
    {1.894e+03f, 6.450e+16f, 6.005e+16f, 6.660e+03f},
    {1.946e+03f, 5.526e+16f, 4.371e+16f, 6.667e+03f},
    {1.989e+03f, 4.826e+16f, 3.447e+16f, 6.674e+03f}
};

Vec model_c7_ic(Grid& grid) {
    // Phase-1 geometry choices (quiet-sun relaxation):
    //   * Uniform B = 1, ∂(1/B)/∂s = 0 — single straight field line; the
    //     flux-tube expansion term in the source/flux drops out trivially.
    //     A non-uniform B(s) (e.g. Spruit 1981 thin-flux-tube) is the
    //     follow-on; the integrator already supports it.
    //   * Gravity ON: φ_g(s) = g · (h(s) - h_base) populated below.
    //
    // HSE strategy (Option A from the Phase-2 plan): we leave C7's tabulated
    // (n, T) as the IC even though C7 was built with NLTE radiative balance,
    // not ideal-gas HSE. Under the code's EOS the C7 profile is off HSE by a
    // nearly uniform 20–26%: |dp/dh + ρg| / |ρg| ≈ 0.25 across the chromosphere
    // (dp/dh is more negative than -ρg, i.e. C7 is "lighter at the top" than
    // a pure-HSE atmosphere with the same T(h) would be). The run absorbs
    // this imbalance during the first few acoustic crossings; the residual
    // drift away from C7 is dominated by missing radiative cooling, not by
    // the IC choice.
    grid.B_imh.ones();
    grid.B_iph.ones();
    grid.B_i.ones();
    grid.dinvB_ds_i.zeros();

    // Tabulated profile (15 rows from writeup Table 1).
    Vec h_tab(15), ne_tab(15), nn_tab(15), T_tab(15);
    for (arma::uword i = 0; i < 15; ++i) {
        h_tab(i)  = MODEL_C7[i][0];
        ne_tab(i) = MODEL_C7[i][1];
        nn_tab(i) = MODEL_C7[i][2];
        T_tab(i)  = MODEL_C7[i][3];
    }

    // Faces: indices -1..ns+3 (size ns+5), uniformly over [h_tab(0), h_tab(14)].
    const arma::uword nF = grid.ns + 5;
    Vec h_F(nF);
    const float h0 = h_tab(0), h1 = h_tab(14);
    for (arma::uword k = 0; k < nF; ++k) {
        const float frac = (float)k / (float)(nF - 1);
        h_F(k) = h0 + (h1 - h0) * frac;
    }

    // φ_g at each face, measured relative to the C7 base height (h0 = 1003 km)
    // so the absolute potential stays comparable to thermal energy in float32.
    // Only ∂φ_g/∂s matters physically; the constant offset is gauge.
    Vec phi_g_F(nF);
    for (arma::uword k = 0; k < nF; ++k) {
        phi_g_F(k) = grid.g * (h_F(k) - h0) * 1000.0f; // km -> m
    }
    for (arma::uword i = 0; i < grid.ns; ++i) {
        grid.phi_g_imh(i) = phi_g_F(i);
        grid.phi_g_iph(i) = phi_g_F(i + 1);
    }

    Vec ne_F = interp1_spline(h_tab, ne_tab, h_F);
    Vec nn_F = interp1_spline(h_tab, nn_tab, h_F);
    Vec T_F  = interp1_spline(h_tab, T_tab,  h_F);

    // Extended-cell quantities: averages of adjacent faces, indices 0..ns+2.
    const arma::uword nE = grid.ns + 3;
    Vec len_E(nE), ne_E(nE), nn_E(nE), T_E(nE);
    for (arma::uword k = 0; k < nE; ++k) {
        len_E(k) = (h_F(k + 1) - h_F(k)) * 1000.0f; // km -> m
        ne_E(k)  = 0.5f * (ne_F(k + 1) + ne_F(k));
        nn_E(k)  = 0.5f * (nn_F(k + 1) + nn_F(k));
        T_E(k)   = 0.5f * (T_F(k + 1) + T_F(k));
    }

    // Cell-centered interior: first ns of the extended array.
    for (arma::uword i = 0; i < grid.ns; ++i) grid.ds_i(i) = len_E(i);

    // Build the conserved-variable state vector. Total energy includes the
    // ρ·φ_g term so that cons2prim (state.cpp:108) recovers the right pressure.
    Vec xn = arma::zeros<Vec>(grid.n_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float n_e   = ne_E(i);
        const float n_n   = nn_E(i);
        const float T     = T_E(i);
        const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I)) = n_e * grid.m_i;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_N)) = n_n * grid.m_n;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::MOM_I)) = 0.0f;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::MOM_N)) = 0.0f;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_I))   = 1.5f * grid.k_b * n_e * 2.0f * T + n_e * grid.m_i * phi_g;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_N))   = 1.5f * grid.k_b * n_n * T       + n_n * grid.m_n * phi_g;
    }

    // Outer ghost cells use extended-array entries ns+1 and ns+2.
    // Outer T_i is doubled at the outermost boundary to mimic the corona
    // (writeup §4). The ghost's φ_g is taken as the upper face of the last
    // interior cell, matching rhs.cpp:154 and integrators.cpp:188 which use
    // grid.phi_g_iph[ns-1] when reaching past the outer boundary.
    {
        const float phi_g_outer = grid.phi_g_iph(grid.ns - 1);
        for (arma::uword k = 0; k < 2; ++k) {
            const arma::uword idx = grid.ns + 1 + k;
            const float n_e = ne_E(idx);
            const float n_n = nn_E(idx);
            const float T   = T_E(idx);
            Vec& ob = (k == 0 ? grid.outer_boundary0_i : grid.outer_boundary1_i);
            ob(cons::RHO_I) = n_e * grid.m_i;
            ob(cons::RHO_N) = n_n * grid.m_n;
            ob(cons::MOM_I) = 0.0f;
            ob(cons::MOM_N) = 0.0f;
            ob(cons::E_I)   = 1.5f * grid.k_b * n_e * (4.0f * T) + n_e * grid.m_i * phi_g_outer; // T_e doubled at outer
            ob(cons::E_N)   = 1.5f * grid.k_b * n_n * (2.0f * T) + n_n * grid.m_n * phi_g_outer;
        }
    }

    // Inner ghost cells replicate the innermost extended cell. φ_g taken as
    // the lower face of cell 0 (which is 0 by gauge choice), matching the
    // im1 path in rhs.cpp:178 and integrators.cpp:189.
    {
        const float n_e = ne_E(0);
        const float n_n = nn_E(0);
        const float T   = T_E(0);
        const float phi_g_inner = grid.phi_g_imh(0);
        grid.inner_boundary0_i(cons::RHO_I) = n_e * grid.m_i;
        grid.inner_boundary0_i(cons::RHO_N) = n_n * grid.m_n;
        grid.inner_boundary0_i(cons::MOM_I) = 0.0f;
        grid.inner_boundary0_i(cons::MOM_N) = 0.0f;
        grid.inner_boundary0_i(cons::E_I)   = 1.5f * grid.k_b * n_e * 2.0f * T + n_e * grid.m_i * phi_g_inner;
        grid.inner_boundary0_i(cons::E_N)   = 1.5f * grid.k_b * n_n * T       + n_n * grid.m_n * phi_g_inner;
        grid.inner_boundary1_i = grid.inner_boundary0_i;
    }

    grid.broadcast();
    return xn;
}

void model_c7_update_bc(Grid& grid, const Vec& xn) {
    Vec rho_i  = get_scalar(grid, xn, cons::RHO_I);
    Vec rho_n  = get_scalar(grid, xn, cons::RHO_N);
    Vec rhoV_i = get_scalar(grid, xn, cons::MOM_I);
    Vec rhoU_n = get_scalar(grid, xn, cons::MOM_N);
    const float V_ns = rhoV_i(grid.ns - 1) / rho_i(grid.ns - 1);
    const float U_ns = rhoU_n(grid.ns - 1) / rho_n(grid.ns - 1);

    // Pull the pinned (n, T) state from the current outer-ghost storage, but
    // rebuild the energy with the halved velocity. Densities & temperatures
    // stay frozen at Model C7 values.
    auto refresh = [&](Vec& ob, float v_target, float u_target) {
        const float n_i_B = ob(cons::RHO_I) / grid.m_i;
        const float n_n_B = ob(cons::RHO_N) / grid.m_n;
        const float V_B_old = ob(cons::MOM_I) / (grid.m_i * n_i_B);
        const float U_B_old = ob(cons::MOM_N) / (grid.m_n * n_n_B);
        // Match the convention used in model_c7_ic for outer-ghost energies:
        // the ghost cell's φ_g is the upper face of the last interior cell.
        const float phi_g = grid.phi_g_iph(grid.ns - 1);
        const float T_i_B = (2.0f/3.0f * ob(cons::E_I) - 1.0f/3.0f * grid.m_i * n_i_B * V_B_old * V_B_old
                            - 2.0f/3.0f * grid.m_i * n_i_B * phi_g) / (2.0f * n_i_B * grid.k_b);
        const float T_n_B = (2.0f/3.0f * ob(cons::E_N) - 1.0f/3.0f * grid.m_n * n_n_B * U_B_old * U_B_old
                            - 2.0f/3.0f * grid.m_n * n_n_B * phi_g) / (n_n_B * grid.k_b);
        ob(cons::MOM_I) = grid.m_i * n_i_B * v_target;
        ob(cons::MOM_N) = grid.m_n * n_n_B * u_target;
        ob(cons::E_I)   = 1.5f * grid.k_b * n_i_B * 2.0f * T_i_B
                        + 0.5f * grid.m_i * n_i_B * v_target * v_target + grid.m_i * n_i_B * phi_g;
        ob(cons::E_N)   = 1.5f * grid.k_b * n_n_B * T_n_B
                        + 0.5f * grid.m_n * n_n_B * u_target * u_target + grid.m_n * n_n_B * phi_g;
    };

    const float v0 = 0.5f * V_ns;
    const float u0 = 0.5f * U_ns;
    refresh(grid.outer_boundary0_i, v0,        u0);
    refresh(grid.outer_boundary1_i, 0.5f * v0, 0.5f * u0);
}

} // namespace chromosphere
