#include "chromosphere.hpp"
#include <armadillo>

namespace chromosphere {

void chromo_init(uword ns_in, float CFL_in) {
    ns          = ns_in;
    num_of_elem = ns * num_of_eq;
    CFL         = CFL_in;

    gammamono = 5.0f / 3.0f;
    alpha_p   = 0.18f;
    pi        = 3.14159265358979323846f;

    m_i  = 1.6726219e-27f;
    m_n  = m_i;
    m_e  = 9.10938356e-31f;
    g    = 0.27395e3f;
    mu_0 = 4.0f * pi * 1.0e-7f;
    k_b  = 1.380649e-23f;
    e_   = 1.602176634e-19f;

    ds_i.zeros(ns);
    cellX_i.zeros(ns);  cellY_i.zeros(ns);  cellZ_i.zeros(ns);
    B_imh.zeros(ns);    B_iph.zeros(ns);    B_i.zeros(ns);
    dinvB_ds_i.zeros(ns);
    gPotential_imh.zeros(ns);
    gPotential_iph.zeros(ns);

    outer_boundary0_i.zeros(num_of_eq);
    outer_boundary1_i.zeros(num_of_eq);
    inner_boundary0_i.zeros(num_of_eq);
    inner_boundary1_i.zeros(num_of_eq);

    B_iimh.zeros(num_of_elem);
    B_iiph.zeros(num_of_elem);
    B_ii.zeros(num_of_elem);
    ds_ii.zeros(num_of_elem);
    dinvB_ds_ii.zeros(num_of_elem);
    dt_ii.zeros(num_of_elem);
}

// Natural cubic spline interpolation, port of interp1 in ModChromosphereTest.f90
static Vec interp1_spline(const Vec& xData, const Vec& yData, const Vec& xVal) {
    const uword n = xData.n_elem;
    Vec h(n - 1);
    Vec alpha(n, fill::zeros);
    Vec a(n), b(n - 1), c(n, fill::zeros), d(n - 1);
    Vec l(n, fill::zeros), mu(n, fill::zeros), z(n, fill::zeros);

    for (uword i = 0; i < n - 1; ++i) h(i) = xData(i + 1) - xData(i);
    for (uword i = 1; i < n - 1; ++i)
        alpha(i) = 3.0f / h(i) * (yData(i + 1) - yData(i))
                 - 3.0f / h(i - 1) * (yData(i) - yData(i - 1));

    l(0) = 1.0f;
    for (uword i = 1; i < n - 1; ++i) {
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
    for (uword i = 0; i < xVal.n_elem; ++i) {
        uword j = 0;
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

Vec model_c7_ic() {
    // Uniform B, no gravity in this run (matches Fortran defaults).
    B_imh.ones();        B_iph.ones();        B_i.ones();
    dinvB_ds_i.zeros();
    gPotential_imh.zeros();
    gPotential_iph.zeros();

    // Tabulated profile (15 rows from writeup Table 1)
    Vec h_tab(15), ne_tab(15), nn_tab(15), T_tab(15);
    for (uword i = 0; i < 15; ++i) {
        h_tab(i)  = MODEL_C7[i][0];
        ne_tab(i) = MODEL_C7[i][1];
        nn_tab(i) = MODEL_C7[i][2];
        T_tab(i)  = MODEL_C7[i][3];
    }

    // Faces: indices -1..ns+3 (size ns+5), uniformly distributed over [h_tab(0), h_tab(14)].
    const uword nF = ns + 5;
    Vec h_F(nF);
    const float h0 = h_tab(0), h1 = h_tab(14);
    for (uword k = 0; k < nF; ++k) {
        // Fortran index i runs from -1 to ns+3 with step (i+1)/(nF-1)
        const float frac = (float)k / (float)(nF - 1);
        h_F(k) = h0 + (h1 - h0) * frac;
    }

    Vec ne_F = interp1_spline(h_tab, ne_tab, h_F);
    Vec nn_F = interp1_spline(h_tab, nn_tab, h_F);
    Vec T_F  = interp1_spline(h_tab, T_tab,  h_F);

    // Extended-cell quantities: averages of adjacent faces, indices 0..ns+2 (size ns+3).
    const uword nE = ns + 3;
    Vec len_E(nE), ne_E(nE), nn_E(nE), T_E(nE);
    for (uword k = 0; k < nE; ++k) {
        len_E(k) = (h_F(k + 1) - h_F(k)) * 1000.0f; // km -> m
        ne_E(k)  = 0.5f * (ne_F(k + 1) + ne_F(k));
        nn_E(k)  = 0.5f * (nn_F(k + 1) + nn_F(k));
        T_E(k)   = 0.5f * (T_F(k + 1) + T_F(k));
    }

    // Cell-centered (interior) — first ns of extended array
    for (uword i = 0; i < ns; ++i) ds_i(i) = len_E(i);

    // Build conserved-variable state vector
    Vec xn = zeros<Vec>(num_of_elem);
    for (uword i = 0; i < ns; ++i) {
        const float ne = ne_E(i);
        const float nn = nn_E(i);
        const float T  = T_E(i);
        const float phig = 0.0f;
        xn(sub2ind(size(ns, num_of_eq), i, CNI)) = ne * m_i;
        xn(sub2ind(size(ns, num_of_eq), i, CNN)) = nn * m_n;
        xn(sub2ind(size(ns, num_of_eq), i, CNV)) = 0.0f;
        xn(sub2ind(size(ns, num_of_eq), i, CNU)) = 0.0f;
        xn(sub2ind(size(ns, num_of_eq), i, CEI)) = 1.5f * k_b * ne * 2.0f * T + ne * m_i * phig;
        xn(sub2ind(size(ns, num_of_eq), i, CEN)) = 1.5f * k_b * nn * T       + nn * m_n * phig;
    }

    // Outer ghost cells use extended-array entries ns+1 and ns+2 (Fortran indices).
    // Outer T_i is doubled at the outermost boundary to mimic the corona (writeup §4).
    for (uword k = 0; k < 2; ++k) {
        const uword idx = ns + 1 + k;
        const float ne = ne_E(idx);
        const float nn = nn_E(idx);
        const float T  = T_E(idx);
        Vec& ob = (k == 0 ? outer_boundary0_i : outer_boundary1_i);
        ob(CNI) = ne * m_i;
        ob(CNN) = nn * m_n;
        ob(CNV) = 0.0f;
        ob(CNU) = 0.0f;
        ob(CEI) = 1.5f * k_b * ne * (4.0f * T); // T_e doubled at outer
        ob(CEN) = 1.5f * k_b * nn * (2.0f * T);
    }

    // Inner ghost cells use extended-array indices 0 and -1 (negative in Fortran).
    // Here we approximate the "-1" face simply by replicating index 0 (innermost extended cell).
    {
        const float ne = ne_E(0);
        const float nn = nn_E(0);
        const float T  = T_E(0);
        inner_boundary0_i(CNI) = ne * m_i;
        inner_boundary0_i(CNN) = nn * m_n;
        inner_boundary0_i(CNV) = 0.0f;
        inner_boundary0_i(CNU) = 0.0f;
        inner_boundary0_i(CEI) = 1.5f * k_b * ne * 2.0f * T;
        inner_boundary0_i(CEN) = 1.5f * k_b * nn * T;
        inner_boundary1_i = inner_boundary0_i;
    }

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

    return xn;
}

void model_c7_update_bc(const Vec& xn) {
    // Recompute Ti, Tn, vv, uu at the outer-boundary ghost cells, halve velocity each step.
    Vec rhoi = get_scalar(xn, CNI);
    Vec rhon = get_scalar(xn, CNN);
    Vec rhov = get_scalar(xn, CNV);
    Vec rhou = get_scalar(xn, CNU);
    const float vv_ns = rhov(ns - 1) / rhoi(ns - 1);
    const float uu_ns = rhou(ns - 1) / rhon(ns - 1);

    // Pull pinned (n, T) state from current outer-ghost storage, but rebuild energy
    // with halved velocity. Densities & temperatures stay frozen at Model C7 values.
    auto refresh = [&](Vec& ob, float v_target, float u_target) {
        const float ni_B = ob(CNI) / m_i;
        const float nn_B = ob(CNN) / m_n;
        const float vv_B_old = ob(CNV) / (m_i * ni_B);
        const float uu_B_old = ob(CNU) / (m_n * nn_B);
        const float phig = 0.5f * (gPotential_imh(ns - 1) + gPotential_iph(ns - 1));
        const float Ti_B = (2.0f / 3.0f * ob(CEI) - 1.0f / 3.0f * m_i * ni_B * vv_B_old * vv_B_old
                            - 2.0f / 3.0f * m_i * ni_B * phig) / (2.0f * ni_B * k_b);
        const float Tn_B = (2.0f / 3.0f * ob(CEN) - 1.0f / 3.0f * m_n * nn_B * uu_B_old * uu_B_old
                            - 2.0f / 3.0f * m_n * nn_B * phig) / (nn_B * k_b);
        ob(CNV) = m_i * ni_B * v_target;
        ob(CNU) = m_n * nn_B * u_target;
        ob(CEI) = 1.5f * k_b * ni_B * 2.0f * Ti_B
                + 0.5f * m_i * ni_B * v_target * v_target + m_i * ni_B * phig;
        ob(CEN) = 1.5f * k_b * nn_B * Tn_B
                + 0.5f * m_n * nn_B * u_target * u_target + m_n * nn_B * phig;
    };

    const float v0 = 0.5f * vv_ns;
    const float u0 = 0.5f * uu_ns;
    refresh(outer_boundary0_i, v0,        u0);
    refresh(outer_boundary1_i, 0.5f * v0, 0.5f * u0);
}

} // namespace chromosphere
