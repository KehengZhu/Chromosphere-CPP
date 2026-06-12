#include "chromosphere.hpp"

namespace chromosphere {

// ============================================================================
// Packed-state accessors
// ============================================================================

Vec get_scalar(const Grid& grid, const Vec& xn_state, arma::uword index) {
    Vec scalar(grid.ns);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        scalar(i) = xn_state(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, index));
    }
    return scalar;
}

Vec scalar_to(const Grid& grid, const Vec& xn_i, arma::uword index) {
    Vec xn_state = arma::zeros<Vec>(grid.n_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        xn_state(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, index)) = xn_i(i);
    }
    return xn_state;
}

// ============================================================================
// Index shifts. ip1(xn)(i) = xn(i+1) in the interior; the i = ns-1 slot
// is filled from outer_boundary0_i (or kept Neumann if USE_NEUMANN_BC).
// Likewise im1 fills i = 0 from inner_boundary0_i.
// ============================================================================

Vec ip1(const Grid& grid, const Vec& xn, arma::uword nk) {
    Vec xn_ip1(arma::size(xn));
#pragma omp parallel for collapse(2)
    for (arma::uword i = 0; i < grid.ns - 1; ++i) {
        for (arma::uword k = 0; k < nk; ++k) {
            xn_ip1(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, k)) =
                xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i + 1, k));
        }
    }
    for (arma::uword k = 0; k < nk; ++k) {
        xn_ip1(arma::sub2ind(arma::size(grid.ns, num_of_eq), grid.ns - 1, k)) =
            grid.outer_boundary0_i[k];
        if (USE_NEUMANN_BC || nk == SLICE) {
            xn_ip1(arma::sub2ind(arma::size(grid.ns, num_of_eq), grid.ns - 1, k)) =
                xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), grid.ns - 1, k));
        }
    }
    return xn_ip1;
}

Vec im1(const Grid& grid, const Vec& xn, arma::uword nk) {
    Vec xn_im1(arma::size(xn));
#pragma omp parallel for collapse(2)
    for (arma::uword i = 1; i < grid.ns; ++i) {
        for (arma::uword k = 0; k < nk; ++k) {
            xn_im1(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, k)) =
                xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i - 1, k));
        }
    }
    for (arma::uword k = 0; k < nk; ++k) {
        xn_im1(arma::sub2ind(arma::size(grid.ns, num_of_eq), 0, k)) =
            grid.inner_boundary0_i[k];
        if (USE_NEUMANN_BC || nk == SLICE) {
            xn_im1(arma::sub2ind(arma::size(grid.ns, num_of_eq), 0, k)) =
                xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), 0, k));
        }
    }
    return xn_im1;
}

Vec ip2(const Grid& grid, const Vec& xn, arma::uword nk) {
    Vec xn_ip2 = ip1(grid, ip1(grid, xn, nk), nk);
    if (USE_NEUMANN_BC || nk == SLICE) return xn_ip2;
    for (arma::uword k = 0; k < nk; ++k) {
        xn_ip2(arma::sub2ind(arma::size(grid.ns, num_of_eq), grid.ns - 1, k)) =
            grid.outer_boundary1_i[k];
    }
    return xn_ip2;
}

Vec im2(const Grid& grid, const Vec& xn, arma::uword nk) {
    Vec xn_im2 = im1(grid, im1(grid, xn, nk), nk);
    if (USE_NEUMANN_BC || nk == SLICE) return xn_im2;
    for (arma::uword k = 0; k < nk; ++k) {
        xn_im2(arma::sub2ind(arma::size(grid.ns, num_of_eq), 0, k)) =
            grid.inner_boundary1_i[k];
    }
    return xn_im2;
}

// ============================================================================
// Equation of state: cons ↔ prim (writeup eq 38, electron quasi-neutrality
// folds T_e = T_i and n_e = n_i into the ion pressure with an extra factor of 2).
// ============================================================================

Vec cons2prim(const Grid& grid, const Vec& cons_state) {
    Vec prim_state(arma::size(cons_state), arma::fill::zeros);
    const Vec rho_i  = get_scalar(grid, cons_state, cons::RHO_I);
    const Vec rho_n  = get_scalar(grid, cons_state, cons::RHO_N);
    const Vec rhoV_i = get_scalar(grid, cons_state, cons::MOM_I);
    const Vec rhoU_n = get_scalar(grid, cons_state, cons::MOM_N);
    const Vec e_i    = get_scalar(grid, cons_state, cons::E_I);
    const Vec e_n    = get_scalar(grid, cons_state, cons::E_N);
    const Vec e_e    = get_scalar(grid, cons_state, cons::E_E);

    const Vec V     = rhoV_i / rho_i;
    const Vec U     = rhoU_n / rho_n;
    const Vec phi_g = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    // p_i is the TOTAL charged pressure (protons + electrons); see cons::E_I note.
    const Vec p_i   = 2.0/3.0 * e_i - 1.0/3.0 * rho_i % V % V - 2.0/3.0 * rho_i % phi_g;
    const Vec p_n   = 2.0/3.0 * e_n - 1.0/3.0 * rho_n % U % U - 2.0/3.0 * rho_n % phi_g;
    const Vec p_e   = 2.0/3.0 * e_e;   // electron partial pressure (no KE / gravity)

    prim_state += scalar_to(grid, rho_i, prim::RHO_I);
    prim_state += scalar_to(grid, rho_n, prim::RHO_N);
    prim_state += scalar_to(grid, V,     prim::V);
    prim_state += scalar_to(grid, U,     prim::U);
    prim_state += scalar_to(grid, p_i,   prim::P_I);
    prim_state += scalar_to(grid, p_n,   prim::P_N);
    prim_state += scalar_to(grid, p_e,   prim::P_E);
    return prim_state;
}

Vec prim2cons(const Grid& grid, const Vec& prim_state) {
    Vec cons_state(arma::size(prim_state), arma::fill::zeros);
    const Vec rho_i = get_scalar(grid, prim_state, prim::RHO_I);
    const Vec rho_n = get_scalar(grid, prim_state, prim::RHO_N);
    const Vec V     = get_scalar(grid, prim_state, prim::V);
    const Vec U     = get_scalar(grid, prim_state, prim::U);
    const Vec p_i   = get_scalar(grid, prim_state, prim::P_I);
    const Vec p_n   = get_scalar(grid, prim_state, prim::P_N);
    const Vec p_e   = get_scalar(grid, prim_state, prim::P_E);

    const Vec rhoV_i = rho_i % V;
    const Vec rhoU_n = rho_n % U;
    const Vec phi_g  = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    // E_I is the TOTAL charged energy (p_i already includes the electron pressure
    // p_e); E_E carries the electron internal energy alone. See cons::E_I note.
    const Vec e_i    = 3.0/2.0 * p_i + 0.5 * rho_i % V % V + rho_i % phi_g;
    const Vec e_n    = 3.0/2.0 * p_n + 0.5 * rho_n % U % U + rho_n % phi_g;
    const Vec e_e    = 3.0/2.0 * p_e;

    cons_state += scalar_to(grid, rho_i,  cons::RHO_I);
    cons_state += scalar_to(grid, rho_n,  cons::RHO_N);
    cons_state += scalar_to(grid, rhoV_i, cons::MOM_I);
    cons_state += scalar_to(grid, rhoU_n, cons::MOM_N);
    cons_state += scalar_to(grid, e_i,    cons::E_I);
    cons_state += scalar_to(grid, e_n,    cons::E_N);
    cons_state += scalar_to(grid, e_e,    cons::E_E);
    return cons_state;
}

// ============================================================================
// MUSCL flux limiter (minmod; writeup eq 15)
// ============================================================================

Vec flux_lim(const Vec& r) {
    Vec one (arma::size(r), arma::fill::ones);
    Vec zero(arma::size(r), arma::fill::zeros);
    Vec res = arma::max(zero, arma::min(one, r));
    res(arma::find_nan(res)).zeros();
    return res;
}

} // namespace chromosphere
