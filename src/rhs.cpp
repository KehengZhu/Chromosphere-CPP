#include "chromosphere.hpp"
#include "physics.hpp"

#include <cstdlib>
#include <iostream>

namespace chromosphere {

// ============================================================================
// Explicit RHS: TVD-MUSCL reconstruction + predictor-corrector +
// Rusanov flux differencing, minus the explicit source.
// (1.5D field-aligned; z-direction not considered. CoMFi Re_MUSCL analogue.)
// ============================================================================

Vec rhs_explicit_state(const Grid& grid, const Vec& xn_state) {
    // ---- shifted neighbours (conserved) ---------------------------------
    Vec cons_xn_state_ip1 = ip1(grid, xn_state);
    Vec cons_xn_state_im1 = im1(grid, xn_state);
    Vec cons_xn_state_ip2 = ip2(grid, xn_state);
    Vec cons_xn_state_im2 = im2(grid, xn_state);

    // ---- limit on primitive variables -----------------------------------
    Vec prim_xn_state     = cons2prim(grid, xn_state);
    Vec prim_xn_state_ip1 = cons2prim(grid, cons_xn_state_ip1);
    Vec prim_xn_state_im1 = cons2prim(grid, cons_xn_state_im1);
    Vec prim_xn_state_ip2 = cons2prim(grid, cons_xn_state_ip2);
    Vec prim_xn_state_im2 = cons2prim(grid, cons_xn_state_im2);

    const Vec dxn_state_iph = prim_xn_state_ip1 - prim_xn_state;
    const Vec dxn_state_imh = prim_xn_state     - prim_xn_state_im1;

    if ((dxn_state_iph.has_nan() + dxn_state_imh.has_nan()) > 0) {
        std::cout << "XN" << std::endl;
        print_xn(grid, xn_state);
        std::cout << "XNPRIM" << std::endl;
        print_xn(grid, prim_xn_state);
        std::cout << "dxn_iph has nan: " << dxn_state_iph.has_nan() << std::endl
                  << "dxn_imh has nan: " << dxn_state_imh.has_nan() << std::endl;
        std::exit(1);
    }

    // limiter inputs
    Vec r_state     = dxn_state_imh / dxn_state_iph;
    Vec r_state_ip1 = dxn_state_iph / (prim_xn_state_ip2 - prim_xn_state_ip1);
    Vec r_state_im1 = (prim_xn_state_im1 - prim_xn_state_im2) / dxn_state_imh;

    // extrapolated cell-edge variables (MUSCL paper eq 4.5)
    Vec Rxn_state_iph = prim_xn_state_ip1 - 0.5 * flux_lim(r_state_ip1) % (prim_xn_state_ip2 - prim_xn_state_ip1);
    Vec Lxn_state_iph = prim_xn_state     + 0.5 * flux_lim(r_state)     % (prim_xn_state_ip1 - prim_xn_state);
    Vec Rxn_state_imh = prim_xn_state     - 0.5 * flux_lim(r_state)     % (prim_xn_state_ip1 - prim_xn_state);
    Vec Lxn_state_imh = prim_xn_state_im1 + 0.5 * flux_lim(r_state_im1) % (prim_xn_state     - prim_xn_state_im1);

    if ((Lxn_state_iph.has_nan() + Rxn_state_iph.has_nan() +
         Lxn_state_imh.has_nan() + Rxn_state_imh.has_nan()) > 0) {
        std::cout << "Lxn_iph has nan: " << Lxn_state_iph.has_nan() << std::endl
                  << "Rxn_iph has nan: " << Rxn_state_iph.has_nan() << std::endl
                  << "Lxn_imh has nan: " << Lxn_state_imh.has_nan() << std::endl
                  << "Rxn_imh has nan: " << Rxn_state_imh.has_nan() << std::endl;
        std::exit(1);
    }

    Rxn_state_iph = prim2cons(grid, Rxn_state_iph);
    Lxn_state_iph = prim2cons(grid, Lxn_state_iph);
    Rxn_state_imh = prim2cons(grid, Rxn_state_imh);
    Lxn_state_imh = prim2cons(grid, Lxn_state_imh);

    // ---- prediction step (half-step) ------------------------------------
    Vec prim_xt_state     = cons2prim(grid,
        xn_state - grid.dt_state / grid.ds_state %
        (cal_flux_state(grid, Lxn_state_iph) - cal_flux_state(grid, Rxn_state_imh)));
    Vec prim_xt_state_ip1 = cons2prim(grid, ip1(grid, prim2cons(grid, prim_xt_state)));
    Vec prim_xt_state_im1 = cons2prim(grid, im1(grid, prim2cons(grid, prim_xt_state)));

    // ---- 2nd-order-in-time reconstruction -------------------------------
    Rxn_state_iph = 0.5 * (prim_xn_state_ip1 + prim_xt_state_ip1) - 0.5 * flux_lim(r_state_ip1) % (prim_xn_state_ip2 - prim_xn_state_ip1);
    Lxn_state_iph = 0.5 * (prim_xn_state     + prim_xt_state)     + 0.5 * flux_lim(r_state)     % (prim_xn_state_ip1 - prim_xn_state);
    Rxn_state_imh = 0.5 * (prim_xn_state     + prim_xt_state)     - 0.5 * flux_lim(r_state)     % (prim_xn_state_ip1 - prim_xn_state);
    Lxn_state_imh = 0.5 * (prim_xn_state_im1 + prim_xt_state_im1) + 0.5 * flux_lim(r_state_im1) % (prim_xn_state     - prim_xn_state_im1);

    Rxn_state_iph = prim2cons(grid, Rxn_state_iph);
    Lxn_state_iph = prim2cons(grid, Lxn_state_iph);
    Rxn_state_imh = prim2cons(grid, Rxn_state_imh);
    Lxn_state_imh = prim2cons(grid, Lxn_state_imh);

    // ---- Rusanov / local Lax–Friedrichs flux ----------------------------
    const Vec a_state_imh = arma::max(cal_spectral_radius_state(grid, Lxn_state_imh),
                                      cal_spectral_radius_state(grid, Rxn_state_imh));
    const Vec a_state_iph = arma::max(cal_spectral_radius_state(grid, Lxn_state_iph),
                                      cal_spectral_radius_state(grid, Rxn_state_iph));

    const Vec F_state_imh = 0.5 * (cal_flux_state(grid, Lxn_state_imh) + cal_flux_state(grid, Rxn_state_imh)
                                   - a_state_imh % (Rxn_state_imh - Lxn_state_imh));
    const Vec F_state_iph = 0.5 * (cal_flux_state(grid, Lxn_state_iph) + cal_flux_state(grid, Rxn_state_iph)
                                   - a_state_iph % (Rxn_state_iph - Lxn_state_iph));

    if ((a_state_imh.has_nan() + a_state_iph.has_nan() +
         F_state_imh.has_nan() + F_state_iph.has_nan()) > 0) {
        std::cout << "a_imh has nan: " << a_state_imh.has_nan() << std::endl
                  << "a_iph has nan: " << a_state_iph.has_nan() << std::endl
                  << "Fximh has nan: " << F_state_imh.has_nan() << std::endl
                  << "Fxiph has nan: " << F_state_iph.has_nan() << std::endl;
        std::exit(1);
    }

    return -1.0 * grid.B_state % (F_state_iph / grid.B_state_iph - F_state_imh / grid.B_state_imh) / grid.ds_state
           + cal_source_state(grid, xn_state);
}

// ============================================================================
// Implicit RHS (writeup §4.3, R_I): ion-neutral drag, collisional + frictional
// heating, conservative field-aligned heat conduction. Pressure-area and
// gravity live in cal_source_state (R_E).
// ============================================================================

Vec rhs_implicit_state(const Grid& grid, const Vec& xn_state) {
    Vec RI_state(arma::size(xn_state), arma::fill::zeros);

    // ---- cell-centered primitives ---------------------------------------
    const Vec rho_i  = get_scalar(grid, xn_state, cons::RHO_I);
    const Vec rho_n  = get_scalar(grid, xn_state, cons::RHO_N);
    const Vec rhoV_i = get_scalar(grid, xn_state, cons::MOM_I);
    const Vec rhoU_n = get_scalar(grid, xn_state, cons::MOM_N);
    const Vec e_i    = get_scalar(grid, xn_state, cons::E_I);
    const Vec e_n    = get_scalar(grid, xn_state, cons::E_N);

    const Vec n_i   = rho_i / grid.m_i;
    const Vec n_n   = rho_n / grid.m_n;
    const Vec V     = rhoV_i / rho_i;
    const Vec U     = rhoU_n / rho_n;
    const Vec phi_g = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    const Vec p_i   = 2.0/3.0 * e_i - 1.0/3.0 * rho_i % V % V - 2.0/3.0 * rho_i % phi_g;
    const Vec p_n   = 2.0/3.0 * e_n - 1.0/3.0 * rho_n % U % U - 2.0/3.0 * rho_n % phi_g;
    const Vec T_i   = p_i / (2.0 * n_i * grid.k_b);
    const Vec T_n   = p_n / (n_n * grid.k_b);

    const Vec Ke = kappa_e(n_i, n_n, T_i);
    const Vec Ki(arma::size(T_i), arma::fill::zeros);
    const Vec Kn = kappa_n(n_i, n_n, T_i, T_n);

    // ---- i+1 neighbour (uses outer ghost at i = ns-1) -------------------
    const Vec xn_state_ip1 = ip1(grid, xn_state);
    const Vec rho_i_ip1  = get_scalar(grid, xn_state_ip1, cons::RHO_I);
    const Vec rho_n_ip1  = get_scalar(grid, xn_state_ip1, cons::RHO_N);
    const Vec rhoV_i_ip1 = get_scalar(grid, xn_state_ip1, cons::MOM_I);
    const Vec rhoU_n_ip1 = get_scalar(grid, xn_state_ip1, cons::MOM_N);
    const Vec e_i_ip1    = get_scalar(grid, xn_state_ip1, cons::E_I);
    const Vec e_n_ip1    = get_scalar(grid, xn_state_ip1, cons::E_N);

    const Vec n_i_ip1 = rho_i_ip1 / grid.m_i;
    const Vec n_n_ip1 = rho_n_ip1 / grid.m_n;
    const Vec V_ip1   = rhoV_i_ip1 / rho_i_ip1;
    const Vec U_ip1   = rhoU_n_ip1 / rho_n_ip1;
    Vec phi_g_ip1 = ip1(grid, 0.5 * (grid.phi_g_iph + grid.phi_g_imh), SLICE);
    phi_g_ip1[grid.ns - 1] = grid.phi_g_iph[grid.ns - 1];
    const Vec p_i_ip1 = 2.0/3.0 * e_i_ip1 - 1.0/3.0 * rho_i_ip1 % V_ip1 % V_ip1 - 2.0/3.0 * rho_i_ip1 % phi_g_ip1;
    const Vec p_n_ip1 = 2.0/3.0 * e_n_ip1 - 1.0/3.0 * rho_n_ip1 % U_ip1 % U_ip1 - 2.0/3.0 * rho_n_ip1 % phi_g_ip1;
    const Vec T_i_ip1 = p_i_ip1 / (2.0 * n_i_ip1 * grid.k_b);
    const Vec T_n_ip1 = p_n_ip1 / (n_n_ip1 * grid.k_b);

    const Vec Ke_ip1 = kappa_e(n_i_ip1, n_n_ip1, T_i_ip1);
    const Vec Ki_ip1(arma::size(T_i), arma::fill::zeros);
    const Vec Kn_ip1 = kappa_n(n_i_ip1, n_n_ip1, T_i_ip1, T_n_ip1);

    // ---- i-1 neighbour (uses inner ghost at i = 0) ----------------------
    const Vec xn_state_im1 = im1(grid, xn_state);
    const Vec rho_i_im1  = get_scalar(grid, xn_state_im1, cons::RHO_I);
    const Vec rho_n_im1  = get_scalar(grid, xn_state_im1, cons::RHO_N);
    const Vec rhoV_i_im1 = get_scalar(grid, xn_state_im1, cons::MOM_I);
    const Vec rhoU_n_im1 = get_scalar(grid, xn_state_im1, cons::MOM_N);
    const Vec e_i_im1    = get_scalar(grid, xn_state_im1, cons::E_I);
    const Vec e_n_im1    = get_scalar(grid, xn_state_im1, cons::E_N);

    const Vec n_i_im1 = rho_i_im1 / grid.m_i;
    const Vec n_n_im1 = rho_n_im1 / grid.m_n;
    const Vec V_im1   = rhoV_i_im1 / rho_i_im1;
    const Vec U_im1   = rhoU_n_im1 / rho_n_im1;
    Vec phi_g_im1 = im1(grid, 0.5 * (grid.phi_g_iph + grid.phi_g_imh), SLICE);
    phi_g_im1[0] = grid.phi_g_imh[0];
    const Vec p_i_im1 = 2.0/3.0 * e_i_im1 - 1.0/3.0 * rho_i_im1 % V_im1 % V_im1 - 2.0/3.0 * rho_i_im1 % phi_g_im1;
    const Vec p_n_im1 = 2.0/3.0 * e_n_im1 - 1.0/3.0 * rho_n_im1 % U_im1 % U_im1 - 2.0/3.0 * rho_n_im1 % phi_g_im1;
    const Vec T_i_im1 = p_i_im1 / (2.0 * n_i_im1 * grid.k_b);
    const Vec T_n_im1 = p_n_im1 / (n_n_im1 * grid.k_b);

    const Vec Ke_im1 = kappa_e(n_i_im1, n_n_im1, T_i_im1);
    const Vec Ki_im1(arma::size(T_i), arma::fill::zeros);
    const Vec Kn_im1 = kappa_n(n_i_im1, n_n_im1, T_i_im1, T_n_im1);

    // ---- ion-neutral drag (α = ρ_i ν_in) --------------------------------
    const Vec nu_in_v = nu_in(grid, n_n, T_i, T_n);
    const Vec alpha   = rho_i % nu_in_v;
    const Vec w       = V - U;
    RI_state += scalar_to(grid, -alpha % w, cons::MOM_I);
    RI_state += scalar_to(grid,  alpha % w, cons::MOM_N);

    // ---- collisional + frictional heating (R_I energy rows) -------------
    RI_state += scalar_to(grid,
        alpha / (grid.m_i + grid.m_n) % (3.0 * grid.k_b * (T_n - T_i) + grid.m_n * w % w) - alpha % V % w,
        cons::E_I);
    RI_state += scalar_to(grid,
        alpha / (grid.m_i + grid.m_n) % (3.0 * grid.k_b * (T_i - T_n) + grid.m_i * w % w) + alpha % U % w,
        cons::E_N);

    // ---- conservative field-aligned heat conduction (writeup eq 4.4) ---
    //   C = (B/ds) [ K_{i+1/2}/B_{i+1/2} (T_{i+1}-T_i)/ds
    //              - K_{i-1/2}/B_{i-1/2} (T_i-T_{i-1})/ds ]
    const Vec Kei_iph = 0.5 * ((Ke + Ki) + (Ke_ip1 + Ki_ip1));
    const Vec Kei_imh = 0.5 * ((Ke + Ki) + (Ke_im1 + Ki_im1));
    const Vec Kn_iph  = 0.5 * (Kn + Kn_ip1);
    const Vec Kn_imh  = 0.5 * (Kn + Kn_im1);

    const Vec q_i_iph = Kei_iph / grid.B_iph % (T_i_ip1 - T_i)     / grid.ds_i;
    const Vec q_i_imh = Kei_imh / grid.B_imh % (T_i     - T_i_im1) / grid.ds_i;
    const Vec q_n_iph = Kn_iph  / grid.B_iph % (T_n_ip1 - T_n)     / grid.ds_i;
    const Vec q_n_imh = Kn_imh  / grid.B_imh % (T_n     - T_n_im1) / grid.ds_i;

    const Vec C_i = grid.B_i % (q_i_iph - q_i_imh) / grid.ds_i;
    const Vec C_n = grid.B_i % (q_n_iph - q_n_imh) / grid.ds_i;
    RI_state += scalar_to(grid, C_i, cons::E_I);
    RI_state += scalar_to(grid, C_n, cons::E_N);

    return RI_state;
}

} // namespace chromosphere
