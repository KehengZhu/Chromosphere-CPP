#include "chromosphere.hpp"

namespace chromosphere {

Vec cal_flux_state(const Grid& grid, const Vec& xn_state) {
    Vec F_state(arma::size(xn_state), arma::fill::zeros);
    const Vec rho_i  = get_scalar(grid, xn_state, cons::RHO_I);
    const Vec rho_n  = get_scalar(grid, xn_state, cons::RHO_N);
    const Vec rhoV_i = get_scalar(grid, xn_state, cons::MOM_I);
    const Vec rhoU_n = get_scalar(grid, xn_state, cons::MOM_N);
    const Vec e_i    = get_scalar(grid, xn_state, cons::E_I);
    const Vec e_n    = get_scalar(grid, xn_state, cons::E_N);

    const Vec V     = rhoV_i / rho_i;
    const Vec U     = rhoU_n / rho_n;
    const Vec phi_g = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    const Vec p_i   = 2.0/3.0 * e_i - 1.0/3.0 * rho_i % V % V - 2.0/3.0 * rho_i % phi_g;
    const Vec p_n   = 2.0/3.0 * e_n - 1.0/3.0 * rho_n % U % U - 2.0/3.0 * rho_n % phi_g;

    F_state += scalar_to(grid, rhoV_i,              cons::RHO_I);
    F_state += scalar_to(grid, rhoU_n,              cons::RHO_N);
    F_state += scalar_to(grid, rho_i % V % V + p_i, cons::MOM_I);
    F_state += scalar_to(grid, rho_n % U % U + p_n, cons::MOM_N);
    F_state += scalar_to(grid, (e_i + p_i) % V,     cons::E_I);
    F_state += scalar_to(grid, (e_n + p_n) % U,     cons::E_N);
    return F_state;
}

Vec cal_spectral_radius_state(const Grid& grid, const Vec& xn_state) {
    Vec res_state(arma::size(xn_state), arma::fill::zeros);
    const Vec rho_i  = get_scalar(grid, xn_state, cons::RHO_I);
    const Vec rho_n  = get_scalar(grid, xn_state, cons::RHO_N);
    const Vec rhoV_i = get_scalar(grid, xn_state, cons::MOM_I);
    const Vec rhoU_n = get_scalar(grid, xn_state, cons::MOM_N);
    const Vec e_i    = get_scalar(grid, xn_state, cons::E_I);
    const Vec e_n    = get_scalar(grid, xn_state, cons::E_N);

    const Vec V     = rhoV_i / rho_i;
    const Vec U     = rhoU_n / rho_n;
    const Vec phi_g = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    const Vec p_i   = 2.0/3.0 * e_i - 1.0/3.0 * rho_i % V % V - 2.0/3.0 * rho_i % phi_g;
    const Vec p_n   = 2.0/3.0 * e_n - 1.0/3.0 * rho_n % U % U - 2.0/3.0 * rho_n % phi_g;

    const Vec cs_i = arma::sqrt(grid.gamma_mono * p_i / rho_i);
    const Vec cs_n = arma::sqrt(grid.gamma_mono * p_n / rho_n);

    Vec spc_r_i = arma::max(arma::abs(V) + arma::abs(cs_i),
                            arma::abs(U) + arma::abs(cs_n));
    for (arma::uword i = 0; i < num_of_eq; ++i) {
        res_state += scalar_to(grid, spc_r_i, i);
    }
    return res_state;
}

// Explicit source (writeup §4.3):
//   S_E = [0, 0, p_i B D - ρ_i G, p_n B D - ρ_n G, 0, 0]^T
// with D = ∂(1/B)/∂s and G = ∂φ_g/∂s.
Vec cal_source_state(const Grid& grid, const Vec& xn_state) {
    Vec S_state(arma::size(xn_state), arma::fill::zeros);
    const Vec rho_i  = get_scalar(grid, xn_state, cons::RHO_I);
    const Vec rho_n  = get_scalar(grid, xn_state, cons::RHO_N);
    const Vec rhoV_i = get_scalar(grid, xn_state, cons::MOM_I);
    const Vec rhoU_n = get_scalar(grid, xn_state, cons::MOM_N);
    const Vec e_i    = get_scalar(grid, xn_state, cons::E_I);
    const Vec e_n    = get_scalar(grid, xn_state, cons::E_N);

    const Vec V     = rhoV_i / rho_i;
    const Vec U     = rhoU_n / rho_n;
    const Vec phi_g = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    const Vec p_i   = 2.0/3.0 * e_i - 1.0/3.0 * rho_i % V % V - 2.0/3.0 * rho_i % phi_g;
    const Vec p_n   = 2.0/3.0 * e_n - 1.0/3.0 * rho_n % U % U - 2.0/3.0 * rho_n % phi_g;

    // Gravity at cell center from the face potentials.
    const Vec gb = -(grid.phi_g_iph - grid.phi_g_imh) / grid.ds_i;

    S_state += scalar_to(grid, p_i % grid.B_i % grid.dinvB_ds_i + rho_i % gb, cons::MOM_I);
    S_state += scalar_to(grid, p_n % grid.B_i % grid.dinvB_ds_i + rho_n % gb, cons::MOM_N);
    return S_state;
}

} // namespace chromosphere
