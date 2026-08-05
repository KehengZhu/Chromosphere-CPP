#include "chromosphere.hpp"

#include <cmath>
#include <limits>

namespace chromosphere {

Vec cal_flux_state(const Grid& grid, const Vec& xn_state) {
    if (!grid.eos_gamma_table.empty()) {
        Vec flux(arma::size(xn_state), arma::fill::zeros);
        const auto sz = arma::size(grid.ns, num_of_eq);
        for (arma::uword i = 0; i < grid.ns; ++i) {
            auto at = [&](arma::uword k) -> double {
                return static_cast<double>(xn_state(arma::sub2ind(sz, i, k)));
            };
            const double phi = 0.5 * static_cast<double>(grid.phi_g_imh(i)
                                                       + grid.phi_g_iph(i));
            const MixtureThermo th = decode_equilibrium_mixture(
                grid.eos_gamma_table, at(cons::RHO_I), at(cons::RHO_N),
                at(cons::MOM_I), at(cons::MOM_N), at(cons::E_I), at(cons::E_N),
                phi, grid.eos_temperature_hint(i),
                grid.eos_gamma_debug_clamp);
            grid.store_eos_temperature_hint(i, th.T);
            const double v = (at(cons::MOM_I) + at(cons::MOM_N)) / th.rho;
            const double f[7] = {
                at(cons::RHO_I)*v, at(cons::RHO_N)*v,
                at(cons::RHO_I)*v*v + th.p_i,
                at(cons::RHO_N)*v*v + th.p_n,
                (at(cons::E_I)+th.p_i)*v,
                (at(cons::E_N)+th.p_n)*v,
                at(cons::E_E)*v};
            for (arma::uword k = 0; k < num_of_eq; ++k)
                flux(arma::sub2ind(sz, i, k)) = static_cast<float>(f[k]);
        }
        return flux;
    }
    Vec F_state(arma::size(xn_state), arma::fill::zeros);
    const Vec rho_i  = get_scalar(grid, xn_state, cons::RHO_I);
    const Vec rho_n  = get_scalar(grid, xn_state, cons::RHO_N);
    const Vec rhoV_i = get_scalar(grid, xn_state, cons::MOM_I);
    const Vec rhoU_n = get_scalar(grid, xn_state, cons::MOM_N);
    const Vec e_i    = get_scalar(grid, xn_state, cons::E_I);
    const Vec e_n    = get_scalar(grid, xn_state, cons::E_N);
    const Vec e_e    = get_scalar(grid, xn_state, cons::E_E);

    const Vec V     = rhoV_i / rho_i;
    const Vec U     = rhoU_n / rho_n;
    const Vec phi_g = 0.5 * (grid.phi_g_imh + grid.phi_g_iph);
    // p_i = total charged pressure (protons + electrons): the momentum flux and
    // the charged total-energy flux are unchanged by the T_e split.
    const Vec p_i   = grid.gm1() * e_i - grid.half_gm1() * rho_i % V % V - grid.gm1() * rho_i % phi_g;
    const Vec p_n   = grid.gm1() * e_n - grid.half_gm1() * rho_n % U % U - grid.gm1() * rho_n % phi_g;

    F_state += scalar_to(grid, rhoV_i,              cons::RHO_I);
    F_state += scalar_to(grid, rhoU_n,              cons::RHO_N);
    F_state += scalar_to(grid, rho_i % V % V + p_i, cons::MOM_I);
    F_state += scalar_to(grid, rho_n % U % U + p_n, cons::MOM_N);
    F_state += scalar_to(grid, (e_i + p_i) % V,     cons::E_I);
    F_state += scalar_to(grid, (e_n + p_n) % U,     cons::E_N);
    // Electron internal energy advects at the charged-fluid velocity V (electron
    // inertia negligible). The −p_e ∇·V compression work is an explicit source
    // (cal_source_state); here just the advective flux ε_e V.
    F_state += scalar_to(grid, e_e % V,             cons::E_E);
    return F_state;
}

Vec cal_spectral_radius_state(const Grid& grid, const Vec& xn_state) {
    if (!grid.eos_gamma_table.empty()) {
        Vec result(arma::size(xn_state), arma::fill::zeros);
        const auto sz = arma::size(grid.ns, num_of_eq);
        for (arma::uword i = 0; i < grid.ns; ++i) {
            auto at = [&](arma::uword k) -> double {
                return static_cast<double>(xn_state(arma::sub2ind(sz, i, k)));
            };
            const double phi = 0.5 * static_cast<double>(grid.phi_g_imh(i)
                                                       + grid.phi_g_iph(i));
            const MixtureThermo th = decode_equilibrium_mixture(
                grid.eos_gamma_table, at(cons::RHO_I), at(cons::RHO_N),
                at(cons::MOM_I), at(cons::MOM_N), at(cons::E_I), at(cons::E_N),
                phi, grid.eos_temperature_hint(i),
                grid.eos_gamma_debug_clamp);
            grid.store_eos_temperature_hint(i, th.T);
            const double v = (at(cons::MOM_I) + at(cons::MOM_N)) / th.rho;
            const float a = static_cast<float>(std::abs(v) + std::sqrt(
                th.gamma1 * (th.p_i + th.p_n) / th.rho));
            for (arma::uword k = 0; k < num_of_eq; ++k)
                result(arma::sub2ind(sz, i, k)) = a;
        }
        return result;
    }
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
    Vec p_i   = grid.gm1() * e_i - grid.half_gm1() * rho_i % V % V - grid.gm1() * rho_i % phi_g;
    Vec p_n   = grid.gm1() * e_n - grid.half_gm1() * rho_n % U % U - grid.gm1() * rho_n % phi_g;

    // Flare scenario: in the fast (~10²–10³ km/s), hot (~10⁷ K) evaporated column
    // the float32 cancellation p = ⅔E − ⅓ρv² − ⅔ρφ can leave p slightly NEGATIVE,
    // making c_s = √(γp/ρ) a NaN that poisons the Rusanov flux. Clamp p ≥ 0 for
    // the sound speed only. Gated on beam heating so steady runs / tests are unchanged.
    if (grid.enable_beam_heating) {
        p_i = arma::clamp(p_i, 0.0f, arma::datum::inf);
        p_n = arma::clamp(p_n, 0.0f, arma::datum::inf);
    }

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
    const Vec p_i   = grid.gm1() * e_i - grid.half_gm1() * rho_i % V % V - grid.gm1() * rho_i % phi_g;
    const Vec p_n   = grid.gm1() * e_n - grid.half_gm1() * rho_n % U % U - grid.gm1() * rho_n % phi_g;

    // Gravity at cell center from the face potentials.
    const Vec gb = -(grid.phi_g_iph - grid.phi_g_imh) / grid.ds_i;

    S_state += scalar_to(grid, p_i % grid.B_i % grid.dinvB_ds_i + rho_i % gb, cons::MOM_I);
    S_state += scalar_to(grid, p_n % grid.B_i % grid.dinvB_ds_i + rho_n % gb, cons::MOM_N);

    // Electron compression work −p_e ∇·V (three-temperature model only). In the
    // 1.5D field-aligned geometry (cross-section A ∝ 1/B) the velocity divergence
    // consistent with the flux differencing is ∇·V = B [V_{i+½}/B_{i+½} −
    // V_{i-½}/B_{i-½}]/Δs with face velocities V_{i±½} the cell-center averages.
    // Adding this to ε_e (advected as ε_e V in cal_flux_state) gives the correct
    // adiabatic electron response dε_e/dt = −(5/2) p_e ∇·V. The complementary
    // proton compression is automatic: p_proton = p_total − p_e and p_total
    // (from the conservatively-evolved E_I) carries the full −p_total ∇·V work.
    // Gated on enable_Te so the single-temperature baseline is untouched.
    if (grid.enable_Te) {
        const Vec V     = rhoV_i / rho_i;
        const Vec e_e   = get_scalar(grid, xn_state, cons::E_E);
        const Vec p_e   = grid.gm1() * e_e;
        const Vec V_ip1 = ip1(grid, V, SLICE);
        const Vec V_im1 = im1(grid, V, SLICE);
        // Face velocity by interpolation between cell centers. Uniform mesh: the
        // midpoint average. Refined mesh: distance-weighted so the face value is a
        // proper linear interpolant (metric-weighted), keeping ∇·V second-order.
        Vec Vf_iph, Vf_imh;
        if (grid.uniform_mesh) {
            Vf_iph = 0.5 * (V + V_ip1);
            Vf_imh = 0.5 * (V_im1 + V);
        } else {
            const Vec ds_ip1 = ip1(grid, grid.ds_i, SLICE);
            const Vec ds_im1 = im1(grid, grid.ds_i, SLICE);
            Vf_iph = (0.5 * ds_ip1 % V + 0.5 * grid.ds_i % V_ip1) / grid.ds_iph_i;
            Vf_imh = (0.5 * grid.ds_i % V_im1 + 0.5 * ds_im1 % V) / grid.ds_imh_i;
        }
        const Vec divV = grid.B_i % (Vf_iph / grid.B_iph - Vf_imh / grid.B_imh) / grid.ds_i;
        S_state += scalar_to(grid, -p_e % divV, cons::E_E);
    }
    return S_state;
}

} // namespace chromosphere
