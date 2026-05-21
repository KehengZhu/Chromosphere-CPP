#include "pfss_field_line.hpp"
#include "data_file_parser.hpp"
#include "scenario.hpp"

namespace chromosphere {

arma::uword pfss_peek_ns(const std::string& data_path) {
    return peek_ns_from_file(data_path);
}

Vec pfss_ic(Grid& grid, const std::string& data_path) {
    const ScenarioDataFile f = parse_scenario_data_file(data_path);
    if (f.ns != grid.ns)
        throw std::runtime_error("pfss_ic: data file ns=" + std::to_string(f.ns)
                                 + " differs from Grid ns=" + std::to_string(grid.ns));

    // Field-line geometry from the file.
    for (arma::uword i = 0; i < grid.ns; ++i) {
        grid.ds_i  (i) = f.ds_m   (i);
        grid.B_imh (i) = f.B_imh_T(i);
        grid.B_iph (i) = f.B_iph_T(i);
        grid.B_i   (i) = 0.5f * (f.B_imh_T(i) + f.B_iph_T(i));
        grid.phi_g_imh(i) = f.phi_g_imh(i);
        grid.phi_g_iph(i) = f.phi_g_iph(i);
        grid.dinvB_ds_i(i) = (1.0f / f.B_iph_T(i) - 1.0f / f.B_imh_T(i)) / f.ds_m(i);
    }

    // Build the conserved-variable state. Energy includes ρ·φ_g so that
    // cons2prim (state.cpp:108) recovers the right pressure.
    Vec xn = arma::zeros<Vec>(grid.n_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float n_e   = f.ne_im3(i);
        const float n_n   = f.nn_im3(i);
        const float T     = f.T_K  (i);
        const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_I)) = n_e * grid.m_i;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::RHO_N)) = n_n * grid.m_n;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::MOM_I)) = 0.0f;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::MOM_N)) = 0.0f;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_I))
            = 1.5f * grid.k_b * n_e * 2.0f * T + n_e * grid.m_i * phi_g;
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_N))
            = 1.5f * grid.k_b * n_n * T       + n_n * grid.m_n * phi_g;
    }

    // Ghost cells. φ_g convention matches model_c7_ic (and rhs.cpp:154,178):
    // outer ghost uses phi_g_iph[ns-1], inner ghost uses phi_g_imh[0]. The
    // T_e_factor column scales the ghost's effective temperature (model_c7
    // uses 2.0 at the outer ghost to mimic coronal T_e).
    auto fill_ghost = [&](Vec& ob, arma::uword idx, float phi_g_eff) {
        const float n_e    = f.ghost_ne_im3    (idx);
        const float n_n    = f.ghost_nn_im3    (idx);
        const float T_eff  = f.ghost_T_K       (idx) * f.ghost_T_e_factor(idx);
        ob(cons::RHO_I) = n_e * grid.m_i;
        ob(cons::RHO_N) = n_n * grid.m_n;
        ob(cons::MOM_I) = 0.0f;
        ob(cons::MOM_N) = 0.0f;
        ob(cons::E_I)   = 1.5f * grid.k_b * n_e * 2.0f * T_eff + n_e * grid.m_i * phi_g_eff;
        ob(cons::E_N)   = 1.5f * grid.k_b * n_n * T_eff       + n_n * grid.m_n * phi_g_eff;
    };

    const float phi_g_outer = grid.phi_g_iph(grid.ns - 1);
    const float phi_g_inner = grid.phi_g_imh(0);
    fill_ghost(grid.outer_boundary0_i, /*idx=*/0, phi_g_outer);
    fill_ghost(grid.outer_boundary1_i, /*idx=*/1, phi_g_outer);
    fill_ghost(grid.inner_boundary0_i, /*idx=*/2, phi_g_inner);
    fill_ghost(grid.inner_boundary1_i, /*idx=*/3, phi_g_inner);

    grid.broadcast();
    return xn;
}

void pfss_update_bc(Grid& grid, const Vec& xn) {
    apply_open_bcs(grid, xn);
}

} // namespace chromosphere
