#include "pfss_field_line.hpp"
#include "data_file_parser.hpp"
#include "scenario.hpp"

#include <cmath>
#include <cstdlib>
#include <string>

namespace chromosphere {

arma::uword pfss_peek_ns(const std::string& data_path) {
    return peek_ns_from_file(data_path);
}

Vec pfss_ic(Grid& grid, const std::string& data_path) {
    const ScenarioDataFile f = parse_scenario_data_file(data_path);
    if (f.ns != grid.ns)
        throw std::runtime_error("pfss_ic: data file ns=" + std::to_string(f.ns)
                                 + " differs from Grid ns=" + std::to_string(grid.ns));

    // Topology sets the outer-face closure. Both "closed" (half loop: outer face
    // is the apex, a reflecting symmetry plane) and "full" (whole loop: outer face
    // is footpoint B, a reflecting chromospheric wall) want a reflecting outer face;
    // the inner face is always reflecting. An "open" line keeps the coronal outflow.
    grid.outer_reflecting = (f.topology == "closed" || f.topology == "full");

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
        xn(arma::sub2ind(arma::size(grid.ns, num_of_eq), i, cons::E_E))
            = 1.5f * grid.k_b * n_e * T;       // electron internal energy, T_e = T
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
        ob(cons::E_E)   = 1.5f * grid.k_b * n_e * T_eff;       // T_e = T_eff at the ghost
    };

    const float phi_g_outer = grid.phi_g_iph(grid.ns - 1);
    const float phi_g_inner = grid.phi_g_imh(0);
    fill_ghost(grid.outer_boundary0_i, /*idx=*/0, phi_g_outer);
    fill_ghost(grid.outer_boundary1_i, /*idx=*/1, phi_g_outer);
    fill_ghost(grid.inner_boundary0_i, /*idx=*/2, phi_g_inner);
    fill_ghost(grid.inner_boundary1_i, /*idx=*/3, phi_g_inner);

    // --- Optional flare beam-heating overlay (Phase-3 de-risk) -------------
    // Opt-in via PFSS_FLARE=1: run the model_flare physics (Fisher, Canfield &
    // McClymont 1985 nonthermal-beam explosive evaporation) on THIS field line's
    // PFSS geometry + open BC instead of the uniform-B C7 column. Mirrors
    // model_flare_ic: force ionization + radiative cooling (the radiative sink
    // sets the explosive threshold), then add the support physics model_c7_ic
    // configures but pfss_ic otherwise omits — grid-scaled numerical diffusivity
    // (shock stabilization), TRAC (under-resolved TR), and the RTV (1978) coronal
    // conductive flux q(T) — and switch the beam on. All knobs share the FLARE_*
    // env names with model_flare_ic so one build serves both. Default off ⇒ the
    // plain pfss_field_line scenario (and the test suite) is unchanged.
    //
    // NOTE (de-risk shortcut): the C7 background here is NOT re-calibrated to be a
    // Stage-E fixed point (model_c7_ic does that by inverting P_phot over the
    // active network; pfss_ic uses the scalar fallback ~1e-4 s⁻¹). Over the ~15 s
    // beam-driven run the resulting quiet-background ionization drift is ~0.1% and
    // negligible against the beam, but the per-line ensemble (Phase 2) will need
    // the full height-dependent calibration to keep non-flaring lines steady.
    auto env_f = [](const char* key, float fallback) -> float {
        if (const char* e = std::getenv(key)) {
            try { return std::stof(std::string(e)); } catch (...) {}
        }
        return fallback;
    };
    if (env_f("PFSS_FLARE", 0.0f) != 0.0f) {
        grid.enable_ionization        = true;
        grid.enable_radiative_cooling = true;

        const float ds_mean = arma::mean(grid.ds_i);
        grid.numerical_diffusivity = 2.0e3f * ds_mean * env_f("FLARE_DIFF_MULT", 1.0f);
        grid.enable_trac  = true;
        grid.trac_T_chrom = 2.0e4f;
        {   // RTV F_c = (2/7) κ₀ T_cor^{7/2} / L (quiet-Sun defaults, as model_c7_ic).
            // Per-line L from the PFSS trace is a Phase-2 refinement.
            const float kappa0 = 1.0e-11f, T_cor = 8.0e5f, L_cor = 4.5e7f;
            grid.outer_heat_flux        = (2.0f / 7.0f) * kappa0 * std::pow(T_cor, 3.5f) / L_cor;
            grid.impose_outer_heat_flux = true;
        }

        grid.enable_beam_heating = true;
        grid.beam_flux        = env_f("FLARE_BEAM_FLUX", 5.0e7f);  // W/m² (explosive default)
        grid.beam_t_on        = env_f("FLARE_T_ON",   2.0f);
        grid.beam_duration    = env_f("FLARE_DUR",   10.0f);
        grid.beam_ramp        = env_f("FLARE_RAMP",   1.0f);
        grid.beam_h_lo_km     = env_f("FLARE_H_LO", 1600.0f);
        grid.beam_h_hi_km     = env_f("FLARE_H_HI", 2000.0f);
        // Self-consistent thick-target beam (Emslie 1978): inject at the loop apex
        // and attenuate down each leg by collisional stopping, with the stopping
        // depth rising as evaporation fills the loop. On by default for the loop
        // (the fixed window is unphysical here); FLARE_THICK_TARGET=0 reverts to it.
        grid.beam_thick_target = (env_f("FLARE_THICK_TARGET", 1.0f) != 0.0f);
        grid.beam_E_cut_keV    = env_f("FLARE_E_CUT", 20.0f);   // low-energy cutoff [keV]
        grid.beam_delta        = env_f("FLARE_DELTA",  5.0f);   // spectral index δ
        // Informational here: pfss_update_bc's apply_open_bcs is already a
        // transparent (Neumann) outflow, so the supersonic upflow can exit
        // regardless. Set for parity with model_flare and any future BC swap.
        grid.outer_free_outflow = (env_f("FLARE_FREE_OUTFLOW", 1.0f) != 0.0f);
    }

    grid.broadcast();
    return xn;
}

void pfss_update_bc(Grid& grid, const Vec& xn) {
    apply_open_bcs(grid, xn);
}

} // namespace chromosphere
