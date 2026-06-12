#include "scenario.hpp"
#include "analytic_canopy.hpp"
#include "model_c7.hpp"
#include "model_flare.hpp"
#include "pfss_field_line.hpp"

#include <stdexcept>
#include <cstdlib>

namespace chromosphere {

Scenario make_scenario(const std::string& name, const std::string& data_path) {
    Scenario sc;
    sc.name = name;

    if (name == "model_c7") {
        sc.peek_ns   = []() -> arma::uword { return 100; };
        sc.ic        = [](Grid& g) { return model_c7_ic(g); };
        sc.update_bc = [](Grid& g, const Vec& xn) { model_c7_update_bc(g, xn); };
        return sc;
    }
    if (name == "model_flare") {
        // Grid resolution is env-tunable (FLARE_NS, default 100) — a finer grid
        // resolves the steep condensation-front contact discontinuity and shrinks
        // the single-cell T = p/n overheating there. Matches the FLARE_* knobs.
        sc.peek_ns   = []() -> arma::uword {
            arma::uword ns = 100;
            if (const char* e = std::getenv("FLARE_NS")) { try { ns = std::stoul(e); } catch (...) {} }
            return ns;
        };
        sc.ic        = [](Grid& g) { return model_flare_ic(g); };
        sc.update_bc = [](Grid& g, const Vec& xn) { model_flare_update_bc(g, xn); };
        return sc;
    }
    if (name == "analytic_canopy") {
        sc.peek_ns   = []() { return analytic_canopy_peek_ns(); };
        sc.ic        = [](Grid& g) { return analytic_canopy_ic(g); };
        sc.update_bc = [](Grid& g, const Vec& xn) { analytic_canopy_update_bc(g, xn); };
        return sc;
    }
    if (name == "pfss_field_line") {
        if (data_path.empty())
            throw std::runtime_error("make_scenario: pfss_field_line requires a data_path "
                                     "(6th positional argument to chromo_main)");
        sc.peek_ns   = [data_path]() { return pfss_peek_ns(data_path); };
        sc.ic        = [data_path](Grid& g) { return pfss_ic(g, data_path); };
        sc.update_bc = [](Grid& g, const Vec& xn) { pfss_update_bc(g, xn); };
        return sc;
    }

    throw std::runtime_error("make_scenario: unknown scenario name '" + name + "'");
}

// Reflecting wall at the inner face + zero-gradient outflow at the outer face.
// See scenario.hpp for the BC math; the implementation just reads cells 0, 1,
// and ns-1 out of the packed state, mirrors / copies them into the four ghost
// buffers, and re-runs Grid::broadcast so the packed _state caches stay in sync.
void apply_open_bcs(Grid& grid, const Vec& xn) {
    const arma::uword ns = grid.ns;
    const float k_b = grid.k_b;
    const float m_i = grid.m_i;
    const float m_n = grid.m_n;

    auto cell_prim = [&](arma::uword i,
                         float& rho_i, float& rho_n,
                         float& V,     float& U,
                         float& T_i,   float& T_n) {
        const auto sz = arma::size(ns, num_of_eq);
        rho_i = xn(arma::sub2ind(sz, i, cons::RHO_I));
        rho_n = xn(arma::sub2ind(sz, i, cons::RHO_N));
        const float momI = xn(arma::sub2ind(sz, i, cons::MOM_I));
        const float momN = xn(arma::sub2ind(sz, i, cons::MOM_N));
        const float E_i  = xn(arma::sub2ind(sz, i, cons::E_I));
        const float E_n  = xn(arma::sub2ind(sz, i, cons::E_N));
        V = momI / rho_i;
        U = momN / rho_n;
        const float phi_g_cell = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        // cons2prim convention (state.cpp): p_i = 2 n_i k_b T_i with electron
        // quasi-neutrality folded in, p_n = n_n k_b T_n.
        const float p_i = 2.0f/3.0f * E_i - 1.0f/3.0f * rho_i * V * V
                         - 2.0f/3.0f * rho_i * phi_g_cell;
        const float p_n = 2.0f/3.0f * E_n - 1.0f/3.0f * rho_n * U * U
                         - 2.0f/3.0f * rho_n * phi_g_cell;
        const float n_i = rho_i / m_i;
        const float n_n = rho_n / m_n;
        T_i = p_i / (2.0f * n_i * k_b);
        T_n = p_n /        (n_n * k_b);
    };

    auto pack_ghost = [&](Vec& ob,
                          float rho_i, float rho_n,
                          float V_g,   float U_g,
                          float T_i,   float T_n,
                          float phi_g_g) {
        const float n_i = rho_i / m_i;
        const float n_n = rho_n / m_n;
        ob(cons::RHO_I) = rho_i;
        ob(cons::RHO_N) = rho_n;
        ob(cons::MOM_I) = rho_i * V_g;
        ob(cons::MOM_N) = rho_n * U_g;
        ob(cons::E_I)   = 1.5f * k_b * n_i * 2.0f * T_i
                        + 0.5f * rho_i * V_g * V_g + rho_i * phi_g_g;
        ob(cons::E_N)   = 1.5f * k_b * n_n * T_n
                        + 0.5f * rho_n * U_g * U_g + rho_n * phi_g_g;
        ob(cons::E_E)   = 1.5f * k_b * n_i * T_i;   // electron internal energy (T_e = T_i)
    };

    // --- Inner reflecting wall (mirror cells 0 and 1) ----------------------
    {
        float rho_i, rho_n, V, U, T_i, T_n;
        const float phi_g_in = grid.phi_g_imh(0);
        cell_prim(0, rho_i, rho_n, V, U, T_i, T_n);
        pack_ghost(grid.inner_boundary0_i, rho_i, rho_n, -V, -U, T_i, T_n, phi_g_in);
        if (ns >= 2) {
            cell_prim(1, rho_i, rho_n, V, U, T_i, T_n);
            pack_ghost(grid.inner_boundary1_i, rho_i, rho_n, -V, -U, T_i, T_n, phi_g_in);
        } else {
            grid.inner_boundary1_i = grid.inner_boundary0_i;
        }
    }

    // --- Outer face -------------------------------------------------------
    // Default (open line / coronal top): pure Neumann outflow — ρ, T, V, U all
    // zero-gradient at both outer ghosts. The boundary has no impedance —
    // outgoing waves see a transparent face and no externally-specified state on
    // the incoming characteristic. Stable only when the IC is near HSE /
    // steady-state; transient runs that need boundary damping may go unstable.
    //
    // When grid.outer_reflecting (closed-loop apex): mirror cells ns-1, ns-2 with
    // V,U → −V,−U so the face is a zero-flux symmetry plane — the evaporated
    // plasma is confined and fills the loop. g_∥→0 at the apex, so the gravity
    // potential is flat there and reusing phi_g_iph(ns-1) for both ghosts is exact
    // to leading order.
    {
        float rho_i, rho_n, V, U, T_i, T_n;
        const float phi_g_out = grid.phi_g_iph(ns - 1);
        cell_prim(ns - 1, rho_i, rho_n, V, U, T_i, T_n);
        if (grid.outer_reflecting) {
            pack_ghost(grid.outer_boundary0_i, rho_i, rho_n, -V, -U, T_i, T_n, phi_g_out);
            if (ns >= 2) {
                cell_prim(ns - 2, rho_i, rho_n, V, U, T_i, T_n);
                pack_ghost(grid.outer_boundary1_i, rho_i, rho_n, -V, -U, T_i, T_n, phi_g_out);
            } else {
                grid.outer_boundary1_i = grid.outer_boundary0_i;
            }
        } else {
            pack_ghost(grid.outer_boundary0_i, rho_i, rho_n, V, U, T_i, T_n, phi_g_out);
            pack_ghost(grid.outer_boundary1_i, rho_i, rho_n, V, U, T_i, T_n, phi_g_out);
        }
    }

    grid.broadcast();
}

} // namespace chromosphere
