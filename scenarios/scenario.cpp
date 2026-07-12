#include "scenario.hpp"
#include "analytic_canopy.hpp"
#include "model_c7.hpp"
#include "model_column.hpp"
#include "model_flare.hpp"
#include "pfss_field_line.hpp"

#include <stdexcept>
#include <cstdlib>
#include <string>

namespace {
// Set an environment variable only if the user has not already set it — used to
// apply a scenario "preset" (default env knobs) that the user can still override.
void set_env_default(const char* key, const char* value) {
    setenv(key, value, /*overwrite=*/0);
}

// Grid resolution for the unified column scenario: a resolved corona (ISO_CORONA)
// needs more cells than the chromosphere-only default. ISO_NS overrides.
arma::uword column_peek_ns() {
    auto set = [](const char* k) {
        const char* e = std::getenv(k);
        return e && std::string(e) != "0" && !std::string(e).empty();
    };
    arma::uword ns = set("ISO_CORONA") ? 1000 : 600;
    if (const char* e = std::getenv("ISO_NS")) { try { ns = std::stoul(e); } catch (...) {} }
    return ns;
}
} // namespace

namespace chromosphere {

Scenario make_scenario(const std::string& name, const std::string& data_path) {
    Scenario sc;
    sc.name = name;

    if (name == "model_c7") {
        sc.peek_ns   = []() -> arma::uword { return 100; };
        // "New explanation" upper BC (docs/gentle_evaporation_downflow.md):
        // well-balanced reconstruction + hydrostatic ghost pressure + imposed
        // TR temperature jump (downward q(T) via the Stage-D Dirichlet ghost-T)
        // + EOS ghost density. tr_jump_bc=true is exclusive to this bare scenario
        // — model_flare / analytic_canopy / model_gentle keep the RTV reservoir.
        sc.ic        = [](Grid& g) { return model_c7_ic(g, /*extended=*/false, /*tr_jump_bc=*/true); };
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
    // The unified field-aligned chromosphere→corona column scenario (the DEFAULT).
    // "model_isentropic" is a backward-compatibility alias for the old name.
    // "model_gentle" is a preset: the documented stable full-physics resolved-corona
    //   gentle conduction-driven evaporation configuration (iso_corona_full) —
    //   resolved corona + radiative sink + two-fluid + ionization, conduction on. The
    //   base is raised to h = 1003 km (model_c7's validated floor): the Stage-E n²
    //   channel counts and the radiative loss overflow float32 at photospheric density
    //   (n ~ 1e23 ⇒ n² ≫ FLT_MAX). All are env defaults the user can still override.
    if (name == "model_column" || name == "model_isentropic" || name == "model_gentle") {
        if (name == "model_gentle") {
            set_env_default("ISO_CORONA",     "1");
            set_env_default("ISO_H_BASE",     "1003");
            set_env_default("ISO_HEAT_FLUX",  "1");
            set_env_default("ISO_COOLING",    "1");
            set_env_default("ISO_IONIZATION", "1");
            set_env_default("ISO_TWO_FLUID",  "1");
            set_env_default("ISO_NS",         "600");
        }
        sc.peek_ns   = []() { return column_peek_ns(); };
        sc.ic        = [](Grid& g) { return model_column_ic(g); };
        sc.update_bc = [](Grid& g, const Vec& xn) { model_column_update_bc(g, xn); };
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
        const float p_i = grid.gm1() * E_i - grid.half_gm1() * rho_i * V * V
                         - grid.gm1() * rho_i * phi_g_cell;
        const float p_n = grid.gm1() * E_n - grid.half_gm1() * rho_n * U * U
                         - grid.gm1() * rho_n * phi_g_cell;
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
        ob(cons::E_I)   = grid.inv_gm1() * k_b * n_i * 2.0f * T_i
                        + 0.5f * rho_i * V_g * V_g + rho_i * phi_g_g;
        ob(cons::E_N)   = grid.inv_gm1() * k_b * n_n * T_n
                        + 0.5f * rho_n * U_g * U_g + rho_n * phi_g_g;
        ob(cons::E_E)   = grid.inv_gm1() * k_b * n_i * T_i;   // electron internal energy (T_e = T_i)
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
