/*!
 * @file scenarios/model_column.cpp
 * @brief The field-aligned chromosphere-to-corona column: the `model_column`
 *        RELEASE scenario.
 * @ingroup scenarios
 *
 * The IC/BC of the release Saha/Gamma1 column. It loads the
 * production `Gamma1` table unconditionally, builds the Saha-HSE initial column
 * from the C7 temperature profile, sets the validated release numerics
 * (log-space MUSCL with the MC3/Koren limiter, the SWMF exact-Riemann Godunov
 * flux, the pressure-based reconstruction) and the 22,000 K external conductive
 * reservoir at the outer face — and rejects every historical two-fluid knob.
 * The historical fixed-gamma two-fluid column that used to share this file is
 * the legacy `model_gentle` scenario (legacy/scenarios/model_gentle.cpp); the
 * two-temperature twin that also shared it until 2026-11-16 is the legacy
 * `model_column_2t` (legacy/scenarios/model_column_2t.cpp). See
 * scenarios/model_column.hpp for the full environment-variable list.
 */
#include "model_column.hpp"
#include "single_fluid/mixture.hpp"
#include "c7_profile.hpp"  // c7_full_profile + c7_full_temperature_pchip (C7 table)
#include "mesh.hpp"       // shared static local-refinement mesh builder

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

//! @cond INTERNAL
extern char** environ;   // POSIX: the process environment
//! @endcond

namespace chromosphere {

// ---------------------------------------------------------------------------
// Module state captured by model_column_ic and reused by the boundary update.
// The interior is single-fluid neutral-dominated hydrogen with a realistic
// frozen ionization profile (neutral chromosphere, ionized corona).
// ---------------------------------------------------------------------------
namespace {
// ---------------------------------------------------------------------------
// Inner boundary — the LOWER TRUNCATION of the modeled domain at h_base
// (1600 km by default). This is NOT the photosphere: it is a mid-chromospheric
// cut chosen because the Saha/LTE equilibrium closure this solver uses is only
// intended to be trusted from about that height upward. The real atmosphere
// continues below it and is simply outside the model. Consequently the two
// ghost cells are a NUMERICAL closure for the MUSCL stencil, not a claim that
// the LTE model extends downward, and the physical content of the boundary is
// exactly one thing: the lower chromosphere below the cut behaves, on the
// timescales resolved here, as a quasi-static stratified reservoir.
// ---------------------------------------------------------------------------
Real kInnerRhoI   = 0.0f;      // base ionized-carrier density (C7 n_e m_i)
Real kInnerRhoN   = 0.0f;      // base neutral-carrier density (C7 n_HI m_n)
Real kInnerTRef     = 0.0f;     // reservoir temperature at the truncation
// RELEASE lower reservoir, captured once at the IC base. Unlike the historical
// two-fluid first-order ladder (legacy model_gentle), each rung is the
// trapezoidal, EOS-closed hydrostatic step
// p_{k+1} = p_k + Δs·½(ρ_k + ρ_{k+1})·g with ρ_{k+1} = ρ_EOS(p_{k+1}, T₀),
// iterated to a fixed point — the SAME construction the outer face already used
// (hse_rung), so the two boundaries are now discretely consistent. The first-
// order form this replaces carried an O((Δs/H)²) pressure error, i.e. an O(Δs/H)
// error in the hydrostatic force balance of the lowest interior cells.
Real kInnerPRes0   = 0.0f;     // reservoir pressure at ghost 0 (1 Δs below cell 0)
Real kInnerPRes1   = 0.0f;     // reservoir pressure at ghost 1 (2 Δs below)
Real kInnerRhoRes0 = 0.0f;     // EOS density at (kInnerPRes0, kInnerTRef)
Real kInnerRhoRes1 = 0.0f;
Real kTtopRef     = 0.0f;      // FIXED outer ghost temperature [K]: CHROMO_BC_T_TOP_K, else the IC top cell
// FIXED reference outer ghost pressure [Pa] — the coronal reservoir back-pressure,
// captured one-sided-hydrostatically from the IC top cell. Subsonic outflow admits
// exactly ONE incoming characteristic,
// so exactly one condition may come from outside: this back-pressure. Anchoring it to
// a FIXED reference rather than to a live interior cell is what makes the outer face
// well posed — see the outer-face block for why both live-cell anchors fail.
Real kOuterPRef    = 0.0f;
bool  kOuterPRefSet = false;    // captured lazily on the first update_bc call
bool  kHeatFluxOn = false;     // conduction on ⇒ open (Mach-capped) top
Real kVcapMach    = 0.1f;      // outer outflow cap as a fraction of c_s
// CHROMO_BC_OUTER_P=hydrostatic: the outer ghost pressure stops being the FIXED reservoir
// back-pressure kOuterPRef and becomes the same EOS-closed hydrostatic ladder
// (hse_rung) ANCHORED ON THE LIVE TOP CELL — a hydrostatically consistent
// closure that imposes no external pressure datum. Subsonic outflow admits exactly
// ONE incoming characteristic, so this swaps WHICH variable carries it, from the
// back-pressure to the fixed ghost temperature. Diagnostic; measured in
// docs/studies/boundaries/outer_bc_pressure_temperature_swap_recap.md (it held
// 4000 s but suppressed the evaporation signal ~19x). The literal dp/ds = 0
// variant (mode 1) destroyed the model there and was removed.
bool  kOuterPHydrostatic = false;
// The release hydrostatic reconstruction (Grid::hydrostatic_reconstruct)
// balances the discrete relation
// ln p_{i+1} - ln p_i = -(g/2)(ds_i rho_i/p_i + ds_{i+1} rho_{i+1}/p_{i+1}), so the
// IC and both ghost ladders below are built in that same form.

double gamma_density_from_pressure(const Grid& grid, double pressure, double temperature) {
    if (!(pressure > 0.0) || !(temperature > 0.0))
        throw std::domain_error("gamma ghost/HSE pressure and temperature must be positive");
    const double lo = grid.eos_gamma_table.min_n_h();
    const double hi = grid.eos_gamma_table.max_n_h();
    auto pressure_at = [&](double n_h) {
        const double x = saha_ionization_fraction_n_h(n_h, temperature);
        return (1.0+x)*n_h*eos_constants::k_b*temperature;
    };
    // EOS density-domain guard, unchanged: p(n_H) is monotone in n_H at fixed T,
    // so bracketing the request between the two axis endpoints is exactly the
    // in-table test the former bisection performed before iterating.
    if (pressure < pressure_at(lo) || pressure > pressure_at(hi))
        throw std::out_of_range("gamma ghost/HSE pressure implies density outside EOS table");
    // The 100 geometric bisections this replaced (200 Saha evaluations per call)
    // solved the same closure the algebraic inverse below solves in closed form.
    return equilibrium_density_from_pressure(pressure, temperature);
}

// Gravity magnitude the release solver actually applies in cell i [m/s^2]: its
// source differences phi_g, which carries a float-literal km->m factor, so it
// differs from grid.g by ~5e-8 relative. The hydrostatic reconstruction must
// balance against this value, not grid.g.
double solver_g(const Grid& grid, arma::uword i) {
    return static_cast<double>((grid.phi_g_iph(i) - grid.phi_g_imh(i))/grid.ds_i(i));
}

// One rung of an EOS-consistent hydrostatic ladder taken DOWNWARD by `ds`
// metres from an anchor (p_anchor, rho_anchor) at fixed temperature T_gh, in the
// discrete relation the hydrostatic reconstruction balances:
//     ln p = ln p_anchor + ½Δs·g·(ρ_anchor/p_anchor + ρ/p),     ρ = ρ_EOS(p, T_gh),
// iterated to a fixed point.
void gamma_hse_rung_down(const Grid& grid, double p_anchor, double rho_anchor,
                         double T_gh, double ds, double& p_out, double& rho_out) {
    double rho_gh = rho_anchor;
    double p_gh = p_anchor;
    for (int it = 0; it < 40; ++it) {
        const double p =
            p_anchor*std::exp(0.5*ds*solver_g(grid, 0)*(rho_anchor/p_anchor + rho_gh/p_gh));
        p_gh = p;
        const double rho_next = gamma_density_from_pressure(grid, p, T_gh);
        const bool done = std::abs(rho_next - rho_gh) <= 1.0e-14*rho_next;
        p_out = p;
        rho_gh = rho_next;
        if (done) break;
    }
    rho_out = rho_gh;
}

// Per-cell (T, n_i, n_n) at height h_km from the real C7 atmosphere via
// c7_full_profile (Avrett & Loeser Table 26): n_i = n_e (ionized), n_n = n_HI.
void c7_cell(Real h_km, Real& T, Real& n_i, Real& n_n) {
    Real n_e, n_HI;
    c7_full_profile(h_km, T, n_e, n_HI);
    n_i = n_e;
    n_n = n_HI;
}

// Cell-average of the prescribed C7/PCHIP temperature for the finite-volume
// Gamma/Saha IC. Two-point Gauss-Legendre quadrature is enough for the smooth
// PCHIP segments and removes the O(ds) top-cell centre-sampling shift under
// refinement without introducing a general quadrature framework.
double c7_temperature_cell_average(double h_lo_km, double h_hi_km) {
    constexpr double inv_sqrt3 = 0.57735026918962576451;
    const double mid = 0.5*(h_lo_km+h_hi_km);
    const double half = 0.5*(h_hi_km-h_lo_km);
    const double h0 = mid-half*inv_sqrt3;
    const double h1 = mid+half*inv_sqrt3;
    return 0.5*(c7_full_temperature_pchip(h0)+c7_full_temperature_pchip(h1));
}

// The release reads only CHROMO_* names. Until 2026-11-16 its knobs carried
// the ISO_ prefix of a retired experimental model (and GRID_REFINE_* /
// TT_OUTER_TI / GAMMA_TABLE); the frozen legacy scenarios still read those. A stale old name
// would otherwise be ignored silently and the run would quietly use the release
// default instead, so any variable in those families stops a release run with
// the name to use instead.
void reject_retired_env() {
    struct Renamed { const char* old_name; const char* new_name; };
    static const Renamed kRenamed[] = {
        {"ISO_H_BASE", "CHROMO_DOMAIN_BASE_KM"},
        {"ISO_DH", "CHROMO_DOMAIN_SPAN_KM"},
        {"ISO_NS", "CHROMO_MESH_NS"},
        {"ISO_REFINE_FACTOR", "CHROMO_MESH_REFINE_FACTOR"},
        {"ISO_REFINE_PROFILE", "CHROMO_MESH_REFINE_PROFILE"},
        {"ISO_REFINE_S_LO_KM", "CHROMO_MESH_REFINE_START_KM"},
        {"ISO_REFINE_S_HI_KM", "CHROMO_MESH_REFINE_END_KM"},
        {"ISO_REFINE_TRANSITION_KM", "CHROMO_MESH_REFINE_GRADE_KM"},
        {"GRID_REFINE_FACTOR", "CHROMO_MESH_REFINE_FACTOR"},
        {"GRID_REFINE_PROFILE", "CHROMO_MESH_REFINE_PROFILE"},
        {"GRID_REFINE_S_LO_KM", "CHROMO_MESH_REFINE_START_KM"},
        {"GRID_REFINE_S_HI_KM", "CHROMO_MESH_REFINE_END_KM"},
        {"GRID_REFINE_TRANSITION_KM", "CHROMO_MESH_REFINE_GRADE_KM"},
        {"ISO_HEAT_FLUX", "CHROMO_CONDUCTION"},
        {"ISO_P_TOP", "CHROMO_IC_P_TOP_PA"},
        {"ISO_T_TOP", "CHROMO_BC_T_TOP_K"},
        {"ISO_VCAP", "CHROMO_BC_VCAP_MACH"},
        {"ISO_OUTER_P_NEUMANN", "CHROMO_BC_OUTER_P (reservoir | hydrostatic)"},
        {"ISO_COND_KIRCHHOFF", "CHROMO_COND_FACE=kirchhoff"},
        {"ISO_COND_ARITH", "CHROMO_COND_FACE=arithmetic"},
        {"ISO_RIEMANN", "CHROMO_NUM_RIEMANN (godunov | roe)"},
        {"GAMMA_TABLE", "CHROMO_EOS_TABLE"}};
    // Research configuration of the historical two-fluid column.
    static const char* const kTwoFluid[] = {
        "ISO_GAMMA", "ISO_TWO_FLUID", "ISO_IONIZATION", "ISO_COOLING",
        "ISO_TRAC", "ISO_CORONA", "ISO_CHEAT", "ISO_QFLUX", "ISO_TBOOST",
        "ISO_NUMERICAL_DIFFUSIVITY_MULT"};
    for (char** entry = environ; entry && *entry; ++entry) {
        const std::string pair(*entry);
        const std::string name = pair.substr(0, pair.find('='));
        const bool retired_family = name.rfind("ISO_", 0) == 0
            || name.rfind("GRID_REFINE_", 0) == 0 || name == "TT_OUTER_TI"
            || name == "GAMMA_TABLE";
        if (!retired_family) continue;
        for (const Renamed& r : kRenamed)
            if (name == r.old_name)
                throw std::runtime_error(
                    "model_column: " + name + " was renamed to " + r.new_name
                    + "; unset it (the release reads only CHROMO_* names)");
        for (const char* knob : kTwoFluid)
            if (name == knob)
                throw std::runtime_error(
                    "model_column: " + name + " is historical two-fluid research "
                    "configuration and has no meaning for the single-fluid release "
                    "solver; run the legacy model_gentle scenario for that path "
                    "(build with -DCHROMO_ENABLE_LEGACY=ON)");
        throw std::runtime_error(
            "model_column: " + name + " is not a release setting (it was removed, "
            "or it belongs to a legacy scenario); unset it. See the configuration "
            "reference for the CHROMO_* names.");
    }
}
} // namespace

Vec model_column_ic(Grid& grid) {
    auto env_f = [](const char* key, Real fallback) -> Real {
        if (const char* e = std::getenv(key)) {
            try { return std::stof(std::string(e)); } catch (...) {}
        }
        return fallback;
    };

    reject_retired_env();
    if (const char* gamma_path = std::getenv("CHROMO_EOS_TABLE")) {
        if (grid.eos_gamma_table.empty()) grid.eos_gamma_table = EosGammaTable::load(gamma_path);
    }
    // A loaded Gamma1 table is what selects the Saha/Gamma1 solver. `model_column`
    // always supplies one, so an empty table here means the caller asked for the
    // historical fixed-gamma two-fluid column, which is the separate legacy
    // `model_gentle` scenario.
    if (grid.eos_gamma_table.empty())
        throw std::runtime_error(
            "model_column_ic: no Gamma1 table is loaded. model_column always "
            "runs the Saha/Gamma1 release solver; the historical "
            "fixed-gamma two-fluid column is the legacy model_gentle scenario "
            "(build with -DCHROMO_ENABLE_LEGACY=ON)");

    // The historical two-fluid research knobs (a fixed adiabatic index, a neutral
    // fluid, finite-rate ionization, a radiative sink, TRAC, coronal heating, an
    // imposed coronal flux, an IC superheat, artificial conduction) are rejected
    // by reject_retired_env() above, which is what keeps `model_column` a single
    // unambiguous production identity.

    // --- IC mode ----------------------------------------------------------
    // The IC is ALWAYS the real Model C7 profile with the density re-integrated
    // hydrostatically (clean V≈0 start; the analytic-isentrope toy IC is retired).

    // --- parameters -------------------------------------------------------
    // Lower truncation height of the modeled domain, in km measured from the
    // photosphere (that is the HEIGHT GAUGE only -- the domain itself starts
    // here, at 1600 km, well above it). See the inner-boundary block.
    const Real h_base     = env_f("CHROMO_DOMAIN_BASE_KM", 0.0f);
    Real       DH_km      = env_f("CHROMO_DOMAIN_SPAN_KM", 1303.0f); // km
    // Base (T, n_H) from the C7 profile at h_base.
    Real T_base, ni_b, nn_b;
    c7_cell(h_base, T_base, ni_b, nn_b);
    T_base = static_cast<Real>(c7_full_temperature_pchip(h_base));
    const bool  heat_flux_on = (env_f("CHROMO_CONDUCTION", 0.0f) != 0.0f);
    kHeatFluxOn        = heat_flux_on;
    kVcapMach          = env_f("CHROMO_BC_VCAP_MACH", 0.1f);
    // Outer-face closure variant: which ONE condition the subsonic-outflow face
    // imposes. Release (0) = fixed reservoir back-pressure + fixed ghost T;
    // 2 = the swap, a hydrostatic pressure anchored on the live top cell + the
    // same fixed ghost T.
    kOuterPHydrostatic = false;
    if (const char* e = std::getenv("CHROMO_BC_OUTER_P")) {
        const std::string mode(e);
        if (mode == "hydrostatic") kOuterPHydrostatic = true;
        else if (mode != "reservoir")
            throw std::runtime_error(
                "model_column_ic: CHROMO_BC_OUTER_P must be reservoir (fixed "
                "back-pressure, the release) or hydrostatic (dp/ds = -rho g, "
                "anchored on the live top cell)");
    }

    // The release is a single-fluid common-temperature equilibrium mixture.
    grid.single_fluid = true;
    grid.enable_ionization = false;
    grid.enable_Te = false;

    // --- geometry: straight field line, gravity on -----------------------
    // Static local refinement (scenarios/mesh.hpp). peek_ns() sized the grid to the
    // coarse-equivalent count (CHROMO_MESH_NS); CHROMO_MESH_REFINE_* ADD cells inside
    // the refined band and grade smoothly to the coarse spacing outside it.
    // CHROMO_MESH_REFINE_PROFILE picks
    // which end the band is anchored to: `lower` (default, historical: fine [s_lo,s_hi]
    // then grade back up to coarse) or `outer` (fine from s_lo to the domain TOP, with
    // the grade in [s_lo−τ, s_lo]) — the latter targets the upper TR / outer boundary.
    // The conduction solver uses these actual mesh metrics to form the face-local
    // numerical diffusivity; ds_m remains useful only for uniform-grid geometry.
    // Uniform mesh reproduced byte-for-byte.
    const arma::uword ns_coarse = grid.ns;
    const Real ds_m = (DH_km * 1000.0f) / static_cast<Real>(ns_coarse);   // coarse Δs [m]
    const RefineParams rp = refine_params_from_env(kReleaseRefineEnv);

    std::vector<double> face_m;   // canonical faces [m] from the inner face (empty ⇒ uniform)
    if (rp.enabled()) {
        face_m = build_static_mesh_faces(static_cast<double>(DH_km) * 1000.0, ns_coarse, rp);
        grid.resize(static_cast<arma::uword>(face_m.size()) - 1);   // zeros the ns-sized arrays
    }

    // Set the field-line geometry AFTER any resize (resize reallocates/zeros these).
    grid.B_i.ones();
    grid.B_imh.ones();
    grid.B_iph.ones();
    grid.dinvB_ds_i.zeros();

    Vec h_F(grid.ns + 1);                                   // face heights [km]
    if (!face_m.empty()) {
        for (arma::uword k = 0; k <= grid.ns; ++k)
            h_F(k) = h_base + static_cast<Real>(face_m[k]) * 1.0e-3f;   // m -> km
        for (arma::uword i = 0; i < grid.ns; ++i) {
            grid.ds_i(i)      = static_cast<Real>(face_m[i + 1] - face_m[i]);
            grid.phi_g_imh(i) = grid.g * (h_F(i)     - h_base) * 1000.0f;   // J/kg
            grid.phi_g_iph(i) = grid.g * (h_F(i + 1) - h_base) * 1000.0f;
        }
    } else {
        for (arma::uword k = 0; k <= grid.ns; ++k)
            h_F(k) = h_base + DH_km * static_cast<Real>(k) / static_cast<Real>(grid.ns);
        for (arma::uword i = 0; i < grid.ns; ++i) {
            grid.ds_i(i)      = ds_m;
            grid.phi_g_imh(i) = grid.g * (h_F(i)     - h_base) * 1000.0f;   // J/kg
            grid.phi_g_iph(i) = grid.g * (h_F(i + 1) - h_base) * 1000.0f;
        }
    }

    // --- finite-volume profile (T, n_i, n_n) -------------------------------
    // Real C7 T(h) and ionization FRACTION with the total density re-integrated
    // hydrostatically from the base (p(h)=p_base·exp(−∫ g/(R_s T) dh), R_s local),
    // so V≈0 is a clean discrete fixed point and any flow is heat-flux driven.
    // The production Gamma/Saha path stores a cell-average C7/PCHIP temperature
    // rather than a centre point sample, so refinement represents the same
    // continuum upper thermal profile instead of shifting the top-cell datum.
    Vec  T_c(grid.ns), ni_c(grid.ns), nn_c(grid.ns);
    // One hydrostatic integration from the base pressure exp(ln_p_base) with base
    // ionization x_base: fills T_c/ni_c/nn_c and returns ln p of the top cell.
    auto integrate_column = [&](double ln_p_base, double x_base) -> double {
        double ln_p = ln_p_base, prev_h = h_base, prev_ds = 0.0;
        double prev_invH = grid.g*grid.m_i/((1.0+x_base)*grid.k_b*T_base);
        for (arma::uword i = 0; i < grid.ns; ++i) {
            const double h_c = 0.5 * (h_F(i) + h_F(i + 1));        // km
            Real T, n_i, n_n;
            c7_cell(static_cast<Real>(h_c), T, n_i, n_n);
            T = static_cast<Real>(c7_temperature_cell_average(h_F(i), h_F(i+1)));
            const double dh_m = (h_c-prev_h)*1000.0;
            double x = n_i/(n_i+n_n);
            {
                const double previous_ln_p = ln_p;
                double trial = ln_p-prev_invH*dh_m;
                bool converged = false;
                for (int it = 0; it < 60; ++it) {
                    const double rho = gamma_density_from_pressure(grid, std::exp(trial), T);
                    x = saha_ionization_fraction_n_h(rho/eos_constants::m_h, T);
                    const double invH = grid.g*grid.m_i/((1.0+x)*grid.k_b*T);
                    // Cell 0 is anchored half a cell above the base.
                    const double next = (i > 0)
                        ? previous_ln_p-0.5/grid.g*(prev_invH*solver_g(grid, i-1)*prev_ds
                                                    +invH*solver_g(grid, i)*grid.ds_i(i))
                        : previous_ln_p-0.5*(prev_invH+invH)*dh_m;
                    if (std::abs(next-trial) < 2.0e-13) {
                        trial = next;
                        converged = true;
                        break;
                    }
                    trial = next;
                }
                if (!converged) {
                    throw std::runtime_error(
                        "gamma Saha-HSE iteration did not converge at cell "
                        + std::to_string(i)+", h_km="+std::to_string(h_c)
                        +", T="+std::to_string(T)+", p="+std::to_string(std::exp(trial)));
                }
                ln_p = trial;
            }
            const double n_tot = std::exp(ln_p)/((1.0+x)*grid.k_b*T);
            x = saha_ionization_fraction_n_h(n_tot, T);
            n_i = x * n_tot;
            n_n = (1.0-x)*n_tot;
            prev_h = h_c;
            prev_ds = grid.ds_i(i);
            prev_invH = grid.g*grid.m_i/((1.0+x)*grid.k_b*T);
            T_c(i)  = T;
            ni_c(i) = n_i;
            nn_c(i) = n_n;
        }
        return ln_p;
    };
    {
        const double nbase = static_cast<double>(ni_b)+nn_b;
        const double xb = saha_ionization_fraction_n_h(nbase, T_base);
        const double ln_p_c7 = std::log((1.0+xb)*nbase*grid.k_b*T_base);
        const double ln_p_top_c7 = integrate_column(ln_p_c7, xb);

        // CHROMO_IC_P_TOP_PA > 0 [Pa]: rebuild the column with the SAME C7 T(h) but the base
        // pressure chosen so the top cell sits at this pressure. The whole column
        // moves consistently (both reservoirs are captured from the rebuilt IC), so
        // the pressure is a controlled parameter rather than a one-cell mismatch.
        // Saha makes p_top a nonlinear function of p_base (a denser column is less
        // ionized and has a shorter scale height), hence the secant iteration.
        if (const char* e = std::getenv("CHROMO_IC_P_TOP_PA")) {
            double target = std::numeric_limits<double>::quiet_NaN();
            try { target = std::stod(e); } catch (...) {}
            if (!(std::isfinite(target) && target > 0.0))
                throw std::runtime_error(std::string("model_column_ic: CHROMO_IC_P_TOP_PA must be a "
                    "positive pressure in Pa, got '") + e + "'");
            const double ln_target = std::log(target);
            auto base_x = [&](double ln_pb) {
                const double rho_b = gamma_density_from_pressure(grid, std::exp(ln_pb), T_base);
                return saha_ionization_fraction_n_h(rho_b/eos_constants::m_h, T_base);
            };
            double a0 = ln_p_c7, f0 = ln_p_top_c7-ln_target;
            double a1 = a0-f0;
            double f1 = integrate_column(a1, base_x(a1))-ln_target;
            int it = 0;
            for (; it < 40 && std::abs(f1) > 1.0e-12; ++it) {
                const double a2 = a1-f1*(a1-a0)/(f1-f0);
                a0 = a1; f0 = f1; a1 = a2;
                f1 = integrate_column(a1, base_x(a1))-ln_target;
            }
            if (std::abs(f1) > 1.0e-12)
                throw std::runtime_error("model_column_ic: CHROMO_IC_P_TOP_PA secant iteration "
                    "did not converge (target " + std::to_string(target) + " Pa)");
            std::cerr << "[model_column] CHROMO_IC_P_TOP_PA: top-cell p=" << std::exp(f1+ln_target)
                      << " Pa (target " << target << ", C7 " << std::exp(ln_p_top_c7)
                      << ")  base p=" << std::exp(a1) << " Pa (C7 " << std::exp(ln_p_c7)
                      << ")  secant iterations=" << it+1 << std::endl;
        }
    }

    Vec xn = arma::zeros<Vec>(grid.n_mixture_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const Real T     = T_c(i);
        const Real n_i   = ni_c(i);
        const Real n_n   = nn_c(i);
        const Real phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        const Real rho_i = n_i * grid.m_i;
        const Real rho_n = n_n * grid.m_n;
        mixture_pack_cell(grid, xn, i, static_cast<double>(rho_i)+rho_n,
                          0.0, T, phi_g);
    }

    // Outer ghost temperature, FIXED for the whole run: the unresolved TR/corona
    // above sits at ~constant T, so the ghost is a fixed hot reservoir and the
    // downward conductive flux q ∝ (T_ghost − T_top) SELF-LIMITS as the
    // chromosphere heats. CHROMO_BC_T_TOP_K > 0 sets it absolutely [K] (release 22000,
    // C7's TR-base T); otherwise it is the IC top-cell temperature.
    kTtopRef = T_c(grid.ns - 1);
    const Real t_top_abs = env_f("CHROMO_BC_T_TOP_K", 0.0f);
    if (t_top_abs > 0.0f) kTtopRef = t_top_abs;

    // --- lower truncation reservoir (h_base = 1600 km, NOT the photosphere) ---
    // The ghosts sit one and two Δs below cell 0 and stand in for the unmodeled
    // lower chromosphere, which on these timescales is a quasi-static stratified
    // reservoir. The RELEASE ladder is trapezoidal and EOS-closed at the
    // truncation temperature T₀ (the legacy two-fluid column, model_gentle, keeps
    // the historical first-order form p_0 + kρ_0gΔs with ρ ∝ p).
    {
        const Real T0    = T_c(0);
        kInnerTRef = T0;
        const Real n_i0  = ni_c(0);
        const Real n_n0  = nn_c(0);
        kInnerRhoI = n_i0 * grid.m_i;
        kInnerRhoN = n_n0 * grid.m_n;
        // Center-to-ghost distance = the BASE cell's own width (mirrored ghost). On a
        // refined mesh cell 0 is a fine cell, so this is the fine Δs, NOT the coarse ds_m.
        const Real ds_base = grid.ds_i(0);
        const double rho0 = kInnerRhoI+kInnerRhoN;
        const double p0 = (1.0+saha_ionization_fraction_n_h(
            rho0/eos_constants::m_h, T0))*rho0/eos_constants::m_h*eos_constants::k_b*T0;
        double pg0, rg0, pg1, rg1;
        gamma_hse_rung_down(grid, p0,  rho0, T0, ds_base, pg0, rg0);
        gamma_hse_rung_down(grid, pg0, rg0,  T0, ds_base, pg1, rg1);
        kInnerPRes0 = pg0; kInnerRhoRes0 = rg0;
        kInnerPRes1 = pg1; kInnerRhoRes1 = rg1;
    }

    // --- numerical scheme (always on, no toggles) ---
    // MC3/Koren limiter (β=2) on the (ln rho, u, ln p) reconstruction. The
    // release reconstruction is log-space by construction and decodes each cell
    // with its own φ_g, so the two-fluid-only Grid::log_reconstruct and
    // Grid::well_balanced are not set. The scheme is reference-free and exactly
    // balanced at rest through the hydrostatic reconstruction (see the eq_wb note
    // at the end of this IC).
    grid.out_base_km              = h_base;   // label output heights from the true base
    // DIAGNOSTIC-ONLY conduction face-conductivity rule (TR-pulse study).
    // harmonic = the release rule (distance-weighted harmonic mean; arithmetic on
    // a uniform mesh); kirchhoff = Kirchhoff secant for the Spitzer part;
    // arithmetic = plain arithmetic mean of the total conductivity.
    grid.conduction_kirchhoff_face = false;
    grid.conduction_arithmetic_face = false;
    if (const char* e = std::getenv("CHROMO_COND_FACE")) {
        const std::string rule(e);
        if (rule == "kirchhoff") grid.conduction_kirchhoff_face = true;
        else if (rule == "arithmetic") grid.conduction_arithmetic_face = true;
        else if (rule != "harmonic")
            throw std::runtime_error(
                "model_column_ic: CHROMO_COND_FACE must be harmonic (the release), "
                "kirchhoff or arithmetic");
    }
    grid.mc3_limiter              = true;
    grid.limiter_beta             = 2.0f;
    // --- release numerical flux and MUSCL primitive set (Gamma/Saha path) ---
    // The released Gamma/Saha model_column advances the corrector with the
    // SWMF-style exact-Riemann Godunov flux — the local Riemann problem solved at
    // the frozen-composition index gamma = 5/3 with the ionization energy carried
    // as SWMF's passive offset E0 — and limits (log rho, V, log p) with the ln p
    // slot hydrostatic, recovering the face temperature by inverting the same
    // authoritative Saha closure.
    grid.swmf_godunov_flux       = true;
    grid.roe_characteristic_flux = false;
    grid.hydrostatic_reconstruct = true;
    // Reference override. CHROMO_NUM_RIEMANN=roe selects the mixture Roe
    // characteristic flux — the equilibrium (Gamma1) linearization that was the
    // previous release flux, kept as the controlled comparison solver. Never
    // needed for a production run. The flux is then the only thing that changes,
    // including the timestep, which mixture_timestep sizes from the signal speed
    // of whichever flux is in force.
    if (const char* riemann = std::getenv("CHROMO_NUM_RIEMANN")) {
        const std::string choice(riemann);
        if (choice == "roe") {
            grid.roe_characteristic_flux = true;
            grid.swmf_godunov_flux = false;
        } else if (choice != "godunov") {
            throw std::invalid_argument(
                "CHROMO_NUM_RIEMANN must be godunov (the release) or roe");
        }
        std::cerr << "[model_column] Riemann override: CHROMO_NUM_RIEMANN=" << choice
                  << " (release default: godunov)" << std::endl;
    }
    // One line naming the numerical method actually in force, so a run log always
    // records it without having to enable the face-flux sidecar.
    std::cerr << "[model_column] numerics: flux="
              << (grid.swmf_godunov_flux ? "swmf-godunov" : "roe-local")
              << " reconstruction=lnrho-v-lnp hydrostatic=1 limiter=mc3"
              << " beta=" << grid.limiter_beta << std::endl;

    // --- runtime physics toggles ------------------------------------------
    // The RELEASE branch is deliberately a straight line: hydro plus physical
    // conduction, nothing else. The optional two-fluid research stages (finite-
    // rate ionization, radiative sink, TRAC, imposed q(T) flux, volumetric coronal
    // heating) live with the legacy model_gentle scenario, and the release
    // rejects their env knobs above.
    grid.single_fluid             = true;
    grid.enable_ionization        = false;
    grid.enable_radiative_cooling = false;
    grid.enable_beam_heating      = false;
    grid.enable_coronal_heating   = false;
    grid.impose_outer_heat_flux   = false;   // heat enters via the ghost-T wall
    grid.enable_Te                = false;
    grid.enable_trac              = false;   // no TRAC stage in the release solver
    grid.enable_vacuum_floor      = false;

    kOuterPRefSet   = false;        // fresh IC ⇒ recapture the outer back-pressure
    // CHROMO_CONDUCTION=0 runs the column with conduction off (an adiabatic hydrostatic
    // check); otherwise conduction is on from t = 0, driven by the outer wall.
    grid.enable_conduction        = heat_flux_on;
    // The release conduction operator is structurally physical-only and never
    // reads this field; the mesh-scaled artificial diffusivity belongs to the
    // historical two-fluid Stage D only.
    grid.numerical_diffusivity_per_length = 0.0f;

    model_column_update_bc(grid, xn);
    grid.broadcast();

    // The release is REFERENCE-FREE: it does not know a hydrostatic equilibrium
    // in advance and does not subtract a frozen residual. It is nevertheless
    // exactly balanced at rest: the hydrostatic reconstruction (set above) makes
    // any column in the discrete relation ln p_{i+1} - ln p_i = o_i + o_{i+1} a
    // fixed point, built from the live state each step, so the t = 0 face mass
    // flux is round-off (1.5e-19 kg m^-2 s^-1 at N=500). The history below is the
    // residual of the plain (ln rho, u, ln p) scheme it replaced. The release formerly
    // froze this (V=0, HSE) IC as grid.eq_state and subtracted R(U_eq) every
    // step, which made that one state an exact discrete fixed point but bound
    // the solver to a pre-known global equilibrium and kept applying a residual
    // evaluated at t = 0 (and at the first step's Δt) long after the transition
    // region had restructured. Completing the MUSCL-Hancock source treatment and
    // making the lower ghost ladder second order dropped the raw discrete
    // hydrostatic face mass-flux defect of the release mesh by ~560x, so the
    // frozen reference no longer has anything material to remove. Under the
    // release double-precision state that residual is 6.0e-14 kg m^-2 s^-1 at
    // N=500 (1.78e-13 before the MC3 branch fix) -- far above the double
    // round-off floor, i.e. pure discretization, and second-order convergent
    // (order 1.95 to N=1000). The old "~1.3x of the float32 round-off floor"
    // reading was a property of float32 storage.
    // Evidence: docs/studies/numerics/reference_free_release_recap.md,
    // docs/studies/numerics/hydrostatic_reconstruction_and_mc3_fix.md.
    //
    // The legacy two-fluid column (model_gentle) is a different solver that
    // never received the predictor source term, so it keeps the frozen
    // reference. eq_wb is a two-fluid-only mechanism.
    grid.eq_wb = false;
    return xn;
}

void model_column_update_bc(Grid& grid, const Vec& xn) {
    const Real g   = grid.g;

    // Decode total density / pressure / temperature / velocity / ionization
    // fraction of an interior cell straight out of (rho, rho u, E).
    auto cell_state = [&](arma::uword i, Real& rho_tot, Real& p_tot,
                          Real& T, Real& V, Real& x) {
        const auto msz = arma::size(grid.ns, num_of_mixture_eq);
        const MixtureThermo th = decode_equilibrium_mixture(
            grid.eos_gamma_table,
            xn(arma::sub2ind(msz, i, mix::RHO)),
            xn(arma::sub2ind(msz, i, mix::MOM)),
            xn(arma::sub2ind(msz, i, mix::ENERGY)),
            mixture_cell_phi(grid, i),
            std::numeric_limits<double>::quiet_NaN(), grid.eos_gamma_debug_clamp);
        rho_tot = th.rho; p_tot = th.p; T = th.T;
        V = static_cast<double>(xn(arma::sub2ind(msz, i, mix::MOM)))/th.rho;
        x = th.x;
    };

    // ====================================================================
    // Outer face — bottom of the transition region (upper BC):
    //   1. pressure    — hydrostatic, ONE-SIDED from the top cell (see below).
    //   2. temperature — the fixed reservoir temperature kTtopRef in both ghosts.
    //   3. density     — from the EOS at (p_ghost, T_ghost).
    // ====================================================================
    {
        const arma::uword nl = grid.ns - 1;
        const Real ds       = grid.ds_i(nl);
        const Real phi_g_out = grid.phi_g_iph(nl);   // ghost decode uses this gauge

        Real rho_top, p_top, T_top, V_top, x_top;
        cell_state(nl, rho_top, p_top, T_top, V_top, x_top);

        // (2) temperature. One fixed ghost temperature serves both the hydro ghosts
        // (their EOS density and internal energy, hence the outer-face
        // reconstruction and Riemann state) and conduction, which reads the same
        // ghost cell as its Dirichlet datum at the ghost centre. The outer face
        // therefore fixes pressure AND temperature; that is inert as a hydro
        // condition because conduction thermostats the top cell (~0.1 s), see
        // docs/studies/boundaries/outer_temperature_bc_single_t_recap.md.
        const Real T_g0 = kTtopRef;
        const Real T_g1 = T_g0;
        grid.outer_conduction_temperature_override = false;
        grid.outer_conduction_temperature          = kTtopRef;

        // (1)+(3) Hydrostatic ghost ladder: each rung one-sided on the rung below,
        //   p_{k+1} = p_k − Δs·½(ρ_k + ρ_{k+1})·g,
        // with ρ closed by the EOS at the imposed T. The ladder is ANCHORED on the
        // fixed reservoir back-pressure kOuterPRef, NOT on a live interior cell.
        //
        // Both live-cell anchors are wrong, in opposite directions:
        //   * p[ns-2] (the previous CENTERED form, p_g0 = p[ns-2] − 2Δs·ρ_top·g)
        //     extrapolates ACROSS the very cell it bounds. That is well balanced only
        //     while the interior satisfies the same discrete HSE; inside the transition
        //     region it does not — the measured residual (dp/ds+ρg)/(ρg) runs +0.03 at
        //     2140 km, +0.21 at 2145 km and −5.0 in the top cell — so the ghost landed
        //     ABOVE p_top and INVERTED the pressure gradient across the outer face
        //     (+4.1e-5 Pa where it should be −4.3e-6): a standing ~9 g downward push
        //     that reversed the top-cell velocity and seeded a ~25 km boundary layer.
        //   * p_top (a one-sided anchor) makes pressure zero-gradient. Combined with the
        //     already-extrapolated V that leaves the outer face transmissive in EVERY
        //     variable but T, so nothing sets the back-pressure and the column drains
        //     out the top (measured: V → 1.9 km/s, p falling monotonically, crash at
        //     t ≈ 12 s). Subsonic outflow has one incoming characteristic and therefore
        //     needs exactly one externally imposed condition.
        // A fixed reference back-pressure supplies that one condition, cannot invert
        // against the interior, and is exactly hydrostatic at the V=0 IC by
        // construction. Physically: a corona of fixed pressure and T sits above 2153 km.
        // The ghost density inverts the equilibrium Saha closure at the imposed
        // ghost temperature.
        auto ghost_density = [&](Real p_gh, Real T_gh) {
            return static_cast<Real>(gamma_density_from_pressure(grid, p_gh, T_gh));
        };
        // ρ_ghost depends on p_ghost, so take one fixed-point pass seeded with the
        // rung below (the correction is O(Δs·Δρ/ρ) and converges in a single sweep).
        auto hse_rung = [&](Real p_below, Real rho_below, Real T_gh,
                            Real& p_out, Real& rho_out) {
            Real rho_gh = rho_below;
            p_out = p_below;
            for (int it = 0; it < 40; ++it) {
                p_out = p_below*std::exp(-0.5*ds*solver_g(grid, nl)*(rho_below/p_below + rho_gh/p_out));
                if (p_out < 1.0e-12f) p_out = 1.0e-12f;
                rho_gh = ghost_density(p_out, T_gh);
            }
            rho_out = rho_gh;
        };
        // Capture the reservoir back-pressure once, one-sided from the (hydrostatic,
        // V=0) initial top cell, at the ghost temperature T_g0, so the ladder is
        // exactly the hydrostatic one this BC reproduces at the initial state.
        if (!kOuterPRefSet) {
            Real p_ref, rho_ref;
            hse_rung(p_top, rho_top, T_g0, p_ref, rho_ref);
            kOuterPRef    = p_ref;
            kOuterPRefSet = true;
            // The reservoir back-pressure is the one externally imposed condition at
            // the outer face; log it (once per capture) so post-run diagnostics can
            // form the top-cell pressure deficit p_top − p_g0 without re-deriving it.
            std::cerr << "[model_column] outer BC anchor: kOuterPRef=" << kOuterPRef
                      << " Pa  p_top=" << p_top << " Pa  T_top=" << T_top
                      << " K  T_ghost=" << T_g0 << " K  outer_p="
                      << (kOuterPHydrostatic ? "hydrostatic" : "reservoir")
                      << std::endl;
        }
        // The outer ghost pressures. Release: the fixed reservoir back-pressure,
        // then one hydrostatic rung up. Hydrostatic variant (kOuterPHydrostatic):
        // the same ladder re-anchored on the live top cell, so the externally
        // imposed condition is the ghost temperature instead.
        Real p_g0 = kOuterPRef;
        if (kOuterPHydrostatic) {
            Real rho_ladder;                      // ladder off the LIVE top cell
            hse_rung(p_top, rho_top, T_g0, p_g0, rho_ladder);
        }
        if (p_g0 < 1.0e-12f) p_g0 = 1.0e-12f;
        const Real rho_g0 = ghost_density(p_g0, T_g0);
        Real p_g1, rho_g1;
        hse_rung(p_g0, rho_g0, T_g1, p_g1, rho_g1);

        // Velocity. Conduction on: Mach-capped OUTFLOW so the evaporation upflow can
        // leave, while |V| ≤ kVcapMach·c_s prevents the ill-posed-inflow runaway.
        // Conduction off: a static reservoir, V = 0.
        Real V_g = 0.0f;
        if (kHeatFluxOn) {
            const Real c_s = std::sqrt(gamma_state(grid.eos_gamma_table, rho_top, T_top,
                                        grid.eos_gamma_debug_clamp).gamma_sound*p_top/rho_top);
            const Real vcap = kVcapMach * c_s;
            V_g = V_top;
            if (V_g >  vcap) V_g =  vcap;
            if (V_g < -vcap) V_g = -vcap;
        }
        mixture_pack_ghost(grid, grid.mix_outer_boundary0, rho_g0, V_g, T_g0,
                           phi_g_out);
        mixture_pack_ghost(grid, grid.mix_outer_boundary1, rho_g1, V_g, T_g1,
                           phi_g_out);
    }

    // ====================================================================
    // Inner face — the LOWER TRUNCATION of the modeled domain at h_base
    // (1600 km). Not the photosphere: it is the height below which this
    // solver's Saha/LTE equilibrium closure is not intended to be trusted, so
    // the atmosphere underneath is deliberately outside the model and the two
    // ghosts are a numerical closure for the MUSCL stencil.
    //
    // The physical content is one statement: the lower chromosphere below the
    // cut is a quasi-static stratified reservoir at temperature T₀. That fixes
    // the ghost thermal state, the hydrostatic ghost pressures kInnerPRes0/1
    // (trapezoidal, EOS-closed — see gamma_hse_rung_down) and V = 0.
    //
    // A gravity-aware CHARACTERISTIC alternative was prototyped and measured:
    // impose only the incoming (upward, u+c) acoustic invariant δu + δp/Z = 0
    // about the base reference state with impedance Z = ρc, and take the
    // outgoing (u−c) one from the live base cell, giving ghost deviations
    // V_g = −(p_0 − p_0,ref)/Z and p_g = p_g,ref − Z·V_0. It degenerates
    // exactly to this wall at equilibrium. Over 20 s conduction-off, 100 s and
    // 1000 s conduction-driven runs it changed the hydrostatic residual not at
    // all, the evaporation window velocity by ≤0.4 %, and the detrended
    // base-region oscillation content not at all: the disturbance that reaches
    // 1600 km is a slow quasi-static pressure adjustment, not acoustic ringing,
    // so this face is not an important reflector for this model. It was
    // therefore NOT promoted — see main.tex for the measurement.
    //
    // Conduction sees the same reservoir: a Dirichlet datum T₀ at the ghost
    // centre.
    // ====================================================================
    {
        const Real phi_g_in = grid.phi_g_imh(0);
        mixture_pack_ghost(grid, grid.mix_inner_boundary0, kInnerRhoRes0, 0.0,
                           kInnerTRef, phi_g_in);
        mixture_pack_ghost(grid, grid.mix_inner_boundary1, kInnerRhoRes1, 0.0,
                           kInnerTRef, phi_g_in);
    }

    grid.broadcast();
}

} // namespace chromosphere
