/*!
 * @file scenarios/model_column.cpp
 * @brief The field-aligned chromosphere-to-corona column: the `model_column`
 *        RELEASE scenario and its EXPERIMENTAL `model_column_2t` twin.
 * @ingroup scenarios
 *
 * One IC/BC implementation for the two Saha/Gamma1 scenarios. It loads the
 * production `Gamma1` table unconditionally, builds the Saha-HSE initial column
 * from the C7 temperature profile, sets the validated release numerics
 * (log-space MUSCL with the MC3/Koren limiter, the SWMF exact-Riemann Godunov
 * flux, the pressure-based reconstruction) and the 22,000 K external conductive
 * reservoir at the outer face — and rejects every historical two-fluid knob.
 * The historical fixed-gamma two-fluid column that used to share this file is
 * the legacy `model_gentle` scenario (legacy/scenarios/model_gentle.cpp). See
 * scenarios/model_column.hpp for the full environment-variable list.
 */
#include "model_column.hpp"
#include "single_fluid/mixture.hpp"
#include "two_temp/two_temp.hpp"   // EXPERIMENTAL model_column_2t leg only
#include "c7_profile.hpp"  // c7_full_profile + c7_full_temperature_pchip (C7 table)
#include "mesh.hpp"       // shared static local-refinement mesh builder

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace chromosphere {

// ---------------------------------------------------------------------------
// Module state captured by model_column_ic and reused by the boundary update.
// The interior is single-fluid neutral-dominated hydrogen with a realistic
// frozen ionization profile (neutral chromosphere, ionized corona).
// ---------------------------------------------------------------------------
namespace {
Real kAjump       = 1.0f;      // top temperature jump: T_ghost1 = a·T_top
Real kBjump       = 1.0f;      // T_ghost2 = b·T_ghost1
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
Real kTtopRef     = 0.0f;      // FIXED reference top temperature [K] for the jump
bool  kTtopFixedAbs = false;    // ISO_T_TOP>0: kTtopRef is an ABSOLUTE imposed top T
// FIXED reference outer ghost pressure [Pa] — the coronal reservoir back-pressure,
// captured one-sided-hydrostatically from the IC top cell (and recaptured at the
// relax→Stage-2 edge). Subsonic outflow admits exactly ONE incoming characteristic,
// so exactly one condition may come from outside: this back-pressure. Anchoring it to
// a FIXED reference rather than to a live interior cell is what makes the outer face
// well posed — see the outer-face block for why both live-cell anchors fail.
Real kOuterPRef    = 0.0f;
bool  kOuterPRefSet = false;    // captured lazily on the first update_bc call
Real kHeatFluxOn  = false;     // Stage 2 (conduction + jump) ⇒ open (Mach-capped) top
Real kVcapMach    = 0.1f;      // outer outflow cap as a fraction of c_s
bool  kInnerTNeumann = false;   // lower-BC conduction temperature: Dirichlet vs Neumann
// Relaxation phase (ISO_RELAX_TIME > 0): Stage-1 adiabatic relaxation (heat flux off,
// no jump, V=0 reservoir top) for sim_time < kRelaxTime, THEN Stage 2.
Real kRelaxTime     = 0.0f;
bool  kHeatFluxTarget = false;
bool  kStage2Started  = true;   // false during relax; flips to recapture kTtopRef
// ISO_HYDRO_T_DECOUPLE=1: split the two roles the outer ghost temperature plays.
// The HYDRO ghost temperature becomes a zero-gradient extrapolation of the live top
// cell (T_hydro_g0 = T_hydro_g1 = T_top) instead of the fixed a·T_ref wall, while
// CONDUCTION keeps seeing exactly the same fixed wall a·T_ref through
// grid.outer_conduction_temperature. Subsonic outflow admits ONE incoming
// characteristic and the fixed reservoir back-pressure kOuterPRef already supplies
// it; the additional fixed hydro T pins the thermal/entropy state as well, which is
// the suspected source of the ~20-30 km top boundary layer. This flag removes ONLY
// that extra hydro constraint, holding the conductive energy input fixed.
// Code default 0, and the release preset also sets 0: on the N=500/R4 mesh the top
// cell relaxes to the wall in ~0.1 s, so the extra hydro constraint is inert
// (V_top within 1 %) — see outer_temperature_bc_single_t_recap.md.
bool  kHydroTDecouple = false;
// ISO_OUTER_T_HYDRO_WALL=1 (decoupled mode only): the HYDRO ghost temperature is
// PINNED to the same fixed wall the conduction rows see, T_hydro = T_cond = a·T_ref,
// instead of zero-gradient-extrapolating the live top cell. The wall stays a
// PHYSICAL-FACE conduction datum (unlike baseline decouple=0, where the wall reverts
// to a ghost-CENTRE datum), so relative to decouple=1 this moves exactly one thing:
// the thermal state the outer-face reconstruction and Riemann problem carry.
bool  kOuterHydroTWall = false;
// ISO_OUTER_P_NEUMANN: the outer ghost pressure stops being the FIXED reservoir
// back-pressure kOuterPRef and is extrapolated off the LIVE top cell instead — a
// Neumann pressure closure. Subsonic outflow admits exactly ONE incoming
// characteristic, so exactly one condition may be imposed from outside; this swaps
// WHICH variable carries it, from the back-pressure to the temperature, and is
// therefore only well posed with a FIXED hydro ghost temperature
// (ISO_OUTER_T_HYDRO_WALL=1, or baseline decouple=0). Two readings of "no imposed
// pressure", both selectable because they are physically different statements:
//   1 = dp/ds = 0            — p_g1 = p_g0 = p_top, the literal zero-gradient BC.
//       Note this is NOT gravitationally consistent: the interior needs
//       dp/ds = −ρg, so the ghost withholds the hydrostatic support of the top cell.
//   2 = dp/ds = −ρg          — the same trapezoidal EOS-closed ladder the release
//       uses (hse_rung), but ANCHORED ON THE LIVE TOP CELL instead of kOuterPRef.
//       Hydrostatically consistent, and imposes no external pressure datum at all.
int   kOuterPNeumann = 0;
// EXPERIMENTAL two-temperature leg (`model_column_2t`). Set only by
// model_column_2t_ic; leaves the release leg untouched.
bool  kTwoTemperature = false;

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

// One rung of an EOS-consistent hydrostatic ladder taken DOWNWARD by `ds`
// metres from an anchor (p_anchor, rho_anchor) at fixed temperature T_gh:
//     p = p_anchor + Δs·½(ρ_anchor + ρ)·g,     ρ = ρ_EOS(p, T_gh),
// iterated to a fixed point. Trapezoidal in ρ, so the ladder reproduces the
// exponential stratification to O(Δs²) instead of the O(Δs) of p + ρ_anchor gΔs.
void gamma_hse_rung_down(const Grid& grid, double p_anchor, double rho_anchor,
                         double T_gh, double ds, double& p_out, double& rho_out) {
    double rho_gh = rho_anchor;
    for (int it = 0; it < 8; ++it) {
        const double p = p_anchor + ds*0.5*(rho_anchor + rho_gh)*grid.g;
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
} // namespace

Vec model_column_ic(Grid& grid) {
    auto env_f = [](const char* key, Real fallback) -> Real {
        if (const char* e = std::getenv(key)) {
            try { return std::stof(std::string(e)); } catch (...) {}
        }
        return fallback;
    };

    if (const char* gamma_path = std::getenv("GAMMA_TABLE")) {
        if (std::getenv("ISO_GAMMA"))
            throw std::logic_error("GAMMA_TABLE and ISO_GAMMA are mutually exclusive");
        if (grid.eos_gamma_table.empty()) grid.eos_gamma_table = EosGammaTable::load(gamma_path);
    }
    // A loaded Gamma1 table is what selects the Saha/Gamma1 solvers. `model_column`
    // and `model_column_2t` always supply one, so an empty table here means the
    // caller asked for the historical fixed-gamma two-fluid column, which is the
    // separate legacy `model_gentle` scenario.
    if (grid.eos_gamma_table.empty())
        throw std::runtime_error(
            "model_column_ic: no Gamma1 table is loaded. model_column and "
            "model_column_2t always run a Saha/Gamma1 solver; the historical "
            "fixed-gamma two-fluid column is the legacy model_gentle scenario "
            "(build with -DCHROMO_ENABLE_LEGACY=ON)");

    // Release configuration surface. Everything below is historical two-fluid
    // research physics that the single-fluid release timestep does not contain:
    // a fixed adiabatic index, a neutral fluid, finite-rate ionization, a
    // radiative sink, TRAC broadening, volumetric coronal heating, an imposed
    // Neumann coronal flux, an IC coronal superheat, and the mesh-scaled
    // artificial conduction term. Rejecting them outright is what keeps
    // `model_column` a single unambiguous production identity — a release run
    // can never be silently reconfigured into something else.
    {
        static const char* const kNonReleaseKnobs[] = {
            "ISO_GAMMA", "ISO_TWO_FLUID", "ISO_IONIZATION", "ISO_COOLING",
            "ISO_TRAC", "ISO_CORONA", "ISO_CHEAT", "ISO_QFLUX", "ISO_TBOOST",
            "ISO_NUMERICAL_DIFFUSIVITY_MULT"};
        for (const char* knob : kNonReleaseKnobs) {
            if (std::getenv(knob))
                throw std::runtime_error(
                    std::string(knob) + " is historical two-fluid research "
                    "configuration and has no meaning for the single-fluid "
                    "release solver; run the legacy model_gentle scenario for "
                    "that path (build with -DCHROMO_ENABLE_LEGACY=ON)");
        }
    }

    // --- IC mode ----------------------------------------------------------
    // The IC is ALWAYS the real Model C7 profile with the density re-integrated
    // hydrostatically (clean V≈0 start; the analytic-isentrope toy IC is retired).

    // --- parameters -------------------------------------------------------
    // Lower truncation height of the modeled domain, in km measured from the
    // photosphere (that is the HEIGHT GAUGE only -- the domain itself starts
    // here, at 1600 km, well above it). See the inner-boundary block.
    const Real h_base     = env_f("ISO_H_BASE", 0.0f);
    Real       DH_km      = env_f("ISO_DH", 1303.0f); // km
    // Base (T, n_H) from the C7 profile at h_base.
    Real T_base, ni_b, nn_b;
    c7_cell(h_base, T_base, ni_b, nn_b);
    T_base = static_cast<Real>(c7_full_temperature_pchip(h_base));
    const bool  heat_flux_on = (env_f("ISO_HEAT_FLUX", 0.0f) != 0.0f);
    kHeatFluxOn        = heat_flux_on;
    kAjump             = env_f("ISO_TJUMP_A", 1.0f);
    kBjump             = env_f("ISO_TJUMP_B", 1.0f);
    kVcapMach          = env_f("ISO_VCAP", 0.1f);
    kInnerTNeumann     = (env_f("ISO_INNER_T_NEUMANN", 0.0f) != 0.0f);
    // Upper-BC over-specification experiment: hydro ghost T free-floats (zero
    // gradient), conduction wall unchanged. Default 0 ⇒ the release single-T BC.
    kHydroTDecouple    = (env_f("ISO_HYDRO_T_DECOUPLE", 0.0f) != 0.0f);
    // Outer-face closure variant: which ONE condition the subsonic-outflow face
    // imposes. Release = fixed back-pressure + fixed ghost T (decouple=0). Setting
    // both knobs is the swap: fixed hydro/conduction T + a Neumann pressure.
    kOuterHydroTWall   = (env_f("ISO_OUTER_T_HYDRO_WALL", 0.0f) != 0.0f);
    kOuterPNeumann     = static_cast<int>(env_f("ISO_OUTER_P_NEUMANN", 0.0f));
    if (kOuterPNeumann < 0 || kOuterPNeumann > 2)
        throw std::runtime_error(
            "model_column_ic: ISO_OUTER_P_NEUMANN must be 0 (fixed reservoir "
            "back-pressure, the release), 1 (dp/ds = 0) or 2 (dp/ds = -rho g, "
            "anchored on the live top cell)");
    if (kOuterPNeumann && kHydroTDecouple && !kOuterHydroTWall)
        throw std::runtime_error(
            "model_column_ic: a nonzero ISO_OUTER_P_NEUMANN removes the imposed "
            "outer pressure, so the one condition the subsonic-outflow face admits "
            "must come from the temperature. Combine it with "
            "ISO_OUTER_T_HYDRO_WALL=1 (fixed hydro ghost T at the conduction wall) "
            "or with ISO_HYDRO_T_DECOUPLE=0; with the decoupled extrapolated hydro "
            "temperature the face would impose nothing and the column drains.");

    // The release is a single-fluid common-temperature equilibrium mixture.
    grid.single_fluid = true;
    grid.enable_ionization = false;
    grid.enable_Te = false;

    // --- geometry: straight field line, gravity on -----------------------
    // Static local refinement (scenarios/mesh.hpp). peek_ns() sized the grid to the
    // coarse-equivalent count (ISO_NS); ISO_REFINE_* ADD cells inside the refined band
    // and grade smoothly to the coarse spacing outside it. ISO_REFINE_PROFILE picks
    // which end the band is anchored to: `lower` (default, historical: fine [s_lo,s_hi]
    // then grade back up to coarse) or `outer` (fine from s_lo to the domain TOP, with
    // the grade in [s_lo−τ, s_lo]) — the latter targets the upper TR / outer boundary.
    // The conduction solver uses these actual mesh metrics to form the face-local
    // numerical diffusivity; ds_m remains useful only for uniform-grid geometry.
    // Uniform mesh reproduced byte-for-byte.
    const arma::uword ns_coarse = grid.ns;
    const Real ds_m = (DH_km * 1000.0f) / static_cast<Real>(ns_coarse);   // coarse Δs [m]
    const RefineParams rp = refine_params_from_env(/*iso_alias=*/true);

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
        double ln_p = ln_p_base, prev_h = h_base;
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
                    const double next = previous_ln_p-0.5*(prev_invH+invH)*dh_m;
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

        // ISO_P_TOP > 0 [Pa]: rebuild the column with the SAME C7 T(h) but the base
        // pressure chosen so the top cell sits at this pressure. The whole column
        // moves consistently (both reservoirs are captured from the rebuilt IC), so
        // the pressure is a controlled parameter rather than a one-cell mismatch.
        // Saha makes p_top a nonlinear function of p_base (a denser column is less
        // ionized and has a shorter scale height), hence the secant iteration.
        if (const char* e = std::getenv("ISO_P_TOP")) {
            double target = std::numeric_limits<double>::quiet_NaN();
            try { target = std::stod(e); } catch (...) {}
            if (!(std::isfinite(target) && target > 0.0))
                throw std::runtime_error(std::string("model_column_ic: ISO_P_TOP must be a "
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
                throw std::runtime_error("model_column_ic: ISO_P_TOP secant iteration "
                    "did not converge (target " + std::to_string(target) + " Pa)");
            std::cerr << "[model_column] ISO_P_TOP: top-cell p=" << std::exp(f1+ln_target)
                      << " Pa (target " << target << ", C7 " << std::exp(ln_p_top_c7)
                      << ")  base p=" << std::exp(a1) << " Pa (C7 " << std::exp(ln_p_c7)
                      << ")  secant iterations=" << it+1 << std::endl;
        }
    }

    Vec xn = arma::zeros<Vec>(
        kTwoTemperature ? grid.ns*num_of_two_temp_eq : grid.n_mixture_state);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const Real T     = T_c(i);
        const Real n_i   = ni_c(i);
        const Real n_n   = nn_c(i);
        const Real phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        const Real rho_i = n_i * grid.m_i;
        const Real rho_n = n_n * grid.m_n;
        if (kTwoTemperature) {
            // T_e = T_i at t = 0. With tau_eq of 0.09-3.9 ms across this column
            // against a step of tens of ms, a common initial temperature is the
            // physically correct start; any split the run develops is produced by
            // the solver, not imposed by the IC.
            two_temp_pack_cell(grid, xn, i, static_cast<double>(rho_i)+rho_n,
                               0.0, T, T, phi_g);
            continue;
        }
        mixture_pack_cell(grid, xn, i, static_cast<double>(rho_i)+rho_n,
                          0.0, T, phi_g);
    }

    // Reference top temperature for the imposed ghost jump (Stage 2), FIXED at the
    // IC top value: the unresolved TR/corona above sits at ~constant T, so
    // T_ghost = a·T_ref is a fixed hot reservoir and the downward conductive flux
    // q ∝ (T_ghost − T_top) SELF-LIMITS as the chromosphere heats.
    kTtopRef = T_c(grid.ns - 1);
    // ISO_T_TOP>0: override with an ABSOLUTE fixed outer-boundary temperature [K]
    // (e.g. 22000 = C7's TR-base T). Not re-anchored across a relax phase.
    const Real t_top_abs = env_f("ISO_T_TOP", 0.0f);
    kTtopFixedAbs = (t_top_abs > 0.0f);
    if (kTtopFixedAbs) kTtopRef = t_top_abs;

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

    // --- numerical scheme: the validated "bestwb" configuration (always on) ---
    // well-balanced + log-space MUSCL + MC3/Koren (β=2) + equilibrium-reference
    // δ-form well-balancing. These hold a hydrostatic column at V≈0 to round-off
    // for any stratification (docs/studies/evaporation/gentle_evaporation_downflow.md). No toggles.
    grid.out_base_km              = h_base;   // label output heights from the true base
    grid.inner_conduction_neumann = kInnerTNeumann;
    grid.well_balanced            = true;
    grid.log_reconstruct          = true;
    grid.mc3_limiter              = true;
    grid.limiter_beta             = 2.0f;
    // DIAGNOSTIC-ONLY reconstruction override for the scheme-comparison leg of the
    // face-flux study (docs/studies/numerics/top_ripple_face_flux_diagnosis.md). Unset — the default
    // — leaves the "bestwb" MC3(beta=2) configuration above completely untouched.
    //   ISO_LIMITER=mc3    MC3/Koren, beta from ISO_MC3_BETA (default 2) [= default]
    //   ISO_LIMITER=minmod symmetric minmod
    //   ISO_LIMITER=first  beta=0 with MC3 ⇒ phi≡0 ⇒ piecewise-constant (1st order)
    // Never set for a production run; it exists so the ripple can be attributed to
    // (or cleared of) the limiter without hand-editing the scenario.
    if (const char* lim = std::getenv("ISO_LIMITER")) {
        const std::string choice(lim);
        if (choice == "minmod") {
            grid.mc3_limiter = false;
        } else if (choice == "first") {
            grid.mc3_limiter = true;
            grid.limiter_beta = 0.0f;
        } else if (choice == "mc3") {
            grid.limiter_beta = env_f("ISO_MC3_BETA", 2.0f);
        } else {
            throw std::invalid_argument("ISO_LIMITER must be mc3, minmod or first");
        }
        std::cerr << "[model_column] DIAGNOSTIC reconstruction override: ISO_LIMITER="
                  << choice << " mc3=" << grid.mc3_limiter
                  << " beta=" << grid.limiter_beta << std::endl;
    }
    // --- release numerical flux and MUSCL primitive set (Gamma/Saha path) ---
    // The released Gamma/Saha model_column advances the corrector with the
    // SWMF-style exact-Riemann Godunov flux — the local Riemann problem solved at
    // the frozen-composition index gamma = 5/3 with the ionization energy carried
    // as SWMF's passive offset E0 — and limits (log rho, V, log p), recovering the
    // face temperature by inverting the same authoritative Saha closure. Both
    // require the equilibrium-manifold face builder, which every Saha/Gamma1 run
    // has (the legacy two-fluid column keeps Rusanov + (log rho, V, log T)).
    grid.swmf_godunov_flux       = true;
    grid.pressure_reconstruct    = true;
    grid.roe_characteristic_flux = false;
    // Reference/debug overrides. ISO_RIEMANN=roe-local selects the mixture Roe
    // characteristic flux — the equilibrium (Gamma1) linearization that was the
    // previous release flux and is kept as the controlled comparison solver —
    // ISO_RIEMANN=rusanov restores the local Lax-Friedrichs flux, and
    // ISO_RECONSTRUCTION=lnrho-v-lnt the temperature-based primitive set; all are
    // retained for regression and controlled numerical comparison only. None is
    // needed for a production run. The three flux choices are mutually exclusive
    // by construction here, so the face Riemann solver is the only variable they
    // change — including the timestep, which mixture_timestep sizes from the
    // signal speed of whichever flux is in force.
    if (const char* riemann = std::getenv("ISO_RIEMANN")) {
        const std::string choice(riemann);
        if (choice == "rusanov") {
            grid.roe_characteristic_flux = false;
            grid.swmf_godunov_flux = false;
        } else if (choice == "roe-local") {
            grid.roe_characteristic_flux = true;
            grid.swmf_godunov_flux = false;
        } else if (choice == "swmf-godunov") {
            grid.roe_characteristic_flux = false;
            grid.swmf_godunov_flux = true;
        } else {
            throw std::invalid_argument(
                "ISO_RIEMANN must be rusanov, roe-local or swmf-godunov");
        }
        std::cerr << "[model_column] Riemann override: ISO_RIEMANN=" << choice
                  << " (release default: swmf-godunov)" << std::endl;
    }
    if (const char* recon = std::getenv("ISO_RECONSTRUCTION")) {
        const std::string choice(recon);
        if (choice == "lnrho-v-lnt") {
            grid.pressure_reconstruct = false;
        } else if (choice == "lnrho-v-lnp") {
            grid.pressure_reconstruct = true;
        } else {
            throw std::invalid_argument(
                "ISO_RECONSTRUCTION must be lnrho-v-lnt or lnrho-v-lnp");
        }
        std::cerr << "[model_column] reconstruction override: ISO_RECONSTRUCTION="
                  << choice << " (release default: lnrho-v-lnp)" << std::endl;
    }
    // One line naming the numerical method actually in force, so a run log always
    // records it without having to enable the face-flux sidecar.
    std::cerr << "[model_column] numerics: flux="
              << (grid.swmf_godunov_flux ? "swmf-godunov"
                  : grid.roe_characteristic_flux ? "roe-local" : "rusanov")
              << " reconstruction="
              << (grid.pressure_reconstruct ? "lnrho-v-lnp" : "lnrho-v-lnt")
              << " limiter=" << (grid.mc3_limiter ? "mc3" : "minmod")
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

    // Optional relaxation phase: run Stage-1 adiabatic relaxation (conduction off,
    // no jump, V=0 reservoir) for sim_time < kRelaxTime so the C7 IC settles to a
    // near-hydrostatic V≈0 equilibrium, THEN start Stage 2 (switched per step in
    // model_column_update_bc off grid.sim_time).
    kRelaxTime      = env_f("ISO_RELAX_TIME", 0.0f);
    kHeatFluxTarget = heat_flux_on;
    const bool relaxing0 = (kRelaxTime > 0.0f);
    kStage2Started  = !relaxing0;   // recapture kTtopRef at the relax→Stage-2 edge
    kOuterPRefSet   = false;        // fresh IC ⇒ recapture the outer back-pressure
    kHeatFluxOn     = relaxing0 ? false : heat_flux_on;

    // Stage 1: conduction OFF (adiabatic Euler steady state). Stage 2: ON, driven
    // by the top T wall.
    grid.enable_conduction        = relaxing0 ? false : heat_flux_on;
    // The release conduction operator is structurally physical-only and never
    // reads this field; the mesh-scaled artificial diffusivity belongs to the
    // historical two-fluid Stage D only.
    grid.numerical_diffusivity_per_length = 0.0f;

    model_column_update_bc(grid, xn);
    grid.broadcast();

    // The release is REFERENCE-FREE: it does not know a hydrostatic equilibrium
    // in advance and does not subtract a frozen residual. The release formerly
    // froze this (V=0, HSE) IC as grid.eq_state and subtracted R(U_eq) every
    // step, which made that one state an exact discrete fixed point but bound
    // the solver to a pre-known global equilibrium and kept applying a residual
    // evaluated at t = 0 (and at the first step's Δt) long after the transition
    // region had restructured. Completing the MUSCL-Hancock source treatment and
    // making the lower ghost ladder second order dropped the raw discrete
    // hydrostatic face mass-flux defect of the release mesh by ~560x, so the
    // frozen reference no longer has anything material to remove. Under the
    // release double-precision state that residual is 1.78e-13 kg m^-2 s^-1 at
    // N=500 -- ~1e8x above the double round-off floor, i.e. pure discretization,
    // and second-order convergent (order 2.19 to N=1000). The old "~1.3x of the
    // float32 round-off floor" reading was a property of float32 storage.
    // Evidence: docs/studies/numerics/reference_free_release_recap.md.
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
    // fraction of an interior cell.
    // Release state: decode the mixture cell straight out of (rho, rho u, E).
    auto mixture_cell_state = [&](arma::uword i, Real& rho_tot, Real& p_tot,
                                  Real& T, Real& V, Real& x,
                                  Real& T_e, Real& T_i) {
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
        T_e = th.T; T_i = th.T;   // the release is a common-temperature mixture
    };

    // EXPERIMENTAL two-temperature leg. `T` is reported as the heavy-particle
    // temperature so every downstream use that means "the gas temperature"
    // (the Mach cap sound speed, the reference top-T anchor) keeps its meaning.
    auto two_temp_cell_state = [&](arma::uword i, Real& rho_tot, Real& p_tot,
                                   Real& T, Real& V, Real& x,
                                   Real& T_e, Real& T_i) {
        const auto tsz = arma::size(grid.ns, num_of_two_temp_eq);
        const TwoTempThermo th = two_temp_decode(
            grid.eos_gamma_table,
            xn(arma::sub2ind(tsz, i, tt::RHO)),
            xn(arma::sub2ind(tsz, i, tt::MOM)),
            xn(arma::sub2ind(tsz, i, tt::ENERGY)),
            xn(arma::sub2ind(tsz, i, tt::E_ELEC)),
            two_temp_cell_phi(grid, i),
            std::numeric_limits<double>::quiet_NaN());
        rho_tot = th.rho; p_tot = th.p; T = th.T_i;
        V = static_cast<double>(xn(arma::sub2ind(tsz, i, tt::MOM)))/th.rho;
        x = th.x;
        T_e = th.T_e; T_i = th.T_i;
    };

    // The release (common-temperature mixture) unless the EXPERIMENTAL
    // two-temperature leg is active. model_column_ic guarantees a Gamma1 table.
    auto cell_state = [&](arma::uword i, Real& rho_tot, Real& p_tot,
                          Real& T, Real& V, Real& x,
                          Real& T_e, Real& T_i) {
        if (kTwoTemperature) two_temp_cell_state(i, rho_tot, p_tot, T, V, x, T_e, T_i);
        else                 mixture_cell_state(i, rho_tot, p_tot, T, V, x, T_e, T_i);
    };

    // --- relaxation → Stage-2 phase switch (off grid.sim_time) -------------
    // While relaxing: conduction OFF, no jump, V=0 reservoir top. At the transition
    // recapture kTtopRef from the RELAXED top cell (anchor the FIXED jump to the
    // equilibrium), then run Stage 2 (conduction ON + jump + Mach-capped outflow).
    const bool relaxing = (kRelaxTime > 0.0f) && (grid.sim_time < kRelaxTime);
    kHeatFluxOn            = relaxing ? false : kHeatFluxTarget;
    grid.enable_conduction = relaxing ? false : kHeatFluxTarget;
    const Real a_live = relaxing ? 1.0f : kAjump;
    const Real b_live = relaxing ? 1.0f : kBjump;
    if (!relaxing && !kStage2Started) {
        if (!kTtopFixedAbs) {           // an ABSOLUTE imposed top T (ISO_T_TOP) stays put
            Real rho_t, p_t, T_t, V_t, x_t, Te_t, Ti_t;
            cell_state(grid.ns - 1, rho_t, p_t, T_t, V_t, x_t, Te_t, Ti_t);
            kTtopRef = T_t;             // else anchor the jump to the relaxed top T
        }
        kOuterPRefSet = false;          // re-anchor the back-pressure on the relaxed top
        kStage2Started = true;
    }

    // ====================================================================
    // Outer face — bottom of the transition region (new upper BC):
    //   1. pressure    — hydrostatic, ONE-SIDED from the top cell (see below).
    //   2. temperature — imposed jump T_ghost1 = a·T_ref, T_ghost2 = b·T_ghost1.
    //   3. density     — from the EOS at (p_ghost, T_ghost).
    // ====================================================================
    {
        const arma::uword nl = grid.ns - 1;
        const Real ds       = grid.ds_i(nl);
        const Real phi_g_out = grid.phi_g_iph(nl);   // ghost decode uses this gauge

        Real rho_top, p_top, T_top, V_top, x_top, Te_top, Ti_top;
        cell_state(nl, rho_top, p_top, T_top, V_top, x_top, Te_top, Ti_top);

        // (2) temperature. TWO distinct roles, kept explicitly separate:
        //   * T_cond_wall — the fixed hot wall that DRIVES Stage-D conduction, always
        //     the imposed jump a·T_ref off the reference top T. With ISO_CORONA the
        //     corona is resolved and a=b=1 (the ghost just extends it).
        //   * T_hydro_g0/g1 — what the HYDRO ghost carries (its EOS density and
        //     internal energy, hence the outer-face reconstruction and Riemann state).
        // Baseline (default, and the release): the hydro ghost IS the wall,
        // T_hydro = a·T_ref, b·a·T_ref — one ghost state serving both roles, so the
        // outer face fixes pressure AND temperature. ISO_HYDRO_T_DECOUPLE=1: the hydro ghost instead zero-gradient-
        // extrapolates the live top cell (T_hydro_g0 = T_hydro_g1 = T_top) and only the
        // conduction rows keep the wall, via grid.outer_conduction_temperature.
        // ISO_OUTER_T_HYDRO_WALL=1 (decoupled mode): the hydro ghost carries the wall
        // value again, T_hydro = T_cond = a·T_ref, while the wall REMAINS a physical-
        // face conduction datum — the fixed-T half of the pressure/temperature swap.
        const Real T_cond_wall = a_live * kTtopRef;
        const bool hydro_T_free = kHydroTDecouple && !kOuterHydroTWall;
        const Real T_g0 = hydro_T_free ? T_top : T_cond_wall;
        const Real T_g1 = kHydroTDecouple ? T_g0  : b_live * T_g0;
        grid.outer_conduction_temperature_override = kHydroTDecouple;
        grid.outer_conduction_temperature          = T_cond_wall;
        // EXPERIMENTAL two-temperature leg. The hydro ghost carries BOTH
        // temperatures, and in the model_column_2t preset (hydro_T_free) BOTH free-float
        // with the top cell — which is the characteristically correct closure: the
        // 4x4 system's eigenvalues are {u-c, u, u, u+c}, so at this subsonic
        // OUTFLOW face three characteristics leave and only u-c enters, and the
        // single condition it admits is already spent on the reservoir
        // back-pressure kOuterPRef. The 22 000 K wall reaches the gas only through
        // the ELECTRON conduction channel (grid.outer_conduction_temperature).
        const Real Te_g0 = hydro_T_free ? Te_top : T_cond_wall;
        const Real Ti_g0 = hydro_T_free ? Ti_top : T_cond_wall;
        const Real Te_g1 = kHydroTDecouple ? Te_g0 : b_live * Te_g0;
        const Real Ti_g1 = kHydroTDecouple ? Ti_g0 : b_live * Ti_g0;

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
        // The ghost density inverts whichever pressure closure the active solver
        // uses at the imposed ghost temperature(s): the two-temperature closure
        // p = n_H k_B (T_i + x(n_H,T_e) T_e) on the experimental leg, the
        // equilibrium Saha closure on the release leg.
        auto ghost_density = [&](Real p_gh, Real T_e_gh, Real T_i_gh) {
            if (kTwoTemperature)
                return static_cast<Real>(two_temp_density_from_pressure(
                    grid.eos_gamma_table, p_gh, T_e_gh, T_i_gh));
            return static_cast<Real>(gamma_density_from_pressure(grid, p_gh, T_i_gh));
        };
        // ρ_ghost depends on p_ghost, so take one fixed-point pass seeded with the
        // rung below (the correction is O(Δs·Δρ/ρ) and converges in a single sweep).
        auto hse_rung = [&](Real p_below, Real rho_below, Real T_e_gh, Real T_i_gh,
                            Real& p_out, Real& rho_out) {
            Real rho_gh = rho_below;
            for (int it = 0; it < 2; ++it) {
                p_out = p_below - ds * 0.5f * (rho_below + rho_gh) * g;
                if (p_out < 1.0e-12f) p_out = 1.0e-12f;
                rho_gh = ghost_density(p_out, T_e_gh, T_i_gh);
            }
            rho_out = rho_gh;
        };
        // Capture the reservoir back-pressure once, one-sided from the (hydrostatic,
        // V=0) reference top cell. Recaptured at the relax→Stage-2 edge alongside
        // kTtopRef so a Stage-1 relaxation re-anchors on the settled equilibrium.
        // The capture uses the HYDRO ghost temperature T_g0 (not the conduction wall),
        // so the ladder is exactly the hydrostatic one this BC will later reproduce:
        // at the reference state T_g0 = T_top there, and BC(q_ref) = q_g,ref holds in
        // decoupled mode exactly as it does in baseline mode.
        if (!kOuterPRefSet) {
            Real p_ref, rho_ref;
            hse_rung(p_top, rho_top, Te_g0, Ti_g0, p_ref, rho_ref);
            kOuterPRef    = p_ref;
            kOuterPRefSet = true;
            // The reservoir back-pressure is the one externally imposed condition at
            // the outer face; log it (once per capture) so post-run diagnostics can
            // form the top-cell pressure deficit p_top − p_g0 without re-deriving it.
            std::cerr << "[model_column] outer BC anchor: kOuterPRef=" << kOuterPRef
                      << " Pa  p_top=" << p_top << " Pa  T_top=" << T_top
                      << " K  T_hydro_g0=" << T_g0 << " K  T_cond_wall="
                      << T_cond_wall << " K  hydro_T_decouple="
                      << (kHydroTDecouple ? 1 : 0)
                      << "  hydro_T_wall=" << (kOuterHydroTWall ? 1 : 0)
                      << "  p_neumann=" << kOuterPNeumann << std::endl;
        }
        // The outer ghost pressures. Mode 0 (release): the fixed reservoir back-
        // pressure, then one hydrostatic rung up. Mode 1: dp/ds = 0, so the ghost pair
        // is uniform in p and the reconstruction sees no imposed pressure slope either.
        // Mode 2: the same hydrostatic ladder as mode 0, re-anchored on the live top
        // cell. In modes 1 and 2 the externally imposed condition is T_g instead.
        Real p_g0 = kOuterPNeumann == 1 ? p_top : kOuterPRef;
        if (kOuterPNeumann == 2) {
            Real rho_ladder;                      // ladder off the LIVE top cell
            hse_rung(p_top, rho_top, Te_g0, Ti_g0, p_g0, rho_ladder);
        }
        if (p_g0 < 1.0e-12f) p_g0 = 1.0e-12f;
        const Real rho_g0 = ghost_density(p_g0, Te_g0, Ti_g0);
        Real p_g1, rho_g1;
        if (kOuterPNeumann == 1) {
            p_g1   = p_g0;
            rho_g1 = ghost_density(p_g1, Te_g1, Ti_g1);
        } else {
            hse_rung(p_g0, rho_g0, Te_g1, Ti_g1, p_g1, rho_g1);
        }

        // Velocity. Stage 1 (relaxation): static reservoir V=0. Stage 2: Mach-capped
        // OUTFLOW so the evaporation upflow can leave, while |V| ≤ kVcapMach·c_s
        // prevents the ill-posed-inflow runaway.
        Real V_g = 0.0f;
        if (kHeatFluxOn) {
            // The two-temperature leg caps on its own FROZEN signal speed
            // sqrt(5/3 p/rho) — the speed its Godunov flux and its CFL actually
            // use — rather than the equilibrium Gamma1 speed of the release.
            const Real c_s = kTwoTemperature
                ? std::sqrt((5.0/3.0)*p_top/rho_top)
                : std::sqrt(gamma_state(grid.eos_gamma_table, rho_top, T_top,
                                        grid.eos_gamma_debug_clamp).gamma_sound*p_top/rho_top);
            const Real vcap = kVcapMach * c_s;
            V_g = V_top;
            if (V_g >  vcap) V_g =  vcap;
            if (V_g < -vcap) V_g = -vcap;
        }
        if (kTwoTemperature) {
            two_temp_pack_ghost(grid, grid.tt_outer_boundary0, rho_g0, V_g,
                                Te_g0, Ti_g0, phi_g_out);
            two_temp_pack_ghost(grid, grid.tt_outer_boundary1, rho_g1, V_g,
                                Te_g1, Ti_g1, phi_g_out);
        } else {
            mixture_pack_ghost(grid, grid.mix_outer_boundary0, rho_g0, V_g, T_g0,
                               phi_g_out);
            mixture_pack_ghost(grid, grid.mix_outer_boundary1, rho_g1, V_g, T_g1,
                               phi_g_out);
        }
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
    // The lower-BC TEMPERATURE choice (Dirichlet vs Neumann/zero-flux) is a
    // separate CONDUCTION boundary condition imposed inside the Stage-D solver
    // (grid.inner_conduction_neumann); the hydro closure is independent of it.
    // ====================================================================
    {
        const Real phi_g_in = grid.phi_g_imh(0);
        if (kTwoTemperature) {
            // The reservoir is dense, cool and collisionally equilibrated
            // (tau_eq = 0.09 ms at 1600 km), so T_e = T_i = T_0 there is a
            // statement about the reservoir, not an extra modelling assumption.
            two_temp_pack_ghost(grid, grid.tt_inner_boundary0, kInnerRhoRes0, 0.0,
                                kInnerTRef, kInnerTRef, phi_g_in);
            two_temp_pack_ghost(grid, grid.tt_inner_boundary1, kInnerRhoRes1, 0.0,
                                kInnerTRef, kInnerTRef, phi_g_in);
        } else {
            mixture_pack_ghost(grid, grid.mix_inner_boundary0, kInnerRhoRes0, 0.0,
                               kInnerTRef, phi_g_in);
            mixture_pack_ghost(grid, grid.mix_inner_boundary1, kInnerRhoRes1, 0.0,
                               kInnerTRef, phi_g_in);
        }
    }

    grid.broadcast();
}

// ---------------------------------------------------------------------------
// EXPERIMENTAL two-temperature leg — `model_column_2t`
// ---------------------------------------------------------------------------
// The same atmosphere, mesh, gravity, numerics and boundary geometry as the
// release column, advanced by the four-row (T_e != T_i) solver of src/two_temp/.
// Keeping ONE IC/BC implementation is deliberate: a two-temperature run that
// differed from the release in any configuration dimension other than the
// temperature split would not be a controlled experiment.
// ---------------------------------------------------------------------------

Vec model_column_2t_ic(Grid& grid) {
    kTwoTemperature = true;
    grid.two_temperature = true;
    // Outer-face heavy-particle conduction condition. Neumann (the default) is
    // the physically and characteristically motivated one; `dirichlet` is the
    // controlled comparison. See src/two_temp/conduction.cpp.
    grid.tt_outer_Ti_neumann = true;
    if (const char* choice = std::getenv("TT_OUTER_TI")) {
        const std::string value(choice);
        if (value == "neumann") {
            grid.tt_outer_Ti_neumann = true;
        } else if (value == "dirichlet") {
            grid.tt_outer_Ti_neumann = false;
        } else {
            throw std::invalid_argument(
                "TT_OUTER_TI must be neumann (default: no heavy conductive flux "
                "from the fully ionized plasma above the domain) or dirichlet "
                "(the controlled comparison, the electron wall imposed on T_i too)");
        }
    }
    // The two-temperature solver has exactly one numerical flux (the release
    // SWMF exact-Riemann Godunov flux on its three release rows) and exactly one
    // primitive set ((ln rho, u, ln p, ln T_e)). Rejecting the release's
    // reference-only overrides keeps a setting from being silently ignored.
    for (const char* knob : {"ISO_RIEMANN", "ISO_RECONSTRUCTION"}) {
        if (std::getenv(knob))
            throw std::runtime_error(
                std::string(knob) + " is a release-solver reference override and "
                "has no meaning for the experimental two-temperature solver, which "
                "has one flux and one primitive set");
    }
    Vec xn = model_column_ic(grid);
    std::cerr << "[model_column_2t] EXPERIMENTAL two-temperature solver: state="
                 "(rho, rho u, E, E_e)  outer_Te=dirichlet@"
              << grid.outer_conduction_temperature << "K  outer_Ti="
              << (grid.tt_outer_Ti_neumann ? "neumann" : "dirichlet")
              << "  inner=" << (grid.inner_conduction_neumann
                                ? "neumann(both)" : "dirichlet(both)")
              << std::endl;
    return xn;
}

void model_column_2t_update_bc(Grid& grid, const Vec& xn) {
    model_column_update_bc(grid, xn);
}

} // namespace chromosphere
