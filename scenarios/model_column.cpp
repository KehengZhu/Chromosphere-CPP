/*!
 * @file scenarios/model_column.cpp
 * @brief The field-aligned chromosphere-to-corona column: the `model_column`
 *        RELEASE scenario and the `model_gentle` historical two-fluid preset.
 * @ingroup scenarios
 *
 * One IC/BC implementation, two scenarios that never cross solvers. In release
 * mode this file loads the production `Gamma1` table unconditionally, builds the
 * Saha-HSE initial column from the C7 temperature profile, sets the validated
 * release numerics (log-space MUSCL with the MC3/Koren limiter, the mixture Roe
 * flux, the pressure-based reconstruction) and the 22,000 K external conductive
 * reservoir at the physical outer face — and rejects every historical two-fluid
 * knob. See scenarios/model_column.hpp for the full environment-variable list.
 */
#include "model_column.hpp"
#include "single_fluid/mixture.hpp"
#include "model_c7.hpp"   // c7_full_profile + c7_route_b_photoionization (C7 library)
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
float kAjump       = 1.0f;      // top temperature jump: T_ghost1 = a·T_top
float kBjump       = 1.0f;      // T_ghost2 = b·T_ghost1
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
float kInnerRhoI   = 0.0f;      // base ion density (RHO slot ghost 0)
float kInnerRhoN   = 0.0f;
// Two-fluid (non-release) ghost ladder: first-order p_0 + k ρ_0 g Δs with ρ ∝ p.
float kInnerPiGh   = 0.0f;      // ion ghost pressure   p_i,0 + ρ_i,0 g Δs
float kInnerPnGh   = 0.0f;      // neutral ghost pressure p_n,0 + ρ_n,0 g Δs
float kInnerPiGh2  = 0.0f;      // second inner ghost (2 Δs below): p_0 + 2ρgΔs
float kInnerPnGh2  = 0.0f;
float kInnerRhoIGh  = 0.0f;     // isothermal-hydrostatic ghost densities (ρ ∝ p at T₀):
float kInnerRhoNGh  = 0.0f;     //   ρ_{i,n},G0 = ρ_{i,n},0 · p_{i,n},G0 / p_{i,n},0
float kInnerRhoIGh2 = 0.0f;     //   ρ_{i,n},G1 = ρ_{i,n},0 · p_{i,n},G1 / p_{i,n},0
float kInnerRhoNGh2 = 0.0f;
float kInnerTRef     = 0.0f;     // reservoir temperature at the truncation
// RELEASE lower reservoir, captured once at the IC base. Unlike the two-fluid
// ladder above, each rung is the trapezoidal, EOS-closed hydrostatic step
// p_{k+1} = p_k + Δs·½(ρ_k + ρ_{k+1})·g with ρ_{k+1} = ρ_EOS(p_{k+1}, T₀),
// iterated to a fixed point — the SAME construction the outer face already used
// (hse_rung), so the two boundaries are now discretely consistent. The first-
// order form this replaces carried an O((Δs/H)²) pressure error, i.e. an O(Δs/H)
// error in the hydrostatic force balance of the lowest interior cells.
float kInnerPRes0   = 0.0f;     // reservoir pressure at ghost 0 (1 Δs below cell 0)
float kInnerPRes1   = 0.0f;     // reservoir pressure at ghost 1 (2 Δs below)
float kInnerRhoRes0 = 0.0f;     // EOS density at (kInnerPRes0, kInnerTRef)
float kInnerRhoRes1 = 0.0f;
float kTtopRef     = 0.0f;      // FIXED reference top temperature [K] for the jump
bool  kTtopFixedAbs = false;    // ISO_T_TOP>0: kTtopRef is an ABSOLUTE imposed top T
// FIXED reference outer ghost pressure [Pa] — the coronal reservoir back-pressure,
// captured one-sided-hydrostatically from the IC top cell (and recaptured at the
// relax→Stage-2 edge). Subsonic outflow admits exactly ONE incoming characteristic,
// so exactly one condition may come from outside: this back-pressure. Anchoring it to
// a FIXED reference rather than to a live interior cell is what makes the outer face
// well posed — see the outer-face block for why both live-cell anchors fail.
float kOuterPRef    = 0.0f;
bool  kOuterPRefSet = false;    // captured lazily on the first update_bc call
float kHeatFluxOn  = false;     // Stage 2 (conduction + jump) ⇒ open (Mach-capped) top
float kVcapMach    = 0.1f;      // outer outflow cap as a fraction of c_s
bool  kInnerTNeumann = false;   // lower-BC conduction temperature: Dirichlet vs Neumann
// Relaxation phase (ISO_RELAX_TIME > 0): Stage-1 adiabatic relaxation (heat flux off,
// no jump, V=0 reservoir top) for sim_time < kRelaxTime, THEN Stage 2.
float kRelaxTime     = 0.0f;
bool  kHeatFluxTarget = false;
bool  kStage2Started  = true;   // false during relax; flips to recapture kTtopRef
// Corona testbed (ISO_CORONA=1): resolved ~1 MK corona (needs Spitzer κ_e + TRAC).
bool  kCorona = false;
// Imposed Neumann coronal flux q(T) ramp (ISO_QFLUX*), the model_c7/RTV closure.
bool  kQflux = false;
// ISO_HYDRO_T_DECOUPLE=1: split the two roles the outer ghost temperature plays.
// The HYDRO ghost temperature becomes a zero-gradient extrapolation of the live top
// cell (T_hydro_g0 = T_hydro_g1 = T_top) instead of the fixed a·T_ref wall, while
// CONDUCTION keeps seeing exactly the same fixed wall a·T_ref through
// grid.outer_conduction_temperature. Subsonic outflow admits ONE incoming
// characteristic and the fixed reservoir back-pressure kOuterPRef already supplies
// it; the additional fixed hydro T pins the thermal/entropy state as well, which is
// the suspected source of the ~20-30 km top boundary layer. This flag removes ONLY
// that extra hydro constraint, holding the conductive energy input fixed.
// Default 0 ⇒ unchanged behaviour.
bool  kHydroTDecouple = false;

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

// Pack a LEGACY two-fluid ghost cell from explicit ion/neutral densities,
// temperature, velocity and the gravitational potential (p_i = 2 n_i k T incl.
// electrons, p_n = n_n k T, ε_e = 3/2 p_e ⇒ T_e = T_i). Taking (rho_i, rho_n)
// explicitly lets each ghost carry the LOCAL ionization fraction (ionized corona
// ↔ neutral chromosphere). The release solver packs its ghosts through
// mixture_pack_ghost instead, which needs only the total density.
void pack_ghost(Vec& ob, const Grid& grid, float rho_i, float rho_n, float T, float V, float phi_g) {
    const float n_i = rho_i / grid.m_i;
    const float n_n = rho_n / grid.m_n;
    const float p_i = 2.0f * n_i * grid.k_b * T;   // protons + electrons
    const float p_n =        n_n * grid.k_b * T;
    const float p_e =        n_i * grid.k_b * T;   // ε_e = 3/2 p_e ⇒ T_e = T
    ob(cons::RHO_I) = rho_i;
    ob(cons::RHO_N) = rho_n;
    ob(cons::MOM_I) = rho_i * V;
    ob(cons::MOM_N) = rho_n * V;
    ob(cons::E_I)   = grid.inv_gm1() * p_i + 0.5f * rho_i * V * V + rho_i * phi_g;
    ob(cons::E_N)   = grid.inv_gm1() * p_n + 0.5f * rho_n * V * V + rho_n * phi_g;
    ob(cons::E_E)   = grid.inv_gm1() * p_e;
}

// Per-cell (T, n_i, n_n) at height h_km from the real C7 atmosphere via
// c7_full_profile (Avrett & Loeser Table 26): n_i = n_e (ionized), n_n = n_HI.
void c7_cell(float h_km, float& T, float& n_i, float& n_n) {
    float n_e, n_HI;
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
    auto env_f = [](const char* key, float fallback) -> float {
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
    // A loaded Gamma1 table means the single-fluid RELEASE solver; an empty one
    // means the historical two-fluid research solver (scenario `model_gentle`).
    // `model_column` always supplies the table, so it always lands here as a
    // release run.
    const bool gamma_mode = !grid.eos_gamma_table.empty();

    // Release configuration surface. Everything below is historical two-fluid
    // research physics that the single-fluid release timestep does not contain:
    // a fixed adiabatic index, a neutral fluid, finite-rate ionization, a
    // radiative sink, TRAC broadening, volumetric coronal heating, an imposed
    // Neumann coronal flux, an IC coronal superheat, and the mesh-scaled
    // artificial conduction term. Rejecting them outright is what keeps
    // `model_column` a single unambiguous production identity — a release run
    // can never be silently reconfigured into something else.
    if (gamma_mode) {
        static const char* const kNonReleaseKnobs[] = {
            "ISO_GAMMA", "ISO_TWO_FLUID", "ISO_IONIZATION", "ISO_COOLING",
            "ISO_TRAC", "ISO_CORONA", "ISO_CHEAT", "ISO_QFLUX", "ISO_TBOOST",
            "ISO_NUMERICAL_DIFFUSIVITY_MULT"};
        for (const char* knob : kNonReleaseKnobs) {
            if (std::getenv(knob))
                throw std::runtime_error(
                    std::string(knob) + " is historical two-fluid research "
                    "configuration and has no meaning for the single-fluid "
                    "release solver; run the model_gentle scenario for that path");
        }
    }

    // --- IC mode ----------------------------------------------------------
    // The IC is ALWAYS the real Model C7 profile with the density re-integrated
    // hydrostatically (clean V≈0 start; the analytic-isentrope toy IC is retired).
    // ISO_CORONA extends the top up into a resolved ~1 MK corona (to
    // ISO_CORONA_TOP_KM) and turns on TRAC so the ionized corona conducts via
    // Spitzer κ_e and can drive A&S78 conduction-driven evaporation.
    kCorona = (env_f("ISO_CORONA", 0.0f) != 0.0f);

    // --- parameters -------------------------------------------------------
    // Lower truncation height of the modeled domain, in km measured from the
    // photosphere (that is the HEIGHT GAUGE only -- the domain itself starts
    // here, at 1600 km, well above it). See the inner-boundary block.
    const float h_base     = env_f("ISO_H_BASE", 0.0f);
    const float corona_top = env_f("ISO_CORONA_TOP_KM", 10000.0f);          // km, coronal top
    float       DH_km      = env_f("ISO_DH", kCorona ? (corona_top - h_base) : 1303.0f); // km
    // Base (T, n_H) from the C7 profile at h_base.
    float T_base, ni_b, nn_b;
    c7_cell(h_base, T_base, ni_b, nn_b);
    if (gamma_mode) T_base = static_cast<float>(c7_full_temperature_pchip(h_base));
    const bool  heat_flux_on = (env_f("ISO_HEAT_FLUX", 0.0f) != 0.0f);
    // Stage 2b radiative sink. Gentle evaporation is the competition between the
    // downward conductive flux and radiative cooling (Antiochos & Sturrock 1978).
    const bool  cooling_on = heat_flux_on && (env_f("ISO_COOLING", 0.0f) != 0.0f);
    // LEGACY two-fluid research opt-ins (model_gentle). They are read only when
    // no Gamma1 table is loaded: the release solver is a single-fluid equilibrium
    // mixture and has neither a neutral fluid nor a finite-rate ionization stage.
    const bool  two_fluid_on  = !gamma_mode && (env_f("ISO_TWO_FLUID", 0.0f) != 0.0f);
    const bool  ionization_on = !gamma_mode && (env_f("ISO_IONIZATION", 0.0f) != 0.0f);
    kHeatFluxOn        = heat_flux_on;
    kAjump             = env_f("ISO_TJUMP_A", 1.0f);
    kBjump             = env_f("ISO_TJUMP_B", 1.0f);
    kVcapMach          = env_f("ISO_VCAP", 0.1f);
    kInnerTNeumann     = (env_f("ISO_INNER_T_NEUMANN", 0.0f) != 0.0f);
    // Upper-BC over-specification experiment: hydro ghost T free-floats (zero
    // gradient), conduction wall unchanged. Default 0 ⇒ baseline behaviour.
    kHydroTDecouple    = (env_f("ISO_HYDRO_T_DECOUPLE", 0.0f) != 0.0f);

    // Adiabatic index. Default 5/3. ISO_GAMMA near 1 (e.g. 1.05) makes the gas
    // nearly isothermal — a polytropic stand-in for the radiative thermostat the
    // conduction-only experiment leaves out. Set BEFORE the IC/BC build.
    if (!gamma_mode) grid.gamma_mono = env_f("ISO_GAMMA", grid.gamma_mono);
    grid.single_fluid = gamma_mode || !two_fluid_on;
    grid.enable_ionization = ionization_on;
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
    const float ds_m = (DH_km * 1000.0f) / static_cast<float>(ns_coarse);   // coarse Δs [m]
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
            h_F(k) = h_base + static_cast<float>(face_m[k]) * 1.0e-3f;   // m -> km
        for (arma::uword i = 0; i < grid.ns; ++i) {
            grid.ds_i(i)      = static_cast<float>(face_m[i + 1] - face_m[i]);
            grid.phi_g_imh(i) = grid.g * (h_F(i)     - h_base) * 1000.0f;   // J/kg
            grid.phi_g_iph(i) = grid.g * (h_F(i + 1) - h_base) * 1000.0f;
        }
    } else {
        for (arma::uword k = 0; k <= grid.ns; ++k)
            h_F(k) = h_base + DH_km * static_cast<float>(k) / static_cast<float>(grid.ns);
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
    double ln_p = 0.0, prev_h = h_base, prev_invH = 0.0;
    {
        const double nbase = static_cast<double>(ni_b)+nn_b;
        const double xb = gamma_mode ? saha_ionization_fraction_n_h(nbase, T_base)
                                     : ni_b/nbase;
        ln_p = std::log((1.0+xb)*nbase*grid.k_b*T_base);
        prev_invH = grid.g*grid.m_i/((1.0+xb)*grid.k_b*T_base);
    }
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const double h_c = 0.5 * (h_F(i) + h_F(i + 1));        // km
        float T, n_i, n_n;
        c7_cell(static_cast<float>(h_c), T, n_i, n_n);
        if (gamma_mode) {
            T = static_cast<float>(c7_temperature_cell_average(h_F(i), h_F(i+1)));
        }
        const double dh_m = (h_c-prev_h)*1000.0;
        double x = n_i/(n_i+n_n);
        if (gamma_mode) {
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
        } else {
            const double invH = grid.g*grid.m_i/((1.0+x)*grid.k_b*T);
            ln_p -= 0.5*(prev_invH+invH)*dh_m;
        }
        const double n_tot = std::exp(ln_p)/((1.0+x)*grid.k_b*T);
        if (gamma_mode) x = saha_ionization_fraction_n_h(n_tot, T);
        n_i = x * n_tot;
        n_n = (1.0-x)*n_tot;
        prev_h = h_c;
        prev_invH = grid.g*grid.m_i/((1.0+x)*grid.k_b*T);
        T_c(i)  = T;
        ni_c(i) = n_i;
        nn_c(i) = n_n;
    }

    Vec xn = arma::zeros<Vec>(
        gamma_mode ? grid.n_mixture_state : grid.n_state);
    const auto sz = arma::size(grid.ns, num_of_eq);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float T     = T_c(i);
        const float n_i   = ni_c(i);
        const float n_n   = nn_c(i);
        const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        const float rho_i = n_i * grid.m_i;
        const float rho_n = n_n * grid.m_n;
        if (gamma_mode) {
            mixture_pack_cell(grid, xn, i, static_cast<double>(rho_i)+rho_n,
                              0.0, T, phi_g);
            continue;
        }
        const float p_i   = 2.0f * n_i * grid.k_b * T;
        const float p_n   =        n_n * grid.k_b * T;
        const float p_e   =        n_i * grid.k_b * T;
        xn(arma::sub2ind(sz, i, cons::RHO_I)) = rho_i;
        xn(arma::sub2ind(sz, i, cons::RHO_N)) = rho_n;
        xn(arma::sub2ind(sz, i, cons::MOM_I)) = 0.0f;
        xn(arma::sub2ind(sz, i, cons::MOM_N)) = 0.0f;
        xn(arma::sub2ind(sz, i, cons::E_I))   = grid.inv_gm1() * p_i + rho_i * phi_g;
        xn(arma::sub2ind(sz, i, cons::E_N))   = grid.inv_gm1() * p_n + rho_n * phi_g;
        xn(arma::sub2ind(sz, i, cons::E_E))   = grid.inv_gm1() * p_e;
    }

    // Reference top temperature for the imposed ghost jump (Stage 2), FIXED at the
    // IC top value: the unresolved TR/corona above sits at ~constant T, so
    // T_ghost = a·T_ref is a fixed hot reservoir and the downward conductive flux
    // q ∝ (T_ghost − T_top) SELF-LIMITS as the chromosphere heats.
    kTtopRef = T_c(grid.ns - 1);
    // ISO_T_TOP>0: override with an ABSOLUTE fixed outer-boundary temperature [K]
    // (e.g. 22000 = C7's TR-base T). Not re-anchored across a relax phase.
    const float t_top_abs = env_f("ISO_T_TOP", 0.0f);
    kTtopFixedAbs = (t_top_abs > 0.0f);
    if (kTtopFixedAbs) kTtopRef = t_top_abs;

    // --- one-time IC coronal superheat (ISO_TBOOST > 1) -------------------
    // HISTORICAL two-fluid option (rejected in release mode).
    // Heat-flux-driven evaporation with NO heating ramp: raise the IC coronal
    // temperature above its own conduction+radiation balance so its own Spitzer
    // flux conducts down the TR and drives evaporation from t=0 (Antiochos &
    // Sturrock 1978). Switched on smoothly ABOVE the TR (tanh in T) so the
    // chromosphere/TR base is untouched; density held fixed (heat raises p).
    const float T_boost = env_f("ISO_TBOOST", 1.0f);
    if (T_boost > 1.0f) {
        const float Tc = env_f("ISO_TBOOST_TC", 3.0e5f);   // selector center [K]
        const float Tw = env_f("ISO_TBOOST_W",  1.0e5f);   // selector width  [K]
        for (arma::uword i = 0; i < grid.ns; ++i) {
            const float rho_i = ni_c(i) * grid.m_i;
            const float rho_n = nn_c(i) * grid.m_n;
            const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
            const float n_i = rho_i / grid.m_i;
            const float n_n = rho_n / grid.m_n;
            const float T   = T_c(i);
            const float w   = 0.5f * (1.0f + std::tanh((T - Tc) / Tw));   // ~0 chromo → ~1 corona
            const float Tb  = (1.0f + (T_boost - 1.0f) * w) * T;
            if (gamma_mode) {
                mixture_pack_cell(grid, xn, i, static_cast<double>(rho_i)+rho_n,
                                  0.0, Tb, phi_g);
                continue;
            }
            const float p_i_b = 2.0f * n_i * grid.k_b * Tb;               // protons + electrons
            const float p_n_b =        n_n * grid.k_b * Tb;
            xn(arma::sub2ind(sz, i, cons::E_I)) = grid.inv_gm1() * p_i_b + rho_i * phi_g;
            xn(arma::sub2ind(sz, i, cons::E_N)) = grid.inv_gm1() * p_n_b + rho_n * phi_g;
            xn(arma::sub2ind(sz, i, cons::E_E)) = grid.inv_gm1() * n_i * grid.k_b * Tb;   // T_e = Tb
        }
    }

    // --- lower truncation reservoir (h_base = 1600 km, NOT the photosphere) ---
    // The ghosts sit one and two Δs below cell 0 and stand in for the unmodeled
    // lower chromosphere, which on these timescales is a quasi-static stratified
    // reservoir. The RELEASE ladder (gamma_mode) is trapezoidal and EOS-closed at
    // the truncation temperature T₀; the legacy two-fluid ladder below it keeps
    // the historical first-order form p_0 + kρ_0gΔs with ρ ∝ p.
    {
        const float T0    = T_c(0);
        kInnerTRef = T0;
        const float n_i0  = ni_c(0);
        const float n_n0  = nn_c(0);
        kInnerRhoI = n_i0 * grid.m_i;
        kInnerRhoN = n_n0 * grid.m_n;
        const float p_i0 = 2.0f * n_i0 * grid.k_b * T0;
        const float p_n0 =        n_n0 * grid.k_b * T0;
        // Center-to-ghost distance = the BASE cell's own width (mirrored ghost). On a
        // refined mesh cell 0 is a fine cell, so this is the fine Δs, NOT the coarse ds_m.
        const float ds_base = grid.ds_i(0);
        if (gamma_mode) {
            const double rho0 = kInnerRhoI+kInnerRhoN;
            const double p0 = (1.0+saha_ionization_fraction_n_h(
                rho0/eos_constants::m_h, T0))*rho0/eos_constants::m_h*eos_constants::k_b*T0;
            double pg0, rg0, pg1, rg1;
            gamma_hse_rung_down(grid, p0,  rho0, T0, ds_base, pg0, rg0);
            gamma_hse_rung_down(grid, pg0, rg0,  T0, ds_base, pg1, rg1);
            kInnerPRes0 = pg0; kInnerRhoRes0 = rg0;
            kInnerPRes1 = pg1; kInnerRhoRes1 = rg1;
            kInnerPiGh = pg0; kInnerPnGh = 0.0f;
            kInnerPiGh2 = pg1; kInnerPnGh2 = 0.0f;
            kInnerRhoIGh = rg0;
            kInnerRhoNGh = 0.0f;
            kInnerRhoIGh2 = rg1;
            kInnerRhoNGh2 = 0.0f;
        } else {
        kInnerPiGh  = p_i0 + kInnerRhoI * grid.g * ds_base;
        kInnerPnGh  = p_n0 + kInnerRhoN * grid.g * ds_base;
        kInnerPiGh2 = p_i0 + 2.0f * kInnerRhoI * grid.g * ds_base;
        kInnerPnGh2 = p_n0 + 2.0f * kInnerRhoN * grid.g * ds_base;
        // ρ ∝ p at fixed T₀ (isothermal-hydrostatic ghost densities).
        kInnerRhoIGh  = kInnerRhoI * (kInnerPiGh  / p_i0);
        kInnerRhoIGh2 = kInnerRhoI * (kInnerPiGh2 / p_i0);
        kInnerRhoNGh  = kInnerRhoN * (kInnerPnGh  / p_n0);
        kInnerRhoNGh2 = kInnerRhoN * (kInnerPnGh2 / p_n0);
        }
    }

    // --- numerical scheme: the validated "bestwb" configuration (always on) ---
    // well-balanced + log-space MUSCL + MC3/Koren (β=2) + equilibrium-reference
    // δ-form well-balancing. These hold a hydrostatic column at V≈0 to round-off
    // for any stratification (docs/gentle_evaporation_downflow.md). No toggles.
    grid.out_base_km              = h_base;   // label output heights from the true base
    grid.inner_conduction_neumann = kInnerTNeumann;
    grid.well_balanced            = true;
    grid.log_reconstruct          = true;
    grid.mc3_limiter              = true;
    grid.limiter_beta             = 2.0f;
    // DIAGNOSTIC-ONLY reconstruction override for the scheme-comparison leg of the
    // face-flux study (docs/top_ripple_face_flux_diagnosis.md). Unset — the default
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
    // The released Gamma/Saha model_column advances the corrector with the 3x3
    // mixture Roe characteristic flux and limits (log rho, V, log p), recovering
    // the face temperature by inverting the same authoritative Saha closure. Both
    // require the equilibrium-manifold face builder, so they are enabled only in
    // gamma_mode; the legacy non-gamma path keeps Rusanov + (log rho, V, log T).
    grid.roe_characteristic_flux = gamma_mode;
    grid.pressure_reconstruct    = gamma_mode;
    // Reference/debug overrides. ISO_RIEMANN=rusanov restores the local
    // Lax-Friedrichs flux and ISO_RECONSTRUCTION=lnrho-v-lnt the temperature-based
    // primitive set; both are retained for regression and controlled numerical
    // comparison only. Neither is needed for a production run.
    if (const char* riemann = std::getenv("ISO_RIEMANN")) {
        const std::string choice(riemann);
        if (choice == "rusanov") {
            grid.roe_characteristic_flux = false;
        } else if (choice == "roe-local") {
            if (!gamma_mode)
                throw std::invalid_argument("ISO_RIEMANN=roe-local requires Gamma/Saha mode");
            grid.roe_characteristic_flux = true;
        } else {
            throw std::invalid_argument("ISO_RIEMANN must be rusanov or roe-local");
        }
        std::cerr << "[model_column] Riemann override: ISO_RIEMANN=" << choice
                  << " (release default: roe-local)" << std::endl;
    }
    if (const char* recon = std::getenv("ISO_RECONSTRUCTION")) {
        const std::string choice(recon);
        if (choice == "lnrho-v-lnt") {
            grid.pressure_reconstruct = false;
        } else if (choice == "lnrho-v-lnp") {
            if (!gamma_mode)
                throw std::invalid_argument(
                    "ISO_RECONSTRUCTION=lnrho-v-lnp requires Gamma/Saha mode");
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
              << (grid.roe_characteristic_flux ? "roe-local" : "rusanov")
              << " reconstruction="
              << (grid.pressure_reconstruct ? "lnrho-v-lnp" : "lnrho-v-lnt")
              << " limiter=" << (grid.mc3_limiter ? "mc3" : "minmod")
              << " beta=" << grid.limiter_beta << std::endl;

    // --- runtime physics toggles ------------------------------------------
    // The RELEASE branch is deliberately a straight line: hydro plus physical
    // conduction, nothing else. The optional stages below it exist only for the
    // historical two-fluid research path, and the release rejects their env
    // knobs above, so none of them can be reached from a model_column run.
    grid.single_fluid             = gamma_mode || !two_fluid_on;
    grid.enable_ionization        = ionization_on;
    grid.enable_radiative_cooling = cooling_on;
    grid.enable_beam_heating      = false;
    grid.enable_coronal_heating   = false;
    grid.impose_outer_heat_flux   = false;   // default: heat enters via the ghost-T jump
    grid.enable_Te                = false;    // two-fluid (ion/neutral), NOT three-temperature
    if (gamma_mode) {
        grid.enable_trac         = false;    // no TRAC stage in the release solver
        grid.enable_vacuum_floor = false;
    } else {
        // TRAC broadens the unresolved TR. ISO_TRAC force-enables it independent of
        // cooling/corona (auto = cooling || corona). ISO_TRAC_TCHROM sets the region base
        // T_b; ISO_TRAC_TCMAXFRAC lifts Johnston's 0.2·T_peak cap (needed on a capped
        // low-T domain where 0.2·T_peak < T_b would pin T_c to the floor).
        const float iso_trac = env_f("ISO_TRAC", -1.0f);   // -1 = auto (cooling||corona)
        grid.enable_trac         = (iso_trac >= 0.0f) ? (iso_trac != 0.0f)
                                                      : (cooling_on || kCorona);
        grid.trac_T_chrom        = env_f("ISO_TRAC_TCHROM", 2.0e4f);
        grid.trac_Tc_max_frac    = env_f("ISO_TRAC_TCMAXFRAC", 0.2f);
        grid.trac_cutoff_T       = grid.trac_T_chrom;
        grid.enable_vacuum_floor = cooling_on || kCorona || ionization_on;
    }

    // --- Route-B photoionization closure (only when the network is on) ----
    // Invert the local ionization-equilibrium balance of the active Stage-E network
    // on the IC so the C7 chromosphere is a fixed point (shared C7-library helper).
    if (ionization_on) {
        Vec h_cell(grid.ns);
        for (arma::uword i = 0; i < grid.ns; ++i)
            h_cell(i) = 0.5f * (h_F(i) + h_F(i + 1));
        grid.photoionization_rate_i = c7_route_b_photoionization(grid, T_c, ni_c, nn_c, h_cell);
    }

    // --- imposed Neumann coronal flux q(T) + ramp (ISO_QFLUX*) ------------
    // HISTORICAL two-fluid option (rejected in release mode).
    // Alternative to the ghost-T jump: impose the downward coronal conductive flux
    // as a Stage-D Neumann outer BC (RTV 1978 / model_c7 closure), optionally ramped
    // up to drive gentle evaporation (Antiochos & Sturrock 1978). The base is C7's
    // own top-gradient flux q = κ₀ T^{5/2} dT/ds unless ISO_QFLUX0 overrides it.
    kQflux = (env_f("ISO_QFLUX", 0.0f) != 0.0f);
    if (kQflux) {
        const float kappa0 = 1.0e-11f;   // W m^-1 K^-7/2
        const arma::uword nl = grid.ns - 1;
        const float T_top = T_c(nl);
        const float dTds  = (T_c(nl) - T_c(nl - 1))
                          / (0.5f * (grid.ds_i(nl) + grid.ds_i(nl - 1)));   // K/m
        const float q_c7  = kappa0 * std::pow(T_top, 2.5f) * dTds;
        const float q_base = env_f("ISO_QFLUX0", q_c7);
        grid.impose_outer_heat_flux   = true;
        grid.outer_heat_flux_base     = q_base;
        grid.outer_heat_flux          = q_base;
        grid.outer_heat_flux_enhance  = env_f("ISO_QFLUX_ENHANCE", 1.0f);
        grid.outer_heat_flux_t_on     = env_f("ISO_QFLUX_T_ON",  600.0f);
        grid.outer_heat_flux_ramp     = env_f("ISO_QFLUX_RAMP",  120.0f);
    }

    // --- ambient volumetric coronal heating H(s) (ISO_CHEAT=1) ------------
    // HISTORICAL two-fluid option (rejected in release mode).
    // The RTV (1978) static-balance term the resolved corona otherwise lacks:
    //   H(s) = E_H0 · exp(−(s−s0)/s_H)  [W/m³]  (Aschwanden & Schrijver 2002),
    // deposited through the corona so it balances the radiative loss locally and
    // pins V≈0; ramping it up (ISO_CHEAT_ENHANCE>1, slow) drives gentle evaporation.
    // Calibration inverts the RTV/AS2002 static scaling for the E_H0 that relaxes a
    // column of length L to apex temperature T_max. When on, the boundary q(T) is
    // turned OFF (the corona is heated inside the domain; no double-heating).
    if (env_f("ISO_CHEAT", 0.0f) != 0.0f) {
        const float L    = arma::sum(grid.ds_i);                 // column length [m]
        const float Tmax = env_f("ISO_CHEAT_TMAX_MK", 1.5f) * 1.0e6f;
        const float sHf  = env_f("ISO_CHEAT_SH_FRAC", 0.4f);
        const float s_H  = sHf * L;
        const float L_cm = L * 100.0f;
        const float p0   = std::pow(Tmax / 1400.0f, 3.0f) / L_cm;
        const float Eun  = 9.8e4f * std::pow(p0, 7.0f / 6.0f) * std::pow(L_cm, -5.0f / 6.0f);
        const float E0   = 0.1f * Eun * std::exp(0.5f * L / s_H);   // erg/cm³/s → W/m³

        grid.enable_coronal_heating = true;
        grid.coronal_heat_E0        = env_f("ISO_CHEAT_E0", E0);
        grid.coronal_heat_sH        = s_H;
        grid.coronal_heat_s0        = env_f("ISO_CHEAT_S0_MM", 0.0f) * 1.0e3f;
        grid.coronal_heat_two_sided = false;                       // open column, base-anchored
        grid.coronal_heat_enhance   = env_f("ISO_CHEAT_ENHANCE", 1.0f);
        grid.coronal_heat_t_on      = env_f("ISO_CHEAT_T_ON",  0.0f);
        grid.coronal_heat_ramp      = env_f("ISO_CHEAT_RAMP", 30.0f);
        // Corona heated internally ⇒ drop the boundary q(T) crutch.
        grid.impose_outer_heat_flux = false;
        grid.outer_heat_flux_base   = 0.0f;
        grid.outer_heat_flux        = 0.0f;
        kQflux = false;
    }

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
    // by the top T jump (or the imposed q(T) flux / coronal heating).
    grid.enable_conduction        = relaxing0 ? false : heat_flux_on;
    // The release conduction operator is structurally physical-only and never
    // reads this field. The historical two-fluid Stage D still folds a
    // mesh-scaled artificial diffusivity into its conduction operator; keep it
    // there behind an explicit diagnostic multiplier for controlled legacy
    // comparisons, but never enable it implicitly.
    // ISO_NUMERICAL_DIFFUSIVITY_MULT=1 restores the historical C_num=2000 m/s
    // (4x in the resolved-corona testbed); unset/0 is the two-fluid baseline.
    if (gamma_mode) {
        grid.numerical_diffusivity_per_length = 0.0f;
    } else {
        const float diff_mult = (kCorona ? 4.0f : 1.0f)
                              * env_f("ISO_NUMERICAL_DIFFUSIVITY_MULT", 0.0f);
        grid.numerical_diffusivity_per_length = heat_flux_on
            ? (diff_mult * 2.0e3f) : 0.0f;
    }

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
    // hydrostatic face mass-flux defect of the release mesh by ~560x, to within
    // ~1.3x of the float32 round-off floor of the stored state, so the frozen
    // reference no longer has anything material to remove.
    //
    // The non-gamma (two-fluid) leg of this scenario is a different solver that
    // never received the predictor source term, so it keeps the frozen
    // reference exactly as before. eq_wb is now a two-fluid-only mechanism.
    grid.eq_wb = !gamma_mode;
    if (grid.eq_wb) grid.eq_state = xn;
    return xn;
}

void model_column_update_bc(Grid& grid, const Vec& xn) {
    const auto  sz  = arma::size(grid.ns, num_of_eq);
    const float m_i = grid.m_i;
    const float m_n = grid.m_n;
    const float k_b = grid.k_b;
    const float g   = grid.g;

    // Decode total density / pressure / temperature / velocity / ionization fraction
    // x = ρ_i/ρ_tot of an interior cell (cons2prim convention). Single-fluid T uses
    // the LOCAL ionization (p = (2 n_i + n_n) k T).
    // Release state: decode the mixture cell straight out of (rho, rho u, E).
    auto mixture_cell_state = [&](arma::uword i, float& rho_tot, float& p_tot,
                                  float& T, float& V, float& x) {
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

    // Legacy two-fluid state (single-fluid T from the local ionization).
    auto two_fluid_cell_state = [&](arma::uword i, float& rho_tot, float& p_tot,
                                    float& T, float& V, float& x) {
        const float rho_i = xn(arma::sub2ind(sz, i, cons::RHO_I));
        const float rho_n = xn(arma::sub2ind(sz, i, cons::RHO_N));
        const float momI  = xn(arma::sub2ind(sz, i, cons::MOM_I));
        const float momN  = xn(arma::sub2ind(sz, i, cons::MOM_N));
        const float E_i   = xn(arma::sub2ind(sz, i, cons::E_I));
        const float E_n   = xn(arma::sub2ind(sz, i, cons::E_N));
        const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        const float Vi    = momI / rho_i;
        const float Un    = momN / rho_n;
        const float p_i   = grid.gm1() * E_i - grid.half_gm1() * rho_i * Vi * Vi - grid.gm1() * rho_i * phi_g;
        const float p_n   = grid.gm1() * E_n - grid.half_gm1() * rho_n * Un * Un - grid.gm1() * rho_n * phi_g;
        const float n_i   = rho_i / m_i;
        const float n_n   = rho_n / m_n;
        rho_tot = rho_i + rho_n;
        p_tot   = p_i + p_n;
        T       = p_tot / ((2.0f * n_i + n_n) * k_b);   // single-fluid T, local ionization
        V       = Vi;
        x       = rho_i / rho_tot;                      // ionization fraction
    };

    const bool release = !grid.eos_gamma_table.empty();
    auto cell_state = [&](arma::uword i, float& rho_tot, float& p_tot,
                          float& T, float& V, float& x) {
        if (release) mixture_cell_state(i, rho_tot, p_tot, T, V, x);
        else         two_fluid_cell_state(i, rho_tot, p_tot, T, V, x);
    };

    // --- imposed q(T) flux ramp (ISO_QFLUX) --------------------------------
    // q(T) = q_base · enhance(t): holds at 1 until t_on, then cosine-rises to
    // outer_heat_flux_enhance over [t_on, t_on+ramp]. Applied to the Stage-D
    // Neumann outer flux (impose_outer_heat_flux). No-op when kQflux is off.
    if (kQflux && grid.outer_heat_flux_base > 0.0f) {
        float amp = 1.0f;
        if (grid.outer_heat_flux_enhance != 1.0f) {
            const float t   = grid.sim_time;
            const float ton = grid.outer_heat_flux_t_on;
            const float r   = (grid.outer_heat_flux_ramp > 1.0e-6f) ? grid.outer_heat_flux_ramp : 1.0e-6f;
            float w;
            if (t <= ton)          w = 0.0f;
            else if (t >= ton + r) w = 1.0f;
            else w = 0.5f * (1.0f - std::cos(static_cast<float>(arma::datum::pi) * (t - ton) / r));
            amp = 1.0f + (grid.outer_heat_flux_enhance - 1.0f) * w;
        }
        grid.outer_heat_flux = grid.outer_heat_flux_base * amp;
    }

    // --- relaxation → Stage-2 phase switch (off grid.sim_time) -------------
    // While relaxing: conduction OFF, no jump, V=0 reservoir top. At the transition
    // recapture kTtopRef from the RELAXED top cell (anchor the FIXED jump to the
    // equilibrium), then run Stage 2 (conduction ON + jump + Mach-capped outflow).
    const bool relaxing = (kRelaxTime > 0.0f) && (grid.sim_time < kRelaxTime);
    kHeatFluxOn            = relaxing ? false : kHeatFluxTarget;
    grid.enable_conduction = relaxing ? false : kHeatFluxTarget;
    const float a_live = relaxing ? 1.0f : kAjump;
    const float b_live = relaxing ? 1.0f : kBjump;
    if (!relaxing && !kStage2Started) {
        if (!kTtopFixedAbs) {           // an ABSOLUTE imposed top T (ISO_T_TOP) stays put
            float rho_t, p_t, T_t, V_t, x_t;
            cell_state(grid.ns - 1, rho_t, p_t, T_t, V_t, x_t);
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
        const float ds       = grid.ds_i(nl);
        const float phi_g_out = grid.phi_g_iph(nl);   // ghost decode uses this gauge

        float rho_top, p_top, T_top, V_top, x_top;
        cell_state(nl, rho_top, p_top, T_top, V_top, x_top);

        // (2) temperature. TWO distinct roles, kept explicitly separate:
        //   * T_cond_wall — the fixed hot wall that DRIVES Stage-D conduction, always
        //     the imposed jump a·T_ref off the reference top T. With ISO_CORONA the
        //     corona is resolved and a=b=1 (the ghost just extends it).
        //   * T_hydro_g0/g1 — what the HYDRO ghost carries (its EOS density and
        //     internal energy, hence the outer-face reconstruction and Riemann state).
        // Baseline (default): the hydro ghost IS the wall, T_hydro = a·T_ref, b·a·T_ref
        // — one ghost state serving both roles, so the outer face fixes pressure AND
        // temperature. ISO_HYDRO_T_DECOUPLE=1: the hydro ghost instead zero-gradient-
        // extrapolates the live top cell (T_hydro_g0 = T_hydro_g1 = T_top) and only the
        // conduction rows keep the wall, via grid.outer_conduction_temperature.
        const float T_cond_wall = a_live * kTtopRef;
        const float T_g0 = kHydroTDecouple ? T_top : T_cond_wall;
        const float T_g1 = kHydroTDecouple ? T_g0  : b_live * T_g0;
        grid.outer_conduction_temperature_override = kHydroTDecouple;
        grid.outer_conduction_temperature          = T_cond_wall;

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
        auto ghost_density = [&](float p_gh, float T_gh) {
            return release ? gamma_density_from_pressure(grid, p_gh, T_gh)
                           : p_gh*m_i/((1.0f+x_top)*k_b*T_gh);
        };
        // ρ_ghost depends on p_ghost, so take one fixed-point pass seeded with the
        // rung below (the correction is O(Δs·Δρ/ρ) and converges in a single sweep).
        auto hse_rung = [&](float p_below, float rho_below, float T_gh,
                            float& p_out, float& rho_out) {
            float rho_gh = rho_below;
            for (int it = 0; it < 2; ++it) {
                p_out = p_below - ds * 0.5f * (rho_below + rho_gh) * g;
                if (p_out < 1.0e-12f) p_out = 1.0e-12f;
                rho_gh = ghost_density(p_out, T_gh);
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
            float p_ref, rho_ref;
            hse_rung(p_top, rho_top, T_g0, p_ref, rho_ref);
            kOuterPRef    = p_ref;
            kOuterPRefSet = true;
            // The reservoir back-pressure is the one externally imposed condition at
            // the outer face; log it (once per capture) so post-run diagnostics can
            // form the top-cell pressure deficit p_top − p_g0 without re-deriving it.
            std::cerr << "[model_column] outer BC anchor: kOuterPRef=" << kOuterPRef
                      << " Pa  p_top=" << p_top << " Pa  T_top=" << T_top
                      << " K  T_hydro_g0=" << T_g0 << " K  T_cond_wall="
                      << T_cond_wall << " K  hydro_T_decouple="
                      << (kHydroTDecouple ? 1 : 0) << std::endl;
        }
        float p_g0 = kOuterPRef;
        if (p_g0 < 1.0e-12f) p_g0 = 1.0e-12f;
        const float rho_g0 = ghost_density(p_g0, T_g0);
        float p_g1, rho_g1;
        hse_rung(p_g0, rho_g0, T_g1, p_g1, rho_g1);

        // Velocity. Stage 1 (relaxation): static reservoir V=0. Stage 2: Mach-capped
        // OUTFLOW so the evaporation upflow can leave, while |V| ≤ kVcapMach·c_s
        // prevents the ill-posed-inflow runaway.
        float V_g = 0.0f;
        if (kHeatFluxOn) {
            const float c_s = release
                ? std::sqrt(gamma_state(grid.eos_gamma_table, rho_top, T_top,
                                        grid.eos_gamma_debug_clamp).gamma_sound*p_top/rho_top)
                : std::sqrt(2.0f*grid.gamma_mono*k_b*T_top/m_i);
            const float vcap = kVcapMach * c_s;
            V_g = V_top;
            if (V_g >  vcap) V_g =  vcap;
            if (V_g < -vcap) V_g = -vcap;
        }
        if (release) {
            mixture_pack_ghost(grid, grid.mix_outer_boundary0, rho_g0, V_g, T_g0,
                               phi_g_out);
            mixture_pack_ghost(grid, grid.mix_outer_boundary1, rho_g1, V_g, T_g1,
                               phi_g_out);
        } else {
            pack_ghost(grid.outer_boundary0_i, grid, x_top * rho_g0,
                       (1.0f - x_top) * rho_g0, T_g0, V_g, phi_g_out);
            pack_ghost(grid.outer_boundary1_i, grid, x_top * rho_g1,
                       (1.0f - x_top) * rho_g1, T_g1, V_g, phi_g_out);
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
        const float phi_g_in = grid.phi_g_imh(0);
        if (release) {
            mixture_pack_ghost(grid, grid.mix_inner_boundary0, kInnerRhoRes0, 0.0,
                               kInnerTRef, phi_g_in);
            mixture_pack_ghost(grid, grid.mix_inner_boundary1, kInnerRhoRes1, 0.0,
                               kInnerTRef, phi_g_in);
            grid.broadcast();
            return;
        }
        Vec& ob0 = grid.inner_boundary0_i;
        ob0(cons::RHO_I) = kInnerRhoIGh;
        ob0(cons::RHO_N) = kInnerRhoNGh;
        ob0(cons::MOM_I) = 0.0f;
        ob0(cons::MOM_N) = 0.0f;
        ob0(cons::E_I)   = grid.inv_gm1() * kInnerPiGh + kInnerRhoIGh * phi_g_in;
        ob0(cons::E_N)   = grid.inv_gm1() * kInnerPnGh + kInnerRhoNGh * phi_g_in;
        ob0(cons::E_E)   = 0.5f * grid.inv_gm1() * kInnerPiGh;   // p_e = ½ p_i ⇒ T_e = T_i

        Vec& ob1 = grid.inner_boundary1_i;
        ob1(cons::RHO_I) = kInnerRhoIGh2;
        ob1(cons::RHO_N) = kInnerRhoNGh2;
        ob1(cons::MOM_I) = 0.0f;
        ob1(cons::MOM_N) = 0.0f;
        ob1(cons::E_I)   = grid.inv_gm1() * kInnerPiGh2 + kInnerRhoIGh2 * phi_g_in;
        ob1(cons::E_N)   = grid.inv_gm1() * kInnerPnGh2 + kInnerRhoNGh2 * phi_g_in;
        ob1(cons::E_E)   = 0.5f * grid.inv_gm1() * kInnerPiGh2;
    }

    grid.broadcast();
}

} // namespace chromosphere
