#include "model_column.hpp"
#include "model_c7.hpp"   // c7_full_profile + c7_route_b_photoionization (C7 library)
#include "mesh.hpp"       // shared static local-refinement mesh builder

#include <cmath>
#include <cstdlib>
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
// Inner (photospheric) discrete-HSE reservoir, captured at the IC base. Both
// ghosts continue the hydrostatic reservoir and carry ρ ∝ p at the base T₀, so
// the inner-face reconstruction is consistent and V = 0 is an exact fixed point
// (docs/gentle_evaporation_downflow.md; the well-balanced inner BC — always on).
float kInnerRhoI   = 0.0f;      // base ion density (RHO slot ghost 0)
float kInnerRhoN   = 0.0f;
float kInnerPiGh   = 0.0f;      // ion ghost pressure   p_i,0 + ρ_i,0 g Δs
float kInnerPnGh   = 0.0f;      // neutral ghost pressure p_n,0 + ρ_n,0 g Δs
float kInnerPiGh2  = 0.0f;      // second inner ghost (2 Δs below): p_0 + 2ρgΔs
float kInnerPnGh2  = 0.0f;
float kInnerRhoIGh  = 0.0f;     // isothermal-hydrostatic ghost densities (ρ ∝ p at T₀):
float kInnerRhoNGh  = 0.0f;     //   ρ_{i,n},G0 = ρ_{i,n},0 · p_{i,n},G0 / p_{i,n},0
float kInnerRhoIGh2 = 0.0f;     //   ρ_{i,n},G1 = ρ_{i,n},0 · p_{i,n},G1 / p_{i,n},0
float kInnerRhoNGh2 = 0.0f;
float kTtopRef     = 0.0f;      // FIXED reference top temperature [K] for the jump
bool  kTtopFixedAbs = false;    // ISO_T_TOP>0: kTtopRef is an ABSOLUTE imposed top T
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

// Pack a ghost cell from explicit ion/neutral densities, temperature, velocity and
// the gravitational potential (p_i = 2 n_i k T incl. electrons, p_n = n_n k T,
// ε_e = 3/2 p_e ⇒ T_e = T_i). Taking (rho_i, rho_n) explicitly lets each ghost
// carry the LOCAL ionization fraction (ionized corona ↔ neutral chromosphere).
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
} // namespace

Vec model_column_ic(Grid& grid) {
    auto env_f = [](const char* key, float fallback) -> float {
        if (const char* e = std::getenv(key)) {
            try { return std::stof(std::string(e)); } catch (...) {}
        }
        return fallback;
    };

    // --- IC mode ----------------------------------------------------------
    // The IC is ALWAYS the real Model C7 profile with the density re-integrated
    // hydrostatically (clean V≈0 start; the analytic-isentrope toy IC is retired).
    // ISO_CORONA extends the top up into a resolved ~1 MK corona (to
    // ISO_CORONA_TOP_KM) and turns on TRAC so the ionized corona conducts via
    // Spitzer κ_e and can drive A&S78 conduction-driven evaporation.
    kCorona = (env_f("ISO_CORONA", 0.0f) != 0.0f);

    // --- parameters -------------------------------------------------------
    const float h_base     = env_f("ISO_H_BASE", 0.0f);                     // km, gauge (photosphere)
    const float corona_top = env_f("ISO_CORONA_TOP_KM", 10000.0f);          // km, coronal top
    float       DH_km      = env_f("ISO_DH", kCorona ? (corona_top - h_base) : 1303.0f); // km
    // Base (T, n_H) from the C7 profile at h_base.
    float T_base, ni_b, nn_b;
    c7_cell(h_base, T_base, ni_b, nn_b);
    const bool  heat_flux_on = (env_f("ISO_HEAT_FLUX", 0.0f) != 0.0f);
    // Stage 2b radiative sink. Gentle evaporation is the competition between the
    // downward conductive flux and radiative cooling (Antiochos & Sturrock 1978).
    const bool  cooling_on = heat_flux_on && (env_f("ISO_COOLING", 0.0f) != 0.0f);
    // Full-physics opt-ins (default off ⇒ clean single-fluid, frozen-ionization):
    const bool  two_fluid_on  = (env_f("ISO_TWO_FLUID", 0.0f) != 0.0f);
    const bool  ionization_on = (env_f("ISO_IONIZATION", 0.0f) != 0.0f);
    kHeatFluxOn        = heat_flux_on;
    kAjump             = env_f("ISO_TJUMP_A", 1.0f);
    kBjump             = env_f("ISO_TJUMP_B", 1.0f);
    kVcapMach          = env_f("ISO_VCAP", 0.1f);
    kInnerTNeumann     = (env_f("ISO_INNER_T_NEUMANN", 0.0f) != 0.0f);

    // Adiabatic index. Default 5/3. ISO_GAMMA near 1 (e.g. 1.05) makes the gas
    // nearly isothermal — a polytropic stand-in for the radiative thermostat the
    // conduction-only experiment leaves out. Set BEFORE the IC/BC build.
    grid.gamma_mono = env_f("ISO_GAMMA", grid.gamma_mono);

    // --- geometry: straight field line, gravity on -----------------------
    // Static local refinement (scenarios/mesh.hpp). peek_ns() sized the grid to the
    // coarse-equivalent count (ISO_NS); ISO_REFINE_* add lower-domain cells and grade
    // back to the outer coarse spacing. ds_m stays the coarse-equivalent width so the
    // numerical diffusivity keeps its magnitude. Uniform mesh reproduced byte-for-byte.
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

    // --- profile (T, n_i, n_n) sampled at cell centers --------------------
    // Real C7 T(h) and ionization FRACTION with the total density re-integrated
    // hydrostatically from the base (p(h)=p_base·exp(−∫ g/(R_s T) dh), R_s local),
    // so V≈0 is a clean discrete fixed point and any flow is heat-flux driven.
    Vec  T_c(grid.ns), ni_c(grid.ns), nn_c(grid.ns);
    float ln_p = 0.0f, prev_h = h_base, prev_invH = 0.0f;
    {
        const float xb = ni_b / (ni_b + nn_b);
        ln_p      = std::log((2.0f * ni_b + nn_b) * grid.k_b * T_base);   // base total pressure
        prev_invH = grid.g * grid.m_i / ((1.0f + xb) * grid.k_b * T_base);
    }
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float h_c = 0.5f * (h_F(i) + h_F(i + 1));        // km
        float T, n_i, n_n;
        c7_cell(h_c, T, n_i, n_n);
        const float x    = n_i / (n_i + n_n);                  // keep C7's ionization fraction
        const float invH = grid.g * grid.m_i / ((1.0f + x) * grid.k_b * T);
        const float dh_m = (h_c - prev_h) * 1000.0f;
        ln_p     -= 0.5f * (prev_invH + invH) * dh_m;
        const float n_tot = std::exp(ln_p) / ((1.0f + x) * grid.k_b * T);
        n_i = x * n_tot;
        n_n = (1.0f - x) * n_tot;
        prev_h = h_c; prev_invH = invH;
        T_c(i)  = T;
        ni_c(i) = n_i;
        nn_c(i) = n_n;
    }

    Vec xn = arma::zeros<Vec>(grid.n_state);
    const auto sz = arma::size(grid.ns, num_of_eq);
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float T     = T_c(i);
        const float n_i   = ni_c(i);
        const float n_n   = nn_c(i);
        const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
        const float rho_i = n_i * grid.m_i;
        const float rho_n = n_n * grid.m_n;
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
            const float rho_i = xn(arma::sub2ind(sz, i, cons::RHO_I));
            const float rho_n = xn(arma::sub2ind(sz, i, cons::RHO_N));
            const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
            const float n_i = rho_i / grid.m_i;
            const float n_n = rho_n / grid.m_n;
            const float T   = T_c(i);
            const float w   = 0.5f * (1.0f + std::tanh((T - Tc) / Tw));   // ~0 chromo → ~1 corona
            const float Tb  = (1.0f + (T_boost - 1.0f) * w) * T;
            const float p_i_b = 2.0f * n_i * grid.k_b * Tb;               // protons + electrons
            const float p_n_b =        n_n * grid.k_b * Tb;
            xn(arma::sub2ind(sz, i, cons::E_I)) = grid.inv_gm1() * p_i_b + rho_i * phi_g;
            xn(arma::sub2ind(sz, i, cons::E_N)) = grid.inv_gm1() * p_n_b + rho_n * phi_g;
            xn(arma::sub2ind(sz, i, cons::E_E)) = grid.inv_gm1() * n_i * grid.k_b * Tb;   // T_e = Tb
        }
    }

    // --- inner (photospheric) discrete-HSE V=0 reservoir ------------------
    // The ghost sits one Δs below cell 0, so a hydrostatic reservoir carries the
    // extra weight ρg·Δs ⇒ p_ghost = p_0 + ρ_0 g Δs. BOTH ghosts continue the
    // reservoir (p_0 + 2ρgΔs for ghost 2) and carry ρ ∝ p at the base T₀, so the
    // inner-face MUSCL reconstruction is consistent and V = U = 0 is an exact
    // discrete fixed point for momentum AND mass (always-on well-balanced inner BC).
    {
        const float T0    = T_c(0);
        const float n_i0  = ni_c(0);
        const float n_n0  = nn_c(0);
        kInnerRhoI = n_i0 * grid.m_i;
        kInnerRhoN = n_n0 * grid.m_n;
        const float p_i0 = 2.0f * n_i0 * grid.k_b * T0;
        const float p_n0 =        n_n0 * grid.k_b * T0;
        // Center-to-ghost distance = the BASE cell's own width (mirrored ghost). On a
        // refined mesh cell 0 is a fine cell, so this is the fine Δs, NOT the coarse ds_m.
        const float ds_base = grid.ds_i(0);
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

    // --- runtime physics toggles ------------------------------------------
    grid.single_fluid             = !two_fluid_on;       // ISO_TWO_FLUID: separate ion/neutral fluids
    grid.enable_ionization        = ionization_on;       // ISO_IONIZATION: Stage-E network
    grid.enable_radiative_cooling = cooling_on;
    grid.enable_beam_heating      = false;
    grid.enable_coronal_heating   = false;
    // TRAC broadens the unresolved TR. ISO_TRAC force-enables it independent of
    // cooling/corona (auto = cooling || corona). ISO_TRAC_TCHROM sets the region base
    // T_b; ISO_TRAC_TCMAXFRAC lifts Johnston's 0.2·T_peak cap (needed on a capped
    // low-T domain where 0.2·T_peak < T_b would pin T_c to the floor).
    const float iso_trac = env_f("ISO_TRAC", -1.0f);   // -1 = auto (cooling||corona)
    grid.enable_trac              = (iso_trac >= 0.0f) ? (iso_trac != 0.0f)
                                                       : (cooling_on || kCorona);
    grid.trac_T_chrom             = env_f("ISO_TRAC_TCHROM", 2.0e4f);
    grid.trac_Tc_max_frac         = env_f("ISO_TRAC_TCMAXFRAC", 0.2f);
    grid.trac_cutoff_T            = grid.trac_T_chrom;
    grid.enable_vacuum_floor      = cooling_on || kCorona || ionization_on;
    grid.impose_outer_heat_flux   = false;   // default: heat enters via the ghost-T jump
    grid.enable_Te                = false;    // two-fluid (ion/neutral), NOT three-temperature

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
    kHeatFluxOn     = relaxing0 ? false : heat_flux_on;

    // Stage 1: conduction OFF (adiabatic Euler steady state). Stage 2: ON, driven
    // by the top T jump (or the imposed q(T) flux / coronal heating).
    grid.enable_conduction        = relaxing0 ? false : heat_flux_on;
    // Numerical diffusivity is folded into Stage D, so it acts only when conduction
    // is on. A grid-scaled value damps TR-gradient ringing; corona needs more (the
    // resolved-corona drainage front is near-transonic on a coarse grid).
    const float diff_mult = kCorona ? 4.0f : 1.0f;
    grid.numerical_diffusivity    = heat_flux_on ? (diff_mult * 2.0e3f * ds_m) : 0.0f;

    model_column_update_bc(grid, xn);
    grid.broadcast();

    // Equilibrium-reference well-balancing (always on): freeze this (V=0, HSE) IC as
    // the reference equilibrium. The integrator caches the explicit-RHS residual here
    // on the first step and subtracts it every step, so this hydrostatic state is an
    // exact discrete fixed point (V=0 held to round-off for any stratification).
    grid.eq_wb = true;
    grid.eq_state = xn;
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
    auto cell_state = [&](arma::uword i, float& rho_tot, float& p_tot,
                          float& T, float& V, float& x) {
        const float rho_i = xn(arma::sub2ind(sz, i, cons::RHO_I));
        const float rho_n = xn(arma::sub2ind(sz, i, cons::RHO_N));
        const float momI  = xn(arma::sub2ind(sz, i, cons::MOM_I));
        const float momN  = xn(arma::sub2ind(sz, i, cons::MOM_N));
        const float E_i   = xn(arma::sub2ind(sz, i, cons::E_I));
        const float E_n   = xn(arma::sub2ind(sz, i, cons::E_N));
        const float Vi    = momI / rho_i;
        const float Un    = momN / rho_n;
        const float phi_g = 0.5f * (grid.phi_g_imh(i) + grid.phi_g_iph(i));
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
        kStage2Started = true;
    }

    // ====================================================================
    // Outer face — bottom of the transition region (new upper BC):
    //   1. pressure    — hydrostatic (centered HSE ⇒ V=0 fixed point).
    //   2. temperature — imposed jump T_ghost1 = a·T_ref, T_ghost2 = b·T_ghost1.
    //   3. density     — from the EOS at (p_ghost, T_ghost).
    // ====================================================================
    {
        const arma::uword nl = grid.ns - 1;
        const float ds       = grid.ds_i(nl);
        const float phi_g_out = grid.phi_g_iph(nl);   // ghost decode uses this gauge

        float rho_top, p_top, T_top, V_top, x_top;
        cell_state(nl, rho_top, p_top, T_top, V_top, x_top);
        float rho_2, p_2, T_2, V_2, x_2;
        cell_state(grid.ns >= 2 ? nl - 1 : nl, rho_2, p_2, T_2, V_2, x_2);

        // (1) hydrostatic ghost pressures (centered HSE at the boundary cell).
        float p_g0 = p_2 - 2.0f * ds * rho_top * g;
        // (2) imposed temperature jump, FIXED relative to the reference top T. With
        //     ISO_CORONA the corona is resolved and a=b=1 (the ghost just extends it).
        const float T_g0 = a_live * kTtopRef;
        const float T_g1 = b_live * T_g0;
        // (3) EOS ghost density at the imposed (p, T), at the top-cell ionization
        //     fraction x_top (R_s = (1+x) k/m). Split into ion/neutral by x_top.
        if (p_g0 < 1.0e-12f) p_g0 = 1.0e-12f;
        const float rho_g0 = p_g0 * m_i / ((1.0f + x_top) * k_b * T_g0);
        float p_g1 = p_top - 2.0f * ds * rho_g0 * g;
        if (p_g1 < 1.0e-12f) p_g1 = 1.0e-12f;
        const float rho_g1 = p_g1 * m_i / ((1.0f + x_top) * k_b * T_g1);

        // Velocity. Stage 1 (relaxation): static reservoir V=0. Stage 2: Mach-capped
        // OUTFLOW so the evaporation upflow can leave, while |V| ≤ kVcapMach·c_s
        // prevents the ill-posed-inflow runaway.
        float V_g = 0.0f;
        if (kHeatFluxOn) {
            const float c_s = std::sqrt(2.0f * grid.gamma_mono * k_b * T_top / m_i);
            const float vcap = kVcapMach * c_s;
            V_g = V_top;
            if (V_g >  vcap) V_g =  vcap;
            if (V_g < -vcap) V_g = -vcap;
        }
        pack_ghost(grid.outer_boundary0_i, grid, x_top * rho_g0, (1.0f - x_top) * rho_g0,
                   T_g0, V_g, phi_g_out);
        pack_ghost(grid.outer_boundary1_i, grid, x_top * rho_g1, (1.0f - x_top) * rho_g1,
                   T_g1, V_g, phi_g_out);
    }

    // ====================================================================
    // Inner face — photosphere, discrete-HSE V=0 reservoir (Pandey 2024),
    // well-balanced through BOTH ghosts (pressure p_0+ρgΔs / p_0+2ρgΔs and
    // density ρ ∝ p at T₀). The lower-BC TEMPERATURE choice (Dirichlet vs
    // Neumann/zero-flux) is a CONDUCTION boundary condition imposed inside the
    // Stage-D solver (grid.inner_conduction_neumann); the hydro wall stays a
    // stable Dirichlet V=0 anchor either way.
    // ====================================================================
    {
        const float phi_g_in = grid.phi_g_imh(0);
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
