/*!
 * Inline physics helpers — collision frequencies, heat conductivities, and
 * hydrogen ionization / recombination rate coefficients.
 *
 * kappa_e / kappa_n are pure formulas in their inputs (no Grid dependency).
 * nu_in needs k_b and m_i so it takes the Grid; the ionization/recombination
 * rates take the Grid for χ_H and k_b.
 */

#pragma once

#include "chromosphere.hpp"

#include <cmath>

namespace chromosphere {

/// TRAC broadening factor ε(T) (Johnston et al. 2020 §2.2): ε = (T_c/T)^{5/2}
/// inside the TRAC region T_b ≤ T < T_c, and ε = 1 elsewhere (chromosphere
/// below T_b and corona at/above T_c are untouched). The conduction stage
/// multiplies κ_e by ε (κ' = κ(T_c), broadening the TR) and the cooling stage
/// divides the optically-thin loss by ε (Λ' = Λ(T/T_c)^{5/2}), so κ'Λ' = κΛ and
/// the TR-integrated radiation is preserved. Returns all-ones when TRAC is off.
inline Vec trac_broadening_factor(const Grid& grid, const Vec& T) {
    Vec eps(T.n_elem, arma::fill::ones);
    if (!grid.enable_trac) return eps;
    const float Tc = grid.trac_cutoff_T;
    const float Tb = grid.trac_T_chrom;
    if (!(Tc > Tb)) return eps;                 // degenerate cutoff → no broadening
    for (arma::uword i = 0; i < T.n_elem; ++i) {
        const float t = T(i);
        if (t >= Tb && t < Tc) eps(i) = std::pow(Tc / t, 2.5f);
    }
    return eps;
}

/// Electron heat conductivity in SI units (writeup eq 53).
inline Vec kappa_e(const Vec& n_e, const Vec& n_n, const Vec& T_e) {
    return (9.2048e-12f * n_e % arma::pow(T_e, 2.5)
            / (n_e + 3.5609e-12f * n_n % T_e % T_e));
}

/// Neutral heat conductivity in SI units (writeup eq 59).
inline Vec kappa_n(const Vec& n_i, const Vec& n_n, const Vec& T_i, const Vec& T_n) {
    return (0.0342006f * n_n % T_n)
         / (1.20613f * n_i % arma::sqrt(T_n + T_i)
          + 1.70573f * n_n % arma::sqrt(T_n));
}

/// Ion-neutral collision frequency, target-density form (writeup eq 57):
///   ν_in = (2 a₀)² n_n √(8 π k_b (T_i + T_n) / m_i)
inline Vec nu_in(const Grid& grid, const Vec& n_n, const Vec& T_i, const Vec& T_n) {
    const float bohr_r = 53e-12f;
    return (2.0f * bohr_r) * (2.0f * bohr_r) * n_n
         % arma::sqrt(arma::abs(8.0f * static_cast<float>(arma::datum::pi)
                                * grid.k_b * (T_i + T_n) / grid.m_i));
}

/// Hydrogen electron-impact ionization rate coefficient S_i(T_e) [m^3/s].
/// Voronov (1997) fit in SI; cgs original is
///   S_i = A (1+P√U) U^K e^{-U} / (X+U),
/// with A = 0.291×10^{-7} cm^3/s, P = 0, X = 0.232, K = 0.39, and
/// U = χ_H / (k_b T_e). The cm^3 → m^3 conversion gives 2.91×10^{-14}.
inline Vec ionization_rate_S(const Grid& grid, const Vec& T_e) {
    const Vec U = grid.chi_H_J / (grid.k_b * arma::clamp(T_e, 1.0f, arma::datum::inf));
    return 2.91e-14f * arma::pow(U, 0.39f) % arma::exp(-U) / (0.232f + U);
}

/// Hydrogen case-B radiative recombination rate coefficient α_r(T_e) [m^3/s].
/// Low-T power-law fit to Hummer (1994); accurate to ~10% for 10^3–10^5 K.
/// Case B is the consistent choice with the optically-thin χ_H bookkeeping
/// used by Stage E (writeup §6.5).
inline Vec recombination_rate_alpha(const Grid& /*grid*/, const Vec& T_e) {
    return 2.7e-19f * arma::pow(arma::clamp(T_e, 1.0f, arma::datum::inf) / 1.0e4f, -0.75f);
}

/// Saha ionization function Φ(T_e) [m^-3] for hydrogen (writeup eq:Seff):
///   Φ = (2π m_e k_B T_e / h^2)^{3/2} e^{-χ_H / k_B T_e}.
/// At LTE the ground-state Saha relation is n_e n_i / n_n = Φ, so Φ ties the
/// multilevel collisional ionization S_CR to its three-body inverse α_c.
/// NOTE: the prefactor coefficient 2π m_e k_B / h² ≈ 1.79985×10^{14} m^-2 K^-1
/// is lumped into a single constant (computed in double precision) because the
/// intermediates m_e k_B (~8×10^-53) and h² (~4.4×10^-67) both underflow to
/// zero in float32 — forming them separately would divide by zero.
inline Vec saha_phi(const Grid& grid, const Vec& T_e) {
    constexpr float SAHA_PREF_C = 1.79985e14f;   // 2π m_e k_B / h² [m^-2 K^-1]
    const Vec Tc = arma::clamp(T_e, 1.0f, arma::datum::inf);
    const Vec pref = arma::pow(SAHA_PREF_C * Tc, 1.5f);
    return pref % arma::exp(-grid.chi_H_J / (grid.k_b * Tc));
}

/// Three-body (collisional–radiative) recombination, density coefficient
/// κ_c(T_e) [m^6/s] such that α_c = κ_c · n_e. Hinnov & Hirschberg (1962) /
/// Stevefelt (1975) low-T form; cgs original
///   α_c = 5.6×10^{-27} (k_B[eV] T_e)^{-4.5} n_e[cm^-3]  [cm^3/s].
/// The cm^3→m^3 (×1e-6) and cm^-3→m^-3 (×1e-6 on n_e) conversions give the
/// 5.6×10^{-39} SI prefactor below. Keeping κ_c (the n_e-independent part)
/// separate lets Stage E carry the three-body channel as the cubic f^3 term.
inline Vec recombination_rate_kappa_c(const Grid& /*grid*/, const Vec& T_e) {
    constexpr float k_b_eV = 8.617333262e-5f;   // Boltzmann constant [eV/K]
    const Vec kT_eV = k_b_eV * arma::clamp(T_e, 1.0f, arma::datum::inf);
    // κ_c = 5.6e-27 (k_B[eV] T)^{-4.5} cm^3/s · 1e-12 (cm^3·cm^-3 → m^6) [m^6/s].
    // Apply the 1e-12 conversion last so the 5.6e-39 product never forms as a
    // (low-precision) float32 denormal.
    return (5.6e-27f * arma::pow(kT_eV, -4.5f)) * 1.0e-12f;   // [m^6/s]
}

/// Three-body recombination coefficient α_c(T_e, n_e) = κ_c · n_e [m^3/s].
inline Vec recombination_rate_alpha_c(const Grid& grid, const Vec& T_e, const Vec& n_e) {
    return recombination_rate_kappa_c(grid, T_e) % n_e;
}

/// Multilevel (Rydberg-ladder) collisional ionization coefficient
/// S_CR(T_e) [m^3/s], the detailed-balance inverse of three-body
/// recombination: S_CR = α_c Φ / n_e = κ_c Φ. Because α_c ∝ n_e, S_CR is
/// independent of n_e (a function of T_e alone). This is the dominant
/// collisional ionization channel in the chromosphere — orders of magnitude
/// larger than the direct ground-state Voronov S_i, since it proceeds through
/// high-n levels (writeup §3.1 eq:Seff).
inline Vec ionization_rate_S_CR(const Grid& grid, const Vec& T_e) {
    return recombination_rate_kappa_c(grid, T_e) % saha_phi(grid, T_e);
}

/// Effective hydrogen photoionization rate P_phot [s^-1] (per neutral atom).
/// Single-rate closure for the Balmer-continuum (n=2) channel that sets
/// chromospheric H ionization (Carlsson & Stein 2002 ApJ 572, 626; the Lyman
/// continuum is negligible). For the 2-state reduction adopted here this is
/// the radiative complement to the collisional Voronov rate:
///   dn_i/dt |_phot = +P_phot · n_n   (does not deplete electron thermal pool
///                                     — energy comes from the radiation field,
///                                     not from local kinetic temperature).
/// Returns `grid.photoionization_rate_i` if the scenario IC populated a
/// per-cell profile. model_c7_ic populates it with the height-structured
/// frozen-radiation-field rate R_ik(h) from Chae (2021) Table 1 (FAL-C),
/// via photoionization_rate_chae() below; otherwise this falls back to the
/// uniform scalar `grid.photoionization_rate`.
inline Vec photoionization_rate_P(const Grid& grid) {
    if (grid.photoionization_rate_i.n_elem == grid.ns) {
        return grid.photoionization_rate_i;
    }
    return Vec(grid.ns, arma::fill::value(grid.photoionization_rate));
}

// ============================================================================
// Optically-thick chromospheric radiative cooling
//   Carlsson & Leenaarts 2012 A&A 539 A39, Eq. 1:
//     Q_X = -L_X(T) E_X(tau or m_c) (N_X*/N_X)(T) A_X (N_H/rho) n_e rho
//   for X in {H I, Ca II, Mg II}. Tables digitized from CL2012 Figs 3-11
//   (see writeup §6 Tables 1-3). Returns Q_rad in W/m^3, positive for
//   net cooling — caller subtracts to add to e_i.
// ============================================================================

namespace cl2012 {

// Table 1: log10 L_X(T) [erg s^-1 per e- per X-particle], T in kK.
// Below the first T entry, L_X is held at the first value (a large negative
// log -> effectively zero), then linear-interpolated. CL2012 Figs 3-5.
constexpr int N_L = 13;
static constexpr float L_T_kK[N_L]   = { 4.f,  5.f, 6.f,  7.f,  8.f, 10.f, 12.f, 15.f, 20.f, 25.f, 30.f, 40.f, 50.f};
static constexpr float L_logL_H[N_L] = {-30.f, -25.0f,-24.5f,-24.0f,-23.6f,-22.8f,-22.2f,-21.85f,-21.75f,-21.85f,-22.0f,-22.4f,-22.8f};
static constexpr float L_logL_Ca[N_L]= {-22.5f,-21.3f,-20.5f,-19.9f,-19.5f,-18.9f,-18.6f,-18.35f,-18.15f,-18.05f,-18.0f,-18.0f,-18.0f};
static constexpr float L_logL_Mg[N_L]= {-23.0f,-21.7f,-20.5f,-19.7f,-19.2f,-18.7f,-18.4f,-18.2f, -18.05f,-18.0f, -18.0f,-18.0f,-18.0f};

// Table 2a: E_H(log10 tau_Lya), Fig 6.
constexpr int N_EH = 11;
static constexpr float EH_logTau[N_EH] = {-2.f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 6.0f};
static constexpr float EH_E[N_EH]      = { 1.0f,1.00f,0.99f,0.95f,0.85f,0.70f,0.55f,0.40f,0.25f,0.15f,0.03f};

// Table 2b: E_Ca(log10 m_c [g/cm^2]), Fig 7.
constexpr int N_ECa = 11;
static constexpr float ECa_logMc[N_ECa] = {-8.f, -7.0f, -6.5f, -6.0f, -5.5f, -5.0f, -4.5f, -4.0f, -3.5f, -3.0f, -2.0f};
static constexpr float ECa_E[N_ECa]     = { 1.0f, 0.99f, 0.97f, 0.90f, 0.75f, 0.58f, 0.43f, 0.30f, 0.20f, 0.12f, 0.02f};

// Table 2c: E_Mg(log10 m_c [g/cm^2]), Fig 8.
constexpr int N_EMg = 10;
static constexpr float EMg_logMc[N_EMg] = {-8.f, -6.0f, -5.5f, -5.0f, -4.5f, -4.0f, -3.5f, -3.0f, -2.5f, -2.0f};
static constexpr float EMg_E[N_EMg]     = { 1.0f, 1.00f, 0.97f, 0.85f, 0.65f, 0.45f, 0.30f, 0.18f, 0.10f, 0.04f};

// Table 3: N(X*)/N(X) (T in kK), Figs 9-11 red curves.
constexpr int N_F = 11;
static constexpr float F_T_kK[N_F] = { 5.f, 8.f, 10.f, 12.f, 14.f, 15.f, 16.f, 18.f, 20.f, 25.f, 30.f};
static constexpr float F_HI[N_F]   = { 1.0f, 0.99f, 0.85f, 0.42f, 0.17f, 0.10f, 0.07f, 0.04f, 0.025f,0.010f,0.005f};
static constexpr float F_Ca[N_F]   = { 1.0f, 0.97f, 0.90f, 0.75f, 0.55f, 0.42f, 0.30f, 0.15f, 0.08f, 0.025f,0.010f};
static constexpr float F_Mg[N_F]   = { 0.95f,0.95f, 0.90f, 0.75f, 0.55f, 0.42f, 0.28f, 0.10f, 0.05f, 0.010f,0.005f};

// Asplund et al. 2009 abundances and N_H/rho [g^-1].
constexpr float A_H        = 1.0f;
constexpr float A_Ca       = 2.19e-6f;   // 10^(6.34-12)
constexpr float A_Mg       = 3.98e-5f;   // 10^(7.60-12)
constexpr float NH_over_rho_cgs = 4.407e23f; // per gram

/// 1-D linear interpolation: piecewise linear, constant extrapolation.
inline float lin_interp(float x, const float* xs, const float* ys, int n) {
    if (x <= xs[0])   return ys[0];
    if (x >= xs[n-1]) return ys[n-1];
    int k = 0;
    while (k < n-1 && xs[k+1] < x) ++k;
    const float t = (x - xs[k]) / (xs[k+1] - xs[k]);
    return ys[k] + t * (ys[k+1] - ys[k]);
}

/// Vectorized look-up wrapper. Returns one Vec same size as q.
inline Vec lookup_vec(const Vec& q, const float* xs, const float* ys, int n) {
    Vec out(q.n_elem);
    for (arma::uword i = 0; i < q.n_elem; ++i) out[i] = lin_interp(q[i], xs, ys, n);
    return out;
}

} // namespace cl2012

// ============================================================================
// Frozen-radiation-field hydrogen photoionization rate R_ik(z)
//   Chae (2021) J. Astron. Space Sci. 38, 83, Table 1 — the FAL-C
//   (Fontenla, Avrett & Loeser 1993) "frozen radiation field" photoionization
//   rate tabulated vs height z [km] on the tau_500=1 scale, as log10 R_ik
//   [s^-1]. Replaces the earlier C7-calibrated P_phot inversion (which the
//   inversion P = α_r n_e²/n_n − S_i n_e made a Stage-E fixed point of C7 but
//   produced an unphysical, ~exact-balance profile). R_ik dips to ~7e-8 s^-1
//   at the temperature minimum (weak radiation field) and climbs to ~7e-3
//   s^-1 at the chromosphere top. Carlsson & Stein 2002 show the photoionizing
//   (Balmer-continuum, n=2) field varies only ~1–2% even through shocks, so a
//   height-only fixed rate is valid (Sollum 1999; Leenaarts 2007, 2020).
// ============================================================================
namespace chae2021 {
constexpr int N_R = 43;
// z [km] on the tau_500=1 scale (h=0 = photospheric reference level).
static constexpr float Z_km[N_R] = {
    -100.f, -80.f, -60.f, -40.f, -20.f,    0.f,   50.f,  100.f,  150.f,  200.f,
     250.f, 300.f, 350.f, 400.f, 450.f,  490.f,  525.f,  560.f,  600.f,  650.f,
     705.f, 755.f, 805.f, 855.f, 905.f,  980.f, 1065.f, 1180.f, 1278.f, 1378.f,
    1475.f,1580.f,1670.f,1775.f,1860.f, 1915.f, 1980.f, 2017.f, 2043.f, 2062.f,
    2075.f,2087.f,2110.f};
// log10 R_ik [s^-1].
static constexpr float logR[N_R] = {
     1.25f, 0.78f, 0.15f,-0.58f,-1.40f,-2.12f,-3.50f,-4.38f,-5.07f,-5.56f,
    -5.91f,-6.22f,-6.51f,-6.76f,-7.01f,-7.13f,-7.15f,-7.07f,-6.79f,-6.33f,
    -5.74f,-5.26f,-4.88f,-4.55f,-4.38f,-4.17f,-3.97f,-3.73f,-3.56f,-3.40f,
    -3.26f,-3.11f,-2.99f,-2.87f,-2.76f,-2.67f,-2.54f,-2.44f,-2.37f,-2.30f,
    -2.26f,-2.22f,-2.15f};
} // namespace chae2021

/// Frozen-field hydrogen photoionization rate P_phot = R_ik(h) [s^-1] from the
/// Chae (2021) FAL-C Table 1 lookup. `h_km` is height on the tau_500=1 scale.
/// log10 R_ik is linearly interpolated in height (constant extrapolation past
/// the table ends — Chae's geometric-dilution correction Eq. 26 is negligible
/// over the chromospheric range used here), then base-10 exponentiated.
inline Vec photoionization_rate_chae(const Vec& h_km) {
    constexpr float LN10 = 2.302585092994046f;
    const Vec logr = cl2012::lookup_vec(h_km, chae2021::Z_km, chae2021::logR,
                                        chae2021::N_R);
    return arma::exp(LN10 * logr);
}

/// Optically-thick chromospheric radiative cooling Q_rad [W/m^3].
/// Positive = net cooling (caller subtracts from energy row).
/// CL2012 Eq. 1; tables in namespace cl2012 above.
/// `s_up_is_outer` matches the grid convention: index 0 = inner (photosphere),
/// index ns-1 = outer (TR). Column integrals are accumulated downward from
/// the TR (top of column = zero column mass).
inline Vec radiative_loss_thick(const Grid& grid,
                                const Vec& n_i, const Vec& n_n,
                                const Vec& T_e) {
    using namespace cl2012;
    const arma::uword ns = grid.ns;

    // ---- column integrals (CGS) ----
    // m_c [g/cm^2] = 0.1 * integral_{s}^{s_top} rho_SI ds_SI (SI -> cgs surface)
    // N_HI [/cm^2] = 1e-4 * integral n_n[/m^3] ds[m]
    const Vec rho = (n_i * grid.m_i + n_n * grid.m_n);   // SI kg/m^3
    Vec m_c(ns, arma::fill::zeros);                       // g/cm^2
    Vec N_HI(ns, arma::fill::zeros);                      // /cm^2
    // Top cell: half-cell column above the midpoint.
    m_c[ns-1]  = 0.5f * rho[ns-1]    * grid.ds_i[ns-1] * 0.1f;
    N_HI[ns-1] = 0.5f * n_n[ns-1]    * grid.ds_i[ns-1] * 1.0e-4f;
    for (int i = static_cast<int>(ns) - 2; i >= 0; --i) {
        // Trapezoid: mean density times cell-center separation.
        const float ds_face   = 0.5f * (grid.ds_i[i] + grid.ds_i[i+1]);
        const float rho_dl    = 0.5f * (rho[i] + rho[i+1]) * ds_face;
        const float nn_dl     = 0.5f * (n_n[i] + n_n[i+1]) * ds_face;
        m_c[i]  = m_c[i+1]  + rho_dl * 0.1f;
        N_HI[i] = N_HI[i+1] + nn_dl  * 1.0e-4f;
    }

    // ---- log coordinates for table lookup ----
    const Vec T_kK    = T_e / 1000.0f;
    const Vec log_mc  = arma::log10(arma::clamp(m_c,  1.0e-12f, 1.0e6f));
    const Vec tau_Lya = 4.0e-14f * N_HI;
    const Vec log_tau = arma::log10(arma::clamp(tau_Lya, 1.0e-6f, 1.0e10f));

    // ---- per-species factors ----
    constexpr float LN10 = 2.302585092994046f;
    const Vec L_H   = arma::exp(LN10 * lookup_vec(T_kK, L_T_kK, L_logL_H,  N_L));
    const Vec L_Ca  = arma::exp(LN10 * lookup_vec(T_kK, L_T_kK, L_logL_Ca, N_L));
    const Vec L_Mg  = arma::exp(LN10 * lookup_vec(T_kK, L_T_kK, L_logL_Mg, N_L));

    const Vec E_H_  = lookup_vec(log_tau, EH_logTau,  EH_E,  N_EH);
    const Vec E_Ca_ = lookup_vec(log_mc,  ECa_logMc, ECa_E, N_ECa);
    const Vec E_Mg_ = lookup_vec(log_mc,  EMg_logMc, EMg_E, N_EMg);

    const Vec f_HI  = lookup_vec(T_kK, F_T_kK, F_HI, N_F);
    const Vec f_Ca  = lookup_vec(T_kK, F_T_kK, F_Ca, N_F);
    const Vec f_Mg  = lookup_vec(T_kK, F_T_kK, F_Mg, N_F);

    // ---- assemble Q_X [erg/cm^3/s] ----
    // n_e_cgs = n_i [SI /m^3] * 1e-6 ; rho_cgs = rho [kg/m^3] * 1e-3
    const Vec n_e_cgs = 1.0e-6f * n_i;
    const Vec rho_cgs = 1.0e-3f * rho;
    const Vec NH_n_e_rho_cgs = NH_over_rho_cgs * n_e_cgs % rho_cgs;  // = N_H * n_e [/cm^6]

    const Vec Q_H  = L_H  % E_H_  % f_HI % NH_n_e_rho_cgs * A_H;
    const Vec Q_Ca = L_Ca % E_Ca_ % f_Ca % NH_n_e_rho_cgs * A_Ca;
    const Vec Q_Mg = L_Mg % E_Mg_ % f_Mg % NH_n_e_rho_cgs * A_Mg;

    // cgs (erg/cm^3/s) -> SI (W/m^3): factor 0.1
    return 0.1f * (Q_H + Q_Ca + Q_Mg);
}

// ============================================================================
// Optically-thin transition-region / coronal radiative losses
//   Q_thin = n_e n_H Λ(T)   [W/m^3],  positive = net cooling.
//   Λ(T) is the standard coronal-abundance optically-thin radiative loss
//   function (CHIANTI-class; peak ≈ 4×10^-35 W m^3 ≈ 10^-21.35 erg cm^3 s^-1
//   near log T ≈ 5.1). It covers the 10^4–10^7 K regime that the optically-THICK
//   chromospheric CL2012 recipe (H I/Ca II/Mg II, effective only for T ≲ 3×10^4 K)
//   does not. This is the radiative sink the lower TR needs to re-radiate the
//   imposed downward coronal conductive flux q(T) (writeup §Boundary conditions,
//   RTV 1978); without it the imposed flux had no TR-temperature sink and ran
//   the top up to coronal temperatures. The enthalpy/advective flux — the other
//   dominant TR energy term (Bradshaw & Cargill 2010; Klimchuk et al. 2008) —
//   is already carried by the hydrodynamic advection, so adding this radiative
//   term completes the leading-order TR energy balance.
//
//   A smooth low-T switch-on (~2×10^4 K) stitches the two recipes: CL2012
//   governs the chromosphere, the thin loss governs the TR/corona, with a tanh
//   blend across the overlap so neither double-counts the other. Λ(T) tabulated
//   from standard sources (Cook et al. 1989; Klimchuk et al. 2008; CHIANTI,
//   Dere et al. 1997+); refine against CHIANTI v11 (Dufresne 2024) if needed.
//   NOTE: equilibrium Λ(T) underestimates lower-TR losses by up to ~×3 under
//   non-equilibrium ionization (Landi et al. 2012) — a known conservative bias.
// ============================================================================
namespace optthin {
constexpr int N_L = 17;
// log10 T [K]
static constexpr float logT_K[N_L] = {
    4.0f, 4.1f, 4.2f, 4.3f, 4.5f, 4.7f, 4.9f, 5.1f, 5.3f,
    5.5f, 5.7f, 5.9f, 6.1f, 6.3f, 6.5f, 6.7f, 7.0f};
// log10 Λ [erg cm^3 s^-1]
static constexpr float logLambda[N_L] = {
    -25.0f, -23.8f, -22.9f, -22.35f, -21.85f, -21.55f, -21.4f, -21.35f, -21.45f,
    -21.6f, -21.75f, -21.9f, -22.05f, -22.2f, -22.35f, -22.45f, -22.5f};
} // namespace optthin

/// Optically-thin radiative loss Q_thin [W/m^3] (positive = cooling). `n_i` is
/// the electron/proton density (n_e = n_i), `n_n` the neutral density so that
/// n_H = n_i + n_n is total hydrogen; Q = n_e n_H Λ(T). The 1e-13 factor is the
/// cgs→SI conversion for Λ (erg cm^3 s^-1 → W m^3).
inline Vec radiative_loss_thin(const Grid& /*grid*/, const Vec& n_i,
                               const Vec& n_n, const Vec& T_e) {
    constexpr float LN10 = 2.302585092994046f;
    const Vec Tc    = arma::clamp(T_e, 1.0f, arma::datum::inf);
    const Vec logTc = arma::log10(Tc);
    const Vec logLam = cl2012::lookup_vec(logTc, optthin::logT_K,
                                          optthin::logLambda, optthin::N_L);
    const Vec Lambda_SI = arma::exp(LN10 * logLam) * 1.0e-13f;   // W m^3
    const Vec n_H = n_i + n_n;                                   // total hydrogen
    // Smooth switch-on across ~2e4 K (stitch to optically-thick CL2012 below).
    const Vec w = 0.5f * (1.0f + arma::tanh((Tc - 2.0e4f) / 5.0e3f));
    return w % (n_i % n_H % Lambda_SI);                         // W m^-3
}

// ============================================================================
// Flare nonthermal-electron-beam heating (chromospheric evaporation driver)
//   Q_beam(s,t) = beam_flux · g(t) · φ(s)   [W/m³],  positive = heating.
//   The beam is injected from the corona at the outer boundary and stops in the
//   upper chromosphere. Its true thick-target stopping column (~10²⁰ cm⁻² for
//   tens of keV) is unresolved on the coarse C7 grid — a single top cell already
//   holds ~10²⁶ m⁻² — so we deposit the flux over a finite layer of thickness
//   beam_deposition_height [km] below the top, weighted by the local total
//   density n_tot (the densest reachable layer absorbs the most, as a stopping
//   beam does). φ(s) is normalized so ∫φ ds = 1, hence the column-integrated
//   heating equals beam_flux [W/m²]. g(t) is a flat-top window over
//   [t_on, t_on+τ] with cosine ramps of half-width beam_ramp.
//
//   This is the standard Fisher/RADYN explosive-evaporation driver expressed as
//   an interior source term (see docs/explosive_evaporation_plan.md). Explosive
//   evaporation requires beam_flux above the Fisher (1985) threshold
//   F_crit ≈ 7×10⁶ W m⁻² together with an active radiative sink.
// ============================================================================
inline Vec beam_heating_rate(const Grid& grid, const Vec& n_i, const Vec& n_n) {
    Vec Q = arma::zeros<Vec>(grid.ns);
    if (!grid.enable_beam_heating || grid.beam_flux <= 0.0f) return Q;

    // Temporal window g(t): cosine ramp up over [t_on, t_on+ramp], flat top,
    // cosine ramp down over [t_off-ramp, t_off]. Zero outside.
    const float t    = grid.sim_time;
    const float t_on = grid.beam_t_on;
    const float toff = grid.beam_t_on + grid.beam_duration;
    const float r    = (grid.beam_ramp > 1.0e-6f) ? grid.beam_ramp : 1.0e-6f;
    float g;
    if (t <= t_on || t >= toff) {
        return Q;                                  // beam off
    } else if (t < t_on + r) {
        const float x = (t - t_on) / r;            // 0→1
        g = 0.5f * (1.0f - std::cos(static_cast<float>(arma::datum::pi) * x));
    } else if (t > toff - r) {
        const float x = (toff - t) / r;            // 1→0
        g = 0.5f * (1.0f - std::cos(static_cast<float>(arma::datum::pi) * x));
    } else {
        g = 1.0f;                                  // flat top
    }
    if (g <= 0.0f) return Q;

    // Deposition window: cells whose center height (absolute km, C7 base =
    // 1003 km) lies in [beam_h_lo_km, beam_h_hi_km] — the upper chromosphere
    // below the TR. Weight ∝ n_tot (thick-target: densest reachable layer
    // absorbs the most).
    const Vec  ds   = grid.ds_i;                    // [m]
    const Vec  ntot = n_i + n_n;                    // [m^-3]
    const float H_BASE_KM = 1.003e3f;               // C7 base offset (chromo_main)
    Vec  w = arma::zeros<Vec>(grid.ns);
    float h_lo = 0.0f;
    for (arma::uword i = 0; i < grid.ns; ++i) {
        const float h_ctr_km = H_BASE_KM + (h_lo + 0.5f * ds(i)) * 1.0e-3f;
        if (h_ctr_km >= grid.beam_h_lo_km && h_ctr_km <= grid.beam_h_hi_km)
            w(i) = ntot(i);
        h_lo += ds(i);
    }
    // Normalize so ∫ φ ds = 1  ⇒  Σ w_i ds_i = 1 (φ = w / Σ w ds).
    const float norm = arma::dot(w, ds);
    if (norm <= 0.0f) return Q;                     // layer empty (shouldn't happen)

    Q = (grid.beam_flux * g / norm) * w;            // W m^-3
    return Q;
}

} // namespace chromosphere
