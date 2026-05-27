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

namespace chromosphere {

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

/// Effective hydrogen photoionization rate P_phot [s^-1] (per neutral atom).
/// Single-rate closure for the Lyα + Balmer-continuum two-step channel
/// (Carlsson & Stein 2002 ApJ 572, 626). For the 2-state reduction adopted
/// here this is the radiative complement to the collisional Voronov rate:
///   dn_i/dt |_phot = +P_phot · n_n   (does not deplete electron thermal pool
///                                     — energy comes from the radiation field,
///                                     not from local kinetic temperature).
/// Returns `grid.photoionization_rate_i` if the scenario IC populated a
/// per-cell profile (e.g. model_c7_ic calibrates a height-dependent P_phot
/// to make C7's tabulated (n_e, n_n, T) a fixed point of the Stage E
/// quadratic); otherwise falls back to the uniform scalar
/// `grid.photoionization_rate`.
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

} // namespace chromosphere
