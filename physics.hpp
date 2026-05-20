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

} // namespace chromosphere
