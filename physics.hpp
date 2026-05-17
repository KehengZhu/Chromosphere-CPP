/*!
 * Inline physics helpers — collision frequencies and heat conductivities.
 *
 * kappa_e / kappa_n are pure formulas in their inputs (no Grid dependency).
 * nu_in needs k_b and m_i so it takes the Grid.
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

} // namespace chromosphere
