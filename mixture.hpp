/*!
 * mixture.hpp — the RELEASE solver: a single-fluid field-aligned
 * equilibrium-mixture model.
 *
 * ## Physical model
 *
 * One fluid, one velocity, one temperature, on a prescribed flux tube of
 * cross-section A(s) with A(s) B(s) = const:
 *
 *   d_t(rho)     + (1/A) d_s(A rho u)            = 0
 *   d_t(rho u)   + (1/A) d_s(A (rho u^2 + p))    = p d_s(ln A) - rho d_s(phi_g)
 *   d_t(E)       + (1/A) d_s(A (E + p) u)        = -(1/A) d_s(A q)
 *
 * with the total energy, the equilibrium closure, and the conductive flux
 *
 *   E     = e_int(rho,T) + 1/2 rho u^2 + rho phi_g
 *   e_int = 3/2 p + x n_H chi_H,   p = (1 + x) n_H k_B T,   n_H = rho / m_H
 *   x     = x_Saha(n_H, T)   (pure-hydrogen Saha equilibrium)
 *   q     = -kappa(n_e, n_HI, T) d_s T
 *   c_s   = sqrt(Gamma1(T, n_H) p / rho)   (CRASH Gamma1 table)
 *
 * The gravitational potential is carried inside E, so `phi_g` appears in the
 * decode of every cell, ghost and reconstructed face.
 *
 * ## Software contract
 *
 * The conserved state is EXACTLY the three rows of `mix::` — there is no
 * carrier-row storage, no equilibrium projection stage, and no electron
 * compatibility row. The ionization fraction x, the carrier densities
 * rho_i = x rho and rho_n = (1-x) rho, n_e and n_HI, and the electron partial
 * pressure p_e are DERIVED diagnostics of (rho, e_int) produced by the EOS.
 *
 * The historical two-fluid solver lives in `chromosphere.hpp` and owns its own,
 * separate seven-row state; nothing here depends on it.
 */
#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

// ============================================================================
// State decode
// ============================================================================

/// Decode one immutable release state plus its two existing ghost layers.
/// An optional previous field supplies same-cell temperature guesses; failed
/// guesses always fall back to the full bracketed inversion, so the result is
/// independent of them.
MixtureField mixture_decode(const Grid& grid, const Vec& state,
                            std::uint64_t state_generation = 0,
                            const MixtureField* previous = nullptr);
void mixture_decode_into(const Grid& grid, const Vec& state,
                         MixtureField& output,
                         std::uint64_t state_generation = 0,
                         const MixtureField* previous = nullptr);

/// Pack one release cell/ghost from primitive data. E is built through the
/// equilibrium caloric closure, so the packed state is on the Saha manifold by
/// construction.
void mixture_pack_cell(const Grid& grid, Vec& state, arma::uword cell,
                       double rho, double velocity, double temperature,
                       double phi);
/// Same, into a length-num_of_mixture_eq ghost buffer.
void mixture_pack_ghost(const Grid& grid, Vec& ghost, double rho,
                        double velocity, double temperature, double phi);

/// Cell-centre gravitational potential, the gauge every decode uses.
inline double mixture_cell_phi(const Grid& grid, arma::uword i) {
    return 0.5 * (static_cast<double>(grid.phi_g_imh(i))
                + static_cast<double>(grid.phi_g_iph(i)));
}

// ============================================================================
// Hydrodynamics
// ============================================================================

/// Explicit field-aligned RHS: MUSCL reconstruction of (log rho, u, log p) with
/// the MC3/Koren limiter, a predictor--corrector half step, the mixture Roe
/// characteristic flux (face-local Rusanov fallback), the geometric
/// pressure/gravity source, and the frozen equilibrium-reference correction.
/// `dt_predictor` is the timestep the internal half step uses; it never leaves
/// this call.
Vec mixture_rhs_explicit(const Grid& grid, const Vec& state,
                         const MixtureField& decoded, double dt_predictor);

/// CFL-limited timestep from the acoustic signal speed |u| + c_s, plus the
/// optional heating-timescale caps.
Vec mixture_timestep(const Grid& grid, const Vec& state,
                     const MixtureField& decoded);

// ============================================================================
// Implicit physical conduction and the optional energy source stages
// ============================================================================

/// Nonlinear backward-Euler physical conduction on the mixture energy row.
/// Mass and momentum are untouched; the accepted energy is the flux-updated
/// internal energy plus the unchanged kinetic and gravitational parts, so the
/// stage is conservative by construction.
Vec mixture_apply_conduction(const Grid& grid, const Vec& state, double dt);

/// Maximum cellwise normalized residual of the backward-Euler conduction
/// equation, for diagnostics and regression tests.
double mixture_conduction_residual_max(const Grid& grid, const Vec& before,
                                       const Vec& after, double dt);

/// TRAC adaptive cutoff temperature of a release state (Johnston et al. 2020).
float mixture_trac_cutoff_T(const Grid& grid, const Vec& state);

/// Replace the internal energy of every cell, keeping rho and rho u. Used by the
/// optional volumetric heating and radiative-cooling stages. Throws if a target
/// leaves the EOS caloric domain.
Vec mixture_set_internal_energy(const Grid& grid, const Vec& state,
                                const std::vector<double>& target_internal_energy);

// ============================================================================
// Time integration
// ============================================================================

/// One release timestep:
///     U^n -> MUSCL/Roe hydro -> U* -> implicit physical conduction -> U^{n+1}
/// with the optional beam / coronal-heating / radiative-cooling energy stages
/// applied after conduction when enabled. `decoded` must describe `state`.
Vec mixture_advance(Grid& grid, const Vec& state, const Vec& dt_i,
                    const MixtureField& decoded);

} // namespace chromosphere
