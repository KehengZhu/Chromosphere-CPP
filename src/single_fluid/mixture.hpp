/*!
 * single_fluid/mixture.hpp — the RELEASE solver: a single-fluid field-aligned
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
 * The release timestep is exactly
 *
 *     U^n -> MUSCL/Roe hydro -> U* -> implicit physical conduction -> U^{n+1}
 *
 * and nothing else. TRAC, nonthermal beam heating, volumetric coronal heating,
 * radiative cooling and artificial/numerical conduction are NOT release physics;
 * they exist only in the historical two-fluid research solver.
 *
 * The historical two-fluid solver lives in `two_fluid/two_fluid.hpp` and owns
 * its own, separate seven-row state. Nothing here includes it, calls it, or
 * depends on it; the only shared code is the Grid in `chromosphere.hpp` and the
 * equilibrium EOS in `eos.hpp`.
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

/// CFL-limited timestep from the acoustic signal speed |u| + c_s. That is the
/// only constraint: the release has no volumetric source term, and the implicit
/// conduction stage is unconditionally stable.
Vec mixture_timestep(const Grid& grid, const Vec& state,
                     const MixtureField& decoded);

// ============================================================================
// Implicit physical conduction
// ============================================================================

/// Nonlinear backward-Euler PHYSICAL conduction on the mixture energy row:
/// Spitzer electron conduction plus neutral-hydrogen conduction, with no TRAC
/// broadening and no artificial/numerical diffusivity anywhere in the operator.
/// Mass and momentum are untouched; the accepted energy is the flux-updated
/// internal energy plus the unchanged kinetic and gravitational parts, so the
/// stage is conservative by construction. The outer thermal datum, when a
/// scenario imposes one, sits at the PHYSICAL boundary face, half a top cell
/// from the top cell centre — not at a ghost centre.
Vec mixture_apply_conduction(const Grid& grid, const Vec& state, double dt);

/// Maximum cellwise normalized residual of the backward-Euler conduction
/// equation, for diagnostics and regression tests.
double mixture_conduction_residual_max(const Grid& grid, const Vec& before,
                                       const Vec& after, double dt);

// ============================================================================
// Time integration
// ============================================================================

/// The complete release timestep:
///     U^n -> MUSCL/Roe hydro -> U* -> implicit physical conduction -> U^{n+1}
/// These two stages are the whole update; there is no third stage and no
/// equilibrium projection. `decoded` must describe `state`.
Vec mixture_advance(Grid& grid, const Vec& state, const Vec& dt_i,
                    const MixtureField& decoded);

} // namespace chromosphere
