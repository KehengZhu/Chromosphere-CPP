/*!
 * model_column — the unified field-aligned chromosphere→corona column scenario
 * (formerly model_isentropic; absorbs model_gentle). This is the DEFAULT scenario.
 *
 * A single straight field line (B = 1, gravity on) initialized from the real
 * Model C7 atmosphere (Avrett & Loeser 2008, Table 26; via c7_full_profile in
 * model_c7.cpp) with a realistic FROZEN ionization profile — neutral chromosphere,
 * ionized corona — and driven to conduction-driven ("gentle", Antiochos & Sturrock
 * 1978) or, with a beam elsewhere, explosive evaporation.
 *
 * The numerical scheme is the well-balanced, high-order configuration validated in
 * docs/gentle_evaporation_downflow.md (the "bestwb" set) — ALWAYS ON here, no toggles:
 *   - well-balanced explicit reconstruction (φ_g-consistent shifted pressures),
 *   - log-space MUSCL reconstruction of the exponentially-stratified slots,
 *   - the asymmetric MC3/Koren limiter (β = 2),
 *   - equilibrium-reference ("δ-form") well-balancing (V = 0 exact for any strat.),
 *   - inner (photospheric) discrete-HSE reservoir well-balanced through BOTH ghosts,
 *     density included (ρ ∝ p at T₀) so V = 0 is an exact fixed point for mass too.
 * Together these hold a hydrostatic column at V ≈ 0 to round-off, so any residual
 * flow is physical (conduction-driven evaporation), not truncation drainage.
 *
 * IC: the C7 temperature and ionization FRACTION with the total density re-integrated
 * hydrostatically from the base (clean V ≈ 0 start). ISO_CORONA extends the domain up
 * into a resolved ~1 MK corona (needed for A&S78 evaporation: only an ionized corona
 * conducts via Spitzer κ_e, and only a resolved coronal VOLUME gives the evaporated
 * mass somewhere to flow into) with TRAC broadening the under-resolved TR.
 *
 * Conduction drivers (physical options):
 *   - top ghost-T jump (default): T_ghost = a·T_ref, b·a·T_ref sets the downward
 *     conductive flux through the Stage-D Dirichlet ghost-T (ISO_TJUMP_A/B, ISO_T_TOP).
 *   - imposed Neumann coronal flux q(T) with a cosine ramp (ISO_QFLUX*), the
 *     model_c7/RTV closure — an alternative to the ghost-T jump.
 *   - ambient volumetric coronal heating H(s) = E₀·exp(−d/s_H) (ISO_CHEAT*), the
 *     Aschwanden & Schrijver (2002) footpoint-anchored RTV static-balance term.
 *   - one-time IC coronal superheat above the TR (ISO_TBOOST*), a hot-corona IC that
 *     conducts down from t = 0.
 *
 * Physical env knobs (numerics are fixed; only physics is tunable):
 *   ISO_NS            grid cells (default 600 chromosphere / 1000 resolved corona)
 *   ISO_CORONA        0 (default) / 1 = extend into a resolved corona
 *   ISO_CORONA_TOP_KM resolved-corona domain top [km] (default 10000)
 *   ISO_H_BASE        domain base height [km] (default 0, photosphere)
 *   ISO_DH            domain height span [km] (default 1303 chromosphere→TR base)
 *   ISO_GAMMA         adiabatic index (default 5/3; ~1.05 ⇒ near-isothermal thermostat)
 *   ISO_HEAT_FLUX     0 = Stage-1 relaxation (conduction off) / 1 = conduction on
 *   ISO_RELAX_TIME    Stage-1 adiabatic relaxation phase length [s] (default 0)
 *   ISO_TJUMP_A/_B    top ghost-T jump factors a, b (default 1)
 *   ISO_T_TOP         absolute imposed top-boundary temperature [K] (default 0 = C7 top)
 *   ISO_VCAP          outer outflow Mach cap (default 0.1)
 *   ISO_INNER_T_NEUMANN  conduction lower BC: 0 = Dirichlet (default) / 1 = Neumann
 *   ISO_COOLING       optically-thick + thin radiative sink (Stage R)
 *   ISO_TWO_FLUID     evolve ion & neutral as separate fluids (default single-fluid)
 *   ISO_IONIZATION    Stage-E ionization/recombination network (default frozen)
 *   ISO_TRAC / ISO_TRAC_TCHROM / ISO_TRAC_TCMAXFRAC   TRAC TR-broadening
 *   ISO_QFLUX / _ENHANCE / _T_ON / _RAMP / _0   imposed coronal q(T) flux + ramp
 *   ISO_CHEAT / _E0 / _TMAX_MK / _SH_FRAC / _S0_MM / _ENHANCE / _T_ON / _RAMP
 *                     ambient volumetric coronal heating H(s) + ramp
 *   ISO_TBOOST / _TC / _W    one-time IC coronal superheat (tanh-in-T selector)
 *   ISO_REFINE_*      static local mesh refinement (scenarios/mesh.hpp)
 */

#pragma once

#include "chromosphere.hpp"

namespace chromosphere {

Vec  model_column_ic(Grid& grid);
void model_column_update_bc(Grid& grid, const Vec& xn);

} // namespace chromosphere
