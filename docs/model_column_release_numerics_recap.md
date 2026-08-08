# `model_column` release numerics

**Status: released.** This document describes the numerical method a normal `model_column`
run uses *today*. It is the reference for the active release configuration; the studies that
produced the decision are listed under "Evidence".

## Release configuration

A normal run needs **no** `ISO_RIEMANN` / `ISO_RECONSTRUCTION` setting:

```bash
scripts/run_chromo_realtime.sh outputs/model_column/<run>.txt \
    full no-ionization model_column - 20.0 no-cooling
```

(`CHROMO_T_END` sets the stop time; add `CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=<s>`
for snapshots. `GAMMA_TABLE` selects the production Gamma/Saha table; the `model_column`
scenario supplies it and the `ISO_*` model/grid defaults itself.)

| Element | Release setting |
| --- | --- |
| Model | 1D field-aligned hydrodynamics on a straight field line, gravity on, `B ≡ 1` |
| Thermodynamics | Gamma/Saha equilibrium closure with ionization energy, production `Gamma1` table |
| Reconstruction | MUSCL on `(ln ρ, V, ln p)`, MC3/Koren limiter, β = 2 |
| Face temperature | exact inversion of the same Saha closure, `T(ρ, p)` |
| Numerical flux | mixture **Roe characteristic** flux (3×3 local linearization) |
| Well balancing | equilibrium-reference (δ-form), plus inner discrete-HSE ghosts (p and ρ) |
| Conduction | physical only (no artificial mesh-scaled diffusivity) |
| Upper boundary | 22,000 K external conductive reservoir at the physical face |
| Mesh | coarse-equivalent N = 500 with R4 outer refinement → 661 actual cells |
| CFL | 0.50 (validated for this configuration; the hard-coded solver default remains 0.25) |

Both release choices are set in `model_column_ic` (`scenarios/model_column.cpp`) from
`gamma_mode`, so the code is the single source of truth. The `Grid` struct defaults stay
false, which leaves the non-Gamma scenarios on Rusanov + `(ln ρ, V, ln T)` unchanged.

## Reference (non-production) overrides

Retained for regression and controlled numerical comparison only. Neither is used in
production; both log a line naming the release default when set.

| Override | Effect | Why it is not the release |
| --- | --- | --- |
| `ISO_RIEMANN=rusanov` | local Lax–Friedrichs flux | its acoustic-scale dissipation `−½ a ΔU` with `a ~ c_s` is far too large where `|V| ≪ c_s`, leaving persistent cell-centred velocity/momentum ripple in the upper chromosphere/TR |
| `ISO_RECONSTRUCTION=lnrho-v-lnt` | limits `(ln ρ, V, ln T)` | density and temperature are limited independently and only then pushed through the nonlinear EOS, manufacturing face-pressure mismatch across the Saha transition |

Invalid values, and the Gamma-only choices requested off the Gamma/Saha path, throw.

## Evidence

- **Pressure reconstruction.** `(ln ρ, V, ln p)` cuts the reconstructed face-pressure
  mismatch by 2–3 orders of magnitude, and does so under *both* solvers — at 4000 s under
  Rusanov, `eps_p max` went 0.2326 % → 0.0015 % (N500) and 0.0875 % → 0.0008 % (N1000). It
  is a reconstruction improvement, not a Roe-specific patch.
  `docs/pressure_reconstruction_recap.md`; runs `outputs/model_column/acc_*`.
- **Solver.** With lnP, Rusanov at 4000 s still gives `Q_V ~ 0.081 / Q_M ~ 0.053` (N500) and
  `0.018 / 0.012` (N1000). Roe + lnP at N500/4000 s gives `eps_p max ~ 0.0023 %`,
  `Q_V ~ 0.0031`, `Q_M ~ 0.0019`, `V ∈ [+3, 11.6] m/s` (no reversal), `F_eff ~ 3.57e-10`.
  Run `outputs/model_column/lnp_roe_N500_4000s*`.
- **Cost.** The lnP path uses the optimized exact Saha pressure inversion (parent-cell
  temperature as Newton hint only — no table, no cache, no tolerance change), ~6.2 → ~3.1 EOS
  evaluations per face state. `docs/pressure_inversion_performance_recap.md`.
- **Promotion.** Running the new defaults with no environment variables reproduces an
  explicit `ISO_RIEMANN=roe-local ISO_RECONSTRUCTION=lnrho-v-lnp` run bitwise.
  Release smoke (N500, 100 s, defaults only): `eps_p max` 0.0003 %, `Q_V` 1.35e-3,
  `Q_M` 7.3e-4, `Q_eff` 7.3e-4, `F_eff` 1.09e-9 kg m⁻² s⁻¹, `T_top` 21838 K,
  `p_top` 0.010288 Pa, `n(V<0) = 0`, mass drift −3.0e-4.
  Run `outputs/model_column/release_roe_lnp_N500_100s*`.

Re-run both checks with `scripts/release_validation.sh`.

## Known limitation

The release fixes the **numerical method**, not the **resolution needed for a given
quantitative claim**. Long-duration evaporation mass flux still shows meaningful resolution
dependence: the Rusanov acceptance matrix showed a ~10–11 % N500 → N1000 gap at 4000 s. Roe
converged much better at 1000 s, but a complete Roe N1000 4000 s comparison has not been run
(`outputs/model_column/lnp_roe_N1000_partial2080s*` stopped at 2080 s).

**N500 long-duration evaporation mass flux is not established as < 10 % grid-converged.**
Do not describe N500 as grid-converged, and do not quote sub-10 % long-time mass-flux
accuracy at N500 without new evidence. This does not affect the numerical-method decision.
