# `model_column` release numerics

> **Partly superseded — the numerical flux row is now historical.** The release face Riemann solver is the SWMF-style exact-Riemann Godunov flux at frozen composition (`docs/studies/numerics/swmf_godunov_flux_experiment.md` section 8); the mixture Roe characteristic flux described below is now the reference/comparison solver, reachable with `ISO_RIEMANN=roe-local`. Everything this document says about the Roe flux itself remains accurate — it is no longer the *release* flux. The reconstruction, conduction, mesh, boundary and limitation content is unaffected.
>
> **Also partly superseded (see `docs/studies/numerics/reference_free_release_recap.md`).** The release hydrodynamic operator has since gained the MUSCL-Hancock predictor momentum source and a trapezoidal EOS-closed lower ghost ladder, and **equilibrium-reference well balancing (`eq_wb` / `ISO_EQ_WB`) has been retired from the release** — it survives only in the two-fluid research solver. Statements below that describe `eq_wb` as active release method, or that quantify the discrete hydrostatic defect, are historical. Everything else in this document still stands.


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
| Solver / state | single-fluid equilibrium mixture, `U = (rho, rho u, E)` (`src/single_fluid/`); carrier quantities derived, never stored |
| Update path | `U^n -> MUSCL-Hancock/Godunov hydro -> U* -> implicit physical conduction -> U^{n+1}` (no projection stage) |
| Thermodynamics | Gamma/Saha equilibrium closure with ionization energy, production `Gamma1` table |
| Reconstruction | MUSCL on `(ln ρ, V, ln p)`, MC3/Koren limiter, β = 2 |
| Face temperature | exact inversion of the same Saha closure, `T(ρ, p)` |
| Numerical flux | SWMF-style **exact-Riemann Godunov** flux at frozen composition (`gamma = 5/3`, ionization energy in the passive offset `E0`), face-local Rusanov fallback. *(This row was `mixture Roe characteristic flux` when this document was written.)* |
| Well balancing | equilibrium-reference (δ-form), plus inner discrete-HSE ghosts (p and ρ) |
| Conduction | physical only, **structurally**: the operator has no TRAC factor and no artificial mesh-scaled diffusivity term at all |
| Upper boundary | 22,000 K external conductive reservoir at the physical face |
| Mesh | coarse-equivalent N = 500 with R4 outer refinement → 661 actual cells |
| CFL | 0.50 (validated for this configuration; the hard-coded solver default remains 0.25) |

Both release choices are set in `model_column_ic` (`scenarios/model_column.cpp`) from
`gamma_mode`, so the code is the single source of truth. The `Grid` struct defaults stay
false, which leaves the non-release two-fluid scenarios on Rusanov + `(ln ρ, V, ln T)` unchanged.

## Architecture

The release solver lives in `src/single_fluid/` (`mixture.hpp`, `mixture.cpp`,
`integrator.cpp`) and advances **three** conserved rows per cell. `x`, `rho_i = x·rho`,
`rho_n = (1-x)·rho`, `n_e`, `n_HI` and `p_e` are derived EOS diagnostics recomputed wherever
needed. `mixture_advance` contains exactly two stages — the explicit MUSCL-Hancock/Godunov hydro update
and the implicit physical conduction solve. There is no TRAC stage, no beam heating, no
volumetric coronal heating, no radiative cooling and no artificial conduction anywhere in
that path; those are two-fluid research physics, not release physics disabled by a flag.

The historical seven-row two-fluid / three-temperature / finite-rate-ionization solver is a
separate non-release path (`src/two_fluid/`: `two_fluid.hpp`, `state.cpp`, `flux.cpp`,
`rhs.cpp`, `integrators.cpp`) with its own state width, packing, boundaries and source
stages. The two directories share only `chromosphere.hpp` (the Grid) and `eos.hpp`; neither
includes or calls the other. A `Gamma1` table selects the release solver; each solver rejects
the other's `Grid`, and a release run rejects `SINGLE_FLUID` / `ENABLE_TE` / `ISO_TWO_FLUID`
/ `ISO_IONIZATION`.

## Solver identity of `model_column`

`model_column` names exactly one production model. `make_scenario` supplies the production
`Gamma1` table unconditionally, so the scenario always selects the release solver, and it
refuses `ISO_GAMMA` outright. `model_column_ic` additionally rejects every historical
two-fluid research knob — `ISO_GAMMA`, `ISO_TWO_FLUID`, `ISO_IONIZATION`, `ISO_COOLING`,
`ISO_TRAC`, `ISO_CORONA`, `ISO_CHEAT`, `ISO_QFLUX`, `ISO_TBOOST`,
`ISO_NUMERICAL_DIFFUSIVITY_MULT` — rather than silently accepting them. The retired
`model_isentropic` alias has been removed. The historical two-fluid column is reached through
the separate `model_gentle` scenario, which never loads a `Gamma1` table.

The former conservative equilibrium projection has been **removed**: it was the identity on
`(rho, rho·u, E)` and existed only to repair the four redundant degrees of freedom of the
carrier-row layout. Removing it also removes one caloric inversion per cell per step and the
small spurious per-step drift heating it used to thermalize.

## Reference (non-production) overrides

Retained for regression and controlled numerical comparison only. None is used in
production; each logs a line naming the release default when set.

| Override | Effect | Why it is not the release |
| --- | --- | --- |
| `ISO_RIEMANN=roe-local` | mixture Roe characteristic flux (3×3 equilibrium linearization, acoustic eigenvalues carrying `Gamma1`) | it is the *equilibrium* acoustic closure; the release takes the frozen-composition limit instead (`docs/studies/numerics/swmf_godunov_flux_experiment.md` §8.1). Selecting it also reverts the CFL to the equilibrium signal speed, so the comparison mode stays self-consistent |
| `ISO_RIEMANN=rusanov` | local Lax–Friedrichs flux | its acoustic-scale dissipation `−½ a ΔU` with `a ~ c_s` is far too large where `|V| ≪ c_s`, leaving persistent cell-centred velocity/momentum ripple in the upper chromosphere/TR |
| `ISO_RECONSTRUCTION=lnrho-v-lnt` | limits `(ln ρ, V, ln T)` | density and temperature are limited independently and only then pushed through the nonlinear EOS, manufacturing face-pressure mismatch across the Saha transition |

Invalid values, and the Gamma-only choices requested off the Gamma/Saha path, throw.

## Evidence

- **Pressure reconstruction.** `(ln ρ, V, ln p)` cuts the reconstructed face-pressure
  mismatch by 2–3 orders of magnitude, and does so under *both* solvers — at 4000 s under
  Rusanov, `eps_p max` went 0.2326 % → 0.0015 % (N500) and 0.0875 % → 0.0008 % (N1000). It
  is a reconstruction improvement, not a Roe-specific patch.
  `docs/studies/numerics/pressure_reconstruction_recap.md`; runs `outputs/model_column/acc_*`.
- **Solver.** With lnP, Rusanov at 4000 s still gives `Q_V ~ 0.081 / Q_M ~ 0.053` (N500) and
  `0.018 / 0.012` (N1000). Roe + lnP at N500/4000 s gives `eps_p max ~ 0.0023 %`,
  `Q_V ~ 0.0031`, `Q_M ~ 0.0019`, `V ∈ [+3, 11.6] m/s` (no reversal), `F_eff ~ 3.57e-10`.
  Run `outputs/model_column/lnp_roe_N500_4000s*`.
- **Cost.** The lnP path uses the optimized exact Saha pressure inversion (parent-cell
  temperature as Newton hint only — no table, no cache, no tolerance change), ~6.2 → ~3.1 EOS
  evaluations per face state. `docs/studies/performance/pressure_inversion_performance_recap.md`.
- **Promotion.** Running the new defaults with no environment variables reproduces an
  explicit `ISO_RIEMANN=roe-local ISO_RECONSTRUCTION=lnrho-v-lnp` run bitwise (the
  equivalence leg now pins `ISO_RIEMANN=swmf-godunov` instead).
  Release smoke (N500, 100 s, defaults only): `eps_p max` 0.0003 %, `Q_V` 1.35e-3,
  `Q_M` 7.3e-4, `Q_eff` 7.3e-4, `F_eff` 1.09e-9 kg m⁻² s⁻¹, `T_top` 21838 K,
  `p_top` 0.010288 Pa, `n(V<0) = 0`, mass drift −3.0e-4.
  Run `outputs/model_column/release_roe_lnp_N500_100s*` (the release smoke artefact is
  now `release_godunov_lnp_N500_100s*`).
- **Three-variable architecture.** After the state refactor the same smoke test gives
  `eps_p max` 0.0002 %, `Q_V` 1.32e-3, `Q_M` 7.0e-4, `Q_eff` 6.88e-4,
  `F_eff` 1.11e-9 kg m⁻² s⁻¹, `T_top` 21838.4 K, `p_top` 0.010288 Pa, `n(V<0) = 0`,
  mass drift −2.84e-4, and the defaults-vs-explicit-override run is still bitwise identical.
  A 60 s N500 comparison against the pre-refactor binary takes the identical step count
  (10668) and end time; T and p agree to a median 1e-5 / max 1.4e-3 relative, the peak
  evaporation velocity moves 40.97 → 40.51 m/s, and column mass agrees to 1.7e-5. That
  residual is accumulated float32 representation round-off, not a systematic error: it does
  not shrink with `dt` (it grows with step count at CFL 0.25) and it is **smaller** than the
  solver's own CFL 0.50 → 0.25 discretization sensitivity (median |Δv| 0.26 vs 0.49 m/s).

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
