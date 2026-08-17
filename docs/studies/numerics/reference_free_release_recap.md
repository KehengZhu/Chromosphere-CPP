# Reference-free `model_column` release: MUSCL-Hancock source term, lower truncation, retirement of `eq_wb`

**Status: GO — reference-free release.** "Reference-free" means precisely: *the interior hydrodynamic operator needs no frozen global hydrostatic reference state.* The boundary closures still capture fixed reservoir data once from the IC (`kOuterPRef`, `kInnerTRef`/`kInnerPRes0/1`), as any truncated-domain problem must. Supersedes the well-balancing sections of `docs/studies/numerics/model_column_release_numerics_recap.md`, `docs/studies/evaporation/gentle_evaporation_downflow.md` and `docs/studies/boundaries/boundary_conditions_plan.md` for current release behaviour.

## What changed

Three release changes, in dependency order.

1. **`src/single_fluid/mixture.cpp` — the MUSCL-Hancock predictor now applies the momentum source.** `build_source` is hoisted above the predictor loop and its result is added as `dt_predictor*source(k)`. The predictor differences cell *i*'s own two extrapolated fluxes over `ds_i` — that is cell *i*'s pressure gradient — so it must apply the term that balances it. Omitting it left the half state of a hydrostatic column with a spurious `dt*|g|` velocity that the corrector then reconstructed into every face. For `model_column` (`B ≡ 1`) the source is gravity alone; the flux-tube term is retained for general `B(s)`.

2. **`scenarios/model_column.cpp` — the lower ghost ladder is second order and EOS-closed.** `gamma_hse_rung_down` replaces `p_0 + k·ρ_0·g·Δs` with the trapezoidal step `p_{k+1} = p_k + Δs·½(ρ_k+ρ_{k+1})·g`, `ρ_{k+1} = ρ_EOS(p_{k+1}, T_0)`, iterated to a fixed point — the same construction the outer face already used. The old form carried an O((Δs/H)²) ghost pressure error, i.e. an O(Δs/H) force-balance error in the lowest interior cells.

3. **`eq_wb` retired from the release.** `model_column_ic` sets `grid.eq_wb = false` on the gamma (release) path; `ISO_EQ_WB` is gone; `mixture_rhs_explicit` no longer subtracts `grid.eq_residual`; `ensure_eq_residual` is removed from `src/single_fluid/integrator.cpp`; the dead `eq_residual_mass` face-flux sidecar column and its consumers are removed. **`Grid::eq_wb/eq_state/eq_residual` remain for the two-fluid research solver**, which never received change 1 and still needs them — `model_column_ic` keeps `eq_wb = true` on its non-gamma leg so legacy two-fluid behaviour is byte-preserved.

Terminology, in code comments and both writeups: **1600 km is the lower mid-chromospheric truncation of the modeled domain, not the photosphere.** It is the height above which the Saha/LTE closure is intended to be used; the two inner ghosts are a MUSCL stencil closure, not a claim that LTE extends below.

## Evidence

Runs under `outputs/model_column/lowbc/`; the pre-change baselines are `outputs/model_column/mtr_proto/`. Release config throughout unless stated: N=500/R4, 661 cells, CFL 0.50 via `scripts/run_chromo_realtime.sh` (CFL 0.25 for the step-capped scans).

### Discrete hydrostatic defect, no reference subtraction (N=500, 20 s, `ISO_HEAT_FLUX=0`)

| Configuration | max\|F^ρ\| at t=0 [kg m⁻² s⁻¹] | max\|V\| at 20 s [m/s] | location of max\|V\| |
|---|---|---|---|
| previous release | 7.47e-10 | 0.710 | 1647 km (base) |
| + predictor source | 1.43e-11 | 0.0325 | 1608 km (base) |
| + trapezoidal lower ghosts | **1.33e-12** | **0.0302** | 2090 km (upper TR) |
| float32 floor ε₃₂·ρ·(\|v\|+c_s) | 1.01e-12 | — | — |

Net 560× reduction. The endpoint is comparable to the estimated float32 round-off scale of the stored state (1.32× at N=500), and doubling the resolution does not reduce it (N=1000 gives 2.66e-12, i.e. 2.63× the same estimate), so **no clean truncation-error trend survives over the resolutions tested**. That is deliberately weaker than claiming the residual *is* the representation limit: reconstruction, the EOS inversion, the nonuniform mesh and the boundary closures all contribute at this level, and two resolutions cannot separate them. The residual velocity moved off the lower boundary and onto the upper transition region.

### Reference-free vs frozen reference, conduction-driven

Analysis window 2130–2150 km, mean V:

| Duration / mesh | eq_wb on | reference-free | difference |
|---|---|---|---|
| 100 s, N=500 | 14.961 m/s | 14.957 m/s | 0.03 % |
| 100 s, N=1000 | 16.314 m/s | 16.338 m/s | 0.15 % |
| 1000 s, N=500 | 6.900 m/s | 6.925 m/s | 0.37 % |

1000 s peak velocity 16.256 vs 16.333 m/s (0.47 %); column mass drift −1.955e-4 vs −1.972e-4 (0.9 %). All runs `termination=end_time`. For contrast, dropping `eq_wb` **without** the predictor fix cost 17 % of the 100 s window velocity (12.240 m/s) — the frozen reference had been masking the defect, not fixing it.

The 1000 s comparison is the substantive reference-free acceptance test: by then the transition region has restructured well away from the frozen C7 initial condition, so `R(U_eq)` no longer describes the state it is subtracted from, yet the two solutions agree to 0.4 %.

### Equilibrium not identical to the frozen C7 reference

Off-canonical domain (`ISO_H_BASE=1700 ISO_DH=453`, a different stratification and mesh), 20 s conduction-off, reference-free: max\|V\| = 0.071 m/s, relative density drift 5.2e-5, t=0 defect 3.95e-12. The discretization maintains a quasi-hydrostatic chromosphere for a stratification it was not tuned on.

### Resolution

5-step reference-free max\|V\| after `model_column_ic`, CFL 0.25:

| `ISO_NS` | cells | max\|V\| [m/s] |
|---|---|---|
| 48 | 70 | 7.58e-1 |
| 125 | 169 | 7.53e-2 |
| 250 | 331 | 1.77e-2 |
| 500 | 661 | **5.82e-3** |
| 1000 | 1320 | 3.38e-3 |

The bound is resolution dependent — that is the real cost of retiring `eq_wb`, which made any chosen state exact at any resolution. It is acceptable because the release ships one mesh, and because `eq_wb` only ever made the *initial* state exact; it never removed the same truncation error from the evolving solution.

### Lower-boundary acoustic reflection (investigated, not promoted)

A gravity-aware characteristic lower boundary was prototyped: impose only the incoming (upward, `u+c`) acoustic invariant `δu + δp/Z = 0` about the base reference state with `Z = ρc`, take the outgoing invariant from the live base cell, giving `V_g = −(p_0 − p_0,ref)/Z`, `p_g = p_g,ref − Z·V_0`. It degenerates exactly to the fixed reservoir at equilibrium.

Measured against the reservoir wall: identical t=0 defect (1.33e-12) and identical 20 s residual profile; 100 s window velocity 14.951 vs 14.957 m/s; 1000 s 6.896 vs 6.925 m/s (0.4 %). Reflection probe (300 s conduction-on, 1 s cadence, mean V over 1600–1700 km, detrended for t > 150 s): 3 sign changes and 0.31 m/s peak-to-peak **either way**. The disturbance the TR sends down arrives as a slow quasi-static pressure adjustment, not acoustic ringing, so this face is not an important reflector for this model. **Not promoted; the code was removed rather than left as a default-off switch.**

**Standing limitation — do not read the null result too broadly.** The reservoir wall means mass cannot cross 1600 km, so the column cannot be replenished from the unresolved chromosphere below. The characteristic variant recovered **~17 % of the 1000 s column mass loss** that way while changing the evaporation velocity by <0.4 %. So the lower closure is nearly irrelevant to the velocity field and to the oscillation content, but it is *not* irrelevant to the column mass budget — which matters directly for a coronal mass-loading interpretation. This is a modeling choice of the truncation, not a numerical defect, and it must be revisited if base mass supply becomes a quantity of interest.

### Rejected: full `ModTransitionRegion` pressure/gravity split

Taking pressure out of the momentum flux and differencing the acoustic star pressure explicitly alongside gravity reproduced the predictor-source result and nothing beyond it (t=0 defect 1.43e-11, identical to predictor-source-only; 100 s window velocity 14.647 vs 14.938 m/s). Rejected and removed on that evidence alone.

A conservation caveat is worth stating precisely rather than overstating: for the constant-area release geometry (`B ≡ 1`) the shared face pressure telescopes and the momentum update remains conservative, so conservation is **not** a reason to reject the split for `model_column`. The concern is that the source-split form is not generally equivalent to the conservative area-weighted operator on a variable-area flux tube, which matters only if this solver is later run on non-constant `B(s)`.

## Validation

- `build_omp/chromo_tests`: **11362/11362**, 0 failed.
- `ctest`: **9/9**.
- `scripts/release_validation.sh`: default-path equivalence identical (snapshots and `.gamma_diag`); 100 s smoke `termination=end_time`, mean V 14.956 m/s, `n(V<0)=0`, mass drift −2.858e-4, `q_phys` 0.7783 W/m².
- `test_stage8_gamma_model_column_saha_hse_and_ghosts` now runs its 5-step hydrostatic regression through `make_scenario("model_column")` on the **actual release preset** (661 cells, 1600–2153 km, 22 kK wall) at the **production CFL 0.50** with `ISO_HEAT_FLUX=0`, asserts `eq_wb == false` and an empty `eq_residual`, prints the measured values, and bounds them at 3–5× headroom: max\|V\| = 5.9e-3 m/s, max\|ρv\| = 5.4e-13, ‖ΔU‖/‖U‖ = 9.2e-8 (5.8e-3 / 4.2e-13 / 5.1e-8 at CFL 0.25 — the residual is a quasi-steady balance, not simply ∝ dt). It previously ran on a 48-cell grid with library defaults and was satisfiable only because `eq_wb` made it identically zero.

## Performance

RHS stage 4.75 → 4.91 s per 100 physical seconds (+3.4 %); decode and conduction unchanged; conduction is ~57 % of the step, so total step cost moves within run-to-run noise. The one-time `R(U_eq)` evaluation is saved. Per-step scaling is unchanged at ~0.114 s wall per physical second at 12 threads.

## Not changed in this work

Saha EOS, the face Riemann solver, reconstruction variables, MC3 limiter (β=2), mesh strategy, conduction physics, CFL, floors, and the hydro/conduction operator splitting. (The Riemann solver was the mixture Roe flux throughout this work; the release flux has since become the SWMF exact-Riemann Godunov flux — `docs/studies/numerics/swmf_godunov_flux_experiment.md`. Every hydrostatic-defect number below is a Roe measurement; the Godunov counterpart of the `1.33e-12` figure is `1.08e-12`, at the same float32 floor.) The splitting remains **first-order Lie** (explicit hydro → implicit backward-Euler conduction).

**No second-order temporal accuracy is claimed, for the hydrodynamic stage or for the composition.** MUSCL-Hancock gives second-order space/time centering for the *homogeneous* flux evolution, and the predictor now carries the source so hydrostatic balance is consistent across stages — but the corrector still evaluates the momentum source at `U^n`, so in the source-only limit `dU/dt = S(U)` the update degenerates to forward Euler. Two separate future numerical questions, neither addressed here: raising the order of the source quadrature (e.g. evaluating `S` at the half-time state), and symmetrizing the hydro/conduction composition (Strang).

## Writeup

`docs/writeup-overleaf/paper.tex` and `main.tex` were aligned: MUSCL-Hancock named and its predictor equation given with provenance (van Leer 1979; Toro 2009 ch. 14; Berthon 2006; Koren 1993 for the limiter), the well-balancing section replaced by the reference-free account with the a-priori/no-a-priori framing (Käppeli & Mishra 2014; Chandrashekar & Klingenberg 2015; Berberich et al. 2018), the 1600 km boundary re-described as a mid-chromospheric truncation, and the Lie splitting and its first-order composite accuracy stated explicitly.
