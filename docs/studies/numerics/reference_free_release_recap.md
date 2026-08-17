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

Net 560× reduction. **All four rows are Roe-flux, float32-storage measurements** — see the re-measurement below, which supersedes the floor interpretation.

At the time this table was made the endpoint was comparable to the estimated float32 round-off scale of the stored state (1.32× at N=500), and doubling the resolution did not reduce it (N=1000 gave 2.66e-12, i.e. 2.63× the same estimate), so **no clean truncation-error trend survived over the resolutions tested**. That was deliberately weaker than claiming the residual *was* the representation limit: reconstruction, the EOS inversion, the nonuniform mesh and the boundary closures all contribute, and two resolutions could not separate them. The residual velocity moved off the lower boundary and onto the upper transition region.

### Re-measured after the double-precision release cutover (2026-08-17)

The float32-floor reading above was correct and is now obsolete. Re-measured in both precisions with the release Godunov flux, same configuration (N=500/R4, 661 cells, 20 s, `ISO_HEAT_FLUX=0`, CFL 0.50, `godunov.fallbacks=0`, `termination=end_time`). Runs under `outputs/model_column/precision_cutover/`. The extraction was first validated by reproducing the `1.33e-12 / 0.0302 / 2090 km / 1.01e-12` row above digit-for-digit from its own archived run `outputs/model_column/lowbc/H1_fix_wall.txt`; that also pins two things the original table did not state — the series is the conduction-off `H*` runs, and its "ε₃₂" is `FLT_EPSILON = 1.1921e-7` (machine epsilon, 2⁻²³), not the 5.96e-8 unit roundoff. The convention below is the original one, for continuity.

| Quantity, N=500/R4, Godunov, 20 s | float32 | **double (release)** | ratio |
|---|---|---|---|
| t=0 max\|F^ρ\| [kg m⁻² s⁻¹] | 1.0758e-12 | **1.7782e-13** | 6.05× smaller |
| — location of that max | 1601.1 km (base face) | **2149.1 km (upper TR)** | — |
| — t=0 base-face \|F^ρ\| | 1.0758e-12 | 1.2457e-13 | 8.6× smaller |
| — t=0 rms / median \|F^ρ\| | 3.163e-13 / 2.033e-13 | 1.991e-14 / 3.271e-16 | 16× / 620× |
| max\|V\| at 20 s [m/s] | 3.1943e-02 | **5.7274e-03** | 5.6× smaller |
| — location | 2093.9 km (upper TR) | 2153.0 km (**top cell**) | — |
| — max\|V\| excluding the top cell | 3.1943e-02 @ 2093.9 km | 4.9651e-03 @ 2149.1 km | 6.4× |
| floor ε·ρ·(\|v\|+c_s) | 1.0093e-12 | 1.8798e-21 | — |
| **defect / floor** | **1.07** | **9.5e+7** | — |

**Two results, and they matter differently.**

**(i) The float32 residual was the representation limit; the double residual is not.** float32 sat at 1.07× its own round-off floor (2.13× on the unit-roundoff convention). Double sits ~10⁸× above its floor, so the surviving residual has **no representation component at all** — it is entirely discretization and boundary closure.

**(ii) A clean second-order truncation trend now exists, and did not before.** Same refinement, N=500 → N=1000:

| Quantity | float32 | **double** |
|---|---|---|
| t=0 max\|F^ρ\| | 1.0758e-12 → 2.1518e-12 (apparent order **−1.0**) | 1.7782e-13 → **3.9112e-14** (order **2.19**) |
| t=0 base-face \|F^ρ\| | ratio 0.73 | ratio 3.86 (order **1.95**) |
| t=0 rms \|F^ρ\| | ratio 0.49 | ratio 6.64 (order **2.73**) |
| max\|V\| at 20 s | 2.4885e-02 (order 0.36) | **1.5291e-03** (order **1.90**) |

In float32 refining the mesh made the defect *worse*, which is the signature of per-face ULP noise, not truncation. In double the reference-free discretization converges at second order in all four measures. **This upgrades the reference-free claim**: the sentence "no clean truncation-error trend survives over the resolutions tested" was a statement about float32 storage, not about the scheme, and it no longer holds. The scheme is second-order convergent on the hydrostatic fixed point.

**What survives, and is real.** In double the base face is still dominated by the dissipative part of the flux — `f_central = -4.78e-15`, `f_diff = +1.2935e-13`, `f_total = 1.2457e-13` — so the remaining base defect is a Godunov dissipation response to a genuinely nonzero reconstructed left/right jump. That is a real discretization / inner-ghost-closure defect, not bit noise, and it is the same feature as the open first-interior-face artifact (`state_precision_release_cutover.md` §4). **Do not read "6× better" as "solved."** Likewise the double max\|V\| now sits on the *topmost* cell, and the top-cell velocity is nearly identical in both precisions (5.27e-3 vs 5.73e-3 m/s): that one is an outer-boundary-closure residual that precision does not touch, while the interior 2090 km peak collapsed.

**Not measured:** Roe-flux counterparts in double (the release flux only was re-run), and the 20 s relative density drift, which is not quotable from the main snapshot — that file carries six significant digits, giving a ~1e-6 relative measurement floor, and the double drift (9.5e-6) is only one decade above it.

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

5-step reference-free max\|V\| after `model_column_ic`, CFL 0.25, `ISO_HEAT_FLUX=0`. **The original ladder was float32 storage under the Roe flux**; re-measured after the double-precision cutover with a float32 control at the identical configuration, and extended one rung to `ISO_NS=2000` to test for flattening. `godunov.fallbacks=0` on every run.

| `ISO_NS` | cells | Δs [km] | **double** max\|V\| [m/s] | order | float32 max\|V\| [m/s] | order |
|---|---|---|---|---|---|---|
| 48 | 70 | 11.410 | 8.109e-1 | — | 8.098e-1 | — |
| 125 | 169 | 4.394 | 7.617e-2 | 2.48 | 7.702e-2 | 2.47 |
| 250 | 331 | 2.209 | 1.878e-2 | 2.04 | 1.790e-2 | 2.12 |
| 500 | 661 | 1.105 | **4.986e-3** | 1.91 | 5.845e-3 | 1.62 |
| 1000 | 1320 | 0.553 | 1.317e-3 | 1.92 | 3.381e-3 | **0.79** |
| 2000 | 2638 | 0.276 | 3.321e-4 | 1.99 | 4.590e-3 | **-0.44** |

The archived Roe/float32 ladder above reproduces digit-for-digit in the float32 build with `ISO_RIEMANN=roe-local` (`7.5762e-1, 7.5271e-2, 1.7748e-2, 5.8225e-3, 3.3808e-3`, and t=0 `max|F^ρ| = 1.3341e-12`), which is what makes the comparison trustworthy. The t=0 `max|F^ρ|` in double falls `3.32e-11 -> 1.31e-14` over the same ladder at orders 1.96–2.41, and the rms falls `8.62e-12 -> 5.46e-16` at orders 2.45–2.73. A double-precision **Roe** ladder gives the same orders (1.91–2.41), so second order is a property of the discretization, not of the Godunov flux.

**In double the bound is a clean, unbroken order ≈ 2 across the whole ladder, 48 -> 2000, with no flattening**, and max\|V\| falls monotonically `8.11e-1 -> 3.32e-4 m/s` over a 41x refinement. float32 tracks it only to `ISO_NS = 250`, degrades to order 0.79 by 1000, and **reverses** at 2000, bottoming out near `3.4e-3 m/s`. So the earlier reading — "the bound is resolution dependent, and that is the real cost of retiring `eq_wb`" — was measuring a **round-off floor, not the scheme**. The corrected statement: retiring `eq_wb` leaves an ordinary second-order truncation error, which is the expected and acceptable cost, and the release mesh's `4.99e-3 m/s` sits on that convergent trend rather than on a floor. `eq_wb` only ever made the *initial* state exact; it never removed the same truncation error from the evolving solution.

*(Do not cross-substitute: `4.986e-3 m/s` is the 5-step value at `ISO_NS=500`; the `5.727e-3 m/s` in the table above is the 20 s value. Different measurements.)*

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
