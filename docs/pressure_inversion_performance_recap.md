# Optimizing the Saha pressure-face temperature inversion

> **Superseded by the model_column release-numerics promotion.** Roe (`roe-local`) and
> `(ln rho, V, ln p)` are now the `model_column` release defaults and need no environment
> override; `ISO_RIEMANN=rusanov` and `ISO_RECONSTRUCTION=lnrho-v-lnt` are reference-only.
> Statements below describe the configuration in force **when this study was run** and are
> kept for chronology. Current release description: `docs/model_column_release_numerics_recap.md`.

**Status: GO.** Pure Newton initial-guess optimization. No physics, no reconstruction variables, no tolerance and no release default changed.

## Problem

The `(ln rho, V, ln p)` reconstruction (`ISO_RECONSTRUCTION=lnrho-v-lnp`, commit `93c32d8`) removed the late coarse-grid Roe ripple, but it made the explicit RHS about 70 % more expensive than the `(ln rho, V, ln T)` path.

The cost is the face closure. The temperature builder evaluates the Saha closure once per face state; the pressure builder must first recover `T` from the reconstructed `(rho, p)` by a safeguarded Newton/bisection iteration, and **every Newton step costs one Saha evaluation**. (`invert_pressure_eval` returns the converged `CaloricEval`, so no evaluation is repeated *after* convergence — but that only removes the last one, not the iteration.) The builder passed no initial guess, so the solver started from the generic bracket midpoint `0.75 T_neutral` and consumed ~6.2 evaluations per face.

## What was changed

Every reconstructed one-sided face state is a small MUSCL extrapolation from one specific cell, whose decoded temperature the RHS has already computed. That temperature is now passed to the inversion as its starting point.

- `DecodedMixtureField::extended_temperature` (length `ns+4`, ghosts included) stores the decoded `T` of every extended cell, written by the same `put()` lambda that writes `extended_primitive` in both decode sites (`decode_mixture_field_into` and the predictor `decode_predicted_caloric_primitives_into`), and by the same ghost path. It is **not** a reconstruction variable: it never enters the limiter, the extrapolation, or any flux. The MUSCL state stays exactly `(ln rho, V, ln p)`.
- `temperature_stencil_view()` produces the shifted parent-cell hint vectors, using the same offset arithmetic as `mixture_stencil_view`. Parentage matches the reconstructions in `rhs_explicit_mixture` exactly:

  | face state | reconstructed from | hint |
  |---|---|---|
  | predictor `wl_iph[i]`, `wr_imh[i]` | cell `i` | `T_i` |
  | corrector `wr_iph[i]` | cell `i+1` | `T_{i+1}` |
  | corrector `wl_iph[i]`, `wr_imh[i]` | cell `i` | `T_i` |
  | corrector `wl_imh[i]` | cell `i-1` | `T_{i-1}` |

  Hints for `i=0` / `i=ns-1` come from the decoded ghost layers, through the same extended array.
- `equilibrium_mixture_face_state_from_log_pressure` gained an optional trailing `temperature_guess` (default NaN, so every other caller is unchanged). The hints are only built when `grid.pressure_reconstruct` is set; the `ln T` path is untouched.
- In `invert_pressure_eval`, a finite hint is **clamped** into the exact bracket `(T_neutral/2, T_neutral]` rather than discarded. Both ends matter: where `x -> 0` the root sits at `T_neutral` and where `x -> 1` it sits at `T_neutral/2`, so a neighbouring cell's temperature falls marginally outside the face's own bracket precisely where it is the most accurate estimate. Clamping raised the accepted-hint rate from 97.1 % to 100.0 %.

The bracket, the safeguarded bisection fallback, the `4 eps p` residual tolerance, the fail-loud non-convergence throw and the Saha closure are all unchanged. The hint can only change *where Newton starts*.

## Rejected

| idea | measured | why rejected |
|---|---|---|
| relax residual tolerance to `1e-11 p` | 3.22 -> 3.06 evaluations, RHS 4.87 -> 4.75 s | ~2 % of RHS for a change that requires full numerical revalidation of an already-validated result. Not worth it. |
| second-order Taylor hint from parent `x` and the analytic Saha derivatives | not implemented; bounded above by ~2.2 evaluations from the measured marginal cost of 0.56 s per evaluation-per-face | best case ~4 % of total runtime, in exchange for a new per-cell `x_eq` array, three more shifted views and a derivative expression duplicated outside the EOS. Fails the cost/benefit test once the acceptance target is already met. |
| predictor-face / previous-step / previous-timestep temperature cache (Phase 3) | not implemented | the parent-cell hint already reaches ~3.2 evaluations against a hard floor of ~3 at machine-precision tolerance (one evaluation is unavoidably spent *confirming* convergence). Temporal caching adds invalidation, OpenMP and restart-determinism surface for at most ~1 evaluation. |

## Pressure-inversion statistics

N500 Roe lnP, 4000 s, 712,314 steps (`CHROMO_PROFILE=1`):

| | before | after |
|---|---|---|
| calls | — (uninstrumented) | 2,825,753,187 |
| mean evaluations per call | ~6.2 | **3.156** |
| maximum evaluations | 9 | 10 |
| bisection fallbacks | — | **0** |
| hint accepted | n/a | 2,825,041,290 (99.97 %) |

The remaining unhinted calls are the small number where the reconstructed `(rho, p)` is far enough from the parent cell that the parent temperature is not merely at the bracket edge.

Roughly **3.0 Saha evaluations per face state were removed**, i.e. about 8.6e9 evaluations over the 4000 s run.

New profile lines: `pressure_inversion.calls`, `.hinted_calls`, `.average_evaluations`, `.maximum_evaluations`, `.bisection_fallbacks`. They are separate from `inversion.*`, which counts the `(rho,E) -> T` energy inversion of the decode/projection stages, and cost one predictable false branch when profiling is off.

## Performance

Alternating 3-repetition benchmark, N500 Roe, 100 s, low-I/O (`scripts/run_chromo_realtime.sh` defaults), 12 threads, mean seconds:

| | `ln T` | `ln p` before | `ln p` after |
|---|---|---|---|
| `time.rhs_s` | 3.680 | 6.211 | **4.910** |
| `time.decode_s` | 0.884 | 0.881 | 0.880 |
| `time.projection_s` | 0.792 | 0.787 | 0.775 |
| `time.conduction_s` | 7.934 | 7.950 | 7.899 |
| total profiled | 13.377 | 15.916 | **14.551** |

- RHS: **-21.0 %**; overhead relative to `ln T` **68.8 % -> 33.4 %**.
- Total profiled: **-8.6 %**; overhead relative to `ln T` **19.0 % -> 8.8 %**.

Production-scale confirmation, N500 Roe lnP 4000 s with the documented diagnostic I/O:

| | before | after |
|---|---|---|
| `time.rhs_s` | 257.43 | **191.80** (-25.5 %) |
| `time.decode_s` | 33.32 | 33.47 |
| `time.projection_s` | 31.47 | 31.28 |
| `time.conduction_s` | 252.48 | 261.98 |

Only RHS moved; the conduction difference is run-to-run variance in a stage the change does not touch.

## Numerical parity

**Tests.** `build` and `build_omp`: 11,663/11,663 unit checks and 8/8 ctest cases pass. New `test_pressure_face_hint_does_not_change_face_state` sweeps `rho` in 1e-12..1e-7 kg/m^3 and `T` in 4.5e3..2.2e6 K (and the scalar inversion over `n_H` 1e12..1e26 m^-3, `T` 3.2e3..1e7 K) against 16 hints — exact, perturbed, both bracket ends, out-of-bracket, zero, negative, NaN, +/-Inf — and requires `T`, `x_eq`, `p`, `e_int`, `Gamma1`, `c`, the conserved face rows and the physical face flux to agree with the unhinted result to `< 1e-13` relative.

**Bit-level run parity.** Baseline (`93c32d8`, built in a clean worktree) and optimized binaries, identical configuration, both individually bit-reproducible run-to-run: outputs are **bit-identical through t = 20 s (3,556 steps)** and diverge only afterwards. At 100 s the divergence is 1.5e-5 of the column maximum in density and energy and 5e-3 in momentum — the smallest, most sensitive field, where N500-to-N1000 discretization error is ~3 %.

**Production diagnostics**, N500 Roe lnP (`util/pressure_reconstruction_diag.py`), validated baseline vs optimized:

| | 100 s base / opt | 1000 s base / opt | 4000 s base / opt |
|---|---|---|---|
| `eps_p` max | 0.0003 % / 0.0002 % | 0.0020 % / 0.0020 % | 0.0023 % / 0.0023 % |
| `Q_V` | 0.00133 / 0.00135 | 0.00285 / 0.00285 | 0.00311 / 0.00315 |
| `Q_M` | 0.00070 / 0.00073 | 0.00152 / 0.00150 | 0.00185 / 0.00192 |
| `Q_eff` | 7.17e-04 / 7.63e-04 | 6.62e-04 / 5.63e-04 | 5.73e-04 / 8.46e-04 |
| min V [m/s] | 11.071 / 11.053 | 3.676 / 3.670 | 2.969 / 2.998 |
| mean V [m/s] | 14.669 / 14.652 | 6.639 / 6.621 | 5.693 / 5.714 |
| n(V<0) | 0 / 0 | 0 / 0 | 0 / 0 |
| `F_eff` [kg m^-2 s^-1] | 1.0926e-09 / 1.0917e-09 | 4.2762e-10 / 4.2667e-10 | 3.5708e-10 / 3.5868e-10 |
| `T_top` [K] | 21838.2 / 21838.2 | 21904.7 / 21904.6 | 21912.8 / 21912.8 |
| `p_top` [Pa] | 0.010288 / 0.010288 | 0.0102878 / 0.0102878 | 0.0102878 / 0.0102878 |
| BL thickness [km] | 11.870 / 11.870 | 11.870 / 11.870 | 12.140 / 12.140 |
| `q_phys` [W m^-2] | 0.7747 / 0.7747 | 0.4563 / 0.4566 | 0.4173 / 0.4173 |
| cumulative [J m^-2] | 83.14 / 83.14 | 569.21 / 569.20 | 1853.83 / 1853.73 |
| mass drift | -2.987e-04 / -2.988e-04 | -1.756e-04 / -1.757e-04 | -1.632e-04 / -1.723e-04 |

Every conservative, thermal and conductive observable matches the validated result to 4 significant figures. `Q_V`/`Q_M`/`Q_eff` differ by a few percent of themselves in both directions, at the level expected from a round-off-seeded trajectory difference over 712k steps, and remain two orders below the `ln T` values the reconstruction was introduced to remove.

## Reproduction

```bash
ISO_NS=500 ISO_RIEMANN=roe-local ISO_RECONSTRUCTION=lnrho-v-lnp \
CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FACE_FLUX_DIAG=1 \
CHROMO_OUTER_COND_DIAG=1 CHROMO_FRAME_DT=20 CHROMO_FACE_FLUX_STRIDE=2000 \
CHROMO_OUTER_COND_STRIDE=400 \
scripts/run_chromo_realtime.sh outputs/model_column/lnp_roe_N500_4000s_opt.txt \
    full no-ionization model_column - 20.0 no-cooling

.venv/bin/python util/pressure_reconstruction_diag.py --targets 100,500,1000,4000 \
  --run "lnP-base=outputs/model_column/lnp_roe_N500_4000s.txt" \
  --run "lnP-opt=outputs/model_column/lnp_roe_N500_4000s_opt.txt"
```

## Limitations

The pressure reconstruction remains diagnostic (`ISO_RECONSTRUCTION` default is `lnrho-v-lnt`, `pressure_reconstruct=false`), and this task does not change that. A residual ~9 % total-runtime overhead against `ln T` remains and is close to the achievable floor at machine-precision inversion tolerance: about one of the ~3.2 evaluations is spent confirming convergence and cannot be removed by a better guess. Going materially below it would require relaxing the inversion tolerance or a table-seeded guess, neither of which is justified by the measured gain.
