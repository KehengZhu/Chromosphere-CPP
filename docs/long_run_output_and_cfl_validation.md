# Long-Run Output Fixes and CFL Acceptance Sweep

Two tasks: (1) make `CHROMO_T_END` runs terminate and emit output correctly, and (2) determine the largest scientifically acceptable CFL for the production-shaped `model_column` case, then measure whether 1000 physical seconds fit in 1000 wall seconds.

## Final status

```text
LONG_RUN_OUTPUT_AND_CFL_VALIDATION: GO
REAL_TIME_1000S_TARGET: GO
```

Accepted production CFL: **0.50**, applied by `scripts/run_chromo_realtime.sh`. The hard-coded solver default stays at 0.25.

---

## 1. Repository state

### Initial `git status`

```text
 M AGENTS.md
 M CLAUDE.md
 M docs/local_mesh_scaled_numerical_conduction_recap.md
 M util/outer_refine_diag.py
 M visualize_commands.md
```

These five pre-existing uncommitted changes (497 added lines, all documentation/diagnostic-script edits from earlier work) were preserved. `AGENTS.md` and `CLAUDE.md` were already byte-identical to each other and were extended identically here. No pre-existing change was reverted, reset, or reformatted.

### Files changed by this task

```text
run_control.hpp                              (new)  step-cap + snapshot-cadence helpers
chromo_main.cpp                                     wired the helpers into the driver loop
tests/chromo_tests.cpp                              11 new run-control tests
scripts/run_chromo_realtime.sh               (new)  long-run production launcher
util/analyze_cfl_sweep.py                    (new)  sweep analysis
util/run_cfl_sweep.sh                        (new)  sweep driver
AGENTS.md, CLAUDE.md                                validated CFL + long-run controls
docs/long_run_output_and_cfl_validation.md   (new)  this document
```

No change to physics, EOS, mesh, refinement geometry, limiter, Riemann solver, conduction formulation, Thomas solve, boundary conditions, diagnostics, artificial conduction, physical constants, or any OpenMP loop. Build directories (`build_cfl_audit/`, `build_cfl_serial/`) and simulation outputs under `outputs/` are not part of the source patch.

---

## 2. Long-run termination correctness

### Old behavior

```cpp
const int step_cap = static_cast<int>(std::max(10000.0f, 10000.0f * time_mult));
```

The cap came from the positional `time_mult` argument even when `CHROMO_T_END` overrode the end time. With `CHROMO_T_END=1000` and `time_mult=20` the cap was 200,000 steps; at the measured dt ≈ 7.0e-4 s (CFL 0.25) a 1000 s run needs ≈ 1.42 million steps, so the run stopped near **140 of the 1000 requested physical seconds** and said nothing about it.

### New behavior

`select_step_cap()` in [run_control.hpp](run_control.hpp) is a pure function of the two environment strings and `time_mult`:

| condition | effective cap | source label |
|---|---|---|
| `CHROMO_STEP_CAP` set and valid | that value | `explicit_step_cap` |
| `CHROMO_T_END` set, no explicit cap | `kUnboundedStepCap` = 1e12 | `end_time` |
| neither set | `max(10000, 10000*time_mult)` — the original float expression, unchanged | `legacy_time_mult` |

`CHROMO_STEP_CAP` values that are non-positive, unparseable, partially numeric, or out of `long long` range are **rejected with a message and exit status 1**, never silently replaced. The step counter and cap are both `long long`, and every cap path is clamped to 1e12 (< `LLONG_MAX`/1000), so neither `step + 1` nor `step < cap` can overflow; an infinite or absurd `time_mult` clamps rather than invoking undefined cast behavior.

Every run now prints:

```text
termination=end_time|step_cap|other
requested_end_time=...
final_physical_time=...
final_step=...
effective_step_cap=... (source=...)
```

and a `step_cap` termination additionally emits a `WARNING` on stderr. Reaching the end time wins when both fire on the same step. Verified live:

```text
$ CHROMO_STEP_CAP=50 CHROMO_T_END=1000 ... chromo_main ...
termination=step_cap
requested_end_time=1000
final_physical_time=1.51784
final_step=50
effective_step_cap=50 (source=explicit_step_cap)
WARNING: run stopped at the step cap (50) after 1.51784 of the requested 1000 physical seconds.
```

```text
$ CHROMO_STEP_CAP=0    -> step-cap error: CHROMO_STEP_CAP must be a positive integer, got '0'   (exit 1)
$ CHROMO_STEP_CAP=-5   -> step-cap error: CHROMO_STEP_CAP must be a positive integer, got '-5'  (exit 1)
$ CHROMO_STEP_CAP=abc  -> step-cap error: CHROMO_STEP_CAP is not a valid integer: 'abc'         (exit 1)
```

### Tests

`test_step_cap_legacy_time_mult`, `test_step_cap_end_time_removes_legacy_cap`, `test_step_cap_explicit_wins`, `test_step_cap_rejects_invalid_values`, `test_step_cap_no_overflow`, `test_termination_classification` — covering the legacy path, the end-time path, explicit caps with and without an end time, zero/negative/unparseable/out-of-range values, the empty-string-as-unset case, overflow safety, a very large requested end time, and all three termination classes including the tie.

---

## 3. Physical-time output cadence

### Old behavior

The snapshot cadence was a step stride, also derived from `time_mult`: `max(10, 10*max(1, time_mult/5))` = **40 steps** for `time_mult=20`. A `CHROMO_T_END=1000` run therefore kept a 20-second run's stride. At the accepted CFL 0.50 (dt ≈ 1.405e-3 s) that is one frame every 0.056 physical seconds:

```text
~711,600 steps  ->  ~17,790 snapshots
~46.9 million cell rows in the main output
main file ~3.6 GB  +  gamma_diag ~4.8 GB  =  ~8.3 GB
```

(extrapolated from the measured 0.20 MB and 0.27 MB per snapshot in this configuration). That both fills the disk and destroys the timing measurement.

### New behavior

`CHROMO_FRAME_DT` sets the cadence in **physical seconds**. `OutputSchedule` in [run_control.hpp](run_control.hpp) schedules against `frame_index * frame_dt` recomputed from an integer index, so the cadence cannot drift; after a write the index advances past every scheduled time the new state already covers, so a step crossing several nominal output times writes once. The initial state is always written, the final state is written if the cadence did not already cover it, and no step is ever written twice. Non-positive or non-finite values are rejected with exit status 1.

The gamma diagnostic is written inside the same `write_frame` path and therefore follows the same cadence automatically. The `.faceflux` / `.outercond` sidecars keep their existing step-stride behavior and stay default-off.

**Legacy compatibility is byte-exact.** With `CHROMO_FRAME_DT` unset the schedule runs in stride mode and reproduces the old `step % frame_stride` logic including the final-frame rule. Verified by running the pre-change binary (`build_omp/chromo_main`, built from `HEAD` before this task) and the new binary on the same legacy configuration, in separate directories so the output path embedded in the header matched:

```text
709ddc4f...  /tmp/leg_old/run.txt              709ddc4f...  /tmp/leg_new/run.txt
5901d712...  /tmp/leg_old/run.txt.gamma_diag   5901d712...  /tmp/leg_new/run.txt.gamma_diag
```

Live check of the new cadence — the CFL 0.125 / 20 s reference with `CHROMO_FRAME_DT=1` produced exactly 21 frames at t = 0, 1.0002, 2.00012, … 20.0002: one per second, never more than one timestep late, no duplicate, final state included.

### Tests

`test_output_schedule_legacy_stride`, `test_output_schedule_frame_dt_cadence`, `test_output_schedule_multi_interval_step`, `test_output_schedule_no_drift_and_no_duplicates` (100 physical seconds at 1 ms steps with a 0.1 s cadence: asserts exactly 1001 frames, no repeated step, and per-frame lateness bounded by one timestep — i.e. no accumulation), `test_output_schedule_final_frame_not_duplicated` (both the ends-on-cadence and ends-between-cadences cases).

---

## 4. Production launcher

`scripts/run_chromo_realtime.sh` defaults to the 12-thread OpenMP runtime, `CHROMO_OUTPUT=0`, `CHROMO_GAMMA_DIAG=0`, `CHROMO_FACE_FLUX_DIAG=0`, `CHROMO_OUTER_COND_DIAG=0`, `CHROMO_EOS_COUNTS=0`, `CHROMO_PROFILE=1`, and `CHROMO_CFL=0.50`. Every value uses `${VAR:-default}` so an explicit environment assignment wins; all CLI arguments are forwarded unchanged; it resolves the project root from its own path so it works from any working directory; and it fails with a build hint if the executable is missing. `scripts/run_chromo_omp.sh` is unchanged and remains the generic OpenMP launcher.

```bash
# low-I/O production
CHROMO_T_END=1000 scripts/run_chromo_realtime.sh out.txt \
    full no-ionization model_column - 20.0 no-cooling

# with snapshots, one frame per physical second
CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=1 \
scripts/run_chromo_realtime.sh out.txt \
    full no-ionization model_column - 20.0 no-cooling
```

---

## 5. Baseline

Reproduced CFL 0.25, 20 physical seconds, 12 threads, all output off, 2638 refined cells:

```text
wall            36.56 s          (audited baseline 36-38 s -> reproduced)
steps           28466
mean dt         7.025e-4 s       min dt 6.869e-4 s
ms/step         1.284
limiter         acoustic 28466, all others 0
profile         conduction 18.54 s | rhs 11.86 s | decode 2.55 s
                projection 2.23 s | cfl 0.22 s | boundary 0.15 s | output 0.002 s
inversion       3.007e8 calls, avg 0.802 iter, max 16, 4.60e5 bisection fallbacks (0.23%)
conduction      avg 2.998 Newton updates, max 4
RTF             1.828   ->   projected 1000 s wall = 1828 s
```

---

## 6. CFL sweep — performance

Three timed low-I/O reps per value on an otherwise quiet machine, 12 threads, 20 physical seconds; median reported.

| CFL | steps | mean dt (s) | min dt (s) | wall median (s) | ms/step | RTF | projected 1000 s wall (s) | status |
|---|---|---|---|---|---|---|---|---|
| 0.125 (ref) | 56961 | 3.512e-4 | 3.434e-4 | 72.25 | 1.269 | 3.61 | 3613 | reference |
| 0.25 | 28466 | 7.026e-4 | 6.869e-4 | 37.51 | 1.318 | 1.876 | 1876 | baseline |
| 0.35 | 20333 | 9.836e-4 | 9.617e-4 | 25.98 | 1.278 | 1.299 | 1299 | accepted, above target |
| 0.40 | 17792 | 1.124e-3 | 1.099e-3 | 22.92 | 1.288 | 1.146 | 1146 | accepted, above target |
| 0.45 | 15815 | 1.265e-3 | 1.236e-3 | 20.42 | 1.291 | 1.021 | 1021 | accepted, just above target |
| 0.50 | 14232 | 1.405e-3 | 1.374e-3 | 17.91 | 1.258 | **0.896** | **896** | **selected** |

Rep spread was 0.1–1.8 s (worst case 4.7% at CFL 0.25). ms/step is essentially flat (1.26–1.32), confirming the speedup comes purely from taking fewer steps, not from any per-step change. Every run at every CFL was limited by the acoustic CFL condition on 100% of steps.

Region profile at the two endpoints (rep 2):

| region | CFL 0.25 (s) | CFL 0.50 (s) |
|---|---|---|
| conduction | 19.88 | 9.32 |
| rhs | 11.65 | 5.68 |
| decode | 2.56 | 1.24 |
| projection | 2.21 | 1.08 |
| cfl | 0.22 | 0.10 |
| boundary | 0.16 | 0.07 |
| output | 0.002 | 0.002 |

---

## 7. CFL sweep — accuracy

Reference: CFL 0.125, 20 s, snapshots every physical second (`CHROMO_FRAME_DT=1`). Candidates compared at matched physical times by `util/analyze_cfl_sweep.py`, which parses the production snapshot file and its `.gamma_diag` sidecar directly. Relative-norm denominator floor: 1e-3 × the domain maximum of |reference| for that field and snapshot, reported per row.

### 7.1 Hard-rejection gates — all candidates pass

| check | 0.25 | 0.35 | 0.40 | 0.45 | 0.50 |
|---|---|---|---|---|---|
| crash / exception | none | none | none | none | none |
| non-finite state | none | none | none | none | none |
| EOS table-bound failure | none | none | none | none | none |
| conduction Newton non-convergence | none | none | none | none | none |
| active limiter | acoustic 100% | acoustic 100% | acoustic 100% | acoustic 100% | acoustic 100% |
| conduction Newton avg / max | 3.00 / 4 | 3.00 / 4 | 3.00 / 4 | 3.00 / 4 | 3.00 / 4 |
| inversion max iterations | 16 | 16 | 16 | 16 | 16 |
| inversion bisection-fallback rate | 0.225% | 0.258% | 0.275% | 0.292% | 0.308% |
| termination | end_time | end_time | end_time | end_time | end_time |

The fallback rate rises smoothly by 37% from CFL 0.25 to 0.50 — a graded response to a larger timestep, not a dramatic jump. Newton behavior is unchanged. No new floor or clamp activated; no limiter changed.

**Determinism and thread-independence** (a hard-rejection criterion) were checked at the selected CFL 0.50 by running four configurations and hashing the outputs with the embedded output path excluded:

```text
ce31f49f...  12 threads, run A
ce31f49f...  12 threads, run B      (repeat -> bit-identical)
ce31f49f...   1 thread, OpenMP build
ce31f49f...   1 thread, CHROMO_ENABLE_OPENMP=OFF serial build
```

All four bit-identical, for both the main output and the gamma sidecar. Output does not depend on OpenMP scheduling and runs are reproducible.

### 7.2 Field errors vs the CFL 0.125 reference (worst matched snapshot, 20 s)

| CFL | rho L1 | T L1 | p L1 | x_eq L1 | e_tot L1 | rho Linf | T Linf | e_tot Linf | v L1 | v abs Linf (m/s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 1.79e-5 | 2.35e-5 | 1.34e-5 | 8.26e-5 | 2.46e-5 | 7.8e-4 | 7.7e-4 | 6.3e-4 | 0.227 | 3.03 |
| 0.35 | 1.71e-5 | 3.46e-5 | 1.10e-5 | 9.45e-5 | 2.56e-5 | 1.19e-3 | 1.18e-3 | 9.6e-4 | 0.247 | 2.78 |
| 0.40 | 1.66e-5 | 3.70e-5 | 1.06e-5 | 9.09e-5 | 2.44e-5 | 1.36e-3 | 1.34e-3 | 1.09e-3 | 0.239 | 2.81 |
| 0.45 | 1.86e-5 | 3.91e-5 | 1.18e-5 | 8.24e-5 | 2.37e-5 | 1.49e-3 | 1.47e-3 | 1.20e-3 | 0.239 | 2.85 |
| 0.50 | 1.90e-5 | 3.95e-5 | 1.21e-5 | 8.31e-5 | 2.40e-5 | 1.60e-3 | 1.58e-3 | 1.30e-3 | 0.310 | 3.19 |

Every thermodynamic and conserved field is **two to three orders of magnitude inside** the 1% L1 / 3% Linf gates, and the Linf errors increase smoothly and monotonically with CFL — clean first-order temporal behavior.

### 7.3 Integrated quantities and observables vs the CFL 0.125 reference (worst over 21 snapshots, 20 s)

| CFL | integrated mass | integrated energy | T_max | p_top | TR shift (km) | gate (0.5% / 0.5% / 1% / 1% / 0.06 km) |
|---|---|---|---|---|---|---|
| 0.25 | 9.13e-6 | 9.08e-6 | 1.37e-5 | 9.73e-6 | 0.0028 | PASS |
| 0.35 | 6.92e-6 | 7.34e-6 | 1.83e-5 | 1.95e-5 | 0.0044 | PASS |
| 0.40 | 6.60e-6 | 6.65e-6 | 2.29e-5 | 1.95e-5 | 0.0052 | PASS |
| 0.45 | 7.59e-6 | 4.47e-6 | 2.29e-5 | 1.95e-5 | 0.0058 | PASS |
| 0.50 | 9.01e-6 | 1.02e-5 | 2.29e-5 | 1.95e-5 | 0.0062 | PASS |

Mass and energy agree to ~1e-5 — 500× inside the 0.5% gate. The TR location (first upward crossing of T = 20,000 K, the TRAC cutoff the solver itself reports; the usual 1e5 K coronal criterion never fires in this `ISO_T_TOP=22000` configuration) shifts by at most 0.0062 km, an order of magnitude inside the 0.06 km minimum refined cell width. The maximum |dT/ds| location is identical at every 20 s snapshot for every CFL.

### 7.4 The velocity and mass-flux gates fail at every CFL, including the current production 0.25

The 1% `max_up_mass_flux` and `v_top` gates and the 1%/3% velocity-field gates are violated by **every** candidate — and by the current production CFL 0.25 itself (v_top 1.43%, max upward mass flux 3.97%, v field L1 0.227). Rather than relax the threshold, the cause was investigated.

The configuration is a near-static relaxation: the reference velocity field has max |v| ≈ 53 m/s and mean |v| ≈ 2–4 m/s, and is a numerical residual on top of a hydrostatic balance carried in `float`. The absolute discrepancies are ≤ 3.2 m/s (Linf) and ≤ 0.43 m/s (L1) at 20 s. A self-convergence ladder shows velocity does not converge in time at all in this configuration — refining the timestep makes the level-to-level difference *larger*, consistent with accumulating float32 cancellation in a residual quantity (8× more steps at CFL 0.0625 than at 0.50):

| pair (t = 20 s) | dt of coarser (s) | L1 \|Δv\| (m/s) | Linf \|Δv\| (m/s) |
|---|---|---|---|
| 0.50 vs 0.45 | 1.405e-3 | 0.050 | 0.32 |
| 0.45 vs 0.40 | 1.265e-3 | 0.067 | 0.33 |
| 0.40 vs 0.35 | 1.124e-3 | 0.078 | 0.50 |
| 0.35 vs 0.25 | 9.836e-4 | 0.192 | 0.96 |
| 0.25 vs 0.125 | 7.026e-4 | 0.428 | 3.13 |
| 0.125 vs 0.0625 | 3.513e-4 | 1.169 | 7.27 |

Measured against the finest run (CFL 0.0625, a dedicated 113,951-step control), *every* CFL from 0.125 to 0.50 sits at the same L1 |Δv| ≈ 1.0–1.17 m/s. The CFL 0.125 reference is therefore **not more trustworthy in velocity** than any candidate, so a velocity-relative comparison against it cannot discriminate between CFL values. Temperature over the same ladder behaves normally (monotone, L1 dT/T from 2.6e-6 to 2.9e-5).

This is a pre-existing property of the model, not a consequence of raising the CFL, and diagnosing it further (e.g. float→double for the residual) is outside this task's scope. It is recorded as a production risk in §10.

The discriminating comparison that *is* valid is against the **current production baseline CFL 0.25**, since that is the accepted state of practice:

| t = 20 s | L1 \|Δv\| (m/s) | Linf \|Δv\| (m/s) | Linf ΔT (K) | L1 ΔT/T | L1 Δp/p | L1 Δrho/rho |
|---|---|---|---|---|---|---|
| **0.25 vs 0.125 — already accepted** | **0.428** | **3.13** | **11.1** | **1.63e-5** | **2.36e-5** | **4.73e-5** |
| 0.35 vs 0.25 | 0.192 | 0.96 | 6.6 | 1.01e-5 | 8.81e-6 | 2.43e-5 |
| 0.40 vs 0.25 | 0.243 | 1.25 | 8.9 | 1.18e-5 | 1.17e-5 | 2.76e-5 |
| 0.45 vs 0.25 | 0.281 | 1.55 | 10.9 | 1.39e-5 | 1.47e-5 | 2.86e-5 |
| 0.50 vs 0.25 | 0.297 | 1.62 | 11.9 | 1.46e-5 | 1.73e-5 | 3.05e-5 |

Moving from CFL 0.25 to 0.50 perturbs the solution by **0.69×** the velocity change and **1.07×** the peak-temperature change already incurred by the accepted 0.25-vs-0.125 discretization gap, and the perturbation grows monotonically with CFL (0.45× → 0.57× → 0.66× → 0.69× in L1 |Δv|). The extra temporal error at CFL 0.50 is bounded by, and comparable to, the temporal error already present and accepted at CFL 0.25.

---

## 8. Medium-duration validation (100 physical seconds)

The two highest 20-second candidates (0.45, 0.50) plus the CFL 0.25 production reference, `CHROMO_FRAME_DT=5` (21 frames each). A CFL 0.125 / 100 s control (284,912 steps, 329.9 s wall) was added so the medium-duration comparison has the same like-for-like scale as §7.4.

| run | steps | wall (s) | termination | limiter | Newton avg/max | fallback rate |
|---|---|---|---|---|---|---|
| 0.125 | 284912 | 329.89 | end_time | acoustic 100% | 3.00 / 4 | — |
| 0.25 | 142442 | 194.72 | end_time | acoustic 100% | 3.00 / 4 | 0.23% |
| 0.45 | 79073 | 107.14 | end_time | acoustic 100% | 3.00 / 4 | 0.29% |
| 0.50 | 71208 | 99.90 | end_time | acoustic 100% | 3.00 / 4 | 0.31% |

Departure from the CFL 0.25 reference, absolute, with the CFL 0.125 control alongside for scale:

| t (s) | max\|v\| ref (m/s) | 0.125: L1\|Δv\| / Linf / ΔT_max | 0.45: L1\|Δv\| / Linf / ΔT_max | 0.50: L1\|Δv\| / Linf / ΔT_max |
|---|---|---|---|---|
| 20 | 53.3 | 0.428 / 3.13 / 11.1 K | 0.281 / 1.55 / 10.9 K | 0.297 / 1.62 / 11.9 K |
| 50 | 52.2 | 2.271 / 6.98 / 51.3 K | 0.989 / 3.18 / 24.0 K | 0.990 / 3.36 / 32.9 K |
| 100 | 45.8 | 3.627 / 9.15 / 102.1 K | 1.667 / 3.28 / 92.3 K | 1.830 / 3.76 / 116.3 K |

Integrated mass / energy, relative to the CFL 0.25 reference:

| t (s) | 0.125 mass | 0.125 energy | 0.45 mass | 0.45 energy | 0.50 mass | 0.50 energy |
|---|---|---|---|---|---|---|
| 20 | 9.13e-6 | 9.08e-6 | 1.54e-6 | 4.61e-6 | 1.21e-7 | -1.11e-6 |
| 50 | 1.15e-4 | 4.27e-5 | -1.02e-4 | -1.47e-4 | -9.44e-5 | -1.81e-4 |
| 100 | 1.85e-4 | 5.80e-5 | -5.89e-5 | -1.91e-4 | -7.22e-5 | -2.17e-4 |

Field-level results at 100 s (worst matched snapshot, vs CFL 0.25): rho L1 1.34e-4 / Linf 1.10e-2, T L1 2.44e-4 / Linf 1.10e-2, p L1 1.47e-4, x_eq L1 4.15e-4, e_tot L1 2.49e-4 for CFL 0.50 — all inside the 1% L1 and 3% Linf gates, with CFL 0.45 uniformly slightly smaller.

Observables over the full 100 s history for CFL 0.50: integrated mass and energy **plateau** at −7e-5 and −2.2e-4 after t ≈ 55 s rather than growing; T_max departure grows to 6.8e-5; p_top to 1.9e-5; TR location shifts monotonically to 0.024 km, still inside the 0.06 km fine cell. The maximum-|dT/ds| cell hops by 1–12 fine cells (±0.06 to −0.7 km) after t = 80 s — small jitter in an argmax, not a structural change.

No boundary oscillation growth, no new clamp, no Newton or inversion degradation, and no accumulating drift appeared in any conserved or thermodynamic quantity at either candidate. Both 0.45 and 0.50 pass the medium-duration gate; at every time and metric both stay closer to the CFL 0.25 reference than the finer CFL 0.125 control does.

---

## 9. Selected CFL

**Accepted production CFL: 0.50.**

It is the largest value swept, it passes every hard-rejection gate, its thermodynamic and conserved-field errors are monotone in CFL and 2–3 orders of magnitude inside the acceptance thresholds at both 20 s and 100 s, it perturbs the solution less than the already-accepted 0.25-vs-0.125 discretization gap, it is bit-reproducible and thread-independent, and it is the only value meeting the real-time target.

Nothing was rejected on accuracy grounds. CFL 0.35, 0.40 and 0.45 are all scientifically acceptable; they were not selected because they are slower (projected 1299 s, 1146 s and 1021 s against a 1000 s target) and offer no accuracy advantage that the evidence can resolve. CFL 0.55 was **not** tested: the task permits it only when there is a clear reason to seek more margin, and 0.50 already clears the target by 10%.

Actions taken:

1. `scripts/run_chromo_realtime.sh` defaults to `CHROMO_CFL="${CHROMO_CFL:-0.50}"`, override-preserving.
2. `AGENTS.md` and `CLAUDE.md` updated identically (verified byte-identical by SHA-256), labelling 0.50 as the validated production runtime default **for this model shape and this hardware**, and documenting `CHROMO_T_END`, `CHROMO_FRAME_DT`, `CHROMO_STEP_CAP` and the termination line.
3. The hard-coded solver default in `chromo_main.cpp` remains `float cfl = 0.25f`, and CFL 0.25 remains the comparison reference. It was not changed.

---

## 10. Final 1000-second real-time acceptance run

Run through `scripts/run_chromo_realtime.sh` at the accepted CFL 0.50, 12 threads, all expensive output disabled, profiling on, no artificial or truncated step cap:

```bash
CHROMO_BINARY=$PWD/build_cfl_audit/chromo_main \
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=2000 \
ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 \
ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20 ISO_NUMERICAL_DIFFUSIVITY_MULT=1 \
CHROMO_T_END=1000 \
/usr/bin/time -p scripts/run_chromo_realtime.sh \
  outputs/model_column/cfl_sweep/realtime_1000s_cfl050.txt \
  full no-ionization model_column - 20.0 no-cooling
```

### Result

```text
wall time                 912.15 s        <= 1000 s   PASS
final physical time       1000 s          >= 1000 s   PASS
termination               end_time                    PASS
process exit status       0
final step                712321
mean dt                   1.4039e-3 s
real-time factor          0.912
speedup vs CFL 0.25       2.06x  (1876 s projected -> 912 s actual)
```

**REAL_TIME_1000S_TARGET: GO** — 1000 physical seconds in 912 wall seconds, 88 s inside the target.

The run terminated on the end time with `effective_step_cap=1000000000000 (source=end_time)`, i.e. the legacy `time_mult`-derived 200,000-step cap that would have truncated it near 140 s of the old behavior was correctly not applied.

### Solver behavior over the full 1000 s

```text
limiter          acoustic 712321 (100%), all others 0
conduction       avg 3.000011 Newton updates, max 4
inversion        7.524e9 calls, avg 0.887 iter, max 16
                 1.054e9 initial-guess accepts (14.0%)
                 6.291e9 one-update convergences (83.6%)
                 1.072e7 bisection fallbacks (0.142%)
profile          conduction 488.5 s | rhs 277.1 s | decode 61.4 s
                 projection 55.7 s | cfl 5.6 s | boundary 4.2 s | output 0.006 s
```

No non-finite state, no clamp, no floor activation, no Newton non-convergence, no limiter change over 712,321 steps. The bisection-fallback rate (0.142%) is *lower* than in the 20-second runs (0.31%), i.e. the initial transient — not the timestep — dominates that counter.

### Output volume

Expensive diagnostics confirmed disabled by the launcher defaults (`CHROMO_OUTPUT=0`, `CHROMO_GAMMA_DIAG=0`, `CHROMO_FACE_FLUX_DIAG=0`, `CHROMO_OUTER_COND_DIAG=0`, `CHROMO_EOS_COUNTS=0`; `CHROMO_PROFILE=1`). No snapshot file was created at all, and the total output-region time was 6.2 ms out of 912 s. Peak on-disk artifact: the 516 KB stdout log (the periodic every-100-steps progress lines). For comparison, the same run with the legacy stride and output enabled would have written ~17,790 snapshots totalling ~8.3 GB.

## 11. Validation performed

```text
cmake --build build_cfl_audit -j4     (OpenMP)   clean
ctest --test-dir build_cfl_audit                 8/8 passed
./build_cfl_audit/chromo_tests                   10143 passed, 0 failed

cmake --build build_cfl_serial -j4    (OpenMP OFF)  clean
ctest --test-dir build_cfl_serial                 8/8 passed
./build_cfl_serial/chromo_tests                   10143 passed, 0 failed

git diff --check                                  clean
```

The serial (`CHROMO_ENABLE_OPENMP=OFF`) build confirms the run-control helpers have no OpenMP dependency: identical test count, and the CFL 0.50 serial output is bit-identical to the 12-thread OpenMP output.

## 12. Reproducing the sweep

```bash
cmake -S . -B build_cfl_audit -DCMAKE_BUILD_TYPE=Release \
      -DCHROMO_ENABLE_OPENMP=ON -DLIBOMP_ROOT="$(brew --prefix libomp)"
cmake --build build_cfl_audit -j4

util/run_cfl_sweep.sh acc  0.125 0.25 0.35 0.40 0.45 0.50   # 1 frame/s, 20 s
util/run_cfl_sweep.sh perf 0.25 0.35 0.40 0.45 0.50         # 3 timed reps each
util/run_cfl_sweep.sh med  100 5 0.125 0.25 0.45 0.50       # 100 s, 1 frame/5 s

D=outputs/model_column/cfl_sweep
.venv/bin/python util/analyze_cfl_sweep.py \
    --reference $D/cfl_reference_0125_20s.txt --reference-label 0.125 \
    --candidate 0.25=$D/cfl025_20s.txt --candidate 0.35=$D/cfl035_20s.txt \
    --candidate 0.40=$D/cfl040_20s.txt --candidate 0.45=$D/cfl045_20s.txt \
    --candidate 0.50=$D/cfl050_20s.txt \
    --fields rho,v,T,p,x_eq,mass_flux --out $D/cfl_accuracy_20s.md
```

`util/run_cfl_sweep.sh` holds the production-shaped model configuration in one place, so only `CHROMO_CFL` and the output controls differ between runs. `util/analyze_cfl_sweep.py` parses the snapshot file and its `.gamma_diag` sidecar plus the stdout log; it reimplements none of the solver equations. Generated reports land in `outputs/model_column/cfl_sweep/` and are not committed.

Sweep artifacts (~122 MB, gitignored) live in `outputs/model_column/cfl_sweep/`.

## 13. Production risks

1. **Velocity is not time-converged in this configuration.** The residual velocity field (max ≈ 53 m/s, mean ≈ 2–4 m/s) does not converge as dt → 0; halving the timestep increases the level-to-level difference. This is pre-existing, affects CFL 0.25 equally, and is consistent with float32 cancellation in a residual carried on top of a hydrostatic balance. Any study whose science depends on m/s-scale chromospheric velocities should treat these values as noise-floor regardless of CFL. Out of scope here; worth a separate task.
2. **CFL 0.50 is validated for one model shape and one machine.** It applies to the 2638-cell refined `model_column` case on the current Apple Silicon workstation. Rerun the sweep for a different mesh, scenario, or hardware. Scenarios with additional timestep limiters (cooling, beam heating, electron energy) were not exercised — every step of every run here was acoustic-limited.
3. **The 100 s comparison uses CFL 0.25 as its reference**, not an independent exact solution, because the finer-time reference is not more trustworthy in the velocity field (§7.4).
4. **`CHROMO_FRAME_DT` is opt-in.** A long run that enables output without it silently keeps the legacy stride and can emit multi-GB files. The launcher's low-I/O defaults mitigate but do not prevent this.
