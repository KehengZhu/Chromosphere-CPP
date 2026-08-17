# Double-precision control experiment: the lower-chromosphere mass-flux gradient is float32 storage round-off

**Status: current.** Diagnostic experiment, 2026-08-16/17. **It changes no release default.** The release conserved state is and remains `float`. The double-precision build exists only to isolate one question and is not a supported configuration.

## 1. The question

In a quasi-steady state on a constant-area flux tube, continuity forces the conservative face mass flux `f_total` to be independent of height. The 4000 s release run does not show that: `f_total` at the 1600 km base is 2.8x its value at the top, and the excess is concentrated in the lower chromosphere. Two candidate explanations had to be separated:

1. a physical or boundary-condition effect — an unrelaxed transient, or the inner reservoir closure driving a spurious base flux; or
2. float32 representation round-off in the conserved state, which `docs/roe_n1000_lower_chromosphere_sloshing_diagnosis.md` had already flagged indirectly through an unexplained inconsistency between the `.faceflux` divergence and the density evolution recorded in the same file.

The measured mechanism for (2) was identified first, from the release run alone: `mixture_advance` performs `next = state + dt*rhs` on a float32 state, and in the lower chromosphere `dt*(-df/ds)` is about `2.6e-17 kg m^-3 s^-1` while one float32 ULP of `rho ~ 9e-10` is about `5.6e-17`. **619 of 660 cells sit below half a ULP** (median ratio 0.30), so the addition rounds straight back to `rho`. There is no compensated summation, so this is permanent stagnation rather than slow drift. That is a mechanism, not proof; the control experiment supplies the proof.

## 2. Design

One variable changed. `chromosphere::Vec` was made `arma::Col<double>` under a new default-`OFF` CMake option, `CHROMO_STATE_FLOAT64`, and the two legs were run from the same source tree:

```bash
cmake -S . -B build_omp_f64 -DCMAKE_BUILD_TYPE=Release \
      -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_STATE_FLOAT64=ON
```

Physics, grid, boundary conditions, reconstruction, the SWMF Godunov flux, the CFL rule, the EOS table and every diagnostic are byte-for-byte the same code. The typedef also promotes the static mesh metric arrays, which share the type; see the caveat in §5.

That the two legs are the same *scheme* is confirmed by the run tallies, which are identical: **712,314 steps**, final physical time 4000 s, `termination=end_time`, **941,679,108 face Riemann problems, all solved exactly, 0 Rusanov fallbacks** in both. Wall time 502 s (float32) versus 507 s (float64) on 12 threads.

## 3. Results

| Quantity, N=500 / R4, 4000 s | float32 (release) | float64 (control) |
| --- | --- | --- |
| `f_total` base / top at 2000 s | 2.682 | **1.008** |
| `f_total` base / top at 4000 s | 2.811 | **0.996** |
| `f_total` spread over 1600–1850 km, 4000 s | 119.5 % | **2.0 %** |
| Cells with per-step mass increment below 0.5 ULP(`rho`) | 619 / 660 | **0 / 660** |
| Base cell velocity, 4000 s | 1.392 m/s | 0.506 m/s |
| Top cell velocity, 4000 s | 13.18 m/s | 13.58 m/s |
| `f_total` at the top face, 4000 s | 3.753e-10 | 3.869e-10 |
| `\|f_diff\|` / `f_total` at the first interior face | 16.3 % | 15.9 % |
| Cell-0 `rho V` / `f_total[0]` | 1.280 | 1.275 |

**The gradient is float32 round-off.** In double precision the column reaches a genuinely constant mass flux from base to top — flat to 1.1 % already at 2000 s and 2.0 % at 4000 s — which is what continuity demands and what the float32 run never achieves at any time in the run.

**The column mass budget closes in double precision and is off by two orders of magnitude in float32.** Integrating the *recorded* `f_total` divergence in time over all 217 sidecar frames and comparing with the *recorded* density change over the same run:

| Mean fractional density change, 1605–1850 km | predicted by the recorded flux divergence | observed | pred / obs |
| --- | --- | --- | --- |
| float32 | `+1.29e-2` | `-5.98e-5` | **-215** |
| float64 | `+3.95e-5` | `+4.20e-5` | **+0.94** |

In float32 the recorded fluxes imply the lower chromosphere should be *filling* by about a percent, while the stored density in fact *fell* by `6e-5` — wrong by a factor of order 200 **and wrong in sign**. In float64 the budget closes to 6 %, which is at or below the accuracy a 217-frame trapezoid over an evolving state can resolve. **This identifies the previously unexplained `.faceflux` divergence anomaly**: the captured flux values were always correct; it was the float32 density *response* that could not exist. It is a storage-precision artifact, not a capture-stage bug.

**The first-interior-face artifact is precision independent and therefore real.** `|f_diff|/f_total` at face 0 is 16.3 % versus 15.9 %, and the cell-0 `rho V` overshoot is 1.280 versus 1.275. The control does not touch it, so that localized feature is a genuine scheme/inner-ghost-closure effect and must be diagnosed separately from the broad gradient.

**Top-of-domain observables move by about 3 %.** Top velocity 13.18 → 13.58 m/s and top `f_total` 3.753e-10 → 3.869e-10. The release evaporation numbers are therefore not badly wrong: the error is concentrated where `rho` is large and its ULP coarse, which is exactly where the release already declines to quote velocities (see the m/s noise-floor limitation in @ref validation).

**Incidental finding.** `inversion.initial_guess_accepts` is 754,709,637 of 1,420,486,978 EOS energy inversions in float32 but only 48,151 in float64, and the mean iteration count rises from 0.50 to 1.03. Over half of all float32 inversions "converge" on the initial guess only because float32 cannot represent a better answer. The pressure inversion is unaffected (mean 3.14 evaluations in both). This costs nothing in accuracy that matters here — the inversion tolerance is far below the discretization error — but it means float32 inversion-iteration statistics are not a measure of the inverter's quality.

## 4. Conclusion

The broad 1600–1850 km `f_total` gradient in the release run is **entirely an artifact of float32 conserved-state storage**, produced by the per-step mass-density increment falling below half a ULP of `rho` through most of the lower chromosphere. It is not a transient, not the inner boundary condition, and not the Godunov flux. Consequences for practice:

- **Do not interpret a lower-chromosphere `f_total` height gradient as a steady-state statement**, and do not chase it as a boundary-condition bug.
- **Any quantitative chromospheric mass-budget or long-time mass-loading claim requires a double-precision conserved state or a perturbation-variable formulation first.** This is the concrete, measured cost of the float32 packing already flagged as a limitation.
- The localized first-face artifact remains open and is a separate, real defect.
- Release evaporation observables at the top of the domain are affected at the ~3 % level.

## 5. Caveats

**The typedef promotes the mesh metric arrays too**, because they share `chromosphere::Vec`. Those are geometry computed once, with a float32 relative error of order `1e-7`; they cannot produce a factor-2.8 flux gradient, and the ULP census points directly at the state update. But this single experiment does not formally separate the two contributions. A follow-up that promotes only the state rows would.

**`CHROMO_STATE_FLOAT64` is not a supported configuration.** It is not exercised by `scripts/release_validation.sh`, and no release number may be quoted from it. Under the option the test binary reports 11602 passed / 29 failed out of 11631 checks. Those failures were inspected and none were loosened: 28 are `expected exactly 0` assertions on residuals that are zero in float32 *only because float32 quantizes them away* (they are `1e-18`–`1e-13` in double), and one is a `Grid::dinvB_ds_i` geometry tolerance reflecting the metric promotion. One case, `test_release_production_table_domain_corners`, is skipped under the option: it walks the packed energy by single ULPs of the storage type to build a state that rounds outside the EOS table, so both its premise and its walk are float32 specific. **The float32 build passes 11643 / 11643 and is unchanged by this work.**

## 6. Reproduce

```bash
cmake -S . -B build_omp_f64 -DCMAKE_BUILD_TYPE=Release \
      -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_STATE_FLOAT64=ON
cmake --build build_omp_f64 -j 12

CHROMO_BINARY=$PWD/build_omp_f64/chromo_main \
CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_STRIDE=3300 \
scripts/run_chromo_realtime.sh \
    outputs/model_column/lnp_godunov_N500_4000s_f64.txt \
    full no-ionization model_column - 20.0 no-cooling

.venv/bin/python util/plot_precision_control.py \
    outputs/model_column/lnp_godunov_N500_4000s.txt \
    outputs/model_column/lnp_godunov_N500_4000s_f64.txt
```

The six-panel comparison figure and the animation of the double-precision run are catalogued in `visualize_commands.md` §14. The float32 mechanism panel alone is also produced by `util/plot_base_mass_flux_diagnosis.py`.
