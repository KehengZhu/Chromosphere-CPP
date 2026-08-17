# Resolution convergence scan: does the top-region `rho V` ripple vanish as `Delta h -> 0`?

> **Partly superseded (see `docs/studies/numerics/reference_free_release_recap.md`).** The release hydrodynamic operator has since gained the MUSCL-Hancock predictor momentum source and a trapezoidal EOS-closed lower ghost ladder, and **equilibrium-reference well balancing (`eq_wb` / `ISO_EQ_WB`) has been retired from the release** — it survives only in the two-fluid research solver. Statements below that describe `eq_wb` as active release method, or that quantify the discrete hydrostatic defect, are historical. Everything else in this document still stands.


Follow-up to `docs/studies/numerics/top_ripple_face_flux_diagnosis.md` (which established Case A: the ripple lives in cell-centred `rho V`, not in the conserved face flux) and `docs/studies/boundaries/upper_bc_hydro_temperature_decoupling_recap.md` (which established that the fixed hydro ghost temperature explains the *last cell* only). This stage ran the one experiment those two recaps identified as the cheapest decisive next test, and nothing else. No production numerical scheme was modified.

## Headline

**Two results hold simultaneously, and they must not be merged.**

- **Result A — the ripple converges.** The cell-centred ripple amplitude falls monotonically with `Delta h` at every metric, the reference-subtracted diffusive contribution falls with it, `f_eff` stays smooth at every resolution, and the physical/numerical mismatch `A_mismatch` goes to zero at the same rate as the ripple itself. The answer to the stage's boxed core question — *does `rho V - F_eff` tend to zero as `Delta h -> 0`* — is **yes**, at approximately first order.

- **Result D — the mean evaporation mass flux does not converge.** `mean(F_eff)` in the window falls roughly *linearly* in `Delta h` (3.69 -> 2.40 -> 1.84 x10^-9 kg m^-2 s^-1 at 1000 s) with no sign of flattening, and the realized conductive boundary input falls with it (`E_cond(1000 s)` = 5513 / 3329 / 2120 J m^-2). **No resolution in this scan produces a grid-independent evaporation rate.** The ripple being a truncation artefact does *not* rescue the evaporation number.

So `N = 1000` is adequate for statements about *ripple structure* being numerical, and is **not** adequate for any quantitative statement about the evaporation rate, the top-cell velocity, or the conductive input.

## 1. Configuration parity

All three runs share the identical command apart from `ISO_NS` and the two step-based diagnostic cadences. From the `.faceflux` headers:

| | ns=500 | ns=1000 | ns=2000 |
|---|---:|---:|---:|
| `Delta h` [km] | 1.1060 | 0.5530 | 0.2765 |
| `uniform_mesh` | 1 | 1 | 1 |
| `mc3` / `beta` | 1 / 2 | 1 / 2 | 1 / 2 |
| `eq_wb` | 1 | 1 | 1 |
| faceflux records | 169 | 174 | 177 |
| cells in the 2130-2150 km window | 18 | 37 | 72 |

Held fixed by construction: physical domain 1600-2153 km (`ISO_H_BASE=1600 ISO_DH=553`); the C7/pchip IC with hydrostatically re-integrated density; `data/eos/gamma1_hydrogen_v1.dat`; `ISO_HYDRO_T_DECOUPLE=1` (the decoupled hydro-T upper BC, i.e. the primary configuration this stage was told to use); fixed reservoir back-pressure semantics (`kOuterPRef` captured from the IC in the same way); conduction wall 22 000 K (`ISO_T_TOP`); cooling off (`no-cooling`, `ISO_COOLING` unset); TRAC off; MC3 `beta=2` (`ISO_LIMITER` unset = production); CFL 0.25; Mach cap `ISO_VCAP=0.1`; single-fluid, frozen-ionization gamma-table path; `eq_wb` reference built the same way (lazily, on the first step after `model_column_update_bc`, from a fresh IC in every run).

**Numerical diffusivity was deliberately left to scale.** `scenarios/model_column.cpp:544` sets `numerical_diffusivity = 2.0e3 * ds_m`, i.e. proportional to `Delta h`. It was *not* pinned to any one resolution's value, per the stage instruction.

**Diagnostic cadence matched in physical, not step, units.** `CHROMO_FACE_FLUX_STRIDE` = 500 / 1000 / 2000 steps and `CHROMO_FRAME_STRIDE` = 2500 / 5000 / 10000, so all three capture at ~5.8 s physical cadence. `CHROMO_FACE_FLUX_TOP` = 75 / 150 / 300 cells = the same 83 km physical extent, so the comparison window is strictly interior to the capture at every resolution.

**The existing `ns=1000` run was NOT reused.** Its `.faceflux` held only two records (`t=0` and `t=1000`), so it could not supply a `t=500` face record, and its capture cadence was set in steps against a different frame stride. Per the stage rule ("若存在任何配置差异，则重新运行 N=1000"), `ns=1000` was re-run from scratch inside this matched set.

### One unavoidable, first-order configuration difference — recorded, not hidden

The IC is sampled at cell centres, and the top cell centre sits at a different height at each resolution (2152.45 / 2152.72 / 2152.86 km). The C7 pchip `T(h)` is steep there, so the *initial* top-cell temperature differs:

| | ns=500 | ns=1000 | ns=2000 |
|---|---:|---:|---:|
| IC `T_top` [K] | 21 573.8 | 22 275.3 | 22 667.5 |
| IC `p_top` [Pa] | 0.0102909 | 0.0102888 | 0.0102877 |
| captured `kOuterPRef` [Pa] | 0.0102822 | 0.0102845 | 0.0102856 |

The back-pressure anchor agrees to **0.03 %** across the scan, so the pressure boundary really is the same. The top-cell IC *temperature* however differs by ~1100 K between ns=500 and ns=2000 — a genuine O(`Delta h`) property of point-sampling a steep profile, not a configuration error, but it does mean the initial wall-to-top-cell temperature gap (and hence the initial conductive drive) is itself resolution-dependent. **This is a plausible contributor to Result D and is flagged as not-yet-separated.**

`time_mult` (argv[6]) was raised 30.0 -> 60.0 for all three runs because `ns=2000` needs more than the `10000*time_mult` step cap (it took 350 540 steps; the first attempt truncated at 300 000 steps / t=855 s and was discarded). `CHROMO_T_END` overrides `total_time`, so `time_mult` touches only `step_cap`; verified by re-running ns=500 and ns=1000 at 60.0 and comparing to the 30.0 runs — the `.gamma_diag` and `.faceflux` sidecars are **byte-identical** (1 849 237 and 4 345 986 bytes for ns=500; 3 801 628 and 8 960 788 for ns=1000).

## 2. Diagnostic provenance

The production face-flux capture was reused unchanged: `Grid::capture_face_flux` / `GammaFaceFluxCapture` (`chromosphere.hpp:123`, `double` arrays throughout), filled by `capture_gamma_face_flux` in `src/rhs.cpp`, emitted by `chromo_main` under `CHROMO_FACE_FLUX_DIAG=1`. **No production C++ was touched this stage.** The only new code is the Python analysis script `util/resolution_scan_diag.py`, which imports its readers and the `roughness`/`relstd` definitions from the existing `util/face_flux_diag.py` and `util/upper_bc_decouple_diag.py` so the metric definitions are literally shared, not re-typed.

Per-record sanity checks, printed for every run and time:

- `f_central + f_diff == f_total` to `<= 1.9e-6` relative (float print rounding);
- the `t=0` capture really is the frozen `eq_wb` reference: its divergence reproduces the captured `eq_residual_mass` row to `<= 1.7e-6` relative, at all three resolutions.

## 3. Metric definitions

`M_i = rho_i V_i`; `F_eff = F_total - F_ref` with `F_ref = F_total` at `t=0`; `delta F_diff = F_diff(U) - F_diff(U_ref)`.
`A_x = std(x)/|mean(x)|` for `x = M, F_eff`; `A_diff = std(delta F_diff)/|mean(F_eff)|`; `Q_x = mean|D2 x| / mean|x|`; `R_eff[i] = -(F_eff[i]-F_eff[i-1])/Delta h`.

**Face-to-cell mapping for the mismatch metric.** `F_eff` lives on faces, `M` on cell centres. The face flux is mapped *down* to cell centres with the symmetric two-point average `F_eff_cell[i] = 0.5(F_eff[i-1] + F_eff[i])` (the cell's own lower and upper face). This is second-order accurate for a smooth field, so it cannot manufacture a first-order signal; and the analysis window is strictly interior to the capture, so `i-1` always exists. Then

`A_mismatch = RMS[ M_i - F_eff_cell[i] ] / |mean(F_eff)|`.

All metrics use the **fixed physical window 2130-2150 km**, never a fixed cell count.

## 4. Core results

### t = 1000 s (records at t = 999.996 / 999.998 / 999.999 s)

| metric | ns=500 | ns=1000 | ns=2000 | coarse/mid | mid/fine | p (3-pt fit) |
|---|---:|---:|---:|---:|---:|---:|
| `mean(rho V)` | 3.6438e-09 | 2.4088e-09 | 1.8425e-09 | 1.513 | 1.307 | 0.49 |
| **`A_M`** | **0.15422** | **0.12626** | **0.06648** | 1.221 | **1.899** | **+0.61** |
| **`Q_M`** | **0.15962** | **0.06821** | **0.02140** | 2.340 | **3.188** | **+1.45** |
| `mean(F_eff)` | 3.6922e-09 | 2.3999e-09 | 1.8417e-09 | 1.539 | 1.303 | 0.50 |
| `A_eff` | 0.00935 | 0.01141 | 0.00341 | 0.820 | 3.347 | +0.73 |
| `Q_eff` | 1.879e-04 | 1.440e-04 | 1.405e-04 | 1.304 | 1.025 | +0.21 |
| **`A_diff`** | **0.17882** | **0.12825** | **0.06863** | 1.394 | **1.869** | **+0.69** |
| `A_central` | 0.17859 | 0.12841 | 0.06875 | 1.391 | 1.868 | +0.69 |
| `corr(dcen, ddiff)` | -0.99863 | -0.99605 | -0.99877 | | | |
| median `a/|V|` | 289.2 | 441.9 | 581.6 | | | |
| **`A_mismatch`** | **0.15449** | **0.12592** | **0.06645** | 1.227 | **1.895** | **+0.61** |
| mismatch rms | 5.704e-10 | 3.022e-10 | 1.224e-10 | 1.888 | 2.469 | |
| mismatch absmax | 1.441e-09 | 1.016e-09 | 4.941e-10 | | | |
| `R_eff` rms | 6.026e-15 | 5.107e-15 | 1.446e-15 | 1.180 | 3.533 | +1.03 |
| `R_eff` rms / (`|F|/Delta h`) | 1.81e-03 | 1.18e-03 | 2.17e-04 | 1.534 | 5.422 | +1.53 |
| `R_total` rms / (`|F|/Delta h`) | 2.067 | 0.498 | 0.172 | | | |
| `mean(F_diff)` | 3.279e-10 | 8.028e-11 | 5.696e-12 | 4.08 | 14.1 | |
| `A_total` | 1.279 | 0.478 | 0.143 | | | |

### t = 500 s (records at t = 497.53 / 497.48 / 499.44 s)

| metric | ns=500 | ns=1000 | ns=2000 | coarse/mid | mid/fine | p (3-pt fit) |
|---|---:|---:|---:|---:|---:|---:|
| `mean(rho V)` | 4.2995e-09 | 2.7069e-09 | 1.8230e-09 | 1.588 | 1.485 | 0.62 |
| **`A_M`** | **0.15001** | **0.13652** | **0.06728** | 1.099 | **2.029** | **+0.58** |
| **`Q_M`** | **0.16776** | **0.08447** | **0.02171** | 1.986 | **3.890** | **+1.47** |
| `mean(F_eff)` | 4.3458e-09 | 2.6978e-09 | 1.8222e-09 | 1.611 | 1.481 | 0.63 |
| `A_eff` | 0.01114 | 0.01047 | 0.00868 | 1.064 | 1.206 | +0.18 |
| `Q_eff` | 1.722e-04 | 1.321e-04 | 1.460e-04 | 1.304 | 0.905 | +0.12 |
| **`A_diff`** | **0.16387** | **0.13884** | **0.06911** | 1.180 | **2.009** | **+0.62** |
| `corr(dcen, ddiff)` | -0.99769 | -0.99716 | -0.99216 | | | |
| median `a/|V|` | 244.9 | 396.9 | 583.6 | | | |
| **`A_mismatch`** | **0.15061** | **0.13636** | **0.06690** | 1.105 | **2.038** | **+0.59** |
| `R_eff` rms / (`|F|/Delta h`) | 2.18e-03 | 1.10e-03 | 6.89e-04 | 1.980 | 1.597 | +0.83 |
| `R_total` rms / (`|F|/Delta h`) | 1.756 | 0.443 | 0.174 | | | |

Full tables, per-metric fits, extremum lists and the Richardson block: `outputs/model_column/resconv_report.txt`, `outputs/model_column/resconv_metrics.json`.

### Reading

1. **Every ripple metric falls monotonically.** `A_M` 0.154 -> 0.126 -> 0.066; `Q_M` 0.160 -> 0.068 -> 0.021; `A_diff` 0.179 -> 0.128 -> 0.069; `A_mismatch` 0.154 -> 0.126 -> 0.066. Same direction and nearly the same factors at 500 s.
2. **The clean rate is the mid->fine pair, and it is first order.** `A_M` and `A_mismatch` both halve (1.90 / 1.90 at 1000 s; 2.03 / 2.04 at 500 s) when `Delta h` halves, i.e. local `p ~ 0.93-1.02`. The coarse->mid step is much weaker (1.22 / 1.23) — see the caveat below — which is why the 3-point fit reports the lower `p ~ 0.6`. **The honest statement is `A_M ~ O(Delta h)` between ns=1000 and ns=2000, with the ns=500 point unreliable.**
3. **`Q_M` is faster than first order** (ratios 2.34 then 3.19 at 1000 s, 1.99 then 3.89 at 500 s; fit `p ~ 1.45`). Curvature is the metric most sensitive to grid-scale structure, and it is the one falling fastest — consistent with a grid-scale numerical mode being damped away, not with a physical feature being resolved.
4. **`A_mismatch` tracks `A_M` to within 2 % at every resolution** (0.15449 vs 0.15422, 0.12592 vs 0.12626, 0.06645 vs 0.06648). That is the cleanest single statement of the whole stage: `F_eff_cell` is so smooth that the mismatch *is* the ripple, and it is going to zero.
5. **`F_eff` is smooth at every resolution.** `Q_eff` sits at 1.3-1.9e-4 throughout, i.e. `Q_M/Q_eff` = 850 (ns=500), 474 (ns=1000), 152 (ns=2000) — the gap narrows only because `Q_M` is collapsing toward `Q_eff`, exactly what convergence to the conserved flux looks like. `A_eff` stays at the 0.3-1.1 % level, and it is a monotone rise across the window (the slow mass loading), not an oscillation.
6. **`R_eff` normalised falls at ~first order** (2.18e-03 -> 1.10e-03 -> 6.89e-04 at 500 s). Discrete mass conservation of the effective flux improves with refinement; `R_total` normalised (1.76 -> 0.44 -> 0.17) remains the quantity to ignore, being dominated by the reference divergence `eq_wb` removes.
7. **The compensation mechanism is unchanged, only smaller.** `corr(delta F_central, delta F_diff)` stays at -0.992 to -0.999 at all three resolutions, and `A_central ~ A_diff` to 4 digits. So refinement does not change *what* the ripple is; it shrinks the amplitude of both halves of the anti-correlated pair together. Meanwhile `mean(F_diff)` collapses by 4.1x then 14.1x, i.e. the *net* Rusanov dissipative mass flux is vanishing fast even though `a/|V|` grows from 289 to 582 (the spectral radius ratio is not the whole story; the reconstructed `rho_R - rho_L` shrinks faster).

### The ns=500 caveat

At `Delta h = 1.106 km` the 2130-2150 km window holds only **18 cells** and shows **10 extrema** — a mean extremum spacing of 1.7 cells, i.e. essentially at the 2-cell Nyquist limit. `relstd` over 18 samples of a barely-resolved pattern is a noisy estimator, and a grid-limited pattern cannot get *worse* than Nyquist, which is why `A_M` is nearly flat from ns=500 to ns=1000 while `Q_M` (a local operator, less affected) still halves. **The ns=500 point is retained for completeness but the convergence claim rests on the ns=1000 -> ns=2000 pair.** This also means the reported 3-point `p` values are biased low and should not be quoted as the order of the scheme.

## 5. Ripple wavelength and peak positions

The answer is **both**, and the distinction matters.

**The dominant extrema are pinned to fixed physical heights.** Automated matching against the finest mesh (tolerance `1.5 Delta h` of the coarser run):

| | extrema matched to an ns=2000 counterpart | max offset | median offset | extra extrema on the fine mesh |
|---|---|---:|---:|---:|
| ns=500, t=1000 s | **10 / 10** | 0.70 km | 0.28 km | +7 |
| ns=1000, t=1000 s | **11 / 11** | 0.43 km | 0.15 km | +6 |
| ns=500, t=500 s | **10 / 10** | 0.70 km | 0.15 km | +9 |
| ns=1000, t=500 s | **11 / 11** | 0.43 km | 0.40 km | +8 |

Every extremum present on a coarse mesh survives at a fixed height on the fine mesh (~2131.7, 2133.9, 2135.6, 2138.4, 2141.1, 2141.7, 2143.9, 2145.0, 2145.5, 2147.8, 2149.4 km at ns=1000, all matched within 0.43 km). So the *envelope* is set by the physical TR structure, not by the mesh — consistent with the earlier finding that MC3 and minmod share peak positions.

**But refinement adds new, finer extrema between them, so the mean spacing is fixed in cells, not in km.** `lambda_extrema` (2 x mean extremum spacing):

| | ns=500 | ns=1000 | ns=2000 |
|---|---:|---:|---:|
| `lambda` [km], t=1000 s | 3.69 | 3.54 | 2.18 |
| `lambda` [cells], t=1000 s | 3.33 | 6.40 | 7.87 |
| `lambda` [km], t=500 s | 3.69 | 3.54 | 1.94 |
| `lambda` [cells], t=500 s | 3.33 | 6.40 | 7.00 |

Reading: ns=500 is Nyquist-limited (3.3 cells) so its `lambda` is meaningless as a physical scale. Between ns=1000 and ns=2000 `lambda` shrinks in km (3.54 -> 1.94-2.18) while staying at 6.4-7.9 cells. **The fine-structure component is a fixed-cell-count mode whose physical wavelength shrinks with the mesh — a grid-scale numerical mode — superposed on a fixed-physical-height envelope of dominant peaks.** That combination is exactly what a truncation error driven by a fixed physical gradient looks like.

A periodogram-based `lambda_peak` was also computed and is reported in `resconv_report.txt`, but it is **not trustworthy here** (18-72 samples give ~3 usable Fourier bins and it locks onto the slow trend rather than the ripple; its values jump between 5 and 20 km non-monotonically). It should be ignored in favour of the extremum statistics and the peak matching above.

## 6. Mean-flux non-convergence (Result D), in detail

`mean(F_eff)` and `mean(rho V)` agree to 0.3-1.7 % at every resolution — the *mean* is not in dispute between the two diagnostics. What is in dispute is its value.

Richardson with mesh ratio 2, from three successive values:

| quantity | ns=500 | ns=1000 | ns=2000 | `p_obs` | `Delta h -> 0` limit | rel. error at ns=2000 |
|---|---:|---:|---:|---:|---:|---:|
| `mean(F_eff)`, t=1000 s | 3.6922e-09 | 2.3999e-09 | 1.8417e-09 | +1.21 | 1.417e-09 | **+29.9 %** |
| `mean(rho V)`, t=1000 s | 3.6438e-09 | 2.4088e-09 | 1.8425e-09 | +1.13 | 1.363e-09 | **+35.2 %** |
| `q_top` (diag.), t=1000 s | 5.012 | 3.285 | 2.339 | +0.87 | 1.195 | +95.7 % |
| `mean(F_eff)`, t=500 s | 4.3458e-09 | 2.6978e-09 | 1.8222e-09 | +0.91 | 8.295e-10 | +119.7 % |
| `q_top` (diag.), t=500 s | 5.751 | 3.649 | 2.311 | +0.65 | **-0.039** | (meaningless) |

**These extrapolations are themselves unreliable and are quoted only to bound the error.** The `t=500 s` limit for `mean(F_eff)` (8.3e-10) and the `t=1000 s` limit (1.42e-9) differ by 70 %, and the `t=500 s` `q_top` limit comes out *negative* — both symptoms of applying Richardson to a state that is still relaxing rather than converged in time. The defensible statement is the raw one: **`mean(F_eff)` falls by 1.54x then 1.30x per halving of `Delta h`, roughly `O(Delta h)`, with no flattening. The evaporation rate has no grid-independent value in this scan, and the ns=2000 value is plausibly still tens of percent above the continuum limit.**

The driver is visible in the boundary state, and it is the *thermal* boundary, not the mass boundary:

| quantity, t=1000 s | ns=500 | ns=1000 | ns=2000 |
|---|---:|---:|---:|
| `T_top` [K] | 19 526.0 | 20 695.4 | 21 361.8 |
| `T_wall - T_top` [K] | 2474 | 1305 | 638 |
| `q_top` [W m^-2] (diagnostic) | 5.012 | 3.285 | 2.339 |
| `E_cond(1000 s)` [J m^-2] | 5513.3 | 3328.8 | 2120.4 |
| `V_top` [m/s] | 118.43 | 81.91 | 62.98 |
| `V_top / V_{N-2}` | 1.269 | 1.187 | **1.004** |
| `p_top` [Pa] | 0.010300 | 0.010293 | 0.010289 |
| `M_column` [kg m^-2] | 2.0974e-04 | 2.0978e-04 | 2.0980e-04 |
| `M/M(0) - 1` | -6.91e-04 | -4.78e-04 | -3.51e-04 |

- **`T_top` converges toward the 22 000 K wall** (19 526 -> 20 695 -> 21 362) and the wall-to-top-cell gap collapses (2474 -> 1305 -> 638 K, i.e. ~`O(Delta h)`). That is the expected behaviour of a first-order Dirichlet boundary layer.
- **The realized conductive input is therefore not resolution-converged either**: `E_cond` falls by 1.66x then 1.57x per halving. Since this is the *drive* for the evaporation, an unconverged drive guarantees an unconverged evaporation rate. Result D is a boundary-drive problem, not a continuity problem.
  Caveat: `q_top` here uses the same first-order estimator as `util/upper_bc_decouple_diag.py` (`kappa_solver(top) * (T_wall - T_top)/Delta h`, cell conductivity because the sidecar lacks the ghost density), so its *absolute* value is diagnostic-grade; the resolution *trend* is what matters and it is unambiguous.
- **Reassuringly, the top-cell anomaly converges away.** `V_top/V_{N-2}` = 1.269 -> 1.187 -> **1.004**: at ns=2000 the last cell simply continues the interior, so the residual last-cell excess left by the hydro-T decoupling stage was itself an `O(Delta h)` artefact. No velocity reversal at any resolution.
- **No drainage, no pressure runaway, at any resolution.** `p_top` is constant to 4 digits across the scan, `M/M(0)-1` is `<= 7e-04` and *shrinks* with refinement, and all three runs are in the same physical regime (subsonic upward outflow at the top, positive mass flux throughout the window). Nothing fell into a different regime.

## 7. Resolution required for a given ripple amplitude

Extrapolated from the finest point (`A_M = 0.0665` at `Delta h = 0.2765` km), using both the 3-point global fit and the (more trustworthy) local mid->fine rate:

| target | global fit (`p = 0.61`) | local rate (`p = 0.93`) |
|---|---|---|
| `A_M < 5 %` | `Delta h` = 0.173 km, ns ~ 3200 | `Delta h` = 0.203 km, **ns ~ 2700** |
| `A_M < 2 %` | `Delta h` = 0.038 km, ns ~ 14 500 | `Delta h` = 0.076 km, **ns ~ 7300** |

So a 5 % ripple is within reach (ns ~ 2700-3200, ~1.5x the cost of the ns=2000 run already performed); 2 % is not, at ~ns 7000-14 000. `Q_M` already reaches 2.1 % at ns=2000.

## 8. 500 s vs 1000 s

**The two times give the same trend at every resolution and for every core metric** — see `visualization/model_column/resconv_time_compare.png`. `A_M` fit `p` = +0.58 (500 s) vs +0.61 (1000 s); `Q_M` +1.47 vs +1.45; `A_diff` +0.62 vs +0.69; `A_mismatch` +0.59 vs +0.61. The mid->fine ratios agree to within 7 % between the two times. Nothing about the convergence conclusion is time-specific.

The mean flux differs between the two times at fixed resolution (e.g. ns=1000: 2.698e-09 -> 2.400e-09, -11 %), i.e. the runs are still slowly relaxing at 1000 s. That is the same drift already recorded in the previous stage and is the reason the Richardson extrapolation of the mean flux cannot be trusted (§6).

## 9. Answers to the stage's 14 required questions

1. **Identical physical configuration?** Yes for everything controllable: domain, IC construction, gamma/Saha table, decoupled hydro-T BC, back-pressure semantics (anchor agrees to 0.03 %), conduction wall, cooling off, TRAC off, MC3 `beta=2`, CFL 0.25, Mach cap, `eq_wb` construction, and the production `Delta h`-scaling numerical diffusivity. Two documented exceptions: (a) the top-cell IC temperature is resolution-dependent at `O(Delta h)` because the IC is point-sampled from a steep pchip profile (21 574 / 22 275 / 22 668 K) — unavoidable, and a candidate contributor to Result D; (b) `time_mult` 30 -> 60 to lift the step cap, proven byte-neutral.
2. **Production-identical face-flux diagnostic reused?** Yes, unchanged. No production C++ was modified this stage. `f_central + f_diff = f_total` to 1.9e-06 and the `t=0` reference reproduces `eq_residual_mass` to 1.7e-06 at all three resolutions.
3. **`A_M` vs `Delta h`?** Monotone decreasing: 0.154 -> 0.126 -> 0.066 (1000 s), 0.150 -> 0.137 -> 0.067 (500 s). Local mid->fine ratio 1.90-2.04 for a 2x refinement, i.e. **first order**. The 3-point fit gives `p ~ 0.6` only because the ns=500 point is Nyquist-limited.
4. **`A_diff` vs `Delta h`?** Monotone decreasing at the same rate: 0.179 -> 0.128 -> 0.069 (1000 s), ratio 1.87 on the clean pair. `mean(F_diff)` collapses much faster (4.1x then 14.1x).
5. **Does `A_mismatch` go to zero?** Yes. 0.154 -> 0.126 -> 0.066, tracking `A_M` to within 2 % at every resolution, at ~first order. **This is the boxed core question of the stage and the answer is yes.**
6. **Is `F_eff` smooth at all resolutions?** Yes. `Q_eff` = 1.3-1.9e-4 at every resolution and both times, against `Q_M` = 2.1e-2 to 1.6e-1; `A_eff` stays at 0.3-1.1 % and is a monotone rise, not an oscillation; `R_eff` normalised falls at ~first order. No sign of Result C.
7. **Does the mean evaporation mass flux converge?** **No.** It falls ~linearly in `Delta h` with no flattening (1.54x then 1.30x per halving at 1000 s), and the realized conductive drive falls with it. Result D holds.
8. **Ripple wavelength: fixed physical or fixed cell?** **Both components are present.** The dominant extrema are pinned to fixed physical heights (10/10 and 11/11 coarse extrema matched to fine-mesh counterparts within 0.43-0.70 km), but refinement inserts additional extrema so that the mean spacing is fixed at 6.4-7.9 **cells** and the wavelength shrinks in km (3.54 -> ~2.0 km). A fixed-physical envelope carrying a grid-scale mode.
9. **Do 500 s and 1000 s agree?** Yes, on every core metric and every fitted order, to within 7 % on the mid->fine ratios.
10. **Classification** — see §10 below.
11. **Is `N = 1000` enough for the current science?** **Split answer.** For the *structural* claim (the top-region `rho V` ripple is spatial truncation error and not a boundary or continuity defect) `N = 1000` is sufficient, and the scan is the evidence. For any *quantitative* claim — evaporation mass flux, top-cell velocity, conductive input, top-cell temperature — `N = 1000` is **not** sufficient: those quantities move by 30-80 % between ns=1000 and ns=2000 and are still moving. They should not be quoted from an ns=1000 run without a resolution error bar.
12. **Proceed to the Riemann-solver comparison next stage?** **Not for the ripple** — the ripple is now explained and converging, so replacing Rusanov would be solving a solved problem. **The next stage should target Result D instead:** the resolution-dependence of the realized conductive boundary input and of the top-cell Dirichlet layer. That is where the unconverged physics is.
13. **Paths** — §11.
14. **Tests / neutrality evidence** — §12.

## 10. Classification

### Confirmed
1. The cell-centred `rho V` ripple in 2130-2150 km **converges** under mesh refinement: `A_M` 0.154 -> 0.126 -> 0.066 and `Q_M` 0.160 -> 0.068 -> 0.021 at 1000 s, monotone, same trend at 500 s.
2. `A_mismatch = RMS[rho V - F_eff_cell]/|mean(F_eff)|` **goes to zero** at ~first order, tracking `A_M` to within 2 %.
3. `F_eff` is smooth at every resolution (`Q_eff ~ 1.4e-4`, flat) and its continuity residual `R_eff` normalised **improves** with refinement (2.2e-3 -> 6.9e-4). **Result C is excluded.**
4. The anti-correlated central/diffusive compensation persists unchanged in *character* (`corr` -0.992 to -0.999) and shrinks in *amplitude* with the mesh.
5. The dominant ripple extrema are pinned to fixed physical heights across all three meshes; the additional fine-scale extrema are a fixed-cell-count mode.
6. **`mean(F_eff)` does not converge** — ~`O(Delta h)` with no flattening — and neither does the realized conductive input `E_cond` (5513 / 3329 / 2120 J m^-2) or `T_wall - T_top` (2474 / 1305 / 638 K). **Result D.**
7. The residual last-cell velocity excess left over from the hydro-T decoupling stage is itself `O(Delta h)`: `V_top/V_{N-2}` -> 1.004 at ns=2000.
8. No drainage, no pressure runaway, no regime change, no velocity reversal at any resolution; mass conservation *improves* with refinement.
9. Enabling the face-flux capture is numerically neutral at ns=500 and ns=2000 as well as ns=1000 (§12).

### Most likely explanation
The `rho V` ripple is **spatial truncation error of the Rusanov dissipative mass flux against the steep TR density/temperature gradient**, as hypothesised at the end of the previous stage. The scan supports the hypothesis on all four of its predictions: the amplitude falls with `Delta h` at ~first order; the diffusive half falls with it; the conserved flux the update actually differences stays smooth throughout; and the fine structure is a grid-scale mode riding on a fixed-physical envelope. Per the pre-registered rules this is **Result A**, and **no boundary change is indicated by the ripple**.

Separately and independently, the **unconverged evaporation rate is most likely a first-order Dirichlet conduction boundary layer at the top face**: `T_top` marches toward the 22 000 K wall as `O(Delta h)`, so the wall-to-top-cell gap — and hence the realized conductive drive — is mesh-dependent. This is a hypothesis about the *cause*; only the *fact* of non-convergence is confirmed.

### Not verified
- **Why the mean flux does not converge.** The Dirichlet-layer explanation above is untested. Competing/compounding candidates, all untested: the resolution-dependent IC top-cell temperature (§1); `numerical_diffusivity ∝ Delta h` changing the Stage-D operator between runs (kept scaling by instruction, but it is a real physical difference between the three runs); the `Delta h` in the conductive-flux stencil at the wall; domain truncation inside the steep gradient.
- **Whether any of these three runs is converged in time.** All are still relaxing at 1000 s (`mean(F_eff)` drops 11 % from 500 s to 1000 s at ns=1000). This is what makes the Richardson extrapolations of §6 untrustworthy and is why only the raw ratios are claimed.
- **`Q_eff`'s floor.** `Q_eff ~ 1.4e-4` does not decrease with refinement. It is far too large to be double round-off (the capture is `double` throughout), so something sets a fixed normalised-curvature floor in `F_eff`. Harmless at this magnitude but unexplained.
- **Whether `A_M` continues at first order below `Delta h = 0.2765` km.** The claim rests on a single clean pair (ns=1000 -> ns=2000). An ns=4000 point would settle it and would also test the ns~2700 prediction for `A_M < 5 %`.
- **Energy and momentum rows.** Still only the mass row has been decomposed; whether the same compensation degrades the heat-flux and velocity diagnostics is untested (carried over from the previous stage).
- **Whether a less dissipative Riemann solver shrinks the compensation.** Untested, and now de-prioritised (§9.12).

### Explicitly not done
No change to the Rusanov flux, the limiter, the EOS, the well-balanced subtraction, any boundary condition, the conduction operator, the numerical-diffusivity formula, or the integration order. No HLL/HLLC. No NSCBC, sponge or pressure relaxation. No baseline fixed-hydro-`T` companion runs at the new resolutions (the stage instruction said to skip them unless an unexplained result appeared; none did). No re-run of the limiter or CFL legs. No production default changed.

## 11. Files

New:
- `util/resolution_scan_diag.py` — the analysis script (metric and mapping definitions in its docstring)
- `docs/studies/validation/resolution_convergence_scan_recap.md` — this document
- `outputs/model_column/resconv_ns{500,1000,2000}_dec_1000s.txt{,.faceflux,.gamma_diag}` + `.console.log`
- `outputs/model_column/resconv_report.txt`, `outputs/model_column/resconv_metrics.json`
- `visualization/model_column/resconv_profiles_{500,1000}s.png` — 4 panels vs physical height: cell `rho V`; `rho V / mean(F_eff)` (ripple convergence independent of the mean shift); `F_eff`; `delta F_diff`
- `visualization/model_column/resconv_convergence_loglog.png` — `A_M`, `Q_M`, `A_diff`, `A_mismatch`, `A_eff` vs `Delta h`, with `Delta h` and `Delta h^2` guides, both times
- `visualization/model_column/resconv_mean_flux.png` — mean mass flux and normalised `R_eff` vs resolution
- `visualization/model_column/resconv_time_compare.png` — 500 s vs 1000 s convergence side by side

Modified: `visualize_commands.md` (new "Resolution convergence of the top-region ripple" section with the exact commands).

**No source file under `src/`, `scenarios/`, `chromo_main.cpp`, `chromosphere.hpp` or `tests/` was modified.** Exact commands: the new section of `visualize_commands.md`.

All figures use physical height as the x-axis; none uses cell index.

## 12. Tests and diagnostic neutrality

- **Full suite: 7682 passed, 0 failed** — identical to the previous stage's baseline, as expected since no production code changed. Includes `test_face_flux_capture_matches_production_continuity`, which is what certifies the capture is the production flux.
- **Diagnostic neutrality at the new resolutions.** 20 s runs at ns=500 and ns=2000 with `CHROMO_FACE_FLUX_DIAG=0` vs `=1`, everything else identical:

| | `.gamma_diag` | state `.txt` |
|---|---|---|
| ns=500 | **byte-identical** | identical except the embedded output filename in a header comment (line 3: `diagnostics=...diag0.txt.gamma_diag` vs `...diag1.txt.gamma_diag`); all 1006 lines of data identical |
| ns=2000 | **byte-identical** | same, all 4006 lines of data identical |

  So the capture is numerically neutral at 0.2765 km and 1.106 km as well as at the 0.553 km already covered by the unit test.
- **`time_mult` neutrality.** ns=500 and ns=1000 re-run at `time_mult=60.0` reproduce the `time_mult=30.0` runs **byte-for-byte** in both the `.gamma_diag` and the `.faceflux` sidecar, confirming that `step_cap` is the only thing that argument touches.

## 13. Suggested next stage (not executed)

The ripple question is closed. The open question is Result D.

1. **Isolate the top conduction boundary layer.** Repeat the scan with the *realized* conductive input held fixed instead of the wall temperature (e.g. `impose_outer_heat_flux` with a matched `q`), and see whether `mean(F_eff)` then converges. This directly tests the Dirichlet-layer hypothesis and reuses machinery that already exists.
2. **Separate the IC sampling effect.** Build the IC by cell *averaging* rather than centre point-sampling of the C7 pchip profile, so the initial top-cell temperature is not `O(Delta h)`-dependent, and re-measure. Cheap, and removes one confound from Result D.
3. **Add ns=4000** to confirm `A_M ~ O(Delta h)` continues and to test the ns~2700 prediction for `A_M < 5 %`.
4. **Reach a genuine steady state** (or state explicitly that these are relaxing solutions) before quoting any evaporation rate; the 500 s -> 1000 s drift currently invalidates Richardson extrapolation.
5. Riemann-solver comparison: **de-prioritised**, since the mechanism it would test is now confirmed to be converging.
