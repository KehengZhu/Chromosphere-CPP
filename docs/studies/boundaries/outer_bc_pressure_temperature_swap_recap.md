# Outer face: swapping which variable carries the one imposed condition

**Question.** The release outer face at 2153 km imposes a **fixed reservoir back-pressure** `kOuterPRef` and lets the hydro ghost temperature free-float with the live top cell (`ISO_HYDRO_T_DECOUPLE=1`), with the 22 kK conduction wall a separate physical-face datum. Subsonic outflow admits exactly one incoming characteristic, so exactly one condition may come from outside. What happens if the *other* variable carries it — fixed `T_hydro = T_cond = 22` kK and a **Neumann (loose) pressure**?

**Answer.** It depends entirely on *which* Neumann pressure, and neither variant is usable.

* **`dp/ds = 0` (the literal zero-gradient BC) destroys the model.** It withholds the hydrostatic pressure support the top cell requires — the interior needs `dp/ds = -rho g` — so the boundary applies a standing net downward force. The column collapses inward, the whole domain runs at the `ISO_VCAP` ghost Mach cap (`V_top = -2462.5` m/s against `0.1 c_s = 2460.3` m/s, i.e. **the steady velocity is the clamp, not physics**), **half the column mass is gone** by 1000 s, and the 6.7 kK chromosphere is replaced by a 23.7 kK fully ionized (`x = 1.000`) pass-through downdraft with `f_top = f_base = -4.13e-07` kg m⁻² s⁻¹, 1530x the release flux and pointing the wrong way. The run does not crash: it reaches a perfectly steady, perfectly wrong state and reports `termination=end_time` with `godunov.fallbacks=0`. **A clean termination is not evidence that a boundary closure is sound.**
* **`dp/ds = -rho g` (hydrostatically consistent, anchored on the live top cell) survives 4000 s but suppresses the physics of interest.** Nothing pins the pressure *level*, so `p_top` drifts up **+5.8 %** over 4000 s and is still rising, column mass drains −1.03 % (release: −0.006 %), and the conduction-driven evaporation signal collapses by a factor of ~19: `V_top = +0.49` m/s against the release `+9.50` m/s, top-face `f_total = +1.47e-11` against `+2.70e-10` kg m⁻² s⁻¹, with the top and base fluxes still disagreeing in sign at 4000 s (no steady state).

**Both failures are entirely the pressure closure.** The control that pins the hydro ghost temperature and keeps the release back-pressure is indistinguishable from the release: `V_top` 9.5008 → 9.5148 m/s, top-face `f_total` 2.7023e−10 → 2.7025e−10, column mass identical to five digits, and the largest velocity difference anywhere in the column is 0.14 m/s at 4000 s (0.12 m/s at its worst, near 100 s). Fixing `T_hydro = T_cond` costs nothing here; loosening the pressure costs everything.

The release closure is therefore retained. The value of the experiment is the mechanism: **the outer pressure closure has to be hydrostatically consistent first, and externally anchored second, and the two requirements fail differently.**

Evidence: `outputs/model_column/outer_bc_swap/`, movie `visualization/model_column/outer_bc_closure_swap_4000s.mp4` (`util/animate_outer_bc_compare.py`, catalogued in `visualize_commands.md`).

## What was varied, and what was not

Two new default-off knobs in `scenarios/model_column.cpp`, both diagnostic only:

| knob | effect |
|---|---|
| `ISO_OUTER_T_HYDRO_WALL=1` | In decoupled mode, pins the hydro ghost temperature to the conduction wall (`T_hydro = T_cond = a·T_ref`) instead of zero-gradient-extrapolating the live top cell. The wall **stays a physical-face conduction datum**, unlike baseline `ISO_HYDRO_T_DECOUPLE=0`, where it reverts to a ghost-centre datum. |
| `ISO_OUTER_P_NEUMANN=1` | Outer ghost pressure `p_g1 = p_g0 = p_top`: `dp/ds = 0` through the face and the ghost region. |
| `ISO_OUTER_P_NEUMANN=2` | The same trapezoidal EOS-closed hydrostatic ladder (`hse_rung`) the release uses, re-anchored on the **live top cell** instead of `kOuterPRef`: `dp/ds = -rho g` with no external pressure datum. |

`ISO_OUTER_P_NEUMANN` is rejected outright when the hydro ghost temperature is still free (release mode), because the face would then impose *nothing* — pressure, velocity and temperature all extrapolated — and the problem is ill posed rather than merely wrong.

Everything else is the canonical release configuration: N = 500/R4 (661 cells), `h_base = 1600` km, `dh = 553` km, production `Gamma1` table, SWMF exact-Riemann Godunov flux on `(ln rho, V, ln p)` with the MC3 limiter, physical conduction only, 22 kK wall, `CHROMO_CFL=0.50`, 4000 s, **double-precision state** (`build_omp`). The release baseline is the 22 kK run of the wall-temperature sweep, `outputs/model_column/twall_sweep/twall_22kK_cfl0.50.txt`, so the comparison shares its IC, mesh, flux and cadence exactly. All four runs reported `termination=end_time` and `godunov.fallbacks=0`; `kOuterPRef = 0.0102856` Pa in every case.

## Result at 4000 s

| | release | T pinned only | T pinned + `dp/ds=0` | T pinned + `dp/ds=-rho g` |
|---|---:|---:|---:|---:|
| `ISO_OUTER_T_HYDRO_WALL` / `ISO_OUTER_P_NEUMANN` | 0 / 0 | 1 / 0 | 1 / 1 | 1 / 2 |
| column mass `M/M(0)` | 99.994 % | 99.994 % | **49.943 %** | 98.973 % |
| `V_top` [m/s] | +9.501 | +9.515 | **−2462.46** | +0.49 |
| top-face `f_total` [kg m⁻² s⁻¹] | +2.7023e−10 | +2.7025e−10 | **−4.1324e−07** | +1.4737e−11 |
| base-face `f_total` | +2.7096e−10 | +2.7098e−10 | −4.1324e−07 | −2.4234e−10 |
| `T` at 1900 km [K] | 6654 | 6654 | **23694** | 6670 |
| `x_eq` at 1900 km | 0.458 | 0.458 | **1.000** | 0.468 |
| `p_top` drift over 4000 s | ~0 % | ~0 % | +493 % | **+5.77 %, still rising** |
| `f_top / f_base` (mass balance) | 0.997 | 0.997 | 1.000 | **−0.061** |

### The `dp/ds = 0` collapse, in time

| t [s] | `M/M(0)` | `V_top` [m/s] | `f_top` | `f_base` | `T`@1900 km | `x`@1900 km |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 100.00 % | 0.0 | −1.06e−11 | +1.25e−13 | 6661 | 0.462 |
| 20 | 100.06 % | −614.2 | −1.83e−08 | +1.09e−12 | 6661 | 0.462 |
| 60 | 101.39 % | **−2462.0** | −9.76e−08 | +1.53e−09 | 6723 | 0.484 |
| 100 | 97.00 % | −2461.7 | −1.15e−07 | −7.04e−07 | 6820 | 0.558 |
| 200 | 79.95 % | −2412.0 | −1.84e−07 | −6.00e−07 | 30242 | 1.000 |
| 400 | 38.23 % | −2447.1 | −3.84e−07 | −7.86e−07 | 27856 | 1.000 |
| 1000 | 49.94 % | −2462.5 | −4.13e−07 | −4.13e−07 | 23694 | 1.000 |
| 4000 | 49.94 % | −2462.5 | −4.13e−07 | −4.13e−07 | 23694 | 1.000 |

Reading it: the ghost velocity saturates at the `ISO_VCAP` cap within ~60 s, mass then leaves through the **base** faster than it enters at the top (the column is being pushed out of the bottom of the domain), the transition region marches down from 2153 km to below 1900 km between 100 s and 200 s, and from ~700 s the whole thing is stationary — a fully ionized 23.7 kK conduit carrying a uniform downward flux. From 1000 s to 4000 s not one digit of the state moves.

## Attribution

`ISO_VCAP` sets the collapse speed, and only the speed: `0.1 c_s = 2460.3` m/s at the final top-cell state against the observed `|V_top| = 2462.5` m/s. Removing the cap would not rescue the closure, it would make the collapse faster; the cap is what makes the failure *look* like a steady state.

Because `dp/ds = -rho g` (mode 2) is stable and `dp/ds = 0` (mode 1) is catastrophic under an otherwise identical configuration, **the mode-1 failure is the missing hydrostatic slope, not the missing external anchor.** `g = 273.95` m/s² and the ghost withholds the half-cell hydrostatic support `rho g Δs/2` at the top face, permanently.

The `T`-pinned-only control (`ISO_OUTER_T_HYDRO_WALL=1` with the release back-pressure, `outputs/model_column/outer_bc_swap/tfix_pref_4000s.txt`) assigns **all** of both failures to the pressure closure. It tracks the release throughout, not just at the end:

| t [s] | `V_top` release | `V_top` T pinned | max column-wide `|ΔV|` |
|---:|---:|---:|---:|
| 20 | 40.6405 | 40.7157 | 1.0e−01 |
| 100 | 39.5134 | 39.6342 | 1.2e−01 |
| 1000 | 13.9340 | 13.9616 | 2.8e−02 |
| 4000 | 9.5008 | 9.5148 | 1.4e−02 |

This does **not** retire the earlier decoupling study (`upper_bc_hydro_temperature_decoupling_recap.md`), which measured a real last-cell pathology from the fixed hydro temperature at N = 1000 without outer refinement. It is configuration-specific: in the current release configuration the imposed wall (22 000 K at the physical face) and the live top cell (22 675 K at t = 0) agree to 3 %, so pinning the ghost to the wall is a small perturbation. Read it as "the hydro-T half of this swap is inert *here*", not as "fixing the hydro temperature is harmless".

## The historical "live top-cell pressure anchor" result is not reproducible

`main.tex` and the outer-face comment block in `scenarios/model_column.cpp` record a retired experiment — "use the current top pressure and one-sided HSE outside it" — whose measured outcome was that the column *drained upward*, accelerated to about 1.9 km s⁻¹ and failed near 12 s. That configuration is exactly `ISO_OUTER_P_NEUMANN=2`, and in the current release numerics it does not fail: 4000 s, `termination=end_time`, no fallbacks, mass drift −1.03 %.

Four cheap 200 s controls narrow it, none of them explaining it:

| control | result |
|---|---|
| release Godunov flux | stable, `max|V|` decaying 227 → 11 m/s |
| previous release flux (`ISO_RIEMANN=roe-local`) | stable, indistinguishable (`p_top` differs in the 5th digit) |
| `float32` state (`build_omp_f32`) | stable, indistinguishable |
| baseline hydro-T mode (`ISO_HYDRO_T_DECOUPLE=0`, wall at the ghost centre) | stable |

So neither the face flux, nor the storage precision, nor the conduction-datum geometry is responsible. The remaining candidates are the R4 outer refinement, the reference-free predictor momentum source plus the second-order lower ghost ladder, and the fact that `hse_rung` was first-order at the time of the original experiment. **Not resolved here.** The practical consequence is narrow and should be stated as such: the *stability* claim attached to that retired row is stale, while the *reason it stays retired* — an unanchored pressure level drifts secularly and the evaporation signal it is meant to measure collapses — is now measured directly and is stronger than the original.

## Wall-temperature study under the unanchored closure: 22 kK vs 100 kK

The obvious follow-up is whether the surviving closure (`ISO_OUTER_P_NEUMANN=2`) can still be used for the one measurement `model_column` is normally used for — the response to the conduction-wall temperature. It cannot, and the way it fails is more interesting than a simple offset: **the wall-temperature *scaling* survives almost intact while the *calibration*, the mass budget and the chromosphere itself do not.**

Both walls were run at **CFL 0.25**, because 100 kK fails the implicit conduction Newton solve at CFL 0.50 under this closure exactly as it does under the release closure. The release-closure counterparts are the 22 kK (CFL 0.50) and 100 kK (CFL 0.25) runs of the wall sweep, joined the same way the sweep joins its two CFL groups. The `ISO_OUTER_P_NEUMANN=2` 22 kK run also exists at both CFL values, so the join is measured rather than assumed: `V_top` 0.500 (CFL 0.50) against 0.542 m s⁻¹ (CFL 0.25), i.e. **8 %**, which bounds the CFL contribution to everything below.

Quasi-steady means the mean over the last 500 s; `f` is the conservative face flux at 4000 s.

| | release 22 kK | release 100 kK | unanchored 22 kK | unanchored 100 kK |
|---|---:|---:|---:|---:|
| `V_top` quasi-steady [m/s] | +9.645 | +791.86 | +0.542 | +37.73 |
| drift across the window | −1.6 % | −2.9 % | −1.6 % | −2.8 % |
| top-face `f_total` [kg m⁻² s⁻¹] | +2.702e−10 | +4.942e−09 | +1.599e−11 | +5.422e−10 |
| base-face `f_total` | +2.710e−10 | +4.552e−09 | **−2.388e−10** | **−3.758e−09** |
| `f_top / f_base` | +0.997 | +1.086 | **−0.067** | **−0.144** |
| column mass at 4000 s | 99.994 % | 99.887 % | 98.978 % | **77.272 %** |
| `p_top(4000)/p_top(0)` | 1.000 | 1.000 | 1.057 | **2.276, still rising** |
| `T` at 1900 km [K] | 6654 | 6546 | 6670 | **21931** |

Three readings:

1. **The exponent survives; the amplitude does not.** The two-point exponent of `V ~ T_wall^n` is **2.911** under the release closure and **2.802** under the unanchored one — against the sweep's fitted 2.91 over 15–100 kK. But the suppression factor is 17.8x at 22 kK and 21.0x at 100 kK. Anyone who measured only the *shape* of the wall response under this closure would conclude, correctly, that it is a `T_wall^~3` power law, and would be wrong about every absolute number by a factor of ~20. Two points cannot establish a power law, so read `n` as a consistency check against the sweep, not as a new fit.
2. **The mass budget inverts.** Under the release closure the top and base face fluxes agree to 0.3 % (22 kK) and 8.6 % (100 kK): one coherent upward column. Under the unanchored closure they have *opposite signs* — a small evaporative flux leaves the top while a flux 15x to 7x larger drains out of the base. The column is not evaporating; it is emptying downward with an evaporative trickle on top, and at 100 kK it has lost **22.7 %** of its mass by 4000 s with no sign of stopping.
3. **At 100 kK the chromosphere is destroyed, and this is a wall-temperature-dependent failure, not a constant offset.** With the pressure anchored, a 100 kK wall drives 792 m/s of evaporation while 1900 km stays at 6546 K — the thermal response is confined to the resolved TR, exactly as the sweep reported. With it unanchored, the transition region migrates *down* to about 1900 km and the temperature there reaches **21 931 K**. The pressure level runs away first (+128 % and still rising linearly at 4000 s, against +5.7 % at 22 kK), and the hot wall is what makes the runaway fast enough to matter over 4000 s.

Figure: `visualization/model_column/outer_bc_twall_response.png` (`util/outer_bc_twall_diag.py`; metrics in `outputs/model_column/outer_bc_swap/outer_bc_twall_metrics.json`). Movie: `visualization/model_column/outer_bc_twall_22kK_vs_100kK.mp4` (`util/animate_outer_bc_compare.py --preset twall`), colour = wall temperature, dashed = unanchored pressure.

**Conclusion for anyone tempted to use this closure because it is "less imposed":** it is stable, it reproduces the qualitative wall-temperature scaling, and it will silently give you evaporation numbers ~20x too small against a pressure reference that is itself moving. The externally anchored back-pressure is not a convenience; it is what makes a wall-temperature measurement mean anything.

## Status

* Release closure unchanged. Both knobs default off, are diagnostic only, and no release number may be quoted from them.
* The unanchored closure is unusable for parameter studies, and its failure grows with the wall temperature — 17.8x evaporation suppression and −1.0 % mass at 22 kK, 21.0x and −22.7 % with the chromosphere gone at 100 kK.
* `ISO_OUTER_P_NEUMANN=1` is a **negative result to keep**: it is the cleanest demonstration in this repository that the outer face must be hydrostatically consistent, and that a run can terminate normally, conserve nothing, and look steady.
* Open: why the historical live-anchor drainage does not reproduce. Cost to resolve is one bisection over the mesh/predictor-source/ladder changes; nobody needs the answer for a release decision.
