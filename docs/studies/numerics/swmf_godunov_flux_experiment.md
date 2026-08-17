# SWMF-style Godunov flux — implementation, Roe comparison, and release cutover

**Status: RELEASED.** The SWMF-style exact-Riemann Godunov flux is the `model_column` release numerical flux. The mixture Roe characteristic flux it replaced is retained as the reference/comparison solver, reachable with `ISO_RIEMANN=roe-local`.

Sections 1–4 are the original controlled experiment, written while the flux was still a comparison mode, and are preserved because they are the evidence base. Section 5 records the limitations as they stand for the release. **Section 8 is the release cutover** — the rationale, the CFL audit, the acceptance runs, and the corrections this document had to make to its own earlier conclusions. Where sections 1–4 and section 8 disagree, section 8 is current.

## 1. What was implemented, and from what

The reference implementation is SWMF, read directly:

- `SWMF/util/CRASH/src/test_godunov.f90` — `get_godunov_flux`, the flux assembly and the energy-offset device (this is the routine the advisor pointed to).
- `SWMF/share/Library/src/ModExactRS.f90` — `exact_rs_pu_star`, `exact_rs_sample`, `pressure_function`, `guess_p`: Toro's exact ideal-gas Riemann solver (Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 2nd ed., ch. 4) in the SWMF revision that permits `gamma_L != gamma_R` and uses P. Voinovich's impedance-weighted pressure guess.

`ModExactRS.f90` is ported line for line to `src/single_fluid/exact_rs.{hpp,cpp}`. The only deviations are (a) the Fortran's OpenMP-threadprivate module variables become a returned `ExactRsSolution`, so the solver is reentrant per face, and (b) the port always behaves as the SWMF caller that passes `UseAnotherRS`: it reports a vacuum or non-positive star pressure instead of calling `CON_stop` or silently retrying at `0.1 p_old`. Arithmetic, branch structure, iteration budget (`nIterMax = 10`) and tolerance (`TolP = 1e-3`) are unchanged.

## 2. How the mode works

`build_swmf_godunov_flux` in `src/single_fluid/mixture.cpp` transcribes `get_godunov_flux` onto the release face states. Per face:

1. Take `(rho, u, p)` from the left and right face states the Roe flux would have received — the same MUSCL-reconstructed, Hancock-predicted, equilibrium-face-builder output. Nothing upstream changes.
2. Solve the **ideal-gas** Riemann problem at a fixed `gamma = 5/3`. This is SWMF's choice: `test_godunov.f90` passes `GammaMax = 5/3` on both sides.
3. Carry the non-ideal EOS content as a passive **specific energy offset** taken from the upwind side (`u* > 0` selects left, SWMF's tie-break at exactly zero selects right):

   `E0 = (E - 1/2 rho u^2 - p/(gamma - 1)) / rho = (e_int - p/(gamma - 1))/rho + phi`

   Here `e_int` already contains the Saha ionization energy, and this model's total energy additionally carries the gravitational potential, so both ride in `E0`. Both sides of a face share `phi`, so the gravitational part cancels from the left/right offset difference and only the EOS content is actually upwinded.
4. Sample the self-similar solution at `x/t = 0` and assemble `F = (rho* u*, rho* u*^2 + p*, (E* + p*) u*)` with `E* = p*/(gamma - 1) + 1/2 rho* u*^2 + rho* E0`. When the sampled state equals a side state this reproduces that side's `E` exactly.
5. On any failure — inadmissible input, vacuum guard, non-positive star pressure, iteration budget exhausted, inadmissible sample — fall back to the **same face-local Rusanov flux the Roe path uses**, and count it.

The controlled variable is exactly the face Riemann solver. Reconstruction, the Hancock predictor, the gravity/pressure-area source, boundary conditions, the Saha EOS and the conduction operator split are untouched, and the three flux choices are mutually exclusive in the scenario.

Files changed: `src/single_fluid/exact_rs.hpp` + `.cpp` (new), `src/single_fluid/mixture.cpp` (the new flux, the dispatch branch, the capture split), `chromosphere.hpp` (`Grid::swmf_godunov_flux`, `GodunovFluxStats`, `Grid::godunov_stats`), `scenarios/model_column.cpp` (`ISO_RIEMANN=swmf-godunov`), `chromo_main.cpp` (the `godunov.*` tally, the sidecar header tag), `CMakeLists.txt`, `tests/chromo_tests.cpp`.

## 3. Verification of the solver itself

Done before any chromosphere run, as the standing method requires.

**Toro's five standard tests** (`gamma = 1.4`), star state versus the published values:

| Test | `p*` (ref) | rel. err | `u*` (ref) | abs. err | iters |
| --- | --- | --- | --- | --- | --- |
| 1 Sod | 0.30313 | 5.9e-7 | 0.92745 | 2.5e-6 | 4 |
| 2 "123" (near vacuum) | 0.00189 | 2.1e-3 † | 0.0 | 0.0 | 9 |
| 3 left blast | 460.894 | 5.5e-7 | 19.5975 | 3.3e-3 | 6 |
| 4 right blast | 46.0950 | 9.6e-7 | −6.19633 | 3.4e-5 | 6 |
| 5 collision | 1691.64 | 4.1e-6 | 8.68975 | 1.8e-4 | 4 |

† the reference is quoted to three digits; our 0.00189387 rounds to it.

**`x/t = 0` sampling.** Sod gives `rho = 0.426319` against Toro's tabulated `rho*_L = 0.42632`. Test 5 correctly returns the left *data* state, because both its waves move right (`wl = 0.79 > 0`). Sampling beyond `wl`/`wr` returns the data states unchanged.

**Degenerate cases, at chromospheric magnitudes and `gamma = 5/3`.** A stationary contact (equal `p`, equal `u`, factor-2 density jump) gives **identically zero** mass flux and `p* = p` to 1e-14. A uniform state in uniform motion is reproduced to round-off for both signs of `u`. The vacuum guard and inadmissible input are reported, not returned.

**Finite-volume shock tube with the full flux assembly.** A first-order Godunov solve of Sod to `t = 0.2`, using exactly the flux construction of `build_swmf_godunov_flux` including the offset:

| `E0` | n=100 | n=200 | n=400 | n=800 | L1 rate | max \|E0 drift\| | fallbacks |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1.464e-2 | 9.383e-3 | 6.029e-3 | 3.847e-3 | ~0.64 | 2.6e-16 | 0 |
| 3.7e2 | 1.464e-2 | 9.383e-3 | 6.029e-3 | 3.847e-3 | ~0.64 | 5.7e-14 (1.5e-16 rel.) | 0 |

`~0.64` is the textbook first-order-Godunov L1 rate on a solution containing a shock and a contact. The two rows agree to all printed digits: **the energy offset is exactly inert to the hydrodynamics**, which is the property the SWMF device claims.

Both checks are now permanent suite cases: `test_swmf_exact_riemann_solver` and `test_swmf_godunov_flux_mode`. Suite: **11633/11633 checks, 0 failures**; `ctest` **9/9**.

## 4. Roe versus `swmf-godunov` on the model

Both legs use identical settings and the same launcher; only `ISO_RIEMANN` differs. Both reported `termination=end_time` and **exactly the same step count**, so the comparison is step-for-step.

### 4a. Conduction-off hydrostatic leg (the Hancock+gravity test)

`ISO_HEAT_FLUX=0 CHROMO_T_END=20`, N500/R4, 661 cells, CFL 0.50, 3620 steps each.

| | Roe (then release) | `swmf-godunov` |
| --- | --- | --- |
| `max\|F^rho\|` at `t = 0` | 1.3345e-12 | 1.0758e-12 |
| `max\|F^rho\|` at `t = 20 s` | 1.1713e-11 | 1.0503e-11 |
| `max\|V\|` at 20 s | 0.03022 m/s @ 2090 km | 0.03194 m/s @ 2094 km |
| `time.rhs_s` | 0.946 | 1.097 (+16 %) |
| fallbacks | — | 0 of 4,785,640 |

The Roe column reproduces the documented release value of 1.335e-12 exactly. Both solvers sit at the same float32 round-off floor; the differences are not resolved by this measurement.

> **Precision context, added at the 2026-08-17 double-precision release cutover.** Every hydrostatic-defect number in this document was measured with the pre-cutover **float32** conserved state, so "the float32 round-off floor" was a real and binding floor at the time. The release now stores the state in `double`, which lowers that floor by ~6 orders of magnitude, so a re-measured defect is no longer floor-limited in the same way and the flux comparison above should not be re-read as "the two fluxes are indistinguishable" without redoing it. See `state_precision_release_cutover.md`. The flux *choice* is unaffected: it was made on frozen-composition thermodynamics, not on a measured defect advantage.

### 4b. Production-like release smoke

`CHROMO_T_END=100`, frame cadence 10 s, sidecars on, 17781 steps each.

| | Roe (then release) | `swmf-godunov` | Δ |
| --- | --- | --- | --- |
| mean V [m/s] | 14.956 | 14.909 | −0.31 % |
| min / max V [m/s] | 11.370 / 32.718 | 11.319 / 32.677 | — |
| `n(V < 0)` | 0 | 0 | — |
| `F_eff` [kg m⁻² s⁻¹] | 1.1165e-9 | 1.1128e-9 | −0.33 % |
| `Q_V` (velocity roughness) | 1.287e-3 | 1.298e-3 | +0.8 % |
| `Q_eff` | 9.614e-4 | 8.824e-4 | −8.2 % |
| `eps_p` median | 1.66e-7 | 2.42e-7 | +46 % |
| `T_top` [K] | 21837.4 | 21837.4 | 0 |
| `p_top` [Pa] | 0.010288 | 0.010288 | 0 |
| BL thickness [km] | 12.140 | 12.140 | 0 |
| `q_phys` [W m⁻²] | 0.7783 | 0.7783 | ~0 |
| column mass drift | −2.858e-4 | −2.963e-4 | +3.7 % |
| `time.rhs_s` | 5.119 | 5.785 | **+13.0 %** |
| whole-step measured work | ~13.9 s | ~14.6 s | +4.7 % |
| fallbacks | — | **0 of 23,506,482** | — |

Face-level difference, from the `.faceflux` sidecars: `max|F_god − F_roe| / max|F_roe|` grows from 1 % at `t = 10 s` to about 9 % at `t = 90 s`, so this is a genuine solver difference and not a no-op. The Godunov **dissipative** part of the mass flux is consistently 8–12 % *smaller* than Roe's during the evaporation phase. **No causal attribution is claimed for that.** The obvious hypothesis is contact resolution — the exact solver resolves an isolated contact exactly, whereas the Roe flux here is a general-EOS linearization that is not constructed to satisfy Property U and would retain a small residual contact dissipation — but the measurement is a magnitude comparison of `f_total - f_central` with no isolation of the contact field, and it sits against a competing effect that pushes the other way (the Godunov flux's frozen-`gamma` waves are ~1.24x faster, which would *raise* dissipation). Separating the two needs a dedicated test, which has not been run.

### 4c. Does SWMF's loose `TolP = 1e-3` matter here?

`godunov.iterations_max = 1` for all 23.5 million face solves in the 100 s run: every face in this model is near-linear, so the Voinovich acoustic guess plus one Newton step already meets the tolerance. Replaying 7271 real production face states against a Newton solve converged to 1e-15:

- `p*`: max relative error **6.0e-16** (round-off);
- `u*`: max relative error **1.1e-6**;
- resulting face mass flux: max relative error **1.1e-6**.

So the shipped tolerance costs nothing measurable in this regime. It would cost more in a strong-shock problem — Toro test 3 above shows `u*` carrying a 1.7e-4 relative error from the one-iteration lag.

### 4d. Release default at the time of the experiment

At this point the release default was still Roe: `scripts/release_validation.sh` passed with the no-override run identical to an explicit `ISO_RIEMANN=roe-local ISO_RECONSTRUCTION=lnrho-v-lnp` run, and `swmf_godunov=0` in the release sidecar header. Both statements are now inverted — see section 8.

## 5. Limitations observed

1. **The fixed `gamma = 5/3` is the defining choice of the method.** *(This item was written as a disqualifying limitation and is superseded by section 8.1, which reinterprets it as the intended frozen-composition closure.)* The equilibrium acoustic speed is `sqrt(Gamma1 p/rho)` with the CRASH-tabulated `Gamma1`, which in this column has median **1.090** (minimum 1.089 at 1885 km) and is below 1.5 over **94.9 %** of the cells. The frozen ideal-gas wave speeds this flux uses are therefore about **1.24x faster almost everywhere**. What remains a genuine limitation is that the two closures have not been distinguished *by measurement* on this model: at these Mach numbers no observable separates them.
2. ~~**The timestep is deliberately *not* taken from this flux.**~~ **Fixed at the cutover, and the arithmetic in the original item was wrong.** `mixture_timestep` now uses the frozen signal speed whenever the Godunov flux is in force. The original claim that `CHROMO_CFL=0.50` gave an effective CFL of about 0.62 assumed the CFL-limiting cell carried `Gamma1 ~ 1.09`; it does not. See section 8.2 for the measurement.
3. **Only smooth, low-Mach behaviour has been exercised on the model.** *(Still current, and now a release limitation.)* `swmf-godunov` is reachable only from `model_column`, which is smooth and essentially incompressible; the flare/shock scenarios run the two-fluid solver, which this mode does not touch. The regime where the `gamma` mismatch would matter most — a strong shock inside the partial-ionization zone — has been tested only in the ideal-gas shock tube (where `gamma` is correct by construction), not in the mixture.
4. **The energy offset is upwinded by the sign of `u*` alone.** That is SWMF's construction and it is exact for a contact, but it is a piecewise-constant switch: at a face where `u*` crosses zero the offset flips side discontinuously. No pathology was observed here because `u* ~ 0` also makes the energy flux `~ 0`, but it is not a smooth function of the state.
5. **`Q_V` and `eps_p` are marginally worse, `Q_eff` and the dissipation marginally better.** The two solvers trade small amounts in opposite directions; nothing in these runs identifies a winner. **This remains true after the cutover** — the release choice rests on the thermodynamic argument of section 8.1 and on the SWMF methodology, not on a measured numerical advantage.

## 6. Assessment *(as written at the end of the controlled experiment — superseded by section 8)*

> The conclusion below is preserved verbatim because it is the honest state of the *measured* evidence, and section 8 does not overturn it: **no measured quantity on this model separates the two solvers.** What changed is the basis of the decision. The promotion rests on the thermodynamic argument of section 8.1 — frozen composition is the physically appropriate limit for chromospheric hydrogen — and on the SWMF methodology, not on a measured numerical advantage. The two specific objections raised below have since been addressed: the fixed `gamma` is reinterpreted in 8.1 as the intended closure rather than an error, and the timestep inconsistency is fixed in 8.2.

The method is correctly implemented and demonstrably correct on the standard tests, runs the production model with **zero fallbacks** over 28 million face solves, and produces results within **0.35 %** of the release Roe solver on every physical observable measured, for **+13 % explicit-RHS cost** (+4.7 % of measured step work).

That is enough to say the SWMF approach **works** on this model. It is not enough to call it better: no measured quantity separates the two beyond the noise of this comparison, and the 24 % characteristic-speed error from the fixed `gamma` is a real structural disadvantage that this smooth, low-Mach test cannot exercise.

Continuing the investigation is justified, and the informative next step is a case where the two solvers should *disagree* — a genuinely compressive, shock-forming event in the partial-ionization zone. Considering it as a release candidate is not justified on this evidence, and would in any case require first deciding what to do about the fixed `gamma`.

## 7. Reproduce

```bash
# conduction-off hydrostatic leg
for FLUX in roe-local swmf-godunov; do
  ISO_RIEMANN=$FLUX ISO_HEAT_FLUX=0 \
  CHROMO_T_END=20 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=5 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=0 CHROMO_FACE_FLUX_STRIDE=40 \
  scripts/run_chromo_realtime.sh outputs/model_column/swmf_godunov/hse_$FLUX.txt \
      full no-ionization model_column - 20.0 no-cooling
done

# production-like leg
for FLUX in roe-local swmf-godunov; do
  ISO_RIEMANN=$FLUX \
  CHROMO_T_END=100 CHROMO_FRAME_DT=10 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_OUTER_COND_DIAG=1 \
  CHROMO_FACE_FLUX_STRIDE=1780 CHROMO_OUTER_COND_STRIDE=89 CHROMO_EOS_COUNTS=1 \
  scripts/run_chromo_realtime.sh outputs/model_column/swmf_godunov/prod_$FLUX.txt \
      full no-ionization model_column - 20.0 no-cooling
done

.venv/bin/python util/pressure_reconstruction_diag.py \
  --run roe=outputs/model_column/swmf_godunov/prod_roe-local.txt \
  --run godunov=outputs/model_column/swmf_godunov/prod_swmf-godunov.txt --targets 100
```

---

## 8. Release cutover

This section supersedes sections 4d, 5.1 and 5.2 and the assessment in section 6.

### 8.1 Why the fixed `gamma = 5/3` is the reason for the promotion, not an obstacle to it

Section 5.1 treated the fixed `gamma` as a 24 % error in the characteristic structure. That framing implicitly assumed the equilibrium index `Gamma1` is the correct acoustic response. It is one of two correct responses, and which one applies is a physical question about relaxation, not a numerical one.

The equilibrium Saha internal energy of this mixture splits exactly:

```
e_int = p/(5/3 - 1) + e_ion
```

a monatomic **translational** part and an **ionization reservoir**. The two limits of the acoustic response are:

| Closure inside the wave | Acoustic index | Value in this column |
| --- | --- | --- |
| Saha equilibrium re-established instantaneously | `Gamma1 = (d ln p/d ln rho)_s` | median **1.090**, below 1.5 over 94.9 % of cells |
| composition frozen (`e_ion` inert) | `gamma = 5/3` | 5/3 everywhere by construction |

Chromospheric hydrogen does not re-equilibrate on a wave timescale. Carlsson & Stein (2002, ApJ 572, 626) compute the hydrogen ionization/recombination relaxation time directly and find it dominated by the slow collisional leakage from the ground state to `n = 2`: about 1 s in the photosphere, rising to `~10^5 s` in the mid-chromosphere and falling to `~10^2 s` at the transition-region base, and still `10`–`10^3 s` inside shocks — "the timescale is still too long for the ionization to reach its equilibrium value in the shock". Leenaarts et al. (2007, A&A 473, 625) reach the same conclusion in 2D radiation-MHD and state the energetic consequence directly: "in non-equilibrium, hydrogen ionization is a less effective internal energy buffer", giving shock temperatures of 10,000–13,000 K against `<= 8000 K` in LTE. Leenaarts (2020, Living Rev. Solar Phys. 17, 3) puts the chromospheric hydrodynamic timescale at about one minute and concludes that with a long relaxation time "an increase in internal energy leads directly to a temperature increase" and "the gas in the chromosphere behaves somewhat like an ideal gas".

So the release interpretation is: **during the local Riemann interaction the composition is frozen, the ionization energy rides passively in `E0`, and the translational gas is monatomic with `gamma = 5/3`.** The exact Godunov solver then solves precisely that frozen-composition local problem. `Gamma1` is not discarded — it remains the closure's equilibrium derivative, is still reported in `.gamma_diag`, still sets the outer-boundary Mach cap, and is still what the Roe reference flux linearizes.

**What this argument does not establish.** It does not show that the release *result* is more accurate; no measured observable on this model separates the two closures (section 4b, and 8.4 below). It also does not resolve the separate question of per-step Saha re-equilibration — see 8.6.

### 8.2 CFL audit: the timestep now follows the flux, and the old "effective CFL 0.62" was wrong

During the controlled experiment `mixture_timestep` deliberately kept the equilibrium speed so both legs took identical steps. That is not acceptable for a release: the release flux propagates the frozen waves, so the requested `CHROMO_CFL` must be the CFL of *those* waves.

`mixture_timestep` now uses

```
dt = CFL * min_i  ds_i / ( |u_i| + sqrt( gamma_signal,i * p_i / rho_i ) )
gamma_signal = max(Gamma1, 5/3)   when swmf_godunov_flux    (release)
gamma_signal = Gamma1             otherwise                 (Roe / Rusanov reference)
```

`max` rather than a bare constant keeps the bound conservative if the table ever returns an index above the monatomic value, and it keeps the Roe comparison mode self-consistent rather than merely unchanged.

**Measured effect on the canonical release configuration: none, and the reason matters.** The CFL-limiting cell of the N=500/R4 column is the topmost cell (2153 km, ~22 kK, hydrogen fully ionized, `ds = 270 m`), where the tabulated `Gamma1` is already exactly `5/3`. The 1.24x speed gap lives in the partial-ionization interior, which is far from the CFL limit because its sound speed is much lower. Measured on the `t = 0` and `t = 100 s` states of the release column:

| Quantity | Value |
| --- | --- |
| CFL-limiting cell, equilibrium rule | cell 660, 2153 km, `Gamma1 = 1.66667`, `x = 1` |
| CFL-limiting cell, frozen rule | the same cell |
| `dt` ratio, equilibrium / frozen | `0.999999` |
| **max effective CFL of the OLD equilibrium-sized step, measured against the frozen wave speeds** | **0.5000** |
| step count, 100 s, either rule | 17781 |
| final physical time, either rule | 100.001 s |

**Section 5.2's claim that this flux ran at an effective CFL of about 0.62 at `CHROMO_CFL=0.50` is therefore withdrawn.** It assumed the limiting cell carried `Gamma1 ~ 1.09`; it does not. The corresponding warning that raising `CHROMO_CFL` toward 0.8 would give an effective 0.99 was wrong for the same reason.

The change is still required. It removes a hidden dependency rather than a measured error: on a different mesh, a deeper domain, a different outer temperature, or under refinement that puts a partial-ionization cell at the CFL limit, the frozen rule would give a step up to 1.24x shorter than the equilibrium rule, and only the frozen rule is safe there.

Because `dt` and the step count are unchanged on this configuration, the CFL 0.50 acceptance sweep of `docs/studies/conduction/coarse_model_column_physical_conduction_recap.md` — run under the Roe flux and the equilibrium signal speed — carries over to the release flux for **this configuration only**. It does not carry over to a materially different model shape.

### 8.3 What changed in the code

| File | Change |
| --- | --- |
| `scenarios/model_column.cpp` | `swmf_godunov_flux = gamma_mode`, `roe_characteristic_flux = false`; `ISO_RIEMANN=roe-local` becomes a reference override; the override log line names `swmf-godunov` as the release default |
| `src/single_fluid/mixture.cpp` | `mixture_timestep` signal speed follows the flux in force; the flux's block comment states the frozen-composition rationale and drops the "correspondingly more dissipative" claim, which the measurements contradict |
| `chromosphere.hpp` | `swmf_godunov_flux` documented as the release flux and as the thing that also sets the timestep; `roe_characteristic_flux` documented as the equilibrium-linearization reference, with the Roe-property claim explicitly not made |
| `src/single_fluid/exact_rs.hpp`, `mixture.hpp`, `integrator.cpp`, `scenarios/model_column.hpp`, `chromo_main.cpp` | "experimental" removed; update-path strings become `MUSCL-Hancock/Godunov` |
| `scripts/release_validation.sh` | equivalence leg pins `ISO_RIEMANN=swmf-godunov`; smoke artefact renamed `release_godunov_lnp_N500_100s`; **new hard gate: the run fails if `godunov.fallbacks` is nonzero** |
| `tests/chromo_tests.cpp` | release-defaults test inverted; `test_roe_local_face_flux_and_projection` now clears the Godunov flag so its Roe-vs-Rusanov comparison is meaningful again; **new `test_release_timestep_uses_frozen_signal_speed`** |

Nothing in the reconstruction, the Hancock predictor, the source treatment, the boundary conditions, the EOS or the conduction operator was touched.

### 8.4 Acceptance runs

All legs: N500/R4, 661 cells, CFL 0.50, 12 OpenMP threads, `termination=end_time`. Outputs under `outputs/model_column/godunov_release/`.

**Roe reference preserved bitwise.** The new binary with `ISO_RIEMANN=roe-local` reproduces the *old* binary's default (Roe release) run snapshot-for-snapshot and byte-for-byte in `.gamma_diag`, at `T_END=20` and again over the full 100 s leg. Demoting Roe changed nothing about Roe.

**Conduction-off hydrostatic leg** (`ISO_HEAT_FLUX=0`, `T_END=20`, 3620 steps both):

| | Roe (reference) | Godunov (release) |
| --- | --- | --- |
| `max\|F^rho\|` at `t = 0` | 1.3345e-12 | **1.0758e-12** |
| `max\|F^rho\|` at `t = 20 s` | 1.1713e-11 | **1.0503e-11** |
| `max\|V\|` at 20 s | 0.03022 m/s @ 2090 km | 0.03194 m/s @ 2094 km |
| `time.rhs_s` | 0.953 | 1.101 (+15.4 %) |
| fallbacks | — | 0 of 4,785,640 |

Both sit at the float32 round-off floor of this scheme; the differences are not resolved by this measurement, and the marginally lower Godunov defect must not be read as an improvement. *(float32 floor — see the precision note in §4a.)*

**Production-like 100 s leg** (17781 steps both):

| | Roe (reference) | Godunov (release) | Δ |
| --- | --- | --- | --- |
| mean V [m/s] | 14.956 | 14.959 | +0.02 % |
| min / max V [m/s] | 11.370 / 32.718 | 11.370 / 32.728 | — |
| `n(V < 0)` | 0 | 0 | — |
| `F_eff` [kg m⁻² s⁻¹] | 1.1165e-9 | 1.1168e-9 | +0.03 % |
| `Q_V` | 1.29e-3 | 1.30e-3 | +0.8 % |
| `Q_M` | 6.7e-4 | 6.6e-4 | −1.5 % |
| `Q_eff` | 9.61e-4 | 9.04e-4 | −5.9 % |
| `eps_p` max | 0.0002 % | 0.0002 % | 0 |
| `T_top` [K] | 21837.4 | 21837.4 | 0 |
| `p_top` [Pa] | 0.010288 | 0.010288 | 0 |
| BL thickness [km] | 12.140 | 12.140 | 0 |
| `q_phys` [W m⁻²] | 0.7783 | 0.7784 | ~0 |
| column mass drift | −2.858e-4 | −2.955e-4 | +3.4 % |
| `time.rhs_s` | 5.147 | 5.809 | **+12.9 %** |
| `time.conduction_s` | 7.757 | 7.729 | −0.4 % |
| fallbacks | — | **0 of 23,506,482** | — |
| `godunov.iterations_max` / mean | — | **1 / 1** | — |

**`scripts/release_validation.sh`: PASS.** Default-path equivalence holds against an explicit `ISO_RIEMANN=swmf-godunov ISO_RECONSTRUCTION=lnrho-v-lnp` run (snapshots identical, `.gamma_diag` byte-identical); the smoke sidecar header reads `roe=0 swmf_godunov=1`; the zero-fallback gate passes.

**Suite:** 11643/11643 checks, 0 failures; `ctest` 9/9.

### 8.5 A run-to-run note worth recording

The 100 s Godunov release run is *not* bit-identical to the 100 s Godunov run of section 4b, even though the flux code is unchanged. The cause is isolated: running the pre-cutover binary with `ISO_RIEMANN=swmf-godunov` reproduces section 4b exactly, so the only difference is the CFL rule, which perturbs `dt` at the `~10^-6` relative level (the top-cell `Gamma1` is `5/3` only to table precision). Both runs are individually reproducible, both take 17781 steps, and the resulting spread in `max V` (40.587 vs 40.649 m/s) sits inside the solver's documented float32 round-off sensitivity — the release recap records a 40.97 → 40.51 m/s spread from representation round-off alone. This is a sensitivity of the model at float32, not a defect of either flux.

### 8.6 The unresolved modelling question this cutover exposes

The release now freezes the composition inside the Riemann problem, but the EOS still **decodes and reconstructs every state through instantaneous Saha equilibrium**, so each timestep is effectively

```
Saha-equilibrated state  ->  frozen-composition hydro step  ->  Saha-equilibrated state
```

These are separate modelling choices, and the first does not imply the second. Read charitably, the algorithm is an operator-split *instantaneous-relaxation* approximation: hydrodynamics at frozen composition, followed by an infinitely fast relaxation of the ionization state back to equilibrium. That split is only justified when the relaxation time is short compared with the timestep — which is the opposite of the regime the frozen Riemann solve is justified by. The literature above puts the chromospheric hydrogen relaxation time at `10^2`–`10^5 s`, against a release timestep of `~5.5e-3 s` and a dynamical time of order a minute.

The comparison that matters is **physical relaxation time versus physical dynamical time**, not relaxation time versus timestep; a short numerical timestep is not by itself an argument against Saha. But on that correct comparison the chromospheric hydrogen case still fails the instantaneous-equilibrium assumption, which means the equilibrium decode is a known and unquantified approximation of the released model — arguably the largest remaining physical-consistency issue in it.

Resolving it is a modelling project, not a bug fix: it means evolving the ionization fraction (or the ion/neutral populations) with finite-rate ionization/recombination instead of enforcing Saha every step. The two-fluid research solver already has such a network. Nothing in this cutover attempts it. Discussed at length in `main.tex`.

### 8.7 Reproduce the cutover evidence

```bash
# Roe-reference bitwise preservation (needs a pre-cutover binary at $OLD)
CHROMO_BINARY=$OLD CHROMO_T_END=20 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=5 \
  scripts/run_chromo_realtime.sh outputs/model_column/godunov_release/oldbin_default_roe.txt \
      full no-ionization model_column - 20.0 no-cooling
ISO_RIEMANN=roe-local CHROMO_T_END=20 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=5 \
  scripts/run_chromo_realtime.sh outputs/model_column/godunov_release/newbin_roe_override.txt \
      full no-ionization model_column - 20.0 no-cooling

# acceptance legs: TAG=godunov is the default (no ISO_RIEMANN), TAG=roe adds ISO_RIEMANN=roe-local
# (see the loops in section 7; only the output directory differs)

scripts/release_validation.sh
```
