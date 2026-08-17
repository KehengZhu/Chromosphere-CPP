# Release cutover: the conserved state is double precision

**Status: current. This changed a release default.** 2026-08-17. The `model_column` release stores the conserved state — and every other stored physical scalar — in `double`. `float` storage is retained only as a default-off diagnostic build. `CHROMO_STATE_FLOAT64` is retired.

Companion documents: `float32_precision_control_experiment.md` is the controlled diagnostic that established the mechanism; this document is the production cutover, the SWMF audit that justified making it the release rather than a rescue of float32, and the post-cutover validation.

## 1. What changed

**Precision is now a single knob.** `chromosphere::Real` is the solver's scalar storage type and `Vec` is `arma::Col<Real>`:

```cpp
#ifdef CHROMO_STATE_FLOAT32
typedef float Real;   // DIAGNOSTIC ONLY
#else
typedef double Real;
#endif
typedef arma::Col<Real> Vec;
```

Every stored physical scalar is declared in terms of them: the conserved state, the mesh metric arrays, the `Grid` physical constants (`m_i`, `m_n`, `m_e`, `k_b`, `q_e`, `chi_H_J`, `g`, `mu_0`), `CFL`, `limiter_beta`, the run controls and `sim_time`. Nothing in the release path names `float` or `double` directly.

| | before | after |
|---|---|---|
| release storage | `float` | **`double`** |
| CMake option | `CHROMO_STATE_FLOAT64=ON` (diagnostic, promoted `Vec` only) | `CHROMO_STATE_FLOAT32=ON` (diagnostic, reverts `Real`) |
| retired name | — | `CHROMO_STATE_FLOAT64` — CMake hard-errors on it |
| diagnostic build tree | `build_omp_f64/` | `build_omp_f32/` |

### 1a. Why a `Vec`-only promotion was not enough

The earlier `CHROMO_STATE_FLOAT64` experiment promoted `arma::Col<float>` to `arma::Col<double>` and nothing else. Six sites in the release path still narrowed through `float`:

| site | what it narrowed |
|---|---|
| `src/single_fluid/mixture.cpp:186-188` (`mixture_pack_ghost`) | every packed state and ghost state, including the IC |
| `src/single_fluid/mixture.cpp:799` | the MUSCL-Hancock predictor row |
| `src/single_fluid/mixture.cpp:871` | **the RHS row — the increment at the centre of this whole problem** |
| `src/single_fluid/mixture.cpp:920` | the per-cell timestep |
| `src/single_fluid/integrator.cpp:158-162` | the packed-energy ULP round into the caloric domain |
| `src/grid.cpp:27-34`, `chromosphere.hpp:527-534` | `m_i`, `m_n`, `m_e`, `k_b`, `q_e`, `chi_H_J` — the constants entering pressure and temperature everywhere |

One narrowing cast anywhere in the update chain reimposes the coarser representation on the whole chain regardless of what `Vec` holds. This is measurable, not theoretical: with `Vec` double but the pack still casting through `float`, a state packed at an EOS table endpoint arrived `4.0e-8` relative below the solver's bracket — exactly one float32 ULP — which is what made `test_release_production_table_domain_corners` unrunnable under the old option and forced it to be skipped.

**Lesson to keep: a precision knob has to be coherent or it does nothing.** That is also how SWMF does it — one `PRECISION` variable in `Makefile.conf`, one default `real` kind, no per-array exceptions.

## 2. Why double, and not a rescue of float32

Two independent reasons. Either alone would be sufficient; together they close the question.

### 2a. Measured necessity

From `float32_precision_control_experiment.md`: in the 4000 s N=500 run the per-step continuity increment `dt*(-dF/ds)` is about `2.6e-17 kg m^-3 s^-1` while one float32 ULP of `rho ~ 9e-10` is about `5.6e-17`. **619 of 660 cells sat below half a ULP** (median ratio 0.30), so `rho += dt*rhs` rounded straight back to `rho` every step. There is no compensated summation, so this is permanent stagnation, not slow drift.

The controlling dimensionless number is the ratio of the local evolution timescale `tau = |U|/|dU/dt|` to the timestep, against the reciprocal unit roundoff `1/u`. Stagnation when `tau/dt > 1/u`:

| | `dt` | `tau` | `tau/dt` | `1/u` | margin |
|---|---|---|---|---|---|
| this column, float32 | 5.6e-3 s | 1.9e5 s | 3.5e7 | 1.7e7 | **fails by 2x** |
| this column, double | 5.6e-3 s | 1.9e5 s | 3.5e7 | 9.0e15 | safe by 2.6e8 |

**Double is not a marginal fix — it is eight orders of magnitude of headroom.** That matters for the decision: it means the choice is robust to a factor of 10 in resolution, in run duration, or in how slow the physics we care about turns out to be, none of which we want to re-litigate later.

**What manufactures the problem is resolving the chromosphere.** A fine mesh and a ~7 km/s sound speed give a tiny `dt`, while mass loading a dense stratified layer gives a long `tau`. AWSoM's coronal grid never has both at once (§3f), which is why nobody upstream has hit this.

### 2b. Consistency with the coupling target

A standalone audit of the local SWMF tree (`/Users/zkeheng/SWMFSoftware/SWMF`, `gitinfo.txt`: SWMF master 2026-03-16 `3313194d`, GM/BATSRUS master 2026-03-18 `1573045`) established the following. All of it is read from source, not assumed.

**Persistent state precision.** BATS-R-US declares the state with no kind parameter:

```fortran
! GM/BATSRUS/src/ModAdvance.f90:127
real, allocatable, target :: State_VGB(:,:,:,:,:)
! :131  StateOld_VGB, StateOld_VG
! :151  Source_VC      (holds dt*(S - divF))
! :204  Flux_VXI, Flux_VYI, Flux_VZI
! :140  DtMax_CB
```

so precision is entirely a compile-time decision, and that decision is double in **every** shipped configuration:

```make
# share/build/Makefile.Darwin.gfortran:19-21   (and all 22 other templates)
DOUBLEPREC = -frecord-marker=4 -fdefault-real-8 -fdefault-double-8
PRECISION  = ${DOUBLEPREC}
```

Verified individually across all 23 `share/build/Makefile.{Darwin,Linux}.*` templates — gfortran, ifort, nagfor, nvfortran, pgf90, crayftn, xlf90, flang, absoft, lf95 — **including the nvfortran/OpenACC GPU builds.** `share/Scripts/Config.pl:124` sets `my $DefaultPrecision = 'double'`. A `-single` switch exists (`Config.pl:168`) but is legacy: `grep -n single Makefile.test` over the 109 KB SWMF regression suite returns one hit, and it is the phrase "single-ion" in a test description. The single-precision path is not validated.

**Update precision.** The conservative update is two statements, both at the same default real kind:

```fortran
! ModUpdateState.f90:372-380
DtLocal = DtFactor*DtMax_CB(i,j,k,iBlock)
Source_VC(iVar,i,j,k) = DtLocal*(Source_VC(iVar,i,j,k) + &
   (Flux_VXI(iVar,i,j,k,iGang) - Flux_VXI(iVar,i+1,j,k,iGang) + ...) &
   /CellVolume_GB(i,j,k,iBlock))
! ModUpdateState.f90:728-731
State_VGB(:,i,j,k,iBlock) = StateOld_VGB(:,i,j,k,iBlock) + Source_VC(1:nVar,i,j,k)
```

The semi-implicit conduction stage — the direct analogue of our backward-Euler stage — is the same plain accumulation (`ModHeatConduction.f90:1743-1762`).

**No special mechanism exists.** `grep -rni "compensated|Kahan|neumaier"` over `GM/BATSRUS/src/` and `srcBATL/` returns **zero hits**. `grep -rn "Real8_|real\*8|double precision|selected_real_kind"` over `GM/BATSRUS/src/*.f90` returns 21 hits and **every one is a time or an I/O buffer** — `ModMain.f90:276,283` (`StartTime`/`EndTime`, absolute epoch seconds), CPU clocks, satellite datetimes, input time series, and the `State8_VC`/`State4_VC` restart records. Zero on state, flux or source. There is no well-balanced or hydrostatic-deviation formulation in the cell update (`grep -rni "well.balanced|hydrostatic"` over `GM/BATSRUS/src/` hits only the threaded-field-line modules); the only deviation formulation is B0/B1 splitting for the magnetic field (`ModB0.f90:80`), a different problem.

The one precision *philosophy* signal in the tree is the contrast at `ModMain.f90:50` `real :: Dt` and `:270` `real :: tSimulation` against `:276` `real(Real8_) :: StartTime`: SWMF selectively promotes the quantity with the large additive offset, and leaves the state to the global switch.

**Conclusion of the audit.** SWMF's answer to sub-ULP secular accumulation is not a mechanism; it is a policy — the conservative state has enough precision, everywhere, by default. **Our float32 state was the SWMF-inconsistent component.** This cutover is therefore not the import of an unusual high-precision requirement into SWMF; it is the correction of a standalone prototype to the precision SWMF production already uses. Coupling is strictly simpler afterwards: `share/Library/src/ModMpiModified.f90:16-19` aliases `MPI_REAL` to `MPI_DOUBLE_PRECISION` whenever `nByteReal == 8`, so a double Chromosphere state matches the framework's exchange type directly.

### 2c. Interpretation, against the four options posed

The audit was framed as a choice between: (A) SWMF has a mechanism we should adopt; (B) SWMF is mostly float32 but its regime does not expose the issue; (C) a dynamic chromosphere module genuinely needs higher-precision persistent conservative state; (D) something else.

**The answer is C.**

- **A is refuted.** There is nothing to adopt — no compensated summation, no mixed-precision accumulator, no deviation formulation for `rho`/`p`.
- **B's premise is false.** SWMF production is *not* float32; every shipped build is double and `-single` is unvalidated. B's *regime* half is true and worth recording (§3f), but it explains why the issue was never encountered upstream, not what precision SWMF uses.
- **C is correct, and it is the same precision SWMF production already uses.** That is the sentence to carry into any writeup or advisor conversation.

### 2d. Selective mixed precision was considered and rejected

A double persistent state with `float` scratch and geometry was the earlier preference. It is not the design:

- the arrays are small at N ~ 500–2000, so the memory saved does not justify architectural complexity;
- SWMF is uniformly double, so a split knob is the *less* consistent choice for coupling;
- a split knob is precisely the failure mode of §1a — the previous experiment was an accidental mixed-precision build, and it silently did nothing.

Revisit only if profiling over many field-line threads demonstrates a real cost, and then as a performance optimization with its own evidence, not as a way to preserve float32.

## 3. Post-cutover validation

### 3a. Test suite — nothing loosened

| build | result |
|---|---|
| `build_omp` (double, RELEASE) | **11647 passed / 0 failed** |
| `build_omp_f32` (float, DIAGNOSTIC) | **11647 passed / 0 failed** |

**No tolerance was loosened anywhere in the suite.** The 29 failures previously documented under `CHROMO_STATE_FLOAT64` do not occur once the knob is coherent — they were a symptom of the incomplete `Vec`-only promotion (§1a), not of double precision exposing loosened assumptions. The exact-zero assertions that were suspected of being float32 quantization artifacts turn out to be structurally exact — a `v = 0` inner boundary row, heating windows evaluated outside their support, and an energy difference between two face states built with identical `phi` — and they hold in both precisions.

`test_release_production_table_domain_corners` is no longer skipped in either precision. Its ULP walk now steps in the storage type and nudges against the **solver's own EOS acceptance** rather than a test-side recomputation of the endpoint energy. The two differ by O(1e-14) relative because the solver inverts through a log-domain `n_h` path, which is coarser than the 64-eps bracket tolerance in `temperature_from_rho_eint`: one float ULP (6e-8) hid that mismatch, one double ULP (1.1e-16) does not. Using the real predicate makes the construction precision independent, and is a stronger assertion than the original.

### 3b. Release validation

`scripts/release_validation.sh` passes in the double release build:

- default-path equivalence: snapshots and `.gamma_diag` **bit-identical** to an explicit `ISO_RIEMANN=swmf-godunov ISO_RECONSTRUCTION=lnrho-v-lnp` run;
- N500 (661 cells) 100 s smoke: `termination=end_time`, `godunov.faces=23507804`, **`godunov.fallbacks=0`** (`bad_input=0 vacuum=0 negative_p=0 no_converge=0 bad_sample=0`);
- `eps_p` max 0.0002 % / median 0.0000 %; `Q_V` 0.00133, `Q_M` 0.00073, `Q_eff` 7.26e-04;
- mean V 14.070 m/s, min 10.473, max 31.790, **`n(V<0) = 0`**;
- `F_eff` 1.0414e-09 kg m^-2 s^-1, `T_top` 21838.6 K, `p_top` 0.010288 Pa, boundary-layer thickness 11.870 km, `q_phys` 0.7728 W m^-2, mass drift -3.454e-04.

### 3c. Hydrostatic fixed point — and a result we did not have before

Re-measured in both precisions, N=500/R4, 20 s, conduction off, release Godunov flux. Full table and provenance in `reference_free_release_recap.md`; the two conclusions:

| N=500/R4, Godunov, 20 s | float32 | **double (release)** |
|---|---|---|
| t=0 max\|F^ρ\| [kg m⁻² s⁻¹] | 1.0758e-12 @ 1601 km (base) | **1.7782e-13 @ 2149 km (upper TR)** |
| max\|V\| at 20 s [m/s] | 3.1943e-02 @ 2094 km | **5.7274e-03 @ 2153 km (top cell)** |
| defect / round-off floor | **1.07** | **9.5e+7** |
| N=500 → N=1000, max\|F^ρ\| | 1.08e-12 → 2.15e-12, order **−1.0** | 1.78e-13 → **3.91e-14**, order **2.19** |

**(i) float32 was at its representation limit; double is nowhere near its own.** 1.07× the float floor versus ~10⁸× the double floor. The surviving double residual is entirely discretization and boundary closure, with no representation component.

**(ii) The scheme is second-order convergent on the hydrostatic fixed point, and that could not be measured before.** In double, refining N=500 → N=1000 reduces the defect at order 2.19 (max), 1.95 (base face) and 1.90 (20 s max\|V\|). In float32 the *same* refinement made the defect worse — apparent order −1.0, the signature of per-face ULP noise. So the standing statement "no clean truncation-error trend survives between N=500 and N=1000" was a property of float32 storage, not of the scheme, and it no longer holds. This retires a real ambiguity in the reference-free release argument.

**What survives is real, and is the open defect.** In double the base face is still dominated by the dissipative part of the flux (`f_central = -4.78e-15`, `f_diff = +1.2935e-13`, `f_total = 1.2457e-13`), so the remaining base defect is Godunov dissipation acting on a genuinely nonzero reconstructed jump — the same first-interior-face artifact of §4, not bit noise. And the double max\|V\| now sits on the topmost cell, where the value is nearly precision-independent (5.27e-3 float vs 5.73e-3 double): an outer-boundary-closure residual precision does not touch, while the interior 2090 km peak collapsed.

### 3d. 4000 s N=500 comparison — and a correction to the earlier experiment

Both legs from the same source tree, N=500/R4 (661 cells), 4000 s, CFL 0.50, 12 threads, `CHROMO_FRAME_DT=20`, 217 sidecar frames, `termination=end_time`, **`godunov.fallbacks=0`** in both.

**First, a regression check that makes the rest trustworthy.** The `CHROMO_STATE_FLOAT32` build reproduces the archived pre-cutover release run **byte-identically** — `.txt`, `.gamma_diag` and `.faceflux` all `diff`-clean against `outputs/model_column/lnp_godunov_N500_4000s.*` (712,314 steps, 941,679,108 exact face solves in both). So the `float` → `Real` sweep is exactly semantics-neutral, and every difference below is a genuine precision effect rather than an artifact of the refactor.

The double leg takes 712,029 steps and 941,302,338 face solves — 0.04 % fewer. **The two legs are no longer step-identical**, because a coherent knob promotes the `Grid` physical constants as well as the state, which perturbs `c_s` at the `1e-8` level and hence `dt`. The earlier experiment could claim identical tallies precisely because it left those constants in `float`.

| Quantity, N=500/R4, 4000 s | float32 (diagnostic) | **double (RELEASE)** | earlier "f64 control" |
|---|---|---|---|
| `f_total` base / top at 2000 s | 2.682 | **1.006** | 1.008 |
| `f_total` base / top at 4000 s | 2.811 | **1.003** | 0.996 |
| `f_total` spread over 1600–1850 km, 4000 s | 119.5 % | **13.0 %** | *2.0 %* |
| Cells with per-step mass increment below 0.5 ULP(`rho`) | 619 / 660 | **0 / 660** | 0 / 660 |
| Column mass budget, predicted / observed | **-214.9** | **0.994** | 0.94 |
| Base cell velocity | 1.392 m/s | 0.356 m/s | 0.506 m/s |
| Top cell velocity | 13.18 m/s | **9.50 m/s** | *13.58 m/s* |
| `f_total` at the top face | 3.753e-10 | **2.702e-10** | *3.869e-10* |
| `\|f_diff\|` / `f_total` at the first interior face | 16.30 % | 16.10 % | 15.9 % |
| Cell-0 `rho V` / `f_total[0]` | 1.280 | 1.275 | 1.275 |

**The gross continuity violation is gone.** Base/top goes from 2.811 to 1.003, and the column mass budget closes (predicted/observed 0.994, against -214.9 in float32 — wrong by two orders of magnitude *and* in sign). No cell is frozen.

**Two entries of the earlier experiment are corrected, and both matter.** The italicised column above is that experiment's "float64 control", which we now know was an incoherent hybrid — `Vec` double, but the pack, predictor, RHS, timestep and every physical constant still `float`. It is not a double run and **its numbers must not be quoted**:

- **The lower-chromosphere spread is 13.0 %, not 2.0 %.** The double column is flat *base-to-top* but still carries a broad ~13 % monotone decline within 1600–1850 km, minimum at the top of that window (1849.7 km), and it is not the first-face artifact — dropping face 0 changes the figure by 0.02 pp. It relaxes, but very slowly: 14.3 % at 500 s, 14.2 % at 1000 s, 13.3 % at 2000 s, 13.1 % at 3000 s, 13.0 % at 4000 s. So the column is still relaxing at 4000 s. That is consistent with, and independent evidence for, the standing long-duration caveats — it is not a precision problem.
- **Top-of-domain evaporation observables move by ~28 %, not ~3 %.** Top velocity 13.18 → 9.50 m/s and top `f_total` 3.753e-10 → 2.702e-10. The earlier experiment's conclusion that "the release evaporation numbers are therefore not badly wrong ... affected at the ~3 % level" was an artifact of its hybrid build and is **withdrawn**. The physical reading is coherent: with the lower chromosphere unable to drain, float32 sustained a mass supply the column does not actually have, so it over-predicted the evaporation flux. **Any previously quoted release evaporation number from a float32 run is high by of order 30 % and must be re-measured.**

The comparison figure is `visualization/model_column/precision_cutover_N500_4000s.png` (catalogued in `visualize_commands.md` §14).

### 3e. A pre-existing open issue that the cutover closed

`docs/studies/validation/long_run_output_and_cfl_validation.md` §7.4 recorded a serious standing finding: the residual velocity field **did not converge in time at all** — refining the timestep made the level-to-level difference *larger*, and every CFL from 0.125 to 0.50 sat at the same `L1 |dv| ~ 1.0-1.2 m/s` against the finest run. It was hypothesised to be float32 cancellation and left as an open production risk.

**It was, and it is fixed.** Re-measured on the release mesh (N=500/R4, 20 s, conduction on, release Godunov flux), with every run interpolated in time to exactly `t = 20.000 s` from a 0.02 s cadence — necessary, because runs overshoot `CHROMO_T_END` by up to one `dt` and at `|dv/dt| ~ 0.14` that misalignment alone injects ~7e-4 m/s, the size of the whole double signal:

| dt-halving pair | **double** L1 | order | **double** Linf | order | float32 L1 | order |
|---|---|---|---|---|---|---|
| 0.50 vs 0.25 | 5.340e-3 | — | 6.398e-2 | — | 1.476e-1 | — |
| 0.25 vs 0.125 | 1.409e-3 | **1.92** | 2.610e-2 | **1.29** | 3.727e-1 | **-1.34** |
| 0.125 vs 0.0625 | 4.537e-4 | **1.64** | 1.200e-2 | **1.12** | 8.284e-1 | **-1.15** |

Against the finest run, `L1 |dv|` in double is monotone `7.157e-3 -> 4.537e-4 m/s` across CFL 0.50 → 0.125, while float32 stays pinned at `1.03-1.20 m/s` at every CFL. The mechanism shows directly in the fields: over an 8x step-count range double is stable (mean `|v|` 2.4015 → 2.4084 m/s) and float32 *degrades* (2.418 → 3.271 m/s, max 40.71 → 44.39). **Error that grows with the number of additions is round-off accumulation, not truncation.** Repeated on §7.4's original 2638-cell mesh the conclusion holds and the float32 runaway is worse (mean `|v|` 2.698 → 6.304 m/s).

The float32 control reproducing the pathology at the identical configuration is what makes this a causal attribution rather than a coincidence.

**Bounded honestly.** §7.4's exact configuration is not reproducible on a release build (`ISO_NUMERICAL_DIFFUSIVITY_MULT` is rejected, the Roe flux is demoted), so this is the current release configuration on two meshes with matched controls, not a bit-level replay. The CFL acceptance matrix has **not** been re-run in double, so CFL 0.50 still stands on its original evidence and no positive noise floor has been established for m/s-scale chromospheric velocities — what is retired is the mechanism, not the caution.

The same pattern appears in the mesh-refinement ladder: at `ISO_NS` = 48…2000 the 5-step residual velocity in double runs `8.11e-1, 7.62e-2, 1.88e-2, 4.99e-3, 1.32e-3, 3.32e-4 m/s` — unbroken order ≈ 2 over a 41x refinement — while float32 tracks it only to `ISO_NS = 250`, degrades to order 0.79 at 1000 and *reverses* at 2000. Full tables: `reference_free_release_recap.md`.

### 3f. Would AWSoM hit this?

Estimated, and clearly labelled as an estimate. Production AWSoM does not resolve a dynamic chromosphere: `ModChromosphere.f90:21-22` and `#CHROMOBC` in `SWMFSOLAR/Param/PARAM.in.awsom` make it a boundary condition at `n = 2e17 m^-3`, `T = 5e4 K`, with the unresolved transition region *broadened* by `extension_factor` (`ModChromosphere.f90:86-89`) rather than resolved; AWSoM-R puts the chromosphere/TR on 1-D threads below `r = 1 R_sun` (`ModTransitionRegion.f90:147` `rChromo = 1.0`) where the density is **algebraically reconstructed** from a barometric/hydrostatic closure (`ModThreadedLC.f90:446,456-467`) and never time-integrated at all. The finest base radial cell is ~300 km (`GM/BATSRUS/Param/CORONA/grid_awsom.dat`, values are `ln r`; the first step is `4.309e-4`), about three orders of magnitude coarser than our refined chromospheric mesh, so `dt` is ~1e3 larger while coronal timescales are comparable or shorter. And the production solar-corona run is `#TIMEACCURATE F` — a local-time-stepping steady-state relaxation over 1000 + 70000 + 80000 iterations — in which `dt*R < u*U` is the *convergence criterion*, not a failure.

Estimated `tau/dt ~ 1e3`–`1e4` for AWSoM against thresholds of 1.7e7 (float) and 9.0e15 (double). **AWSoM would plausibly be unaffected even in float32.** That is why this failure mode has no precedent upstream, and it is exactly the difference a resolved dynamic chromosphere introduces.

## 4. What this does NOT fix

**The localized first-interior-face artifact is precision independent and remains open.** `|f_diff|/f_total` at face 0 is 16.3 % (float32) versus 15.9 % (double) and the cell-0 `rho V` overshoot is 1.280 versus 1.275. It is a ghost / reconstruction / boundary-adjacent scheme effect and must be diagnosed on its own. **Double must never be described as solving all lower-boundary mass-flux problems.**

**The resolution limitation is untouched.** N500 long-duration evaporation mass flux is still not established as < 10 % grid-converged. That is a statement about mesh resolution, not storage precision, and the two must not be merged.

**`sim_time` is now double.** It was `float`, accumulating `dt = 5.6e-3 s` toward 4000 s where one float ULP is `2.4e-4 s` — about 23 ULPs per step, so not yet broken, but with only ~1.5 decades of margin before it would be on a longer run. The promotion removes that as a future hazard rather than fixing an observed defect; recorded here so nobody re-derives it.

**One deliberate behaviour change outside the state.** `run_control.hpp` computes the legacy `time_mult`-derived step cap as `max(10000, 10000*time_mult)`, previously in float and annotated "bit-for-bit the historical expression". It is now `double`. The value is cast to `long long`, so the resulting cap is unchanged for every sane `time_mult`; a pathological many-digit `time_mult` could shift the integer cap by one. `CHROMO_T_END` / `CHROMO_STEP_CAP` bypass this path entirely.

**`run_control.hpp` deliberately does NOT use `Real`.** It is `double` unconditionally. That header is intentionally free of Armadillo and of the solver's precision knob so it can be unit tested without a simulation (its own file comment says so), and nothing in it is stored state — it turns one CLI scalar into an integer step cap and schedules snapshots. Making it name `Real` would have forced a dependency on `chromosphere.hpp`, or a second local typedef of the knob, which is the incoherence of §1a in miniature. Verified self-contained: `#include "run_control.hpp"` compiles on its own. The rule is "nothing in the release path may name `float` or `double` directly" *for stored physical scalars*; run control is not one.

## 5. Reproduce

```bash
# Release (double) — the default; no precision flag
cmake -S . -B build_omp -DCMAKE_BUILD_TYPE=Release -DCHROMO_ENABLE_OPENMP=ON
cmake --build build_omp -j 12
(cd build_omp && OMP_NUM_THREADS=12 ./chromo_tests)
scripts/release_validation.sh

# Diagnostic (float32) — reproduces the pre-cutover storage
cmake -S . -B build_omp_f32 -DCMAKE_BUILD_TYPE=Release \
      -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_STATE_FLOAT32=ON
cmake --build build_omp_f32 -j 12
(cd build_omp_f32 && OMP_NUM_THREADS=12 ./chromo_tests)
```

The 4000 s comparison legs and the figure are catalogued in `visualize_commands.md` §14.
