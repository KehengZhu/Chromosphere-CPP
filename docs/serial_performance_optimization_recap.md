# Serial performance optimization of the gamma-table solver

Pre-OpenMP single-thread optimization of the active `GAMMA_TABLE` solver path. Scope was limited to removing redundant work: the physics, the numerical scheme, the boundary conditions and the fail-closed behaviour are unchanged, and **every production output of the reference run is byte-for-byte identical to the baseline**. No parallel code was introduced.

## Result

| | baseline | final |
|---|---|---|
| wall time (20 s benchmark) | 124.72 s | **78.06 s** |
| ms/step | 4.383 | **2.742** |
| speedup | — | **1.60×** |
| steps | 28466 | 28466 (identical) |
| cells | 2638 | 2638 |

Region breakdown (`CHROMO_PROFILE=1`, seconds):

| region | baseline | final | Δ |
|---|---|---|---|
| rhs | 62.66 | 39.25 | −37% |
| conduction | 33.92 | 21.48 | −37% |
| decode | 22.43 | 17.22 | −23% |
| projection | 12.70 | 8.39 | −34% |
| boundary | 4.82 | 0.13 | −97% |
| cfl | 0.19 | 0.18 | — |

EOS inversion counters — all identical except the one the work targeted:

| counter | baseline | final |
|---|---|---|
| `inversion.calls` | 300 694 281 | 300 694 281 |
| `inversion.initial_guess_accepts` | 67 214 120 | 67 214 120 |
| `inversion.one_update_convergences` | 226 719 603 | 226 719 603 |
| `inversion.average_iterations` | 0.801972 | 0.801972 |
| `inversion.maximum_iterations` | 16 | 16 |
| `inversion.bisection_fallbacks` | 459 564 | 459 564 |
| **`inversion.bracket_evaluations`** | **466 960 322** | **13 521 116** (−97%) |
| `conduction.average_newton_updates` | 2.997681 | 2.997681 |
| `conduction.maximum_newton_updates` | 4 | 4 |
| `limiter.acoustic` | 28466 | 28466 |

The identity of `average_iterations`, `maximum_iterations` and `bisection_fallbacks` is the strongest evidence that the new EOS fast path reproduces the full safeguarded solver's decisions exactly rather than merely approximating them; only the bracket-endpoint evaluations it was designed to skip disappeared.

## Original baseline

```
2638 actual cells, 28 466 steps, 124.72 s wall, 4.38 ms/step
rhs 62.5 s | conduction 33.8 s | projection 12.6 s | decode 22.3 s
EOS inversions 300.7 M, average updates 0.80, one-update 226.7 M, bisection 0.15%
```

Reproduced exactly on this machine before any edit (124.72 s), so all deltas below are measured against a same-machine run, not the quoted history.

## Benchmark and verification commands

Timing run (no output files):

```bash
env GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat ISO_H_BASE=1600 ISO_DH=553 \
    ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=2000 ISO_HYDRO_T_DECOUPLE=1 \
    ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 ISO_REFINE_S_LO_KM=500 \
    ISO_REFINE_TRANSITION_KM=20 ISO_NUMERICAL_DIFFUSIVITY_MULT=1 \
    CHROMO_T_END=20 CHROMO_OUTPUT=0 CHROMO_GAMMA_DIAG=0 CHROMO_PROFILE=1 \
    ./build/chromo_main /tmp/chromo_serial_profile.txt \
    full no-ionization model_column - 20.0 no-cooling
```

Equivalence run — same case with every diagnostic sidecar enabled, so the conserved state, the gamma diagnostics, the face-flux capture and the outer-conduction capture are all compared:

```bash
env <same variables> CHROMO_T_END=20 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 \
    CHROMO_FACE_FLUX_DIAG=1 CHROMO_OUTER_COND_DIAG=1 \
    CHROMO_FRAME_STRIDE=250 CHROMO_FACE_FLUX_STRIDE=250 CHROMO_FACE_FLUX_TOP=40 \
    CHROMO_OUTER_COND_STRIDE=1 CHROMO_PROFILE=1 \
    ./build/chromo_main <dir>/out.txt full no-ionization model_column - 20.0 no-cooling
```

then `cmp` of `out.txt`, `out.txt.gamma_diag`, `out.txt.faceflux`, `out.txt.outercond` and a `diff` of the per-100-step `dt`/`time` log lines and of stderr (the outer-BC anchor line). `out.txt` is compared with its one header line excluded, because that line embeds the output path.

## Files changed

| file | change |
|---|---|
| `src/rhs.cpp` | dead first-pass face bundles removed; face-bundle zero-fills removed; 3-row reconstruction scratch; cached non-uniform reconstruction weights |
| `src/grid.cpp` | `broadcast()` split into a fingerprint check and `force_rebuild_metrics()`; direct packed fill; new conduction scratch rows |
| `chromosphere.hpp` | `Grid` geometry fingerprint + `force_rebuild_metrics()` + `metrics_generation()`; `GammaRhsScratch` reconstruction-weight cache; `GammaConductionScratch` fused-evaluation rows |
| `src/eos.cpp` | `CaloricEval` fused evaluation; algebraic Saha pressure inversion; zero/one-update inversion fast path; `InversionResult` carrying the accepted evaluation; fused face-state / projection / known-temperature packing |
| `eos.hpp` | `equilibrium_density_from_pressure`, `CaloricState`, `equilibrium_caloric_state`; optional caloric out-parameter on `decode_equilibrium_mixture`; optional caloric input on `pack_equilibrium_from_known_temperature` |
| `src/integrators.cpp` | conduction Newton pass on one fused evaluation per cell; decode evaluation carried into the first pass; converged thermodynamics reused by the known-temperature packing |
| `scenarios/model_column.cpp` | `gamma_density_from_pressure` uses the algebraic inverse |
| `tests/chromo_tests.cpp` | two new tests (geometry cache; algebraic-vs-bisection density inversion) |

## Optimizations

### 1. Dead RHS face work — 124.72 s → (with #2) 109.22 s

`rhs_explicit_mixture` built four one-sided face bundles in its predictor pass, but the predictor update differences only `fl_iph.flux - fr_imh.flux`. The `fr_iph` and `fl_imh` bundles — and the `wr_iph` / `wl_imh` reconstructions that feed them — were pure dead work; both are recomputed unconditionally by the corrector pass. Face-bundle construction went from **8 to 6 per timestep**.

`build_mixture_face` also zero-filled its three `n_state` outputs before overwriting every element. The loop writes `sub2ind(size(ns,num_of_eq), i, k) = i + k*ns` for all `i,k`, which enumerates `n_state` exactly, so `set_size` replaces `zeros` with no change in value.

Byte-identical.

### 2. Static mesh/geometry caching — combined with #1: **109.22 s** (rhs 62.66→51.32, boundary 4.82→0.47)

Everything `Grid::broadcast()` derives depends only on the static mesh (`ds_i`) and the magnetic geometry (`B_i`, `B_imh`, `B_iph`, `dinvB_ds_i`) — never on the ghost buffers. `model_column_update_bc` nevertheless called it twice per timestep, re-deriving `s_face`, `s_i`, `ds_iph_i`, `ds_imh_i`, `uniform_mesh` and seven packed `n_state` arrays. The packed arrays were built as `dst.zeros(n_state); for k: dst += scalar_to(src, k)`, i.e. **49 `n_state`-sized allocations and zero-fills per broadcast**.

Two changes:

* `broadcast()` now fingerprints the five source arrays and returns immediately when none changed; `force_rebuild_metrics()` does the unconditional rebuild. Static arrays are therefore rebuilt only on initialization, `resize`, mesh construction and magnetic-geometry edits, and an ordinary boundary refresh costs an `O(ns)` compare. This was preferred over simply deleting the per-step calls because it makes *every* scenario correct by construction rather than relying on each `update_bc` not touching geometry.
* The packed fill is direct instead of accumulated. To stay bit-exact the fill writes `0.0f + src[i]`: the old accumulation started from a zeroed buffer, and `0.0f + (-0.0f) == +0.0f`, which a bare assignment would not reproduce.

`resize()`/`init()` invalidate the fingerprint. `test_broadcast_static_metric_cache` asserts the packed caches equal the old accumulated construction **bitwise**, that a repeat `broadcast()` changes nothing and does not bump the generation counter, and that edits to `ds_i` / `B_iph` and a `resize` are all still picked up.

Byte-identical.

### 3. Algebraic boundary density inversion — no measurable wall-time change

`gamma_density_from_pressure` ran **100 geometric bisections (200 Saha evaluations)** for every pressure→density conversion. The pure-H closure inverts in closed form: with `S(T) = (2π m_e k_B T/h²)^{3/2} exp(−χ_H/k_B T)` and `y = p/(k_B T)`,

```
x²/(1−x) = S/n_H  and  p = (1+x) n_H k_B T   ⇒   x = √(S/(S+y)),  n_H = y/(1+x)
```

`equilibrium_density_from_pressure` (new, in `eos.cpp`) evaluates this in the log domain with two algebraically identical branches, each chosen so its exponential argument is non-positive — no overflow for either `y ≫ S` (deep neutral) or `y ≪ S` (fully ionized). It reuses the same `log_c` constant and term order as `saha_ionization_fraction_n_h`, and the same `eos_constants`. The EOS density-domain guard, the `domain_error` on non-positive/non-finite inputs and the `out_of_range` on an out-of-table pressure are all preserved verbatim in the caller.

`test_equilibrium_density_from_pressure_matches_bisection` sweeps 25 temperatures × 29 densities over the whole production domain (3.16 × 10³ – 10⁷ K, 10¹² – 10²⁶ m⁻³) and checks, per point, agreement with a local copy of the old bisection to `1e-12` relative, recovery of the density the pressure was built from, and pressure reconstruction to better than `1e-13` relative. Worst observed values are inside those bounds. The deep-neutral and fully-ionized limits and the three error paths are checked separately.

This is called only from the IC and the two outer ghost rungs (a handful of calls per step), so it is a robustness/cleanliness win rather than a throughput win — but the run stayed byte-identical, which also confirms the algebraic result rounds to the same `float32` ghost state as the 100-iteration bisection.

### 4. Fused caloric evaluation + one-step EOS fast path — 109.22 s → **89.06 s**

The Saha fraction (a `log`, an `exp` and a `sqrt`) dominates every thermodynamic query, yet the internal energy, pressure, heat capacity, ionization fraction and Gamma-state construction each re-derived it. `CaloricEval` computes all of them from **one** Saha call, each field written with exactly the expression and operation order of the function it replaces.

Fused call sites: `equilibrium_internal_energy`, `equilibrium_heat_capacity_impl`, `gamma_state`, `equilibrium_mixture_face_state` (was `gamma_state` + `equilibrium_internal_energy` = 2 Saha solves, now 1 — six times per timestep per cell), `pack_equilibrium_from_known_temperature` (2 → 1), and the Newton loop of `temperature_from_rho_eint` (energy + capacity were 2 → 1).

`temperature_from_rho_eint` gained a zero/one-update fast path:

1. evaluate energy **and** capacity once at the supplied guess;
2. accept if `|residual| ≤ 2e-13·scale`;
3. otherwise take one Newton step, applying **exactly** the safeguards the full solver applies (bracket narrowed to the guess on the residual's side, the same `4·eps` edge test, the same slow-progress test against the initial bracket width);
4. evaluate the candidate once and accept if converged;
5. otherwise fall through to the untouched safeguarded Newton/bisection solver.

This removes the two table-wide bracket-endpoint energies (~470 M Saha solves per 20 s run) from the ~97% of calls that converge in at most one update. A one-update convergence went from 6 Saha solves (guess + 2 bracket + iteration-0 energy + capacity + iteration-1 energy) to 2.

Equivalence argument for the fast path, beyond the measured byte-identity:

* The accepted candidate is `T_guess − residual/C_V` with `residual` and `C_V` computed by the same expressions at the same `T_guess`, so it is bit-identical to what iteration 0 of the full loop produces.
* The full loop's iteration 0 always fails its residual test whenever the guess test failed (they compare the same quantity against the same tolerance), so no full-path exit is skipped there.
* The fast path is additionally gated on the guess **and** the accepted candidate lying inside a `1e-9`-relative inset of the table temperature range. The bracket block below can return the endpoint temperature itself, or throw `out_of_range`, only when the target energy is within `64·eps` (≈1.4e-14 relative) of an endpoint energy — five orders inside that inset. So no call the full path would have terminated at an endpoint, or rejected as out of range, can be answered by the fast path instead.
* Any guard failure falls through to the unmodified solver, which re-runs its own guess check; behaviour for those calls is unchanged.

`invert_caloric` also returns the caloric evaluation at the accepted temperature (it just computed it), so `decode_equilibrium_mixture` and `project_equilibrium_single_fluid` build their `GammaState` with no further Saha solve. The one exit that has no evaluation at its return value — the bisection-midpoint exit — reports `has_state = false` and the caller falls back to a normal evaluation.

The `2e-13` caloric residual tolerance, the table bounds, the out-of-range errors, the non-convergence error and the safeguarded fallback are all unchanged.

Byte-identical.

### 5. Repeated EOS work between projection and conduction — 89.06 s → **81.85 s**, then 78.27 s with the packing reuse

**Conduction Newton loop (89.06 → 81.90 s, conduction 29.40 → 22.86 s).** Each residual pass computed the Saha fraction three times per cell: once for `n_e`/`n_HI`, once inside `equilibrium_heat_capacity`, and once inside `equilibrium_internal_energy` in the residual assembly — all at the *same* `(rho, T)`, since the temperature is not touched between those loops. They now share one `equilibrium_caloric_state` call, with the residual energy carried in a new `e_at_T` scratch row.

**Decode → first Newton pass (81.90 → 81.85 s, conduction 22.86 → 21.89 s).** `decode_equilibrium_mixture` gained an optional `CaloricState` out-parameter. The conduction stage's decode evaluates the caloric EOS at exactly `(th.rho, th.T)`, which is precisely where the first Newton residual pass starts, so that evaluation is carried over and iteration 0 does no Saha work.

**On the projection→conduction cache specifically.** Passing projection's *pre-rounding double* thermodynamic state into conduction was evaluated and rejected: the state conduction must describe is the packed `Vec<float>` conserved state, and the equivalence of the double state to a decode of its own `float32` projection is not provable across all seven rows. Validating such a cache against the packed state costs one caloric evaluation — exactly what the decode's fast path already costs — so it would buy nothing. What *is* provable, and is implemented, is the decode→first-pass carry above: it uses the same `(rho, T)` by construction, with the packed state as the sole source. The normal decode remains the only path into the stage, so there is no cache-provenance failure mode to fall back from.

The conduction Newton method itself is untouched: same quasi-Newton matrix, same lagged conductivity, same `2e-11` tolerance, same Thomas solve, same 40-iteration cap, same boundary treatment, same stagnation error.

Byte-identical.

### Optional items (attempted after the five main items were complete and measured)

**Three-row reconstruction scratch + cached non-uniform weights — 81.85 s → 78.21 s (rhs 42.74 → 38.96 s).** The gamma reconstruction limits exactly three variables (`log ρ`, `V`, `log T`) but its scratch was packed at the `num_of_eq` conserved-state stride, so every difference, ratio, limiter and extrapolation moved seven rows per cell and discarded four. All those operations are element-wise, so the stride is now 3 and every surviving element is unchanged. Separately, the non-uniform weights `W1/W3/W4` and the three slope-ratio metrics are functions of the static mesh alone; they are now built once per mesh generation (keyed on `Grid::metrics_generation()`) instead of rebuilt from six packed temporaries every timestep.

**Reuse of the converged conduction thermodynamics in the known-temperature packing — 78.21 s → 78.27 s (conduction 21.93 → 21.53 s).** On convergence the scratch rows still describe the final accepted temperature, and `pack_gamma_known_temperature` re-forms `rho` from the identical `RHO_I+RHO_N` rows, so `pack_equilibrium_from_known_temperature` accepts the evaluation instead of re-solving Saha. Real in the conduction region, inside noise in the total.

**Gamma1 uniform-axis index lookup — not done.** The stated precondition is that the axes are *verified* uniform in log space. The production table's axes are written as 16-digit decimals (501 × 57 points), so consecutive log-axis differences are not exactly equal in double precision and a computed index cannot be proven to match `std::upper_bound` at every knot. Replacing the search under those conditions risks an off-by-one at a node and a different interpolation, which the "no unexplained numerical change" requirement rules out.

## Tests

* `ctest`: 8/8 passing (`chromo_tests`, `eos_production_table`, `eos_production_table_checksum`, three `gamma_reject_*` guards, `gamma_output_sidecar`, `conductive_flux_python`).
* `chromo_tests`: **8897 / 8897** checks passing (baseline 8146; +751 from the two new tests).
* New: `test_broadcast_static_metric_cache` (19 checks) and `test_equilibrium_density_from_pressure_matches_bisection` (732 checks).
* No test tolerance was weakened.

`model_c7` and `model_flare` fail to start in non-gamma mode with `std::logic_error: Gamma1 table is empty`. This was verified to reproduce identically on a pristine `HEAD` worktree, so it is pre-existing and unrelated to this work; it is not investigated here.

## Numerical-equivalence evidence

Full 20 s reference run, baseline binary vs final binary:

| artifact | result |
|---|---|
| `out.txt` (all 2638 cells × 7 rows, 114 snapshots) | identical |
| `out.txt.gamma_diag` (ρ, v, T, x, n_e, n_HI, p, Γ₁, κ_phys, κ_solver) | identical |
| `out.txt.faceflux` (reconstruction, limiters, Rusanov split, eq-wb residual) | identical |
| `out.txt.outercond` (28466 rows: T_top, T_wall, κ_phys, χ_num, κ_num, q_phys, q_num, q_total) | identical |
| step count | 28466 = 28466 |
| `dt`/`time` sequence (per-100-step log) | identical |
| active CFL limiter | `acoustic` 28466 / 28466 in both |
| conduction Newton counts | avg 2.997681, max 4 in both |
| EOS inversion counters | identical except `bracket_evaluations` |
| outer-BC anchor (stderr) | identical |

`out.txt` conserves the same E_I+E_N remainder structure and the same ρ, so the conservation tests in `chromo_tests` (mass conservation, uniform fixed point, RK4 fixed point, projection conservation, Stage-7 total-energy sources) are exercising the same numbers as before.

Byte identity holds for **all** changes, including the algebraically-restructured EOS ones, so no maximum-absolute / maximum-relative difference needs reporting: it is exactly zero everywhere.

## Final profile

```
wall 78.06 s | 28466 steps | 2638 cells | 2.742 ms/step
rhs        39.25 s  (50.3%)
conduction 21.48 s  (27.5%)
decode     17.22 s  (22.1%)
projection  8.39 s  (10.8%)
boundary    0.13 s
cfl         0.18 s
```

(Region timers overlap: `decode` is entered from within `rhs`.)

## Remaining hotspots

macOS `sample` over 20 s of steady-state stepping, by top-of-stack:

| symbol | samples |
|---|---|
| `EosGammaTable::gamma1` | 2972 |
| `log` (libm) | 2782 |
| `require_finite_positive` | 2567 |
| `exp` (libm) | 1280 |
| `invert_caloric` | ~840 |
| `build_mixture_face` | ~430 |
| `equilibrium_mixture_face_state` | ~380 |
| `pow` (libm, `κ_e ∝ T^{5/2}`) | 302 |

The remaining cost is now almost entirely **transcendentals and the Gamma1 table lookup**, not redundant Newton iterations. Two specific observations for whoever picks this up next:

1. **Duplicate `log()` per thermodynamic query.** One `equilibrium_mixture_face_state` call currently computes `log(n_H)` three times (`require_n_h_in_bounds`, the Saha solve, `gamma1`) and `log(T)` twice (the Saha solve, `gamma1`) — five logs where two suffice. Hoisting them and threading precomputed `log_t` / `log_n` into log-aware overloads of `saha_ionization_fraction_n_h`, `require_n_h_in_bounds` and `gamma1` is a pure value reuse (`std::log` is deterministic), so it would be bit-identical. `log(n_H)` is additionally constant across the Newton iterations of a single inversion at fixed density and could be hoisted out of the loop. Given `log` is 2782 samples and `gamma1` a further 2972 (its two internal logs counted separately), this is the largest remaining serial win by a wide margin. It was **not** attempted here because it is outside the task's explicit optional list.
2. **`require_finite_positive` at 2567 samples** is disproportionate for a two-comparison guard; it is not being inlined out of the hot EOS paths. Worth checking whether marking it `inline` (or moving it to a header) recovers it.

## Recommended next steps for OpenMP

The serial work above also removed several structures that would have been hazards under threading:

* **Per-step `n_state` allocations are gone** from `broadcast()` and from the non-uniform reconstruction-weight construction. Both were hidden allocation storms that would have serialized on the allocator.
* **`GammaRhsScratch` and `GammaConductionScratch` are now the single mutable working set** of the hot path. They are `mutable` members of a `const Grid&`, so they are *shared* across any threads that enter `rhs_explicit_mixture` or `apply_gamma_conduction_stage` concurrently. Nested/concurrent RHS evaluation must be ruled out, or the scratch must become thread-local, before any parallel region is opened.
* `Grid::store_eos_temperature_hint(i, T)` writes a shared per-cell hint from a `const` method. Under a cell-parallel loop each thread writes only its own `i`, which is a benign but formally racy pattern; confirm the hint array is `ns`-sized and index-disjoint before relying on it.
* `Grid::broadcast()`'s fingerprint check is a read-only `O(ns)` compare and is safe to call from a serial region only; keep boundary updates outside parallel regions.

Where to parallelize, in expected-value order:

1. **`build_mixture_face`** — six independent per-cell loops per step, no cross-cell coupling, ~40% of `rhs`. Each iteration is a pure function of `(mixture[i], phi_face[i])`. This is the single best target.
2. **`decode_mixture_field_into`** — a per-cell loop over the inversion. Cell-independent apart from `store_eos_temperature_hint` (index-disjoint) and the ghost decodes, which are done after the loop.
3. **`project_equilibrium_single_fluid`** — same shape as the decode.
4. **The conduction stage's per-cell loops** (fused caloric evaluation, coefficient assembly, residual assembly). The **Thomas solve itself is inherently sequential** and must stay serial; only the loops around it parallelize. With ~4 Newton passes per step, the fork/join cost has to be amortized over the whole iteration, not per-loop.
5. The Armadillo element-wise expressions in the reconstruction are memory-bound and, at 3 rows × 2638 cells, likely too small to pay for a parallel region — measure before touching.

At 2638 cells the per-loop work is modest, so use a single parallel region per stage with `nowait` where possible rather than a `parallel for` per loop, and re-verify byte identity after each step: with a fixed schedule and no reductions, the cell loops above should remain bit-exact, which makes the same `cmp`-based verification harness reusable.

---

```
SERIAL_OPTIMIZATION: GO
```

* all tests passing (8/8 ctest, 8897/8897 checks, no tolerance weakened);
* no unexplained numerical change — every production output of the 20 s reference run is byte-for-byte identical;
* material measured per-step speedup: 4.383 → 2.742 ms/step (**1.60×**);
* no parallel code introduced.

---

# Final audit follow-up: Gamma1/log reuse and conduction packing

## Scope

This was the final narrowly scoped serial refinement before independent OpenMP work. It addressed only the three audit findings left above:

1. skip Gamma1 table queries in caloric-only decode, projection, and known-temperature packing;
2. reuse validated `log(T)` and `log(n_H)` coordinates in hot EOS paths;
3. remove the final conduction packing cache whose allocation/copy cost had not produced a repeatable wall-time benefit.

No physics, numerical parameter, mesh, boundary condition, CFL rule, Newton method, conduction method, diagnostic definition, or parallel code was changed.

## Files changed in this follow-up

* `eos.hpp`
* `src/eos.cpp`
* `src/rhs.cpp`
* `src/state.cpp`
* `src/integrators.cpp`
* `chromo_main.cpp`
* `tests/chromo_tests.cpp`
* `docs/serial_performance_optimization_recap.md`

The workspace already contained unrelated uncommitted work. This follow-up preserved it and did not reformat, revert, or absorb it.

## Gamma1 query removal

`MixtureThermo` still guarantees a valid `gamma1`. Callers that do not need sound speed now use explicit lighter result types instead of receiving a partially valid object:

* `CaloricMixtureThermo` for caloric-only decode;
* `ProjectedMixtureRows` for projection and known-temperature packing.

The full public wrappers share the same underlying caloric/row calculation and add exactly one Gamma1 query at the end. Production call sites were switched as follows:

* MUSCL predicted-state decode: caloric-only;
* equilibrium projection: row-only;
* conduction initial and ghost decode: caloric-only;
* conduction final pack: row-only.

CFL decode and all six MUSCL face bundles continue to query Gamma1. The focused RHS test verifies exactly `6*ns` queries, so the predictor contributes none. At `ns=2638`, the removed call sites represent `4*ns+6 = 10,558` avoided Gamma1 lookups per timestep, or 300,544,028 over the 28,466-step reference run.

Optional `CHROMO_EOS_COUNTS=1` instrumentation records Gamma1 and EOS-log operations for focused profiling. It is disabled by default and must remain disabled, or be made thread-local, during future parallel execution.

## Log-space coordinate reuse

A validated internal `EosCoordinates` value now carries physical and logarithmic coordinates together. The public APIs still accept physical `T` and `n_H`, preserve the same validation/exception behavior, and use exact `std::log`/`std::exp` rather than approximations.

The hot-path changes are:

* `invert_caloric` computes `n_H` and `log(n_H)` once outside the Newton loop;
* each candidate temperature computes `log(T)` once and reuses it for Saha and any Gamma1 lookup;
* Gamma1 density bounds and bilinear interpolation consume the same validated `log(n_H)`;
* MUSCL face construction consumes reconstructed `log(rho)` and `log(T)` directly, avoiding the previous `log(T) -> exp(T) -> log(T)` cycle inside Saha/Gamma1;
* Saha, caloric state construction, density guards, and Gamma1 interpolation share the same coordinates without approximate math or `-ffast-math`.

Direct physical-state and log-aware caloric/Gamma1 paths are exactly equal in the representative tests. The face API necessarily reconstructs physical values through `exp(log(x))`; exact equality was attempted first, and the observed difference was only a few ULPs in the isolated API comparison. The production reference evolution nevertheless remains byte-identical, including all saved face-flux diagnostics and timestep decisions.

## Final conduction packing-cache decision

The simple option was retained:

* removed the per-timestep `std::vector<CaloricState> converged_caloric(ns)` allocation and copy;
* removed the optional unchecked cached-caloric pointer from the packing API;
* retained the useful fused caloric evaluation and first-pass seed inside the conduction Newton solve;
* let the final row-only pack perform one authoritative fused caloric evaluation and retain the existing residual gate and conserved neutral-energy remainder.

There is therefore no stale or mismatched `(rho,T)` cache contract left to validate or fail closed.

## Focused tests added

`chromo_tests` increased from 8,897 to 9,037 checks. The new tests cover:

* caloric decode, row-only projection, and row-only packing perform zero Gamma1 queries;
* their full wrappers perform one Gamma1 query and produce the same numerical rows/thermodynamics;
* an actual explicit RHS evaluation performs exactly six face-state Gamma1 lookups per cell;
* log-aware Saha/caloric/Gamma1 equivalence across mostly neutral, partially ionized, mostly ionized, and lower/upper interior table states;
* face-state equivalence at a tight machine-precision bound after the unavoidable `exp(log(x))` round trip;
* inverted-temperature equivalence;
* one non-exact guess converges in exactly one Newton update with zero bracket-endpoint evaluations;
* one rejected fast candidate continues through the safeguarded bracket/fallback path and reaches the same root.

Final validation:

```text
cmake --build build -j4:             PASS
ctest --test-dir build:              8 / 8 PASS
./build/chromo_tests:                9037 / 9037 PASS
```

No production tolerance was weakened.

## Numerical-equivalence result

The diagnostic-enabled 20 s reference was rerun with the exact production configuration. Comparison against the pre-edit baseline showed:

| artifact / decision | result |
|---|---|
| conserved-state output | byte-identical |
| `.gamma_diag` | byte-identical |
| `.faceflux` | byte-identical |
| `.outercond` | byte-identical |
| step count | 28,466 = 28,466 |
| timestep sequence | identical |
| active CFL limiter | acoustic, 28,466 / 28,466 |
| EOS inversion decisions/counts | identical |
| conduction Newton counts | avg 2.997681, max 4 in both |

The main output header embeds its output filename. Comparing runs written under different names initially found that metadata-only difference on line 3; rerunning the final binary at the identical output path made the complete main file byte-identical. No physical or numerical difference remains unexplained.

## Benchmark before and after

Exact benchmark configuration: 20 s R4-local, `ISO_NS=2000`, 2,638 actual cells, outer refinement factor 4, output and diagnostics disabled, profiling enabled.

| run | wall time | ms/step | decode | RHS | projection | conduction |
|---|---:|---:|---:|---:|---:|---:|
| audit-follow-up baseline | 78.765 s | 2.767 | 17.359 s | 39.358 s | 8.420 s | 21.321 s |
| final 1 | 76.072 s | 2.672 | 8.833 s | 37.550 s | 7.703 s | 20.399 s |
| final 2 | 76.187 s | 2.676 | 8.823 s | 37.567 s | 7.823 s | 20.370 s |
| final 3 | 75.976 s | 2.669 | 8.825 s | 37.509 s | 7.675 s | 20.385 s |
| **final median** | **76.072 s** | **2.672** | — | — | — | — |

The median reduction is 2.693 s, or 3.42% (`1.035x`) relative to this follow-up baseline. Region timers overlap because decode is entered from RHS. All three final runs retained the same 28,466 steps, timestep/limiter sequence, EOS inversion counters, and conduction Newton counters.

## Remaining hotspots and OpenMP readiness

The dominant measured serial regions are now approximately:

```text
rhs         37.5 s
conduction  20.4 s
decode       8.83 s   # overlaps with RHS
projection   7.7 s
cfl          0.18 s
boundary     0.13 s
```

Further serial micro-optimization was not pursued because this task was restricted to the three audit findings and the retained implementation already produces a repeatable total wall-time reduction. The code is ready for independent audit before OpenMP implementation. The earlier scratch/hint/thread-safety cautions still apply; no OpenMP, thread, process, MPI, GPU, or SIMD code was added here.

---

```
FINAL_SERIAL_REFINEMENT: GO
```

* all tests pass;
* the complete production reference remains byte-identical;
* measured serial cost is lower by a repeatable 3.42%;
* no parallel code was added;
* no unchecked caloric-cache contract remains.
