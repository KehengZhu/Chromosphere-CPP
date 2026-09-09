# Validation and known limitations {#validation}

This page describes what is checked automatically, what was established by measurement, and what is explicitly *not* established. For the method being validated see @ref numerics; for the scenarios see @ref scenario_reference.

## Running the checks

Configure and build first. The test target is built by the default `all` target:

```bash
cmake -S . -B build
cmake --build build -j
```

**Full registered suite.** From the build directory:

```bash
cd build && ctest --output-on-failure
```

**C++ unit suite alone.** The test binary is standalone and takes no arguments:

```bash
./build/chromo_tests
```

It prints `[RUN  ] <name>` before each case and `[  OK  ] <name>  (<k> checks)` after it, then a summary block. As measured on the current tree:

```text
===== Summary =====
Passed: 11647
Failed: 0
```

That is **107 test cases and 11647 individual checks, 0 failures**. Note that `Passed:` counts individual `EXPECT_NEAR` / `EXPECT_REL` / `EXPECT_TRUE` invocations across all cases, not the number of cases. The binary exits 0 if and only if no check failed. `ctest` currently reports `100% tests passed, 0 tests failed out of 9`.

**The suite passes in both storage precisions.** 11647 / 0 in the `double` release build and 11647 / 0 in the diagnostic `CHROMO_STATE_FLOAT32` build, with **no tolerance loosened anywhere**. `test_release_production_table_domain_corners` is no longer skipped in either precision: its ULP walk steps in the storage type and nudges against the solver's own EOS acceptance rather than a test-side recomputation of the endpoint energy. (The 29 check failures once documented under the retired `CHROMO_STATE_FLOAT64` option were an artifact of promoting `Vec` alone while the update chain still cast through `float`; they do not occur now that precision is a single coherent knob. See Storage precision in @ref architecture.)

**Release validation script.**

```bash
scripts/release_validation.sh [output_dir]
```

`scripts/release_validation.sh` runs two checks, both through the normal release invocation `scripts/run_chromo_realtime.sh <out> full no-ionization model_column - 20.0 no-cooling` with no numerical-method overrides:

1. *Default-path equivalence* (`CHROMO_T_END=20`). It runs `model_column` twice — once with nothing set, once with `ISO_RIEMANN=swmf-godunov ISO_RECONSTRUCTION=lnrho-v-lnp` set explicitly — and requires the two snapshot files to be identical line for line after discarding the `# EOS_MODE=` header line (which embeds the self-referential sidecar path), and additionally compares the `.gamma_diag` sidecars byte for byte. This is the check that `swmf-godunov` and `(ln rho, V, ln p)` really are the defaults and that no environment variable is needed to obtain them. The script exits non-zero if the snapshots differ.
2. *Release smoke* (`CHROMO_T_END=100`, frame cadence 10 s). It runs the N=500 / R4 release mesh (661 actual cells) for 100 physical seconds with the face-flux and outer-conduction sidecars enabled (`CHROMO_FACE_FLUX_DIAG`, `CHROMO_OUTER_COND_DIAG`, `CHROMO_EOS_COUNTS`), echoes the `# ns=` sidecar header (which records the numerical method actually used — a default release run shows `roe=0 swmf_godunov=1`), the `termination=` line and the `godunov.*` fallback tally, **fails the run if any face fell back to Rusanov**, then runs `util/pressure_reconstruction_diag.py` to report `eps_p`, `Q_V`, `Q_M`, `Q_eff`, `F_eff`, `T_top`, `p_top` and the negative-velocity cell count, writing them to a `_metrics.json` sidecar.

Outputs default to `outputs/model_column`; intermediates go to `tmp/release_validation`; the Python interpreter defaults to `.venv/bin/python` and is overridable with `PY`. Check 1 and the zero-fallback gate in check 2 are pass/fail. Check 2 **reports** its metrics rather than asserting thresholds on them, so a human has to read the numbers.

As measured on the current `double`-precision release tree: check 1 bit-identical, and check 2 `godunov.faces = 23507804`, `godunov.fallbacks = 0`, `termination=end_time`, `eps_p` max `0.0002 %`, `Q_eff = 7.26e-04`, mean `V = 14.070 m/s` with no negative-velocity cell, `F_eff = 1.0414e-09 kg m^-2 s^-1`, `T_top = 21838.6 K`, mass drift `-3.454e-04`. The script is run only against the release build; the diagnostic `CHROMO_STATE_FLOAT32` tree is deliberately not exercised by it.

## What the C++ suite covers

`tests/chromo_tests.cpp` is a single self-contained hand-rolled harness — no external test framework. Cases are plain `static void test_*()` functions invoked from `main()` through a `RUN(fn)` macro that stringifies the function name, so the reported name is the function name. There is no auto-registration and no per-test isolation: several cases mutate process environment variables and rely on an explicit cleanup helper, so ordering is load-bearing.

**Helper and mesh-metric invariants.** Algebraic identities of the indexing and limiter helpers and of the derived mesh metrics: `test_scalar_to_get_scalar_inverse`, `test_ip1_im1_interior_shift`, `test_flux_lim_is_minmod`, `test_flux_lim_mc3`, `test_broadcast_static_metric_cache`, `test_metric_caches_center_to_center`, `test_openmp_thread_safety_contracts`.

**EOS, Saha and `Gamma1` table checks.** The log-domain Saha closure, both EOS inversions (`rho(p,T)` against a bisection reference and `T(rho,p)`), the table loader and its interpolation, the production table's domain corners, and the face-state consistency properties that motivated pressure-based reconstruction: `test_saha_ionization_fraction_log_domain`, `test_equilibrium_density_from_pressure_matches_bisection`, `test_equilibrium_temperature_from_density_pressure`, `test_eos_gamma_table_loader_and_interpolation`, `test_release_production_table_domain_corners`, `test_stage3_equilibrium_mixture_closure`, `test_pressure_face_state_matches_temperature_face_state`, `test_pressure_face_hint_does_not_change_face_state`, `test_reconstruction_preserves_constant_pressure`, `test_smooth_saha_gradient_reduces_pressure_mismatch`, `test_final_serial_log_aware_eos_and_gamma1_contracts`, `test_final_serial_one_update_and_fallback_regression`.

**Release-state round trips.** Packing and decoding `U = (rho, rho u, E)` and deriving the carrier quantities that are never stored: `test_mixture_state_round_trip_and_carrier_derivation`, `test_cons_prim_roundtrip`, `test_mixture_field_cache_contract`, `test_model_column_gamma_ic_uses_cell_average_temperature`.

**Roe reference-flux properties.** The mixture Roe characteristic flux is no longer the release flux, but its correctness is still covered because it is the controlled comparison solver. Its characteristic projection, the face-local Rusanov fallback, the equilibrium face states, the acoustic characteristic speed and the spectral radius: `test_mixture_roe_flux_and_rusanov_fallback`, `test_roe_local_face_flux_and_projection`, `test_stage5_6_equilibrium_face_flux_and_sound_speed`, `test_stage9_acoustic_characteristic_speed`, `test_spectral_radius_uniform_at_rest`, `test_cal_max_v_matches_spectral_radius`.

**Release SWMF Godunov flux.** Three cases cover the release face solver. `test_swmf_exact_riemann_solver` checks the ported SWMF exact Riemann solver (`src/single_fluid/exact_rs.hpp`) against **Toro's five standard tests** to the six significant digits their star pressures are published to, checks that sampling outside the wave fan returns the data states, that `x/t = 0` for Sod reproduces the tabulated `rho*_L = 0.42632`, and that the degenerate cases this model actually lives in are exact: a stationary contact carries **identically zero** mass, a uniform state in uniform motion is reproduced to round-off for either sign of the velocity, and the vacuum and inadmissible-input guards report rather than return. It then runs a first-order Godunov solve of Sod with the SWMF flux assembly and confirms both L1 convergence against the exact solution and that a nonzero specific energy offset — the device that carries the non-ideal EOS content — is advected passively and changes the result by nothing beyond round-off. `test_swmf_godunov_flux_mode` exercises the production RHS path: on a healthy stratified column every face must be solved exactly with zero fallbacks, the captured flux must still be the one the continuity row differences, the result must differ measurably from both Roe and Rusanov, a uniform state must reduce to the physical flux, one real `mixture_advance` step must leave every cell decodable, and a state driven past the two-rarefaction escape speed must fall back without producing a NaN **and be counted** in `Grid::godunov_stats`. `test_release_timestep_uses_frozen_signal_speed` pins the CFL rule: `mixture_timestep` must return exactly `CFL * min_i ds_i/(|u_i| + sqrt(max(Gamma1_i, 5/3) p_i/rho_i))` on the release path and the plain `Gamma1` form under the Roe reference flux, the release step can never exceed the reference step, no cell may exceed the requested CFL measured against the frozen wave speeds, and the frozen/equilibrium speed gap must actually be large (>1.15x) somewhere in the column so the rule is not a disguised no-op.

**Conservation and fixed-point properties.** A uniform state must be an exact fixed point, mass must be conserved, the CFL selection must be respected, the conduction operator must conserve energy on an irregular mesh, and the diagnostic captures must reproduce the production face quantities: `test_uniform_state_is_fixed_point`, `test_mass_conservation_on_uniform_state`, `test_rk4_uniform_fixed_point`, `test_cal_dt_respects_cfl`, `test_mixture_physical_conduction`, `test_face_flux_capture_matches_production_continuity`, `test_outer_conduction_capture_matches_solver_face`, `test_outer_conduction_physical_face_is_mesh_independent`, `test_irregular_uniform_state_preserved`, `test_irregular_linear_pressure_gradient`, `test_irregular_cfl_selection`, `test_irregular_conduction_conserves_energy`, `test_irregular_boundary_hse`, `test_numerical_diffusivity_uniform_reduction`, `test_numerical_diffusivity_nonuniform_face_scaling`, `test_uniform_regression_refinement_unset`, and on the two-fluid side `test_stage_e_conserves_rho_tot_per_cell`, `test_stage_e_mass_momentum_chi_H_drain`, `test_stage_e_kinetic_equilibrium`, `test_stage_e_fixed_point_at_kinetic_equilibrium`, `test_stage_e_large_dt_bounded`, `test_advance_euler_no_op_when_ionization_disabled`, `test_trac_broadening_conserves_kappa_lambda`, `test_beam_heating_partitions_by_heat_capacity`, `test_coronal_heating_partitions_by_heat_capacity`.

**Solver-isolation guards.** These are the tests that keep the two solvers from being confused for one another, in both directions:

| Test | What it asserts |
| --- | --- |
| `test_model_column_release_scenario_defaults_and_overrides` | The `model_column` preset sets none of `ISO_TWO_FLUID`, `ISO_IONIZATION`, `ISO_COOLING`, `ISO_TRAC`, `ISO_CORONA`, `ISO_CHEAT`, `ISO_QFLUX`, `ISO_TBOOST`, `ISO_NUMERICAL_DIFFUSIVITY_MULT`; `ISO_GAMMA` makes `make_scenario("model_column", "")` throw; the retired `model_isentropic` alias throws; and `model_column_ic` itself throws on each rejected knob. |
| `test_mixture_integrator_path_and_guards` | On a release Grid the legacy two-fluid entry points `advance_Euler_state`, `advance_Euler_explicit_state`, `advance_RK4`, `cal_dt_i` and `rhs_explicit_state` must all throw. |
| `test_model_gentle_two_fluid_path` | The reverse direction: `model_gentle` never sets `GAMMA_TABLE`, produces a seven-row state with `single_fluid` false, and `mixture_decode` on it must throw. |
| `test_model_column_release_numerics_defaults` | The release defaults really are `swmf_godunov_flux` and `pressure_reconstruct`, and `roe_characteristic_flux` is off unless asked for; `ISO_RIEMANN=roe-local` selects it *and* clears `swmf_godunov_flux`, so the three fluxes are mutually exclusive; neither Gamma-path flux may leak onto the two-fluid path; invalid `ISO_RIEMANN` / `ISO_RECONSTRUCTION` values, and Gamma-only choices requested off the Gamma path, must throw. |

**Scenario initial-condition and boundary checks.** The Saha-HSE column and its ghost ladder, the upper-boundary temperature decoupling, the C7 closures in both upper-BC modes, the open-BC mirror/extrapolation patterns, the PFSS data-file parser and grid invariants, the analytic canopy `B` profile, and the static-refinement mesh builder: `test_stage8_gamma_model_column_saha_hse_and_ghosts`, `test_upper_bc_hydro_conduction_temperature_decoupling`, `test_model_c7_bc_discrete_hse_inner_mach_capped_outer`, `test_model_c7_tr_jump_bc`, `test_apply_open_bcs_mirrors_and_extrapolates`, `test_apply_open_bcs_outer_velocity_halving_pattern`, `test_open_bcs_stability_under_model_c7`, `test_advance_euler_stable_under_model_c7`, `test_advance_euler_stable_under_model_c7_with_ionization`, `test_parser_reads_synthetic_file`, `test_pfss_ic_grid_invariants`, `test_pfss_advance_euler_stable`, `test_analytic_canopy_B_profile`, `test_mesh_disabled_is_uniform`, `test_mesh_refined_diagnostic_profile`, `test_mesh_outer_refined_profile`, `test_mesh_outer_profile_edge_cases`, `test_refine_params_env_aliases`, `test_refined_mesh_ic_runs`, `test_outer_refined_mesh_ic_runs`.

**Two further clusters** do not belong to any of the layers above. *Microphysics formulas checked against the write-up*: `test_pressure_relation_eq38`, `test_nu_in_collision_formula`, `test_kappa_e_eq53`, `test_kappa_n_eq59`, `test_kappa_n_dominates_in_chromosphere`, `test_ionization_rate_S_voronov_at_1e4K`, `test_ionization_rate_S_increases_with_T`, `test_recombination_rate_alpha_hummer`, `test_photoionization_rate_default_uniform`, `test_stage_e_photoionization_drives_low_T_equilibrium`, `test_stage_e_no_chi_H_drain_from_photoionization`, `test_model_c7_photoionization_uses_route_b_closure`, `test_trac_cutoff_detection_and_limiter`, `test_beam_heating_rate_profile`, `test_coronal_heating_rate_profile` — these exercise @ref two_fluid_solver physics, not the release path. *Run control and output cadence*: the five `test_step_cap_*` cases, `test_termination_classification`, and the five `test_output_schedule_*` cases.

## Registered CTest cases

Nine tests are registered in `CMakeLists.txt`.

| Test | What it verifies |
| --- | --- |
| `chromo_tests` | The whole C++ unit suite above; passes if every check passes. |
| `eos_production_table` | Runs `eos_validate_production_table` on `data/eos/gamma1_hydrogen_v1.dat` — the production `Gamma1` table is self-consistent with the C++ Saha implementation. |
| `eos_production_table_checksum` | `cmake/check_eos_checksum.cmake`: the SHA-256 of the production table matches its `.sha256` sidecar, so the table cannot drift silently. |
| `release_reject_enable_te` | `cmake/check_release_guard.cmake` with `ENABLE_TE=1`: a `model_column` run must exit non-zero *and* name `ENABLE_TE is a legacy two-fluid setting`. |
| `release_reject_single_fluid` | Same guard with `SINGLE_FLUID=0`; must fail with `SINGLE_FLUID is a legacy two-fluid setting`. |
| `release_reject_two_fluid` | Same guard with `ISO_TWO_FLUID=1`; must fail with `ISO_TWO_FLUID is a legacy two-fluid setting`. |
| `release_reject_ionization` | Same guard with the `ionization` CLI argument; must fail with `finite-rate ionization was requested`. |
| `release_output_sidecar` | `cmake/check_gamma_output.cmake`: a short `ISO_NS=48` `model_column` run must write a `.gamma_diag` sidecar carrying `EOS_MODE=gamma_table`, the table's SHA-256 as provenance, the gravity-potential convention, and the ten named physical columns; and the raw snapshot header must declare exactly **3** conserved rows, i.e. `(rho, rho u, E)`. |
| `conductive_flux_python` | `tests/test_conductive_flux.py` — the Python conductive-flux diagnostic reproduces the `physics.hpp` electron and neutral conductivity formulas exactly, and gives byte-identical flux for the fixed-gamma and Gamma-table input paths on identical arrays. Registered only when `.venv/bin/python` exists. |

The four `release_reject_*` cases check the *reason* for the failure, not merely that it failed, so a release run cannot be silently reconfigured into something else and cannot fail the guard for an unrelated reason.

## Known limitations

**N500 long-duration evaporation mass flux is not established as better than 10 % grid-converged.** Do not describe the release mesh as grid-converged, and do not quote sub-10 % long-time mass-flux accuracy at N=500 without new evidence. This is a **resolution** limitation and is untouched by the storage-precision cutover; keep the two distinct. The release fixes the *numerical method*; it does not by itself fix the *resolution required for a particular quantitative scientific claim*. Those are separate questions and must be kept separate in any write-up. The measurements behind this: a Rusanov-era acceptance matrix showed a ~10–11 % N500 → N1000 gap in evaporation mass flux at 4000 s, and a complete N=1000 4000 s comparison has never been run under either flux. Closer to the current configuration the agreement is better — roughly 2.5–3.4 % between N=500 and N=1000 at 500–1000 s — but that is a different, shorter measurement and does not retire the limitation.

**The release is not a well-balanced scheme in the constructive sense.** It removes the dominant discretization inconsistencies and leaves a small, resolution-dependent residual. Completing the MUSCL-Hancock momentum source and making the lower ghost ladder second order (trapezoidal, EOS-closed) cut the raw discrete hydrostatic face mass-flux defect by about 560x. At the release `double` storage that residual is `1.78e-13 kg m^-2 s^-1` at N=500 (20 s, conduction off, release Godunov flux), with its maximum in the upper transition region at 2149 km rather than on the base face, and the 20 s maximum residual velocity is `5.73e-3 m/s`. **The residual is now demonstrably second-order truncation error, not a representation floor**: it is about `1e8` times the double round-off scale of the stored state, and refining N=500 → N=1000 reduces it 4.55x (order 2.19), with order 1.95 at the base face and 1.90 in the residual velocity. The earlier statement that no clean truncation-error trend survived between the two resolutions was a property of `float` storage — in `float` the defect was `1.08e-12`, sat at `1.07x` its own round-off floor, and *grew* under refinement (apparent order `-1.0`). The 5-step residual-velocity ladder now behaves the same way: at `ISO_NS` = 48, 125, 250, 500, 1000, 2000 the maximum residual velocity runs `8.11e-1, 7.62e-2, 1.88e-2, 4.99e-3, 1.32e-3, 3.32e-4 m/s`, a clean unbroken order ≈ 2 over a 41x refinement with no flattening. The `float` ladder (`8.10e-1, 7.70e-2, 1.79e-2, 5.84e-3, 3.38e-3, 4.59e-3 m/s`) tracks it only to `ISO_NS = 250`, degrades to order 0.79 by 1000, and reverses at 2000. So the well-balancing bound is resolution dependent in the ordinary second-order sense — the release mesh sits on a convergent trend rather than on a round-off floor. This remains acceptable because the release ships one mesh, and is still not a general claim. What survives at the base face in `double` is a real defect rather than noise — `f_diff = 1.29e-13` of an `f_total = 1.25e-13`, i.e. Godunov dissipation acting on a genuinely nonzero reconstructed jump — and it is the same first-interior-face artifact listed below. Evidence: `docs/studies/numerics/reference_free_release_recap.md`.


**Reference-free does not mean boundary-free.** The interior hydrodynamic operator needs no frozen global hydrostatic reference state, but the boundary closures still capture fixed reservoir data once from the initial condition (`kOuterPRef`, `kInnerTRef`, `kInnerPRes0/1`), as any truncated-domain problem must.

**`CHROMO_CFL=0.50` is validated for this workstation and this model shape only.** It is the validated production runtime default for the reduced release configuration — coarse-equivalent N=500, R4 outer refinement, 661 actual cells, physical-face 22,000 K conductive boundary, physical conduction only. The conservative hard-coded solver default remains 0.25 and 0.25 stays the comparison reference; do not change the C++ default. **Rerun the sweep before trusting 0.50 on different hardware or a materially different model shape.** The acceptance sweep itself was run on one mesh and one machine, every step of every run in it was acoustic-limited, and scenarios with additional timestep limiters (radiative cooling, beam heating, electron energy) were never exercised. At N=500 the 0.50-vs-0.25 differences at 500 s are 1.6 % in face mass flux, 2.0 % in outer conductive flux, 2.7 % in integrated conductive input and 2.1 % in thermal-layer thickness.

**Residual velocities now converge in time, but no m/s-scale noise floor has been established.** In the pre-cutover `float` build the residual velocity field did not converge as `dt -> 0` — halving the timestep *increased* the level-to-level difference. That was float32 cancellation, and the cutover removes it: re-measured on the release mesh at 20 s, `double` converges at observed order ≈ 2.0 in `L1` (level-to-level `5.34e-3, 1.41e-3, 4.54e-4 m/s` for CFL 0.50→0.25→0.125→0.0625) while a `float` control at the identical configuration still shows apparent order `-1.2` to `-1.3` and stays pinned near `1.0-1.2 m/s` at every CFL. The mechanism is unambiguous in the fields: over an 8x step-count range `double` is stable (mean `|v|` 2.4015 → 2.4084 m/s) and `float` degrades (2.418 → 3.271 m/s), i.e. error growing with the number of additions. **What this does and does not license:** a velocity self-convergence ladder is a valid discriminator again, and the ≈1 m/s floor is gone (≈7e-3 m/s at CFL 0.50). But the CFL acceptance matrix has **not** been re-run in `double`, and no positive noise floor has been measured for m/s-scale chromospheric velocities, so quantitative claims at that scale still require their own measurement. Evidence: the RESOLVED block in `docs/studies/validation/long_run_output_and_cfl_validation.md` §7.4.

**The float32 mass-update stagnation is retired, and it is why the release stores `double`.** With a `float` conserved state, `mixture_advance` computing `next = state + dt*rhs` could not move the lower chromosphere at all: in the 4000 s N=500 run `dt*(-df/ds)` is about `2.6e-17 kg m^-3 s^-1` against a float32 ULP of `rho ~ 9e-10` of about `5.6e-17`, so **619 of 660 cells sat below half a ULP** and the addition rounded straight back to `rho`. There is no compensated summation and no deviation formulation, so that was permanent stagnation rather than slow drift, and its visible signature was a conservative face mass flux 2.81x larger at the 1600 km base than at the top that never flattened, although continuity forces a quasi-steady constant-area column to carry a height-independent flux. The controlling ratio is the stratified layer's mass-loading timescale over the CFL step, `tau/dt ~ 3.5e7`, against `1/u = 1.7e7` for `float` and `9.0e15` for `double`: `float` fails by 2x, `double` has about eight orders of magnitude of margin. **Double precision is now the release default and the only supported production configuration** (see Storage precision in @ref architecture). Measured at 4000 s it brings the base/top flux ratio from 2.811 to **1.003**, closes the column mass budget (predicted/observed -214.9 → **0.994**), and leaves **no** cell frozen (619/660 → 0/660). The 1600–1850 km spread falls from 119 % to **13 %** — not to the 2 % once reported, which came from an incoherent hybrid build; the surviving 13 % is a broad monotone decline that relaxes very slowly (14.3 % at 500 s → 13.0 % at 4000 s) and is **not** the first-face artifact. It is evidence that the column is still relaxing at 4000 s, which is a duration question, not a precision one. Top-of-domain evaporation observables also move by about **28 %** (top velocity 13.18 → 9.50 m/s, top `f_total` 3.753e-10 → 2.702e-10), downward and in the physically coherent direction: a lower chromosphere that cannot drain sustains a mass supply the column does not have. **Any evaporation number quoted from a float32 run is high by of order 30 % and must be re-measured.** What double fixes is that broad, spurious height gradient and the masking of the truncation trend — **not** the localized first-interior-face artifact below, and **not** the resolution limitation above. Evidence: `docs/studies/numerics/state_precision_release_cutover.md`, with the original diagnostic in `docs/studies/numerics/float32_precision_control_experiment.md`.

**A localized first-interior-face artifact is real, precision independent, and open.** Distinct from, and much narrower than, the retired float32 height gradient: face 0 carries about 1900x the interior median `|f_diff|`, its `|f_diff|` is about 16 % of the local `f_total`, and the cell-0 `rho V` overshoots `f_total[0]` by about 28 %. Double precision leaves it unchanged (16.3 % versus 15.9 %), so it is a ghost, reconstruction or boundary-adjacent scheme effect rather than round-off, and it has not been diagnosed. **The precision cutover must not be described as fixing the lower-boundary mass flux**; this one survives it.

**No formal order of accuracy is claimed, in space or in time.** The hydro and conduction stages are composed by first-order Lie splitting, not Strang. The corrector still evaluates the momentum source at `U^n`, so in the source-only limit the update degenerates to forward Euler. No second-order temporal accuracy is claimed for the hydrodynamic stage or for the composition, and no formal second-order spatial convergence is claimed for the release mesh.

**The lower boundary is mass-impermeable.** The 1600 km reservoir wall means mass cannot cross it, so the column cannot be replenished from the unresolved chromosphere below. A gravity-aware characteristic variant recovered about 17 % of the 1000 s column mass loss while changing the evaporation velocity by under 0.4 %. That is nearly irrelevant to velocities and oscillations but **not** irrelevant to the column mass budget, which matters directly for any coronal mass-loading interpretation. This must be revisited if base mass supply becomes a quantity of interest.

**The conduction drive has a first-order Dirichlet boundary layer.** The gap between the imposed 22,000 K wall and the realized top-cell temperature closes only as `O(ds)`, so the realized conductive input into the domain is mesh dependent. This is the mechanism behind the evaporation-rate convergence limitation above: an unconverged drive guarantees an unconverged evaporation rate. Relatedly, the N=500 thermal-layer thickness is about 11 % larger than N=1000 at 500 s, although the layer still spans about 50 cells, and the N=1000 comparison was carried to 500 s rather than 1000 s — only the selected N=500 candidate received the 1000 s long-run stability gate.

**Quote the conserved face flux, not the cell-centred `rho V`.** The cell-centred momentum diagnostic still contains substantial reconstruction-scale ripple; the conserved face mass flux is far smoother and is the quantity the production decisions were made on.

**The `.faceflux` divergence anomaly is explained, and it was not a capture bug.** The sidecar's flux *values* were always consistent; its *divergence* disagreed with the density evolution recorded in the same file by factors of order 60x in the lower chromosphere, which is exactly the float32 mass-update stagnation described above — the fluxes were right and the density *response* was forbidden. In double precision the same integrated divergence matches the recorded density change to 6 %, so `R_total` and `R_eff` in `util/face_flux_diag.py` are usable as mass-budget diagnostics on the release build. A mass-budget residual built from a pre-cutover run, or from the diagnostic `CHROMO_STATE_FLOAT32` build, still measures the representation floor rather than the scheme. Separately and still unexplained: a normalized-curvature floor of about `1.4e-4` in `Q_eff` does not decrease with refinement; it is harmless at that magnitude.

**A lower-chromosphere standing mode is present and only numerically damped.** A broad velocity oscillation of the sub-transition-region column is a conduction-excited, boundary-trapped standing mode. Its period is `O(1e3 s)` but is **not** resolution- or scheme-converged, and the restoring mechanism that sets the period has not been identified. The model as posed contains no physical damping for it: radiative cooling is off in the release, there is no viscosity, and conduction acts on `1e5 s`. It was assessed as a non-blocker for the release because it does not change the transition-region evaporation flux, but it should be understood before time-resolved chromospheric velocities are quoted.

**The 1600 km lower boundary is a mid-chromospheric truncation, not the photosphere.** It is the height below which the Saha/LTE closure is not intended to be trusted; the atmosphere below it is outside the model, and the two inner ghosts are a numerical closure for the MUSCL stencil, not a claim that LTE extends downward. Describing it as photospheric misstates the model.

**The source-split momentum form assumes constant area.** The release runs `B = 1`, where this is not an issue. On a variable-area flux tube the split form is not generally equivalent to the conservative area-weighted operator; that path is not part of the release.

## The single-temperature assumption, measured

Not a limitation but a positive result, and the only place in this document where a release
assumption is checked against a *different solver* rather than against a refinement of itself.

The release is a single-fluid equilibrium mixture with one temperature. A separate experimental
two-temperature solver (`src/two_temp/`, scenario `model_column_2t`, @ref two_temp_physics) was
built and run on exactly the release column — same atmosphere, mesh, gravity, table, CFL and
boundary geometry, one shared IC/BC implementation, so the temperature split is the only variable.
At `T_e = T_i` the two closures degenerate exactly: at t = 0 the two states are **bit-identical**
in all three release rows.

Over 4000 s the two temperatures never meaningfully separate. The global maximum over every cell
and all 712,031 steps is `max |T_e - T_i|/T_i = 3.28e-3`, reached at t ~ 0.022 s in the top cell
during the boundary switch-on and decaying as `1/t`; the quasi-steady value is `9.0e-6`
(0.0715 K). The two-temperature run reproduces the release top-of-domain velocity and mass flux to
**0.176 %** and the column pressure to `3.1e-5`. The mechanism is timescale separation: `tau_eq` is
`9.2e-5` to `4.0e-3 s` against a conduction/evaporation timescale of order `1e3 s`.

**What this licenses and what it does not.** It validates the single-temperature closure *for this
configuration* — a 1600–2153 km chromospheric column with physical conduction only and a 22 kK
electron-conducted outer reservoir. It says nothing about a flare, beam-heated, much hotter or much
more rarefied regime, where `tau_eq` rises and the electron heat flux is far larger. Only N=500 was
run, so the decoupling itself is not grid-converged. Evidence:
`docs/studies/conduction/two_temperature_test_study.md`.

**The experimental solver itself has no registered CTest case** and is not exercised by
`scripts/release_validation.sh` — it is a test study, not a shipped path. What *is* verified after
adding it is that the release is unaffected: all 9 CTest cases pass and `release_validation.sh`
completes clean with the two-temperature sources compiled in.

## Rejected approaches

Both of these were prototyped and measured, and **must not be reintroduced without new evidence**.

**The full SWMF `ModTransitionRegion` momentum pressure/gravity split.** Measurably no better than the predictor source fix alone on this constant-area column, and not generally conservative for variable-area flux-tube geometry. The predictor-gravity half of it, which *is* what the release adopted, is the part that carried the benefit.

**A gravity-aware characteristic (non-reflecting) lower boundary.** It changed the hydrostatic residual not at all and the evaporation velocity by at most 0.4 %; the 1600 km face is not an important acoustic reflector for this model, because the disturbance that reaches it is a slow quasi-static pressure adjustment rather than acoustic ringing. Its one measurable benefit — recovering about 17 % of the 1000 s column mass loss — is recorded above as a standing limitation of the wall rather than a reason to adopt the characteristic form.

## Evidence documents

The authoritative records live in `docs/`. Current code is always the source of truth for active behaviour; these recaps record how the decisions were reached.

| Document | Status | Establishes |
| --- | --- | --- |
| `docs/studies/numerics/reference_free_release_recap.md` | **Current** — the authority for the release as shipped | Retirement of the frozen equilibrium reference, the MUSCL-Hancock momentum source, the trapezoidal EOS-closed lower ghost ladder, the ~560x hydrostatic defect reduction, and the resolution-dependent well-balancing bound. |
| `docs/studies/numerics/model_column_release_numerics_recap.md` | Current for the reconstruction and the limitation; its flux and well-balancing rows are historical | `(ln rho, V, ln p)` reconstruction, physical-only conduction, and the boxed N500 evaporation-flux limitation. Its Roe-as-release-flux rows were superseded by the Godunov cutover. |
| `docs/studies/numerics/swmf_godunov_flux_experiment.md` | **Current** — the authority for the release numerical flux | The SWMF exact-Riemann Godunov flux, its verification against Toro's five standard tests and the passive-offset property, the measured Roe-versus-Godunov comparison on the release model, the frozen-composition rationale, and the flux-consistent CFL. |
| `docs/studies/conduction/coarse_model_column_physical_conduction_recap.md` | Current for the mesh, boundary and CFL; its `eq_wb` and hydro-flux rows are historical | The N=500 / R4 661-cell mesh, the 22,000 K physical-face conductive boundary, zero artificial conduction, the cell-averaged initial condition, and the CFL 0.50 acceptance at N=500. |
| `docs/studies/validation/long_run_output_and_cfl_validation.md` | Sweep configuration retired; conclusions current, except that its velocity noise-floor finding was measured in the pre-cutover float32 build | The CFL 0.50 acceptance sweep, the `CHROMO_T_END` / `CHROMO_FRAME_DT` / `CHROMO_STEP_CAP` run-control semantics, and the m/s velocity noise-floor finding. |
| `docs/studies/numerics/pressure_reconstruction_recap.md` | Superseded by the promotion it recommended | Why `(ln rho, V, ln p)` replaced `(ln rho, V, ln T)`: face-pressure mismatch reduced by ~350x at N=500. |
| `docs/studies/numerics/state_precision_release_cutover.md` | **Current** — the authority for the release storage precision | The production cutover to a `double` conserved state: the single `chromosphere::Real` knob, the `tau/dt` versus `1/u` necessity argument, the SWMF/BATS-R-US `PRECISION = ${DOUBLEPREC}` consistency argument, the retirement of `CHROMO_STATE_FLOAT64`, the diagnostic-only `CHROMO_STATE_FLOAT32` build, and the post-cutover test and release-validation status. |
| `docs/studies/numerics/float32_precision_control_experiment.md` | Current as the original diagnostic; its `CHROMO_STATE_FLOAT64` build no longer exists | That the broad lower-chromosphere mass-flux gradient was float32 conserved-state round-off and not physics, the below-half-a-ULP mass-update mechanism, the closure of the `.faceflux` divergence anomaly, and the precision-independence of the first-interior-face artifact. It is the experiment the cutover acted on. |
| `docs/studies/numerics/roe_n1000_lower_chromosphere_sloshing_diagnosis.md` | Current for the mode; its `.faceflux` divergence anomaly is explained as pre-cutover float32 round-off; measured under the Roe reference flux | The lower-chromosphere standing mode, the reflecting inner boundary, and the mass-reservoir behaviour. |
| `docs/studies/validation/resolution_convergence_scan_recap.md` | Retired configuration | That the cell-centred `rho V` ripple is spatial truncation error, and that no resolution in a uniform ns = 500/1000/2000 scan produced a grid-independent evaporation rate. |
| `docs/studies/numerics/outer_tr_refinement_recap.md` | Retired configuration | The `outer` refinement profile and the `MESH_MAX_RATIO` grading fix; its evaporation numbers predate `q_num = 0` and must not be quoted. |
| `docs/studies/conduction/physical_conductive_flux_recap.md` | Diagnostic definition current, numbers historical | The plotted conductive flux is `q_phys = -(kappa_e + kappa_n) dT/ds` and excludes every solver-only term. |
| `docs/studies/numerics/static_local_refinement.md` | Current | The static local-refinement mesh builder — not AMR. |
| `docs/studies/performance/openmp_parallelization_recap.md` | Current | Thread-count scaling and the bit-identical serial/OpenMP result. |
| `docs/studies/conduction/two_temperature_test_study.md` | Current; **EXPERIMENTAL, non-release** | That the release's single-temperature assumption holds on the release column to about 0.2 % in the evaporation observables, measured against a separate `T_e != T_i` solver. Also the characteristic derivation of the two-temperature boundary conditions, and two bugs found and fixed in the experimental path. |
