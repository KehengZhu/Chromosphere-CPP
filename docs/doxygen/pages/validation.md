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
Passed: 11362
Failed: 0
```

That is **104 test cases and 11362 individual checks, 0 failures**. Note that `Passed:` counts individual `EXPECT_NEAR` / `EXPECT_REL` / `EXPECT_TRUE` invocations across all cases, not the number of cases. The binary exits 0 if and only if no check failed. `ctest` currently reports `100% tests passed, 0 tests failed out of 9`.

**Release validation script.**

```bash
scripts/release_validation.sh [output_dir]
```

`scripts/release_validation.sh` runs two checks, both through the normal release invocation `scripts/run_chromo_realtime.sh <out> full no-ionization model_column - 20.0 no-cooling` with no numerical-method overrides:

1. *Default-path equivalence* (`CHROMO_T_END=20`). It runs `model_column` twice — once with nothing set, once with `ISO_RIEMANN=roe-local ISO_RECONSTRUCTION=lnrho-v-lnp` set explicitly — and requires the two snapshot files to be identical line for line after discarding the `# EOS_MODE=` header line (which embeds the self-referential sidecar path), and additionally compares the `.gamma_diag` sidecars byte for byte. This is the check that `roe-local` and `(ln rho, V, ln p)` really are the defaults and that no environment variable is needed to obtain them. The script exits non-zero if the snapshots differ.
2. *Release smoke* (`CHROMO_T_END=100`, frame cadence 10 s). It runs the N=500 / R4 release mesh (661 actual cells) for 100 physical seconds with the face-flux and outer-conduction sidecars enabled (`CHROMO_FACE_FLUX_DIAG`, `CHROMO_OUTER_COND_DIAG`, `CHROMO_EOS_COUNTS`), echoes the `# ns=` sidecar header (which records the numerical method actually used) and the `termination=` line, then runs `util/pressure_reconstruction_diag.py` to report `eps_p`, `Q_V`, `Q_M`, `Q_eff`, `F_eff`, `T_top`, `p_top` and the negative-velocity cell count, writing them to a `_metrics.json` sidecar.

Outputs default to `outputs/model_column`; intermediates go to `tmp/release_validation`; the Python interpreter defaults to `.venv/bin/python` and is overridable with `PY`. Only check 1 is a pass/fail gate. Check 2 **reports** its metrics rather than asserting thresholds on them, so a human has to read the numbers.

## What the C++ suite covers

`tests/chromo_tests.cpp` is a single self-contained hand-rolled harness — no external test framework. Cases are plain `static void test_*()` functions invoked from `main()` through a `RUN(fn)` macro that stringifies the function name, so the reported name is the function name. There is no auto-registration and no per-test isolation: several cases mutate process environment variables and rely on an explicit cleanup helper, so ordering is load-bearing.

**Helper and mesh-metric invariants.** Algebraic identities of the indexing and limiter helpers and of the derived mesh metrics: `test_scalar_to_get_scalar_inverse`, `test_ip1_im1_interior_shift`, `test_flux_lim_is_minmod`, `test_flux_lim_mc3`, `test_broadcast_static_metric_cache`, `test_metric_caches_center_to_center`, `test_openmp_thread_safety_contracts`.

**EOS, Saha and `Gamma1` table checks.** The log-domain Saha closure, both EOS inversions (`rho(p,T)` against a bisection reference and `T(rho,p)`), the table loader and its interpolation, the production table's domain corners, and the face-state consistency properties that motivated pressure-based reconstruction: `test_saha_ionization_fraction_log_domain`, `test_equilibrium_density_from_pressure_matches_bisection`, `test_equilibrium_temperature_from_density_pressure`, `test_eos_gamma_table_loader_and_interpolation`, `test_release_production_table_domain_corners`, `test_stage3_equilibrium_mixture_closure`, `test_pressure_face_state_matches_temperature_face_state`, `test_pressure_face_hint_does_not_change_face_state`, `test_reconstruction_preserves_constant_pressure`, `test_smooth_saha_gradient_reduces_pressure_mismatch`, `test_final_serial_log_aware_eos_and_gamma1_contracts`, `test_final_serial_one_update_and_fallback_regression`.

**Release-state round trips.** Packing and decoding `U = (rho, rho u, E)` and deriving the carrier quantities that are never stored: `test_mixture_state_round_trip_and_carrier_derivation`, `test_cons_prim_roundtrip`, `test_mixture_field_cache_contract`, `test_model_column_gamma_ic_uses_cell_average_temperature`.

**Roe flux properties.** The mixture Roe characteristic flux, its characteristic projection, the face-local Rusanov fallback, the equilibrium face states, the acoustic characteristic speed and the spectral radius: `test_mixture_roe_flux_and_rusanov_fallback`, `test_roe_local_face_flux_and_projection`, `test_stage5_6_equilibrium_face_flux_and_sound_speed`, `test_stage9_acoustic_characteristic_speed`, `test_spectral_radius_uniform_at_rest`, `test_cal_max_v_matches_spectral_radius`.

**Conservation and fixed-point properties.** A uniform state must be an exact fixed point, mass must be conserved, the CFL selection must be respected, the conduction operator must conserve energy on an irregular mesh, and the diagnostic captures must reproduce the production face quantities: `test_uniform_state_is_fixed_point`, `test_mass_conservation_on_uniform_state`, `test_rk4_uniform_fixed_point`, `test_cal_dt_respects_cfl`, `test_mixture_physical_conduction`, `test_face_flux_capture_matches_production_continuity`, `test_outer_conduction_capture_matches_solver_face`, `test_outer_conduction_physical_face_is_mesh_independent`, `test_irregular_uniform_state_preserved`, `test_irregular_linear_pressure_gradient`, `test_irregular_cfl_selection`, `test_irregular_conduction_conserves_energy`, `test_irregular_boundary_hse`, `test_numerical_diffusivity_uniform_reduction`, `test_numerical_diffusivity_nonuniform_face_scaling`, `test_uniform_regression_refinement_unset`, and on the two-fluid side `test_stage_e_conserves_rho_tot_per_cell`, `test_stage_e_mass_momentum_chi_H_drain`, `test_stage_e_kinetic_equilibrium`, `test_stage_e_fixed_point_at_kinetic_equilibrium`, `test_stage_e_large_dt_bounded`, `test_advance_euler_no_op_when_ionization_disabled`, `test_trac_broadening_conserves_kappa_lambda`, `test_beam_heating_partitions_by_heat_capacity`, `test_coronal_heating_partitions_by_heat_capacity`.

**Solver-isolation guards.** These are the tests that keep the two solvers from being confused for one another, in both directions:

| Test | What it asserts |
| --- | --- |
| `test_model_column_release_scenario_defaults_and_overrides` | The `model_column` preset sets none of `ISO_TWO_FLUID`, `ISO_IONIZATION`, `ISO_COOLING`, `ISO_TRAC`, `ISO_CORONA`, `ISO_CHEAT`, `ISO_QFLUX`, `ISO_TBOOST`, `ISO_NUMERICAL_DIFFUSIVITY_MULT`; `ISO_GAMMA` makes `make_scenario("model_column", "")` throw; the retired `model_isentropic` alias throws; and `model_column_ic` itself throws on each rejected knob. |
| `test_mixture_integrator_path_and_guards` | On a release Grid the legacy two-fluid entry points `advance_Euler_state`, `advance_Euler_explicit_state`, `advance_RK4`, `cal_dt_i` and `rhs_explicit_state` must all throw. |
| `test_model_gentle_two_fluid_path` | The reverse direction: `model_gentle` never sets `GAMMA_TABLE`, produces a seven-row state with `single_fluid` false, and `mixture_decode` on it must throw. |
| `test_model_column_release_numerics_defaults` | The release defaults really are `roe_characteristic_flux` and `pressure_reconstruct`; invalid `ISO_RIEMANN` / `ISO_RECONSTRUCTION` values, and Gamma-only choices requested off the Gamma path, must throw. |

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

**N500 long-duration evaporation mass flux is not established as better than 10 % grid-converged.** Do not describe the release mesh as grid-converged, and do not quote sub-10 % long-time mass-flux accuracy at N=500 without new evidence. The release fixes the *numerical method*; it does not by itself fix the *resolution required for a particular quantitative scientific claim*. Those are separate questions and must be kept separate in any write-up. The measurements behind this: a Rusanov-era acceptance matrix showed a ~10–11 % N500 → N1000 gap in evaporation mass flux at 4000 s, and a complete Roe N=1000 4000 s comparison has never been run. Closer to the current configuration the agreement is better — roughly 2.5–3.4 % between N=500 and N=1000 at 500–1000 s — but that is a different, shorter measurement and does not retire the limitation.

**The release is not a well-balanced scheme in the constructive sense.** It removes the dominant discretization inconsistencies and leaves a small, resolution-dependent residual. Completing the MUSCL-Hancock momentum source and making the lower ghost ladder second order (trapezoidal, EOS-closed) cut the raw discrete hydrostatic face mass-flux defect by about 560x, to `1.33e-12 kg m^-2 s^-1` at N=500. That is comparable to the estimated float32 round-off scale of the stored state — about 1.32x it at N=500 and 2.63x at N=1000 — so **no clean truncation-error trend survives between the two resolutions tested**, and the residual is not demonstrated to *be* the representation floor. Reconstruction, the EOS inversion, the non-uniform mesh and the boundary closures all contribute at that level and two resolutions cannot separate them. The practical consequence is that the well-balancing bound is now resolution dependent: the five-step maximum residual velocity runs `7.58e-1, 7.53e-2, 1.77e-2, 5.82e-3, 3.38e-3 m/s` at `ISO_NS` = 48, 125, 250, 500, 1000. This is acceptable because the release ships one mesh; it would not be acceptable as a general claim.

**Reference-free does not mean boundary-free.** The interior hydrodynamic operator needs no frozen global hydrostatic reference state, but the boundary closures still capture fixed reservoir data once from the initial condition (`kOuterPRef`, `kInnerTRef`, `kInnerPRes0/1`), as any truncated-domain problem must.

**`CHROMO_CFL=0.50` is validated for this workstation and this model shape only.** It is the validated production runtime default for the reduced release configuration — coarse-equivalent N=500, R4 outer refinement, 661 actual cells, physical-face 22,000 K conductive boundary, physical conduction only. The conservative hard-coded solver default remains 0.25 and 0.25 stays the comparison reference; do not change the C++ default. **Rerun the sweep before trusting 0.50 on different hardware or a materially different model shape.** The acceptance sweep itself was run on one mesh and one machine, every step of every run in it was acoustic-limited, and scenarios with additional timestep limiters (radiative cooling, beam heating, electron energy) were never exercised. At N=500 the 0.50-vs-0.25 differences at 500 s are 1.6 % in face mass flux, 2.0 % in outer conductive flux, 2.7 % in integrated conductive input and 2.1 % in thermal-layer thickness.

**Chromospheric velocities at the m/s scale are a noise floor, not a signal.** The residual velocity field does not converge as `dt -> 0`; halving the timestep *increases* the level-to-level difference, which is consistent with float32 cancellation in a residual carried on top of a hydrostatic balance. Any study whose conclusion depends on m/s-scale chromospheric velocities must treat those values as noise regardless of CFL. The stored state is float32 throughout.

**No formal order of accuracy is claimed, in space or in time.** The hydro and conduction stages are composed by first-order Lie splitting, not Strang. The corrector still evaluates the momentum source at `U^n`, so in the source-only limit the update degenerates to forward Euler. No second-order temporal accuracy is claimed for the hydrodynamic stage or for the composition, and no formal second-order spatial convergence is claimed for the release mesh.

**The lower boundary is mass-impermeable.** The 1600 km reservoir wall means mass cannot cross it, so the column cannot be replenished from the unresolved chromosphere below. A gravity-aware characteristic variant recovered about 17 % of the 1000 s column mass loss while changing the evaporation velocity by under 0.4 %. That is nearly irrelevant to velocities and oscillations but **not** irrelevant to the column mass budget, which matters directly for any coronal mass-loading interpretation. This must be revisited if base mass supply becomes a quantity of interest.

**The conduction drive has a first-order Dirichlet boundary layer.** The gap between the imposed 22,000 K wall and the realized top-cell temperature closes only as `O(ds)`, so the realized conductive input into the domain is mesh dependent. This is the mechanism behind the evaporation-rate convergence limitation above: an unconverged drive guarantees an unconverged evaporation rate. Relatedly, the N=500 thermal-layer thickness is about 11 % larger than N=1000 at 500 s, although the layer still spans about 50 cells, and the N=1000 comparison was carried to 500 s rather than 1000 s — only the selected N=500 candidate received the 1000 s long-run stability gate.

**Quote the conserved face flux, not the cell-centred `rho V`.** The cell-centred momentum diagnostic still contains substantial reconstruction-scale ripple; the conserved face mass flux is far smoother and is the quantity the production decisions were made on.

**Two diagnostics are not trustworthy quantitatively.** The `.faceflux` sidecar's flux *values* are consistent, but its *divergence* is inconsistent with the density evolution recorded in the same file, by factors of order 60x in the lower chromosphere down to 0.6x near the top. `R_total` and `R_eff` in `util/face_flux_diag.py` should not be quoted until that capture's stage consistency is checked. Separately, a normalized-curvature floor of about `1.4e-4` in `Q_eff` does not decrease with refinement and is unexplained; it is harmless at that magnitude.

**A lower-chromosphere standing mode is present and only numerically damped.** A broad velocity oscillation of the sub-transition-region column is a conduction-excited, boundary-trapped standing mode. Its period is `O(1e3 s)` but is **not** resolution- or scheme-converged, and the restoring mechanism that sets the period has not been identified. The model as posed contains no physical damping for it: radiative cooling is off in the release, there is no viscosity, and conduction acts on `1e5 s`. It was assessed as a non-blocker for the release because it does not change the transition-region evaporation flux, but it should be understood before time-resolved chromospheric velocities are quoted.

**The 1600 km lower boundary is a mid-chromospheric truncation, not the photosphere.** It is the height below which the Saha/LTE closure is not intended to be trusted; the atmosphere below it is outside the model, and the two inner ghosts are a numerical closure for the MUSCL stencil, not a claim that LTE extends downward. Describing it as photospheric misstates the model.

**The source-split momentum form assumes constant area.** The release runs `B = 1`, where this is not an issue. On a variable-area flux tube the split form is not generally equivalent to the conservative area-weighted operator; that path is not part of the release.

## Rejected approaches

Both of these were prototyped and measured, and **must not be reintroduced without new evidence**.

**The full SWMF `ModTransitionRegion` momentum pressure/gravity split.** Measurably no better than the predictor source fix alone on this constant-area column, and not generally conservative for variable-area flux-tube geometry. The predictor-gravity half of it, which *is* what the release adopted, is the part that carried the benefit.

**A gravity-aware characteristic (non-reflecting) lower boundary.** It changed the hydrostatic residual not at all and the evaporation velocity by at most 0.4 %; the 1600 km face is not an important acoustic reflector for this model, because the disturbance that reaches it is a slow quasi-static pressure adjustment rather than acoustic ringing. Its one measurable benefit — recovering about 17 % of the 1000 s column mass loss — is recorded above as a standing limitation of the wall rather than a reason to adopt the characteristic form.

## Evidence documents

The authoritative records live in `docs/`. Current code is always the source of truth for active behaviour; these recaps record how the decisions were reached.

| Document | Status | Establishes |
| --- | --- | --- |
| `docs/reference_free_release_recap.md` | **Current** — the authority for the release as shipped | Retirement of the frozen equilibrium reference, the MUSCL-Hancock momentum source, the trapezoidal EOS-closed lower ghost ladder, the ~560x hydrostatic defect reduction, and the resolution-dependent well-balancing bound. |
| `docs/model_column_release_numerics_recap.md` | Current for the method; its well-balancing rows are historical | The Roe-local mixture flux, `(ln rho, V, ln p)` reconstruction, physical-only conduction, and the boxed N500 evaporation-flux limitation. |
| `docs/coarse_model_column_physical_conduction_recap.md` | Current for the mesh, boundary and CFL; its `eq_wb` and hydro-flux rows are historical | The N=500 / R4 661-cell mesh, the 22,000 K physical-face conductive boundary, zero artificial conduction, the cell-averaged initial condition, and the CFL 0.50 acceptance at N=500. |
| `docs/long_run_output_and_cfl_validation.md` | Sweep configuration retired; conclusions current | The CFL 0.50 acceptance sweep, the `CHROMO_T_END` / `CHROMO_FRAME_DT` / `CHROMO_STEP_CAP` run-control semantics, and the float32 velocity noise-floor finding. |
| `docs/pressure_reconstruction_recap.md` | Superseded by the promotion it recommended | Why `(ln rho, V, ln p)` replaced `(ln rho, V, ln T)`: face-pressure mismatch reduced by ~350x at N=500. |
| `docs/roe_n1000_lower_chromosphere_sloshing_diagnosis.md` | Current for the mode and the `.faceflux` bug | The lower-chromosphere standing mode, the reflecting inner boundary, the mass-reservoir behaviour, and the untrusted `.faceflux` divergence. |
| `docs/resolution_convergence_scan_recap.md` | Retired configuration | That the cell-centred `rho V` ripple is spatial truncation error, and that no resolution in a uniform ns = 500/1000/2000 scan produced a grid-independent evaporation rate. |
| `docs/outer_tr_refinement_recap.md` | Retired configuration | The `outer` refinement profile and the `MESH_MAX_RATIO` grading fix; its evaporation numbers predate `q_num = 0` and must not be quoted. |
| `docs/physical_conductive_flux_recap.md` | Diagnostic definition current, numbers historical | The plotted conductive flux is `q_phys = -(kappa_e + kappa_n) dT/ds` and excludes every solver-only term. |
| `docs/static_local_refinement.md` | Current | The static local-refinement mesh builder — not AMR. |
| `docs/openmp_parallelization_recap.md` | Current | Thread-count scaling and the bit-identical serial/OpenMP result. |
