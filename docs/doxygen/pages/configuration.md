# Configuration reference {#configuration}

Everything a run is configured with lives in exactly two places: the seven positional command-line arguments of `chromo_main`, and the environment. There is no input deck, no configuration file, and no runtime parser. The authoritative sources for this page are `src/main.cpp`, `src/core/run_control.hpp`, `scenarios/scenario.cpp`, `scenarios/model_column.cpp`, `scenarios/mesh.cpp`, and, for the legacy two-fluid scenarios, `legacy/scenarios/legacy_scenarios.cpp`, `legacy/scenarios/model_gentle.cpp`, `legacy/scenarios/model_c7.cpp`, `legacy/scenarios/model_flare.cpp` and `legacy/scenarios/pfss_field_line.cpp`; no environment variable is read anywhere in `src/` outside the driver `src/main.cpp`.

## Build options

CMake options, set at configure time (`cmake -S . -B <dir> -D<option>=ON`). Build each non-default combination into its own tree so the release build is never shadowed. See @ref architecture for the build layout.

| Option | Default | Effect |
| --- | --- | --- |
| `CHROMO_ENABLE_OPENMP` | `OFF` | Links OpenMP and defines `CHROMO_USE_OPENMP=1`. The production configuration is `ON` (tree `build_omp/`, which the run scripts use); the thread count comes from `OMP_NUM_THREADS` at run time. |
| `CHROMO_ENABLE_LEGACY` | `OFF` | Also builds the frozen two-fluid solver and its scenarios under `legacy/` (`model_gentle`, `model_c7`, `model_flare`, `analytic_canopy`, `pfss_field_line`) and their tests, and defines `CHROMO_ENABLE_LEGACY=1`. With it `OFF` none of that code is compiled or linked, and requesting a legacy scenario fails immediately with a message naming `-DCHROMO_ENABLE_LEGACY=ON`. See `legacy/README.md`. |
| `CHROMO_STATE_FLOAT32` | `OFF` | DIAGNOSTIC ONLY: stores the conserved state in `float`. No release number may be quoted from it. |
| `CHROMO_STATE_FLOAT64` | — | Retired; configuring with it is a hard CMake error (double is the default). |

## Command line

Usage: `chromo_main [output_path] [mode] [ionization] [scenario] [data_path] [time_mult] [cooling]`. All seven are positional and optional; to set a later one you must supply the earlier ones (pass `-` or `""` for an unused `data_path`).

| # | Argument | Accepted values | Default | Meaning |
| --- | --- | --- | --- | --- |
| 1 | `output_path` | any path | `outputs/output.txt` | Snapshot file. Sidecars are derived from it: `<out>.gamma_diag`, `<out>.faceflux`, `<out>.outercond`. See @ref io_formats. A one-line summary is always appended to `outputs/output.log`. |
| 2 | `mode` | `full`, `explicit` | `full` | `explicit` selects the two-fluid explicit-only path (`advance_Euler_explicit_state`, R_I identically zero; legacy build only). Any string other than `explicit` means `full`. Combining `explicit` with a loaded `Gamma1` table is a hard error: "the release solver has no explicit-only mode". |
| 3 | `ionization` | `ionization`, `no-ionization` | `no-ionization` for `model_column`, `ionization` for every other scenario | Sets `Grid::enable_ionization` (finite-rate Stage-E network). Any string other than the literal `no-ionization` counts as on. In release mode `ionization` is a hard error: the release closure is Saha equilibrium. |
| 4 | `scenario` | `model_column`, `model_column_2t`; with `-DCHROMO_ENABLE_LEGACY=ON` also `model_gentle`, `model_c7`, `model_flare`, `analytic_canopy`, `pfss_field_line` | `model_column` | Selects the IC/BC pair via `make_scenario`. An unknown name is a hard error, and so is a legacy name in a build without the legacy code (the message names `-DCHROMO_ENABLE_LEGACY=ON`). See @ref scenario_reference. |
| 5 | `data_path` | path to a field-line `.dat` | `""` | Required by `pfss_field_line` (hard error if empty). Ignored by every other scenario. |
| 6 | `time_mult` | float | `1.0` | Multiplier on the derived total time `10 * sum(ds_i) / 2e4 m/s`. Also drives the legacy step cap `max(10000, 10000*time_mult)` and the legacy frame stride `max(10, 10*max(1, time_mult/5))`. Both legacy rules are superseded when `CHROMO_T_END` / `CHROMO_STEP_CAP` / `CHROMO_FRAME_DT` are set. |
| 7 | `cooling` | `cooling`, `no-cooling` | `no-cooling` for `model_column`, `cooling` for every other scenario | Sets `Grid::enable_radiative_cooling`. Any string other than the literal `no-cooling` counts as on. |

The differing defaults for arguments 3 and 7 are decided by `scenario_name == "model_column"` alone, evaluated before either argument is read, so `chromo_main out.txt` and `chromo_main out.txt full` both run the release column with ionization and cooling off.

Note on argument 7 in release mode: `model_column_ic` sets `Grid::enable_radiative_cooling = false` (and rejects `ISO_COOLING` outright), so passing `cooling` to `model_column` is silently overridden to off rather than rejected. The release timestep contains no radiative-cooling stage in any case — nothing in `src/single_fluid/` reads `enable_radiative_cooling`, `enable_ionization`, `enable_trac` or `enable_beam_heating`.

## Environment variables

Two parsing conventions are used and they behave differently.

- Driver booleans (`CHROMO_*` in `src/main.cpp`) use `env_enabled()`: unset takes the stated default; the values `0`, `false` and `off` disable; **every other value enables, including an empty string**.
- Scenario knobs (`ISO_*`, `FLARE_*`, `GENTLE_*`, `C7_*`, `PFSS_FLARE`) use a local `env_f()`: the value is parsed with `std::stof` and an unparseable value silently falls back to the default. A knob is "on" when its float value is non-zero.
- The release rejection lists test **presence only** (`std::getenv(name) != nullptr`), so `ISO_COOLING=0` is just as fatal as `ISO_COOLING=1`.

### Run control and output cadence

| Variable | Type / values | Default | Effect | Applies to |
| --- | --- | --- | --- | --- |
| `CHROMO_T_END` | float, seconds | unset | Absolute stop time in physical seconds. Replaces the `time_mult`-derived total time when positive and finite; a non-positive or unparseable value is silently ignored. When set (and `CHROMO_STEP_CAP` is not), the step cap becomes unbounded (1e12) so the run cannot be truncated. | both solvers |
| `CHROMO_STEP_CAP` | positive integer | unset | Explicit step cap; always wins over both other rules. A non-integer, trailing-garbage or non-positive value is **rejected** with a `step-cap error:` message and exit code 1, never silently replaced. Values >= 1e12 are clamped to 1e12. | both solvers |
| `CHROMO_FRAME_DT` | positive float, seconds | unset | Switches the snapshot cadence from a step stride to physical seconds, on an integer frame index so the cadence cannot drift. A non-positive, non-finite or unparseable value is rejected with exit code 1. Required for any long `CHROMO_T_END` run with output enabled. | both solvers |
| `CHROMO_FRAME_STRIDE` | integer >= 1 | `max(10, 10*max(1, time_mult/5))` | Legacy step-stride snapshot cadence. Values below 1 are clamped to 1; unparseable values keep the derived default. Overridden by `CHROMO_FRAME_DT` when that is set. | both solvers |
| `CHROMO_PROGRESS_STRIDE` | integer >= 1 | `100` | Steps between progress lines on stdout. | both solvers |
| `CHROMO_OUTPUT` | boolean | **on** | Writes the main snapshot file. `CHROMO_OUTPUT=0` opens no snapshot file at all. | both solvers |
| `CHROMO_GAMMA_DIAG` | boolean | **on** | Writes the `<out>.gamma_diag` sidecar (10 derived EOS columns per cell). Only has an effect in release mode; the two-fluid solver never writes it. | release |
| `CHROMO_CFL` | float | `0.25` | CFL number passed to `Grid::init`. An unparseable value silently keeps 0.25. See "Validated production defaults" below. | both solvers |
| `CHROMO_BINARY` | path | `build_omp/chromo_main` | Read by `scripts/run_chromo_omp.sh` and `scripts/run_chromo_realtime.sh` only, not by the program. | launchers |

Every run prints `termination=end_time|step_cap|other` on exit, together with `requested_end_time`, `final_physical_time`, `final_step` and `effective_step_cap` with its source (`explicit_step_cap`, `legacy_time_mult` or `end_time`). A `step_cap` termination also emits a warning on stderr. Always check that a long production run ended with `termination=end_time`.

### Diagnostics and profiling

All of these are read-only captures: enabling them changes no numerical result.

| Variable | Type / values | Default | Effect | Applies to |
| --- | --- | --- | --- | --- |
| `CHROMO_PROFILE` | boolean | **off** | Enables region timers (CFL, decode, RHS, conduction, boundary, output), EOS-inversion and conduction-iteration statistics and the active timestep limiter, printed at shutdown by `print_runtime_profile`. | both solvers |
| `CHROMO_EOS_COUNTS` | boolean | **off** | Counts EOS table operations and prints `eos.gamma1_queries`, `eos.temperature_logs`, `eos.n_h_logs` at shutdown. | both solvers |
| `CHROMO_FACE_FLUX_DIAG` | boolean | **off** | Writes the `<out>.faceflux` sidecar: the production total-mass face flux (the `mix::RHO` continuity row) and its central/diffusive split, copied out of `mixture_rhs_explicit`. **Requires the release solver** — with the two-fluid solver it throws "CHROMO_FACE_FLUX_DIAG requires the release solver" and the process aborts. | release only |
| `CHROMO_FACE_FLUX_TOP` | integer >= 0 | `0` | Limits each face-flux record to the top N cells; `0` means all cells. Values >= `ns` also mean all cells. | release only |
| `CHROMO_FACE_FLUX_STRIDE` | integer >= 1 | the effective frame stride | Capture cadence in steps for the face-flux sidecar. The final step of the run is always captured regardless. | release only |
| `CHROMO_OUTER_COND_DIAG` | boolean | **off** | Writes the `<out>.outercond` sidecar: outer-face `T_top`, `T_wall`, face conductivity, `ds_face`, area ratio and `q_face` from the final converged Newton iteration of each captured conduction solve. **Requires the release solver** (same abort behaviour as above). | release only |
| `CHROMO_OUTER_COND_STRIDE` | integer >= 1 | the face-flux stride | Capture cadence in steps for the outer-conduction sidecar. | release only |

Sidecar layouts are documented on @ref io_formats.

### Solver selection

| Variable | Type / values | Default | Effect |
| --- | --- | --- | --- |
| `GAMMA_TABLE` | path to a `Gamma1` table | supplied unconditionally by `model_column` and `model_column_2t` as `data/eos/gamma1_hydrogen_v1.dat` (only the path is overridable) | **A loaded table is what selects a Saha/`Gamma1` solver** — the release single-fluid one, or the experimental two-temperature one when the scenario sets `Grid::two_temperature`. The driver SHA-256s the file and records the hash in the `gamma_diag` header. Setting it together with `ISO_GAMMA` is a hard error; setting it with `mode=explicit` is a hard error; setting it with any scenario other than `model_column` or `model_column_2t` is a hard error. Without a table the driver runs the historical two-fluid solver, which exists only in a `-DCHROMO_ENABLE_LEGACY=ON` build. |
| `SINGLE_FLUID` | float, non-zero = on | unset | Two-fluid research knob: slaves neutrals to the ion fluid. Read only when no `Gamma1` table is loaded. **Rejected (exit 1) in release mode if present at all.** |
| `ENABLE_TE` | float, non-zero = on | unset | Two-fluid research knob: separate electron temperature `T_e != T_i`. Read only when no `Gamma1` table is loaded. **Rejected (exit 1) in release mode if present at all.** |

### `model_column` / `model_gentle` model knobs (the `ISO_*` family)

`model_column` and `model_column_2t` share one IC/BC implementation (`scenarios/model_column.cpp`); the legacy `model_gentle` has its own frozen copy of the same column (`legacy/scenarios/model_gentle.cpp`, legacy build only). The three never share a solver. All three install override-preserving defaults with `setenv(..., overwrite=0)` before the IC runs, so an explicit assignment in the environment always wins.

`model_column` presets: `GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat`, `ISO_H_BASE=1600`, `ISO_DH=553`, `ISO_NS=500`, `ISO_HEAT_FLUX=1`, `ISO_T_TOP=22000`, `ISO_HYDRO_T_DECOUPLE=0`, `ISO_REFINE_PROFILE=outer`, `ISO_REFINE_FACTOR=4`, `ISO_REFINE_S_LO_KM=500`, `ISO_REFINE_TRANSITION_KM=20`.

`model_column_2t` presets: **exactly the same list**, deliberately — the experimental two-temperature leg must differ from the release in nothing but the temperature split, or the comparison is not controlled. It adds only `TT_OUTER_TI` (default `neumann`) and rejects `ISO_RIEMANN` / `ISO_RECONSTRUCTION`, which have no meaning for a solver with one flux and one primitive set.

`model_gentle` presets (in `legacy/scenarios/legacy_scenarios.cpp`): `ISO_CORONA=1`, `ISO_H_BASE=1003`, `ISO_HEAT_FLUX=1`, `ISO_COOLING=1`, `ISO_IONIZATION=1`, `ISO_TWO_FLUID=1`, `ISO_NS=600`. It never loads a `Gamma1` table.

Knobs valid in **both** modes:

| Variable | Type | Default (bare, before presets) | Effect |
| --- | --- | --- | --- |
| `ISO_NS` | unsigned integer | `1000` if `ISO_CORONA` is set and non-empty and not `0`, else `600` | Coarse-equivalent cell count read by `peek_ns` before the IC. With refinement on this is the coarse-grid equivalent, not the allocated cell count. |
| `ISO_H_BASE` | float, km | `0` | Height of the lower truncation face above the photosphere. It is a truncation of the modeled domain, not the photosphere. Also sets `Grid::out_base_km`, the height label offset in the output file. |
| `ISO_DH` | float, km | `1303` (or `ISO_CORONA_TOP_KM - ISO_H_BASE` when `ISO_CORONA` is on) | Domain span above the base. |
| `ISO_HEAT_FLUX` | float, non-zero = on | `0` | `0` = Stage-1 adiabatic relaxation (conduction off); non-zero = conduction on. |
| `ISO_RELAX_TIME` | float, seconds | `0` | Length of an initial Stage-1 adiabatic relaxation phase before conduction and the top jump switch on. |
| `ISO_T_TOP` | float, K | `0` (= use the IC top-cell temperature) | Any positive value is an **absolute** imposed outer-boundary temperature and is not re-anchored across a relax phase. Release preset 22000. The column is *very* sensitive to it: the quasi-steady top-of-domain velocity goes as roughly `T_wall^3` and reverses sign near 10 kK — see [Sensitivity to the wall temperature](#config-twall-sensitivity). |
| `ISO_P_TOP` | float, Pa | unset (the C7 column, top-cell `p` = 0.0102877 Pa) | Release only. Rebuilds the initial column with the **same** C7 `T(h)` but a base pressure chosen by secant iteration so the IC top cell sits at this pressure; both boundary reservoirs (`kOuterPRef`, `kInnerPRes0/1`) are then captured from the rebuilt column, so the pressure changes consistently through the domain instead of as a one-cell mismatch at the outer face. `p_top` is not proportional to `p_base` under Saha (a denser column is less ionized and has a shorter scale height): 0.003 Pa needs 0.22x the C7 base pressure, 0.08 Pa needs 12.1x. The achieved pressures are printed as `[model_column] ISO_P_TOP: ...`. Non-positive or unparsable values are rejected, and so is any use outside the release solver. Unset reproduces the C7 column byte for byte. The TR stays at the C7 height, so across a pressure sweep the initial `L_TR` is fixed by construction and the chromospheric temperature `T_c` stays at the C7 plateau; CFL 0.50 is validated only at the C7 pressure. A 5000 s sweep over `ISO_T_TOP` 18–70 kK × `ISO_P_TOP` 0.003–0.08 Pa ends cleanly with zero Godunov fallbacks. At `ISO_T_TOP >= 50000` with `ISO_P_TOP >= 0.04` (and at 70 kK already at 0.02), however, the N=500 upper chromosphere carries intermittent grid-scale flux noise, with top-face reversals. It is a resolution artefact of the sub-cell TR foot: it is locked to the foot crossing cell faces, is unchanged by CFL 0.25 or by moving the mesh grade, and is 4x weaker at N=1000. The 500 s mean flux is unaffected to 1.5 %. At 18–25 kK and 0.04–0.08 Pa the initial TR has still not relaxed by 5000 s. See `docs/studies/boundaries/tp_sweep_evaporation_closure_validation.md`. |
| `ISO_TJUMP_A` | float | `1.0` | Outer ghost temperature jump `T_ghost1 = a * T_ref`. |
| `ISO_TJUMP_B` | float | `1.0` | Second ghost jump `T_ghost2 = b * T_ghost1`. |
| `ISO_HYDRO_T_DECOUPLE` | float, non-zero = on | `0` | Release preset `0`: one fixed outer temperature `ISO_T_TOP` in both hydro ghosts, which conduction also reads as its Dirichlet datum at the ghost centre. `1` lets the **hydro** ghost temperature zero-gradient-extrapolate the live top cell while conduction keeps a fixed wall at the physical outer face (the release until 2026-09-22; `model_column_2t` preset). See [Outer temperature](#config-outer-temperature). |
| `ISO_OUTER_T_HYDRO_WALL` | float, non-zero = on | `0` | Diagnostic only. In decoupled mode, pins the outer **hydro** ghost temperature back to the conduction wall (`T_hydro = T_cond = ISO_T_TOP`) while the wall stays a physical-face conduction datum. Only useful together with `ISO_OUTER_P_NEUMANN` — see [Which variable carries the outer condition](#config-outer-closure-swap). |
| `ISO_OUTER_P_NEUMANN` | integer, `0`\|`1`\|`2` | `0` (fixed reservoir back-pressure) | Diagnostic only. Replaces the fixed outer back-pressure `kOuterPRef` with a Neumann pressure closure extrapolated off the live top cell: `1` = `dp/ds = 0`, `2` = `dp/ds = -rho g` (the release ladder re-anchored on the live top cell). Rejected unless the hydro ghost temperature is fixed, since otherwise the face would impose nothing. Mode `1` **destroys the model** — see [Which variable carries the outer condition](#config-outer-closure-swap). |
| `ISO_VCAP` | float, Mach | `0.1` | Outer outflow velocity cap as a fraction of the sound speed. |
| `ISO_COND_KIRCHHOFF` | float, non-zero = on | `0` | Diagnostic only. Replaces the harmonic-mean face conductivity of the Spitzer electron part with the Kirchhoff secant `(2/7)(kappa_e,R T_R - kappa_e,L T_L)/(T_R - T_L)`, clipped to the two cell values. That secant is exact for a steady `kappa ∝ T^{5/2}` profile however steep it is, whereas the harmonic mean follows the cold side of a TR foot narrower than a cell. Tests whether the grid-locked TR pulses at hot walls (≥ 45 kK) come from the face conductivity. The neutral part and every other path are unchanged; unset reproduces the release byte for byte. |
| `ISO_COND_ARITH` | float, non-zero = on | `0` | Diagnostic only. Arithmetic mean of the two total cell conductivities at every conduction face, on any mesh. At the base of the transition region neutral conduction dominates the cold side and Spitzer conduction the hot side, and the harmonic mean is pinned to the smaller one; the arithmetic mean lies near the Kirchhoff value of the total conductivity there. Cannot be combined with `ISO_COND_KIRCHHOFF` (which applies the Kirchhoff secant to the Spitzer part only and keeps the harmonic mean for the neutral part, so it *lowers* that face conductivity). Unset reproduces the release byte for byte. |
| `ISO_INNER_T_NEUMANN` | float, non-zero = on | `0` (Dirichlet) | Conduction lower BC: Dirichlet fixed reservoir temperature vs Neumann (insulating). On the two-temperature leg it applies to **both** channels. |
| `TT_OUTER_TI` | `neumann` \| `dirichlet` | `neumann` | **`model_column_2t` only.** Outer-face conduction BC on the **heavy-particle** temperature. `neumann` (the physically and characteristically motivated default) makes the top insulating for the heavy channel, because the fully ionized plasma above the domain has no neutral-hydrogen conductive flux to supply. `dirichlet` imposes the same 22,000 K wall on `T_i` that the electron channel gets, as the controlled comparison. The **electron** channel keeps the Dirichlet reservoir either way. Any other value throws. See @ref numerics-two-temp-bc. |

Knobs read only by the **two-fluid** path (`model_gentle`, legacy build only); each is a hard error in release mode unless noted:

| Variable | Type | Default | Effect |
| --- | --- | --- | --- |
| `ISO_GAMMA` | float | `Grid::gamma_mono` (5/3) | Fixed adiabatic index of the two-fluid solver. Rejected three separate times in release mode: by `make_scenario` for `model_column`, by the `GAMMA_TABLE` branch of the driver, and by `model_column_ic`. Not to be confused with the fixed `gamma = 5/3` of the release Godunov flux: that is the frozen-composition index of the *face Riemann problem only*, is not user-settable, and does not replace the Saha closure anywhere else in the solver. |
| `ISO_TWO_FLUID` | float, non-zero = on | `0` | Evolve ion and neutral as separate fluids. Rejected in release mode by the driver. |
| `ISO_IONIZATION` | float, non-zero = on | `0` | Stage-E finite-rate ionization/recombination network. Rejected in release mode by the driver. |
| `ISO_COOLING` | float, non-zero = on | `0` | Radiative sink; only takes effect when `ISO_HEAT_FLUX` is also on. Rejected in release mode. |
| `ISO_CORONA` | float, non-zero = on | `0` | Extend the domain into a resolved ~1 MK corona. Rejected in release mode. |
| `ISO_CORONA_TOP_KM` | float, km | `10000` | Coronal top used to derive the default `ISO_DH` when `ISO_CORONA` is on. Read but inert in release mode (`ISO_CORONA` cannot be set there, and `ISO_DH` is preset). |
| `ISO_TRAC` | float | `-1` = auto (`ISO_COOLING || ISO_CORONA`) | Force TRAC TR-broadening on or off. Rejected in release mode. |
| `ISO_TRAC_TCHROM` | float, K | `2.0e4` | TRAC chromospheric base temperature. Read only on the two-fluid branch. |
| `ISO_TRAC_TCMAXFRAC` | float | `0.2` | TRAC cutoff cap as a fraction of the peak temperature. Read only on the two-fluid branch. |
| `ISO_QFLUX` | float, non-zero = on | `0` | Impose the coronal conductive flux `q(T)` as a Neumann outer BC instead of the ghost-T jump. Rejected in release mode. |
| `ISO_QFLUX0` | float, W/m^2 | C7's own top-gradient flux `kappa0 T^{5/2} dT/ds` | Base flux. Read only when `ISO_QFLUX` is on. |
| `ISO_QFLUX_ENHANCE` | float | `1.0` | Ramp amplitude on the imposed flux. |
| `ISO_QFLUX_T_ON` | float, s | `600.0` | Ramp start time. |
| `ISO_QFLUX_RAMP` | float, s | `120.0` | Ramp duration. |
| `ISO_CHEAT` | float, non-zero = on | `0` | Ambient volumetric coronal heating `H(s) = E0 exp(-(s-s0)/s_H)`. Turns the boundary `q(T)` off when enabled. Rejected in release mode. |
| `ISO_CHEAT_TMAX_MK` | float, MK | `1.5` | Target apex temperature used to calibrate `E0` from the RTV/AS2002 static scaling. |
| `ISO_CHEAT_SH_FRAC` | float | `0.4` | Heating scale height as a fraction of the column length. |
| `ISO_CHEAT_E0` | float, W/m^3 | the calibrated value | Direct override of the peak heating rate. |
| `ISO_CHEAT_S0_MM` | float, Mm | `0.0` | Heating anchor offset. |
| `ISO_CHEAT_ENHANCE` | float | `1.0` | Ramp amplitude (`>1` drives gentle evaporation). |
| `ISO_CHEAT_T_ON` | float, s | `0.0` | Ramp start time. |
| `ISO_CHEAT_RAMP` | float, s | `30.0` | Ramp duration. |
| `ISO_TBOOST` | float | `1.0` | One-time IC coronal superheat factor, applied above the TR through a tanh-in-T selector; only active when `> 1`. Rejected in release mode. |
| `ISO_TBOOST_TC` | float, K | `3.0e5` | Selector centre temperature. |
| `ISO_TBOOST_W` | float, K | `1.0e5` | Selector width. |
| `ISO_NUMERICAL_DIFFUSIVITY_MULT` | float | `0` | Multiplier on the historical mesh-scaled artificial conduction (`1` restores the legacy `C_num = 2000 m/s`, 4x in the resolved-corona testbed). Rejected in release mode; the release conduction operator is structurally physical-only and sets `numerical_diffusivity_per_length = 0`. |

### Mesh refinement

Static local refinement, not AMR: the mesh is generated once at IC time and never regridded. Read by `refine_params_from_env` in `scenarios/mesh.cpp`, and used **only** by `model_column` / `model_gentle` — no other scenario calls the builder, and PFSS field lines carry their own externally specified mesh. The `ISO_REFINE_*` aliases are applied after the `GRID_REFINE_*` values and therefore override them.

| Variable | Alias | Type | Default | Effect |
| --- | --- | --- | --- | --- |
| `GRID_REFINE_FACTOR` | `ISO_REFINE_FACTOR` | float | `1.0` | Fine/coarse resolution ratio. `<= 1` disables refinement, and the builder then returns the exact uniform `ns_coarse` grid. |
| `GRID_REFINE_S_LO_KM` | `ISO_REFINE_S_LO_KM` | float, km from the inner face | `0.0` | Inner edge of the refined band. |
| `GRID_REFINE_S_HI_KM` | `ISO_REFINE_S_HI_KM` | float, km | `700.0` | Outer edge of the refined band. Used by the `lower` profile only; the `outer` profile always runs fine to the domain top. |
| `GRID_REFINE_TRANSITION_KM` | `ISO_REFINE_TRANSITION_KM` | float, km | `100.0` | **Minimum** width of the geometric fine/coarse grade. The actual width is set by the adjacent-cell-width ratio cap `MESH_MAX_RATIO = 1.1`, so a requested width the cap cannot support is widened, never narrowed. |
| `GRID_REFINE_PROFILE` | `ISO_REFINE_PROFILE` | `lower` or `outer`, case-insensitive | `lower` | `lower`: coarse, fine over `[s_lo, s_hi]`, grade back to coarse, coarse. `outer`: coarse, grade over `[s_lo - tau, s_lo]`, fine from `s_lo` to the domain top. Any other string keeps the current value. |

Refinement only **adds** cells: `ISO_NS` remains the coarse-grid-equivalent count and `Grid::resize` grows the allocation to match the generated face grid. The release preset (`outer`, factor 4, `s_lo = 500 km`, transition 20 km, `ISO_NS=500`, `ISO_DH=553`) yields 661 actual cells. If the requested window does not fit strictly inside the domain, the builder silently falls back to the uniform grid.

### Flare, PFSS and C7 knobs

All of these belong to legacy two-fluid scenarios and are read only in a `-DCHROMO_ENABLE_LEGACY=ON` build.

`FLARE_*` are read by `model_flare_ic` and, under the `PFSS_FLARE` overlay, by `pfss_ic`; both use the same names so one build serves both.

| Variable | Type | Default | Effect | Applies to |
| --- | --- | --- | --- | --- |
| `FLARE_NS` | unsigned integer | `100` | Cell count, read by `peek_ns`. | `model_flare` |
| `FLARE_BEAM_FLUX` | float, W/m^2 | `5.0e7` | Nonthermal electron beam energy flux (Fisher et al. 1985 middle explosive run). | `model_flare`, `PFSS_FLARE` |
| `FLARE_T_ON` | float, s | `2.0` | Beam switch-on time. | `model_flare`, `PFSS_FLARE` |
| `FLARE_DUR` | float, s | `10.0` | Beam flat-top duration. | `model_flare`, `PFSS_FLARE` |
| `FLARE_RAMP` | float, s | `1.0` | Cosine ramp width at both ends of the beam window. | `model_flare`, `PFSS_FLARE` |
| `FLARE_H_LO` | float, km | `1600.0` | Lower edge of the fixed deposition window. | `model_flare`, `PFSS_FLARE` |
| `FLARE_H_HI` | float, km | `2000.0` | Upper edge of the fixed deposition window. | `model_flare`, `PFSS_FLARE` |
| `FLARE_FREE_OUTFLOW` | float, non-zero = on | `1.0` (on) | Opens the outer face so the supersonic evaporation upflow can leave. Informational under `PFSS_FLARE`, whose BC is already a transparent Neumann outflow. | `model_flare`, `PFSS_FLARE` |
| `FLARE_DIFF_MULT` | float | `1.0` | Multiplier on the grid-scaled numerical diffusivity that stabilizes the evaporation shock. Under `PFSS_FLARE` it multiplies a fixed `2.0e3`. | `model_flare`, `PFSS_FLARE` |
| `FLARE_THICK_TARGET` | float, non-zero = on | `1.0` (on) | Self-consistent Emslie (1978) thick-target beam injected at the apex instead of the fixed height window. Read **only** on the `PFSS_FLARE` path. | `PFSS_FLARE` |
| `FLARE_E_CUT` | float, keV | `20.0` | Beam low-energy cutoff. Read only on the `PFSS_FLARE` path. | `PFSS_FLARE` |
| `FLARE_DELTA` | float | `5.0` | Beam spectral index. Read only on the `PFSS_FLARE` path. | `PFSS_FLARE` |
| `PFSS_FLARE` | float, non-zero = on | `0` (off) | Overlays the `model_flare` physics (ionization, cooling, TRAC, grid-scaled diffusivity, RTV `q(T)`, beam) on a PFSS field line's geometry and BC. | `pfss_field_line` |
| `C7_TJUMP_A` | float | `1.0` | Ghost-T multiplier on the live top cell temperature; `1` means continuous extrapolation. Read only in the `tr_jump_bc` branch, i.e. only by the bare `model_c7` scenario. | `model_c7` |
| `C7_TJUMP_B` | float | `1.0` | Second ghost-T multiplier, same gating. | `model_c7` |

`GENTLE_*` is not a remnant: it is the live ambient-coronal-heating overlay of `pfss_field_line` (the alternative to `PFSS_FLARE`), read in `pfss_ic`. It is unrelated to the `model_gentle` scenario, whose knobs are the `ISO_*` family above. The former `GENTLE_Q_*` / `GENTLE_HEAT_*` / `GENTLE_T_BOOST*` names of the column scenario were renamed to `ISO_QFLUX*` / `ISO_CHEAT*` / `ISO_TBOOST*` and no longer exist.

| Variable | Type | Default | Effect |
| --- | --- | --- | --- |
| `GENTLE` | float, non-zero = on | `0` (off) | Enables the overlay: ionization, radiative cooling, vacuum floor, TRAC, grid-scaled diffusivity and footpoint-anchored coronal heating `H(s)`. No beam. |
| `GENTLE_WELL_BALANCED` | float, non-zero = on | `1.0` (on) | `0` reverts to the pre-fix explicit reconstruction. |
| `GENTLE_DIFF_MULT` | float | `1.0` | Multiplier on the `2.0e3` grid-scaled numerical diffusivity. |
| `GENTLE_TMAX_MK` | float, MK | `2.0` | Target apex temperature for the RTV/AS2002 heating calibration. |
| `GENTLE_SH_FRAC` | float | `0.4` | Heating scale height as a fraction of the loop half-length. |
| `GENTLE_E0` | float, W/m^3 | the calibrated `E_H0` | Direct override of the peak footpoint heating rate. |
| `GENTLE_S0_MM` | float, Mm | `0.0` | Heating anchor offset. |
| `GENTLE_ENHANCE` | float | `1.0` | Slow ramp amplitude; `>1` drives gentle evaporation. |
| `GENTLE_T_ON` | float, s | `0.0` | Ramp start time. |
| `GENTLE_RAMP` | float, s | `30.0` | Ramp duration. |
| `GENTLE_FREE_OUTFLOW` | float, non-zero = on | `1.0` (on) | Open outer face for the gentle upflow. |

### OpenMP runtime

The program reads no OpenMP variable itself except to echo it: it prints `[openmp] compiled=... max_threads=... conduction_team=... requested=<OMP_NUM_THREADS>` at startup, and it aborts if `omp_get_max_threads()` is outside `[1, 256]`. The conduction team is capped at `min(4, max_threads)`.

| Variable | Recommended value | Effect |
| --- | --- | --- |
| `OMP_NUM_THREADS` | `12` on the current Apple Silicon workstation | Thread count. Set to `1` for serial/OpenMP numerical-equivalence checks and `16` only when measuring maximum throughput. |
| `OMP_DYNAMIC` | `FALSE` | Prevents the runtime from silently varying the team size. |
| `OMP_MAX_ACTIVE_LEVELS` | `1` | No nested parallelism. |

`scripts/run_chromo_omp.sh` applies these three as override-preserving defaults; `scripts/run_chromo_realtime.sh` adds low-I/O output defaults (`CHROMO_OUTPUT=0`, `CHROMO_GAMMA_DIAG=0`, `CHROMO_FACE_FLUX_DIAG=0`, `CHROMO_OUTER_COND_DIAG=0`, `CHROMO_EOS_COUNTS=0`, `CHROMO_PROFILE=1`, `CHROMO_PROGRESS_STRIDE=10000`) and, when the fourth positional argument is exactly `model_column`, `CHROMO_CFL=0.50`. Explicit environment assignments always win over all of these.

## Reference-only overrides — never use in production

These exist so a numerical result can be attributed to (or cleared of) a specific scheme choice without hand-editing the scenario. They are not needed for any production run, they print a loud override line to stderr when set, and `scripts/release_validation.sh` checks that the default path reproduces the explicit `swmf-godunov` / `lnrho-v-lnp` run bitwise.

| Variable | Values | Release default | Effect |
| --- | --- | --- | --- |
| `ISO_RIEMANN` | `swmf-godunov`, `roe-local`, `rusanov` | `swmf-godunov` | The release value is the SWMF-style exact-Riemann Godunov flux at frozen composition described under [Numerics](@ref numerics). `roe-local` selects the mixture Roe characteristic flux — the equilibrium-`Gamma1` linearization that was the previous release flux, kept as the controlled comparison solver. `rusanov` restores the local Lax-Friedrichs flux (excessively dissipative at the low Mach numbers of this problem), which is also the automatic per-face fallback of the other two. The three are mutually exclusive, so the face Riemann solver is the only thing this variable changes — **including the timestep**, which `mixture_timestep` sizes from the signal speed of whichever flux is in force (frozen `sqrt(5/3 p/rho)` for `swmf-godunov`, equilibrium `sqrt(Gamma1 p/rho)` otherwise). `swmf-godunov` and `roe-local` require Gamma/Saha mode. Any other value throws `std::invalid_argument`. |
| `ISO_RECONSTRUCTION` | `lnrho-v-lnt`, `lnrho-v-lnp` | `lnrho-v-lnp` | `lnrho-v-lnt` restores the temperature-based primitive set, which closes pressure through the nonlinear EOS and manufactures face-pressure mismatch across the partial-ionization transition. It also turns off the hydrostatic reconstruction, which lives in the `ln p` slot, so the initial condition and ghosts revert to the trapezoidal hydrostatic form. `lnrho-v-lnp` requires Gamma/Saha mode. Any other value throws. |
| `ISO_LIMITER` | `mc3`, `minmod`, `first` | `mc3` with beta = 2 | `minmod` selects the symmetric minmod limiter; `first` sets beta = 0 with MC3, giving phi identically 0, i.e. piecewise-constant first order. Any other value throws. |
| `ISO_MC3_BETA` | float | `2.0` | Read **only** when `ISO_LIMITER=mc3` is explicitly set. Setting `ISO_MC3_BETA` alone has no effect. |

The numeric-method toggles that used to select the "bestwb" path — `ISO_LOG_RECON`, `ISO_MC3`, `ISO_EQ_WB`, `ISO_INNER_WB`, `ISO_INNER_WB_RHO`, `ISO_C7_IC`, `ISO_C7_HSE` — are now **unconditional** in `model_column` / `model_gentle` and are simply ignored if set, so old command lines still run. Log-space `(ln rho, u, ln p)` MUSCL, the MC3/Koren limiter at beta = 2, the C7 IC and its hydrostatic density rebuild are always on (the two-fluid `well_balanced` φ_g pressure correction has no counterpart in the release, whose reconstruction decodes each cell with its own potential; the release's balance at rest comes from its hydrostatic reconstruction instead — see @ref numerics); the equilibrium-reference well-balancing that `ISO_EQ_WB` used to select was **retired** from the release, which is reference-free (`Grid::eq_wb` survives for the two-fluid research solver only and must never be set on a release Grid). The diagnostic toggles `ISO_DIAG_*`, the base sponge `ISO_SPONGE_*`, `ISO_TRAC_TC_FIXED` and the analytic-isentrope toy IC knobs (`ISO_T_BASE`, `ISO_N_BASE`, `ISO_F_ION`) were removed outright and are not read anywhere.

## Rejected in release mode

"Release mode" means a loaded `Gamma1` table *and* a scenario that did not select the experimental two-temperature solver, which in practice means the `model_column` scenario. The guards fire in three places and every one of them stops the run.

`make_scenario` (before anything is allocated), for `model_column` only:

- `ISO_GAMMA` present -> `std::runtime_error`: "model_column is the single-fluid release scenario and derives every thermodynamic index from the Saha/Gamma1 closure; ISO_GAMMA ... has no meaning here. Use the legacy model_gentle scenario for the historical fixed-gamma column."

For `model_column_2t` only:

- `ISO_GAMMA` present -> `std::runtime_error`: "model_column_2t is built on the Saha/Gamma1 closure; ISO_GAMMA ... has no meaning here."

`chromo_main`, in the `GAMMA_TABLE` branch and immediately after (all print `release-solver error:` and exit 1):

- `ISO_GAMMA` present together with `GAMMA_TABLE` -> "GAMMA_TABLE and explicitly supplied ISO_GAMMA are mutually exclusive".
- `mode == "explicit"` -> "the release solver has no explicit-only mode".
- any scenario other than `model_column` or `model_column_2t` -> "the Gamma1/Saha solvers currently support only model_column (release) and model_column_2t (experimental two-temperature)".
- the positional `ionization` argument set to anything but `no-ionization` -> "finite-rate ionization was requested, but the release closure is Saha equilibrium; pass no-ionization".
- `SINGLE_FLUID`, `ENABLE_TE`, `ISO_TWO_FLUID` or `ISO_IONIZATION` **present at all**, regardless of value -> "<name> is a legacy two-fluid setting and has no meaning for the single-fluid release solver".
- a scenario IC returning a state whose element count is not `grid.n_mixture_state` -> "the scenario returned a state of N elements, expected M (rho, rho u, E per cell)". The two-temperature leg has the matching check against `grid.ns*num_of_two_temp_eq` -> "(rho, rho u, E, E_e per cell)".

`model_column_ic`, for every knob in `kNonReleaseKnobs` **present at all**, regardless of value, throwing `std::runtime_error` "<name> is historical two-fluid research configuration and has no meaning for the single-fluid release solver; run the legacy model_gentle scenario for that path (build with -DCHROMO_ENABLE_LEGACY=ON)": `ISO_GAMMA`, `ISO_TWO_FLUID`, `ISO_IONIZATION`, `ISO_COOLING`, `ISO_TRAC`, `ISO_CORONA`, `ISO_CHEAT`, `ISO_QFLUX`, `ISO_TBOOST`, `ISO_NUMERICAL_DIFFUSIVITY_MULT`.

`model_column_ic` also rejects `ISO_HSE_RECON` **present at all** ("ISO_HSE_RECON was removed; the hydrostatic reconstruction is the release scheme and has no switch"): it selected that reconstruction while it was a prototype, and a stale `ISO_HSE_RECON=0` must not silently run the release.

`model_column_ic` also refuses a Grid with no `Gamma1` table loaded ("no Gamma1 table is loaded ... the historical fixed-gamma two-fluid column is the legacy model_gentle scenario"); the two-fluid leg it used to contain is `model_gentle_ic` in `legacy/scenarios/model_gentle.cpp`.

Additionally, `CHROMO_FACE_FLUX_DIAG` and `CHROMO_OUTER_COND_DIAG` are rejected the other way round: both require the release solver and throw if enabled on a two-fluid run.

## Validated production defaults

`CHROMO_CFL=0.50` is the validated production runtime default for the **reduced `model_column` release configuration only**: coarse-equivalent N = 500, R4 outer refinement, 661 actual cells, a fixed 22,000 K outer ghost temperature, and physical conduction only. Evidence: `docs/studies/conduction/coarse_model_column_physical_conduction_recap.md`. That sweep ran with the previous decoupled outer temperature closure. It carries over because the step is limited by the top cell's sound speed and is unchanged (`dt = 5.63e-3 s` under both closures), and the 4000 s single-temperature run ends with `termination=end_time` and zero Godunov fallbacks.

The hard-coded C++ solver default remains `CHROMO_CFL=0.25` (`src/main.cpp` and `Grid::CFL`), 0.25 stays the comparison reference, and **the C++ default must not be changed**. `scripts/run_chromo_realtime.sh` supplies 0.50 as an override-preserving default and only when the scenario argument is exactly `model_column`; every other scenario keeps 0.25 unless `CHROMO_CFL` is given explicitly. Rerun the CFL sweep before trusting 0.50 on different hardware or a materially different model shape.

12 OpenMP threads is the workstation default (performance cores only, measured within about 1 % of the 16-thread maximum-throughput result). Do not assume it is optimal on other hardware, and do not hard-code a thread count in C++, CMake, or the numerical algorithm.

A low-I/O production run is therefore:

```bash
CHROMO_T_END=1000 scripts/run_chromo_realtime.sh out.txt \
    full no-ionization model_column - 20.0 no-cooling
```

and the same run with snapshots at one frame per physical second:

```bash
CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=1 \
scripts/run_chromo_realtime.sh out.txt \
    full no-ionization model_column - 20.0 no-cooling
```

No environment variable is required to select the release solver, its Riemann solver, or its reconstruction: the scenario supplies all of it. See @ref numerics for what those defaults mean numerically and @ref scenario_reference for the scenarios themselves.

Known limitation to keep in mind when choosing a configuration: N500 long-duration evaporation mass flux is **not** established as better than 10 % grid-converged. That is a statement about the resolution required for a particular quantitative claim, not about the release numerical method. See @ref validation.

## Sensitivity to the wall temperature {#config-twall-sensitivity}

`ISO_T_TOP` is the most consequential single knob in the release configuration, so it deserves a number rather than a warning. With everything else fixed at the release configuration of the time (C7 initial condition, N = 500/R4, `ISO_HYDRO_T_DECOUPLE = 1` — the previous decoupled closure, not re-run under the current single-temperature release — physical conduction only, 4000 s, double-precision state) a one-parameter sweep of the wall over 8–100 K gives a quasi-steady top-of-domain velocity of

| `ISO_T_TOP` [kK] | 8 | 10 | 12 | 15 | 18 | **22** | 30 | 45 | 70 | 100 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `V_top` [m/s] | −0.44 | −0.10 | +0.54 | +2.1 | +4.6 | **+9.6** | +27.6 | +93.5 | +310 | +792 |

Four decades of velocity for one decade of wall temperature. The response is a smooth power law on the upflow branch — a fit over `ISO_T_TOP >= 15000` gives `V ~ T_wall^3.09`, with the local exponent falling from 3.7 at the 22 kK release point to 2.6 at 100 kK — and it reverses sign near 10 kK, below which the column drains instead of evaporating. There is no saturation anywhere in the tested range and the `ISO_VCAP` outflow cap never binds (at most 15 % of the cap, at 100 kK).

The mechanism is entirely conductive: the wall thermostats the top cell to within 0.3 %, Spitzer conductivity supplies the `T^2.5`, and a column with no radiative sink must export the resulting conductive flux as the enthalpy flux of the evaporating material, so `V ~ q_wall / [rho_top (h + x chi_H/m_H)]` to within a factor of two over the whole range. Two practical consequences: any quantitative evaporation result is a statement about the chosen wall temperature as much as about the chromosphere, and a wall below about 8 kK drives the top of the column below the `Gamma1` table temperature floor and aborts the run.

These are means over 3500–4000 s, **not steady values**, and the flux keeps falling after them. The release has no chromospheric energy source, so the gas lifted by the evaporation cools adiabatically, the chromosphere contracts, and the transition region lengthens under the fixed-height wall. The mass flux therefore falls as `F = J/L_TR`, drifting 2–6 % per 500 s at 4000 s. The time-invariant quantity is `J = F·L_TR`, which a steady conduction–enthalpy TR integral predicts from the boundary state alone. Its exponent in `T_wall` is a nearly constant 2.6–2.7 over 22–100 kK, so most of the softening of the local `V_top` exponent quoted above comes from reading every wall at the same time. Before the drift, a run needs ~520–840 s for the acoustic transient and, below 30 kK, up to ~2300 s for the initial C7 transition region to relax. Derivations and a run-length table: `docs/studies/boundaries/twall_mass_flux_relaxation_timescales.md`.

Full study, including the time-convergence and CFL caveats: `docs/studies/boundaries/upper_wall_temperature_sensitivity_recap.md`.

## Outer temperature: one fixed ghost-cell temperature {#config-outer-temperature}

The release (`ISO_HYDRO_T_DECOUPLE=0`) puts one fixed temperature, `T_c = ISO_T_TOP = 22 000 K`, in both upper hydro ghosts. Together with the fixed back-pressure `kOuterPRef`, the ghost state is a coronal reservoir at `(p_ref, T_c)`. Conduction reads its Dirichlet datum from the same first ghost, at the ghost centre, one top-cell width (276 m) from the top-cell centre: `q_top = kappa_face (T_c - T_top) / ds_top`, with `kappa_face` the harmonic mean of the top-cell and ghost conductivities.

Formally that is two conditions on a subsonic outflow face, which admits one. It makes no measurable difference here because conduction resets the top-cell temperature far faster than the flow renews it. The top cell's thermal relaxation time `tau_th = rho c_v ds_top d / kappa_face` is 0.08 s, against a top-cell advection time `ds_top / V` of 7–29 s. So `T_top` is set by the conductive balance whatever the ghost carries, and at an outflow face the Riemann solution samples the interior entropy anyway. Against the previous decoupled closure (hydro ghost T from the top cell, conduction datum at the physical face) on the release mesh at 4000 s:

| | decoupled (previous) | pinned hydro T, face datum | **single T (release)** |
|---|---:|---:|---:|
| `V_top` [m/s] | 9.501 | 9.515 | **9.419** |
| `T_c - T_top` [K] | 69.7 | — | 139.4 |
| top-face `f_total` [kg m⁻² s⁻¹] | 2.702e−10 | — | 2.683e−10 |

Pinning the hydro ghost temperature is inert (+0.15 %). Nearly all of the −0.86 % comes from moving the conduction datum half a cell outward, a first-order offset in the top-cell width that converges away under refinement. **The equivalence needs conduction.** With `ISO_HEAT_FLUX=0` the fixed 22 kK ghost is inconsistent with the 22.67 kK initial top cell and nothing relaxes it: the 20 s conduction-off residual rises from 5.7e−3 to 1.4e−2 m/s and its refinement order falls from ~2 to ~1. With the ghost matched to the initial top cell it returns to the decoupled values to four digits. Full study: `docs/studies/boundaries/outer_temperature_bc_single_t_recap.md`.

## Which variable carries the outer condition {#config-outer-closure-swap}

Subsonic outflow at the 2153 km face admits exactly one incoming characteristic. This study was run under the previous decoupled closure, which put that one condition in the **pressure** (fixed reservoir back-pressure `kOuterPRef`) and left the hydro ghost temperature free. `ISO_OUTER_T_HYDRO_WALL` and `ISO_OUTER_P_NEUMANN` exist to put it in the **temperature** instead, and both outcomes are worth knowing before anyone tries it again. The current release fixes both the pressure and the ghost temperature; the temperature is inert on the release mesh (see [Outer temperature](#config-outer-temperature)), and the pressure is the condition that cannot be dropped.

| | decoupled (previous release) | `T_HYDRO_WALL=1`, `P_NEUMANN=1` (`dp/ds = 0`) | `T_HYDRO_WALL=1`, `P_NEUMANN=2` (`dp/ds = -rho g`) |
|---|---:|---:|---:|
| column mass at 4000 s | 99.994 % | **49.94 %** | 98.97 % |
| `V_top` [m/s] | +9.50 | **−2462.5** (= the `ISO_VCAP` cap) | +0.49 |
| top-face `f_total` [kg m⁻² s⁻¹] | +2.70e−10 | **−4.13e−07** | +1.47e−11 |
| `T` at 1900 km [K] | 6654 | **23694** (fully ionized) | 6670 |
| `p_top` drift over 4000 s | ~0 % | +493 % | +5.8 %, still rising |
| termination | `end_time` | `end_time` | `end_time` |

`dp/ds = 0` **destroys the model.** The interior needs `dp/ds = -rho g`, so a zero-gradient face permanently withholds the top cell's hydrostatic support; the column collapses inward, the transition region marches down past 1900 km within 200 s, half the column mass is pushed out through the base, and by 1000 s the run is a perfectly steady 23.7 kK fully ionized downdraft whose velocity is set by `ISO_VCAP` and nothing else. It never crashes and it reports `godunov.fallbacks=0`: **a clean termination is not evidence that a boundary closure is sound.**

`dp/ds = -rho g` is hydrostatically consistent and survives 4000 s, but with the pressure *level* unpinned it drifts secularly and the conduction-driven evaporation signal — the thing the model exists to measure — falls by a factor of ~19. Stability was not the reason to prefer the fixed back-pressure; measurability is.

That last point was then tested where it matters, on a wall-temperature comparison at 22 kK and 100 kK (both CFL 0.25; 100 kK fails the conduction Newton at 0.50 under either closure):

| quasi-steady, 4000 s | release 22 kK | release 100 kK | unanchored 22 kK | unanchored 100 kK |
|---|---:|---:|---:|---:|
| `V_top` [m/s] | +9.645 | +791.9 | +0.542 | +37.7 |
| `f_top / f_base` | +0.997 | +1.086 | −0.067 | −0.144 |
| column mass | 99.994 % | 99.887 % | 98.978 % | **77.3 %** |
| `p_top(t)/p_top(0)` | 1.000 | 1.000 | 1.057 | **2.276, rising** |
| `T` at 1900 km [K] | 6654 | 6546 | 6670 | **21931** |

The two-point exponent of `V ~ T_wall^n` is 2.911 anchored and 2.802 unanchored, so the **shape** of the wall response survives while every absolute value is ~20x too small (17.8x at 22 kK, 21.0x at 100 kK). Meanwhile the top and base mass fluxes take opposite signs — the column drains downward with an evaporative trickle on top — and at 100 kK the transition region migrates down to 1900 km and destroys the chromosphere there. A closure that reproduces the right power law and the wrong everything else is the worst kind to leave available without this warning.

Both knobs are diagnostic only and no release number may be quoted from them. Full study, including the controls that place the failure in the missing hydrostatic slope rather than the missing anchor: `docs/studies/boundaries/outer_bc_pressure_temperature_swap_recap.md`.
