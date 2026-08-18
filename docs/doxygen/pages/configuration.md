# Configuration reference {#configuration}

Everything a run is configured with lives in exactly two places: the seven positional command-line arguments of `chromo_main`, and the environment. There is no input deck, no configuration file, and no runtime parser. The authoritative sources for this page are `chromo_main.cpp`, `run_control.hpp`, `scenarios/scenario.cpp`, `scenarios/model_column.cpp`, `scenarios/mesh.cpp`, `scenarios/model_c7.cpp`, `scenarios/model_flare.cpp` and `scenarios/pfss_field_line.cpp`; no environment variable is read anywhere in `src/`.

## Command line

Usage: `chromo_main [output_path] [mode] [ionization] [scenario] [data_path] [time_mult] [cooling]`. All seven are positional and optional; to set a later one you must supply the earlier ones (pass `-` or `""` for an unused `data_path`).

| # | Argument | Accepted values | Default | Meaning |
| --- | --- | --- | --- | --- |
| 1 | `output_path` | any path | `outputs/output.txt` | Snapshot file. Sidecars are derived from it: `<out>.gamma_diag`, `<out>.faceflux`, `<out>.outercond`. See @ref io_formats. A one-line summary is always appended to `outputs/output.log`. |
| 2 | `mode` | `full`, `explicit` | `full` | `explicit` selects the two-fluid explicit-only path (`advance_Euler_explicit_state`, R_I identically zero). Any string other than `explicit` means `full`. Combining `explicit` with a loaded `Gamma1` table is a hard error: "the release solver has no explicit-only mode". |
| 3 | `ionization` | `ionization`, `no-ionization` | `no-ionization` for `model_column`, `ionization` for every other scenario | Sets `Grid::enable_ionization` (finite-rate Stage-E network). Any string other than the literal `no-ionization` counts as on. In release mode `ionization` is a hard error: the release closure is Saha equilibrium. |
| 4 | `scenario` | `model_column`, `model_gentle`, `model_c7`, `model_flare`, `analytic_canopy`, `pfss_field_line` | `model_column` | Selects the IC/BC pair via `make_scenario`. An unknown name is a hard error. See @ref scenario_reference. |
| 5 | `data_path` | path to a field-line `.dat` | `""` | Required by `pfss_field_line` (hard error if empty). Ignored by every other scenario. |
| 6 | `time_mult` | float | `1.0` | Multiplier on the derived total time `10 * sum(ds_i) / 2e4 m/s`. Also drives the legacy step cap `max(10000, 10000*time_mult)` and the legacy frame stride `max(10, 10*max(1, time_mult/5))`. Both legacy rules are superseded when `CHROMO_T_END` / `CHROMO_STEP_CAP` / `CHROMO_FRAME_DT` are set. |
| 7 | `cooling` | `cooling`, `no-cooling` | `no-cooling` for `model_column`, `cooling` for every other scenario | Sets `Grid::enable_radiative_cooling`. Any string other than the literal `no-cooling` counts as on. |

The differing defaults for arguments 3 and 7 are decided by `scenario_name == "model_column"` alone, evaluated before either argument is read, so `chromo_main out.txt` and `chromo_main out.txt full` both run the release column with ionization and cooling off.

Note on argument 7 in release mode: `model_column_ic` recomputes `Grid::enable_radiative_cooling` from `ISO_HEAT_FLUX && ISO_COOLING`, and `ISO_COOLING` is rejected outright in release mode, so passing `cooling` to `model_column` is silently overridden to off rather than rejected. The release timestep contains no radiative-cooling stage in any case — nothing in `src/single_fluid/` reads `enable_radiative_cooling`, `enable_ionization`, `enable_trac` or `enable_beam_heating`.

## Environment variables

Two parsing conventions are used and they behave differently.

- Driver booleans (`CHROMO_*` in `chromo_main.cpp`) use `env_enabled()`: unset takes the stated default; the values `0`, `false` and `off` disable; **every other value enables, including an empty string**.
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
| `GAMMA_TABLE` | path to a `Gamma1` table | supplied unconditionally by `model_column` as `data/eos/gamma1_hydrogen_v1.dat` (only the path is overridable) | **A loaded table is what selects the release single-fluid solver.** The driver SHA-256s the file and records the hash in the `gamma_diag` header. Setting it together with `ISO_GAMMA` is a hard error; setting it with `mode=explicit` is a hard error; setting it with any scenario other than `model_column` is a hard error ("the release solver currently supports only model_column"). Without a table the driver runs the historical two-fluid solver. |
| `SINGLE_FLUID` | float, non-zero = on | unset | Two-fluid research knob: slaves neutrals to the ion fluid. Read only when no `Gamma1` table is loaded. **Rejected (exit 1) in release mode if present at all.** |
| `ENABLE_TE` | float, non-zero = on | unset | Two-fluid research knob: separate electron temperature `T_e != T_i`. Read only when no `Gamma1` table is loaded. **Rejected (exit 1) in release mode if present at all.** |

### `model_column` / `model_gentle` model knobs (the `ISO_*` family)

`model_column` and `model_gentle` share one IC/BC implementation (`scenarios/model_column.cpp`) but never share a solver. Both scenarios install override-preserving defaults with `setenv(..., overwrite=0)` before the IC runs, so an explicit assignment in the environment always wins.

`model_column` presets: `GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat`, `ISO_H_BASE=1600`, `ISO_DH=553`, `ISO_NS=500`, `ISO_HEAT_FLUX=1`, `ISO_T_TOP=22000`, `ISO_HYDRO_T_DECOUPLE=1`, `ISO_REFINE_PROFILE=outer`, `ISO_REFINE_FACTOR=4`, `ISO_REFINE_S_LO_KM=500`, `ISO_REFINE_TRANSITION_KM=20`.

`model_gentle` presets: `ISO_CORONA=1`, `ISO_H_BASE=1003`, `ISO_HEAT_FLUX=1`, `ISO_COOLING=1`, `ISO_IONIZATION=1`, `ISO_TWO_FLUID=1`, `ISO_NS=600`. It never loads a `Gamma1` table.

Knobs valid in **both** modes:

| Variable | Type | Default (bare, before presets) | Effect |
| --- | --- | --- | --- |
| `ISO_NS` | unsigned integer | `1000` if `ISO_CORONA` is set and non-empty and not `0`, else `600` | Coarse-equivalent cell count read by `peek_ns` before the IC. With refinement on this is the coarse-grid equivalent, not the allocated cell count. |
| `ISO_H_BASE` | float, km | `0` | Height of the lower truncation face above the photosphere. It is a truncation of the modeled domain, not the photosphere. Also sets `Grid::out_base_km`, the height label offset in the output file. |
| `ISO_DH` | float, km | `1303` (or `ISO_CORONA_TOP_KM - ISO_H_BASE` when `ISO_CORONA` is on) | Domain span above the base. |
| `ISO_HEAT_FLUX` | float, non-zero = on | `0` | `0` = Stage-1 adiabatic relaxation (conduction off); non-zero = conduction on. |
| `ISO_RELAX_TIME` | float, seconds | `0` | Length of an initial Stage-1 adiabatic relaxation phase before conduction and the top jump switch on. |
| `ISO_T_TOP` | float, K | `0` (= use the IC top-cell temperature) | Any positive value is an **absolute** imposed outer-boundary temperature and is not re-anchored across a relax phase. Release preset 22000. The column is *very* sensitive to it: the quasi-steady top-of-domain velocity goes as roughly `T_wall^3` and reverses sign near 10 kK — see [Sensitivity to the wall temperature](#config-twall-sensitivity). |
| `ISO_TJUMP_A` | float | `1.0` | Outer ghost temperature jump `T_ghost1 = a * T_ref`. |
| `ISO_TJUMP_B` | float | `1.0` | Second ghost jump `T_ghost2 = b * T_ghost1`. |
| `ISO_HYDRO_T_DECOUPLE` | float, non-zero = on | `0` | Lets the outer **hydro** ghost temperature zero-gradient-extrapolate the live top cell while the conduction solver keeps the fixed wall at the physical outer face. Release preset 1. |
| `ISO_VCAP` | float, Mach | `0.1` | Outer outflow velocity cap as a fraction of the sound speed. |
| `ISO_INNER_T_NEUMANN` | float, non-zero = on | `0` (Dirichlet) | Conduction lower BC: Dirichlet fixed reservoir temperature vs Neumann (insulating). |

Knobs read only by the **two-fluid** path (`model_gentle`); each is a hard error in release mode unless noted:

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
| `ISO_RECONSTRUCTION` | `lnrho-v-lnt`, `lnrho-v-lnp` | `lnrho-v-lnp` | `lnrho-v-lnt` restores the temperature-based primitive set, which closes pressure through the nonlinear EOS and manufactures face-pressure mismatch across the partial-ionization transition. `lnrho-v-lnp` requires Gamma/Saha mode. Any other value throws. |
| `ISO_LIMITER` | `mc3`, `minmod`, `first` | `mc3` with beta = 2 | `minmod` selects the symmetric minmod limiter; `first` sets beta = 0 with MC3, giving phi identically 0, i.e. piecewise-constant first order. Any other value throws. |
| `ISO_MC3_BETA` | float | `2.0` | Read **only** when `ISO_LIMITER=mc3` is explicitly set. Setting `ISO_MC3_BETA` alone has no effect. |

The numeric-method toggles that used to select the "bestwb" path — `ISO_LOG_RECON`, `ISO_MC3`, `ISO_EQ_WB`, `ISO_INNER_WB`, `ISO_INNER_WB_RHO`, `ISO_C7_IC`, `ISO_C7_HSE` — are now **unconditional** in `model_column` / `model_gentle` and are simply ignored if set, so old command lines still run. Well-balanced reconstruction, log-space MUSCL, the MC3/Koren limiter at beta = 2, the C7 IC and its hydrostatic density rebuild are always on; the equilibrium-reference well-balancing that `ISO_EQ_WB` used to select was **retired** from the release, which is reference-free (`Grid::eq_wb` survives for the two-fluid research solver only and must never be set on a release Grid). The diagnostic toggles `ISO_DIAG_*`, the base sponge `ISO_SPONGE_*`, `ISO_TRAC_TC_FIXED` and the analytic-isentrope toy IC knobs (`ISO_T_BASE`, `ISO_N_BASE`, `ISO_F_ION`) were removed outright and are not read anywhere.

## Rejected in release mode

"Release mode" means a loaded `Gamma1` table, which in practice means the `model_column` scenario. The guards fire in three places and every one of them stops the run.

`make_scenario` (before anything is allocated), for `model_column` only:

- `ISO_GAMMA` present -> `std::runtime_error`: "model_column is the single-fluid release scenario and derives every thermodynamic index from the Saha/Gamma1 closure; ISO_GAMMA ... has no meaning here. Use model_gentle for the historical fixed-gamma column."

`chromo_main`, in the `GAMMA_TABLE` branch and immediately after (all print `release-solver error:` and exit 1):

- `ISO_GAMMA` present together with `GAMMA_TABLE` -> "GAMMA_TABLE and explicitly supplied ISO_GAMMA are mutually exclusive".
- `mode == "explicit"` -> "the release solver has no explicit-only mode".
- any scenario other than `model_column` -> "the release solver currently supports only model_column".
- the positional `ionization` argument set to anything but `no-ionization` -> "finite-rate ionization was requested, but the release closure is Saha equilibrium; pass no-ionization".
- `SINGLE_FLUID`, `ENABLE_TE`, `ISO_TWO_FLUID` or `ISO_IONIZATION` **present at all**, regardless of value -> "<name> is a legacy two-fluid setting and has no meaning for the single-fluid release solver".
- a scenario IC returning a state whose element count is not `grid.n_mixture_state` -> "the scenario returned a state of N elements, expected M (rho, rho u, E per cell)".

`model_column_ic`, for every knob in `kNonReleaseKnobs` **present at all**, regardless of value, throwing `std::runtime_error` "<name> is historical two-fluid research configuration and has no meaning for the single-fluid release solver; run the model_gentle scenario for that path": `ISO_GAMMA`, `ISO_TWO_FLUID`, `ISO_IONIZATION`, `ISO_COOLING`, `ISO_TRAC`, `ISO_CORONA`, `ISO_CHEAT`, `ISO_QFLUX`, `ISO_TBOOST`, `ISO_NUMERICAL_DIFFUSIVITY_MULT`.

Additionally, `CHROMO_FACE_FLUX_DIAG` and `CHROMO_OUTER_COND_DIAG` are rejected the other way round: both require the release solver and throw if enabled on a two-fluid run.

## Validated production defaults

`CHROMO_CFL=0.50` is the validated production runtime default for the **reduced `model_column` release configuration only**: coarse-equivalent N = 500, R4 outer refinement, 661 actual cells, the 22,000 K conductive boundary at the physical outer face, and physical conduction only. Evidence: `docs/studies/conduction/coarse_model_column_physical_conduction_recap.md`.

The hard-coded C++ solver default remains `CHROMO_CFL=0.25` (`chromo_main.cpp` and `Grid::CFL`), 0.25 stays the comparison reference, and **the C++ default must not be changed**. `scripts/run_chromo_realtime.sh` supplies 0.50 as an override-preserving default and only when the scenario argument is exactly `model_column`; every other scenario keeps 0.25 unless `CHROMO_CFL` is given explicitly. Rerun the CFL sweep before trusting 0.50 on different hardware or a materially different model shape.

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

`ISO_T_TOP` is the most consequential single knob in the release configuration, so it deserves a number rather than a warning. With everything else fixed at the release configuration (C7 initial condition, N = 500/R4, `ISO_HYDRO_T_DECOUPLE = 1`, physical conduction only, 4000 s, double-precision state) a one-parameter sweep of the wall over 8–100 K gives a quasi-steady top-of-domain velocity of

| `ISO_T_TOP` [kK] | 8 | 10 | 12 | 15 | 18 | **22** | 30 | 45 | 70 | 100 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `V_top` [m/s] | −0.44 | −0.10 | +0.54 | +2.1 | +4.6 | **+9.6** | +27.6 | +93.5 | +310 | +792 |

Four decades of velocity for one decade of wall temperature. The response is a smooth power law on the upflow branch — a fit over `ISO_T_TOP >= 15000` gives `V ~ T_wall^3.09`, with the local exponent falling from 3.7 at the 22 kK release point to 2.6 at 100 kK — and it reverses sign near 10 kK, below which the column drains instead of evaporating. There is no saturation anywhere in the tested range and the `ISO_VCAP` outflow cap never binds (at most 15 % of the cap, at 100 kK).

The mechanism is entirely conductive: the wall thermostats the top cell to within 0.3 %, Spitzer conductivity supplies the `T^2.5`, and a column with no radiative sink must export the resulting conductive flux as the enthalpy flux of the evaporating material, so `V ~ q_wall / [rho_top (h + x chi_H/m_H)]` to within a factor of two over the whole range. Two practical consequences: any quantitative evaporation result is a statement about the chosen wall temperature as much as about the chromosphere, and a wall below about 8 kK drives the top of the column below the `Gamma1` table temperature floor and aborts the run.

Full study, including the time-convergence and CFL caveats: `docs/studies/boundaries/upper_wall_temperature_sensitivity_recap.md`.
