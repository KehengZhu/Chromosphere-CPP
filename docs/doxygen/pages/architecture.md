# Source tree and architecture {#architecture}

This page describes where code lives, what depends on what, who owns the state, how a solver gets selected, and what the build produces.

## The dependency rule

There are two solvers over one shared mesh. `chromosphere.hpp` (the Grid, the two conserved-row index sets, the scratch and diagnostic-capture structures) and `eos.hpp` (the pure-hydrogen Saha closure and the CRASH `Gamma1` table) form the shared layer. Both solver directories depend on that layer and on nothing else of each other's.

@verbatim
                 chromosphere.hpp  +  eos.hpp
                    (Grid)            (closure)
                        ^                ^
                        |                |
            +-----------+                +-----------+
            |                                        |
   src/single_fluid/                         src/two_fluid/
   RELEASE                                   HISTORICAL research
   U = (rho, rho u, E)   -- no edge -->       7-row carrier state
   3 conserved rows        <-- no edge --     ion/neutral/electron
@endverbatim

`src/single_fluid/mixture.hpp` does not include `two_fluid/two_fluid.hpp`, and `two_fluid/two_fluid.hpp` does not include `single_fluid/mixture.hpp`. The two do not share a conserved-state width, a packing convention, a timestep function, or a boundary-ghost buffer. The only translation unit that includes both headers is `chromo_main.cpp`, which selects between them at run time. See @ref shared, @ref release_solver and @ref two_fluid_solver.

## Annotated file tree

### Shared headers (repository root)

| Path | What it holds |
| --- | --- |
| `chromosphere.hpp` | chromosphere::Grid — the object that owns all solver state — plus the storage scalar `Real` and its vector `Vec` (see Storage precision below), `mix::` (3 release rows), `cons::`/`prim::` (7 two-fluid rows), MixtureField, MixtureFaceArrays, MixtureRhsScratch, MixtureConductionScratch, MixtureFaceFluxCapture, OuterConductionCapture. See @ref grid. |
| `eos.hpp` | Pure-hydrogen Saha equilibrium closure, `EosGammaTable` loader, MixtureThermo / CaloricMixtureThermo, the caloric inversions and equilibrium face builders. See @ref eos. |
| `physics.hpp` | Inline closed-form physics: collision frequencies, `physical_kappa_e` / `physical_kappa_n` / `physical_conductivity`, radiative losses, heating rates, TRAC, hydrogen ionization/recombination rate coefficients. |
| `parallel.hpp` | Thin OpenMP abstraction: thread queries, per-thread exception collector, `parallel_for_cells`. Degrades to serial loops without `CHROMO_ENABLE_OPENMP`; rethrows the exception from the lowest physical cell index so a parallel failure matches a serial one. |
| `profiling.hpp` | Opt-in region timers, EOS-inversion and conduction-iteration statistics, active timestep limiter. Inert unless `CHROMO_PROFILE=1`. |
| `run_control.hpp` | Driver run-control helpers: `select_step_cap` (`CHROMO_STEP_CAP` / `CHROMO_T_END`), `OutputSchedule` (stride or `CHROMO_FRAME_DT`), `classify_termination`. Header-only, free of Armadillo and OpenMP, unit-testable without running a simulation. |
| `chromo_main.cpp` | The command-line driver. Owns no physics. See @ref driver and @ref io_formats. |

### Shared implementation, `src/*.cpp`

| Path | What it holds |
| --- | --- |
| `src/grid.cpp` | `Grid::init`, `Grid::resize`, `Grid::broadcast` / `force_rebuild_metrics`, and the fingerprinted static-mesh metric cache (`ds_i`, `B_i`, `B_imh`, `B_iph`, `dinvB_ds_i`). |
| `src/eos.cpp` | Saha closure, `Gamma1` table loader, caloric inversions and equilibrium face builders. Everything in double precision regardless of the storage type, log-domain where the state is known only through logarithms; safeguarded Newton with bisection fallback. |
| `src/profiling.cpp` | Profile counters in cache-line-padded per-thread slots, summed only when the report prints. |

### Release solver, `src/single_fluid/`

| Path | What it holds |
| --- | --- |
| `src/single_fluid/mixture.hpp` | Release API and governing equations: `mixture_decode` / `mixture_decode_into`, `mixture_rhs_explicit`, `mixture_timestep`, `mixture_advance`, `mixture_cell_phi`. |
| `src/single_fluid/mixture.cpp` | State decode, MUSCL reconstruction of `(ln rho, u, ln p)` with the MC3/Koren limiter, the release SWMF exact-Riemann Godunov flux with face-local Rusanov fallback, the Roe and Rusanov reference fluxes, the flux-tube pressure + gravity source, and the flux-consistent CFL timestep. |
| `src/single_fluid/exact_rs.hpp` / `.cpp` | Line-for-line port of the SWMF exact ideal-gas Riemann solver (`share/Library/src/ModExactRS.f90`): the star-state pressure/velocity iteration and the self-similar sampler the release flux calls at every face. |
| `src/single_fluid/integrator.cpp` | Implicit backward-Euler physical conduction, `mixture_conduction_residual_max`, and `mixture_advance` — the whole two-stage release timestep. |

### Historical two-fluid solver, `src/two_fluid/`

| Path | What it holds |
| --- | --- |
| `src/two_fluid/two_fluid.hpp` | Seven-row research-solver API. Every entry point rejects a Grid that carries a `Gamma1` table. |
| `src/two_fluid/state.cpp` | Packed-state accessors, index shifts (`ip1`/`im1`/`ip2`/`im2`), `cons2prim` / `prim2cons`. |
| `src/two_fluid/flux.cpp` | Cell-centred flux, spectral radius, pressure-area + gravity source. |
| `src/two_fluid/rhs.cpp` | Explicit TVD-MUSCL + Rusanov RHS and the implicit drag / frictional-heating / conduction right-hand sides. |
| `src/two_fluid/integrators.cpp` | `cal_dt_i`, `advance_Euler_state`, `advance_Euler_explicit_state`, `advance_RK4`, and the operator-split source Stages A–E plus beam, coronal heating, cooling and floors. |

### Scenarios, `scenarios/`

| Path | What it holds |
| --- | --- |
| `scenarios/scenario.hpp` | chromosphere::Scenario (`peek_ns`, `ic`, `update_bc`, `name`), `make_scenario`, and the shared open-hyperbolic `apply_open_bcs`. |
| `scenarios/scenario.cpp` | Name dispatch and the per-scenario override-preserving environment defaults, including the unconditional `GAMMA_TABLE` default that makes `model_column` the release scenario. |
| `scenarios/mesh.hpp` | Shared STATIC local-refinement face-grid construction (not AMR): generated once at IC time, never regridded. |
| `scenarios/mesh.cpp` | Face-grid construction and its `ISO_REFINE_*` / `GRID_REFINE_*` environment parsing. |
| `scenarios/model_column.hpp` | `model_column` (RELEASE) and `model_gentle` (historical two-fluid): one IC/BC implementation, two scenarios that never cross solvers. |
| `scenarios/model_column.cpp` | The column IC and boundary closures, including the rejection of every historical two-fluid knob on a release run. |
| `scenarios/model_c7.hpp` | Model C7 quiet-Sun chromosphere scenario and the shared C7 atmosphere library. |
| `scenarios/model_c7.cpp` | The C7 table, `c7_full_profile`, the IC and the `model_c7` boundary closures. |
| `scenarios/model_flare.hpp` | Explosive chromospheric evaporation scenario (Fisher, Canfield & McClymont 1985). |
| `scenarios/model_flare.cpp` | Reuses the C7 column and adds the transient nonthermal-electron beam and its `FLARE_*` controls. |
| `scenarios/analytic_canopy.hpp` | Exponential magnetic-canopy `B(z)` overlay on the C7 thermodynamic profile. |
| `scenarios/analytic_canopy.cpp` | The canopy recipe and its IC/BC. |
| `scenarios/data_file_parser.hpp` | Parser interface for the tabulated `[META]` / `[CELLS]` / `[GHOSTS]` scenario data-file format. |
| `scenarios/data_file_parser.cpp` | The parser. `[CELLS]` rows are `i ds_m B_imh_T B_iph_T phi_g_imh phi_g_iph ne_im3 nn_im3 T_K`; `[GHOSTS]` is exactly four rows tagged `outer_0`, `outer_1`, `inner_0`, `inner_1`. |
| `scenarios/pfss_field_line.hpp` | PFSS field-line scenario driven by a tabulated data file from `util/extract_field_line.py`. |
| `scenarios/pfss_field_line.cpp` | Grid population from the parsed field-line table and its boundary update. |
| `scenarios/data/` | Tracked field-line and loop input tables (PFSS lines, event ensembles, closed loops, SOL2014). Excluded from Doxygen. |

### Everything else

| Path | What it holds |
| --- | --- |
| `tests/chromo_tests.cpp` | The C++ test binary. Physical correctness against the write-up (pressure relation, sound speeds, spectral radius, conductivities, collision rate) plus numerical regressions. Each test constructs its own Grid. |
| `tests/test_conductive_flux.py` | Python regression tests for the shared conductive-flux diagnostic; imports `util/`. |
| `include/armadillo`, `include/armadillo_bits/` | Vendored header-only Armadillo, built with `ARMA_DONT_USE_LAPACK` / `ARMA_DONT_USE_BLAS`. Excluded from Doxygen; not ours. |
| `data/eos/` | Versioned EOS runtime data: `gamma1_hydrogen_v1.dat` (the production `Gamma1` table) and its committed `.sha256`. |
| `util/` | Python analysis and plotting scripts, the `util/eos/` table generator and validators, `util/extract_field_line.py`, and the venv setup. Figures default to `visualization/`. |
| `fortran/` | Optional Fortran/SWMF coupling shim: `chromo_interface.cpp` (C interface over the Grid and solver), `ModChromoParams.f90`, `ModChromosphereCPP.f90`, `ModChromosphereTest.f90`, `chromo_main.f90`. Not part of the CMake build. |
| `docs/` | Engineering recaps and plans, `supporting-papers/`, this Doxygen tree, and the separate Overleaf write-up worktree `docs/writeup-overleaf/`. |
| `outputs/` | Code outputs: snapshots, sidecars, console logs, `output.log`. Organized into per-scenario subdirectories. |
| `visualization/` | Generated figures and movies, per-scenario subdirectories. Gitignored. |
| `scripts/` | `run_chromo_omp.sh` (generic OpenMP launcher), `run_chromo_realtime.sh` (long production runs), `release_validation.sh` (release checks). |
| `cmake/` | CTest driver scripts: `check_eos_checksum.cmake`, `check_gamma_output.cmake`, `check_release_guard.cmake`. |

## State ownership

All solver state lives in one chromosphere::Grid struct that the **caller** owns. There is no global mutable solver state, no file-scope solver singleton, and no function-static cache in a numerical path.

- Solver functions take `const Grid&` when they only read the mesh, the prescribed field, gravity and the runtime toggles.
- They take `Grid&` only when they write a Grid-owned scratch or diagnostic buffer (`MixtureRhsScratch`, `MixtureConductionScratch`, `MixtureFaceFluxCapture`, `OuterConductionCapture`, the broadcast metric caches, the `eos_T_hint` accelerator).
- Scratch is scratch: `MixtureConductionScratch` is explicitly not a thermodynamic cache, and `eos_T_hint` only seeds the Newton inversion — a stale or absent hint changes no converged root, only the iteration count.
- The conserved state itself is a separate `Vec` the caller passes in and receives back; `MixtureField` binds a decode to one immutable `(Grid, Vec)` pair and `MixtureField::require_matches` throws if either the state buffer or a boundary ghost changed underneath it.

Multiple Grids can therefore coexist in one process, and `tests/chromo_tests.cpp` constructs a fresh Grid per test rather than sharing one.

## Solver selection

Selection is by **data**, not by a flag: a non-empty `Grid::eos_gamma_table` selects the release single-fluid solver; an empty one selects the historical two-fluid solver.

The chain is:

1. `make_scenario("model_column", ...)` in `scenarios/scenario.cpp` calls `set_env_default("GAMMA_TABLE", "data/eos/gamma1_hydrogen_v1.dat")` unconditionally, and rejects `ISO_GAMMA` outright. `model_gentle` and the other two-fluid scenarios never set it.
2. `chromo_main.cpp` reads `GAMMA_TABLE`, hashes the file and assigns `grid.eos_gamma_table = EosGammaTable::load(path)`. `release_mode` is then simply `!grid.eos_gamma_table.empty()`.
3. Both solvers guard their own entry points against the other's Grid, so a mismatch throws rather than silently computing nonsense:

| Guard site | Condition | Result |
| --- | --- | --- |
| `mixture_decode_into` (`src/single_fluid/mixture.cpp`) | `eos_gamma_table.empty()` | throws `"the mixture solver requires a Gamma1 table"` |
| `mixture_advance` (`src/single_fluid/integrator.cpp`) | `eos_gamma_table.empty()` | throws `"the mixture solver requires a Gamma1 table"` |
| `mixture_conduction_residual_max` (`src/single_fluid/integrator.cpp`) | `eos_gamma_table.empty()` | throws `"mixture conduction residual requires a Gamma1 table"` |
| `rhs_explicit_state` (`src/two_fluid/rhs.cpp`) | `!eos_gamma_table.empty()` | throws; points at the release solver in `mixture.hpp` |
| `cal_dt_i` (`src/two_fluid/integrators.cpp`) | `!eos_gamma_table.empty()` | throws; points at `mixture_timestep` |
| `advance_Euler_state` (`src/two_fluid/integrators.cpp`) | `!eos_gamma_table.empty()` | throws; points at `mixture_advance` |
| `advance_Euler_explicit_state` (`src/two_fluid/integrators.cpp`) | `!eos_gamma_table.empty()` | throws |
| `advance_RK4` (`src/two_fluid/integrators.cpp`) | `!eos_gamma_table.empty()` | throws |

The driver adds its own release-mode preconditions before the first step: it rejects `explicit` mode, any scenario other than `model_column`, a simultaneously supplied `ISO_GAMMA`, finite-rate `ionization`, and the legacy variables `SINGLE_FLUID`, `ENABLE_TE`, `ISO_TWO_FLUID` and `ISO_IONIZATION`. It also checks that the scenario returned exactly `grid.n_mixture_state` elements.

## Lifecycle of one run

@verbatim
  make_scenario(name, data_path)          scenarios/scenario.cpp
        |
  Grid::init(sc.peek_ns(), cfl)           src/grid.cpp  -- allocate ns-sized fields
        |
  [GAMMA_TABLE loaded?] -> release_mode   chromo_main.cpp
        |
  xn = sc.ic(grid)                        scenario IC: geometry, B, gravity,
        |                                 ghost buffers, initial conserved state;
        |                                 may call Grid::resize(ns_new) for a
        |                                 locally refined mesh
        |
  +---- per step ------------------------------------------------+
  |  sc.update_bc(grid, xn)               refresh ghost buffers   |
  |  dt = mixture_timestep(...)  |  cal_dt_i(...)                 |
  |  xn = mixture_advance(...)   |  advance_Euler_state(...)      |
  |  write_frame / sidecar records when the schedule says so      |
  +---------------------------------------------------------------+
        |
  final frame + termination report        termination=end_time|step_cap|other
@endverbatim

`Grid::resize` is a pure reallocation of the `ns`-sized fields: physical constants, CFL, ghost buffers and every runtime toggle are deliberately preserved, so a scenario can size a coarse grid, set its flags, then resize to the refined cell count. Metric caches are invalidated (`metrics_valid = false`) and rebuilt by the next `Grid::broadcast`, which fingerprints the five geometry arrays and skips the rebuild when none changed. See @ref driver for the driver's own control flow and @ref scenarios for what each IC/BC pair models.

## Build layout

One static library, five executables, and a CTest suite.

| Target | Sources | Purpose |
| --- | --- | --- |
| `chromosphere` (STATIC) | `src/*.cpp`, `src/single_fluid/*.cpp`, `src/two_fluid/*.cpp`, `scenarios/*.cpp` | Both solvers, the shared layer, and all scenarios. `PUBLIC` include dirs: repository root, `src`, `include`. |
| `chromo_main` | `chromo_main.cpp` | The simulation driver. |
| `chromo_tests` | `tests/chromo_tests.cpp` | The C++ test binary. |
| `eos_validate_saha` | `util/eos/validate_saha.cpp` | Offline cross-check of every CRASH `zAv` row against the actual double/log-domain C++ Saha implementation. Built but not registered as a CTest. |
| `eos_validate_production_table` | `util/eos/validate_production_table.cpp` | Validates the tracked production `Gamma1` table against the v1 loader contract. |
| `eos_validate_refined_table` | `util/eos/validate_refined_table.cpp` | Same check for a refined table. |

C++14, `-O3` in the default `Release` build type, `CMAKE_EXPORT_COMPILE_COMMANDS ON`.

### Registered CTest cases

| Test name | What it runs |
| --- | --- |
| `chromo_tests` | The `chromo_tests` executable. |
| `eos_production_table` | `eos_validate_production_table data/eos/gamma1_hydrogen_v1.dat`. |
| `eos_production_table_checksum` | `cmake/check_eos_checksum.cmake` — CMake's built-in SHA-256 against the committed `.sha256`, no external checksum program. |
| `release_reject_enable_te` | `cmake/check_release_guard.cmake` with `GUARD=enable_te`. |
| `release_reject_single_fluid` | Same script, `GUARD=single_fluid`. |
| `release_reject_two_fluid` | Same script, `GUARD=two_fluid`. |
| `release_reject_ionization` | Same script, `GUARD=ionization`. |
| `release_output_sidecar` | `cmake/check_gamma_output.cmake` — a release run must produce the snapshot and its `.gamma_diag` sidecar. See @ref io_formats. |
| `conductive_flux_python` | `tests/test_conductive_flux.py` under `.venv/bin/python`. Registered **only** when that interpreter exists. |

The four `release_reject_*` cases are generated by a `foreach` over `enable_te single_fluid two_fluid ionization`.

### Storage precision

Precision is **one knob**. `chromosphere::Real` is the scalar storage type and `chromosphere::Vec` is `arma::Col<Real>`, and every stored physical scalar is declared in terms of them: the conserved state, the static mesh metric arrays, the `Grid` physical constants (`m_i`, `m_n`, `m_e`, `k_b`, `q_e`, `chi_H_J`, `g`, `mu_0`), the run controls, `CFL`, `limiter_beta` and `Grid::sim_time`. **Nothing in the release path names `float` or `double` directly.** That coherence is the point: a single narrowing cast anywhere in the update chain reimposes the coarser representation on the whole chain no matter what `Vec` holds, which is exactly how the pre-cutover float32 storage survived an earlier `Vec`-only promotion — `mixture_pack_ghost`, the MUSCL-Hancock predictor row, the RHS row, the timestep row, the packed-energy round in the integrator and all the `Grid` constants still cast through `float`.

**The release storage type is `double`.** It is the default and the only supported production configuration, for two independent reasons.

*Measured necessity.* The column resolves a stratified layer whose mass-loading timescale exceeds the CFL step by `tau/dt ~ 3.5e7`. One relative ULP of `float` is `6e-8`, so `1/u = 1.7e7` — smaller than that ratio, and the per-step continuity increment therefore falls below half a ULP of the stored `rho`. In the 4000 s `N = 500` gentle-evaporation run `dt*(-dF/ds)` is about `2.6e-17 kg m^-3 s^-1` against a float32 ULP of `rho ~ 9e-10` of about `5.6e-17`, and **619 of the 660 cells sat below half a ULP**, so the slow density adjustment could not accumulate at all. For `double`, `1/u = 9.0e15` leaves about eight orders of magnitude of margin.

*Consistency with the target coupling framework.* Every shipped SWMF/BATS-R-US build sets `PRECISION = ${DOUBLEPREC}` in `share/build/Makefile.<OS>.<compiler>` (all 23 templates, the nvfortran/GPU ones included), `Config.pl` defaults to double, and BATS-R-US declares `State_VGB`, `StateOld_VGB`, `Source_VC` and the fluxes as default `real` — so `U^{n+1} = U^n + dt*R` is a float64 accumulation there, with no compensated summation and no deviation formulation. A double-precision conserved state is the SWMF-consistent choice; the float32 state was the inconsistent one.

`CHROMO_STATE_FLOAT32` (option, default `OFF`) reverts to the pre-cutover `float` storage. It is **DIAGNOSTIC ONLY** — kept so the round-off sensitivity stays reproducible — is not exercised by `scripts/release_validation.sh`, and no release number may be quoted from it. Build it into a separate tree so the release build is never shadowed:

```bash
cmake -S . -B build_omp_f32 -DCMAKE_BUILD_TYPE=Release \
      -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_STATE_FLOAT32=ON
```

`CHROMO_STATE_FLOAT64` was **retired**: double is the default, so CMake hard-errors on the old flag rather than producing a build that looks configured and is not. There is no `build_omp_f64` tree any more.

Both precisions pass the whole unit suite with no loosened tolerance. Evidence: `docs/studies/numerics/state_precision_release_cutover.md` for the cutover, `docs/studies/numerics/float32_precision_control_experiment.md` for the original diagnostic that motivated it, and Known limitations in @ref validation for what the cutover did and did not fix.

### OpenMP

`CHROMO_ENABLE_OPENMP` (option, default `OFF`) links OpenMP into the `chromosphere` library `PUBLIC`ly and defines `CHROMO_USE_OPENMP=1`. It first tries `find_package(OpenMP)`; on Apple Clang it falls back to locating `libomp` under `LIBOMP_ROOT`, `/opt/homebrew/opt/libomp` or `/usr/local/opt/libomp` and adds `-Xpreprocessor -fopenmp`. If neither path works the configure step fails hard rather than silently producing a serial build. With the option `OFF`, `parallel.hpp` compiles to plain serial loops and nothing else in the solver branches on it. The thread count is never hard-coded in C++ or CMake; it comes from `OMP_NUM_THREADS` at run time.
