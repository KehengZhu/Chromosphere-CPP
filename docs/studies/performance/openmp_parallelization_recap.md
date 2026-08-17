# OpenMP Parallelization Recap

## Final status

```text
OPENMP_PARALLELIZATION: GO
```

The active gamma-table chromosphere solver now has an explicit, deterministic CPU OpenMP path and a clean serial fallback. The physical model, numerical method, mesh, timestep sequence, limiter, Riemann solver, EOS, nonlinear conduction formulation, Thomas solver, boundary conditions, diagnostics, and fail-closed behavior were preserved.

The best stable result on the audited machine is:

```text
OMP_NUM_THREADS=16
20 physical seconds
2638 actual refined cells
28466 steps
36.346 s median wall time
1.277 ms/step
2.170x versus the same-session serial median
```

The 2.170x result is below the task's approximate 2.5x expectation, but it is material, stable, and achieved without any numerical or physical change. The remaining wall time is dominated by the conduction stage, especially the serial Thomas solve and the barriers required by the nonlinear Newton sequence.

---

## Repository and build state

### Initial repository state

The workspace was initially clean:

```text
## main...origin/main
```

No pre-existing uncommitted changes had to be preserved. Separate build directories were used throughout validation and removed after the final audit so that generated build products are not part of the source patch.

### Source files changed

```text
CMakeLists.txt
parallel.hpp
chromo_main.cpp
chromosphere.hpp
src/eos.cpp
src/grid.cpp
src/integrators.cpp
src/profiling.cpp
src/rhs.cpp
src/state.cpp
tests/chromo_tests.cpp
docs/studies/performance/openmp_parallelization_recap.md
```

No scenario, mesh, boundary-condition, EOS-table, physical-constant, limiter, Riemann-solver, or numerical-conduction parameter file was changed.

### Compiler, runtime, and CPU topology

```text
Compiler:       Apple clang 21.0.0 (clang-2100.1.1.101)
Target:         arm64-apple-darwin25.6.0
OpenMP runtime: Homebrew libomp 22.1.8
Runtime path:   /opt/homebrew/opt/libomp/lib/libomp.dylib
CPU:            16 physical/logical cores
Performance:    12 cores
Efficiency:      4 cores
```

The serial executable was verified not to link `libomp`. The OpenMP executable was verified to link:

```text
/opt/homebrew/opt/libomp/lib/libomp.dylib
```

### Explicit CMake configuration

OpenMP is now opt-in:

```cmake
option(CHROMO_ENABLE_OPENMP "Enable OpenMP parallel execution" OFF)
```

`OFF` has no OpenMP dependency. `ON` first uses standard `find_package(OpenMP)`. AppleClang additionally supports `LIBOMP_ROOT`, the `LIBOMP_ROOT` environment variable, and the conventional Homebrew roots:

```text
/opt/homebrew/opt/libomp
/usr/local/opt/libomp
```

If `CHROMO_ENABLE_OPENMP=ON` is requested and no working runtime is found, configuration fails instead of silently building serial code. This fail-closed behavior was exercised before `libomp` was installed.

Exact builds:

```bash
cmake -S . -B build_serial \
  -DCMAKE_BUILD_TYPE=Release \
  -DCHROMO_ENABLE_OPENMP=OFF
cmake --build build_serial -j4

cmake -S . -B build_omp \
  -DCMAKE_BUILD_TYPE=Release \
  -DCHROMO_ENABLE_OPENMP=ON \
  -DLIBOMP_ROOT=/opt/homebrew/opt/libomp
cmake --build build_omp -j4
```

The runtime banner reports whether OpenMP is compiled, the maximum requested team, the conduction team, and `OMP_NUM_THREADS` when present. `OMP_NUM_THREADS=1` and `OMP_NUM_THREADS=4` were independently verified to report different runtime thread availability, confirming that the OpenMP runtime was genuinely active.

---

## Thread-safety changes

### Temperature-hint lifecycle

`Grid::eos_T_hint` is now fully allocated and initialized to NaN in both `Grid::init` and `Grid::resize`.

`Grid::store_eos_temperature_hint` no longer performs lazy whole-vector allocation. It fails closed if the array was not pre-sized or the cell index is invalid. Parallel loops write only `eos_T_hint[i]` for their own physical cell.

The existing provenance rule remains unchanged: a caller only reads a hint when that call site already had the right cell-matched provenance.

### Scratch ownership

`GammaRhsScratch` and `GammaConductionScratch` remain shared members of one `Grid`, but all parallel writes are index-disjoint and every shared array is sized before entering a parallel region.

The supported ownership contract is explicit:

```text
One Grid instance may have one active top-level solve.
One Grid instance is not concurrently reentrant.
Nested OpenMP execution is disabled.
```

No per-cell vector is resized inside an OpenMP loop.

### Runtime profiling and EOS counters

Runtime profiling and `CHROMO_EOS_COUNTS` now use fixed-capacity, cache-line-aligned per-thread slots:

```text
maximum slots: 256
slot ownership: OpenMP thread ID
merge order:   slot 0, 1, 2, ...
```

Each thread updates only its own slot. Queries, reset operations, and final printing merge the slots deterministically. There are no hot-path atomic increments for Saha, Gamma1, inversion, limiter, or conduction counters.

Stage-level wall timers remain outside parallel cell loops and continue to report elapsed wall time rather than summed thread CPU time.

Serial builds use slot zero and retain the same public profiling interface.

### Deterministic exception propagation

`parallel.hpp` defines a reusable `ParallelFailure` collector. Every physical-cell worker catches exceptions locally and records:

```text
physical cell index
std::exception_ptr
```

After the parallel region, the caller deterministically selects and rethrows the exception associated with the lowest physical cell index. The successful path uses no atomics. Scalar/ghost failures are assigned an index after all physical cells, so a simultaneous physical-cell failure remains the deterministic winner.

A focused test injects failures at cells 11 and 3 and verifies that `cell=3` is rethrown.

### Legacy OpenMP cleanup

Two pre-existing `collapse(2)` pragmas in the small `ip1` and `im1` shift helpers became active once OpenMP was explicitly enabled. They were removed because:

* the loops are too small to amortize a parallel region;
* the task explicitly leaves tiny reconstruction/shift expressions serial;
* ThreadSanitizer first reported them while examining the new OpenMP build.

The helpers now remain serial in both build modes.

---

## Parallelized regions

### Physical-cell mixture decode

Functions:

```text
decode_mixture_field_into
decode_predicted_caloric_primitives_into
```

Parallel work:

* reads seven conserved rows, cell potential, and the index-matched prior temperature guess;
* performs the existing caloric/Gamma1 decode in the same cell-local order;
* writes `out.cells[i]` where applicable;
* writes the three physical primitive entries for cell `i`;
* writes `eos_T_hint[i]`.

Serial work:

* metadata and boundary-signature setup;
* all four ghost decodes.

Scheduling and failure behavior:

```text
schedule(static)
index-disjoint writes
deterministic lowest-cell exception rethrow
```

### Equilibrium projection

Function:

```text
project_equilibrium_single_fluid(const Grid&, const Vec&)
```

Each worker independently reads the seven conserved rows for one cell, runs the existing caloric inversion and row-only equilibrium packing, writes the seven projected rows for that cell, and updates only `eos_T_hint[i]`.

The row-only path from the final serial refinement is preserved. No Gamma1 query was reintroduced.

### Six MUSCL face-state constructions

The original six face bundles remain exactly:

```text
predictor: 2
corrector: 4
```

They are now built through an exact cell-level helper and two batched parallel loops:

```text
predictor batch: one parallel region, two bundles per cell
corrector batch: one parallel region, four bundles per cell
```

All destination arrays are sized before entering the region. For each face bundle and cell, the following order is unchanged:

1. read reconstructed `log(rho)`, velocity, and `log(T)`;
2. perform the log-aware EOS and Gamma1 face construction;
3. pack conserved rows;
4. calculate physical flux rows;
5. calculate spectral radius.

Gamma1 remains evaluated for every face state. Each worker writes only its own cell positions in each output bundle.

The limiter, stencil expressions, Rusanov combination, well-balanced residual, and source addition remain serial.

### Nonlinear conduction around the Thomas solve

The gamma-table conduction stage uses one persistent OpenMP region per timestep. The physical-cell loops cover:

1. initial caloric decode;
2. per-Newton fused caloric evaluation;
3. conductivity construction;
4. face coefficient and heat-capacity construction;
5. residual and tridiagonal matrix assembly;
6. temperature update;
7. final known-temperature row packing.

The following remain `single` or serial:

* ghost thermodynamic evaluation;
* scalar boundary quantities;
* deterministic max merges;
* convergence and stagnation decisions;
* `thomas_solve_double_inplace`;
* scalar outer-conduction diagnostic finalization.

The Newton matrix, lagged conductivity, residual normalization, `2e-11` tolerance, 40-update cap, stagnation behavior, boundary-row modifications, face-local numerical conduction, and final conservative energy remainder are unchanged.

Per-thread residual and step maxima are written into fixed slots, followed by an explicit barrier and a serial fixed-order max merge. All threads observe the same stop/convergence decision after the `single` barrier.

#### Conduction team cap

Short production-shaped probes showed that the many necessary Newton barriers make large conduction teams counterproductive at 2638 cells:

```text
2-thread conduction:  approximately 0.092 s per 0.1-physical-s probe
4-thread conduction:  approximately 0.087 s
8-thread conduction:  approximately 0.132 s
12-thread conduction: approximately 0.204 s before capping
```

The persistent conduction region therefore uses:

```text
min(4, OMP maximum threads)
```

The RHS, decode, projection, face construction, and final packing still use the full requested team. With a global 12-thread configuration, the capped conduction probe returned to approximately 0.092 s.

---

## Validation

### Test suites

Final source state:

```text
build_serial ctest:          8 / 8 PASS
build_omp ctest, 12 threads: 8 / 8 PASS
build_serial/chromo_tests:   9070 / 9070 PASS
build_omp/chromo_tests:      9070 / 9070 PASS
```

New focused checks cover:

* OpenMP compile-time ON/OFF behavior;
* runtime thread bounds;
* temperature-hint pre-sizing and fail-closed behavior;
* deterministic lowest-cell exception propagation;
* exact per-thread EOS-counter aggregation;
* exact per-thread inversion-profile aggregation.

Existing production tests exercise decode, projection, face construction, and conduction under the OpenMP build.

### Gate A: serial build versus OpenMP one thread

The complete 20-physical-second diagnostic run compared:

```text
build_serial
build_omp with OMP_NUM_THREADS=1
```

Results:

```text
steps                         28466 = 28466
timestep sequence             identical
active limiter                acoustic for every step
EOS inversion decisions       identical
Gamma1/EOS counters           identical
conduction Newton counts      identical
conserved output              byte-identical
.gamma_diag                   byte-identical
.faceflux                     byte-identical
.outercond                    byte-identical
```

### Gate B: serial build versus OpenMP twelve threads

The same full comparison passed with `OMP_NUM_THREADS=12`. All four output files were byte-identical to the serial build.

### Gate C: repeated multi-thread determinism

The final source state was run three times with:

```text
OMP_DYNAMIC=FALSE
OMP_MAX_ACTIVE_LEVELS=1
OMP_NUM_THREADS=12
```

Each run was byte-identical to serial and to the other two OpenMP runs. Final SHA-256 hashes:

```text
conserved output  47794ea5da5507261ba9ec55e33bc84ccb4e391ef44cc17335c8303d0f64cbc6
.gamma_diag       f2ff63533394d18d4121a4c1631360c7d720df800e5a17a2c2bf96db932a4821
.faceflux         31becf7190f42e6744e90909a38a055d5a6d8e2b5f654dfcd563c068918c0148
.outercond        dbddffbb04d66092d686646d49fc6ce21a400b264842ad1239a243c34da2a86b
```

Diagnostic counter identity in every repeated run:

```text
inversion.calls                         300995013
inversion.bisection_fallbacks             1699711
inversion.bracket_evaluations             14122580
conduction.calls                              28466
conduction.average_newton_updates          2.997681
conduction.maximum_newton_updates                 4
eos.gamma1_queries                       526262300
eos.temperature_logs                     879186665
eos.n_h_logs                            1052442140
```

Performance runs, where output and EOS diagnostic counters were disabled, retained the same values at every thread count:

```text
steps                                     28466
inversion.calls                       300694281
inversion.initial_guess_accepts        67214120
inversion.one_update_convergences     226719603
inversion.maximum_iterations                  16
inversion.bisection_fallbacks            459564
inversion.bracket_evaluations          13521116
conduction.calls                           28466
conduction.average_newton_updates       2.997681
conduction.maximum_newton_updates              4
```

### ThreadSanitizer result and platform limitation

An OpenMP+TSan build compiled successfully. AppleClang ThreadSanitizer then reported reads after standard OpenMP implicit barriers as races.

A minimal independent program containing only:

```cpp
#pragma omp parallel for
for (...) array[i] = i;
// serial read of array after the implicit barrier
```

produced the same report with AppleClang 21 and Homebrew libomp 22.1.8. Therefore this platform combination does not model libomp's barrier happens-before relation correctly and cannot provide a clean OpenMP TSan verdict.

This limitation is supplementary rather than a substitute for thread-safety review. Race closure is instead supported by:

* strict index-disjoint ownership;
* no shared resizing in regions;
* no hot atomics;
* explicit barriers before per-thread max merges;
* deterministic failure tests;
* exact counter identity;
* serial/one-thread/multi-thread byte identity;
* three repeated multi-thread hash-identical runs;
* removal of the two legacy small-loop pragmas found during the TSan investigation.

---

## Performance

### Protocol

Exact production case:

```text
20 physical seconds
2638 actual refined cells
28466 steps
output disabled
Gamma diagnostics disabled
EOS operation counters disabled
runtime profiling enabled
```

Every OpenMP thread count was run three times after the final source cleanup. The same-session serial build was run four times because a connector retry produced one additional complete run; its median is used below.

Historical audited serial median:

```text
76.072 s
```

Same-session serial median:

```text
78.856 s
```

The current machine was approximately 3.7% slower during this benchmark session, so speedups below use 78.856 s. Both baselines are retained for transparency.

### Scaling table

| threads | wall median s | run range s | ms/step | speedup | efficiency | RHS s | decode s | projection s | conduction s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| serial | 78.856 | 78.554–78.907 | 2.770 | 1.000x | 100.0% | 40.775 | 8.951 | 7.684 | 20.359 |
| 1 | 80.477 | 80.284–80.772 | 2.827 | 0.980x | 98.0% | 41.255 | 9.509 | 7.856 | 20.702 |
| 4 | 40.724 | 40.647–40.785 | 1.431 | 1.936x | 48.4% | 14.189 | 3.879 | 2.941 | 18.514 |
| 8 | 37.537 | 37.456–37.568 | 1.319 | 2.101x | 26.3% | 11.206 | 3.002 | 2.498 | 19.581 |
| 12 | 36.684 | 36.476–36.699 | 1.289 | 2.150x | 17.9% | 11.643 | 2.487 | 2.242 | 19.098 |
| **16** | **36.346** | **36.335–36.551** | **1.277** | **2.170x** | **13.6%** | **11.143** | **2.540** | **2.293** | **19.212** |

Against the historical 76.072-second baseline, the final 16-thread speedup is approximately 2.093x. Against the original pre-serial-optimization runtime of approximately 124.7 seconds, the combined serial+OpenMP improvement is approximately 3.43x; that comparison is context only and was not used as the OpenMP acceptance baseline.

### Interpretation

The largest gains occur in genuinely independent per-cell work:

```text
RHS:        about 40.8 s -> 11.1 s
decode:      about 9.0 s ->  2.5 s
projection:  about 7.7 s ->  2.3 s
```

Conduction remains approximately 19.2 seconds because:

* the Thomas solve is intentionally serial;
* each Newton pass requires ordered coefficient, residual, solve, and update phases;
* barriers are required for correctness;
* the 2638-cell workload is too small for large conduction teams.

The final speedup is therefore Amdahl-limited. Increasing threads beyond 16 was not attempted because the machine has 16 physical cores and the task excludes redesigning or parallelizing the Thomas solver.

---

## Remaining serial fraction

The remaining material serial work is:

1. Thomas tridiagonal solve inside every conduction Newton pass;
2. limiter and small Armadillo reconstruction expressions;
3. Rusanov combination and source addition around the parallel face batches;
4. boundary updates and four ghost decodes;
5. CFL calculation;
6. output and diagnostic serialization;
7. deterministic scalar convergence decisions and counter merges.

These regions were deliberately retained because they are sequential by method, too small to amortize another region, or explicitly out of scope.

No MPI, multiprocessing, `std::thread`, pthread, GPU, SIMD intrinsic, fast-math, local time stepping, super-time-stepping, cyclic reduction, new EOS solver, new Newton formulation, or new tridiagonal solver was introduced.

---

## Production recommendation

### Recommended runtime

For maximum throughput on this audited 12-performance-core + 4-efficiency-core machine:

```bash
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_NUM_THREADS=16
```

No affinity variable was required for correctness or for the reported result. The code explicitly uses `schedule(static)` for cell loops.

`OMP_NUM_THREADS=12` is a reasonable lower-power alternative and is only about 0.9% slower in the measured median, but 16 threads was the best stable measured configuration and is therefore the primary production recommendation.

The conduction region automatically caps itself at four threads and reports the selected team in the runtime banner.

### Expected runtime

```text
20 physical seconds:   approximately 36.35 s wall time
500 physical seconds:  approximately 908.65 s
                       approximately 15.14 minutes
```

The 500-second estimate is a linear projection from the measured final 20-second benchmark and assumes the same mesh, timestep behavior, output settings, and machine load.

### AppleClang setup

```bash
brew install libomp

cmake -S . -B build_omp \
  -DCMAKE_BUILD_TYPE=Release \
  -DCHROMO_ENABLE_OPENMP=ON \
  -DLIBOMP_ROOT="$(brew --prefix libomp)"
cmake --build build_omp -j4
```

Linux GCC and Clang builds continue to use standard CMake OpenMP discovery.

### Serial fallback

```bash
cmake -S . -B build_serial \
  -DCMAKE_BUILD_TYPE=Release \
  -DCHROMO_ENABLE_OPENMP=OFF
cmake --build build_serial -j4
```

This build has no `libomp` dependency and preserves the original serial numerical path.

---

```text
OPENMP_PARALLELIZATION: GO
```

Acceptance basis:

* OpenMP runtime definitely active;
* explicit fail-closed ON/OFF build behavior;
* serial fallback clean and independent of libomp;
* serial and OpenMP tests all pass;
* OpenMP one-thread output matches serial byte-for-byte;
* OpenMP multi-thread output matches serial byte-for-byte;
* three final 12-thread runs are hash-identical;
* profiling and EOS counters are thread-safe and deterministic;
* no unresolved ownership, resize, exception, or synchronization issue remains;
* measured 2.170x same-session speedup is material and stable;
* no physics or numerical method changed.

Stop here for independent audit. Do not begin MPI, GPU, local-time-stepping, or further numerical optimization.
