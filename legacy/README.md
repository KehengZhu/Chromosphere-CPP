# legacy/ — the frozen two-fluid solver and its scenarios

**Status: FROZEN.** Nothing here is part of the release. No new features go in
this directory; fix only what breaks the legacy build.

## What is here

| Path | What it is |
| --- | --- |
| `two_fluid/` | The historical seven-row two-fluid / three-temperature / finite-rate-ionization research solver (`two_fluid.hpp`, `state.cpp`, `flux.cpp`, `rhs.cpp`, `integrators.cpp`). Every entry point rejects a Grid that carries a `Gamma1` table. |
| `scenarios/legacy_scenarios.hpp` / `.cpp` | Name dispatch for the scenarios below (`make_legacy_scenario`, with each scenario's environment preset) and the shared open boundary update `apply_open_bcs`. |
| `scenarios/model_gentle.hpp` / `.cpp` | `model_gentle`, the historical two-fluid C7 column. It used to share `scenarios/model_column.cpp` with the release; it is now its own frozen copy. |
| `scenarios/model_c7.hpp` / `.cpp` | `model_c7`, the quiet-Sun C7 baseline, and the Route-B photoionization closure `c7_route_b_photoionization`. |
| `scenarios/model_flare.hpp` / `.cpp` | `model_flare`, explosive evaporation driven by a nonthermal electron beam. |
| `scenarios/analytic_canopy.hpp` / `.cpp` | `analytic_canopy`, the exponential magnetic-canopy `B(z)` overlay on the C7 column. |
| `scenarios/pfss_field_line.hpp` / `.cpp` | `pfss_field_line`, a tabulated PFSS field line. |

They still use the shared layer: `src/core/` (Grid, EOS, physics formulas),
`scenarios/c7_profile.*` (the C7 table), `scenarios/mesh.*` (static refinement)
and `scenarios/data_file_parser.*` (the field-line table reader). Nothing in
`src/` or `scenarios/` includes anything from here except behind
`#ifdef CHROMO_ENABLE_LEGACY` (`src/main.cpp`, `scenarios/scenario.cpp` and the
legacy tests in `tests/chromo_tests.cpp`). `src/core/chromosphere.hpp` still
carries the two-fluid `Grid` fields (`eq_wb`, `cons::`, the two-fluid ghost
buffers); splitting `Grid` per solver is a separate task.

## How to build it

Not built by default. A default build neither compiles nor links any file
here, and asking it for a legacy scenario fails at once with a message naming
the option. Configure a separate tree:

```sh
cmake -S . -B build_legacy -DCMAKE_BUILD_TYPE=Release \
      -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_ENABLE_LEGACY=ON
cmake --build build_legacy -j 8
(cd build_legacy && ctest --output-on-failure)   # includes the legacy tests

CHROMO_BINARY=build_legacy/chromo_main scripts/run_chromo_omp.sh \
    outputs/model_c7/out.txt full ionization model_c7
```

The option adds the files above to the `chromosphere` library and defines
`CHROMO_ENABLE_LEGACY=1` for every target that links it. Combining it with
`-DCHROMO_STATE_FLOAT32=ON` is allowed.
