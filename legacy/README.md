# legacy/ — the frozen two-fluid and two-temperature solvers and their scenarios

**Status: FROZEN.** Nothing here is part of the release. No new features go in
this directory; fix only what breaks the legacy build.

## What is here

| Path | What it is |
| --- | --- |
| `two_fluid/` | The historical seven-row two-fluid / three-temperature / finite-rate-ionization research solver (`two_fluid.hpp`, `state.cpp`, `flux.cpp`, `rhs.cpp`, `integrators.cpp`). Every entry point rejects a Grid that carries a `Gamma1` table. |
| `two_temp/` | The two-temperature (`T_e != T_i`) test-study solver, state `U = (rho, rho u, E, E_e)` (`two_temp.hpp`, `closure.cpp`, `hydro.cpp`, `conduction.cpp`). Moved from `src/two_temp/` on 2026-11-16 when its study closed (`docs/studies/conduction/two_temperature_test_study.md`). Its one edge outside `legacy/` and `src/core/` is `src/single_fluid/exact_rs.hpp`, the standalone SWMF ideal-gas Riemann solver port. |
| `scenarios/legacy_scenarios.hpp` / `.cpp` | Name dispatch for the scenarios below (`make_legacy_scenario`, with each scenario's environment preset) and the shared open boundary update `apply_open_bcs`. |
| `scenarios/model_gentle.hpp` / `.cpp` | `model_gentle`, the historical two-fluid C7 column. It used to share `scenarios/model_column.cpp` with the release; it is now its own frozen copy. |
| `scenarios/model_c7.hpp` / `.cpp` | `model_c7`, the quiet-Sun C7 baseline, and the Route-B photoionization closure `c7_route_b_photoionization`. |
| `scenarios/model_flare.hpp` / `.cpp` | `model_flare`, explosive evaporation driven by a nonthermal electron beam. |
| `scenarios/analytic_canopy.hpp` / `.cpp` | `analytic_canopy`, the exponential magnetic-canopy `B(z)` overlay on the C7 column. |
| `scenarios/pfss_field_line.hpp` / `.cpp` | `pfss_field_line`, a tabulated PFSS field line. |
| `scenarios/model_column_2t.hpp` / `.cpp` | `model_column_2t`, the two-temperature twin of the release column on `two_temp/`. A frozen copy of the column IC/BC at commit 7e06466: it still reads the historical `ISO_*` names and `TT_OUTER_TI`, so the study's runs reproduce bit for bit. |

They still use the shared layer: `src/core/` (Grid, EOS, physics formulas),
`scenarios/c7_profile.*` (the C7 table), `scenarios/mesh.*` (static refinement)
and `scenarios/data_file_parser.*` (the field-line table reader). Nothing in
`src/` or `scenarios/` includes anything from here except behind
`#ifdef CHROMO_ENABLE_LEGACY` (`src/main.cpp`, `scenarios/scenario.cpp` and the
legacy tests in `tests/chromo_tests.cpp`). `src/core/chromosphere.hpp` still
carries the two-fluid `Grid` fields (`eq_wb`, `cons::`, the two-fluid ghost
buffers) and the two-temperature ones (`two_temperature`, `tt_outer_Ti_neumann`,
the `tt_*` ghost buffers), none read by the release solver; splitting `Grid` per
solver is a separate task.

## Environment names

The legacy scenarios keep the historical `ISO_*` environment names (and
`ISO_REFINE_*` for the mesh, `TT_OUTER_TI` for `model_column_2t`) exactly as
they stood at commit 7e06466. The release `model_column` was renamed to
`CHROMO_*` names on 2026-11-16 and stops on any `ISO_*` variable, so the same
shell cannot silently drive both. The one shared exception is the table path:
the driver reads `CHROMO_EOS_TABLE` for every scenario, and the
`model_column_2t` preset sets that name. The old-to-new table is in
`docs/doxygen/pages/configuration.md`.

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
CHROMO_BINARY=build_legacy/chromo_main CHROMO_T_END=1000 scripts/run_chromo_realtime.sh \
    outputs/model_column_2t/out.txt full no-ionization model_column_2t - 20.0 no-cooling
```

The option adds the files above to the `chromosphere` library and defines
`CHROMO_ENABLE_LEGACY=1` for every target that links it.
