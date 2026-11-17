# Chromosphere2026

A 1.5D field-aligned solver for the solar chromosphere and the chromosphere-to-corona transition. C++14, with header-only Armadillo vendored in `include/` — no LAPACK, no BLAS, no other dependency.

The production configuration (`model_column`) is a **single-fluid field-aligned equilibrium-mixture** solver: Gamma/Saha equilibrium thermodynamics with the ionization energy included, MUSCL-Hancock hydrodynamics with the SWMF-style exact-Riemann Godunov flux, and implicit physical heat conduction. A separate historical two-fluid research solver is also in the tree and drives the other scenarios.

**This README is a quick start.** Everything technical — equations, numerical method, scenarios, the full command-line and environment reference, output formats, validation — lives in the Doxygen reference and the write-up:

| | |
| --- | --- |
| Code reference | **[kehengphysics.site/docs/chromosphere/](https://kehengphysics.site/docs/chromosphere/)** — physics model, numerics, scenarios, configuration, I/O formats, validation |
| Write-up (current release) | [docs/chromosphere_paper.pdf](docs/chromosphere_paper.pdf) |
| Write-up (engineering record) | [docs/chromosphere_writeup.pdf](docs/chromosphere_writeup.pdf) |
| Reference papers | `docs/supporting-papers/` |

---

## 1. Configure

Nothing to configure for a default build. Three optional pieces:

**OpenMP** (recommended for production runs) needs a separate build directory:

```sh
cmake -S . -B build_omp -DCMAKE_BUILD_TYPE=Release -DCHROMO_ENABLE_OPENMP=ON
cmake --build build_omp -j
```

On macOS with Apple Clang, install `libomp` first (`brew install libomp`); CMake finds it automatically, or point it at another install with `-DLIBOMP_ROOT=<path>`.

**Legacy code** (the two-fluid scenarios `model_gentle`, `model_c7`, `model_flare`, `analytic_canopy`, `pfss_field_line` and the solver they run on, plus the closed two-temperature test study `model_column_2t` and its solver) is frozen under `legacy/` and not built by default; a default build refuses those scenarios with a message naming the option. Build it into its own tree:

```sh
cmake -S . -B build_legacy -DCMAKE_BUILD_TYPE=Release -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_ENABLE_LEGACY=ON
cmake --build build_legacy -j
```

**Python** for the plotting and field-line-extraction scripts:

```sh
bash util/setup_venv.sh          # creates .venv/ and installs the deps
```

## 2. Compile

```sh
cmake -S . -B build              # Release by default
cmake --build build -j
```

This produces `build/chromo_main` (the solver), `build/chromo_tests` (the test suite), and three `eos_validate_*` cross-check tools. Use `-DCMAKE_BUILD_TYPE=Debug` for a debug build; re-run `cmake` after editing `CMakeLists.txt`.

Check the build:

```sh
cd build && ctest --output-on-failure
```

## 3. Choose parameters

A run is `chromo_main [output_path] [mode] [ionization] [scenario] [data_path] [time_mult] [cooling]`, plus environment variables. The three things you normally set:

**Scenario** — decides the physical model, the geometry, and which solver runs.

| Scenario | Solver | What it is |
| --- | --- | --- |
| `model_column` | release single-fluid | **The default and the release configuration.** Field-aligned chromosphere-to-corona column, conduction-driven evaporation. |
| `model_column_2t` | two-temperature (legacy build) | Closed test study: the same column with `T_e != T_i`. Not release. |
| `model_gentle` | two-fluid (legacy build) | Historical research preset on the same column. |
| `model_c7` | two-fluid (legacy build) | Quiet-Sun Model C7 baseline. |
| `model_flare` | two-fluid (legacy build) | Explosive evaporation driven by a nonthermal electron beam. |
| `analytic_canopy` | two-fluid (legacy build) | C7 plus an exponential magnetic-canopy `B(z)`. |
| `pfss_field_line` | two-fluid (legacy build) | A real field line traced from a synoptic magnetogram. |

**Stop time and snapshot cadence** — always set both on a long run:

```sh
CHROMO_T_END=<seconds>     # absolute stop time; without it a run can be truncated by the step cap
CHROMO_FRAME_DT=<seconds>  # snapshot cadence in physical seconds, not steps
```

**CFL** — `CHROMO_CFL=0.50` is the validated production value for `model_column`, and `scripts/run_chromo_realtime.sh` applies it for you. The hard-coded solver default stays at the conservative 0.25.

`model_column` needs **no** environment variable to select its physics or its numerical method — the scenario supplies all of it, and it rejects the historical two-fluid knobs outright. For every other variable, see the *Configuration reference* page in the Doxygen docs.

## 4. Run

The release run, low I/O:

```sh
CHROMO_T_END=1000 scripts/run_chromo_realtime.sh \
    outputs/model_column/run.txt \
    full no-ionization model_column - 20.0 no-cooling
```

The same run with snapshots once per physical second:

```sh
CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=1 \
scripts/run_chromo_realtime.sh outputs/model_column/run.txt \
    full no-ionization model_column - 20.0 no-cooling
```

`run_chromo_realtime.sh` adds the OpenMP runtime defaults (12 threads), CFL 0.50 and low-I/O settings; explicit environment assignments still win. `scripts/run_chromo_omp.sh` is the generic OpenMP launcher, and `build/chromo_main` runs directly if you want no defaults at all.

Always check the run ended with `termination=end_time`, which `chromo_main` prints alongside the requested end time, final time and final step.

The legacy scenarios run the same way from a legacy build:

```sh
build_legacy/chromo_main outputs/model_c7/out.txt      full ionization model_c7
build_legacy/chromo_main outputs/model_flare/out.txt   full ionization model_flare
build_legacy/chromo_main outputs/pfss/out.txt          full ionization pfss_field_line \
    scenarios/data/pfss_qs_20190801.dat
CHROMO_BINARY=build_legacy/chromo_main CHROMO_T_END=1000 scripts/run_chromo_realtime.sh \
    outputs/model_column_2t/out.txt full no-ionization model_column_2t - 20.0 no-cooling
```

Simulation output goes under `outputs/`, organized per scenario.

## 5. Visualize

Plotting scripts live in `util/release/` (retired scenarios: `util/legacy/`) and default their output to `visualization/` (gitignored), organized in the same per-scenario subdirectories as `outputs/`.

```sh
.venv/bin/python util/release/animate_isentropic.py \
    outputs/model_column/run.txt visualization/model_column/run_evolution.mp4 20
```

[`docs/viz_catalog/`](docs/viz_catalog/README.md) is the catalog of the exact command behind every file in `visualization/` — add an entry whenever you generate a new figure or movie.

## 6. Validate

```sh
cd build && ctest --output-on-failure   # full suite, or ./chromo_tests directly
scripts/release_validation.sh           # release-configuration acceptance checks
```

## Repository layout

```
src/core/                                 shared: Grid, EOS closure, physics formulas, OpenMP/profiling
src/single_fluid/                         RELEASE solver, U = (rho, rho u, E)
src/main.cpp                              CLI driver (builds chromo_main)
scenarios/                                model_column IC/BC, static mesh, C7 profile, data tables
legacy/                                   frozen two-fluid and two-temperature solvers + their scenarios (-DCHROMO_ENABLE_LEGACY=ON)
tests/                                    test suite
util/                                     Python diagnostics: release/ legacy/ run/ lib/ eos/
docs/doxygen/                             code reference (config + narrative pages)
docs/writeup-overleaf/                    write-up sources (separate Overleaf repo)
docs/supporting-papers/                   reference PDFs
outputs/  visualization/                  run output and figures, per scenario
```

## Documentation

The code reference is published at **https://kehengphysics.site/docs/chromosphere/** — no sign-in, nothing to install. The Doxyfile sets `SOURCE_BROWSER = NO` and `VERBATIM_HEADERS = NO`, so the site carries the API reference and the narrative pages but **no verbatim source** from this private repository.

`docs/doxygen/html/` is a gitignored build artifact — the site is the only published copy, so it has to be refreshed when the code changes. Build it locally, or build and publish:

```sh
brew install doxygen                 # or: sudo apt-get install doxygen
scripts/build_docs.sh --open         # build and open it in a browser
scripts/publish_docs.sh              # build if needed, then push to the site
scripts/publish_docs.sh --dry-run    # show what would change
```

**When you change documented behavior, update the narrative pages in `docs/doxygen/pages/` too** — Doxygen regenerates the API listings from the sources, but it cannot update prose, so that is the part that goes stale. Then rebuild and publish. `scripts/install_hooks.sh` installs a pre-commit hook that rebuilds and rejects the commit on any Doxygen warning; it deliberately does not publish for you.

The server tree is `/srv/www/kehengphysics.site/docs/`, outside `/var/www/wordpress`, and documents itself in a `README.md` at its root: `docs/index.html` is the hand-edited landing page, `docs/chromosphere/` is the generated reference, and the serving config is the single file `/etc/nginx/snippets/kehengphysics-docs.conf`.

The write-up PDFs in `docs/` are refreshed the same way, by a `latexmk` hook that calls `scripts/sync_writeup_pdf.sh` whenever `main.tex` or `paper.tex` is recompiled.

## References

- Al Shidi, Q., et al. (2019). *Time-Dependent Two-Fluid Magnetohydrodynamic Model of the Solar Chromosphere.* `docs/supporting-papers/`
- Avrett, E. H. & Loeser, R. (2008). *ApJS* **175**, 229 — Model C7.
- Upstream CoMFi (MIT license): https://github.com/qalshidi/phd-dissertation
