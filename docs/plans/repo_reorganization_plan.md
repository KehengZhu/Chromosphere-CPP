# Repository reorganization plan

**Status: PLAN** — in progress, started 2026-10-07. Mark DONE, with the merge commit, when phase 5 is closed.

## Why

The skeleton was sound (`src/<solver>/`, `docs/studies/<subject>/`, Doxygen pages), but three things made the tree hard to navigate:

- **Fat entry points.** Root `CLAUDE.md` (258 lines) carried study conclusions and their arguments; `visualize_commands.md` (2789 lines) mixed current and retired runs, had duplicate section numbers (§12, §13 twice), still listed retired `model_isentropic` runs, and catalogued the withdrawn `_f64` hybrid-build run as current.
- **History mixed with the release.** The two-fluid solver and its five scenarios compiled into every build; its C7 profile library was entangled with the release IC (`scenarios/model_column.cpp` includes `model_c7.hpp`), and the `model_gentle` two-fluid leg lived inside `model_column.cpp`. `docs/design/` held an authoritative standard, forward plans, background notes and a tool profile side by side. Study documents had no uniform current/superseded marker.
- **Flat directories.** `util/` held 99 scripts with no structure; release diagnostics sat next to scripts for retired scenarios. Six core headers and `chromo_main.cpp` sat in the repository root although `src/` existed.

## Target layout

```
CMakeLists.txt  README.md  CLAUDE.md  .gitignore        <- the only root files
src/
  main.cpp                     (was chromo_main.cpp)
  core/                        chromosphere.hpp grid.cpp eos.hpp eos.cpp parallel.hpp
                               profiling.hpp profiling.cpp run_control.hpp
                               (physics.hpp too, unless only legacy code still needs it)
  single_fluid/                release solver, unchanged
  two_temp/                    experimental two-temperature solver, unchanged
scenarios/                     release + experimental scenarios only
  scenario.* mesh.* data_file_parser.* model_column.* c7_profile.*
  data/                        input data, path unchanged (runtime paths depend on it)
legacy/                        NOT built by default: -DCHROMO_ENABLE_LEGACY=ON
  README.md                    what is here, how to build it, that it is frozen
  two_fluid/                   was src/two_fluid/
  scenarios/                   model_c7, model_flare, analytic_canopy, pfss_field_line,
                               the model_gentle leg
fortran/  tests/  cmake/  scripts/  data/  include/   unchanged roles
util/
  lib/        shared Python helpers imported by other scripts
  release/    model_column diagnostics, sweeps, plots
  two_temp/   model_column_2t scripts
  legacy/     two-fluid scenarios, iso_*, gentle_*, loop/pfss/flare/event scripts
  run/        sweep and ensemble drivers
  eos/        unchanged (CMake builds the validators from here)
docs/
  README.md           map of docs/
  standards/          authoritative rules: writeup_standards.md, release_facts.md
  plans/              forward-looking and in-progress plans (this file)
  studies/<subject>/  measured evidence; every file starts with a Status line
  viz_catalog/        command catalog for every file in visualization/ (was visualize_commands.md)
  doxygen/ presentations/ posters/ supporting-papers/ archive/ writeup-overleaf/
```

## Decisions

- **Two-fluid code: option (b).** Move it to `legacy/`, off by default behind `CHROMO_ENABLE_LEGACY`. With the option off, asking for a legacy scenario fails with a message naming the option. `chromosphere.hpp` keeps its two-fluid `Grid` fields for now; splitting `Grid` per solver is a refactor and a separate task.
- **Code move: yes.** It is a pure relocation. The release output must stay **bitwise identical**, which is the one place byte comparison is the right standard.
- **`model_column_2t` / `src/two_temp/` stay built.** They are experimental, not legacy.
- **`docs/design/codexpro.md`** is a tool workspace profile, not project documentation. Remove it; git history keeps it.

## Phases

| Phase | Content | Owner |
|---|---|---|
| 0 | Commit the pending T,p sweep work so the reorganization starts clean | done, `425e37c` |
| A | Code move + legacy split (`src/`, `scenarios/`, `legacy/`, CMake, tests, `fortran/` include paths, `util/eos` include paths, Doxyfile, `architecture.md`, `configuration.md` build options, `scenario_reference.md`, README build section) | subagent, worktree |
| B | Docs + util: split `docs/design/` into `standards/` and `plans/`, Status lines on every study, `visualize_commands.md` → `docs/viz_catalog/` (current vs `_historical.md`), `util/` reorganization with every path reference updated (docs, catalog, `scripts/release_validation.sh`, presentations' READMEs), `docs/README.md` | subagent, worktree |
| C | Merge A and B; rewrite `CLAUDE.md` to rules + pointers, moving its release facts to `docs/standards/release_facts.md`; Doxygen build warning-free; publish; check the writeup for moved paths (report only, the Overleaf repo is not committed) | main session |
| 5 | Local scratch: list `tmp/` (46 MB) and the largest `outputs/` runs (13 GB) for the user to decide. **Nothing is deleted without confirmation.** | main session |

A and B run concurrently (at most two subagents). Their file sets are disjoint except `docs/doxygen/pages/*.md` and `README.md`, where A edits code paths and B edits util/doc paths; those are merged in phase C.

## Acceptance

**A**
- Default build (legacy OFF, OpenMP) from a fresh build directory: `ctest` passes.
- Legacy ON build: `ctest` passes, including the legacy tests.
- A 20 s `model_column` run and a short `model_column_2t` run give byte-identical snapshots to a baseline binary built from the pre-move commit, at the same thread count. So does a short `model_c7` run under legacy ON.
- `scripts/release_validation.sh` passes.
- The `CHROMO_STATE_FLOAT32=ON` build compiles.
- Doxygen is warning-free.

**B**
- No dangling references. Every `util/...`, `docs/design/...` and `visualize_commands.md` path cited anywhere in the tracked tree resolves.
- Every moved Python script byte-compiles, and its local imports resolve.
- Two or three release scripts actually run, including the one `release_validation.sh` uses.
- Doxygen is warning-free.

**C**
- Merged tree passes A's checks once more.
- `CLAUDE.md` ≤ ~120 lines.
- Docs published.

## Keeping it tidy: rules for `CLAUDE.md`

These go into `CLAUDE.md` as a short "Where things go" section. Its long-form reasoning goes into `docs/standards/`.

1. **Placement table is binding.** New code, scripts, docs and outputs go where the table says. Nothing new goes in the repository root.
2. **Legacy is frozen.** No features in `legacy/` or `util/legacy/`; fix only what breaks the legacy build. Retiring something means moving its code, scripts and catalog entries to the legacy locations and marking its studies `SUPERSEDED`/`HISTORICAL` **in the same change**.
3. **Every study and plan starts with a Status line**, one of `CURRENT`, `SUPERSEDED by <file>`, `HISTORICAL`, `PLAN`, or `DONE → <recap>`. Whoever supersedes a document updates the old header and the Evidence table in `validation.md`.
4. **No orphan scripts.** A diagnostic script is either kept, which means it lives in the right `util/` subdirectory and its figures are in `docs/viz_catalog/`, or deleted before commit. Throwaway work belongs in the session scratchpad, never in the tree.
5. **`CLAUDE.md` holds rules and pointers only.** Measurements, arguments and study conclusions go into `docs/`, linked from one line. Keep it under ~120 lines; when adding a rule, look for one to merge or drop.
