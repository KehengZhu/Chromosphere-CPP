# `docs/` — what lives where

Several kinds of document live here, and they are not interchangeable.

| If you want… | Go to |
| --- | --- |
| **The API reference** — what a class or function is, with units | `doxygen/` → build with `scripts/build_docs.sh`, published to <https://kehengphysics.site/docs/chromosphere/> |
| **The release write-ups** | `chromosphere_paper.pdf` (short) and `chromosphere_writeup.pdf` (full engineering record) |
| **Why a decision was made, with measurements** | `studies/<subject>/` |
| **A talk or poster** | `presentations/` and `posters/`, one dated folder each |
| **The rules a change must follow** (writeup catch-up, …) | `standards/` |
| **What is planned or in progress** | `plans/` |
| **The command behind a figure or movie in `visualization/`** | `viz_catalog/` |

**Current code is always the source of truth for active behaviour.** The documents under `studies/` record how a decision was reached; several describe configurations that have since been retired. When a study and the code disagree, the code is right and the study is history.

**Every study and plan opens with a Status line**, directly under its title, in one of five forms:

| Status | Meaning |
| --- | --- |
| `**Status: CURRENT**` | Describes the release as it is, or evidence it still rests on. An optional clause says which part is current. |
| `**Status: SUPERSEDED by <link>**` | Replaced by the linked document; read that one instead. |
| `**Status: HISTORICAL** — <why>` | Records a retired configuration or a non-release solver (two-fluid, `model_gentle`, `model_isentropic`, `iso_*`, float32-era numbers). Kept for provenance. |
| `**Status: PLAN**` | Forward-looking; not yet carried out. |
| `**Status: DONE → <recap link>**` | A plan that was carried out; the recap records the result. |

Whoever supersedes or retires a document updates its Status line and the **Evidence documents** table in `doxygen/pages/validation.md` in the same change.

## `doxygen/` — the API reference

Config in `doxygen/Doxyfile`, hand-written narrative pages in `doxygen/pages/`. The narrative pages are the ones that go stale, because Doxygen regenerates API listings but cannot update prose:

`configuration.md` · `numerics.md` · `physics_model.md` · `scenario_reference.md` · `io_formats.md` · `validation.md` · `architecture.md`

`doxygen/html/` is a gitignored build artifact. Build with `scripts/build_docs.sh` (must be warning-free), publish with `scripts/publish_docs.sh`. See the root `CLAUDE.md` for the rule that publishing is part of the change, not a follow-up.

## The two write-up PDFs stay at this level on purpose

`chromosphere_writeup.pdf` and `chromosphere_paper.pdf` are **mirrored automatically** from `writeup-overleaf/` by a `$success_cmd` hook in `writeup-overleaf/.latexmkrc` that calls `scripts/sync_writeup_pdf.sh`. Those two exact paths are hard-coded in that script and are what `README.md` and the Doxygen pages link. **Do not move or rename them.**

## `writeup-overleaf/` — a separate Git repository

Overleaf-linked, with its own history; it is not part of this repository's commits. Read `standards/writeup_standards.md` before editing `paper.tex` or `main.tex`, and inspect the two worktrees separately.

## `presentations/` and `posters/` — one dated folder each

```
presentations/2026-09-29-tp-sweep/         (T_wall, p_top) evaporation sweep against the analytic closures
presentations/2026-08-17-group-meeting/    Saha closure, the 1600 km truncation, the SWMF exact-Riemann flux
presentations/2026-07-11-group-meeting/    two-fluid ionization & evaporation talk
posters/2026-07-11-shine/                  SHINE 2026 poster (beamerposter, UMich theme) + web export
posters/2026-06-23-agu26/                  AGU 2026 poster draft and abstract
```

Name new ones `YYYY-MM-DD-<event-or-topic>/`. Keep each talk's sources, its own figures and its build directory inside its own folder so an old deck stays buildable after the repository moves on around it.

## `studies/` — measured evidence, by subject

| Directory | Subject | Start here |
| --- | --- | --- |
| `studies/numerics/` | Flux, reconstruction, well-balancing, mesh refinement, storage precision | `swmf_godunov_flux_experiment.md` — the release numerical flux; `reference_free_release_recap.md` — the release as shipped; `state_precision_release_cutover.md` — why the state is `double` |
| `studies/boundaries/` | Inner truncation reservoir, outer conductive boundary, and the evaporation response to the boundary state | `outer_temperature_bc_single_t_recap.md` — the release outer boundary; `twall_mass_flux_relaxation_timescales.md`, `analytic_evaporation_J_table.md`, `tp_sweep_evaporation_closure_validation.md` — the evaporation flux and its analytic closures |
| `studies/eos-ionization/` | Saha closure, the CRASH `Gamma1` table, ionization/recombination rates | `crash_gamma_eos_plan.md` then the `stage*_recap.md` sequence (stages 0–3 current; 4–9 describe the seven-row path the single-fluid solver replaced) |
| `studies/conduction/` | Physical and numerical conduction, the release mesh and its CFL | `coarse_model_column_physical_conduction_recap.md` |
| `studies/validation/` | Convergence scans and long-run acceptance | `long_run_output_and_cfl_validation.md` |
| `studies/evaporation/` | Gentle (conduction-driven) and explosive (beam-driven) evaporation in the two-fluid solver — all HISTORICAL | `gentle_evaporation_plan.md`, `explosive_evaporation_plan.md`. Release evaporation work is under `studies/boundaries/`. |
| `studies/performance/` | OpenMP, serial optimization, EOS-inversion speedups | `openmp_parallelization_recap.md` |

The authoritative index of which study is current and which is superseded is the **Evidence documents** table in `doxygen/pages/validation.md`, not this file.

## `standards/` — authoritative rules

| File | What it is |
| --- | --- |
| `writeup_standards.md` | **Authoritative.** How to catch the write-ups up after a code change. Referenced by the root `CLAUDE.md`. |
| `release_facts.md` | **Authoritative.** What is true about the `model_column` release and how to state it: numerical method, frozen-composition flux, storage precision, validated CFL and run controls, known limitations. Read before writing about the release. |
| `engineering_discipline.md` | Priority order, validation-by-risk, staged validation, protecting existing work, and the reporting standard. |

## `plans/` — forward-looking and in-progress plans

| File | What it is |
| --- | --- |
| `repo_reorganization_plan.md` | The 2026-10 repository reorganization (layout, phases, acceptance) |
| `publication_strategy_plan.md` | The AWSoM-R gap this work fills, the publication plan built on it, measured cost against the host's real-time budget, and the instantaneous-Saha justification plan |
| `project_context_zh.md` | **中文.** Plain-language explanation of AWSoM-R, the opening in its chromospheric closure, and what this work delivers — for a reader who has not read Sokolov et al. (2021). Derived from `publication_strategy_plan.md`, which stays authoritative for every planning decision. |
| `scenario_reconciliation_plan.md` | Proposed scenario consolidation; partly implemented |
| `electron_temperature_plan.md` | Three-temperature (`T_e != T_i`) design, two-fluid solver only. **Note:** the *release* `T_e != T_i` question is scoped separately — see `publication_strategy_plan.md` §2.7, which measures τ_eq = 0.09–3.9 ms over the release domain and defers the split to Tier 2. |
| `event-data-to-code.md` | How observations of the 2024-08-01 M8.2 flare become model inputs (two-fluid event campaigns) |

## `viz_catalog/` — how every figure was made

The exact command behind every file in `visualization/`, split by scenario area: `model_column.md` (release), `model_column_2t.md`, `eos.md`, `reference.md`, `tutorials.md`, `presentations.md`, and `_historical.md` for retired scenarios and configurations. **Every file in `visualization/` has an entry**; add or update it in the same change that generates the figure. Index and conventions: `viz_catalog/README.md`. The scripts themselves live in `util/release/`, `util/two_temp/` and `util/legacy/`.

## `supporting-papers/` — the literature collection

A stable collection of reference PDFs. Check here before downloading a paper. Left as-is by design: filenames match how they are cited elsewhere in the repository.

## `archive/`

`Keheng_s_chromosphere.pdf` — an early write-up PDF, superseded by `chromosphere_writeup.pdf`. Kept for provenance.

## Build directories

`latex-build/`, `presentation-build/` and `doxygen/html/` are gitignored build output. Nothing source-controlled should ever point into them.

## Paths in older documents

Studies and plans are historical records, so their plain-text file citations describe the tree as it was when they were written; only their Markdown links were rewritten after the 2026-10-07 reorganization (`plans/repo_reorganization_plan.md`). To follow an old citation:

| Old path | Now |
| --- | --- |
| `chromo_main.cpp` | `src/main.cpp` |
| `chromosphere.hpp`, `eos.hpp`, `physics.hpp`, `parallel.hpp`, `profiling.hpp`, `run_control.hpp`, `src/{grid,eos,profiling}.cpp` | `src/core/` |
| `src/two_fluid/`, and the older `src/{rhs,flux,state,integrators}.cpp` | `legacy/two_fluid/` |
| `scenarios/{model_c7,model_flare,analytic_canopy,pfss_field_line}.*`, the `model_gentle` leg of `model_column` | `legacy/scenarios/` (C7 table: `scenarios/c7_profile.*`) |
| `util/<script>` | `util/{release,two_temp,legacy,run,lib}/<script>` |
| `visualize_commands.md` | `viz_catalog/` |
| `design/` | `standards/` and `plans/` |
