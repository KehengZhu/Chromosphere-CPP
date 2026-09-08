# `docs/` — what lives where

Four kinds of document live here, and they are not interchangeable.

| If you want… | Go to |
| --- | --- |
| **The API reference** — what a class or function is, with units | `doxygen/` → build with `scripts/build_docs.sh`, published to <https://kehengphysics.site/docs/chromosphere/> |
| **The release write-ups** | `chromosphere_paper.pdf` (short) and `chromosphere_writeup.pdf` (full engineering record) |
| **Why a decision was made, with measurements** | `studies/<subject>/` |
| **A talk or poster** | `presentations/` and `posters/`, one dated folder each |

**Current code is always the source of truth for active behaviour.** The documents under `studies/` record how a decision was reached; several describe configurations that have since been retired, and each says so at the top. When a study and the code disagree, the code is right and the study is history.

## `doxygen/` — the API reference

Config in `doxygen/Doxyfile`, hand-written narrative pages in `doxygen/pages/`. The narrative pages are the ones that go stale, because Doxygen regenerates API listings but cannot update prose:

`configuration.md` · `numerics.md` · `physics_model.md` · `scenario_reference.md` · `io_formats.md` · `validation.md` · `architecture.md`

`doxygen/html/` is a gitignored build artifact. Build with `scripts/build_docs.sh` (must be warning-free), publish with `scripts/publish_docs.sh`. See the root `CLAUDE.md` for the rule that publishing is part of the change, not a follow-up.

## The two write-up PDFs stay at this level on purpose

`chromosphere_writeup.pdf` and `chromosphere_paper.pdf` are **mirrored automatically** from `writeup-overleaf/` by a `$success_cmd` hook in `writeup-overleaf/.latexmkrc` that calls `scripts/sync_writeup_pdf.sh`. Those two exact paths are hard-coded in that script and are what `README.md` and the Doxygen pages link. **Do not move or rename them.**

## `writeup-overleaf/` — a separate Git repository

Overleaf-linked, with its own history; it is not part of this repository's commits. Read `design/writeup_standards.md` before editing `paper.tex` or `main.tex`, and inspect the two worktrees separately.

## `presentations/` and `posters/` — one dated folder each

```
presentations/2026-08-17-group-meeting/    Saha closure, the 1600 km truncation, the SWMF exact-Riemann flux
presentations/2026-07-11-group-meeting/    two-fluid ionization & evaporation talk
posters/2026-07-11-shine/                  SHINE 2026 poster (beamerposter, UMich theme) + web export
posters/2026-06-23-agu26/                  AGU 2026 poster draft and abstract
```

Name new ones `YYYY-MM-DD-<event-or-topic>/`. Keep each talk's sources, its own figures and its build directory inside its own folder so an old deck stays buildable after the repository moves on around it.

## `studies/` — measured evidence, by subject

| Directory | Subject | Start here |
| --- | --- | --- |
| `studies/numerics/` | Flux, reconstruction, well-balancing, mesh refinement, storage precision | `swmf_godunov_flux_experiment.md` — the release numerical flux; `reference_free_release_recap.md` — the release as shipped; `float32_precision_control_experiment.md` — why lower-chromosphere mass flux is not a steady-state statement |
| `studies/boundaries/` | Inner truncation reservoir and outer conductive boundary | `boundary_conditions_plan.md`, `q_of_T_boundary_condition.md` |
| `studies/eos-ionization/` | Saha closure, the CRASH `Gamma1` table, ionization/recombination rates | `crash_gamma_eos_plan.md` then the `stage*_recap.md` sequence |
| `studies/conduction/` | Physical and numerical conduction, the release mesh and its CFL | `coarse_model_column_physical_conduction_recap.md` |
| `studies/validation/` | Convergence scans and long-run acceptance | `long_run_output_and_cfl_validation.md` |
| `studies/evaporation/` | Gentle (conduction-driven) and explosive (beam-driven) evaporation | `gentle_evaporation_plan.md`, `explosive_evaporation_plan.md` |
| `studies/performance/` | OpenMP, serial optimization, EOS-inversion speedups | `openmp_parallelization_recap.md` |

The authoritative index of which study is current and which is superseded is the **Evidence documents** table in `doxygen/pages/validation.md`, not this file.

## `design/` — cross-cutting standards and forward-looking plans

| File | What it is |
| --- | --- |
| `writeup_standards.md` | **Authoritative.** How to catch the write-ups up after a code change. Referenced by the root `CLAUDE.md`. |
| `publication_strategy_plan.md` | The AWSoM-R gap this work fills, the publication plan built on it, measured cost against the host's real-time budget, and the instantaneous-Saha justification plan |
| `project_context_zh.md` | **中文.** Plain-language explanation of AWSoM-R, the opening in its chromospheric closure, and what this work delivers — for a reader who has not read Sokolov et al. (2021). Derived from `publication_strategy_plan.md`, which stays authoritative for every planning decision. |
| `scenario_reconciliation_plan.md` | Proposed scenario consolidation; partly implemented |
| `electron_temperature_plan.md` | Three-temperature (`T_e != T_i`) design, two-fluid solver only. **Note:** the *release* `T_e != T_i` question is scoped separately — see `publication_strategy_plan.md` §2.7, which measures τ_eq = 0.09–3.9 ms over the release domain and defers the split to Tier 2. |
| `event-data-to-code.md` | How observations of the 2024-08-01 M8.2 flare become model inputs |
| `codexpro.md` | Saved workspace-profile notes |

## `supporting-papers/` — the literature collection

A stable collection of reference PDFs. Check here before downloading a paper. Left as-is by design: filenames match how they are cited elsewhere in the repository.

## `archive/`

`Keheng_s_chromosphere.pdf` — an early write-up PDF, superseded by `chromosphere_writeup.pdf`. Kept for provenance.

## Build directories

`latex-build/`, `presentation-build/` and `doxygen/html/` are gitignored build output. Nothing source-controlled should ever point into them.
