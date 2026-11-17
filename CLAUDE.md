# Chromosphere2026 — Agent Instructions

This file holds rules and pointers only. Measurements, arguments and study conclusions live in `docs/`; link them from one line here. Keep this file under ~120 lines: when adding a rule, merge or drop one.

## Where things go (binding)

| What | Where |
|---|---|
| Core shared code (Grid, EOS, run control, OpenMP, profiling) | `src/core/` |
| Release solver | `src/single_fluid/` |
| Release scenario (`model_column`); scenario input data | `scenarios/`; `scenarios/data/` |
| Frozen two-fluid solver and its scenarios (`model_c7`, `model_flare`, `analytic_canopy`, `pfss_field_line`, `model_gentle`); frozen two-temperature solver (`legacy/two_temp/`) and `model_column_2t` | `legacy/` — built only with `-DCHROMO_ENABLE_LEGACY=ON` |
| Python/shell tools | `util/release/`, `util/legacy/`, `util/run/`, shared helpers `util/lib/`, EOS validators `util/eos/` |
| Project scripts (build docs, run, validate, publish) | `scripts/` |
| Measured evidence | `docs/studies/<subject>/` |
| Forward-looking / in-progress plans | `docs/plans/` |
| Authoritative standards | `docs/standards/` |
| Command for every file in `visualization/` | `docs/viz_catalog/` |
| Run outputs / figures (both gitignored) | `outputs/<scenario>/` / `visualization/<scenario>/`; retired material in `_archive/` |

- **Nothing new in the repository root.** Throwaway work goes in the session scratchpad, never in the tree.
- **No orphan scripts.** A diagnostic script is either kept (right `util/` subdir, its figures catalogued) or deleted before commit.
- **Legacy is frozen.** No features in `legacy/` or `util/legacy/`; fix only what breaks its build. Retiring something moves its code, scripts and catalog entries to the legacy locations and marks its studies in the same change.
- **Every study and plan starts with a Status line**: `CURRENT`, `SUPERSEDED by <file>`, `HISTORICAL`, `PLAN`, or `DONE → <recap>`. Whoever supersedes a document updates the old header and the Evidence table in `docs/doxygen/pages/validation.md`.
- **Every generated figure or movie** gets an entry (filename, exact command from the project root, one-line description) in the matching `docs/viz_catalog/` file, in the same change.

## Release invariants

Full statements, numbers and evidence: **`docs/standards/release_facts.md` — read it before writing about the release anywhere.**

- `model_column` means exactly one thing: the release single-fluid equilibrium-mixture solver, state `U = (rho, rho u, E)`; ionization and carrier densities are derived, never stored. One step = explicit MUSCL-Hancock/Godunov hydro, then implicit physical conduction (Lie splitting). No TRAC, beam heating, coronal heating, radiative cooling or artificial conduction in the release path, not even behind flags.
- Face flux: SWMF exact-Riemann Godunov at **frozen composition** (`gamma = 5/3`, ionization energy in `E0`), face-local Rusanov fallback counted in `godunov.*`. `Gamma1` is the *other* (equilibrium) limit, not "wrong". Saha re-equilibration after each step is an open modelling question — do not redesign ionization without asking.
- Release configuration: `(ln rho, u, ln p)` MC3/Koren reconstruction with the `ln p` slot hydrostatic (`Grid::hydrostatic_reconstruct`, no switch), N=500/R4 mesh (661 cells), CFL 0.50, fixed outer reservoir with one 22 kK ghost temperature. The C++ default CFL stays 0.25.
- The 1600 km lower boundary is a **mid-chromospheric truncation, never "photospheric"**.
- The release is reference-free **and** exactly balanced at rest for `B = 1` (hydrostatic reconstruction, since 2026-10-12); "not well-balanced" describes the earlier plain `ln p` scheme only. `eq_wb` is legacy-only and never set on a release Grid.
- Storage is **double**; nothing in the release path names `float`/`double` (use `Real`/`Vec`). `CHROMO_STATE_FLOAT32=ON` is diagnostic only. Any float32-era evaporation number is ~30 % high.
- Do not overclaim double: it did not remove the first-interior-face artifact (the hydrostatic reconstruction did) or the ~13 % lower-chromosphere flux spread (physical: adiabatic contraction). N500 long-run evaporation flux is not established < 10 % grid-converged — a resolution limit, separate from precision.
- `ISO_RIEMANN`, `ISO_RECONSTRUCTION` are reference-only overrides; never production. Rejected (need new evidence to revisit): the `ModTransitionRegion` momentum split, a characteristic lower boundary, selective mixed precision.

## Running

- OpenMP build, 12 threads (`OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 OMP_NUM_THREADS=12`) via `scripts/run_chromo_omp.sh`; long production runs via `scripts/run_chromo_realtime.sh`. Exceptions: 1 thread for serial-equivalence checks, 16 only for max-throughput measurement, `CHROMO_ENABLE_OPENMP=OFF` for a true serial baseline. Never hard-code a thread count; rerun scaling on new hardware.
- Long runs: set `CHROMO_T_END` (physical s) and, with output on, `CHROMO_FRAME_DT`; confirm the run ended `termination=end_time`. Release checks: `scripts/release_validation.sh`.

```bash
CHROMO_T_END=1000 scripts/run_chromo_realtime.sh out.txt full no-ionization model_column - 20.0 no-cooling
```

## Documentation duties (part of every change, not a follow-up)

- **API reference** (Doxygen; config `docs/doxygen/Doxyfile`, prose in `docs/doxygen/pages/`). Any change to code, behaviour, configuration, numerics, physics, scenarios, I/O, validation or the writeup requires: (1) update the affected narrative pages — rebuilding does not update prose; (2) `scripts/build_docs.sh`, warning-free; (3) `scripts/publish_docs.sh`. The pre-commit hook checks (2) only, and only where `scripts/install_hooks.sh` was run.
- The site `https://kehengphysics.site/docs/chromosphere/` is **public** and this repository is private: never enable `SOURCE_BROWSER`, `VERBATIM_HEADERS`, `INLINE_SOURCES` or `REFERENCES_LINK_SOURCE`. Every public entity gets a real description with units in `[...]`.
- `README.md` is a quick start only; technical detail belongs in the Doxygen pages.
- **Writeup** (`docs/writeup-overleaf/`, a separate Overleaf-linked repo): physics, numerics, BC, defaults, validation or interpretation changes trigger a catch-up per `docs/standards/writeup_standards.md`. Code is the source of truth. `paper.tex` = current release only, concise; `main.tex` = engineering record with inactive/historical material labelled. Never commit or push either repo unless asked. Successful `latexmk` runs mirror the PDFs to `docs/chromosphere_writeup.pdf` / `docs/chromosphere_paper.pdf` via the `.latexmkrc` hook — keep the hook; run `scripts/sync_writeup_pdf.sh` by hand otherwise.

## Visualization

- Generate movie frames in parallel.
- Chromosphere-to-corona profiles: split x-axis (zoom chromosphere+TR, compress corona; `set_xscale("function", ...)`), region boundaries as translucent dashed lines. Full closed loops: plot vs arc length with both footpoints zoomed and the apex dotted. Reference implementation: `util/legacy/animate_loop_flare.py`.

## Engineering discipline

Full text: `docs/standards/engineering_discipline.md`. Priority: correct → scientifically accurate → reproducible evidence → simple → fast.

- Smallest change that solves the problem; no frameworks or speculative generality; no unrelated refactors.
- Validate according to risk; byte identity only where identity itself matters; staged funnel for expensive runs.
- Protect existing work: inspect status/diffs first, preserve unrelated changes, keep patches narrow, keep build products out.
- Stop when the decision is supported; recaps are compact and decision-oriented (what changed, why, validation, results, risks, status).

## Reference papers

Check `docs/supporting-papers/` before asking the user to download a paper.
