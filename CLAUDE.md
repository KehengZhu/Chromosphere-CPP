# Chromosphere2026 — Agent Instructions

## Writeup standards

The project writeup lives in `docs/writeup-overleaf/`. The primary sources are `main.tex` and `paper.tex`, with bibliography in `reference.bib`. They have deliberately different roles. Do **not** treat `paper.tex` as a shortened copy of `main.tex`, and do not automatically edit both files unless the task explicitly asks for both.

### `main.tex` — engineering and development record

`main.tex` is the authoritative technical record of Chromosphere2026. It should explain both the current implementation and how the project reached it. It may contain substantially more detail than a paper, including:

- the current release architecture and numerical workflow;
- governing models and important derivations;
- implementation details that materially affect correctness, stability, convergence, performance, or reproducibility;
- validation and convergence evidence;
- important design decisions and alternatives;
- major failed approaches and what they taught us;
- historical, experimental, inactive, and retired paths;
- remaining limitations and future work.

Preserve useful technical history, but do not turn `main.tex` into a raw change log, command transcript, code dump, or concatenation of recap files. Organize development history by the scientific or numerical question being answered rather than merely by date. When useful, label material clearly as **Active release path**, **Implemented but inactive**, **Experimental**, **Historical**, **Retired**, or **Planned** so a reader cannot confuse old capability with the current solver.

### `paper.tex` — release scientific paper

`paper.tex` is a concise academic paper describing only the latest implemented release capability, currently centered on `model_column`. It should read like a scientific article, not an engineering report or development diary.

The paper should include only material needed to understand and assess the released calculation:

- physical motivation and scientific context;
- the active physical model and assumptions;
- the important thermodynamic closure;
- the numerical method at paper-appropriate depth;
- initial and boundary conditions;
- computational domain, grid, and material runtime settings;
- focused verification needed to establish credibility;
- the latest defensible preliminary results;
- physical interpretation, limitations, and future coupling goals.

Do not include obsolete implementation history, failed boundary experiments, full test inventories, low-level performance work, inactive model extensions, or detailed derivations that are not needed to understand the release result. In particular, do not present older two-fluid, kinetic-ionization, flare, beam-heating, three-temperature, loop, or multidimensional capabilities as active unless the current code path being reported actually uses them.

### Code is the source of truth

Before writing or revising **each substantive section**, inspect the current implementation that section describes. Do not write a section from memory and do not infer current behavior solely from an existing `.tex` file, an old recap, a plan, or source-code comments.

For `model_column` work, inspect the relevant current code as needed, including:

- `scenarios/model_column.cpp` / `.hpp`;
- `scenarios/model_c7.cpp` / `.hpp`;
- `scenarios/mesh.cpp` / `.hpp`;
- `src/eos.cpp` and `eos.hpp`;
- `src/rhs.cpp`, `src/flux.cpp`, `src/state.cpp`, `src/grid.cpp`, and `src/integrators.cpp`;
- `chromosphere.hpp` and `chromo_main.cpp`;
- the production run scripts and relevant tests.

Use documents under `docs/` to recover motivation, validation evidence, measured results, and development history, but cross-check any statement about **current behavior, defaults, equations actually advanced, boundary semantics, grid settings, or runtime configuration** against the code. If code and documentation disagree, describe the code and fix the stale documentation rather than propagating the discrepancy.

For quantitative results, identify the exact run/configuration that produced them. Never combine a numerical value from one configuration with prose describing another. Distinguish clearly between production defaults, controlled comparison settings, and historical settings.

### Writing quality

Write concise, precise academic English. Technical depth should be sufficient to make the model understandable and reproducible, but not exhaustive for its own sake.

- Prefer direct sentences and compact paragraphs.
- Give each paragraph one clear purpose.
- Do not repeat the same explanation in several sections.
- Avoid promotional or vague language such as "robust", "realistic", or "accurate" unless the surrounding text states what was actually validated.
- Use equations when they define the physical model or clarify an important numerical method; do not reproduce algebra that adds little explanatory value.
- Spend prose on choices that materially affect physics, correctness, stability, convergence, interpretation, performance, or reproducibility. Compress minor implementation details.
- Use terminology consistently and distinguish concepts that are easy to conflate, such as physical versus numerical conduction, `gamma_E` versus `Gamma1`, hydrodynamic versus conductive boundary data, and coarse-equivalent versus actual refined cell count.
- Distinguish a stable result from a converged result, and an apparently approaching-convergence trend from formal convergence.
- Do not overclaim. If a result remains provisional because an input, boundary discretization, or observable is not converged, say so explicitly.

The target is not maximum detail. The target is the **smallest amount of technical detail that accurately explains the important physics, numerics, evidence, and limitations**.

### Citations and figures

Physical claims should be supported by appropriate literature. Reference papers already stored in `docs/supporting-papers/` should be checked first. Do not invent citations.

Before retaining or adding a numerical figure, verify its source output, run configuration, plotted quantity, and relevance to the current narrative. A figure should support a material conclusion rather than merely document that a diagnostic exists. New visualizations must follow the repository visualization rules below and be entered in `visualize_commands.md`.

### Compiling the writeup

When compiling LaTeX, always build inside `docs/writeup-overleaf/latex-build/` — never let auxiliary files (`.aux`, `.log`, `.fls`, `.pdf`, etc.) land next to the `.tex` sources.

Examples:
```bash
cd docs/writeup-overleaf
latexmk -xelatex -output-directory=latex-build main.tex
latexmk -xelatex -output-directory=latex-build paper.tex
```

Resolve compilation errors, undefined citations, broken references, duplicate labels, and obsolete figure references. Do not spend substantial effort eliminating harmless warnings that do not affect the document.

## Output directories

- **Code outputs** (simulation results, data dumps, logs from `chromo_main` / scenarios): save under `outputs/`.
- **Visualizations** (plots, animations, figures generated from outputs): save under `visualization/` (project root). This directory is gitignored. Plotting scripts live in `util/` and already default their output here — keep new/edited scripts defaulting to `visualization/` too. Never write figures or movies to the repo root, into `util/`, or alongside source code.
- **Both `outputs/` and `visualization/` are organized into per-scenario subdirs** mirroring `scenarios/`: `model_c7/`, `model_column/` (the unified column, which absorbed the old `iso_*`/`gentle_*` runs), `model_flare/`, `analytic_canopy/`, `loops/`, `pfss/`, `events/` (real-event campaigns — 2024-08-01 AR13768, SOL2014, ensembles; keeps `event_ensemble/`, `event_te/`, `sol2014/` nested inside), and `two_fluid/`. `visualization/` additionally has `reference/` (physics-reference plots + their scripts). `_archive/` holds dumps/figures from now-deleted diagnostic code and retired analytic-isentrope stages. **Place new runs/figures in the matching scenario subdir** — never loose in `outputs/` or `visualization/` root.

Do not scatter output files in the repo root or alongside source code.

## visualize_commands.md

`visualize_commands.md` (project root) is the catalog of the exact command used to generate every file in `visualization/`. **Whenever you generate a visualization, add or update its entry** in `visualize_commands.md`: the output filename, the exact command run (assume the project root as the working directory), and a one-line description. Keep it in sync — if you change a script's invocation or add a new figure/movie, update the catalog in the same change.

## Visualization
When asked to generate a movie, generate frames in parallel if possible.

When plotting field-aligned profiles that span the chromosphere through the corona (e.g. loop flare runs), use a **split x-axis**: zoom in on the thin chromosphere+TR (sub-Mm) and compress the extended corona (tens of Mm) using different x-scales (a continuous piecewise transform via matplotlib's `set_xscale("function", ...)` works well — see `util/animate_loop_flare.py`). Mark the chromosphere/TR/corona region boundaries with translucent vertical dashed lines.

For a **full closed loop** (footpoint A → apex → footpoint B, `topology=full`), plot vs **arc length** `s`, not height — on a closed loop height folds the two legs onto each other. Apply the split symmetrically: zoom **both** footpoints and compress the coronal middle, laying the loop out flat `[footpoint A | apex | footpoint B]`, with the apex marked by a dotted line and region boundaries on both legs. `util/animate_loop_flare.py` auto-detects topology from the `.dat` and switches between height (half loop) and arc-length (full loop) accordingly.

## Default OpenMP runtime

For production-shaped simulations on the current Apple Silicon workstation, use the OpenMP build and default to 12 threads:

```bash
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_NUM_THREADS=12
```

Prefer `scripts/run_chromo_omp.sh` for agent-run simulations. It applies these defaults while allowing explicit environment overrides. The 12-thread default uses the performance cores and was measured within about 1% of the 16-thread maximum-throughput result while leaving the efficiency cores available for the operating system, compilation, Python analysis, and other agent work.

Exceptions:

- Use `OMP_NUM_THREADS=1` for serial/OpenMP numerical-equivalence checks.
- Use `OMP_NUM_THREADS=16` only when measuring maximum throughput on this workstation.
- Use a build configured with `CHROMO_ENABLE_OPENMP=OFF` for a true serial baseline.
- Do not assume 12 threads is optimal on different hardware; rerun thread scaling after changing machines.
- Do not hard-code the thread count in C++, CMake, or the numerical algorithm.

## `model_column` release numerical method

A normal `model_column` run uses 1D field-aligned hydrodynamics with Gamma/Saha equilibrium thermodynamics (ionization energy included, production `Gamma1` table), MUSCL reconstruction of **`(ln ρ, V, ln p)`** with the MC3/Koren limiter (β = 2), the **mixture Roe characteristic numerical flux**, equilibrium-reference well balancing, physical conduction only, and a 22,000 K external conductive reservoir at the physical outer face, on the coarse-equivalent N = 500 / R4 mesh (661 actual cells) at CFL 0.50.

**No environment variable is required to select the release solver or reconstruction.** `ISO_RIEMANN=rusanov` and `ISO_RECONSTRUCTION=lnrho-v-lnt` are reference-only overrides kept for regression and controlled numerical comparison; never use them in production. Full description and evidence: `docs/model_column_release_numerics_recap.md`. Release checks: `scripts/release_validation.sh`.

Known limitation to preserve in any writeup: **N500 long-duration evaporation mass flux is not established as < 10 % grid-converged.** Distinguish the release numerical method from the resolution required for a particular quantitative scientific claim.

## Validated production CFL and long-run controls

`CHROMO_CFL=0.50` is the **validated production runtime default for the reduced `model_column` release configuration: coarse-equivalent N=500, R4 outer refinement, 661 actual cells, physical-face 22,000 K conductive boundary, and physical conduction only**. Evidence: `docs/coarse_model_column_physical_conduction_recap.md`.

The conservative hard-coded solver default remains `CHROMO_CFL=0.25`, and 0.25 stays the comparison reference. Do not change the C++ default. Rerun the sweep before trusting 0.50 on different hardware or a materially different model shape.

The `model_column` scenario itself supplies the override-preserving defaults for the validated Gamma/Saha N=500/R4 physical-conduction-only release model, including the production Gamma table and model/grid `ISO_*` settings. Use `scripts/run_chromo_realtime.sh` for long production runs; it adds CFL 0.50, the 12-thread OpenMP runtime, and low-I/O runtime defaults but does not independently define the physical model. Explicit environment assignments still win. `scripts/run_chromo_omp.sh` remains the generic OpenMP launcher.

Two run-control variables matter for long runs:

- `CHROMO_T_END` — absolute stop time in physical seconds. When it is set, the legacy `time_mult`-derived step cap no longer applies, so a long run cannot be silently truncated.
- `CHROMO_FRAME_DT` — snapshot cadence in **physical seconds** instead of the legacy `time_mult`-derived step stride. Required whenever output is enabled on a long run; without it a 1000 s run emits tens of thousands of ASCII snapshots and the wall time is dominated by formatting.

`CHROMO_STEP_CAP` sets an explicit step cap and overrides both rules; invalid or non-positive values are rejected rather than silently replaced. Every run prints `termination=end_time|step_cap|other` together with the requested end time, final time, final step, and effective cap — always check that a long run ended with `termination=end_time`.

Low-I/O production run:

```bash
CHROMO_T_END=1000 scripts/run_chromo_realtime.sh out.txt \
    full no-ionization model_column - 20.0 no-cooling
```

Same run with snapshots at one frame per physical second:

```bash
CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=1 \
scripts/run_chromo_realtime.sh out.txt \
    full no-ionization model_column - 20.0 no-cooling
```

## Repository-wide engineering discipline

Treat every task as a production-quality engineering task. The default priority order is:

1. Correct and reliable behavior.
2. Scientifically and numerically accurate results.
3. Reproducible evidence that supports the decision being made.
4. Simple, maintainable implementation.
5. Performance and convenience improvements.

Do not trade correctness for speed, but also do not confuse more code, more abstraction, or more validation with higher quality.

### Avoid overengineering

- Prefer the smallest change that clearly solves the actual problem.
- Use focused helper functions instead of introducing a framework, class hierarchy, policy system, or generalized subsystem without a demonstrated need.
- Do not refactor unrelated code while implementing a narrow task.
- Do not design for hypothetical future requirements unless the current task explicitly depends on them.
- Do not add layers of wrappers, configuration objects, templates, or indirection when a direct implementation is easier to review and equally safe.
- Preserve existing interfaces and behavior unless changing them is necessary for correctness or explicitly requested.

### Validate according to risk

Verification should answer material questions about correctness, reliability, numerical behavior, or performance. Do enough to establish trustworthy evidence, then stop.

- Prioritize tests that exercise the changed behavior and realistic production paths.
- Re-run broad test suites after changes that could affect them, but do not repeatedly rerun unchanged checks without a concrete reason.
- Use hashes and byte comparisons only when exact identity is itself important, such as regression preservation, deterministic output, or legacy compatibility.
- Do not generate or compare hashes for every intermediate artifact when the result cannot affect the decision.
- Do not require byte identity where numerical or physical equivalence is the correct standard.
- Use absolute errors, relative errors with defensible denominator floors, and physically meaningful observables where appropriate.
- Treat small timing variation, harmless formatting changes, and differences far below discretization or physical uncertainty as non-blocking unless they indicate a real bug.

### Focus on material issues

Investigate a discrepancy when it is scientifically meaningful, numerically suspicious, persistent, or large enough to change the production decision. Do not spend substantial effort on cosmetic details or artificial edge cases that cannot affect real use.

A newly discovered issue is a blocker only if it can materially affect one or more of:

- correctness;
- scientific interpretation;
- numerical stability or convergence;
- data integrity;
- requested output or termination behavior;
- production performance;
- reproducibility;
- the final GO/NO_GO decision.

Record unrelated minor issues briefly and leave them for a separate task.

### Use staged validation

For expensive experiments, simulations, or parameter sweeps, use a funnel:

1. Run inexpensive smoke and stability checks.
2. Reject clearly invalid candidates early.
3. Perform detailed comparisons only for plausible candidates.
4. Run long or expensive acceptance cases only for the selected candidate or the smallest set needed to resolve uncertainty.

Do not run several costly production-scale cases when one controlled run is sufficient.

### Protect existing work

- Inspect repository status and existing diffs before editing.
- Preserve unrelated uncommitted changes.
- Do not reset, revert, overwrite, or broadly reformat work outside the task scope.
- Keep patches narrow, reviewable, and attributable to the requested task.
- Keep generated outputs and build products out of source patches unless they are explicitly required deliverables.

### Stop when the decision is supported

Once the requested behavior is correct, the important tests pass, the evidence supports the conclusion, and material risks are documented, stop. Do not continue polishing, benchmarking, or redesigning merely because additional work is possible.

If a target is not achieved, report the measured gap and the most important remaining blocker. Do not silently expand into a new optimization or redesign phase.

### Reporting standard

Final recaps should be compact and decision-oriented. Emphasize:

- what changed;
- why it was necessary;
- what was validated;
- the material results;
- remaining production risks;
- the final status.

Avoid exhaustive command transcripts, repeated hashes, low-value implementation trivia, and long lists of insignificant differences.

The governing principle is:

```text
Reliable and accurate production behavior is mandatory.
Simple and sufficient evidence is preferred.
More validation is not automatically better validation.
```

## Reference Papers
Reference Papers that you may need to understand a physical process or code implementation are in docs/supporting-papers. Check them first before you prompt the user to download the papers.
