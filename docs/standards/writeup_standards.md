# Chromosphere2026 Writeup Maintenance Standard

This is the authoritative project-level standard for catching the release writeups up after code changes, auditing them for stale claims, and cleaning obsolete active-release language. The operating principle is:

> **Update forward, then sweep backward.**

Use this standard when asked to perform a writeup catch-up, audit, cleanup, or substantial edit.

## 1. Repository roles

### Main project repository

`/Users/zkeheng/SWMFSoftware/Chromosphere2026` is the source repository for production code, the current numerical implementation, release defaults, tests, validation evidence, release recaps, and engineering history.

### Writeup repository

`docs/writeup-overleaf/` is a separate Overleaf-linked Git working tree. It contains primarily `paper.tex`, `main.tex`, the bibliography, and writeup figures/assets. It is a documentation target, not the source of truth for current solver behavior.

Do not assume that the main repository's status, diff, or commit history includes changes under the writeup repository. For writeup maintenance:

- inspect the two Git roots, statuses, and diffs separately when relevant;
- keep project standards and agent/workflow policy in the main repository, never in `docs/writeup-overleaf/`;
- do not commit or push either repository unless explicitly requested.

## 2. Source-of-truth hierarchy

Resolve current release facts in this order:

1. current production code;
2. active release configuration and defaults;
3. current tests and release-validation evidence;
4. current release recap documents;
5. existing `paper.tex` and `main.tex`;
6. older recap and history documents.

Existing writeup text is never authoritative for current implementation behavior when it conflicts with current code. Inspect code directly for equations actually advanced, active numerical methods, boundary semantics, update ordering, runtime defaults, and the canonical release configuration. Use historical documents for motivation, experimental evidence, design decisions, and development history.

## 3. What triggers a writeup catch-up

Review the writeups when a code or release change affects any of the following:

- governing physical equations;
- EOS or thermodynamic closure;
- conserved, primitive, or reconstructed variables;
- Riemann solver, numerical flux, limiter, or reconstruction;
- well balancing;
- time integration or update ordering;
- physical source terms active in the release;
- physical or numerical conduction;
- boundary conditions, ghost construction, or boundary independence;
- mesh or refinement strategy;
- canonical grid, runtime configuration, or production defaults;
- validation or convergence conclusions;
- known limitations;
- scientifically interpreted results.

Pure internal refactors, diagnostic-only tools, visualization-only changes, and performance optimizations that leave the mathematical behavior unchanged usually do not require a scientific-paper update. They may belong in `main.tex` when they materially affect reproducibility, reliability, or engineering behavior, but they do not automatically belong in `paper.tex`.

## 4. Decide which writeup needs updating

### `paper.tex`: concise scientific release paper

Update `paper.tex` when a change materially affects the current physical model, governing equations, thermodynamic closure, active numerical method, boundary conditions, canonical computational configuration, validation, scientific result, or limitation.

Describe only the current released calculation and the evidence needed to understand and assess it. Obsolete release implementations should normally disappear rather than remain for historical completeness.

### `main.tex`: engineering, numerical, and development record

Update `main.tex` for the current implementation in greater technical depth, important derivations, correctness- or reproducibility-relevant implementation details, validation, significant design decisions, and useful experimental or failed paths.

Historical material may remain when informative, but label it unambiguously as **Active release path**, **Implemented but inactive**, **Experimental**, **Historical**, **Retired**, or **Planned**. Historical material must never read as though it is still the active release.

Do not assume that every change requires both files. Update `paper.tex`, `main.tex`, both, or neither according to these roles.

## 5. Catch-up workflow

### A. Determine what changed

Inspect the relevant code changes, current branch and worktree state, release configuration, tests, release recaps, and validation evidence. Build a small impact map rather than rereading every file blindly.

Examples:

- Riemann-solver change: numerical-method sections, configuration tables, validation, conclusions, and limitations;
- boundary change: boundary sections, ghost semantics, conduction formulation, configuration tables, and boundary verification;
- grid/default change: every canonical grid number, cell count, runtime statement, result, caption, and limitation tied to that configuration.

### B. Update affected current-release sections

Bring the current release description forward to match the implementation. Write the smallest coherent patch that accurately describes the changed physics, numerics, evidence, and limitations.

### C. Sweep backward for stale truth

Search both TeX sources for superseded solver names, reconstruction choices, grid sizes, cell counts, CFLs, boundary distances, source-stage descriptions, retired physical/numerical terms, and obsolete validation conclusions.

Remove stale release material from `paper.tex`. In `main.tex`, preserve useful history but classify it correctly.

### D. Check internal consistency

Equations, prose, tables, captions, configuration summaries, verification, results, conclusions, and limitations must all describe the same release. A new correct paragraph is insufficient if an old table or caption still states the previous configuration.

## 6. Release-model contract

Both writeups must satisfy the following contract. It exists because a release refactor once left the documented model describing a storage layout the solver no longer used.

1. **Self-contained closed statement.** `paper.tex` and `main.tex` must each contain a complete, closed statement of the current release mathematical model — governing equations, the definition of every conserved variable, and every closure needed to make the system mathematically determinate (equation of state, ionization closure, conductive flux, acoustic index). A reader must not have to reconstruct the governing equations from the numerical-method sections.
2. **State variables must match the solver.** The conserved state written in the document must be exactly the state the production solver stores and advances — same variables, same count. Quantities the code derives from that state must be presented as derived, never as independently evolved.
3. **Numerical-method sections are not the source of the model.** Discretization sections may refine or specialize the governing equations, but the equations themselves must be stated first, in a physics section.
4. **Prescribed field-aligned geometry must be explicit.** For a field-aligned model, the general governing equations must show the prescribed `B(s)` (or equivalently `A(s)`, with `A(s)B(s) = const`) and the resulting geometric source terms, *even when the reported run uses constant `B`*. State the constant-`B` specialization explicitly and separately. Do not present prescribed geometry as though the field were evolved, and do not add induction or Lorentz-force terms the release does not solve.
5. **Historical architecture must not read as active.** A superseded state representation, stage, or update path may be documented, but only under an explicit non-active label (`main.tex`) or not at all (`paper.tex`). Removing a solver stage means removing it from the active update sequence, the abstract, the configuration tables, the verification section, and the conclusions — not only from the section that defined it.
6. **Only structurally present physics may be documented as release physics.** Describe a stage as part of the release update path only if it is actually implemented inside the canonical release integrator. A stage that exists elsewhere, or that would run only behind a default-off flag, is not release physics and must carry a non-active label (`main.tex`) or be absent (`paper.tex`). This mirrors the code-side rule that the canonical release integrator should contain only physics the production timestep actually performs; inactive experimental stages belong in the research solver, not embedded in the release path behind flags.
7. **Current code is the source of truth.** Verify every one of the above directly against the production implementation, not against an earlier revision of the writeup or a recap document.

## 7. Numerical and physical accuracy

Verify every substantive current-release claim directly in the implementation before writing. In particular, inspect:

- the actual authoritative variables advanced, and their exact count;
- the exact Roe, Rusanov, reconstruction, and limiter formulation;
- EOS derivatives and acoustic speed;
- active, inactive, and removed source stages;
- which quantities are stored and which are derived from the stored state;
- which boundary data are independent, extrapolated, or EOS-derived;
- ghost-state construction;
- physical versus numerical conduction;
- physical-face versus ghost-center thermal boundary semantics;
- update ordering;
- scenario and launcher production defaults.

Do not infer these solely from comments, old recaps, plans, or existing TeX when current code can answer them. Write equations for the implemented mathematical model, not merely transcriptions of code variable names.

Use concise, precise academic English. Support physical claims with real literature; check `docs/supporting-papers/` first. Before retaining or adding a numerical figure, verify its source run, configuration, plotted quantity, and relevance.

## 8. Avoid configuration mixing

Never combine a numerical value from one configuration with prose describing another. Configuration identity includes at least:

- resolution and refinement;
- solver and reconstruction;
- boundary condition;
- conduction model;
- CFL and end time;
- active source set.

Distinguish canonical production defaults, controlled comparison runs, and historical runs. Identify the exact run/configuration behind quantitative evidence. A stable result is not converged unless convergence evidence supports that claim.

## 9. Cleanup standard

Every catch-up must remove or reclassify superseded active-release statements. Common targets are:

- old canonical resolution or actual cell count;
- old Riemann solver or reconstruction;
- obsolete boundary placement or distance;
- artificial or numerical terms no longer active, or removed from the operator entirely;
- optional source modules no longer in the release flow, or no longer implemented in the release solver at all;
- a solver-source layout or module list that no longer matches the source tree;
- a scenario name documented as selecting a solver it can no longer select;
- limitations that later work resolved;
- convergence claims invalidated by later evidence;
- a superseded conserved-state representation or a removed solver stage still described as active.

Obsolete release details should normally disappear from `paper.tex`. They may remain in `main.tex` only when historically useful and explicitly separated from the active release. Never allow several incompatible configurations to be described as current in different chapters.

## 10. Validation after writeup changes

At minimum:

1. inspect the main project repository state/diff relevant to the code or release change;
2. inspect the Overleaf writeup repository state/diff separately;
3. search both TeX files for superseded terms and configuration values related to the update, including the names of removed stages and retired state variables;
4. confirm the release-model contract in section 6 still holds in both documents;
5. compile both documents in the required build directory and confirm there are no undefined references;
6. inspect the resulting writeup diff for accidental unrelated rewrites;
7. report what was updated, retired, and intentionally preserved as history.

Compile with:

```bash
cd docs/writeup-overleaf
latexmk -xelatex -output-directory=latex-build paper.tex
latexmk -xelatex -output-directory=latex-build main.tex
```

Keep all LaTeX auxiliary and generated files under `docs/writeup-overleaf/latex-build/`.

## 11. Scope discipline

A catch-up is not permission to rewrite the whole writeup after every code change. Prefer the smallest coherent update that restores accuracy and internal consistency. Do not refactor unrelated prose merely because the document is open.

When an active release fact changes, however, search globally enough to remove or reclassify every contradictory active statement. Stop once the current description is correct, the backward sweep is clean, both documents compile, and material limitations are stated accurately.
