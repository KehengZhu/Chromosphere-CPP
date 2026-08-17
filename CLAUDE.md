# Chromosphere2026 — Agent Instructions

## Writeup maintenance

`docs/writeup-overleaf/` is a separate Overleaf-linked Git repository, not part of the main repository's Git history. Before auditing, catching up, cleaning, or substantially editing `paper.tex` or `main.tex`, read the authoritative main-repository standard in `docs/design/writeup_standards.md` and inspect the two Git worktrees separately.

Changes to release physics, governing equations, thermodynamics, numerics, boundary conditions, active update paths, canonical defaults/configuration, validation, limitations, or scientific interpretation trigger a writeup catch-up review. Current production code is the source of truth for active implementation behavior; existing TeX and older recaps are not.

Keep `paper.tex` concise and current-release-only. Use `main.tex` for the deeper engineering record and clearly label inactive, experimental, historical, retired, and planned material. Update only the document or documents whose role is affected, avoid configuration mixing, and do not commit or push either repository unless explicitly requested.

The compiled PDFs are mirrored into the main repository automatically: `main.pdf` → `docs/chromosphere_writeup.pdf` and `paper.pdf` → `docs/chromosphere_paper.pdf`. `docs/writeup-overleaf/latex-build/` is gitignored in both repositories, so those two copies are the only tracked PDFs of the writeup, and they are the ones to link from `README.md` and the Doxygen pages. The sync is a `$success_cmd` hook in `docs/writeup-overleaf/.latexmkrc` that calls `scripts/sync_writeup_pdf.sh`, so every successful `latexmk` run — VS Code LaTeX Workshop recipe or command line — refreshes them; the hook no-ops where the main repository is absent, e.g. when Overleaf compiles the same sources. Do not delete the hook when editing `.latexmkrc`, and run `scripts/sync_writeup_pdf.sh` by hand if a PDF was produced some other way.

## Code documentation

The API reference is Doxygen: config in `docs/doxygen/Doxyfile`, narrative pages in `docs/doxygen/pages/`, built with `scripts/build_docs.sh`. `README.md` is deliberately a quick start only — configure, build, choose parameters, run, visualize — and technical detail belongs in the Doxygen pages, not in the README. When you change public behavior, update the relevant page (`configuration.md`, `numerics.md`, `physics_model.md`, `scenario_reference.md`, `io_formats.md`, `validation.md`, `architecture.md`) in the same change.

**`docs/doxygen/html/` is a gitignored build artifact, not tracked.** The browsable copy is published **publicly, with no authentication** to `https://kehengphysics.site/docs/chromosphere/` by `scripts/publish_docs.sh` (rsync over ssh into `/srv/www/kehengphysics.site/docs/chromosphere/`, served by `/etc/nginx/snippets/kehengphysics-docs.conf`, documented by a `README.md` at `/srv/www/kehengphysics.site/`). That is the only published copy and the one `README.md` links, so **a stale site is a stale reference for everyone** — nothing in the repository compensates for forgetting to publish. `publish_docs.sh` builds first if the artifact is missing, so it works from a fresh clone.

### Refresh the reference whenever you change what it documents

This is a required step of the change, not a follow-up. Any change to code, public behavior, configuration, numerics, physics, scenarios, I/O formats, validation, or the writeup triggers all three of:

1. **Update the hand-written narrative pages** in `docs/doxygen/pages/` — `configuration.md`, `numerics.md`, `physics_model.md`, `scenario_reference.md`, `io_formats.md`, `validation.md`, `architecture.md`. Doxygen regenerates the API listings from the sources automatically, but it cannot update prose, so **this is the step that actually goes stale.** Rebuilding alone is not enough. If a change alters the writeup (`paper.tex` / `main.tex`), check whether the same fact is asserted in a narrative page and in the two mirrored PDFs (`docs/chromosphere_paper.pdf`, `docs/chromosphere_writeup.pdf`) that the pages link.
2. **Rebuild** with `scripts/build_docs.sh`, warning-free.
3. **Publish** with `scripts/publish_docs.sh`.

`.githooks/pre-commit` rebuilds and fails the commit on any Doxygen warning, but it deliberately does **not** publish — pushing to a public website is not something a commit should do silently — and it cannot tell whether the prose is current. Publishing and prose are on you. The hook also only exists in clones that have run `scripts/install_hooks.sh` once, so never assume it is active.

**That public site is only safe because the Doxyfile sets `SOURCE_BROWSER = NO` and `VERBATIM_HEADERS = NO`.** This repository is private; those two settings are what keep the generated tree free of verbatim source. Do not turn either on — and do not add any other setting that emits code (`INLINE_SOURCES`, `REFERENCES_LINK_SOURCE`) — without first deciding that publishing the entire source tree to the open internet is acceptable. `scripts/build_docs.sh` wipes `docs/doxygen/html/` before each run for the same reason: Doxygen overwrites but never deletes, so output it no longer generates would otherwise linger and still be published.

Keep the build warning-free. `EXTRACT_ALL` is off and `WARN_IF_UNDOCUMENTED` is on, so a new public entity with no documentation is a warning in `docs/doxygen/doxygen-warnings.log` and the pre-commit hook rejects the commit. Every public entity gets a real statement of what it is, with physical units in square brackets where it has them — not a restatement of its name.

`scripts/build_docs.sh` wipes `docs/doxygen/html/` before every run. This matters: Doxygen overwrites files but never deletes them, so output it no longer generates would otherwise linger and still be published. That is exactly how disabling `SOURCE_BROWSER` once left 34 stale source listings in place, ready to ship.

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

## `model_column` release architecture and numerical method

The release is a **single-fluid field-aligned equilibrium-mixture solver** whose conserved state is exactly `U = (rho, rho u, E)` (`mix::RHO/MOM/ENERGY`, `src/single_fluid/`). The ionization fraction `x`, the carrier densities `x*rho` / `(1-x)*rho`, `n_e`, `n_HI` and `p_e` are **derived** EOS diagnostics, never stored or advanced, so the release has no carrier-row storage and no equilibrium-projection stage. One step is
`U^n -> explicit MUSCL-Hancock/Godunov hydro -> U* -> implicit backward-Euler physical conduction -> U^{n+1}` — those two stages are the whole `mixture_advance`, composed by first-order Lie splitting (not Strang). TRAC, beam heating, volumetric coronal heating, radiative cooling and artificial/numerical conduction are **not** release physics and are not present in that path, not even behind default-off flags; they live only in the two-fluid solver.

A normal `model_column` run uses Gamma/Saha equilibrium thermodynamics (ionization energy included, production `Gamma1` table), MUSCL reconstruction of **`(ln rho, u, ln p)`** with the MC3/Koren limiter (β = 2), a **MUSCL-Hancock** half-step predictor that applies the same momentum source (flux-tube pressure term + gravity) as the corrector, the **SWMF-style exact-Riemann Godunov numerical flux** (face-local Rusanov fallback) evaluated on the time-centred face states, physical conduction only, and a 22,000 K external conductive reservoir at the physical outer face, on the coarse-equivalent N = 500 / R4 mesh (661 actual cells) at CFL 0.50.

**The release is reference-free: the interior hydrodynamic operator needs no frozen global hydrostatic reference state.** (The boundary closures still capture fixed reservoir data once from the IC — `kOuterPRef`, `kInnerTRef`/`kInnerPRes` — as any truncated-domain problem must. The release is *not* a well-balanced scheme in the constructive sense; it removes the dominant discretization inconsistencies and leaves a small, resolution-dependent residual.) `Grid::eq_wb` / `eq_state` / `eq_residual` and `ISO_EQ_WB` were retired from `model_column` once the predictor source term and a second-order (trapezoidal, EOS-closed) lower ghost ladder cut the raw discrete hydrostatic face mass-flux defect by ~560x, down to a level comparable to the estimated float32 round-off scale of the stored state at release resolution (no clean truncation-error trend survives between N=500 and N=1000). `eq_wb` still exists in `chromosphere.hpp` for the **two-fluid research solver only** — never set it on a release Grid. Evidence: `docs/studies/numerics/reference_free_release_recap.md`.

The **lower boundary at 1600 km is a mid-chromospheric truncation of the modeled domain, not the photosphere.** It is the height below which the Saha/LTE closure is not intended to be trusted; the atmosphere below it is outside the model, and the two inner ghosts are a numerical closure for the MUSCL stencil, not a claim that LTE extends downward. Its physical content is a quasi-static stratified reservoir at the truncation temperature with `v = 0`. Never describe it as photospheric.

The historical **two-fluid / three-temperature / finite-rate-ionization** solver is a genuinely separate non-release path (`src/two_fluid/`: `two_fluid.hpp`, `state.cpp`, `flux.cpp`, `rhs.cpp`, `integrators.cpp`) that keeps its own seven-row state and source stages. The source tree mirrors this: `src/single_fluid/` and `src/two_fluid/` share only `chromosphere.hpp` (the Grid) and `eos.hpp`, and neither includes the other. Loading a `Gamma1` table selects the release solver; each solver rejects the other's `Grid`. A release run also rejects `SINGLE_FLUID`, `ENABLE_TE`, `ISO_TWO_FLUID` and `ISO_IONIZATION` outright — a release user never needs them. Scenarios `model_gentle`, `model_c7`, `model_flare`, `analytic_canopy` and `pfss_field_line` use the two-fluid solver.

**`model_column` means exactly one thing: the release solver.** The scenario supplies the production `Gamma1` table unconditionally and rejects `ISO_GAMMA`; `model_column_ic` rejects every historical two-fluid knob (`ISO_GAMMA`, `ISO_TWO_FLUID`, `ISO_IONIZATION`, `ISO_COOLING`, `ISO_TRAC`, `ISO_CORONA`, `ISO_CHEAT`, `ISO_QFLUX`, `ISO_TBOOST`, `ISO_NUMERICAL_DIFFUSIVITY_MULT`). The `model_isentropic` alias has been retired; the historical two-fluid column is `model_gentle`.

Rejected and **not** to be reintroduced without new evidence: the full SWMF `ModTransitionRegion` momentum pressure/gravity split (measurably no better than the predictor source fix alone on this constant-area column, and not generally conservative for variable-area flux-tube geometry), and a gravity-aware characteristic (non-reflecting) lower boundary (changed the hydrostatic residual not at all and the evaporation velocity by <= 0.4 %; the 1600 km face is not an important acoustic reflector for this model).

**The release face Riemann solver is the SWMF-style exact-Riemann Godunov flux at FROZEN COMPOSITION.** It solves the local Riemann problem exactly as an ideal gas at a fixed `gamma = 5/3`, carrying the Saha ionization energy (and this model's gravitational potential) in SWMF's passive specific-energy offset `E0`, and samples the self-similar solution at `x/t = 0`. The physical basis is that the equilibrium Saha internal energy splits exactly as `e_int = p/(5/3-1) + e_ion`, so the translational gas is monatomic and the ionization reservoir is inert *across the wave interaction* — the right limit when the ionization/recombination relaxation time is long compared with the dynamical time, which is the chromospheric case for hydrogen. Every face that cannot be solved exactly falls back to face-local Rusanov and is counted in `Grid::godunov_stats`; the driver prints a `godunov.*` tally at the end of every release run and `scripts/release_validation.sh` fails on a nonzero fallback count.

**The equilibrium acoustic index `Gamma1` is NOT "wrong" — it is the other limit.** `Gamma1` (median 1.090 in this column) is the correct acoustic response if Saha equilibrium is re-established instantaneously *within* the wave; frozen `5/3` is the correct response when ionization cannot relax on the wave timescale. The release takes the frozen limit, on that physical argument and on the SWMF methodology. Never write that `Gamma1` is universally wrong. `Gamma1` remains the closure's equilibrium derivative and is still used by the `.gamma_diag` diagnostics, the outer-boundary Mach cap, and the Roe reference flux.

**The timestep follows the flux.** `mixture_timestep` sizes the step from `|u| + sqrt(max(Gamma1, 5/3) p/rho)` when the Godunov flux is in force and from `|u| + sqrt(Gamma1 p/rho)` under the Roe/Rusanov reference fluxes, so `CHROMO_CFL` always means the CFL of the scheme actually running. On the canonical N=500/R4 column this is a measured no-op — the CFL-limiting cell is the fully ionized top cell where `Gamma1 = 5/3` — so `dt`, the step count and the end time are unchanged. A previously documented "effective CFL 0.62" for this flux was wrong; the measured value is 0.500.

**Open modelling question, not to be quietly dropped:** the Riemann solve freezes composition, but the EOS decode re-imposes instantaneous Saha equilibrium after every timestep. Those are separate choices; the algorithm is effectively `Saha-equilibrated state -> frozen-composition hydro step -> Saha-equilibrated state`, an operator-split *instantaneous-relaxation* approximation that is not justified when the physical relaxation time is long. Documented in `main.tex`; do not redesign the ionization model without asking.

**No environment variable is required to select the release solver, its flux, or its reconstruction.** `ISO_RIEMANN=roe-local` (the previous release flux, the equilibrium-`Gamma1` mixture Roe linearization), `ISO_RIEMANN=rusanov` and `ISO_RECONSTRUCTION=lnrho-v-lnt` are reference-only overrides kept for regression and controlled numerical comparison; never use them in production. Full description and evidence: `docs/studies/numerics/swmf_godunov_flux_experiment.md` (flux and CFL) and `docs/studies/numerics/model_column_release_numerics_recap.md` (reconstruction). Release checks: `scripts/release_validation.sh`.

Known limitation to preserve in any writeup: **N500 long-duration evaporation mass flux is not established as < 10 % grid-converged.** Distinguish the release numerical method from the resolution required for a particular quantitative scientific claim.

## Validated production CFL and long-run controls

`CHROMO_CFL=0.50` is the **validated production runtime default for the reduced `model_column` release configuration: coarse-equivalent N=500, R4 outer refinement, 661 actual cells, physical-face 22,000 K conductive boundary, and physical conduction only**. Evidence: `docs/studies/conduction/coarse_model_column_physical_conduction_recap.md`. That sweep was run under the Roe flux and the equilibrium signal speed; it carries over to the release Godunov flux because the flux-consistent CFL rule leaves `dt` and the step count unchanged on this exact configuration (measured, `docs/studies/numerics/swmf_godunov_flux_experiment.md`). It does not carry over automatically to a different mesh or model shape, where a partial-ionization cell could become CFL-limiting and the step would then be up to 1.24x shorter.

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
