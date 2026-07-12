# Chromosphere2026 — Agent Instructions

## Writeup location

The project writeup lives in `docs/writeup-overleaf/`. The primary source files are `main.tex` and `paper.tex`, with bibliography in `reference.bib`. `main.tex` is for operational use and contains technical details, while `paper.tex` is a cleaner version with publication style. When asked to update the writeup, you will need to update both `main.tex` and `paper.tex`.

## Compiling the writeup

When compiling LaTeX, always build inside `docs/writeup-overleaf/latex-build/` — never let auxiliary files (`.aux`, `.log`, `.fls`, `.pdf`, etc.) land next to the `.tex` sources.

Example:
```bash
cd docs/writeup-overleaf
latexmk -xelatex -output-directory=latex-build main.tex
```

## Output directories

- **Code outputs** (simulation results, data dumps, logs from `chromo_main` / scenarios): save under `outputs/`.
- **Visualizations** (plots, animations, figures generated from outputs): save under `visualization/` (project root). This directory is gitignored. Plotting scripts live in `util/` and already default their output here — keep new/edited scripts defaulting to `visualization/` too. Never write figures or movies to the repo root, into `util/`, or alongside source code.

Do not scatter output files in the repo root or alongside source code.

## visualize_commands.md

`visualize_commands.md` (project root) is the catalog of the exact command used to generate every file in `visualization/`. **Whenever you generate a visualization, add or update its entry** in `visualize_commands.md`: the output filename, the exact command run (assume the project root as the working directory), and a one-line description. Keep it in sync — if you change a script's invocation or add a new figure/movie, update the catalog in the same change.

## Visualization
When asked to generate a movie, generate frames in parallel if possible.

When plotting field-aligned profiles that span the chromosphere through the corona (e.g. loop flare runs), use a **split x-axis**: zoom in on the thin chromosphere+TR (sub-Mm) and compress the extended corona (tens of Mm) using different x-scales (a continuous piecewise transform via matplotlib's `set_xscale("function", ...)` works well — see `util/animate_loop_flare.py`). Mark the chromosphere/TR/corona region boundaries with translucent vertical dashed lines.

For a **full closed loop** (footpoint A → apex → footpoint B, `topology=full`), plot vs **arc length** `s`, not height — on a closed loop height folds the two legs onto each other. Apply the split symmetrically: zoom **both** footpoints and compress the coronal middle, laying the loop out flat `[footpoint A | apex | footpoint B]`, with the apex marked by a dotted line and region boundaries on both legs. `util/animate_loop_flare.py` auto-detects topology from the `.dat` and switches between height (half loop) and arc-length (full loop) accordingly.

## Reference Papers
Reference Papers that you may need to understand a physical process or code implementation are in docs/supporting-papers. Check them first before you prompt the user to download the papers.