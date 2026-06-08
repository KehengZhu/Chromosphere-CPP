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
- **Visualizations** (plots, animations, figures generated from outputs): save under `util/visualization/`.

Do not scatter output files in the repo root or alongside source code.

## Visualization
When asked to generate a movie, generate frames in parallel if possible.
