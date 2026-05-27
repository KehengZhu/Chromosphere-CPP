# Chromosphere2026 — Agent Instructions

## Writeup location

The project writeup lives in `docs/668746597ab9c90ccebae1a0/`. The primary source files are `main.tex` and `paper.tex`, with bibliography in `reference.bib`.

## Compiling the writeup

When compiling LaTeX, always build inside `docs/668746597ab9c90ccebae1a0/latex-build/` — never let auxiliary files (`.aux`, `.log`, `.fls`, `.pdf`, etc.) land next to the `.tex` sources.

Example:
```bash
cd docs/668746597ab9c90ccebae1a0
latexmk -xelatex -output-directory=latex-build main.tex
```

## Output directories

- **Code outputs** (simulation results, data dumps, logs from `chromo_main` / scenarios): save under `outputs/`.
- **Visualizations** (plots, animations, figures generated from outputs): save under `util/visualization/`.

Do not scatter output files in the repo root or alongside source code.
