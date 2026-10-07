# Visualization catalog — educational animations

Self-contained teaching animations of the face Riemann problem. See [README.md](README.md) for conventions.

---

## Educational / tutorial animations (`visualization/godunov_tutorial/`, `visualization/roe_tutorial/`)

Self-contained teaching material, not production output. Plain 1D ideal-gas Euler
(`gamma = 1.4`, Sod data) — no gravity, ionization, conduction, source terms or well
balancing. Two companion pieces in the same visual language, on the same initial data,
meant to be watched back to back: `godunov_tutorial/` (the **exact** Riemann solver) and
`roe_tutorial/` (the **linearised** one). Scripts and outputs live together in each
directory, alongside a `README.md` and a printed numerical summary.

### `godunov_face_flux.mp4`, `godunov_face_flux.gif`
Nine-scene explainer of `U_L, U_R -> exact Riemann solution -> U(x/t = 0) -> F_face`:
two cells / one face, MUSCL-reconstructed face states, the local Riemann problem,
the unfolding fan (with tracer particles advected by the exact `u(x/t)`), the `x`-`t`
self-similarity picture, sampling at `x/t = 0`, the physical Euler flux, the
finite-volume update, and a closing exact-Godunov vs Roe comparison.
1920x1080, 25 fps, 116 s; frames rendered in parallel then encoded with ffmpeg.
```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
  visualization/godunov_tutorial/make_godunov_animation.py --jobs 12
```
Options: `--probe` (still frames per scene for layout checks, into `probe/`),
`--no-gif`, `--fps N`, `--out PATH`, `--gif PATH`, `--keep-frames`.

### `riemann_summary.txt`
Numerical summary that drives the animation: `p*`, `u*`, wave types and speeds, the
state sampled at `x/t = 0`, and the resulting mass/momentum/energy fluxes, with
`F(U_L)`, `F(U_R)`, both averages and the Roe flux for contrast.
```bash
.venv/bin/python visualization/godunov_tutorial/exact_riemann.py \
  > visualization/godunov_tutorial/riemann_summary.txt
```

### `roe_face_flux.mp4`, `roe_face_flux.gif`
Twelve-scene explainer of `U_L, U_R -> Roe average -> linearised characteristic waves
-> F_Roe`: the same face, why the nonlinear flux is linearised, the **secant-not-tangent**
idea leading to `F_R - F_L = A~ (U_R - U_L)`, the `sqrt(rho)`-weighted Roe average (with
`(U_L+U_R)/2` explicitly crossed out), the three linear waves, the waterfall decomposition
`dU = sum_k alpha_k r_k`, the `x`-`t` diagram and eigenvalue-sign upwinding, the two halves
of the flux formula, the equivalent upwind walk sampled at `x/t = 0`, a side-by-side with
exact Godunov, when the linearisation is accurate, and where it sits in the solver step.
1920x1080, 25 fps, 144 s; frames rendered in parallel then encoded with ffmpeg. Imports
`solve_riemann` from `../godunov_tutorial/exact_riemann.py` for the contrast scenes.
```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
  visualization/roe_tutorial/make_roe_animation.py --jobs 12
```
Options: same as above (`--probe`, `--no-gif`, `--fps N`, `--out PATH`, `--gif PATH`,
`--keep-frames`).

### `roe_summary.txt`
Numerical summary and verification residuals that drive the animation: `U_L/U_R`,
`F_L/F_R`, `u~`, `H~`, `a~`, the three eigenvalues, right eigenvectors, wave strengths
`alpha_k`, each `alpha_k r_k` and `|lam_k| alpha_k r_k`, the intermediate states, the Roe
matrix, the flux in three equivalent forms, and residuals for `dU = sum alpha_k r_k`,
`dF = A~ dU`, `A~ = A(u~,H~)`, and consistency. Includes the two non-degenerate checks
(moving-gas variant and a 200-case random sweep) that the symmetric Sod problem cannot
provide, plus an exact-Godunov comparison.
```bash
.venv/bin/python visualization/roe_tutorial/roe_solver.py \
  > visualization/roe_tutorial/roe_summary.txt
```

---
