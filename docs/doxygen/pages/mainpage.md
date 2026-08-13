# Chromosphere2026 {#mainpage}

A 1.5D field-aligned solver for the solar chromosphere and the chromosphere-to-corona transition. C++14, header-only vendored Armadillo, no LAPACK/BLAS.

This is the code reference. For getting the model built and running, see the repository `README.md`; for the scientific write-up see `docs/writeup-overleaf/paper.tex` (current release) and `main.tex` (engineering record).

## Two solvers, one mesh

The project carries two solvers over one shared field-aligned mesh. They share only the @ref grid "Grid" and the @ref eos "equation of state"; neither includes the other, and they do not share a conserved-state width, a packing convention, or a timestep path.

| | @ref release_solver "Release (single fluid)" | @ref two_fluid_solver "Historical two-fluid" |
| --- | --- | --- |
| Source | `src/single_fluid/` | `src/two_fluid/` |
| Conserved state | 3 rows: `(rho, rho u, E)` | 7 rows: ion / neutral / electron carriers |
| Thermodynamics | Saha equilibrium + CRASH `Gamma1` | fixed adiabatic index |
| Ionization | equilibrium (derived) | finite-rate network |
| Temperatures | one | up to three (`T_e`, `T_i`, `T_n`) |
| Timestep | hydro, then implicit physical conduction | operator-split Stages A–E, R, beam, heating |
| Scenarios | `model_column` | `model_gentle`, `model_c7`, `model_flare`, `analytic_canopy`, `pfss_field_line` |
| Status | production | research / historical |

Which solver runs is decided by data, not by a flag: **a loaded `Gamma1` table selects the release solver**, and the `model_column` scenario supplies the production table unconditionally. Each solver rejects the other's Grid, so the choice cannot be made by accident.

## Where to start

| If you want to | Read |
| --- | --- |
| understand the equations being solved | @ref physics_model |
| understand the discretization | @ref numerics |
| pick a scenario | @ref scenario_reference |
| set a command-line argument or environment variable | @ref configuration |
| parse a snapshot or sidecar file | @ref io_formats |
| find your way around the source tree | @ref architecture |
| read a symbol off the write-up | @ref naming |
| know what is verified and what is not | @ref validation |

## Module index

- @ref release_solver — the production solver.
- @ref two_fluid_solver — the historical research solver.
- @ref shared — @ref grid, @ref eos, @ref physics_formulas, @ref runtime.
- @ref scenarios — initial conditions, boundary conditions, meshes, geometry.
- @ref driver — the CLI program and its run-control helpers.

## The release configuration in one paragraph

`model_column` is a field-aligned chromosphere-to-corona column on a single straight field line, initialized from the Model C7 atmosphere with the density re-integrated hydrostatically through the Saha closure. It runs the single-fluid equilibrium-mixture solver with Gamma/Saha thermodynamics (ionization energy included), MUSCL reconstruction of `(ln rho, u, ln p)` with the MC3/Koren limiter, a MUSCL-Hancock predictor, the mixture Roe characteristic flux, and implicit physical conduction driven by a 22,000 K external conductive reservoir at the physical outer face — on a coarse-equivalent `N = 500` / R4 mesh (661 actual cells) at `CFL = 0.50`. No environment variable is needed to select any of it; the scenario does. See @ref numerics.

## Primary references

- Al Shidi, Q., et al. (2019), *Time-Dependent Two-Fluid Magnetohydrodynamic Model of the Solar Chromosphere* — the 2D two-fluid model this code is the field-aligned reduction of.
- Avrett, E. H. & Loeser, R. (2008), *ApJS* **175**, 229 — Model C7, Table 26.
- Rosner, Tucker & Vaiana (1978); Antiochos & Sturrock (1978) — static loop balance and gentle evaporation.
- Fisher, Canfield & McClymont (1985), *ApJ* **289**, 425 — explosive evaporation.
- Johnston et al. (2019, 2020) — TRAC.
- Voronov (1997); Hummer (1994); Carlsson & Leenaarts (2012); Chae (2021) — rate coefficients and radiative losses.

Full PDFs are in `docs/supporting-papers/`.
