# Face-local mesh-scaled numerical conduction

## Decision

Implemented and validated. Artificial conduction now follows the actual local
centre-to-centre face spacing without changing physical conductivity or boundary
conditions.

## Structural change

- `Grid::numerical_diffusivity_per_length` stores the mesh-independent
  `C_num` in m/s. For `model_column`,
  `C_num = (ISO_CORONA ? 4 : 1) * 2000 * ISO_NUMERICAL_DIFFUSIVITY_MULT`.
- `physics.hpp::numerical_diffusivity_at_face` is the authoritative definition:
  `chi_num,f = C_num * ds_face` and `K_num,f = chi_num,f * C_V,f`.
- The gamma-table nonlinear solve, its residual checker, and the legacy charged
  and neutral implicit conduction rows all use `ds_imh_i(i)` / `ds_iph_i(i)`.
  Mirrored boundary geometry therefore needs no special case.
- The default-off outer capture records the accepted `chi_num_face`,
  `kappa_num_face`, `q_num`, and `q_total` from the final converged gamma solve.
- Other scenarios were migrated without changing their uniform-grid intent:
  their former `2000*ds_mean` assignments are now `C_num=2000 m/s`, with existing
  scenario multipliers applied to `C_num`.

## Targeted 500 s result

Same 1600--2153 km gamma/Saha, C7/pchip, hydro-temperature-decoupled wall,
fixed back pressure, 22 kK conduction wall, cooling/TRAC off, MC3 beta 2,
CFL 0.25, and equilibrium-reference well balancing as the prior study.

| metric at 500 s | U, local | old R4 coarse-scaled | old R4 global 0.25 | R4-local |
|---|---:|---:|---:|---:|
| cells | 1000 | 1320 | 1320 | 1320 |
| coarse / fine width [m] | 553 / 553 | 553 / 138 | 553 / 138 | 553 / 138 |
| outer `chi_num` [m2/s] | 1.106e6 | 1.106e6 | 2.765e5 | **2.760e5** |
| outer `kappa_phys` [W m-1 K-1] | 0.6090 | -- | -- | **0.6517** |
| outer `kappa_num` [W m-1 K-1] | 0.8314 | -- | -- | **0.1956** |
| `q_num/q_total` | 0.5772 | 0.5516 | 0.2311 | **0.2308** |
| mean effective mass flux [kg m-2 s-1] | 2.694e-9 | 3.640e-9 | 1.072e-9 | **1.057e-9** |
| `A_M` | 0.13658 | 0.00943 | 0.03071 | **0.02954** |
| `Q_M` | 0.08438 | 0.00308 | 0.00677 | **0.00684** |
| face-flux `A_eff` / `Q_eff` | 0.01088 / 1.47e-4 | 0.00089 / 7.00e-5 | 0.02031 / 1.92e-4 | **0.01804 / 2.00e-4** |
| steps; cells x steps | 86,434; 8.64e7 | 355,731; 4.70e8 | 356,051; 4.70e8 | **356,052; 4.70e8** |

The new uniform run reproduces the saved old U metrics exactly, including step
count, conduction split, mass flux, and ripple metrics. At the R4 outer face the
local diffusivity is 0.2496 of the uniform/coarse value, as required. The old
refinement-induced evaporation increase is gone: R4-local agrees within 1.5% with
the independent global-0.25 diagnostic and is 71% below old coarse-scaled R4.

The graded band has no new interface-localized artifact. With reduced artificial
conduction its small mass-flux variation is not strictly monotone, but its
normalized curvature is `3.79e-4`, essentially the same as the old global-0.25 R4
run (`3.92e-4`), with fewer slope-sign changes. The accepted face flux in the TR
window remains smooth (`A_eff=0.0180`, `Q_eff=2.00e-4`).

## Verification

- `chromo_tests`: 8,145 checks passed, 0 failed.
- Full CTest: 8/8 passed.
- U and R4-local both completed to 500 s with finite state and no conduction
  Newton or mesh-transition failure.

## Previous next step (completed below)

Run one higher coarse-equivalent resolution with the same local formulation to
measure whether the remaining mean evaporation flux is grid-converging; do not
tune `C_num` first.

## Higher-resolution follow-up: R4-local N=2000

The requested single higher-resolution run kept the same 1600--2153 km domain,
outer refined band, gamma/Saha physics, boundary conditions, limiter, CFL,
well-balancing, physical conductivity, and `C_num=2000 m/s`. It completed to
500 s without a nonfinite state or conduction-Newton failure.

| metric at 500 s | R4-local N=1000 | R4-local N=2000 | relative change |
|---|---:|---:|---:|
| total cells | 1320 | 2638 | -- |
| coarse / fine / outer width [m] | 552.84 / 138.02 / 138.02 | 276.45 / 69.10 / 69.10 | about 1/2 |
| completed steps | 356,052 | 712,320 | 2.00x |
| cells x steps | 4.700e8 | 1.879e9 | 4.00x |
| outer `chi_num` [m2/s] | 2.7604e5 | 1.3820e5 | -49.9% |
| outer `kappa_phys` [W m-1 K-1] | 0.6517 | 0.6574 | +0.9% |
| outer `kappa_num` [W m-1 K-1] | 0.1956 | 0.09725 | -50.3% |
| `q_phys` [W m-2] | 1.1359 | 0.8544 | -24.8% |
| `q_num` [W m-2] | 0.3409 | 0.1264 | -62.9% |
| `q_total` [W m-2] | 1.4768 | 0.9808 | -33.6% |
| `q_num/q_total` | 0.2308 | 0.1289 | -44.2% |
| mean cell `rho V` [kg m-2 s-1] | 1.0566e-9 | 6.0980e-10 | -42.3% |
| mean `F_eff` [kg m-2 s-1] | 1.0566e-9 | 6.1003e-10 | **-42.3%** |
| `A_M` | 0.02954 | 0.01237 | -58.1% |
| `Q_M` | 0.00684 | 0.00198 | -71.0% |
| `A_mismatch` | 0.02504 | 0.01010 | -59.6% |
| `A_eff` | 0.01804 | 0.00842 | -53.3% |
| `Q_eff` | 2.00e-4 | 4.21e-4 | +110% |

**Decision: Case B.** The incorrect coarse-scaled artificial-conduction trend
remains removed, and the run is stable, but a 42% change in mean evaporation flux
is not evidence of convergence. The physical/conduction solution remains
under-resolved. The outer physical conductivity itself changes by only 0.9%; most
of the `q_phys` change comes from the resolved wall-temperature gradient, while
halving `chi_num` also halves `kappa_num` and reduces `q_num` more strongly.

Ripple and mismatch amplitudes improve substantially. Face-flux amplitude
smoothness (`A_eff`) also improves, but normalized curvature (`Q_eff`) worsens.
No distinct refinement-interface failure appears: the run stays finite, and the
grade-band curvature remains the same order as the adjacent bands, though its
fine-scale roughness is higher and should not be described as strictly monotone.

The conservation coverage gap was closed in the existing
`test_irregular_conduction_conserves_energy`: it now verifies that nonzero
face-local artificial conduction changes the graded-mesh RHS and that the total
zero-boundary-flux energy contribution still telescopes. `chromo_tests` and all
8 CTest targets pass.

**Next step:** if another convergence point is warranted, run one R4-local
N=4000 case with the formulation and coefficient unchanged; do not tune `C_num`.
