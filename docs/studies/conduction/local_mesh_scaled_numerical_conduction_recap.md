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

## Third convergence point: R4-local N=4000 (OpenMP)

The requested N=4000 point was run with the formulation, coefficient
(`C_num = 2000 m/s`), mesh profile, physics, boundary conditions, limiter, CFL,
and well balancing unchanged. Only `ISO_NS` and the diagnostic strides differ
from N=2000. It used the committed OpenMP build at 12 threads; that path is
already established byte-identical to serial, so no equivalence run was repeated.
`ctest --test-dir build_omp` with `OMP_NUM_THREADS=12`: **8/8 passed**.

The run completed normally: `END! step=1397913 time=500`, no nonfinite state, no
conduction-Newton failure, no EOS-inversion or mesh-transition exception. The
timestep was flat at `dt ~ 3.5097e-4 s` throughout.

### Mesh and cost

| | N=1000 | N=2000 | N=4000 |
|---|---:|---:|---:|
| total cells | 1320 | 2638 | **5275** |
| coarse width [m] | 552.8 | 276.4 | **138.2** |
| fine width [m] | 138.02 | 69.10 | **34.55** |
| outer face centre-to-ghost `ds_face` [m] | 138.02 | 69.10 | **34.55** |
| cells in 2130--2150 km | 145 | 290 | **579** |
| completed steps | 356,052 | 712,320 | **1,397,913** |
| cells x steps | 4.700e8 | 1.879e9 | **7.374e9** |
| relative cost | 1.00x | 4.00x | **15.69x** |

### Outer conduction

| metric at 500 s | N=1000 | N=2000 | N=4000 | 2000->4000 |
|---|---:|---:|---:|---:|
| `T_top` [K] | 21 759.4 | 21 910.2 | **21 961.3** | +0.23% |
| `T_wall - T_top` [K] | 240.56 | 89.81 | **38.74** | -56.9% |
| `G_wall = (T_wall-T_top)/ds_face` [K/m] | 1.7430 | 1.2996 | **1.1213** | **-13.7%** |
| `kappa_phys` [W m-1 K-1] | 0.65172 | 0.65742 | **0.65935** | +0.29% |
| `chi_num` [m2/s] | 2.7604e5 | 1.3820e5 | **6.9100e4** | **-50.0%** |
| `kappa_num` [W m-1 K-1] | 0.19559 | 0.097249 | **0.048523** | -50.1% |
| `q_phys` [W m-2] | 1.1359 | 0.85441 | **0.73935** | **-13.5%** |
| `q_num` [W m-2] | 0.34090 | 0.12639 | **0.054411** | -57.0% |
| `q_total` [W m-2] | 1.4768 | 0.98080 | **0.79376** | -19.1% |
| `q_num/q_total` | 0.23083 | 0.12886 | **0.068548** | **-46.8%** |

`chi_num` halves exactly with the fine spacing (-49.93%, -50.00%), which is the
face-local formulation behaving as specified. The numerical share of the outer
conductive input is now under 7%.

### Evaporation and ripple

| metric at 500 s | N=1000 | N=2000 | N=4000 | 2000->4000 |
|---|---:|---:|---:|---:|
| mean cell `rho V` [kg m-2 s-1] | 1.0566e-9 | 6.0980e-10 | **6.2240e-10** | +2.07% |
| mean `F_eff` [kg m-2 s-1] | 1.05663e-9 | 6.10035e-10 | **6.22495e-10** | **+2.04%** |
| `A_M` | 0.029538 | 0.012375 | **0.015186** | +22.7% |
| `Q_M` | 0.0068407 | 0.0019825 | **0.0020046** | +1.1% |
| `A_mismatch` | 0.025036 | 0.010102 | **0.0025292** | **-75.0%** |
| `A_eff` | 0.018039 | 0.0084246 | **0.015223** | +80.7% |
| `Q_eff` | 2.0007e-4 | 4.2089e-4 | **3.8650e-4** | -8.2% |
| `corr(dcen,ddiff)` | -0.8124 | -0.7563 | **-0.1820** | -- |

### Convergence

```text
delta_12 = (F_2000 - F_1000)/F_1000 = -42.27%
delta_23 = (F_4000 - F_2000)/F_2000 =  +2.04%

D_12 = |F_1000 - F_2000| = 4.4660e-10 kg m-2 s-1
D_23 = |F_2000 - F_4000| = 1.2460e-11 kg m-2 s-1
D_12 / D_23 = 35.8
```

`D_23 < D_12` by a factor of **35.8**. The sequence is however **not monotone**
(down 42%, then up 2%), so no apparent order and no Richardson limit is quoted
for `F_eff`. Forcing one here would be meaningless.

The outer boundary quantities *are* monotone, with mutually consistent apparent
orders on the 2x sequence:

```text
G_wall   ratio 2.486   p_app = 1.31
q_phys   ratio 2.447   p_app = 1.29
q_total  ratio 2.652   p_app = 1.41
```

Provisional Richardson limits from those three only (diagnostic estimates, not
validated physical values): `G_wall -> 1.00 K/m`, `q_phys -> 0.66 W m-2`,
`q_total -> 0.68 W m-2`. So the conductive drive is converging at roughly first
order but is **not yet converged**: `q_total` still moves -19% over the last
refinement while `F_eff` moves only +2%.

### Outer thermal boundary layer

One definition, applied identically at all three resolutions (new
`outer_boundary_layer` helper in `util/outer_refine_diag.py`): with
`G(h) = |dT/ds|` evaluated on the actual non-uniform cell centres, the layer is
the contiguous run of top cells with `G >= 0.1 * max(G)`, and its thickness is
measured from that edge to the outer cell centre.

| | N=1000 | N=2000 | N=4000 |
|---|---:|---:|---:|
| layer thickness [km] | 8.69 | 9.88 | **10.33** |
| cells across the layer | 64 | 144 | **300** |
| layer edge / `T_edge` | 2144.30 km / 7619 K | 2143.12 km / 7507 K | **2142.64 km / 7462 K** |
| `max|dT/ds|` [K/m] | 2.2637 | 2.1410 | **2.3917** |
| height of `max|dT/ds|` [km] | 2149.40 | 2146.72 | **2145.30** |
| `|dT/ds|` in the top cell [K/m] | 2.1486 | 1.7100 | **1.2425** |
| peak gradient at the outer cell? | no | no | **no** |

The layer is **not** a one-cell or few-cell feature: it is a ~10 km structure
resolved by 300 fine cells at N=4000, and its thickness is itself converging
(+13.7% then +4.55%). The steepest gradient sits **inside** the resolved domain
and moves progressively away from the wall (2149.4 -> 2146.7 -> 2145.3 km) while
the top-cell gradient falls monotonically (2.15 -> 1.71 -> 1.24 K/m). The
Dirichlet wall is therefore no longer the location of the sharpest structure.

### Refinement-interface check

No artifact. Grade-band (2079--2101 km) normalized curvature is
`3.79e-4 / 6.91e-4 / 6.11e-4` -- same order at all three meshes and not growing.
Maximum adjacent width ratio in the capture region is 1.0071 at N=4000, inside
the 1.1 cap. The flux split self-consistency `f_central + f_diff = f_total`
holds to `9.9e-8`. The captured `rho V` profile has the same shape and magnitude
at N=2000 and N=4000 across the whole capture region (5.6--6.3e-10 vs
5.98--6.37e-10 kg m-2 s-1); N=4000 is simply flatter with height.

### Reading of the ripple metrics

`A_mismatch` falls 75% and `corr(dcen,ddiff)` collapses from -0.756 to -0.182,
so the cell-centred `rho V` has now essentially converged onto the conserved face
flux and the anti-correlated central/diffusive decomposition artifact is largely
gone. `A_M`, `A_eff` and `Q_eff` are normalized amplitudes: their mean stopped
shrinking between N=2000 and N=4000 while 579 (rather than 290) cells now sample
shorter resolved scales in the same physical window, so their small rise is a
sensitivity effect of the metric, not new error. It is not evidence that OpenMP
or local refinement introduced anything.

### Decision: Outcome A

> The local face-scaled formulation shows evidence of approaching grid
> convergence, but the limiting evaporation rate remains provisional.

Formal convergence is not claimed from three points, and the mean evaporation
flux must still not be quoted as a physical result. The remaining caveat is that
`q_total` is still moving -19% per refinement while `F_eff` is flat to 2%, so the
near-flatness of `F_eff` at this pair is not yet backed by a converged drive.

**Next step:** do not run N=8000. The cheapest remaining lever is the outer
conductive discretization itself -- the apparent order of `G_wall`, `q_phys` and
`q_total` is ~1.3, consistent with a first-order-dominated ghost-centre Dirichlet
stencil, and `T_wall - T_top` is still shrinking by ~57% per refinement rather
than converging to a fixed boundary-layer value. Audit whether 22 000 K is
intended at the physical boundary face or at the ghost-cell centre before
spending another 4x on resolution. No boundary change was made in this stage.
