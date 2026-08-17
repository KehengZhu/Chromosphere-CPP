# Outer (upper-TR) static local refinement: does it reduce the numerical ripple with less artificial conduction?

> Historical pre-local-scaling result. Its globally coarse-scaled R4 and R4
> multiplier-0.25 runs remain the comparison baselines. The follow-up implementation
> and validation are in `docs/studies/conduction/local_mesh_scaled_numerical_conduction_recap.md`.

Follow-up to `docs/studies/validation/resolution_convergence_scan_recap.md` (the `rho V` ripple in 2130–2150 km is spatial truncation error and converges at ~first order under *global* refinement) and `docs/studies/conduction/numerical_conduction_convergence_recap.md` (`numerical_diffusivity = 2000 Δh` supplies 58 % of the outer conductive drive at `ns=1000`). This stage adds a **static local** refinement profile that puts the fine cells only where the ripple lives — the upper TR against the outer boundary — and measures the same metrics.

No AMR, no regridding, no remapping. The mesh is built once at IC time by the existing shared builder.

## Headline

**Yes on both counts.** Outer refinement at fixed coarse-equivalent `ISO_NS = 1000` drops `A_M` from **0.1366 → 0.0514 → 0.0094** for uniform → ×2 → ×4, i.e. a **14.5× reduction** using 1320 cells instead of 1000. That beats *global* refinement on both axes: the uniform `ns=2000` run of the previous stage reached only `A_M = 0.0673` with 2000 cells and the same number of steps. And the ×4 mesh runs stably with `ISO_NUMERICAL_DIFFUSIVITY_MULT=0.25` for the full 500 s, with `A_M = 0.0307` — still 4.4× below the uniform run — while the numerical share of the outer conductive input falls from 55 % to 23 %.

**Caveat that must travel with the result:** refinement does **not** fix Result D. With `D_num` left at its coarse-equivalent value the mean evaporation flux *rises* (2.69 → 2.78 → 3.64 ×10⁻⁹ kg m⁻² s⁻¹), because the artificial conductivity is a single global scalar that local refinement does not reduce while the wall stencil `ds_face` shrinks with the local cell.

## 1. What was added

A second first-class profile in the shared static-mesh builder, selected by `ISO_REFINE_PROFILE` / `GRID_REFINE_PROFILE`:

| profile | layout (s measured from the inner face) |
|---|---|
| `lower` (default, historical) | coarse `[0,s_lo]` → fine `[s_lo,s_hi]` → grade fine→coarse `[s_hi,s_hi+τ]` → coarse `[s_hi+τ,L]` |
| **`outer` (new)** | coarse `[0,s_lo−τ]` → **grade coarse→fine** `[s_lo−τ,s_lo]` → **fine `[s_lo,L]` up to and including the outer face** |

For `outer`, `s_lo` is the foot of the fully refined region, the grade sits immediately *below* it, and `s_hi` is unused. Full semantics: `docs/studies/numerics/static_local_refinement.md`.

### A real bug the tight transition exposed, and the fix

The old grade chose its cell count as `ceil(max(n_continuity, n_ratio_cap))` and then stretched an exponential map to fill exactly the requested τ. Rounding `n` up destroys the endpoint width match, so the *junction* ratios are not the internal ratio. With the target profile (τ = 20 km, factor 4, `ISO_NS=1000` on a 553 km domain) that produced an adjacent width ratio of **1.172** at the coarse/grade junction — a silent violation of `MESH_MAX_RATIO = 1.1`. The generous `lower` profile (τ = 100 km) was far from the limit and never showed it.

The grade is now anchored the other way round: widths are `fine_ds·r^m` with `r = factor^{1/n} ≤ 1.1`, so the narrow end matches the fine spacing *exactly* and **every** adjacent ratio — internal and both junctions — is `r`. That fixes the width to `fine_ds·(factor−1)/(r−1)`, so `*_TRANSITION_KM` becomes a **minimum**: `n` is the smallest count whose grade is at least that wide. A width the cap cannot support is widened to the nearest one it can, instead of exceeding the cap. Every uniform segment now uses `ceil` so no segment is ever coarser than its target, which is what keeps the junction bounded. Side effect on the existing profile: the documented `lower` transition is 700 → **800.21** km instead of exactly 800 km. Nothing else about it changed.

## 2. Files changed

- `scenarios/mesh.hpp` — `RefineProfile` enum, `RefineParams::profile`, documented `outer` layout.
- `scenarios/mesh.cpp` — `env_profile`; `make_grade` / `append_grade` replace the old `append_transition`; `append_uniform` now takes a target spacing and `ceil`s; `build_static_mesh_faces` dispatches on the profile.
- `scenarios/model_column.cpp` — **comment only**. It already routed through `refine_params_from_env` + `build_static_mesh_faces` + `Grid::resize`, so selecting the profile needed no code.
- `tests/chromo_tests.cpp` — 3 new tests, 1 existing assertion updated (below).
- `docs/studies/numerics/static_local_refinement.md` — `outer` profile, the transition-as-minimum semantics, and an explicit note that refinement does *not* reduce `numerical_diffusivity` inside the fine band.
- `util/outer_refine_diag.py`, `visualize_commands.md`, this recap.

**No hydro, Rusanov, EOS, well-balancing, conduction stencil, boundary-condition or time-integration code was touched.**

## 3. Environment controls

```text
ISO_REFINE_PROFILE=outer        # or `lower` (default); GRID_REFINE_PROFILE is the generic name
ISO_REFINE_FACTOR=4             # fine/coarse ratio; 1 disables
ISO_REFINE_S_LO_KM=500          # foot of the fully refined region, from the inner face
ISO_REFINE_TRANSITION_KM=20     # MINIMUM grade width, immediately below s_lo
```

On the 1600–2153 km domain (`ISO_H_BASE=1600 ISO_DH=553`) that is: coarse 1600–2080, grade 2080–2100, fine 2100–2153 km.

## 4. Mesh actually produced (`ISO_NS = 1000`, L = 553 km, coarse Δs = 553.0 m)

| | U | R2 | R4 |
|---|---:|---:|---:|
| total cells | 1000 | **1111** | **1320** |
| fine Δs [m] | — | **276.04** | **138.02** |
| coarse Δs [m] | 553.0 | 552.76 | 552.84 |
| cells: coarse / grade / fine | 1000 / – / – | 868 / 51 / 192 | 868 / 68 / 384 |
| grade ratio `r` | — | **1.0137** | **1.0206** |
| **max adjacent width ratio** | 1.0 | **1.0137** | **1.0206** |
| actual grade span [km] | — | 2079.79–2100 (20.21) | 2079.86–2100 (20.14) |
| cells in the 2130–2150 km window | 37 | 72 | 145 |

All widths positive and finite; inner face exactly at `s=0`, outer face exactly at `s=L`; `s_lo` an exact face. The ratio is bounded well inside 1.1 in both cases, and the fine cells run through the outer boundary (`ds_i(ns−1) == min(ds_i)`).

## 5. Physical validation

Identical to the previous two stages apart from the mesh: domain 1600–2153 km, C7/pchip IC, `data/eos/gamma1_hydrogen_v1.dat`, `ISO_HYDRO_T_DECOUPLE=1`, fixed reservoir back-pressure, 22 000 K conduction wall, cooling off, TRAC off, MC3 β=2, CFL 0.25, `eq_wb` reference, `ISO_NUMERICAL_DIFFUSIVITY_MULT=1` for the primary comparison, `CHROMO_T_END=500`. Diagnostic strides scaled by the refinement factor so all runs capture at the same *physical* cadence.

Metrics are the **existing** definitions, imported (not re-typed) from `resolution_scan_diag.face_metrics` and `numerical_conduction_diag.read_outercond`, over the fixed physical window 2130–2150 km. Full table: `outputs/model_column/outref_report.txt`; JSON: `outputs/model_column/outref_metrics.json`; figure: `visualization/model_column/outref_profiles_500s.png`.

### t = 500 s

| metric | U | R2 | R4 | R4, `mult=0.25` |
|---|---:|---:|---:|---:|
| window cell width [km] | 0.5530 | 0.2761 | 0.1379 | 0.1379 |
| **`A_M`** | **0.13658** | **0.05142** | **0.00943** | **0.03071** |
| **`Q_M`** | **0.08438** | **0.01710** | **0.00308** | **0.00677** |
| `A_mismatch` | 0.13631 | 0.04964 | 0.00951 | 0.02477 |
| `A_diff` | 0.13881 | 0.05210 | 0.01005 | 0.02585 |
| `A_eff` (face flux) | 0.01088 | 0.01182 | 0.00089 | 0.02031 |
| `Q_eff` (face flux) | 1.467e-04 | 9.60e-05 | 6.99e-05 | 1.92e-04 |
| `corr(dcen, ddiff)` | −0.9969 | −0.9759 | −0.9961 | −0.7729 |
| `mean(rho V)` | 2.703e-09 | 2.785e-09 | 3.640e-09 | 1.072e-09 |
| `mean(F_eff)` | 2.694e-09 | 2.781e-09 | 3.640e-09 | 1.072e-09 |
| `T_top` [K] | 20 548.9 | 21 240.7 | 21 552.3 | 21 758.5 |
| `q_phys` [W m⁻²] | 1.598 | 1.738 | 2.088 | 1.140 |
| `q_num` [W m⁻²] | 2.182 | 2.211 | 2.569 | 0.343 |
| **`q_num/q_total`** | **0.577** | **0.560** | **0.552** | **0.231** |
| steps to 500 s | 86 433 | 176 190 | 355 731 | 356 051 |

`f_central + f_diff = f_total` to ≤1.9e-06 in every run, so the capture is self-consistent.

### Reading

1. **The primary success criterion is met with margin:** `A_M(R4) = 0.0094 < A_M(R2) = 0.0514 < A_M(U) = 0.1366`. Ratios 2.66× then 5.45×, i.e. **faster than the first-order rate global refinement gave** (1.10× then 2.03× over the same nominal Δh halvings in the previous scan). Putting the cells only where the gradient is beats spreading them over the whole column, which is the whole point of the profile.
2. **Same story on every ripple metric.** `Q_M` 0.0844 → 0.0171 → 0.0031; `A_mismatch` tracks `A_M` to within 4 % at every mesh (0.136/0.050/0.0095), so `rho V` is converging onto the conserved face flux exactly as the global scan concluded. `A_eff` and `Q_eff` stay small and *improve*, so the conserved flux never roughened.
3. **Better ripple per unit cost than global refinement.** Uniform `ns=2000` reached `A_M = 0.0673` at 500 s with 2000 cells and ~175 k steps (cells×steps = 3.5e8). R2 gets `A_M = 0.0514` — 1.31× smaller — for **0.56× that cost**; R4 gets 0.0094 — **7.1× smaller** — for 1.34× that cost. Global refinement cannot reach `A_M < 0.01` at any resolution in the previous scan's extrapolation (`ns ~ 7000–14 000` was needed merely for 2 %).
4. **No interface artifact.** In the graded band 2079–2101 km, `rho V` is monotone in every refined run (R2: 2.482 → 2.568e-09 across the grade; R4: 3.589 → 3.598e-09) with normalised curvature 1.0–1.1e-04 — the *same* order as the uniform run's in that same physical band (6.6e-05) and ~1000× smaller than the TR ripple it is being compared against. R2's face-flux capture reaches down to 2037 km and shows `max|D2(rho V)|` = 9.4e-13 in the grade versus a window ripple amplitude of ~1.4e-10. No damping of any kind was added at the interface.
5. **Reduced artificial conduction is now viable.** R4 with `ISO_NUMERICAL_DIFFUSIVITY_MULT=0.25` completes the full 500 s (356 051 steps, all state finite, top-100 km velocities 6–47 m/s, no reversal), with `q_num/q_total` down from 0.55 to **0.23** and `A_M = 0.0307` — 4.4× below the uniform production run, though 3.3× above R4 at full diffusivity. So the trade is real and favourable: the refined mesh buys either a much smaller ripple at the same artificial conductivity, or a much more physical drive at a ripple still well below production. (Context: `D_num=0` at uniform `ns=2000` *aborted* at t = 921 s in the previous stage. This run is a 4× reduction rather than removal, so it is not a direct refutation of that, but it is the first configuration in which the coefficient was cut and the ripple stayed below the production level.)
6. **Result D is untouched, and its sign flips.** With `D_num` fixed at the coarse-equivalent value, `mean(F_eff)` *rises* 2.69 → 2.78 → 3.64e-09 as the mesh is refined, because the wall stencil `ds_face` shrinks with the local cell while `numerical_diffusivity` does not: `q_num` goes 2.18 → 2.57 W m⁻² and `q_phys` 1.60 → 2.09 W m⁻². `T_top` marches toward the 22 000 K wall (20 549 → 21 552 K) exactly as in the global scan, so the first-order Dirichlet boundary layer is being resolved — but the evaporation rate still has no grid-independent value. **The evaporation number must not be quoted from any of these runs.**

## 6. Tests

Full suite: **8122 checks passed, 0 failed.**

New:
- `test_mesh_outer_refined_profile` — the 553 km / `ns=1000` outer profile: exact endpoints, exact `s_lo` face, strictly increasing faces, positive finite widths, adjacent ratio ≤ 1.1, outer spacing ≈ coarse/4, *every* cell above `s_lo` fine (so the fine band really reaches the boundary), first cell ≈ coarse, grade foot at least the requested 20 km below `s_lo`, and cell count > `ns_coarse`.
- `test_mesh_outer_profile_edge_cases` — uniform fallback when the grade does not fit below `s_lo` or `s_lo` is outside the domain; and that `τ = 0` still produces a bounded grade rather than an abrupt jump.
- `test_outer_refined_mesh_ic_runs` — end-to-end on the 1600–2153 km domain: `model_column` builds the mesh, `Grid::resize` grows the state, one full semi-implicit step stays finite with `V < 1 cm/s` (the well-balanced hydrostatic fixed point survives refinement, i.e. no boundary-scale blow-up), the finest cell is the last one and the coarsest the first, the adjacent ratio is bounded, and the non-uniform CFL scales with the smallest *local* cell.

Updated: `test_mesh_refined_diagnostic_profile` now checks that the `lower` grade's coarse spacing resumes in [800, 805] km instead of asserting an exact 800 km face (the transition width is now an output — §1). `test_refine_params_env_aliases` gained the `*_PROFILE` parsing/override checks.

Unchanged and still passing: `test_mesh_disabled_is_uniform`, `test_uniform_regression_refinement_unset` (bitwise-reproducible uniform step), `test_irregular_*` (metric caches, manufactured operators, CFL, conduction energy conservation, boundary HSE), `test_refined_mesh_ic_runs` (the `lower` profile end-to-end).

## 7. Cost

Wall-clock is not cleanly separable (runs overlapped), so the honest cost proxy is cells × steps: **1.00 / 2.26 / 5.43** for U / R2 / R4. The step count is set by the acoustic CFL at the smallest cell (86 433 / 176 190 / 355 731 — almost exactly ×2 and ×4), and the cell count grows only 1.11× / 1.32×. Against uniform `ns=2000` at the same 500 s (2000 cells, ~175 k steps ⇒ 3.5e8), R2 costs **0.56×** and R4 **1.34×**.

## 8. Not done / out of scope

No AMR, regridding, moving mesh, conservative remapping, multiple patches, local time stepping, **local diffusivity sensor**, HLL/HLLC, TRAC or BC redesign, cell-averaged IC, or fixed-flux experiment. Only one reduced-diffusivity case (`0.25`, on R4) was run — no sweep. No 1000 s runs. No `lower`-profile physics re-run (its mesh changed by 0.2 km at the transition top; only mesh-property tests cover it).

## 9. Single highest-priority next step

**Make the artificial conductivity local.** Every remaining problem in this stage is the same one: `numerical_diffusivity` is a global scalar tied to the *coarse-equivalent* Δs, so (a) it is not reduced where the mesh is refined, (b) it therefore *increases* the outer conductive drive under refinement (`ds_face` shrinks, `κ_num` does not), and (c) it is why the evaporation rate still does not converge. Replacing `2000·ds_m` with `2000·ds_i(i)` — a one-line change that is exactly grid-consistent and reduces to the current behaviour on a uniform mesh — would let refinement reduce the artificial term automatically and is the direct test of whether `mean(F_eff)` then converges. That is a numerical-scheme change and was deliberately left out of this stage.
