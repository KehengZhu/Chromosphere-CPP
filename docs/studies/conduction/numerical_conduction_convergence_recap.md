# Is `numerical_diffusivity = 2000 Δh` why the evaporation flux does not converge?

> Historical global-scaling study. The production implementation now stores
> `C_num = 2000 m/s` (before scenario/diagnostic multipliers) and evaluates
> `chi_num,f = C_num*ds_face` locally. See
> `docs/studies/conduction/local_mesh_scaled_numerical_conduction_recap.md` for the implementation and
> targeted validation.

Targeted follow-up to Result D of `docs/studies/validation/resolution_convergence_scan_recap.md`. Only the boxed question was tested; nothing was extended.

## Headline

**Partially — it is the single largest term, but it is not the whole story.**

Numerical conduction supplies **57.7 %** of the outer conductive input at `ns=1000` and **38.5 %** at `ns=2000`, and switching it off collapses the outer flux ~7× and the mean evaporation flux ~6×. So most of what the production runs call "evaporation drive" is the artificial conductivity. But with `D_num=0` the evaporation flux is **still not grid-converged** — it moves **+34 % / +16 %** from `ns=1000` to `ns=2000` (500 s / ~900 s), the *opposite sign* to production's −33 % / −25 % and of comparable magnitude. And `D_num=0` at `ns=2000` **aborts at t = 921.5 s**.

## 1. The production outer conductive flux, exactly

From `src/integrators.cpp::apply_gamma_conduction_stage` (top row `i = ns-1`), with `outer_conduction_ghost` supplying the wall state and `scenarios/model_column.cpp:544` the coefficient:

```
q_total = ( κ_phys,face + κ_num,face ) · (T_wall − T_top) / ds_face     [W m⁻²]

κ_phys,face = ½ ( κ(top) + κ(wall) )                      (uniform_mesh ⇒ arithmetic mean;
                                                           non-uniform ⇒ distance-weighted harmonic)
κ(·)        = κ_e(n_e, n_HI, T)·TRAC(T) + κ_n(n_e, n_HI, T)
κ_num,face  = numerical_diffusivity · C_V^eff(ρ_top, T_top)
numerical_diffusivity = (kCorona ? 4 : 1) · 2.0e3 · Δh     ⇒ 1.106e6 at ns=1000, 5.53e5 at ns=2000
ds_face     = ds_iph_i(ns-1) = Δh                          (top-cell CENTRE → ghost CENTRE, NOT Δh/2)
```

Sign convention: **positive `q_total` heats the top cell**. It enters the discrete energy equation as the flux-tube divergence `g_right·(T_wall − T_top) = (B_i/B_iph) · q_total / ds_i`, i.e. `area_ratio = B_i(ns-1)/B_iph(ns-1)` (= 1 for this straight column).

Two points that a cell-centred coefficient cannot represent, and which the `.gamma_diag` sidecar therefore got wrong:

- `κ(wall)` uses the **ghost density with Saha re-evaluated at `T_wall`**, not the top cell's ionization. At the wall's 22 000 K that is a materially different `κ_e ∝ T^{5/2}`.
- `κ_num,face` at the top uses the **cell** `C_V` (no face averaging, since there is no `i+1` cell), while every interior face averages two cells.

`impose_outer_heat_flux` (`ISO_QFLUX`) would replace this Dirichlet term entirely with `grid.outer_heat_flux/ds_i`; it is **off** in every run here, and the capture records that flag.

## 2. Diagnostic added

A minimal, default-off, write-only capture — `Grid::capture_outer_conduction` → `OuterConductionCapture` (`chromosphere.hpp`), filled at the end of `apply_gamma_conduction_stage` from the coefficients of the **final converged Newton iteration** (`diag_kr_top` / `diag_cv_r_top` are rewritten every iteration, so after the accepted solve they are the accepted values). No solver stage reads it back. `chromo_main` writes one row per captured step to `<out>.outercond` under `CHROMO_OUTER_COND_DIAG=1` (`CHROMO_OUTER_COND_STRIDE`, default = the face-flux stride).

Columns: `t step T_top T_wall kappa_phys_face kappa_num_face ds_face area_ratio q_phys q_num q_total`. No per-cell output, no energy-budget infrastructure, no generic conduction framework.

Override: `ISO_NUMERICAL_DIFFUSIVITY_MULT` multiplies **only** `grid.numerical_diffusivity`. Unset ⇒ production. `=0` ⇒ `κ_num ≡ 0` with physical conductivity, hydro, boundaries, EOS, limiter, CFL and well-balancing untouched.

## 3. The four runs

Domain 1600–2153 km, decoupled hydro-T upper BC, 22 000 K conduction wall, fixed reservoir back-pressure, cooling off, TRAC off, MC3 β=2, CFL 0.25, `data/eos/gamma1_hydrogen_v1.dat`, C7/pchip IC, `CHROMO_T_END=1000`. `ns=500` not run.

`ISO_NUMERICAL_DIFFUSIVITY_MULT=1` is verified equivalent to unset: the production runs reproduce the previous stage's numbers to every printed digit — `mean(F_eff)` = `2.6978e-09` / `2.3999e-09` (ns=1000, 500 s / 1000 s) and `1.8222e-09` / `1.8417e-09` (ns=2000), and `T_top(1000 s)` = 20 695.4 / 21 361.8 K.

## 4. Results

`q` [W m⁻²] at the outer face, into the top cell. `mean(F_eff)` over 2130–2150 km [kg m⁻² s⁻¹]. Full table: `outputs/model_column/ncond_report.txt`; JSON: `outputs/model_column/ncond_metrics.json`.

### t = 500 s (all four runs matched)

| run | Δh [km] | `D_num` | `T_top` [K] | `T_wall−T_top` [K] | `q_phys` | `q_num` | `q_total` | `q_num/q_total` | `mean(F_eff)` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| prod, ns=1000 | 0.5530 | 1.106e6 | 20 548.2 | 1451.8 | 1.5987 | 2.1827 | 3.7814 | **0.577** | 2.6978e-09 |
| prod, ns=2000 | 0.2765 | 5.53e5 | 21 371.5 | 628.5 | 1.4496 | 0.9074 | 2.3570 | **0.385** | 1.8222e-09 |
| `D_num=0`, ns=1000 | 0.5530 | 0 | 21 550.0 | 450.0 | 0.5242 | 0 | 0.5242 | 0 | 4.6400e-10 |
| `D_num=0`, ns=2000 | 0.2765 | 0 | 21 721.5 | 278.5 | 0.6553 | 0 | 0.6553 | 0 | 6.2312e-10 |

### t ≈ 900 s (the latest time the aborted run reaches)

| run | `T_top` [K] | `q_phys` | `q_num` | `q_total` | `q_num/q_total` | `mean(F_eff)` |
|---|---:|---:|---:|---:|---:|---:|
| prod, ns=1000 | 20 667.0 | 1.4776 | 1.9918 | 3.4694 | **0.574** | 2.4584e-09 |
| prod, ns=2000 | 21 360.4 | 1.4743 | 0.9240 | 2.3983 | **0.385** | 1.8361e-09 |
| `D_num=0`, ns=1000 | 21 590.6 | 0.4780 | 0 | 0.4780 | 0 | 3.7411e-10 |
| `D_num=0`, ns=2000 | 21 745.2 | 0.6003 | 0 | 0.6003 | 0 | 4.3378e-10 |

`D_num=0`, ns=2000 is quoted from `t_F = 897.3 s` for `mean(F_eff)` but `t_q = 762.5 s` for `q` — its `.outercond` write buffer was lost in the abort. `T_top` and `q_phys` are drifting by <1 % per 100 s there, so this does not affect any conclusion, but it is not a matched time and is flagged in the report and on the figure.

### `ns = 1000 → 2000`

| branch | time | `q_total` | `q_phys` | `mean(F_eff)` |
|---|---|---|---|---|
| production | 500 s | 3.7814 → 2.3570 (**−37.7 %**) | 1.5987 → 1.4496 (−9.3 %) | 2.6978e-09 → 1.8222e-09 (**−32.5 %**) |
| production | 900 s | 3.4694 → 2.3983 (**−30.9 %**) | 1.4776 → 1.4743 (**−0.2 %**) | 2.4584e-09 → 1.8361e-09 (**−25.3 %**) |
| `D_num=0` | 500 s | 0.5242 → 0.6553 (**+25.0 %**) | same | 4.6400e-10 → 6.2312e-10 (**+34.3 %**) |
| `D_num=0` | 900 s | 0.4780 → 0.6003 (**+25.6 %**) | same | 3.7411e-10 → 4.3378e-10 (**+15.9 %**) |

Figure: `visualization/model_column/ncond_numerical_conduction.png`.

## 5. Reading

1. **Numerical conduction dominates the production drive, and its share is exactly what scales.** `q_num` falls 1.99 → 0.92 W m⁻² (−54 %) between `ns=1000` and `ns=2000` at 900 s, while `q_phys` is **flat to 0.2 %** (1.4776 → 1.4743). The entire −31 % change in the production `q_total` is `q_num`. On the drive side, the hypothesis is fully confirmed.
2. **Absolutely, the artificial term is not a correction.** At `ns=1000` it *exceeds* the physical flux (2.18 vs 1.60 W m⁻² at 500 s). Removing it drops `q_total` 3.47 → 0.48 W m⁻² and `mean(F_eff)` 2.46e-9 → 3.74e-10 at `ns=1000` — a factor 7 and 6.6. **The production evaporation rate is mostly artificial-conductivity-driven.** That is the more consequential finding than the convergence question.
3. **But `D_num=0` does not converge.** `mean(F_eff)` moves +34 % (500 s) / +16 % (900 s) and `q_phys` +25 % across the same refinement. Removing the `O(Δh)` term did not leave a grid-independent remainder; it left a comparable one of the opposite sign — consistent with the earlier scan's two *unseparated* confounds cancelling against `κ_num` in the production runs: the `O(Δh)`-dependent IC top-cell temperature (21 574 / 22 275 / 22 668 K) and the `Δh` in the Dirichlet stencil at the wall. Note `T_wall − T_top` is 450 → 278 K in the `D_num=0` branch, i.e. **not** halving with `Δh`, so `(T_wall−T_top)/Δh` *rises* by 25 % — the wall gradient itself is what moves now.
4. **Note the sign flip in `q_phys`.** Production `q_phys` is flat with resolution, `D_num=0` `q_phys` rises 25 %. Same physical conductivity formula; the difference is the state it is evaluated on. With `κ_num` present the top-cell temperature is dragged much further from the wall (`T_top` = 20 667 vs 21 591 K at `ns=1000`), which is what makes production `q_phys` look converged. That coincidence is not evidence of a converged physical boundary layer.

## 6. `D_num = 0` stability

`ns=1000` completed cleanly (176 095 steps, t = 1000 s). **`ns=2000` aborted at t = 921.5 s (step ≈ 326 700 of 350 540 needed) with an uncaught `std::domain_error: rho_n must be positive and finite`.** Reported as-is; no limiter, CFL or boundary was touched to rescue it.

What the state looked like on the way in, from the surviving sidecars:

- The **top 300 cells stayed benign to the last capture (t = 919.7 s)**: `T_max` 21 755 K (rising 21 641 → 21 755 over 800 s, flattening), `ρ_min` 2.860e-11 kg m⁻³ (flat to 3 digits), `|V|_max` 37 m/s and *falling*. No shock, no drainage, no thermal runaway, no CFL collapse.
- The one precursor that distinguishes it from the production run: in the last full snapshots `MOM_I` and `MOM_N` go **negative near 2140 km** (−1.37e-10 / −3.26e-11 at t = 874.8 s), first appearing between t ≈ 706 and 847 s. The production `ns=2000` run has strictly positive momentum everywhere at t = 1000 s. So a small localized flow reversal at ~2140 km develops **only without numerical conduction** — which is precisely the ringing the coefficient was added to damp.
- The row that actually threw is the **trace neutral carrier**: `RHO_N` is pinned at `eos_trace_fraction_floor = 1e-8 × ρ` throughout the top of the domain, i.e. ≈2.9e-19 kg m⁻³ in float32. Both runs sit there (`RHO_N/RHO_I = 1.000e-08` exactly at the top cell in both). It is a float32 trace-row failure triggered by the local reversal, not a hydrodynamic blow-up of the physical state.

Consequence for this stage: the `D_num=0` non-convergence claim rests on 500 s and ~900 s, both of which the run reaches. No 1000 s `D_num=0` `ns=2000` number is quoted.

## 7. Verdict

**Only partially explains the result.**

- Confirmed: numerical conduction is 57.7 % / 38.5 % of the outer conductive input at `ns=1000` / `ns=2000`; it is the *entire* resolution-dependence of the production `q_total` (`q_phys` flat to 0.2 % at 900 s); and its presence accounts for 85 % (ns=1000) / 76 % (ns=2000) of the production evaporation mass flux in absolute terms.
- Not confirmed: that it is the *reason* the evaporation flux fails to converge. With it removed, `mean(F_eff)` still moves +34 % / +16 % across the same refinement, so a second, comparable, opposite-signed mechanism remains.
- New, unasked-for but load-bearing: the production evaporation rate is dominated by an artificial coefficient, so its *magnitude* — not only its convergence — is not physically defensible as currently configured.

## 8. Files

Modified (4):
- `chromosphere.hpp` — `OuterConductionCapture` struct; `Grid::capture_outer_conduction` + `outer_conduction_capture`
- `src/integrators.cpp` — `apply_gamma_conduction_stage` mirrors the top face's `kr`/`cv_r` and fills the capture after convergence. No coefficient, stencil or BC changed.
- `scenarios/model_column.cpp` — `diff_mult` gains the `ISO_NUMERICAL_DIFFUSIVITY_MULT` factor (default 1.0 = production)
- `chromo_main.cpp` — `CHROMO_OUTER_COND_DIAG` / `CHROMO_OUTER_COND_STRIDE`, `.outercond` sidecar

New:
- `util/numerical_conduction_diag.py` — reads `.outercond` and reuses `resolution_scan_diag.face_metrics` / `face_flux_diag.read_faceflux` for `mean(F_eff)`, so no metric is re-typed. Includes a tolerant face-flux reader for the sidecar the aborted run truncated.
- `tests/chromo_tests.cpp::test_outer_conduction_capture_matches_solver_face`
- `docs/studies/conduction/numerical_conduction_convergence_recap.md` — this file

Runs (`outputs/model_column/`, each with `.faceflux`, `.gamma_diag`, `.outercond`, `.console.log`):
- `ncond_ns1000_prod_1000s.txt`, `ncond_ns2000_prod_1000s.txt`
- `ncond_ns1000_dnum0_1000s.txt`, `ncond_ns2000_dnum0_1000s.txt` *(ns=2000 ends at 921.5 s)*

Results: `outputs/model_column/ncond_report.txt`, `outputs/model_column/ncond_metrics.json`, `visualization/model_column/ncond_numerical_conduction.png`. Commands: the new section of `visualize_commands.md`.

Nothing was changed in the Rusanov flux, the limiter, the EOS, the well-balanced subtraction, any boundary condition, the conduction stencil, the CFL, or the integration order. No production default changed — both new knobs are default-off/identity.

## 9. Tests and diagnostic neutrality

- **Full suite: 7705 passed, 0 failed** (7704 previous + the one new test).
- `test_outer_conduction_capture_matches_solver_face` checks: bit-for-bit identical advanced state with the capture off vs on; `T_wall == grid.outer_conduction_temperature`; `ds_face == ds_iph_i(ns-1) == ds_i(ns-1)` (the full cell width, not half); `area_ratio == 1`; `κ_phys,face == ½(κ(top) + κ(wall))` with `κ(wall)` rebuilt independently from the ghost density and Saha at `T_wall`; `κ_num,face == numerical_diffusivity · C_V(ρ_top,T_top)`; `q_phys + q_num == q_total`; `q_num > 0.2 q_total` at this resolution; and that `ISO_NUMERICAL_DIFFUSIVITY_MULT=0` zeroes `q_num` exactly while leaving `κ_phys,face` within a few tenths of a percent.
- **Run-level neutrality**, ns=1000, 20 s, `CHROMO_OUTER_COND_DIAG=0` vs `=1`, everything else identical: `.gamma_diag` **byte-identical**; state `.txt` identical in all data lines (only the embedded output filename in the header comment differs).

## 10. Single highest-priority next step

Rebuild the IC by **cell averaging** rather than centre point-sampling of the C7 pchip profile, and repeat the `D_num=0` `ns=1000`/`2000` pair. That removes the one remaining `O(Δh)` confound that is cheap to remove and is the most likely source of the residual +25 % `q_phys` / +16–34 % `mean(F_eff)` drift, and it does not require any scheme change.
