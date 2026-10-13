# MC3 branch fix and the hydrostatic reconstruction

**Status: CURRENT** — the MC3 branch fix and the hydrostatic reconstruction are both in the release. The hydrostatic reconstruction was a prototype behind `ISO_HSE_RECON=1` (commit `7a0e1af`) and was promoted to the release, with the switch removed, on 2026-10-12 (see *Promotion to release*). Labels A/B/C below are the configurations as they were measured. Supersedes the hydrostatic-residual numbers of `reference_free_release_recap.md` (measured before the fix).

## Question

The advisor pointed out that limiting `ln p` directly ignores gravity. In a column at rest `ln p` falls by `Δs/H` per cell, with `H = k_B T (1+x)/(m_H g)`, so the cell-to-cell differences the limiter compares are mostly the hydrostatic fall. The remedy is to limit only the departure from discrete hydrostatic balance and to discretize gravity consistently with it. This study asks three things. Does the release limiter actually act on the column? Is the MC3 implementation right? Does a hydrostatic reconstruction help?

## What changed

1. **MC3 branch fix (release).** `limit_mc3` in `src/single_fluid/mixture.cpp` uses the code's ratio `r = Δ₋/Δ₊` with the limited slope applied to `Δ₊`. BATSRUS `ModFaceValue` 'mc3' gives the upper face `u_i + (Δ₋ + 2Δ₊)/6` and the lower face `u_i − (2Δ₋ + Δ₊)/6`, so the `+` branch must be `(r+2)/3` and the `−` branch `(2r+1)/3`. The code had them swapped, which is the `κ = −1/3` scheme rather than the documented third-order `κ = 1/3`. On cell averages of `s²` at `s = 1, 2, 3` the old code gives the face value at 2.5 as 5.917, BATSRUS and the fixed code 6.25 (exact). The frozen legacy helpers `flux_lim_mc3_plus/minus` keep the old form.
2. **Hydrostatic reconstruction (`Grid::hydrostatic_reconstruct`; the prototype `ISO_HSE_RECON=1` when measured, the release since 2026-10-12).**
   - Each cell carries its half-cell step `o_i = ρ_i g_par ds_i/(2p_i)`, from the decoded state at `U^n`. Ghosts take the boundary cell's `g_par` and width.
   - The `ln p` slot of every stencil neighbour is moved into cell `i`'s hydrostatic frame by removing the steps `h_{j+1/2} = o_j + o_{j+1}` between them. MC3 then limits only the departure, and face values are shifted back by `±o_i`. The predictor half-step centres use the same frame.
   - The gravity source becomes `2 p_i sinh(o_i)/ds_i`, which is exactly the pressure difference of cell `i`'s own extrapolations and equals `ρ_i g_par` to `O((ds/H)²)`.
   - The IC and both ghost ladders are built in the same discrete form, `ln p_{i+1} − ln p_i = o_i + o_{i+1}`.
   - `ln ρ` and `u` are unchanged. The flux-tube term is unchanged, so exact balance holds for `B ≡ 1` (the current model) only.
   - `model_column_2t` and the `ISO_RECONSTRUCTION=lnrho-v-lnt` reference run without it and keep the trapezoidal IC and ghost forms.

   One implementation detail mattered: `phi_g` carries a float-literal km→m factor, so the gravity the solver differences differs from `grid.g` by 4.7e-8 relative. Balancing the IC against `grid.g` left a uniform 1.3e-5 m/s² sinking, 2.6e-4 m/s after 20 s at every resolution. The IC and ghosts therefore use the solver's own `(φ_iph − φ_imh)/ds` (`solver_g` in `scenarios/model_column.cpp`). The release path is not changed by this; the 5e-8 mismatch is far below its truncation error.

## Evidence

Runs under `outputs/model_column/hse_recon/`. Labels: **A** = pre-change binary (commit `f01189b`), **B** = MC3 fixed, **C** = B + `ISO_HSE_RECON=1`. Release preset (N=500/R4, 661 cells; N=1000 is 1320 cells), CFL 0.50, `godunov.fallbacks=0` and `termination=end_time` (or `step_cap`) on every run.

### Does the limiter act?

On the release IC and at 20, 40 and 60 s of a conduction-driven run, the `ln p` slope ratio stays within 0.95–1.13 (|r−1| ≤ 3e-3 below 2000 km, ≤ 4.8e-2 in the lower TR). MC3 with β = 2 clips the upper face for `r < 0.4` or `r > 4` and the lower face for `r < 0.25` or `r > 2.5`, so it **never clips** here. Over the whole 0–1000 s run (66 761 face values) r stays within 0.64–1.19. `ds/H` is about 4e-3 in the chromosphere and 1e-3 in the TR. The hydrostatic residual therefore comes from the unlimited reconstruction plus the cell-centred `ρg` source, not from limiting. Its maximum is in the upper TR (2149 km), where `H` changes fastest. Ratios were estimated from the 6-digit `.gamma_diag` output; that noise is ~1e-3 in `r`, far from the clip thresholds.

### Column at rest (20 s, `ISO_HEAT_FLUX=0`, `ISO_T_TOP=0` so the outer ghost matches the IC top cell)

`util/release/hydrostatic_residual_metrics.py`. A reproduces `reference_free_release_recap.md` digit for digit.

| | A | B | C |
|---|---|---|---|
| N=500 t=0 max\|F^ρ\| [kg m⁻² s⁻¹] | 1.7782e-13 @ 2149.1 km | 6.0148e-14 @ base face | 1.5458e-19 |
| N=500 base face f_diff / f_total | 1.294e-13 / 1.246e-13 | 6.467e-14 / 6.015e-14 | 3.7e-21 / 4.7e-21 |
| N=500 20 s max\|V\| [m/s] | 5.7274e-03 (top cell) | 7.4004e-03 (top cell) | 4.3009e-10 |
| — excluding the top cell | 4.97e-3 @ 2149 km | 4.70e-3 @ 1783 km | — |
| N=1000 t=0 max\|F^ρ\| | 3.9112e-14 | 1.5540e-14 | 2.7066e-19 |
| N=1000 20 s max\|V\| | 1.5291e-03 | 1.9337e-03 | 2.7898e-08 |
| order (F, V) | 2.19, 1.90 | 1.95, 1.94 | round-off |

The base-face dissipative flux is exactly halved by the fix (1.294 → 0.647). That is the `κ = 1/3` prediction: on smooth data the left/right face jump is `(1−κ)/4` times the third difference, 1/6 instead of 1/3.

5-step ladder (CFL 0.25, `CHROMO_STEP_CAP=5`, `ISO_T_TOP=0`), max\|V\| [m/s] and t=0 max\|F^ρ\| [kg m⁻² s⁻¹]:

| `ISO_NS` (cells) | A max\|V\| | B max\|V\| | C max\|V\| | A F | B F | C F |
|---|---|---|---|---|---|---|
| 48 (70) | 8.109e-1 | 1.660 | 4.3e-10 | 3.32e-11 | 1.43e-11 | 4.0e-20 |
| 125 (169) | 7.617e-2 | 1.790e-1 | 1.3e-10 | 4.39e-12 | 1.70e-12 | 1.2e-19 |
| 250 (331) | 1.878e-2 | 4.156e-2 | 2.3e-9 | 8.38e-13 | 3.49e-13 | 1.7e-19 |
| 500 (661) | 4.986e-3 | 1.076e-2 | 1.2e-9 | 2.14e-13 | 8.89e-14 | 1.5e-19 |
| 1000 (1320) | 1.317e-3 | 2.863e-3 | 3.4e-9 | 5.08e-14 | 1.96e-14 | 2.7e-19 |
| 2000 (2638) | 3.321e-4 | 7.438e-4 | 3.0e-9 | 1.31e-14 | 5.04e-15 | 3.5e-19 |

B orders 2.33, 2.12, 1.95, 1.91, 1.94; the same ladder with `ISO_RIEMANN=roe-local` gives 1.634, 1.779e-1, 4.145e-2, 1.074e-2, 2.861e-3, 7.435e-4 m/s (orders 2.32–1.91), so second order is still a property of the discretization, not of the flux. At N=500 the 5-step maximum of A and B sits on the top cells (A 4.99e-3, B 1.08e-2 at 2153 km); below 2100 km B is slightly lower (2.9e-5 vs 3.5e-5 m/s). So the fix improves the base face and the bulk, but roughly doubles the short-time residual in the steepest top cells. C is round-off at every rung.

The release regression `test_stage8_gamma_model_column_saha_hse_and_ghosts` (5 steps, CFL 0.50, decoupled outer T) moves from 4.9e-3 / 1.4e-13 / 7.9e-9 to 1.04e-2 / 2.9e-13 / 3.3e-9 (max\|V\|, max\|ρu\|, ‖ΔU‖/‖U‖); its velocity bound was raised from 2e-2 to 4e-2 to keep 4x headroom. The new `test_hydrostatic_reconstruction_exact_fixed_point` gives 9.7e-10 m/s and 7.6e-15.

### Conduction-driven, 1000 s (release configuration)

`util/release/pressure_reconstruction_diag.py`, metrics in `outputs/model_column/hse_recon/prod_metrics.json`:

| t = 1000 s | A500 | B500 | C500 | A1000 | B1000 | C1000 |
|---|---|---|---|---|---|---|
| mean V [m/s] | 5.766 | 5.787 | 5.767 | 5.817 | 5.839 | 5.818 |
| max V [m/s] | 12.387 | 12.409 | 12.392 | 12.492 | 12.512 | 12.494 |
| F_eff [kg m⁻² s⁻¹] | 3.4301e-10 | 3.4453e-10 | 3.4320e-10 | 3.4422e-10 | 3.4575e-10 | 3.4435e-10 |
| q_phys [W m⁻²] | 0.4210 | 0.4208 | 0.4210 | 0.4237 | 0.4235 | 0.4238 |
| column mass drift | −8.816e-5 | −1.196e-4 | −7.438e-5 | −8.828e-5 | −1.198e-4 | −7.426e-5 |
| base face f_diff/f_total | −0.161 | −0.174 | −0.001 | −0.162 | −0.175 | −0.001 |

- **Evaporation is insensitive to both changes.** B moves mean V and F_eff by +0.36–0.44 %, about 40 % of the N500 → N1000 difference (+0.9 %) and in the same direction. C differs from A by ≤ 0.06 %, so the hydrostatic reconstruction undoes the B shift. At 100 s all three agree to 0.03 %.
- **Column mass.** B loses 36 % more mass by 1000 s and C 16 % less. The differences are the same at both resolutions and are set up in the first ~100 s (−3.39e-4, −3.72e-4, −3.21e-4 at 100 s). They are a start-up-transient boundary effect, not a resolution effect.
- **First-interior-face artifact.** The dissipative share of the base-face mass flux, listed in `release_facts.md` as open (~16 %), is 16–17 % in A and B and **0.1 % in C**, at both resolutions. The artifact is the hydrostatic inconsistency of the reconstruction and the lower ghosts.

## Interpretation and status

- The advisor's diagnosis is right in kind: the release does not reconstruct relative to hydrostatic balance, and the residual it leaves is a second-order truncation error largest where `H` varies fastest. The limiter clipping is not the mechanism on this mesh; the unlimited reconstruction and source mismatch are.
- The MC3 fix makes the code do what its documentation claims. Its measured effect is mixed and small: better at the base face and in the bulk, worse in the top cells for short times, and +36 % column mass loss. Keeping it in the release is a correctness choice, not a measured improvement.
- The hydrostatic reconstruction makes the column at rest exact to round-off, removes the first-interior-face artifact, and leaves evaporation unchanged. It was a candidate for the release, with three open items; their resolution is under *Promotion to release*:
  - the variable-`B` flux-tube term in the balanced form;
  - a decision on `ln ρ`, which is still reconstructed raw (harmless for balance at `u = 0`, since the exact Riemann solver returns zero mass flux at equal pressures);
  - a long-run (4000 s) and CFL check.

## Promotion to release (2026-10-12)

The hydrostatic reconstruction is now always on in `model_column` (`Grid::hydrostatic_reconstruct`). Setting `ISO_HSE_RECON` at all throws, so a stale `ISO_HSE_RECON=0` cannot silently run the release. The `ISO_RECONSTRUCTION=lnrho-v-lnt` reference and `model_column_2t` run without it and keep the trapezoidal IC and ghost forms. The three open items above were resolved as follows:

- **Variable `B`.** Not needed for the current model, which has `B ≡ 1`. Recorded as a scope limit in `release_facts.md`: a variable-`B` model would bring back a residual.
- **`ln ρ`.** Left as is. At rest a density jump is a stationary contact, and the exact Riemann solver gives it zero mass flux. That would not hold for a Rusanov-type flux, which is only the per-face fallback and recorded zero fallbacks on every run here.
- **Long run and CFL.** Measured below.

**Promotion is exact.** A 20 s release run of the new default build is byte-identical (`.txt`, `.gamma_diag`) to the prototype binary run with `ISO_HSE_RECON=1`. All ctest tests pass, including the legacy build, and so does `scripts/release_validation.sh` (default and explicit paths bit-identical; 100 s smoke with `godunov.fallbacks=0`).

The column-at-rest regression moved into `test_hydrostatic_reconstruction_exact_fixed_point` (`max|V| < 1e-8 m/s`). Three face-flux tests normalized the continuity check by the face flux itself; at rest that flux is now round-off. They now impose a linear 0–20 m/s upflow on their test column.

### 4000 s, release configuration (B = plain `ln p` with the MC3 fix, C = release)

N=500/R4, CFL 0.50, single-T outer closure. Runs `outputs/model_column/hse_recon/long_{B,C}_N500_4000s`. Every run ends `termination=end_time` with `godunov.fallbacks=0`. Face metrics come from the `.faceflux` capture nearest each time. The spread is `(max − min)/mean` of `f_total` over faces below 1850 km. The mean V, `F_eff` and drift columns come from `pressure_reconstruction_diag.py`.

| t | | base/top `f_total` | spread 1600–1850 km | base-face `\|f_diff/f_total\|` | cell-0 `ρV/f_total[0]` | top V [m/s] | top `f_total` | mean V [m/s] | `F_eff` | mass drift |
|---|---|---|---|---|---|---|---|---|---|---|
| 500 s | B | 1.152 | 20.8 % | 17.4 % | 1.369 | 19.56 | 5.588e-10 | 7.473 | 4.7125e-10 | −1.748e-4 |
| | C | 1.039 | 11.4 % | 0.10 % | 1.000 | 19.50 | 5.578e-10 | 7.455 | 4.7179e-10 | −1.110e-4 |
| 4000 s | B | 1.0036 | 13.0 % | 17.5 % | 1.369 | 9.418 | 2.682e-10 | 4.689 | 2.4690e-10 | −7.98e-5 |
| | C | 1.0022 | 13.0 % | 0.11 % | 1.000 | 9.412 | 2.683e-10 | 4.689 | 2.4694e-10 | −5.07e-5 |

- **Evaporation is unchanged at 4000 s.** Mean V is identical, and `F_eff`, top V and top flux agree to ≤ 0.06 %. The old release (pre-MC3 fix) gave 9.42 m/s and 2.683e-10 at the top (`release_facts.md`), so the published long-run evaporation numbers stand.
- **The first-interior-face artifact stays removed for the whole run:** 0.10–0.11 % at every time, and the base cell's `ρV` matches the face flux (1.000 against 1.369).
- **The ~13 % lower-chromosphere spread is unchanged** (12.96 % against 13.04 % at 4000 s), as expected for the physical adiabatic-contraction signature (`docs/studies/boundaries/twall_mass_flux_relaxation_timescales.md`). Early on the two differ: at 500 s the spread is 11.4 % in C against 20.8 % in B, and base/top is 1.039 against 1.152. That difference was not decomposed further.
- **Column mass loss is 36 % smaller** by 4000 s, consistent with the 1000 s result.

### CFL 0.25 control at 500 s

| 500 s, CFL 0.50 → 0.25 | mean V | `F_eff` | top V | top `f_total` |
|---|---|---|---|---|
| B | 7.473 → 7.472 | +0.47 % | +0.32 % | +0.29 % |
| C | 7.455 → 7.455 | +0.09 % | +0.16 % | +0.16 % |

The release is no more CFL-sensitive than the scheme it replaces, and CFL 0.50 stays validated for this configuration.

## Reproduce

```bash
# A = binary built at f01189b; B = binary built at 7a0e1af; C = that binary with ISO_HSE_RECON=1,
# which is now the default (current builds reject ISO_HSE_RECON)
# column at rest
ISO_NS=500 ISO_T_TOP=0 ISO_HEAT_FLUX=0 CHROMO_T_END=20 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 \
CHROMO_FRAME_DT=5 CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=0 CHROMO_FACE_FLUX_STRIDE=40 \
  scripts/run_chromo_realtime.sh outputs/model_column/hse_recon/static_tt0_B_N500_20s.txt \
  full no-ionization model_column - 20.0 no-cooling
util/release/hydrostatic_residual_metrics.py outputs/model_column/hse_recon/static_tt0_B_N500_20s.txt
# 5-step ladder: add CHROMO_CFL=0.25 CHROMO_STEP_CAP=5 CHROMO_FRAME_STRIDE=5 CHROMO_FACE_FLUX_STRIDE=5, vary ISO_NS
# conduction-driven
ISO_NS=500 CHROMO_T_END=1000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=10 \
CHROMO_FACE_FLUX_DIAG=1 CHROMO_OUTER_COND_DIAG=1 CHROMO_FACE_FLUX_STRIDE=1780 \
CHROMO_OUTER_COND_STRIDE=89 CHROMO_EOS_COUNTS=1 \
  scripts/run_chromo_realtime.sh outputs/model_column/hse_recon/prod_B_N500_1000s.txt \
  full no-ionization model_column - 20.0 no-cooling
util/release/pressure_reconstruction_diag.py --run B=outputs/model_column/hse_recon/prod_B_N500_1000s.txt --targets 100,1000
```
