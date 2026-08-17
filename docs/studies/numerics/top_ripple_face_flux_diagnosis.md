# Top-region ripple: finite-volume face mass-flux diagnosis

> **Partly superseded (see `docs/studies/numerics/reference_free_release_recap.md`).** The release hydrodynamic operator has since gained the MUSCL-Hancock predictor momentum source and a trapezoidal EOS-closed lower ghost ladder, and **equilibrium-reference well balancing (`eq_wb` / `ISO_EQ_WB`) has been retired from the release** — it survives only in the two-fluid research solver. Statements below that describe `eq_wb` as active release method, or that quantify the discrete hydrostatic defect, are historical. Everything else in this document still stands.


Read-only follow-up to `docs/studies/boundaries/upper_bc_hydro_temperature_decoupling_recap.md`. That stage established that the fixed hydro ghost temperature caused the *last-cell* reversal but not the ~20–30 km ripple below the top. This stage asks the sharper question: **does the ripple exist in the finite-volume face mass flux that continuity is actually differenced from, or only in the cell-centred product `rho*V`?**

Answer: **it exists only in `rho*V`. Case A.** The conserved face flux is smooth to ~1e-4 in normalised curvature while `rho*V` ripples at ~7e-2 — a factor of ~830. There is no continuity violation, no boundary instability, and the limiter and the operator splitting are both cleared by direct experiment.

## 1. Active code paths audited

Everything below was read, not assumed. `model_column` + `GAMMA_TABLE` runs the **gamma-mixture** path exclusively; the fixed-gamma path is dead code for these runs.

| concern | resolved in | finding |
|---|---|---|
| which RHS runs | [src/rhs.cpp:409](src/rhs.cpp#L409) `rhs_explicit_state` | non-empty `eos_gamma_table` ⇒ `rhs_explicit_mixture`; the fixed-gamma body below it is never entered |
| reconstruction variables | [src/rhs.cpp:118-123](src/rhs.cpp#L118-L123) `decode_mixture_field_into::put` | **`(ln rho_total, V_cm, ln T)`** — 3 slots, the *total* mixture, `V_cm = (MOM_I+MOM_N)/rho_total`. Not the 7 primitive rows, and there is no separate `log_reconstruct` step: log is baked into the stencil variables. Slots 3–6 of the packed vector stay zero, so their limiter ratios go non-finite and are zeroed. |
| stencil / ghosts | [src/rhs.cpp:190](src/rhs.cpp#L190) `mixture_stencil_view`, [src/rhs.cpp:159-174](src/rhs.cpp#L159-L174) | one `extended_primitive` array of `ns+4` rows; ghost rows 0,1 and `ns+2`, `ns+3` are the decoded `inner_boundary1/0_i`, `outer_boundary0/1_i`, decoded at the **face** potential `phi_g_iph(ns-1)` |
| limiter | [src/rhs.cpp:223](src/rhs.cpp#L223) `mixture_mc3` | `phi_+ = max(0, min(beta r, beta, (2r+1)/3))`, `phi_- = max(0, min(beta r, beta, (r+2)/3))`, `beta = 2`; non-finite `r` ⇒ 0. Asymmetric, so the L and R faces get different limiters. |
| final face states | [src/rhs.cpp:368-378](src/rhs.cpp#L368-L378) | corrector: `wl_iph = 0.5(w + wt) + 0.5 phi_+(r) (w_ip1 - w)` and `wr_iph = 0.5(w_ip1 + wt_ip1) - 0.5 phi_-(r_ip1) (w_ip2 - w_ip1)`, with `wt` the conservatively-updated predictor re-inverted through the EOS ([src/rhs.cpp:359-363](src/rhs.cpp#L359-L363)). The predictor is *not* a cons2prim of the old state — it is an independent caloric inversion. |
| sound speed / spectral radius | [src/rhs.cpp:247-255](src/rhs.cpp#L247-L255), [src/eos.cpp:462-486](src/eos.cpp#L462-L486) | each one-sided face state is pushed through `equilibrium_mixture_face_state`, i.e. a **full Saha + Gamma1 table solve at the reconstructed `(rho, T)`**; `c_s = sqrt(Gamma1 (p_i+p_n)/rho)`; `a = |V| + c_s` written into **all seven** rows, then `a_iph = max(a_L, a_R)` |
| mass flux rows | [src/eos.cpp:489](src/eos.cpp#L489) `equilibrium_mixture_flux`, [src/state.cpp:157](src/state.cpp#L157) | continuity is **split over two rows**, `RHO_I` and `RHO_N`. The equilibrium projection re-partitions `x_row` every step, so only `RHO_I + RHO_N` is a conserved density. All fluxes below are that sum. |
| geometry / sources | [src/rhs.cpp:264-282](src/rhs.cpp#L264-L282) `mixture_source`, [scenarios/model_column.cpp:208](scenarios/model_column.cpp#L208) | `B_i = B_imh = B_iph = 1`; the source touches only `MOM_I`/`MOM_N`. **Continuity has no source term.** Confirmed. |
| well-balanced subtraction | [src/rhs.cpp:397](src/rhs.cpp#L397) | `rhs -= grid.eq_residual` is applied to the **whole packed vector, continuity rows included**. This is the load-bearing fact of the whole diagnosis — see §3. |
| projection conservation | [src/eos.cpp:350-390](src/eos.cpp#L350-L390) | `rho = rho_i + rho_n` and `momentum = momentum_i + momentum_n` are carried through exactly (velocity is common), energy via a remainder. Total mass and total momentum are preserved; the row split is not. |
| telescoping | verified analytically **and** asserted in the new test | `flux_imh[i] == flux_iph[i-1]` exactly: `r_ip1[i-1] == r[i]` and the predictor terms line up, so the interior scheme is conservative face-by-face. |
| sidecar sufficiency | [chromo_main.cpp:257-276](chromo_main.cpp#L257-L276) | the existing `.gamma_diag` is **not** sufficient to rebuild the production flux offline: it lacks the four ghost rows, `eq_residual`, and the predictor half-step, and reproducing the latter needs the tabulated Saha inversion inside the RHS. Hence the in-code capture below rather than a Python reimplementation. |
| other conduction/RHS paths | [src/rhs.cpp:712](src/rhs.cpp#L712) `rhs_implicit_state` | a second, explicit conduction operator referenced only by tests, never by a scenario. Untouched. |

## 2. The diagnostic hook (default OFF)

`Grid::capture_face_flux` (default `false`) makes one `rhs_explicit_mixture` evaluation copy its **final** face arrays into `Grid::face_flux_capture` ([chromosphere.hpp](../../../chromosphere.hpp), `GammaFaceFluxCapture`; filled by `capture_gamma_face_flux` in [src/rhs.cpp](../../../src/two_fluid/rhs.cpp)). It is pure copy-out: no solver stage reads it back. Captured per face `i+1/2`:

- reconstruction inputs `rho_cell, v_cell, T_cell` (the same `w` the limiter saw);
- one-sided reconstructed states `rho_L/R, v_L/R, T_L/R`, and `c_L/R` recovered exactly from the `|V| + c_s` spectral-radius rows;
- `a_face`, and the total-mass split `f_central = 0.5[(rho v)_L + (rho v)_R]`, `f_diff = -0.5 a (rho_R - rho_L)`, `f_total` (read straight out of the production `flux_iph`, not recomputed);
- `eq_residual_mass`, the frozen well-balanced residual's continuity row;
- the limiter ratios and values actually used: `r, phi_+` (L state) and `r_ip1, phi_-` (R state) for each of the three slots.

`chromo_main` writes them to `<out>.faceflux` under `CHROMO_FACE_FLUX_DIAG=1` (`CHROMO_FACE_FLUX_TOP`, `CHROMO_FACE_FLUX_STRIDE`). The record is emitted *after* the step but holds the RHS of the state at the **start** of that step; since every column comes from that one call, a record is self-consistent and needs no cross-matching. The final step of a run is always captured.

**Face flux obtained is the production flux, not a look-alike.** New test `test_face_flux_capture_matches_production_continuity` checks, on the real ns=200 gamma column:

- `rhs` with the flag on is **bit-for-bit** equal to `rhs` with it off (`approx_equal(..., "absdiff", 0.0)`);
- `f_central + f_diff == f_total` to 1e-5 (float32 rounding only);
- `a_face == max(|v_L| + c_L, |v_R| + c_R)` to 1e-12;
- `-(f_total[i] - f_total[i-1])/ds - eq_residual_mass[i]` reproduces the **production `RHO_I + RHO_N` rows** of `rhs_explicit_state` to 1e-5 of `|F|/ds`, for every interior cell.

The last item is the proof: the captured numbers *are* what continuity is differenced from, including the telescoping and the eq_wb subtraction, with nothing else in the row.

## 3. The correct conserved flux under `eq_wb`

`model_column` runs with `grid.eq_wb = true` and the subtraction hits the continuity rows, so the discrete continuity update is

```
d(rho_i)/dt = -( F[i] - F[i-1] )/ds  +  ( F_ref[i] - F_ref[i-1] )/ds
```

with `F_ref` the frozen reference flux. **The discrete steady state therefore requires `F - F_ref = const`, not `F = const`.** Reporting `std(F_total)` alone is meaningless here, and it is exactly what makes `F_total` look violently structured (relstd 0.42) when nothing is wrong.

`F_ref` is recovered from the t = 0 capture, which is taken at the reference equilibrium. Verified: its divergence equals the captured `eq_residual` continuity row to **1.3e-6** relative (`eq_ref_max_relerr` in every report). Define `f_eff = f_total - f_ref`.

Why `F_ref` is large: the IC carries the imposed 22 kK top, hence a steep density gradient, hence a large `-0.5 a (rho_R - rho_L)`. At the ripple peaks `f_ref` reaches `+5.9e-9` and even changes sign (`-1.17e-9`) — comparable to or larger than the physical flux `~2.9e-9`.

## 4. Quantitative decomposition (2130–2150 km)

Full tables (this window, top 90 cells, and top 90 excluding the last cell) in `outputs/model_column/faceflux_{500,1000}s_report.txt` and `..._metrics.json`.

| metric | 500 s base | 500 s dec | 1000 s base | 1000 s dec |
|---|---:|---:|---:|---:|
| `mean(rho V)` | 2.9026e-09 | 2.7031e-09 | 2.6736e-09 | 2.4088e-09 |
| `relstd(rho V)` | 0.12222 | 0.13658 | 0.11355 | 0.12626 |
| **`roughness(rho V)`** | **0.07375** | **0.08438** | **0.06040** | **0.06821** |
| `mean(f_central)` | 2.9667e-09 | 2.7695e-09 | 2.7377e-09 | 2.4747e-09 |
| `mean(f_diff)` | 1.7573e-11 | 1.6820e-11 | 1.7747e-11 | 1.7280e-11 |
| `relstd(f_total)` | 0.41667 | 0.44640 | 0.45132 | 0.49900 |
| `roughness(f_total)` | 0.31332 | 0.33536 | 0.33889 | 0.37257 |
| **`mean(f_eff)`** | **2.8936e-09** | **2.6956e-09** | **2.6647e-09** | **2.4013e-09** |
| `relstd(f_eff)` | 0.00909 | 0.01129 | 0.01007 | 0.01181 |
| **`roughness(f_eff)`** | **8.853e-05** | **1.481e-04** | **1.580e-04** | **1.431e-04** |
| `diff_fraction` mean | 0.51678 | 2.33852 | 1.13917 | 0.41347 |
| `diff_fraction` median | 0.04283 | 0.03746 | 0.04908 | 0.05250 |
| `diff_fraction` max | 13.601 | 80.461 | 35.503 | 6.891 |
| `R_total` absmax | 1.4176e-11 | 1.4176e-11 | 1.4176e-11 | 1.4177e-11 |
| `R_total` rms | 3.1391e-12 | 3.1390e-12 | 3.1390e-12 | 3.1389e-12 |
| `R_total` rms / (\|F\|/ds) | 0.58169 | 0.62267 | 0.62922 | 0.69170 |
| **`R_eff` absmax** | **7.887e-15** | **1.006e-14** | **9.448e-15** | **8.628e-15** |
| **`R_eff` rms** | **4.802e-15** | **5.569e-15** | **4.848e-15** | **5.040e-15** |
| **`R_eff` rms / (\|F\|/ds)** | **8.899e-04** | **1.100e-03** | **9.719e-04** | **1.110e-03** |
| `relstd(d_central)/mean(f_eff)` | 0.12497 | 0.13831 | 0.11646 | 0.12894 |
| `relstd(d_diff)/mean(f_eff)` | 0.12459 | 0.13766 | 0.11591 | 0.12815 |
| **`corr(d_central, d_diff)`** | **-0.99735** | **-0.99666** | **-0.99625** | **-0.99580** |

`roughness(x) = mean|D2 x| / mean|x|`; `d_central = f_central - f_central,ref`, `d_diff = f_diff - f_diff,ref`.

Reading:

1. `mean(f_eff)` and `mean(rho V)` agree to **0.3 %** (2.8936e-09 vs 2.9026e-09). The mean mass flux is right; only its spatial fluctuation differs.
2. `f_eff` is **smooth**: normalised curvature 8.9e-05 against 7.4e-02 for `rho V` — a factor **833**. Its 0.9 % relstd is a *monotone* 8 % rise across 50 km (bottom panel of the figure), i.e. the genuine slow mass loading of a still-relaxing column, not an oscillation.
3. The compensating pair is explicit: `d_central` and `d_diff` each fluctuate by ~0.125 of `mean(f_eff)` — the same size as the `rho V` ripple, 0.122 — and are anti-correlated at **-0.997**.
4. `R_eff` rms is 4.8e-15 kg m^-3 s^-1, i.e. 8.9e-04 of `|F|/ds`. Over 500 s that integrates to ~2.4e-12 kg m^-3, ~3 % of the local `rho` — the residual slow relaxation, not a broken constraint.
5. `R_total` looks 600x worse and is the quantity to *ignore*: it is dominated by the reference divergence that `eq_wb` removes, and it is essentially identical (1.4176e-11) in all four runs and at both times, i.e. a frozen property of the IC, not of the solution.

### Ripple peaks

Per-peak dumps (`rho_L/R, V_L/R, T_L/R, c_L/R, a, f_central, f_diff, f_total, f_ref, f_eff`, limiter `r`, `phi`, and the active MC3 branch) are printed in the reports. Representative baseline 500 s peaks:

| cell km | `f_central` | `f_diff` | `f_total` | `f_ref` | **`f_eff`** | `rho V` cell |
|---|---:|---:|---:|---:|---:|---:|
| 2149.41 | +3.876e-09 | +4.996e-09 | +8.872e-09 | +5.932e-09 | **+2.939e-09** | +3.905e-09 |
| 2145.54 | +3.661e-09 | +8.057e-10 | +4.466e-09 | +1.550e-09 | **+2.916e-09** | +3.672e-09 |
| 2143.88 | +3.066e-09 | -1.420e-10 | +2.924e-09 | +1.866e-11 | **+2.906e-09** | +3.202e-09 |
| 2144.43 | +2.834e-09 | -1.092e-09 | +1.743e-09 | -1.166e-09 | **+2.909e-09** | +2.839e-09 |

`f_total` spans 1.74e-09 to 8.87e-09 (5x) and `rho V` spans 2.84e-09 to 3.91e-09; `f_eff` is **2.906e-09 to 2.939e-09, i.e. flat to 1.1 %.**

## 5. Secondary: limiter / reconstruction — **cleared**

Same domain, resolution, IC, EOS, boundary and time; only the reconstruction changed, via the new diagnostic-only `ISO_LIMITER` override (default unset = the untouched MC3 beta=2 configuration). `first` is `beta = 0` with MC3, i.e. `phi == 0`, i.e. piecewise-constant.

| metric (2130–2150 km, 500 s) | MC3 beta=2 | minmod | first order |
|---|---:|---:|---:|
| `relstd(rho V)` | 0.12222 | 0.15001 | **0.38826** |
| `roughness(rho V)` | 0.07375 | 0.07363 | 0.06587 |
| `mean(f_diff)` | 1.7573e-11 | 9.7340e-10 | **1.4984e-08** |
| `relstd(f_total)` | 0.41667 | 0.42372 | 0.56349 |
| **`relstd(f_eff)`** | **0.00909** | **0.01042** | **0.01121** |
| **`roughness(f_eff)`** | **8.853e-05** | **1.055e-04** | **1.204e-04** |
| `R_eff` rms / (\|F\|/ds) | 8.899e-04 | 7.401e-04 | 1.355e-04 |
| `corr(d_central, d_diff)` | -0.99735 | -0.99782 | **-0.99982** |

Cell-centred `rho V` extrema (km), top region:

```
MC3     2127.8 2130.1 2131.7 2133.9 2135.6 2138.3 2141.1 2141.7 2143.9 2144.4 2145.5 2147.7 2149.4 2151.1 2152.2
minmod         2132.3 2133.4 2135.0 2136.1 2137.2 2140.0        2143.3 2145.0 2145.5 2147.2 2149.4 2151.6 2152.2
dec     2127.8 2130.1 2131.7 2133.9 2135.6 2138.3 2141.1 2141.7 2143.9 2144.4 2145.5 2147.7 2149.4        2152.2
first                                2136.1 2138.9 2141.1               2145.0 2146.1 2147.7 2149.4
```

Conclusions, against the pre-registered rules:

- **First order does not remove the ripple — it triples it** (relstd 0.122 -> 0.388) while `f_eff` stays equally flat. So reconstruction/limiting is *not* the cause; going first-order raises `mean(f_diff)` by 850x (to 5x the physical flux) and the advective part must cancel proportionally more.
- **MC3 and minmod give the same wavelength and the same peak positions** (2135.6, 2141.1, 2143.9/2143.3, 2145.5, 2147.7/2147.2, 2149.4 all shared) and near-identical roughness (0.0738 vs 0.0736). No limiter-branch-switching signature.
- At every peak both faces sit on the **`third_order` branch** of MC3 (`r ~ 0.7–1.2`), for `rho` and `T` alike — no adjacent-cell switching between `beta r`, the `beta` cap, and zero. The only branch flip is in the **velocity** slot, where `r_v` is frequently negative (`-0.58`, `-2.7`, `-37`) and `phi_+ = 0`, i.e. the V reconstruction is first-order at exactly these cells. That is a *consequence* of the `rho V` fluctuation, and it is not enough to structure the conserved flux.
- All three schemes have flat `f_eff` and only differ in cell-centred `rho V` -> by the pre-registered rule, "mainly the decomposition between reconstruction and numerical diffusion".

## 6. Secondary: operator splitting — **cleared**

CFL 0.25 / 0.125 / 0.0625, everything else fixed, compared at the same physical time (500 s):

| metric (2130–2150 km) | CFL 0.25 | CFL 0.125 | CFL 0.0625 |
|---|---:|---:|---:|
| `relstd(rho V)` | 0.12222 | 0.12280 | 0.12678 |
| `roughness(rho V)` | 0.07375 | 0.07517 | 0.07901 |
| `relstd(f_eff)` | 0.00909 | 0.00921 | 0.00610 |
| `R_eff` rms / (\|F\|/ds) | 8.899e-04 | 8.828e-04 | 8.608e-04 |
| `relstd(d_central)/mean(f_eff)` | 0.12497 | 0.12472 | 0.12949 |
| `corr(d_central, d_diff)` | -0.99735 | -0.99727 | -0.99889 |

The ripple amplitude does **not** decrease with `dt` — it is flat to slightly increasing across a 4x reduction, and `R_eff` normalised is flat to 3 %. Neither `O(dt)` nor `O(dt^2)`. **Hydro–conduction operator splitting is not a source of the ripple.** Caveat recorded: `mean(rho V)` at 500 s differs by ~20 % between CFL 0.25 and 0.0625 (2.90e-09 vs 3.50e-09), i.e. the *slow relaxation transient* is `dt`-sensitive even though the relative ripple is not.

Because the CFL scaling is flat and §4 already accounts for the ripple, the per-substage (hydro / projection / conduction) intra-step dump described in the task was **not implemented**. Recorded as deliberately not run, not as a null result.

## 7. Classification

**Case A.** `f_eff` — the finite-volume face mass flux the continuity update actually sees — is constant to 1 % in amplitude and smooth to 1e-4 in normalised curvature; `rho V` ripples at 7e-2; the central and diffusive halves are anti-correlated at -0.997 with individually matching amplitudes. Per the pre-registered rule this must **not** be called a continuity violation or a boundary instability.

Whether the compensation is *too large* to leave the solution physically interpretable is a separate, legitimate concern: `|f_diff|/|f_total|` has a median of only ~4 % but reaches 13.6 (baseline) and 80.5 (decoupled) at individual faces, and the top-cell `rho V` in the baseline swings to -0.9e-09 (a reversal) while the conserved flux there is smoothly positive. **Cell-centred `rho V` is not a trustworthy diagnostic of mass transport in the top ~20 km of this configuration; `f_eff` is.**

### Confirmed
1. The gamma path reconstructs `(ln rho_total, V_cm, ln T)` with asymmetric MC3 beta=2, a caloric-inverted predictor, and Rusanov `a = max(|V|+c_s)` per face; continuity is `RHO_I + RHO_N`, has no source, and `B == 1`.
2. The captured face flux is the production flux (bit-identical RHS; continuity rows reproduced to 1e-5).
3. Under `eq_wb` the conserved flux is `f_total - f_ref`; `f_ref` verified against `eq_residual` to 1.3e-6.
4. `f_eff` is flat and smooth; `mean(f_eff)` = `mean(rho V)` to 0.3 %; `R_eff` = 8.9e-04 of `|F|/ds`.
5. The `rho V` ripple is the advective half of a -0.997-anti-correlated pair.
6. Reconstruction/limiting is not the cause (first order makes it 3x worse; MC3 and minmod share peak positions).
7. Operator splitting is not the cause (no `dt` scaling over 4x).
8. `R_total` (1.4176e-11) is a frozen property of the IC, identical across all runs and times.

### Most likely explanation
The steep IC density gradient at the imposed 22 kK top makes the Rusanov dissipative mass flux `-0.5 a (rho_R - rho_L)` comparable to, and at some faces several times larger than, the physical flux (`a ~ 1.2–2.0e4` m/s against `V ~ 30–90` m/s, so `a/V ~ 200–600`). `eq_wb` removes its reference value exactly, leaving a smooth conserved flux; but the *live* dissipative flux still varies from face to face, and the cell-centred momentum adjusts locally to cancel it. The visible ripple is that adjustment. Its amplitude is set by the ratio of the local dissipative flux variation to `rho`, not by any boundary or scheme defect.

### Not verified
- **What fixes the ripple positions.** `corr(rho V, f_ref)` is only +0.65 (window) / +0.49 (top 90), so the reference-flux imprint is a partial, not sufficient, explanation. The alternative — that the positions are set by the live `T`/ionisation structure through `a` and `Gamma1` — was not tested.
- **Whether `f_eff` is converged.** `R_eff` is small but nonzero (~3 % of `rho` per 500 s), and `mean(f_eff)` still drifts (2.894e-09 at 500 s -> 2.665e-09 at 1000 s). These are relaxing states, not steady ones.
- **Grid convergence.** No resolution scan was run; the whole study is at `ns = 1000`, `ds = 0.553` km.
- Whether the `-0.997` cancellation degrades the *energy* and *momentum* rows the same way (only the mass row was decomposed).
- Whether a less dissipative Riemann solver (HLL/HLLC) would shrink the compensation. Untested and out of scope.

### Explicitly not done
No characteristic/NSCBC boundary, no pressure relaxation, no sponge, no domain extension, no lower-boundary change, no conduction-wall or fixed-flux change, no EOS redesign, no floor tuning, no change to any production default. `ISO_LIMITER` and the CFL sweep are diagnostic legs, unset/default in production. The EOS `rtsafe` fix and the one-sided fixed-back-pressure closure are untouched.

## 8. Suggested next-step ordering (not executed)

1. **Resolution scan** (`ns` = 500 / 1000 / 2000). If `roughness(rho V)` falls while `roughness(f_eff)` stays flat, the ripple is a fixed-`ds` dissipation artefact and needs no further explanation. Cheapest decisive test.
2. **Riemann-solver dissipation.** The mechanism is `a >> |V|`. An HLL/HLLC or a `beta`-limited Rusanov would directly test it. This is a scheme change, so it needs approval.
3. **Reference-flux dependence.** Re-relax the IC so `F_ref` is small, then re-measure. Separates "imprint of the frozen reference" from "live dissipation variation".
4. **Energy/momentum row decomposition** with the same capture, to see whether the same cancellation degrades the heat-flux and velocity diagnostics.
5. Only then revisit the boundary closure — the boundary is now cleared for the ripple twice over (hydro `T`, and now the conserved flux is smooth through the whole region including the outer face).

## 9. Tests and flag-off regression

- `test_face_flux_capture_matches_production_continuity` (new, see §2). Full suite: **7682 passed, 0 failed** (7275 before; the increase is this test's per-cell assertions over ns=200).
- **Flag-off / diagnostic-neutrality regression.** The four new runs were made with `CHROMO_FACE_FLUX_DIAG=1` and compared cell-by-cell against the archived runs from the previous stage over all ten `.gamma_diag` columns and all 1000 cells:

| run | max relative difference |
|---|---|
| 1000 s baseline vs `..._1000s_hydroTbase` | **0.000e+00** |
| 1000 s decoupled vs `..._1000s_hydroTdecoupled` | **0.000e+00** |
| 500 s baseline vs `..._500s_hydroTbase` | **0.000e+00** |

Bit-for-bit identical, including `V_last` (-15.85370 / +81.90920 / -31.60110). Mass, top velocity, pressure and temperature are all inside those columns.

## 10. Files

New:
- [util/face_flux_diag.py](../../../util/face_flux_diag.py) — the analysis script (metric definitions in its docstring)
- `docs/studies/numerics/top_ripple_face_flux_diagnosis.md` — this document
- `outputs/model_column/faceflux_{500,1000}s_hydroT{base,decoupled}.txt{,.faceflux,.gamma_diag}` + `.console.log`
- `outputs/model_column/faceflux_500s_{lim_minmod,lim_first,cfl0125,cfl00625}.txt{,.faceflux,.gamma_diag}` + `.console.log`
- `outputs/model_column/faceflux_{500,1000}s_report.txt`, `faceflux_{500,1000}s_metrics.json`, `faceflux_500s_{limiter,cfl}_metrics.json`
- `visualization/model_column/faceflux_{500,1000}s_decomposition.png`, `faceflux_500s_limiter_compare.png`, `faceflux_500s_cfl_compare.png`

Modified: [chromosphere.hpp](../../../chromosphere.hpp) (`GammaFaceFluxCapture`, `capture_face_flux`), [src/rhs.cpp](../../../src/two_fluid/rhs.cpp) (`capture_gamma_face_flux` + one guarded call), [chromo_main.cpp](../../../chromo_main.cpp) (`CHROMO_FACE_FLUX_*` sidecar), [scenarios/model_column.cpp](../../../scenarios/model_column.cpp) (`ISO_LIMITER` diagnostic override), [tests/chromo_tests.cpp](../../../tests/chromo_tests.cpp), [visualize_commands.md](../../../visualize_commands.md).

Exact commands: see the "Top-region face mass-flux decomposition" section of [visualize_commands.md](../../../visualize_commands.md).
