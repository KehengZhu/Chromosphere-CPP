# Upper-BC experiment: decoupling the hydro ghost temperature from the conduction wall

Stage goal: in `model_column`, stop the outer hydro ghost from carrying a *fixed* temperature (22 000 K) and let it zero-gradient-extrapolate the live top cell, while the Stage-D conduction solver keeps seeing exactly the same fixed hot Dirichlet wall. This isolates the extra thermal/entropy constraint the hydro side was imposing at a subsonic outflow face where the fixed reservoir back-pressure `kOuterPRef` already supplies the one admissible incoming condition. No characteristic/NSCBC redesign, no domain extension, no lower-BC change, no pressure-closure change.

## 1. Files changed and the core diff

### `chromosphere.hpp` — two new `Grid` fields (default off)

```cpp
bool  outer_conduction_temperature_override = false;
float outer_conduction_temperature          = 0.0f;
```

When the override is on, the Stage-D **outer Dirichlet datum and the outer-ghost conductivity temperature** come from `outer_conduction_temperature` instead of the decoded hydro ghost. The ghost *density* is still the hydro ghost's, so `kappa_outer` keeps the ghost `n_e`/`n_HI` (re-evaluated by Saha *at the wall temperature*, see below). The override has no effect when `impose_outer_heat_flux` is true — that branch has no Dirichlet datum.

### `src/integrators.cpp` — conduction reads the wall, not the ghost

New helpers next to `decode_gamma_ghost`:

```cpp
static double outer_conduction_wall_T(const Grid& grid, double ghost_T);

struct OuterConductionGhost { double T, n_e, n_HI; };
static OuterConductionGhost outer_conduction_ghost(const Grid& grid,
                                                   const MixtureThermo& ghost) {
    if (!grid.outer_conduction_temperature_override)
        return {ghost.T, ghost.n_e, ghost.n_HI};       // byte-identical baseline
    const double T = static_cast<double>(grid.outer_conduction_temperature);
    const double x = saha_ionization_fraction_n_h(ghost.n_H, T);
    return {T, x*ghost.n_H, (1.0-x)*ghost.n_H};
}
```

Applied at **every** place an active conduction path reads the outer ghost temperature:

| path | site | change |
|---|---|---|
| gamma-mixture solve (`apply_gamma_conduction_stage`) — **the path this run uses** | `k_outer` | `outer.{n_e,n_HI,T}` → `wall.{n_e,n_HI,T}` |
| same | `tr` at `i+1==ns` | `outer.T` → `wall.T` |
| gamma residual diagnostic (`gamma_conduction_residual_max`) | `k_outer`, `tr` | same substitution (keeps the solve and its verifier consistent) |
| fixed-gamma Stage-D (`apply_conduction_stage`), charged row | `T_ghost_out_i[0]` | overwritten with the wall right after `ghost_extended_T`/`ghost_extended_Te` |
| fixed-gamma Stage-D, **neutral row** | `T_ghost_out_n[0]` | same |
| fixed-gamma Stage-D, outer-face conductivity | `T_i_ip1_full[ns-1]`, `T_n_ip1_full[ns-1]` (the ip1-shifted slot *is* the ghost) | pinned to the wall before `kappa_e`/`kappa_n`/TRAC broadening |

Audited and deliberately **not** changed:

* `impose_outer_heat_flux` (Neumann) branches — no Dirichlet datum to override; left as-is, and the experiment keeps `impose_outer_heat_flux = false`.
* `inner_conduction_neumann` / all inner-boundary code — untouched.
* `rhs_implicit_state` (`src/rhs.cpp:712`) contains a second, *explicit* conduction operator that also reads the outer ghost. It is referenced **only** by `tests/chromo_tests.cpp:3337` and by no scenario/integrator, so it is not an active conduction path; left untouched. Flagged here so a future stage that revives it knows to wire the override in.
* Refined/non-uniform mesh: both conduction paths take the same `face_k`/`face_conductivity_series` branch for the outer face regardless of mesh, and the override substitutes upstream of it, so the non-uniform path is covered by the same edit (no separate outer-face code exists).

### `scenarios/model_column.cpp` — the experiment switch

New module flag `kHydroTDecouple` from `ISO_HYDRO_T_DECOUPLE` (default `0`). In the outer-face block the single `T_g0`/`T_g1` is split into three explicitly named quantities:

```cpp
const float T_cond_wall = a_live * kTtopRef;                      // conduction wall
const float T_g0 = kHydroTDecouple ? T_top : T_cond_wall;         // hydro ghost 1
const float T_g1 = kHydroTDecouple ? T_g0  : b_live * T_g0;       // hydro ghost 2
grid.outer_conduction_temperature_override = kHydroTDecouple;
grid.outer_conduction_temperature          = T_cond_wall;
```

Everything else in the outer face is unchanged: `p_g0 = kOuterPRef` (fixed reservoir back-pressure), the one-sided HSE ladder for `p_g1`, `rho_gk = rho_EOS(p_gk, T_gk)` via `gamma_density_from_pressure`, the live ionization split, and `V_g = clip(V_top, ±M_cap c_s)` shared by both ghosts. The `kOuterPRef` capture uses the **hydro** ghost temperature `T_g0`, so the captured ladder is exactly the one the BC reproduces at the reference state in both modes.

One diagnostic line was added at the `kOuterPRef` capture (stderr, once per capture) recording `kOuterPRef`, `p_top`, `T_top`, `T_hydro_g0`, `T_cond_wall` and the flag, so post-run analysis can form the pressure deficit `p_top − p_g0` without re-deriving the anchor.

### `util/upper_bc_decouple_diag.py` (new)

Reads the `.gamma_diag` sidecar and reports all twelve requested metrics. Metric definitions are documented in its docstring; the two non-obvious ones:

* **roughness (top N cells)** = `mean|Δ²(ρV)| / mean|ρV|` over the top N cells.
* **q_top** = `kappa_solver(top) · (T_wall − T_top)/Δs` — a first-order estimate of the discrete Dirichlet flux `kappa_face (T_wall − T_top)/Δs`, using the top cell's own conductivity because the wall-side `kappa` needs the ghost density, which the sidecar does not carry. Identical definition for both runs.

## 2. Boundary information flow, before → after

| quantity | baseline | decoupled |
|---|---|---|
| hydro `p_g0` | `kOuterPRef` (fixed) | `kOuterPRef` (fixed, **unchanged**) |
| hydro `p_g1` | one-sided HSE rung off `p_g0` | same |
| hydro `T_g0` | `a·T_ref` = 22 000 K (external) | `T_top` (zero-gradient, **from interior**) |
| hydro `T_g1` | `b·T_g0` = 22 000 K | `T_g0` = `T_top` |
| hydro `rho_gk` | `rho_EOS(p_gk, 22000)` | `rho_EOS(p_gk, T_top)` |
| hydro `V_g` | `clip(V_top, ±0.1 c_s)` | same |
| conduction outer Dirichlet `T` | 22 000 K (via the hydro ghost) | 22 000 K (via `grid.outer_conduction_temperature`) |
| conduction `kappa_outer` | `kappa(n_e^ghost, n_HI^ghost, 22000)` | `kappa(Saha(rho_ghost, 22000), 22000)` |
| externally imposed hydro conditions | **two** (p and T) | **one** (p) |

## 3. Tests

`tests/chromo_tests.cpp::test_upper_bc_hydro_conduction_temperature_decoupling` (64 checks) on the actual experiment configuration (`ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1`, gamma table, `ns=200`):

1. **Flag off ⇒ unchanged.** Override disarmed; both hydro ghosts at 22 000 K; moving the live top cell 20 % off the wall does **not** move the ghost.
2. **Decoupled hydro T.** `T_hydro_g0 = T_top` and `T_hydro_g1 = T_hydro_g0`.
3. **Conduction still sees the fixed wall, not the live ghost.** Same state and same ghosts: disarming the override (so Stage D falls back to the ghost temperature, which is cooler than the 22 kK wall) makes the step come out *colder*. If Stage D were still reading the ghost the two results would be identical.
4. **Wall sensitivity, ghost invariance.** Raising only `outer_conduction_temperature` raises the post-step top-cell temperature (> 0.1 % change) while `grid.outer_boundary0_i` is bitwise unchanged.
5. **Live-T sensitivity, wall invariance.** Perturbing the top cell moves `T_hydro_g0` to the new `T_top`; `outer_conduction_temperature` stays 22 000 K and remains > 5 % away from the ghost temperature.
6. **Pressure not reverted.** `p_g0` stays exactly the captured back-pressure across a state change that alters `p_top` — no live extrapolation crept back in.
7. **Second-ghost EOS/HSE consistency.** `p_g1 = p_g0 − Δs·½(ρ_g0+ρ_g1)g`, and both `ρ_gk` reproduce an *independent* Saha `ρ(p,T)` inversion written in the test.
8. **Well-balance / no launch.** 20 semi-implicit steps from the reference IC in both modes; the decoupled `max|V|` stays within 10× of baseline (in practice smaller — see §5).

Full suite: **7275 checks passed, 0 failed** (baseline before this stage: 7211 — i.e. +64 new checks, no regressions, no coverage removed). The EOS `rtsafe` performance fix and the one-sided fixed-back-pressure code are untouched.

## 4. Baseline vs decoupled — quantitative comparison

Identical in both runs: domain 1600–2153 km, `ns=1000` (Δs = 0.553 km), C7/pchip IC, `data/eos/gamma1_hydrogen_v1.dat`, cooling off, TRAC off, CFL 0.25 numerics, `kOuterPRef = 0.0102845 Pa` (**identical to printed precision in both modes**), Mach cap 0.1, no relaxation phase, conduction wall 22 000 K, `CHROMO_FRAME_STRIDE=500`.

### At t = 500 s

| metric | baseline | decoupled | change |
|---|---:|---:|---|
| `V_last` [m/s] | **−31.60** | **+91.29** | reversal removed |
| `V_{N−2}` [m/s] | +99.51 | +69.80 | — |
| `V_last / V_{N−2}` | −0.318 | +1.308 | sign flip removed |
| `(ρV)_last/(ρV)_mid` | **−0.327** | **+1.105** | sign flip removed |
| rel-std `ρV`, 2130–2150 km | 0.1174 | 0.1316 | **+12 % (worse)** |
| roughness, top 90 cells (incl. last cell) | 0.0538 | 0.0440 | −18 % |
| roughness, top 90 cells (**excl. last cell**) | 0.0371 | 0.0416 | **+12 % (worse)** |
| max `|ρV − 15-cell smooth|/smooth`, 2100–2148 km | 0.271 | 0.319 | **+18 % (worse)** |
| `V` at the 2146 km ripple peak [m/s] | 46.24 | 45.04 | −2.6 % |
| `T_top` [K] | 21 892 | 20 549 | −1343 |
| `T_hydro_g0 − T_top` [K] | +107.7 | 0 (by construction) | constraint removed |
| `p_top − p_g0` [Pa] | **−5.31e−5** | **+1.33e−5** | \|deficit\| ↓ 4× |
| `q_top` [W/m²] | 0.278 | 3.644 | **13× (not matched)** |
| `E_cond(500 s)` [J/m²] | 123.1 | 1690.3 | **13.7× (not matched)** |
| `M(t)` [kg/m²] | 2.09804e−4 | 2.09771e−4 | — |
| `M(t)/M(0) − 1` | −3.80e−4 | −5.37e−4 | no drainage either way |
| boundary-layer thickness [km] (5 % plateau criterion) | 6.08 | 1.11 | 1-cell feature only |
| max ripple vs plateau | 1.340 | 0.541 | dominated by the last cell |

### Time behaviour (t = 100 / 300 / 400 / 500 s)

| | baseline `V_last` | decoupled `V_last` | baseline rel-std | decoupled rel-std |
|---|---:|---:|---:|---:|
| 100 s | −51.83 | +44.76 | 0.2151 | 0.2197 |
| 300 s | −38.15 | +89.36 | 0.1328 | 0.1438 |
| 400 s | −34.49 | +91.73 | 0.1220 | 0.1354 |
| 500 s | −31.60 | +91.29 | 0.1174 | 0.1316 |

Both runs are settled by ~300–400 s and neither shows secular pressure drift, mass drainage, or velocity runaway. `p_top` is constant to 4 significant figures from 300 s on in both.

### The 1000 s pair (same config, `CHROMO_T_END=1000`)

Because the decoupled run was stable but still slowly relaxing, the pair was extended to 1000 s.

| metric | baseline | decoupled |
|---|---:|---:|
| `V_last` [m/s] | −15.85 | +81.91 |
| `V_{N−2}` [m/s] | +96.88 | +68.98 |
| `V_last / V_{N−2}` | −0.164 | +1.187 |
| `(ρV)_last/(ρV)_mid` | −0.1746 | +1.0847 |
| rel-std `ρV`, 2130–2150 km | 0.1051 | 0.1166 (**+11 %**) |
| roughness top 90 (incl. last cell) | 0.0483 | 0.0364 (−25 %) |
| roughness top 90 (**excl. last cell**) | 0.0323 | 0.0353 (**+9.5 %**) |
| max `\|ρV − smooth\|/smooth`, 2100–2148 km | 0.2012 | 0.2152 (**+7 %**) |
| `T_top` [K] | 21 923 | 20 695 |
| `p_top − p_g0` [Pa] | −6.08e−5 | +8.8e−6 |
| `q_top` [W/m²] | 0.200 | 3.285 |
| `E_cond(1000 s)` [J/m²] | 243.1 | 3421.0 |
| `M/M(0) − 1` | −1.77e−4 | −4.78e−4 |
| boundary layer [km] | 6.08 | 0.55 (1 cell) |

The 1000 s conclusion is the 500 s conclusion, slightly attenuated: the last-cell reversal stays gone; the ripple stays 7–11 % worse. Neither run drifts (mass change < 5e−4 over 1000 s, `p_top` flat), so nothing here is a slow instability.

**Flag-off regression against the archived run.** The new 1000 s baseline reproduces the pre-existing `_onesidedbc` run to 4–5 significant figures — `V_last` −15.8537 vs −15.8531, `(ρV)_last/(ρV)_mid` −0.17461 vs −0.17437, rel-std 0.10512 vs 0.10525 — confirming that with `ISO_HYDRO_T_DECOUPLE=0` the boundary and conduction code are unchanged in practice as well as in principle.

### Historic-baseline cross-check

Re-running the diagnostic on the two archived 1000 s runs reproduces four of the five previously reported numbers exactly — `V_last` −12.37/−15.85 (reported −12.4/−15.9) and rel-std 0.1011/0.1053 (reported 0.101/0.105), `(ρV)_last/(ρV)_mid` −0.1170/−0.1744 (reported −0.117/−0.174). The previously quoted "top-90-cell roughness 0.045 → 0.040" could **not** be reproduced: its definition was never recorded in the repo, and every natural definition I tried (mean/rms of first or second differences, residual from a linear fit, normalised by mean or max) makes the one-sided run marginally *rougher* than the centered one, consistent with the direction of every other metric. The definition now used is documented in the script docstring and applied identically to all runs.

## 5. Well-balance verification

* **Reference-state consistency.** At the IC, `T_top,ref = 22275.3 K`. In decoupled mode `T_hydro_g0 = T_hydro_g1 = 22275.3 K` exactly — the ghost is the true zero-gradient continuation, so `BC(q_ref) = q_g,ref`. (In *baseline* mode the ghost is 22 000 K, i.e. 275 K *colder* than the reference top cell — the baseline ghost was never the natural HSE continuation.)
* **`kOuterPRef` recaptured under the new semantics**, from the hydro ghost temperature: `0.0102845 Pa` in both modes (the ~1.2 % ghost-density difference acts on a term that is 0.04 % of `p_top`, so the anchor is identical to printed precision). Back-pressure parity between the two runs therefore holds.
* **Clean runs only.** Both runs start from a fresh IC; `eq_residual` is computed lazily on the first step *after* `model_column_update_bc`, i.e. with the new boundary in place. The flag is never flipped mid-run (it is read once in `model_column_ic`), so a stale `eq_residual` cannot be reused across a semantics change.
* **Undriven HSE regression** (`ISO_HEAT_FLUX=0`, conduction off, 50 s, `ns=1000`):

  | | max\|V\| [m/s] | `V_last` | `M/M(0)−1` |
  |---|---:|---:|---:|
  | baseline | 0.0408 | +0.0065 | +2.45e−8 |
  | decoupled | **0.0306** | −0.0285 | +3.08e−8 |

  The decoupled boundary is *slightly better balanced*, with no acoustic launch and mass conserved to 3e−8.

## 6. Conductive-input parity — **not achieved**

The *imposed* thermal boundary is identical by construction and by test: 22 000 K Dirichlet, `impose_outer_heat_flux = false` in both runs, same `kappa` model, same TRAC state (off), same numerical diffusivity.

The *realized* conductive flux is **not** matched, because the interior responded: freed from the fixed hydro ghost, the top cell settles at 20 549 K instead of 21 892 K, so the wall-to-top-cell temperature difference grows from 108 K to 1451 K and `q_top` grows 13× (`E_cond(500 s)`: 123 → 1690 J/m²). This is a genuine non-equivalence and it must temper the attribution below: the decoupled run is not simply "the same heating with one fewer constraint".

Mechanistically this is itself informative. With outflow at the top face the *upwind* state is the interior, so a hot ghost should not advect enthalpy inward; the fact that removing the 22 kK ghost cools the top cell by 1343 K shows the fixed hot ghost was doing real thermodynamic work on the last cell through the reconstruction and the Rusanov dissipation, and through its (lower) EOS density. That is exactly the mechanism an over-specified thermal condition would produce — but it also means the two runs differ in more than the constraint count.

## 7. Does this support the over-specification hypothesis?

**Partially — for the last-cell reversal only, not for the boundary-layer ripple.**

* **Supported:** the fixed hydro temperature is the cause of the *single-cell* pathology. The last-cell velocity reversal is completely eliminated (−31.6 → +91.3 m/s), the mass-flux sign flip is gone (`(ρV)_last/(ρV)_mid` −0.327 → +1.105, i.e. the top cell now continues the interior flux instead of opposing it), and the top-cell pressure deficit against the reservoir drops 4× and changes sign (−5.31e−5 → +1.33e−5 Pa). No drainage, no drift, no instability, and slightly better well-balancing.
* **Not supported:** the ~20–30 km mass-flux ripple over 2130–2150 km is **unchanged in shape and marginally worse in amplitude**. Every ripple metric that excludes the last cell moves the wrong way by 12–18 % (rel-std 0.117 → 0.132; roughness-excl-last 0.0371 → 0.0416; max deviation from a 15-cell smooth 0.271 → 0.319), and the ripple peaks sit at the same heights (~2135, ~2146, ~2150 km) with nearly the same velocities (46.2 → 45.0 m/s at 2146 km). See `visualization/model_column/upper_bc_hydroT_decouple_500s_top90.png`: the two `ρV` curves are the same curve apart from the final point.

So the earlier "boundary layer" metric was conflating two distinct phenomena. The fixed hydro `T` explains the last physical cell; it does not explain the standing ripple in the 20–30 km below it.

## 8. Remaining risk and suggested next steps

1. **The unmatched conductive input is the main caveat.** A cleaner second round would hold the *realized* flux fixed rather than the wall temperature — e.g. re-run the decoupled case with the wall lowered so that `q_top` matches the baseline's 0.278 W/m², and check whether the last-cell reversal stays away. That separates "removed a constraint" from "changed the heating by 13×".
2. **The ripple needs a different explanation.** Ranked candidates, none tested here: (a) ghost-side MUSCL/MC3 reconstruction across the outer face (the ripple is a ~3-cell-wavelength standing pattern, the signature of a limiter/reconstruction artefact, and it sits exactly where the TR temperature gradient steepens); (b) the split conduction–hydro operator against a strong `T` gradient (the known `ΔS`-proportional artefact recorded for the refined-mesh sweep); (c) the rigidity of the *fixed* `p_back` itself, which the decoupled run leaves fully intact; (d) TR truncation height — the top of the domain sits inside the steep gradient; (e) genuine need for characteristic relaxation.
3. **Float32 storage** at the outer face limits the pressure-deficit metric to ~1e−7 Pa resolution; the deficits reported here (1e−5 Pa) are two decades above that, so they are meaningful, but a tighter experiment would want a double-precision ghost diagnostic.
4. **Not done, deliberately:** no characteristic/NSCBC boundary redesign, no pressure-relaxation or sponge, no domain extension, no lower-BC change, no change to `impose_outer_heat_flux`, and no change to the existing default behaviour (`ISO_HYDRO_T_DECOUPLE` defaults to 0 and the full suite is unchanged at the flag-off setting).

## 9. Commands and artefacts

Runs (project root as working directory):

```bash
# baseline (flag off)
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=1000 ISO_HYDRO_T_DECOUPLE=0 \
    CHROMO_T_END=500 CHROMO_FRAME_STRIDE=500 \
    ./build/chromo_main \
    outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_500s_hydroTbase.txt \
    full no-ionization model_column - 30.0 no-cooling \
    > outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_500s_hydroTbase.console.log 2>&1

# decoupled (flag on) — same command with ISO_HYDRO_T_DECOUPLE=1 and
# ..._500s_hydroTdecoupled.txt / .console.log

# 1000 s pair: identical, CHROMO_T_END=1000, ..._1000s_hydroTbase / _1000s_hydroTdecoupled

# undriven HSE well-balance regression (conduction off)
#   ISO_HEAT_FLUX=0 CHROMO_T_END=50 CHROMO_FRAME_STRIDE=2000
#   -> outputs/model_column/model_column_h1600_2153_ns1000_condoff_50s_hydroT{0,1}.txt
```

Diagnostics and figure: see `visualize_commands.md`
(`upper_bc_hydroT_decouple_500s_top90.png`, `upper_bc_hydroT_decouple_500s_metrics.json`,
and the two evolution movies).
