# Visualization catalog — `model_column` (release)

Figures and movies of the release solver: the single-fluid Saha/`Gamma1` column with the SWMF Godunov flux and `double` storage. Outputs land in `visualization/model_column/`. Conventions and the index of the other catalog files are in [README.md](README.md). Since 2026-10-12 the top level holds only outputs of the current release, which includes the hydrostatic reconstruction. Everything made from the earlier code, which covers every entry below unless it says otherwise, is in `visualization/model_column/2026-10-12_7a0e1af/` (see its `INDEX.md`).

---

## Release mesh and face-flux diagnostics

### Release N=500 mesh and initial temperature

`model_column_N500_mesh_temperature_initial.png` and `.pdf` -- exact cell width and
EOS-aware initial (`t=0`) temperature versus height for the **release**
coarse-equivalent N=500 column. The static x4 outer refinement above 500 km gives
661 cells over 1600--2153 km: coarse 1104.768 m grading to fine 276.042 m, with the
fine region starting at 2100 km. Initial temperature spans 6631.35--22674.80 K.

Provenance is the current release smoke run `release_godunov_lnp_N500_100s`, written by
`scripts/release_validation.sh` with no environment overrides (console log:
`flux=swmf-godunov reconstruction=lnrho-v-lnp hydrostatic=1 limiter=mc3 beta=2`). The
figure was regenerated from it on 2026-10-13. The earlier version, drawn from the Roe-era run
`release_roe_lnp_N500_100s`, is in `2026-10-12_7a0e1af/`. The script independently
rebuilds the faces from `scenarios/mesh.cpp` and cross-checks them against the run
header. It accumulates them in double, as the release does, and the maximum face error
is 5.000 m: the six-significant-digit header rounding.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
  util/release/plot_model_column_mesh_temperature.py \
  --run outputs/model_column/release_godunov_lnp_N500_100s.txt \
  --ns-coarse 500 \
  --output visualization/model_column/model_column_N500_mesh_temperature_initial
```

### `model_column/*_mass_flux.mp4` — cell-centred `rho V` vs the face mass flux

Read-only comparison movie for the release finite-volume mass transport:
cell-centred `rho V` against the production face mass flux and its
equilibrium-reference decomposition. Everything is read from the
`<run>.faceflux` sidecar (`CHROMO_FACE_FLUX_DIAG=1`, release solver only), so all
four series come from the SAME `mixture_rhs_explicit` evaluation at the same step
and the panels are synchronous by construction.

Definitions (index `i` = the UPPER face `i+1/2` of cell `i`; see also the
face-flux decomposition section in [_historical.md](_historical.md)): `rhoV_cell = rho_cell * v_cell` plotted at
the cell centres; `f_total` = the production numerical face mass flux (Roe here; the release flux is now the SWMF Godunov flux,
Rusanov if selected); `f_ref` = `f_total` at `t = 0`, a **diagnostic baseline**;
`f_eff = f_total - f_ref` = the change in face mass transport since initialization.
Since the reference-free release (`docs/studies/numerics/reference_free_release_recap.md`) the solver
subtracts nothing — `f_total` is what the update sees — so `f_ref`/`f_eff` are now
purely an analysis convenience for separating the evolving transport from the
(round-off-level) initial hydrostatic flux. Older runs in this catalog that predate
the retirement were produced with `eq_wb = 1`, where `f_eff` *was* what the update
saw; their captions are left as recorded.

Layout: 2x2, left column = full captured column, right column = top zoom
(default 2100 km -> domain top), top row = cell `rho V`, bottom row = `f_total`
(solid blue) / `f_ref` (dashed grey, time-independent) / `f_eff` (solid red). The
bottom row **repeats** the top row's `rho V` array as a thick light-grey halo drawn
behind `f_eff` (same curve, same y-scale): agreement shows as a grey fringe around
the red line, and any cell-vs-face difference shows as separation. **Both rows of
a column share one y-scale** (same units, comparable magnitude) and the scale is
fixed over the whole movie, so time and cell-vs-face comparisons are both honest;
the zoom column carries its own limits. The 2130--2150 km release analysis window
is shaded, and a metrics box reports the window mean `rho V`, mean `f_eff`, and the
`roughness(rho V)/roughness(f_eff)` ratio.

The startup print is the audit trail: `f_central + f_diff == f_total` (~1e-9) and
a `rho V` cross-check against the nearest-in-time frame of the main snapshot file
(`mix::MOM`). The former `-div f_ref == eq_residual_mass` check was removed with the
`eq_residual_mass` sidecar column when `eq_wb` left the release. The main-output and `.faceflux`
cadences are independent, so that check reports its own `dt`; it is skipped with a
message for older dumps whose main file is not a 3-row release state.

Each run also gets a static last-frame `*_mass_flux_final.png` (same layout).

```bash
# 4000 s Roe + (ln rho, V, ln p) N500 production run (358 captured records)
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/animate_face_mass_flux.py \
  outputs/model_column/lnp_roe_N500_4000s.txt \
  visualization/model_column/lnp_roe_N500_4000s_mass_flux.mp4 --fps 20

# current release configuration, 100 s smoke capture (11 records)
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/animate_face_mass_flux.py \
  outputs/model_column/release_roe_lnp_N500_100s.txt \
  visualization/model_column/release_roe_lnp_N500_100s_mass_flux.mp4 --fps 4
```

Options: `--zoom LO HI`, `--window LO HI`, `--stride N` / `--max-frames N`
(or `ANIM_STRIDE` / `ANIM_MAX_FRAMES`), `--faceflux PATH` for a sidecar that does
not sit next to the run, `--final-figure PATH` / `--no-final-figure`,
`--no-check-output`. A sidecar whose first record is not `t = 0` is rejected, since
`f_ref` is undefined without it.

---

### `model_column/roe_n1000_lower_chromosphere_sloshing.png`

Six-panel diagnosis figure for the broad low-chromosphere velocity / mass-flux oscillation
in the Roe + `(ln rho, V, ln p)` N1000 4000 s run
(`docs/studies/numerics/roe_n1000_lower_chromosphere_sloshing_diagnosis.md`). Top row: `V(h, t)` space-time
heatmaps for Roe N1000 and Roe N500 over 1600--2130 km, plus a high-cadence (2.8 s)
`.faceflux` capture of the first 220 s that resolves the conduction-launched downgoing
acoustic front (base arrival ~58 s, matching the integral of ds/c_s) and its reflection off
the frozen `V = 0` inner wall. Bottom row: `V(t)` at 1800 km for Roe/Rusanov x N500/N1000
(same mode, ~2x amplitude at N1000, Rusanov included -- the mode is not Roe-specific), mode
amplitude vs height (near-node at the inner wall, peak at ~1830 km, all four runs collapse
above ~2050 km), and the evaporation-window 2130--2150 km mean `rho*V` (uniformly positive,
resolution-consistent). Generated from existing runs only; no new production run.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
python util/release/plot_column_sloshing.py \
  --run "Roe+lnP N1000=outputs/model_column/lnp_roe_N1000_4000s.txt.gamma_diag" \
  --run "Roe+lnP N500=outputs/model_column/lnp_roe_N500_4000s_opt.txt.gamma_diag" \
  --run "Rusanov+lnP N1000=outputs/model_column/acc_lnp_N1000.txt.gamma_diag" \
  --run "Rusanov+lnP N500=outputs/model_column/acc_lnp_N500.txt.gamma_diag" \
  --faceflux outputs/model_column/lnp_roe_N1000_1000s.txt.faceflux \
  --out visualization/model_column/roe_n1000_lower_chromosphere_sloshing.png
```

Supporting conduction-off control run used by the same diagnosis (no figure):

```bash
ISO_NS=1000 ISO_HEAT_FLUX=0 \
CHROMO_T_END=1000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
scripts/run_chromo_realtime.sh outputs/model_column/eqctl_lnp_roe_N1000_1000s.txt \
  full no-ionization model_column - 20.0 no-cooling \
  > outputs/model_column/eqctl_lnp_roe_N1000_1000s.console.log 2>&1
```

---

## Release Godunov cutover — N500 4000 s evolution and base mass-flux diagnosis

> Measured in the pre-cutover float32 build: top-of-domain evaporation numbers in this section are high by of order 30 %. The `double` re-measurement is the storage-precision cutover comparison below.

### Mass-flux panel default changed

`util/release/animate_isentropic.py` now plots the **conservative numerical face mass flux
`f_total`** from the `<run>.faceflux` sidecar in its mass-flux panel, at the face
heights, whenever that sidecar exists. The cell-centred product `rho*V` is only a
reconstruction of the flux: it carries the full cell-centred ripple and is
systematically misleading in the boundary cells, where the ghost closure fixes the
face flux and not the cell-centred product. Set `ANIM_MASS_FLUX=cell` to force the
old panel; without a sidecar the script falls back to `rho*V` and labels the panel
accordingly. **Every movie catalogued before this section (now in [_historical.md](_historical.md)) predates this change and shows
`rho*V`.** New runs should therefore enable `CHROMO_FACE_FLUX_DIAG=1`; choose
`CHROMO_FACE_FLUX_STRIDE` so the capture cadence is at least as fine as
`CHROMO_FRAME_DT` (the animator picks the nearest capture in time and prints the
worst mismatch).

### `lnp_godunov_N500_4000s_evolution.mp4`

The 4000 s N500 release run under the **new release numerics** — SWMF exact-Riemann
Godunov flux at frozen composition, `(ln rho, V, ln p)` MUSCL reconstruction with
MC3/beta=2, Gamma/Saha EOS, R4 outer refinement, 22,000 K physical-face conductive
boundary, physical conduction only, CFL 0.50. No `ISO_RIEMANN` override: this is the
default. Reached `termination=end_time` at 4000 s, step 712,314 (the identical step
count to the Roe run, as expected from the CFL audit), 502 s wall on 12 threads,
941,679,108 exact face solves with **zero fallbacks**. Snapshot cadence 20 s gives
201 frames; face-flux stride 3300 steps gives 217 captures (worst snapshot/capture
time mismatch 9.2 s). This is the direct successor to
`lnp_roe_N500_4000s_evolution.mp4` and the mass-flux panel is now the face flux.

```bash
CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_STRIDE=3300 \
/usr/bin/time -p scripts/run_chromo_realtime.sh \
  outputs/model_column/lnp_godunov_N500_4000s.txt \
  full no-ionization model_column - 20.0 no-cooling \
  > outputs/model_column/lnp_godunov_N500_4000s.console.log 2>&1

ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/release/animate_isentropic.py \
  outputs/model_column/lnp_godunov_N500_4000s.txt \
  visualization/model_column/lnp_godunov_N500_4000s_evolution.mp4 25
```

### `godunov_N500_4000s_base_mass_flux_diagnosis.png`

Six-panel diagnosis of the elevated mass flux near the lower boundary, comparing the
new release run against the old `lnp_roe_N500_4000s` run the earlier movie was made
from. Top row: `f_total` versus height at t = 0, 2000, 4000 s for each run, plus
`|f_diff|` for the release run. Bottom row: base zoom (1600-1640 km) comparing
`f_total` against cell-centred `rho*V` for each run, plus the per-step mass-density
increment measured in float32 ULPs of `rho`. Findings, in order of size:

1. **Most of the old base anomaly was a stale run, not a boundary condition.** The
   `lnp_roe_N500_4000s` run predates the MUSCL-Hancock predictor momentum source.
   Its t = 0 face velocity is 0.7569 m/s, which is exactly `0.5 * dt * g`
   (0.0055258 * 273.95 / 2 = 0.75690), and its t = 0 face mass flux is 7.5e-10
   rather than zero. That is the documented pre-fix defect. The new run starts at
   `max|f_total| = 1.08e-12` and its `f_total` tracks cell-centred `rho*V` to
   sub-percent everywhere above cell 0, where the old run had a uniform 1.8x offset.
2. **A genuine, localized first-cell artifact remains.** At face 0 the dissipative
   part `f_diff` is -1.7e-10, about **1900x** the interior median, and cell 0's
   centred `rho*V` overshoots the face flux by 28 %. This is the known base-cell /
   inner-ghost artifact and it is confined to the single lowest face.
3. **The broad base-to-top gradient is a slow transient, and float32 prevents it
   from relaxing.** Above ~1850 km the t = 4000 s profile is flat to about 10 %,
   i.e. genuinely quasi-steady; the non-flat part is the lowest ~200 km, where
   `f_total` still runs 2.4x the interior mean. It cannot work itself off: the
   per-step mass-density increment `dt * (-df/ds)` is **below half a float32 ULP of
   `rho` in 618 of 660 cells** (median ratio 0.29), so `rho += dt*rhs` rounds back to
   `rho` every step and the density in the lower chromosphere is frozen at
   representation round-off. This also explains the previously recorded and
   unexplained `.faceflux` divergence anomaly in
   `docs/studies/numerics/roe_n1000_lower_chromosphere_sloshing_diagnosis.md` section "Unrelated issue
   found": the captured flux values are right, but their divergence cannot be
   reconciled with a density evolution that float32 will not let happen.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
  util/release/plot_base_mass_flux_diagnosis.py \
  "release Godunov (current code)=outputs/model_column/lnp_godunov_N500_4000s.txt" \
  "Roe, pre predictor-source fix=outputs/model_column/lnp_roe_N500_4000s.txt" \
  visualization/model_column/godunov_N500_4000s_base_mass_flux_diagnosis.png
```

---

## Storage-precision cutover comparison (CURRENT — double is the release)

The re-measured comparison after the release cutover, with precision as a single
coherent knob (`chromosphere::Real`). The `double` leg is the plain release build
with no precision flag; the `float32` leg is the default-off diagnostic build.
Both legs ran to `termination=end_time` at 4000 s with `godunov.fallbacks=0`:
712,029 steps / 941,302,338 exact face solves (double) and 712,314 / 941,679,108
(float32). They are **not** step-identical, because a coherent knob promotes the
`Grid` physical constants too, which perturbs `c_s` at the `1e-8` level and hence
`dt`; the earlier section could claim identical tallies only because it left those
constants in `float`.

**Regression check:** the float32 leg is **byte-identical** to the archived
pre-cutover release run `outputs/model_column/lnp_godunov_N500_4000s.*` in all
three files, so the refactor is semantics-neutral and every difference below is a
genuine precision effect.

```bash
# release (double) — no precision flag
cmake -S . -B build_omp -DCMAKE_BUILD_TYPE=Release -DCHROMO_ENABLE_OPENMP=ON
cmake --build build_omp -j 12
# diagnostic (float32)
cmake -S . -B build_omp_f32 -DCMAKE_BUILD_TYPE=Release \
      -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_STATE_FLOAT32=ON
cmake --build build_omp_f32 -j 12

for leg in release_double:build_omp diag_float32:build_omp_f32; do
  CHROMO_BINARY=$PWD/${leg##*:}/chromo_main \
  CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_STRIDE=3300 \
  scripts/run_chromo_realtime.sh \
    outputs/model_column/${leg%%:*}_N500_4000s.txt \
    full no-ionization model_column - 20.0 no-cooling
done
```

### `precision_cutover_N500_4000s.png`

Same six panels as the withdrawn hybrid-build comparison in [_historical.md](_historical.md) (`f_total(h)` at 2000 s and 4000 s; base zoom of
`f_total` against cell-centred `rho*V`; final velocity profile; column mass budget;
per-step mass increment in ULPs of each leg's storage type), re-run against the
coherent builds. The legend now labels float32 as the diagnostic leg.

| Quantity, N=500/R4, 4000 s | float32 (diagnostic) | **double (RELEASE)** |
| --- | --- | --- |
| `f_total` base/top, t = 2000 s | 2.682 | **1.006** |
| `f_total` base/top, t = 4000 s | 2.811 | **1.003** |
| 1600–1850 km spread of `f_total`, t = 4000 s | 119.5 % | **13.0 %** |
| cells with per-step mass increment below 0.5 ULP | 619 / 660 | **0 / 660** |
| mass budget over 1605–1850 km, predicted / observed | **-214.9** | **0.994** |
| base velocity, t = 4000 s | 1.392 m/s | 0.356 m/s |
| top velocity, t = 4000 s | 13.18 m/s | **9.50 m/s** |
| `f_total` at the top face, t = 4000 s | 3.753e-10 | **2.702e-10** |
| `f_diff` at face 0, as a fraction of `f_total[0]` | 16.30 % | 16.10 % |
| cell 0 `rho*V` / `f_total[0]` | 1.280 | 1.275 |

The gross continuity violation is gone — base/top 2.811 → 1.003, the mass budget
closes, and no cell is frozen. **Two corrections to that withdrawn comparison, both against
overclaiming:** the lower-chromosphere spread is 13.0 %, not 2 % (a broad monotone
decline with its minimum at 1849.7 km, *not* the first-face artifact, relaxing only
slowly — 14.3 % at 500 s to 13.0 % at 4000 s, so the column is still relaxing at
4000 s); and the top-of-domain evaporation observables move by ~28 %, not ~3 %, and
downward, because a lower chromosphere that cannot drain sustains a mass supply the
column does not have. **Any evaporation number from a float32 run is high by of
order 30 %.** The first-interior-face artifact is unchanged and remains open.
Full record: `docs/studies/numerics/state_precision_release_cutover.md`.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
  util/release/plot_precision_control.py \
  outputs/model_column/diag_float32_N500_4000s.txt \
  outputs/model_column/release_double_N500_4000s.txt \
  visualization/model_column/precision_cutover_N500_4000s.png
```

---

## Upper-boundary conduction-wall temperature sensitivity (`docs/studies/boundaries/upper_wall_temperature_sensitivity_recap.md`)

One-parameter sweep of the imposed conduction-wall temperature `ISO_T_TOP` at the 2153 km outer face, with the C7 IC and every other piece of release physics/numerics fixed and `ISO_HYDRO_T_DECOUPLE=1` (so only `Grid::outer_conduction_temperature` moves). Double-precision release build (`build_omp`), N=500/R4, 4000 s. All four figures come from one invocation of `util/release/twall_sweep_diag.py`.

Prerequisite runs (about 3 h wall on 16 cores; ~260 MB of dumps under `outputs/model_column/twall_sweep/`):

```bash
util/run/run_twall_sweep.sh          # broad log sweep 10,15,18,22,30,45,70,100 kK + 70 kK CFL control
SWEEP_CASES="6000:0.50 8000:0.50 11000:0.50 12000:0.50 13000:0.50 14000:0.50" \
    util/run/run_twall_sweep.sh      # refinement across the flow reversal
```

The overlay movie's mass-flux panel additionally needs the conservative face flux, which is a separate opt-in because it costs ~230 kB per capture (~1.6 GB for the sweep at the default ~10 s cadence). `SWEEP_FACE_FLUX=1` adds `CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_STRIDE=1780` to every case; the 13 accepted walls only (6 kK aborts, and the 70 kK CFL 0.25 control is not animated). About 1 h wall — the 100 kK CFL 0.25 case runs ~3.1 M steps and is the long pole by a wide margin. Only the `.faceflux` sidecars are kept; the rest of the re-run is discarded after the byte-identity check below.

```bash
SWEEP_FACE_FLUX=1 SWEEP_CASES="8000:0.50 10000:0.50 11000:0.50 12000:0.50 \
13000:0.50 14000:0.50 15000:0.50 18000:0.50 22000:0.50 30000:0.50 45000:0.50 \
70000:0.50 100000:0.25" util/run/run_twall_sweep.sh outputs/model_column/twall_sweep_ff
# the diagnostic must be inert -- verify, then keep only the new sidecars
for f in outputs/model_column/twall_sweep_ff/*.txt.gamma_diag; do
    b=$(basename "$f")
    [ "$(md5 -q "$f")" = "$(md5 -q "outputs/model_column/twall_sweep/$b")" ] \
        && echo "IDENTICAL $b" || echo "DIFFERS   $b"
done
mv outputs/model_column/twall_sweep_ff/*.faceflux outputs/model_column/twall_sweep/
rm -rf outputs/model_column/twall_sweep_ff
```

### `twall_sweep_N500_4000s_V_of_t.png`, `..._V_of_h.png`, `..._summary.png`, `..._support.png`
Per-run colours are **rank-spaced** on `turbo` with a luminance cap, not log-value-spaced on `plasma`: the sweep clusters six walls between 10 and 15 kK to resolve the flow reversal, and value spacing collapsed those six into one indistinguishable dark purple. Affects `V_of_t`, `V_of_h` and `support`; the `summary` panels use fixed per-height colours and are unchanged.

`V(t)` at 1800/2000/2100 km and the top cell (symlog, so the sub-m/s cold-wall cases stay visible next to the 3 km/s 100 kK transient); quasi-steady `V(h)` full column, resolved-TR zoom and symlog; the summary panel set (V, `rho V`, `q_wall`, `T_top`, log-log upflow branch with the local exponent `dlnV/dlnT`, and the steady enthalpy-flux budget); and the supporting `rho V(h)`, `T(h)` and column-mass-drift diagnostics. The 6 kK run is excluded automatically because it aborts at t≈380 s below the EOS temperature floor; the 100 kK point is supplied explicitly because it needed CFL 0.25.

```bash
.venv/bin/python util/release/twall_sweep_diag.py --cfl 0.50 \
    --run 100000=outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt \
    --heights 1800 2000 2100 --avg-window 500 \
    --out-dir visualization/model_column --prefix twall_sweep_N500_4000s \
    --json outputs/model_column/twall_sweep/twall_sweep_metrics.json
```

### `twall_sweep_N500_4000s_overlay.mp4`

The whole sweep animated: all 13 accepted wall temperatures overlaid on the **same six panels as `lnp_godunov_N500_4000s_f64_evolution.mp4`** — `T(h)`, `V(h)`, `q_parallel(h)` on the top row, `rho(h)`, `p(h)`, mass flux `(h)` on the bottom — so the sweep reads against the single-run release movie panel for panel. 201 frames (20 s cadence, 0–4000 s), 20 fps, 10 s. Runs are keyed by a **discrete colorbar with one equal-width band per run**, built from the actual line colours: `wall_colors()` spaces the colours by *rank* and caps their luminance, so a continuous `LogNorm` bar drawn from `turbo` would show colours no line has and place them at the wrong `T_wall`. Equal-width bands are the honest rendering of a rank-spaced assignment — the bar is a legend with an axis, and its width is deliberately not linear in `T_wall`. The 22 kK release baseline is drawn thicker and its tick label is bold. Run discovery, run acceptance (6 kK excluded on `termination != end_time`; 100 kK supplied explicitly because it needed CFL 0.25) and the rank-spaced `turbo` colours are imported from `util/release/twall_sweep_diag.py`, so the movie and the four static figures always describe the same set of runs.

**The mass-flux panel is the conservative numerical FACE flux `f_total`**, read from the `.faceflux` sidecars and plotted at the face heights — the quantity the continuity row is actually differenced from, as in `util/release/animate_isentropic.py`. Producing it required re-running the sweep with the read-only face-flux diagnostic (`SWEEP_FACE_FLUX=1`, new opt-in flag on `util/run/run_twall_sweep.sh`; see the prerequisite block above). Enabling the diagnostic changes no numerical result, and this was verified rather than assumed: **all 13 re-run `.gamma_diag` files are byte-identical to the original sweep's**, so the recap's metrics stand unchanged and the re-run doubles as an exact reproducibility check. Set `ANIM_MASS_FLUX=cell` to force the cell-centred `rho*V` panel instead; the panel title always states which one is drawn. The two grids do not coincide — the capture cadence is a fixed number of *steps* while snapshots are in *seconds*, and a hotter wall runs a shorter `dt`, so the per-run capture count ranges 253 (8 kK) to ~1700 (100 kK, CFL 0.25) — so the script picks the nearest capture per run per frame and prints the worst frame-to-capture mismatch (6.1 s, set by the 8 kK run's 15.8 s cadence).

**The pressure panel is `p(h,t)/p_22kK(h,t)`, pressure relative to the 22 kK release run, not `p` in Pa** — a linear axis with the release baseline as the flat unit reference line. Absolute `p` separates the sweep at *no* y scale, and this was measured rather than assumed: `p(h)` falls 6.6x over the column while the wall-to-wall spread is at most 12 % at any height (1.001 at 1600 km, peaking at 1.121 at 2104 km), so all 13 lines sit within a line width of each other on a log-Pa axis. Zooming y does not fix it — `p` falls that same 12 % in ~30 km, so a y zoom crops x with it and leaves the fan packed against the top of the domain with the lower 400 km off-scale. Dividing by the release run removes the shared stratification *and* the shared early acoustic transient, leaving a 0.90–1.10 fan separated over the whole 1600–2153 km: at 4000 s the 100 kK wall reaches 0.90 near 2100 km while the 8 kK wall peaks 1.6 % above it. Set `ANIM_P_PANEL=abs` for the old absolute log-Pa panel; the panel title always states which one is drawn.

One deliberate deviation from the single-run movie, forced by overlaying responses that span four decades: **`V`, `q` and the mass flux are symlog, and `V` is in m/s rather than km/s.** The 8 kK wall drains at −0.4 m/s while the 100 kK wall evaporates at +792 m/s; on the single-run movie's linear km/s axis every wall below 45 kK is a flat line on zero. `T` and `rho` are logarithmic for the same reason (`T_top` spans 6.7 → 100 kK). `q` gets a 10x wider symlog linear band than `V` and `f_total` because its range is set by the TR/boundary spike (~100 W/m²) while the chromospheric interior sits near 1e-3.

**How much the panel choice actually changes, measured at t = 4000 s** (do not eyeball this from two renders with different y-limits — that is what produced an earlier, wrong claim here that the residual ripple was absent from `f_total`; the ripple is present in both and identical to a ratio of 1.00):

| | 8 kK | 22 kK (release) | 70 kK | 100 kK |
|---|---|---|---|---|
| base cell / base face, `rho*V[0] / f_total[0]` | 1.268 | **1.275** | 1.282 | 1.292 |
| interior max `\|1 - rho*V/f_total\|` | 0.96 % | 1.07 % | 6.25 % | 13.7 % |
| top cell / top face | 1.000 | 0.999 | 0.999 | 1.000 |

Two things follow. First, the base-cell overshoot is **1.27–1.29 at every wall temperature**, reproducing the 1.28 recorded for the release run and showing that the known first-interior-face artifact is essentially independent of the boundary forcing — consistent with its being a ghost/reconstruction closure issue rather than anything the wall drives. Second, the interior discrepancy is not negligible for hot walls: at 100 kK the cell-centred product misstates the conserved flux by up to 14 % somewhere in the interior, so plotting `f_total` is a correctness improvement for the sweep and not a cosmetic one.

What the movie shows that the static figures do not: the 6.7 kK chromospheric plateau and the `rho` stratification below ~2100 km are visually identical across all 13 walls for the entire 4000 s — the *thermal* response is confined to the resolved TR above ~2115 km — while the hot-wall cases ring through a large early acoustic transient before settling onto the ordered, still slowly drifting quasi-steady fan. The relative-pressure panel is the one place where that "identical below the TR" reading is shown to be a limit of the log axis and not a fact: once the release run is divided out, the wall signal is visible as a smooth, monotone fan reaching down through the whole column (at 4000 s the 100 kK wall is 10.0 % below the release baseline at 2100 km, still 4.2 % below at 1900 km and 1.0 % below at 1700 km, pinned back to 0.08 % only at the 1600 km truncation face itself, where the reservoir closure fixes the state) — the boundary forcing is felt hydrostatically far below the TR, where the temperature is untouched.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/animate_twall_sweep.py \
    --cfl 0.50 \
    --run 100000=outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt \
    --out visualization/model_column/twall_sweep_N500_4000s_overlay.mp4 --fps 20
```

### `twall_sweep_N500_4000s_overlay_p_abs.mp4`

The same movie with the **absolute log-Pa pressure panel** (`ANIM_P_PANEL=abs`), kept alongside the relative-pressure version above. Identical in every other panel and in every frame. Keep it for two reasons: it is the only rendering that shows the actual pressure stratification of the sweep in physical units, and it is the control that makes the relative panel's case — the 13 walls are indistinguishable here, which is the measurement, not an artifact of a badly chosen y limit (see the numbers in the entry above). Use this one when the question is "what is the pressure", the relative one when the question is "how does the wall change the pressure".

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl ANIM_P_PANEL=abs .venv/bin/python util/release/animate_twall_sweep.py \
    --cfl 0.50 \
    --run 100000=outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt \
    --out visualization/model_column/twall_sweep_N500_4000s_overlay_p_abs.mp4 --fps 20
```

### `twall_sweep_N500_4000s_overlay_p_ic.mp4`

The same movie with the pressure panel **relative to the shared C7 initial condition**, `p(h,t)/p(h,0)` (`ANIM_P_PANEL=ic`), kept alongside the baseline-relative and absolute versions. Every other panel and every frame are identical. Since all runs start from the same IC, the denominator is the same for every run and none is privileged. The 22 kK release run is a curve like the others (about 1 % below its IC under the TR at 4000 s), and each run's bend sits at **its own** transition region. That is why this version exists. In the baseline-relative panel, the sharp turn near the top (0.90 at 2104 km for 100 kK) reads like an upper-boundary inconsistency, but it is a TR-height difference. The hot TR is nearly isobaric (scale height ≳ 10³ km against ~230 km in the chromosphere), and the 100 kK wall has moved its TR base from ~2138 km down to ~2104 km. At fixed height, that compares hot, flat plasma with the release run's still-stratified chromosphere. Compared at equal column mass instead of equal height, every wall from 8 to 100 kK matches the release pressure to within 0.16 %, and `(p − p_top)/(g·m)` is 1.000–1.004: the pressure is hydrostatic everywhere below the TR. The small ~0.2 % ripple on the 100 kK curve at 1650–1950 km also appears in its `V` and `f_total`, so it is a flow feature of that run, not something the normalization creates.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl ANIM_P_PANEL=ic .venv/bin/python util/release/animate_twall_sweep.py \
    --cfl 0.50 \
    --run 100000=outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt \
    --out visualization/model_column/twall_sweep_N500_4000s_overlay_p_ic.mp4 --fps 20
```

### `twall_sweep_N500_p_abs_top_t2000.png`

Absolute pressure in Pa, not normalized, for all 13 accepted walls at t = 2000 s, restricted to h ≥ 2100 km. Over the full column the walls sit within a line width of each other, but the top ~50 km carries the whole spread, so zooming in separates them without dividing by any reference. A dot marks each run's TR base (first cell with T > 7.5 kK). The dots show that the "kink" in the normalized panels is simply where each run's p(h) turns from the steep chromospheric slope (scale height ~230 km) to the nearly flat hot-TR slope: 2116 km for 100 kK, 2127 km for 70 kK, ~2139 km for the walls near the release. Every p(h) decreases monotonically with height apart from one feature: at the 70 and 100 kK TR base there is a 1–2e-5 Pa (~0.2 %) notch over two to three cells, where T rises 2.4x per 0.28 km cell. That is about 3x the run's whole momentum-flux budget `rho V^2` (6.7e-6 Pa at 100 kK), so it is a feature of the under-resolved thermal front, not the acceleration of the evaporating gas.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/plot_twall_sweep_p_top.py \
    --time 2000 --hmin 2100 \
    --run 100000=outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt \
    --out visualization/model_column/twall_sweep_N500_p_abs_top_t2000.png
```

### `twall_relaxation_phases.png`, `twall_relaxation_spacetime.png`, `twall_relaxation_timescales.png`

How the evaporation mass flux relaxes and then decays across the sweep, each stage against its analytic estimate. Study: `docs/studies/boundaries/twall_mass_flux_relaxation_timescales.md`. Post-processing only: the 13 accepted sweep runs with their `.faceflux` sidecars (prerequisites above), plus the current single-T release run `outputs/model_column/release_singleT_N500_4000s` for the flux-profile check. Metrics go to `outputs/model_column/twall_sweep/twall_relaxation_metrics.json`.

**`twall_relaxation_phases.png`** (six panels):
- (a) **wave phase** — `F_base/F_top − 1` with multiples of the one-way acoustic time (dotted) and a ±5 % band;
- (b) **TR relaxation** — `F·L_TR/J`, with `J` the first-principles steady-TR integral;
- (c) **the drift** — `F` keeps falling while `F·L_TR` stays flat, both normalized at 1000 s;
- (d) **the contraction law** — `L_TR` growth against evaporated mass `Φ`, with the linear coefficient `K`;
- (e) **the adiabatic test** — `T` change at fixed column mass, measured (lines) against the parcel-lift adiabat (circles), cut 10 km below the TR base;
- (f) **the flux profile** — `f_total(z)/f_top`, averaged over 3800–4000 s to suppress the hot-wall ripple, against the contraction prediction (dashed). 100 kK is omitted there because its ±30 % spatial ripple dominates even after averaging; the black line is the release run.

**`twall_relaxation_spacetime.png`**: `f_total(z,t)/f_top(t)` for 22 and 100 kK over 0–900 s, with the ±`c_eq` characteristics launched from the top at t = 0.

**`twall_relaxation_timescales.png`**: per wall temperature (≥ 15 kK):
- when the wave phase is over (5 %);
- when the TR is relaxed, with the fitted `τ_TR`;
- the flux-drift e-folding `1/r` at 4000 s, measured and from contraction alone.

Walls ≤ 14 kK sit at the flow reversal, where `F → 0` makes every ratio meaningless, so they are left out of the ratio panels and the timescale figure. Colours are the sweep's rank-spaced `wall_colors()`.

```bash
.venv/bin/python util/release/twall_relaxation_diag.py \
    --run 100000=outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt
```

## Outer-face closure swap: fixed `T_hydro` + Neumann pressure (`docs/studies/boundaries/outer_bc_pressure_temperature_swap_recap.md`)

The release outer face imposes its one admissible external condition as a **fixed reservoir back-pressure** and lets the hydro ghost temperature free-float. These two movies put the opposite choice on the same axes: hydro ghost temperature **pinned** to the 22 kK conduction wall (`ISO_OUTER_T_HYDRO_WALL=1`) and the pressure **Neumann**, in both readings — `ISO_OUTER_P_NEUMANN=1` (`dp/ds = 0`) and `=2` (`dp/ds = -rho g`, the release hydrostatic ladder re-anchored on the live top cell) — plus the control that pins the temperature and keeps the release back-pressure, which isolates the pressure closure as the cause. Canonical release configuration otherwise (N=500/R4, double state, CFL 0.50, physical conduction only); the release baseline is the 22 kK run of the section-16 sweep, so the IC, mesh, flux and snapshot cadence are shared exactly.

Prerequisite runs (~45 min wall total on 12 threads; ~1.3 GB under `outputs/model_column/outer_bc_swap/`):

```bash
# the three 4000 s cases (the release baseline comes from the conduction-wall temperature sweep above)
for CASE in "tfix_pneumann:1" "tfix_phse:2" "tfix_pref:0"; do
    ISO_OUTER_T_HYDRO_WALL=1 ISO_OUTER_P_NEUMANN="${CASE##*:}" \
    CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
    CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=2000 \
    CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_STRIDE=1780 \
    scripts/run_chromo_realtime.sh \
        "outputs/model_column/outer_bc_swap/${CASE%%:*}_4000s.txt" \
        full no-ionization model_column - 20.0 no-cooling \
        > "outputs/model_column/outer_bc_swap/${CASE%%:*}_4000s.console.log" 2>&1
done

# 1 s-cadence 200 s probes for the transient movie (release, mode 1, mode 2)
for CASE in "release_probe200:" "probe200:1" "hse_probe200:2"; do
    N="${CASE##*:}"
    ${N:+ISO_OUTER_T_HYDRO_WALL=1} ${N:+ISO_OUTER_P_NEUMANN=$N} \
    CHROMO_T_END=200 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=1 \
    scripts/run_chromo_realtime.sh \
        "outputs/model_column/outer_bc_swap/${CASE%%:*}.txt" \
        full no-ionization model_column - 20.0 no-cooling \
        > "outputs/model_column/outer_bc_swap/${CASE%%:*}.console.log" 2>&1
done
```

### `outer_bc_closure_swap_4000s.mp4`

Four runs overlaid on the **same six panels as the release evolution movie** — `T(h)`, `V(h)`, `q_parallel(h)` on the top row, `rho(h)`, `p(h)`, `f_total(h)` on the bottom. 201 frames (20 s cadence, 0–4000 s), 20 fps, 10 s. Colours are assigned by position so a subset renders in the same colours as the full set: release dark grey and thick, `dp/ds = 0` red, `dp/ds = -rho g` blue, the `T`-pinned control thin grey drawn last so its coincidence with the release line is visible.

Panels follow `util/release/animate_twall_sweep.py` with three deliberate differences, all forced by comparing *closures* instead of one parameter: runs are keyed by a **legend**, not a colorbar (all four share the 22 kK wall, so `T_wall` no longer identifies a line); **`p` is absolute log Pa**, not divided by the release run, because a closure swap moves `p` by factors of several and a ratio panel would only rescale that; and a run that did **not** reach `CHROMO_T_END` is kept rather than excluded, with its termination shown in the legend, because for a boundary-closure experiment an early abort is the result. Here all four reached 4000 s with zero Godunov fallbacks, so no legend carries a tag.

The mass-flux panel is the conservative face flux `f_total` from the `.faceflux` sidecars (`ANIM_MASS_FLUX=cell` forces the cell-centred `rho*V`; the panel title always says which). Its symlog linear band is 10x wider than the sweep overlay's because the collapse run sits at ~1e-7 while the release sits at ~1e-10, and the tighter threshold collided the +/- decade labels around zero.

What it shows: the `dp/ds = 0` run leaves the release solution within the first 60 s and never returns — by 200 s the transition region has marched down past 1900 km, and from ~700 s the whole domain is a stationary 23.7 kK fully ionized downdraft at the `ISO_VCAP` cap carrying `f_top = f_base = -4.13e-07` kg m⁻² s⁻¹, 1530x the release flux and pointing the wrong way, with half the column mass gone. The `dp/ds = -rho g` run instead tracks the release in `T` and `rho` for the whole 4000 s and fails quietly: its `p` creeps upward, its `f_total` sits an order of magnitude below the release's and disagrees in sign between the top and base faces. The control is the release line.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/animate_outer_bc_compare.py \
    --out visualization/model_column/outer_bc_closure_swap_4000s.mp4 --fps 20
```

### `outer_bc_singleT_vs_release_1500s.mp4`

The release outer temperature BC (single T: one fixed 22 kK temperature in the hydro ghost cells, which is also the conduction Dirichlet datum at the ghost centre, `ISO_HYDRO_T_DECOUPLE=0`) overlaid on the previous decoupled release (free `T_hydro`, 22 kK conduction datum at the physical face). Same six panels as `outer_bc_closure_swap_4000s.mp4`. 0–1500 s from the two 4000 s runs `outputs/model_column/release_singleT_N500_4000s` and `outputs/model_column/twall_sweep/twall_22kK_cfl0.50`, 76 frames (20 s snapshot cadence), 8 fps, 9.5 s. Previous release dark grey and thick, release thin red. The mass-flux panel is the conservative face flux `f_total`. Face-flux captures are step-strided, so that panel can lag the other five by up to 3.9 s.

What it shows: the lines coincide in every panel. V_top is 4 % lower under single T at 20 s, and the gap settles to about 0.9 % from 400 s on (11.98 vs 12.10 m/s at 1500 s). T_top is 80–200 K cooler, because the conduction datum sits a full cell away instead of half a cell. T differences above 0.1 % are confined to the top ~13 km. Analysis: `docs/studies/boundaries/outer_temperature_bc_single_t_recap.md`.

```bash
ANIM_WORKERS=14 MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/animate_outer_bc_compare.py \
    --preset singleT --t-max 1500 \
    --out visualization/model_column/outer_bc_singleT_vs_release_1500s.mp4 --fps 8
```

### `outer_bc_closure_swap_transient_200s.mp4`

The same six panels over **0–200 s at 1 s cadence** (201 frames, 12 fps, 17 s) for the release and the two swap variants. The 20 s cadence of the 4000 s movie gives the collapse only three frames, which is where the entire mechanism lives: both Neumann variants launch a **downward-propagating compression front** from the outer face at t = 0, and this movie resolves it crossing the column (at t = 40 s it sits near 1780 km with everything above it falling). Their fates then separate — `dp/ds = -rho g` decays back to a few m/s while `dp/ds = 0` keeps accelerating to the cap.

The mass-flux panel is the cell-centred `rho*V` here, not `f_total`: the 1 s-cadence probes were run without the face-flux sidecar, and the script's all-or-nothing rule then draws `rho*V` for every run rather than mixing two different quantities on one axis. The panel title states it.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl ANIM_MASS_FLUX=cell .venv/bin/python \
  util/release/animate_outer_bc_compare.py \
  --run 'release: fixed $p_{\rm ref}$, free $T_{\rm hydro}$=outputs/model_column/outer_bc_swap/release_probe200.txt' \
  --run 'swap: $T_{\rm hydro}$ pinned, $dp/ds=0$=outputs/model_column/outer_bc_swap/probe200.txt' \
  --run 'swap: $T_{\rm hydro}$ pinned, $dp/ds=-\rho g$=outputs/model_column/outer_bc_swap/hse_probe200.txt' \
  --out visualization/model_column/outer_bc_closure_swap_transient_200s.mp4 --fps 12
```

### `outer_bc_twall_response.png` and `outer_bc_twall_22kK_vs_100kK.mp4`

Wall-temperature study **under** the unanchored `dp/ds = -rho g` closure: the 2x2 of {release outer closure, `ISO_OUTER_T_HYDRO_WALL=1 ISO_OUTER_P_NEUMANN=2`} x {22 kK, 100 kK}. Both unanchored runs are **CFL 0.25**, because 100 kK fails the implicit conduction Newton solve at CFL 0.50 under this closure exactly as it does under the release one; the release counterparts come from the section-16 sweep (22 kK at CFL 0.50, 100 kK at CFL 0.25) and are joined the same way the sweep joins its two CFL groups. The join is measured, not assumed: the unanchored 22 kK case exists at both CFL values (`tfix_phse_4000s.txt` is CFL 0.50) and differs by 8 % in `V_top`.

```bash
# the two unanchored 4000 s runs (~20 min and ~75 min wall, run concurrently at 6 threads)
for T in 22000 100000; do
  ( OMP_NUM_THREADS=6 CHROMO_CFL=0.25 ISO_T_TOP=$T \
    ISO_OUTER_T_HYDRO_WALL=1 ISO_OUTER_P_NEUMANN=2 \
    CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
    CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=2000 \
    CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_STRIDE=3560 \
    scripts/run_chromo_realtime.sh "outputs/model_column/outer_bc_swap/hse_$((T/1000))kK_cfl0.25.txt" \
      full no-ionization model_column - 20.0 no-cooling \
      > "outputs/model_column/outer_bc_swap/hse_$((T/1000))kK_cfl0.25.console.log" 2>&1 ) &
done; wait
```

**The PNG** is the time-domain half of the answer, which the profile movie cannot show: (a) `V_top(t)` symlog, (b) `p_top(t)/p_top(0)` — the level the release BC pins and this closure does not, (c) `M(t)/M(0)`, (d) quasi-steady `V(h)`. Colour carries `T_wall` and the dash pattern carries the closure, one visual channel per varied quantity, which a colour-only scheme cannot do for a 2x2. The legend is a single figure-level legend under all four panels: inside panel (a) it sits on top of the 100 kK curves the panel exists to show. Panel (b) is the point of the figure — +128 % and still rising linearly at 4000 s for the unanchored 100 kK case, against a flat 1.000 for both release runs.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/outer_bc_twall_diag.py
```

**The movie** puts the same four runs on the six evolution panels, 201 frames (20 s cadence, 0–4000 s), 20 fps, with the same colour/dash encoding. It uses the `twall` preset of `util/release/animate_outer_bc_compare.py`, which differs from the `closure` preset in two catalogued ways: the line styles above, and a 10x wider symlog linear band on the mass-flux panel (`f_linthresh_frac`), because the 100 kK release run reaches ~5e-9 while the 22 kK unanchored run works at ~1e-11 and the narrower band puts eight labelled decades on one short axis.

What it shows: the unanchored 100 kK run (red dashed) rings violently for the first ~250 s (`V_top` swinging ±1.5 km/s), then settles into a *two-cell* column — evaporation above ~1900 km, drainage below it — with the transition region parked at 1900 km instead of 2100 km and 1900 km heated from 6.5 kK to 21.9 kK. The release 100 kK run (red solid) keeps the same chromosphere as its 22 kK counterpart and puts the entire wall response into the TR above 2100 km.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/animate_outer_bc_compare.py \
    --preset twall \
    --out visualization/model_column/outer_bc_twall_22kK_vs_100kK.mp4 --fps 20
```

---

## Analytic steady-evaporation lookup table (`docs/studies/boundaries/analytic_evaporation_J_table.md`)

### `analytic_J_table.png`

The steady-TR invariant `J = F·L_TR = ∫ κ/(h − h_c) dT` of the release equations, with `h_c = h_Saha(T_c, p_top)`, tabulated over `T_top` 15–200 kK, `p_top` 1e-3–10 Pa and `T_c` 6.0–7.0 kK. Four panels:
- (a) `J(T_top)` at five pressures, `T_c` = 6.7 kK, with the sweep's measured `F·L_TR` at 4000 s (triangles, `p_top` = 0.0103 Pa); grey = beyond the validated 100 kK;
- (b) `J(p_top)/J(0.0103 Pa)` at five `T_top`;
- (c) `J(T_c)/J(6.7 kK)` at 0.0103 and 1 Pa;
- (d) the local exponent `∂lnJ/∂lnp_top`, with the one line the sweep validated and the dashed contour beyond which `J` depends on the TR-base marker by more than 20 %.

Post-processing only; reads the accepted sweep runs and `outputs/model_column/twall_sweep/twall_relaxation_metrics.json` for the end-to-end check. The table goes to `outputs/model_column/analytic_J_table/` (`J_table.npz`, `J_table.csv`, `J_table_summary.json`). Runtime ~6 s.

```bash
.venv/bin/python util/release/analytic_evaporation_table.py
```

## (T_wall, p_top) evaporation sweep against the analytic closures (`docs/studies/boundaries/tp_sweep_evaporation_closure_validation.md`)

The runs come from `util/run/run_tp_sweep.sh`: 5 wall temperatures (18–70 kK) × 6 pressures (0.003–0.08 Pa), each 5000 s, plus an N=1000 control. They live on the mac-mini (6.4 GB).
- `util/release/tp_sweep_compare.py` runs there. It reduces them to `outputs/model_column/tp_sweep/tp_sweep_compare.json`, `tp_sweep_series.npz` and `compare.log`. Those three files are copied back.
- `util/release/plot_tp_sweep_vs_analytic.py` redraws its two figures from the compact files alone.

### `tp_sweep_sim_vs_analytic.png`

Simulated against analytic evaporation, as curves over `p_top`, one colour per `T_wall`. The six panels:
- (a) `J = F·L_TR`: simulation as markers, the steady-TR integral with the nominal `T_c` = 6.7 kK as lines;
- (b) parity of `F_sim` against `F_E = q_wall/(h_top − h_c)` and against `F_J = J/L_TR`, both with the nominal `T_c`;
- (c) measured over nominal chromospheric enthalpy `h_c`;
- (d)–(f) `F_sim/F_E` with the measured `h_c` (= ER), `F_sim/F_E` with the nominal `T_c`, and `F L_TR/J` with the nominal `T_c`.

Circles are read over 300 s from `t_ready` and are open when the run was not ready before 5000 s. Crosses are read over 4700–5000 s. The star is the N=1000 control.

### `tp_sweep_time_series.png`

One column per `T_wall`, one colour per `p_top`. Three rows:
- the energy ratio `ER(t)`, with a dot at `t_ready`;
- the top-cell flux divided by its median after 1000 s;
- `std_t F(z)/|mean_t F(z)|` over faces after 1000 s, which locates the upper-chromosphere unsteadiness of the hot, high-pressure corner.

```bash
.venv/bin/python util/release/plot_tp_sweep_vs_analytic.py --dir outputs/model_column/tp_sweep \
    --out-dir visualization/model_column
```

### `tp_corner_test.png`

Follow-up at the unsteady 70 kK, 0.08 Pa corner. Four runs, each over 0–1500 s:
- the sweep run;
- N=1000;
- CFL 0.25;
- the R4 grade moved to 1931–1949 km (`ISO_REFINE_S_LO_KM=350`).

The panels:
- (a) `std_t F(z)` over 500–1500 s, divided by the column flux scale `median_z |mean_t F|`, with each case's grade band shaded and its TR base dotted;
- (b) the top-face flux divided by its median;
- (c)–(f) the space-time maps `F(z,t)` over 1900–2153 km, in the same units.

The runs are in `outputs/model_column/tp_sweep_corner/` on the mac-mini, launched by its `launch.sh`. The script reads the raw runs and also writes `corner_test.json`, so it is run there.

```bash
.venv/bin/python util/release/tp_corner_test_compare.py
```

### `tp_sweep_flux_compare.png`

The same comparison as heat maps over the (T_wall, p_top) grid. Each run is read over 300 s starting at its `t_ready`; `*` marks a run that was not ready before the end. The six panels:
- `F_sim`;
- `F_E = q_wall/(h_top − h_c)` with the nominal 6.7 kK `T_c`;
- `F_sim/F_E`, nominal `T_c`;
- `F_sim/F_E`, measured `h_c`;
- `F L_TR/J`, nominal `T_c`;
- `t_ready`.

This script also writes `<dir>/tp_sweep_compare.json`, `<dir>/tp_sweep_series.npz` and, redirected, `compare.log`. It needs the raw runs, so it is run on the mac-mini.

```bash
.venv/bin/python util/release/tp_sweep_compare.py --dir outputs/model_column/tp_sweep \
    --out-dir visualization/model_column --workers 8 > outputs/model_column/tp_sweep/compare.log
```

## Anatomy of the analytic J integral

The 2026-09-29 tp-sweep deck crops made from these figures are in [presentations.md](presentations.md).

### `analytic_J_integrand.png`

What the steady-TR invariant `J = ∫_{T1}^{T_top} κ/(h − h_c) dT` is made of, at `T_top` = 70 kK, `T_c` = 6.7 kK, and `p_top` = 0.003, 0.0103 and 0.08 Pa. Same Saha EOS, conductivity and quadrature as `twall_relaxation_diag.j_theory`. Three panels:
- (a) the Saha `x(T)` (solid) and the enthalpy jump `(h − h_c)/(h_top − h_c)` (dashed);
- (b) the integrand `κ/(h − h_c)`, with its log divergence at `T_c` and its Spitzer `T^{3/2}` tail;
- (c) the running integral `J(T)/J(T_top)`, which in a steady TR is `(s − s_1)/L_TR`, i.e. the steady TR profile.

It also prints `J/J_Spitzer` and the share of `h_top − h_c` gained between 6.7 and 10 kK (0.13 / 0.21 / 0.32). No run needed; runtime ~2 s.

```bash
.venv/bin/python util/release/plot_J_integrand.py
```

---

## (T_wall, p_top) sweep movies, one per p_top (`visualization/model_column/tp_sweep_movies/`)

### `tp_sweep_p<P>_N500_5000s_overlay_p_abs.mp4`, P = 0.003, 0.005, 0.0103, 0.02, 0.04, 0.08

One movie per top pressure of the (T_wall, p_top) sweep above, in the layout of `twall_sweep_N500_4000s_overlay_p_abs.mp4`: `T`, `V`, `q_parallel` on top, `rho`, absolute log-Pa `p`, and the conservative face flux `f_total` below. Each movie overlays the five walls (18, 25, 35, 50, 70 kK) at that `p_top`, keyed by the same discrete rank-spaced colorbar. 251 frames at a 20 s cadence over 0–5000 s, 20 fps (12.6 s). The y limits are pooled over each movie's own five runs, so they differ between pressures; compare shapes across movies, and read magnitudes off the axes. The sweep has no 22 kK run, so no line is drawn thick, and the pressure panel is absolute (`ANIM_P_PANEL=abs`); the release pressure is the `0.0103` movie.

The raw runs are on the mac-mini, so the movies are rendered there (six in parallel, 2 frame workers each) and copied back:

```bash
printf '%s\n' 0.003 0.005 0.0103 0.02 0.04 0.08 | xargs -P 6 -I{} env \
    MPLCONFIGDIR=/tmp/chromosphere2026-mpl ANIM_P_PANEL=abs ANIM_MAX_FRAMES=251 ANIM_WORKERS=2 \
    .venv/bin/python util/release/animate_twall_sweep.py --dir outputs/model_column/tp_sweep \
    --p-top {} --fps 20 \
    --out visualization/model_column/tp_sweep_movies/tp_sweep_p{}_N500_5000s_overlay_p_abs.mp4
```

---

## MC3 branch fix and experimental hydrostatic reconstruction (`visualization/model_column/`)

Study: `docs/studies/numerics/hydrostatic_reconstruction_and_mc3_fix.md`. Inputs are the A/B/C runs under `outputs/model_column/hse_recon/` (A = binary built at `f01189b`, B = MC3 branch fixed, C = B + `ISO_HSE_RECON=1`); the run commands are in the study's Reproduce section.

### `hse_recon_static_residual.pdf`

Column at rest, conduction off, outer ghost matched to the IC (`ISO_T_TOP=0`), 20 s: log-scale bars of the t = 0 max|F^ρ| and the 20 s max|V| for A, B and C at 500 and 1000 cells. C sits at round-off.

### `hse_recon_evaporation.pdf`

Conduction-driven release runs, 1000 s, 500 (solid) and 1000 (dashed) cells: mean V over 2130–2150 km and relative column-mass change against time for A, B and C.

```bash
.venv/bin/python util/release/plot_hse_recon_comparison.py
```

### `hse_recon_evaporation_bars.pdf`

The same 500-cell runs at t = 1000 s as bars for A, B and C: mean V over 2130–2150 km [m/s] and the dissipative share |f_diff/f_total| of the base-face mass flux [%] (16.1, 17.4, 0.1).

```bash
.venv/bin/python util/release/plot_hse_recon_comparison.py
```
