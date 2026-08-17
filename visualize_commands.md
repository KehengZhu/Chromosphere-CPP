# Visualization commands

Catalog of the exact command that regenerates every file under `visualization/`.

**Conventions**

- Run every command **from the project root** (`/Users/zkeheng/SWMFSoftware/Chromosphere2026`).
- The interpreter is the project venv: **`.venv/bin/python`** (see `README.md`).
- All plotting/animation scripts live in `util/` (two physics-reference scripts live in `visualization/reference/` and write next to themselves). Outputs land in `visualization/` (gitignored).
- **Both `outputs/` and `visualization/` are organized into per-scenario subdirs** mirroring `scenarios/`: `model_c7/`, `model_column/` (the unified column — absorbs the old `iso_*`/`gentle_*` runs), `model_flare/`, `analytic_canopy/`, `loops/`, `pfss/`, `events/` (real-event campaigns: 2024-08-01 AR13768, SOL2014, ensembles, `event_te/`, `sol2014/`, …), and `two_fluid/`. `visualization/` also has `reference/` (physics-reference plots). Dumps/figures produced by now-deleted diagnostic code (`dynamic_hse_*`, `cond_hydro_*`) and retired analytic-isentrope stages live under `_archive/`.
- Input simulation dumps live under `outputs/<scenario>/`; loop/field-line geometry `.dat` files live under `scenarios/data/`; magnetograms under `util/data/`. Some `outputs/` subdirs are gitignored — adjust the input paths below to your actual run dumps where noted.
- **Keep this file in sync** (per `CLAUDE.md`): whenever you generate a visualization, add or update its entry here.

---

## 1. Model C7 chromosphere (single column)

### `model_c7_evolution.mp4`
9-panel time evolution of a Model C7 chromosphere run: a 2×4 block (neutral/ion density, velocity, pressure, temperature) plus a wide bottom panel with the number-density-weighted bulk velocity `(n_n u + n_i v)/(n_n + n_i)`. Input `model_c7_movie.txt` is a 1.5× run (`time_mult=1.5`, t≈829 s, 976 frames); rendered at 45 fps. **As of the "New explanation" hybrid upper BC** (docs/studies/evaporation/gentle_evaporation_downflow.md): `model_c7` now runs with `well_balanced` + a hydrostatic-pressure + EOS-density TR-base ghost (heat still via the imposed RTV flux q(T), not a Dirichlet T-jump). This removes the *spurious* part of the persistent TR downflow (interior V min ≈ −2 km/s → ≈ −1 km/s) and shifts the top from a ~200 kK corona to a ~10–20 kK TR base; the movie shows the rest → relaxation-upflow → mild-residual-downflow transient.
```bash
# regenerate the 1.5x run dump first:
build/chromo_main outputs/model_c7/model_c7_movie.txt full ionization model_c7 "" 1.5
.venv/bin/python util/animate_output.py outputs/model_c7/model_c7_movie.txt visualization/model_c7/model_c7_evolution.mp4 45
```
Positional: `[input.txt] [out.mp4] [fps]` — defaults `build/output.txt`, `visualization/model_c7/model_c7_evolution.mp4`, `15`.

### `model_c7_explicit_evolution.mp4`
Same 9-panel animation, for the explicit-scheme (R_I = 0) run.
```bash
.venv/bin/python util/animate_output.py outputs/model_c7/output_explicit.txt visualization/model_c7/model_c7_explicit_evolution.mp4
```

### `model_c7_snapshot.png`
Single 8-panel snapshot of the last frame of a run.
```bash
.venv/bin/python util/plot_output.py outputs/model_c7/model_c7_movie.txt visualization/model_c7/model_c7_snapshot.png
```
Positional: `[input.txt] [out.png]` — defaults `build/output.txt`, `visualization/model_c7/model_c7_snapshot.png`.

### `model_c7_compare.mp4`
8-panel animation overlaying two runs: semi-implicit (full) vs explicit (R_I = 0).
```bash
.venv/bin/python util/animate_compare.py outputs/model_c7/output_full.txt outputs/model_c7/output_explicit.txt visualization/model_c7/model_c7_compare.mp4
```
Positional: `[full.txt] [explicit.txt] [out.mp4] [fps]` — defaults `build/output_full.txt`, `build/output_explicit.txt`, `visualization/model_c7/model_c7_compare.mp4`, `15`.

### `model_c7_cooling_compare.mp4`
8-panel animation overlaying optically-thick cooling ON (Stage R) vs OFF. Defaults already point at the real files, so no args needed.
```bash
.venv/bin/python util/animate_cooling_compare.py outputs/model_c7/out_c7_cool_on.txt outputs/model_c7/out_c7_cool_off.txt visualization/model_c7/model_c7_cooling_compare.mp4
```
Positional: `[cool_on.txt] [cool_off.txt] [out.mp4] [fps]` — defaults `outputs/model_c7/out_c7_cool_on.txt`, `outputs/model_c7/out_c7_cool_off.txt`, `visualization/model_c7/model_c7_cooling_compare.mp4`, `6`.

### `model_c7_ionization_compare.mp4`
8-panel animation overlaying ionization ON (Stage E) vs OFF.
```bash
.venv/bin/python util/animate_ionization_compare.py outputs/model_c7/out_ioniz_on.txt outputs/model_c7/out_ioniz_off.txt visualization/model_c7/model_c7_ionization_compare.mp4
```
Positional: `[ionz_on.txt] [ionz_off.txt] [out.mp4] [fps]` — script defaults point at repo root (`../out_ioniz_on.txt`); pass the real `outputs/` files explicitly as above. (README has a worked example of this script.)

---

## 2. Open-BC and total-density (ρ_tot) runs

### `c7_openbc_rho_total.png`, `c7_openbc_rho_total.mp4`
Single-run ρ_tot = ρ_i + ρ_n evolution (space-time map + snapshots + column mass) for the open-BC C7 run.
```bash
.venv/bin/python util/visualize_rho_total.py outputs/model_c7/out_c7_openbc.txt --png visualization/model_c7/c7_openbc_rho_total.png --mp4 visualization/model_c7/c7_openbc_rho_total.mp4
```
Positional `input` (default `out_pfss_10x.txt`); flags `--png`, `--mp4`, `--fps` (24). One invocation emits both files.

### `out_c7_10x_on_rho_total.png`, `out_c7_10x_on_rho_total.mp4`
ρ_tot evolution for the C7 10× run; output names derive from the input stem (`out_c7_10x_on` → `<stem>_rho_total.{png,mp4}`).
```bash
.venv/bin/python util/visualize_rho_total.py outputs/model_c7/out_c7_10x_on.txt
```
(default output)

### `out_pfss_10x_rho_total.png`, `out_pfss_10x_rho_total.mp4`
ρ_tot evolution for the PFSS 10× run — the script's default input/output.
```bash
.venv/bin/python util/visualize_rho_total.py outputs/pfss/out_pfss_10x.txt
```
(default output)

### `pfss_10x_rho_total_compare.png`, `pfss_10x_rho_total_compare.mp4`
ρ_tot ionization on/off comparison for the PFSS 10× runs — the compare script's default invocation.
```bash
.venv/bin/python util/visualize_rho_total_compare.py outputs/pfss/out_pfss_10x.txt outputs/pfss/out_pfss_10x_noioniz.txt
```
Positional `runA runB`; `--labelA`/`--labelB` (default "ionization ON"/"OFF"), `--out-png`/`--out-mp4` (default `visualization/pfss_10x_rho_total_compare.{png,mp4}`), `--fps` (24). (default output)

### `c7_openbc_ioniz_compare.png`, `c7_openbc_ioniz_compare.mp4`
ρ_tot ionization on/off comparison for the open-BC C7 case.
```bash
.venv/bin/python util/visualize_rho_total_compare.py outputs/model_c7/out_c7_openbc.txt outputs/model_c7/out_c7_openbc_noioniz.txt --out-png visualization/model_c7/c7_openbc_ioniz_compare.png --out-mp4 visualization/model_c7/c7_openbc_ioniz_compare.mp4
```

### `pfss_openbc_ioniz_compare.png`, `pfss_openbc_ioniz_compare.mp4`
ρ_tot ionization on/off comparison for the open-BC PFSS case.
```bash
.venv/bin/python util/visualize_rho_total_compare.py outputs/pfss/out_pfss_openbc.txt outputs/pfss/out_pfss_openbc_noioniz.txt --out-png visualization/pfss/pfss_openbc_ioniz_compare.png --out-mp4 visualization/pfss/pfss_openbc_ioniz_compare.mp4
```

### `c7_openbc_8panel_compare.mp4`, `pfss_openbc_8panel_compare.mp4`
8-panel run comparison (`animate_compare.py`) for the open-BC C7 / PFSS cases, with overridden output paths.
```bash
.venv/bin/python util/animate_compare.py outputs/model_c7/out_c7_openbc.txt outputs/model_c7/out_c7_oldlike_openbc.txt visualization/model_c7/c7_openbc_8panel_compare.mp4
.venv/bin/python util/animate_compare.py outputs/pfss/out_pfss_openbc.txt outputs/pfss/out_pfss_openbc_noioniz.txt visualization/pfss/pfss_openbc_8panel_compare.mp4
```

---

## 3. Flare and loop runs

### `flare_height_time.png`, `flare_compare_snapshot.png`, `flare_timeseries.png`
Explosive vs gentle chromospheric-evaporation diagnostics (height–time maps, peak-phase profile comparison, time series). All three from one invocation (fixed output names).
```bash
.venv/bin/python util/visualize_flare.py outputs/model_flare/flare_explosive.txt outputs/model_flare/flare_gentle.txt 25.0
```
Positional (defaults): `explosive.txt`=`outputs/model_flare/flare_explosive.txt`, `gentle.txt`=`outputs/model_flare/flare_gentle.txt`, `t_max_s`=`25.0`.

### `flare_evolution_frames.png`, `flare_height_time_500.png`
Refined 500-cell explosive-evaporation frame strip + continuous height–time map. Both from one invocation (fixed output names).
```bash
.venv/bin/python util/visualize_flare_frames.py outputs/model_flare/flare_explosive_500.txt outputs/model_flare/flare_gentle.txt
```

### `loop_flare_evolution.mp4` (half loop) / `loop_full_flare_evolution.mp4` (full loop)
Topology-aware field-aligned flare movie with split x-axis (auto-detected from the `.dat` `[META] topology=`: half → vs height, full → vs arc length). Beam overlay read from `FLARE_*` env vars. Frames rendered in parallel.
```bash
# Full loop (default output)
FLARE_DELTA=5 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 \
  .venv/bin/python util/animate_loop_flare.py scenarios/data/loop_full_fine.dat outputs/loops/loop_full_flare.txt visualization/loops/loop_full_flare_evolution.mp4
# Half (closed) loop
.venv/bin/python util/animate_loop_flare.py scenarios/data/loop_closed_fine.dat outputs/loops/loop_half_flare.txt visualization/loops/loop_flare_evolution.mp4
```
Positional: `[loop.dat] [flare.txt] [out.mp4] [t_max_s]` — defaults `scenarios/data/loop_full_fine.dat`, `outputs/loops/loop_full_flare.txt`, `visualization/loops/loop_full_flare_evolution.mp4`.

### `loop_flare_verify.png`
Full/closed-loop flare run vs control: T_i and V vs true height, coronal-mass loading, final velocity profile (output path hardcoded).
```bash
.venv/bin/python util/verify_loop_flare.py scenarios/data/loop_closed_fine.dat outputs/loops/loop_closed_quiet.txt outputs/loops/loop_closed_flare.txt
```
Positional (defaults): `<loop.dat>`=`scenarios/data/loop_closed_test.dat`, `<control.txt>`=`outputs/loops/loop_closed_quiet.txt`, `<flare.txt>`=`outputs/loops/loop_closed_flare.txt`.

### `loop_cooling_full.png`, `loop_cooling_half.png`
Post-beam cooldown analysis (apex/EM-weighted coronal T & n vs time, e-folding cooling times, apex (n,T) trajectory). Default output `visualization/loops/loop_cooling.png`; pass an explicit name for the full/half variants.
```bash
.venv/bin/python util/loop_cooling.py scenarios/data/loop_full_fine.dat   outputs/loops/loop_full_flare.txt visualization/loops/loop_cooling_full.png
.venv/bin/python util/loop_cooling.py scenarios/data/loop_closed_fine.dat outputs/loops/loop_half_flare.txt visualization/loops/loop_cooling_half.png
```

### `loop_bc_compare.png`
Half-loop (reflecting apex) vs full-loop (interior apex) apex-BC comparison for one-sided beam heating. All paths hardcoded — no args.
```bash
.venv/bin/python util/compare_loop_bc.py
```

### `loop_ic_check.png`
Sanity-check a loop `.dat` initial condition (T, densities, ionization fraction, hydrostatic residual vs height).
```bash
.venv/bin/python util/check_loop_ic.py scenarios/data/loop_full_fine.dat visualization/loops/loop_ic_check.png
```
Positional (defaults): `<loop.dat>`=`scenarios/data/loop_closed_test.dat`, `[out.png]`=`visualization/loops/loop_ic_check.png`.

### `beam_stopping_feedback.png`
Thick-target beam stopping-depth feedback (loop becomes opaque as it fills): deposition-weighted mean height + chromospheric energy fraction over time.
```bash
.venv/bin/python util/beam_stopping_feedback.py scenarios/data/loop_full_fine.dat outputs/loops/loop_full_thick.txt visualization/loops/beam_stopping_feedback.png 20 5
```
Positional (defaults): `[loop.dat]`, `[thick.txt]`, `[out.png]`, `[E_cut_keV]`=`20`, `[delta]`=`5`, `[T_off_s]`=`62`.

---

## 4. PFSS and 3-D event renders

### `pfss_3d_global.png`, `pfss_3d_local.png`, `pfss_3d_rotate.mp4`
Global PFSS sphere (B_LOS + altitude-coloured field lines), a local footpoint box, and a rotating-camera movie — from one PFSS solve. Output names are hardcoded; `--out-dir` sets location. `--no-movie` skips the MP4.
```bash
.venv/bin/python util/visualize_pfss_3d.py \
    --magnetogram util/data/mrzqs190801t0014c2220_229.fits \
    --nseeds 96 --nrho 60 --rss 2.5 \
    --out-dir visualization/pfss
```

### `pfss_lines_evolution.png`, `pfss_lines_evolution.mp4`
2–3 PFSS field lines side-by-side: 3D geometry + BC table (PNG) and a per-line 4-field strip animation (MP4). `len(--dats)` must equal `len(--outs)`.
```bash
.venv/bin/python util/visualize_pfss_lines_evolution.py \
    --magnetogram util/data/mrzqs190801t0014c2220_229.fits \
    --dats scenarios/data/pfss_lineA.dat scenarios/data/pfss_lineB.dat scenarios/data/pfss_lineC.dat \
    --outs outputs/pfss/out_lineA.txt outputs/pfss/out_lineB.txt outputs/pfss/out_lineC.txt \
    --out-png visualization/pfss/pfss_lines_evolution.png \
    --out-mp4 visualization/pfss/pfss_lines_evolution.mp4
```

### `pfss_lineA_flare_evolution.mp4`
4-panel vertical-profile animation (T_i, V, densities, ionization) of a flaring PFSS line (PFSS_FLARE run), beam window shaded.
```bash
.venv/bin/python util/animate_pfss_flare.py outputs/pfss/pfss_lineA_flare.txt visualization/pfss/pfss_lineA_flare_evolution.mp4
```
Positional: `[input.txt] [out.mp4]` — defaults `outputs/pfss/pfss_lineA_flare.txt`, `visualization/pfss/pfss_lineA_flare_evolution.mp4`.

### `pfss_10x_evolution.mp4`
Same 4-panel flare-profile animation applied to the 10×-duration PFSS relaxation run.
```bash
.venv/bin/python util/animate_pfss_flare.py outputs/pfss/out_pfss_10x.txt visualization/pfss/pfss_10x_evolution.mp4
```

### `pfss_lineA_flare_verify.png`
Phase-3 de-risk: control (beam off) vs explosive (beam on) PFSS-line run; height-time T_i and V heatmaps + printed diagnostics.
```bash
.venv/bin/python util/verify_pfss_flare.py outputs/pfss/pfss_lineA_control.txt outputs/pfss/pfss_lineA_flare.txt visualization/pfss/pfss_lineA_flare_verify.png
```

### `canopy_geometry_evolution.png`, `canopy_evolution.mp4`
`analytic_canopy` scenario: static figure (3D "trumpet" flux tube + B(z)/A(z) + BC table + 4-panel evolution) and the matching evolution movie.
```bash
.venv/bin/python util/visualize_canopy.py \
    --sim-output outputs/analytic_canopy/out_canopy.txt \
    --out-png visualization/analytic_canopy/canopy_geometry_evolution.png \
    --out-mp4 visualization/analytic_canopy/canopy_evolution.mp4
```

### `event_ensemble_3d_peak.png`, `event_ensemble_3d.mp4`
Simulated event ensemble on a 3D solar sphere — each PFSS line drawn along its real path, coloured by its own run's temperature. PNG = snapshot at `--snap-time`; MP4 = flare over time. Names hardcoded; `--out-dir` sets location. `--no-anim` skips the MP4.
```bash
.venv/bin/python util/visualize_event_3d.py \
    --lines-json outputs/events/event_ensemble_lines.json \
    --magnetogram util/data/adapt40311_044012_202408010600_i00053600n1.fts.gz \
    --out-dir visualization/events
```

### `event_ensemble_3d_global.png`
Wide framing of the same 3D ensemble — same script; only `event_ensemble_3d_peak.png` is emitted natively, so this is a second run (wider `--zoom-rsun`) whose PNG was renamed.
```bash
.venv/bin/python util/visualize_event_3d.py \
    --lines-json outputs/events/event_ensemble_lines.json \
    --magnetogram util/data/adapt40311_044012_202408010600_i00053600n1.fts.gz \
    --out-dir visualization/events --no-anim
mv visualization/events/event_ensemble_3d_peak.png visualization/events/event_ensemble_3d_global.png
```

### `event_20240801_AR13768_flare.mp4`, `event_20240801_AR13768_Te.mp4`
Event-starter full-loop flare animation (AR 13768): split-x 4-panel profile movie vs arc length, live thick-target beam Q(s) overlay. `_Te` variant uses the 3-temperature (`ENABLE_TE=1`) dump (auto-detected).
```bash
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 \
  .venv/bin/python util/animate_loop_flare.py scenarios/data/event_20240801_AR13768.dat outputs/events/event_20240801_AR13768_flare.txt visualization/events/event_20240801_AR13768_flare.mp4
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 \
  .venv/bin/python util/animate_loop_flare.py scenarios/data/event_20240801_AR13768.dat outputs/events/event_te/event_20240801_AR13768_te.txt visualization/events/event_20240801_AR13768_Te.mp4
```

### `event_20240801_c05_101Mm_Te.mp4`, `event_20240801_o00_1040Mm_Te.mp4`
3-temperature ensemble-member animations (closed 101 Mm; open 1040 Mm).
```bash
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 \
  .venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/events/event_te/c05_101Mm_te.txt  visualization/events/event_20240801_c05_101Mm_Te.mp4
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 \
  .venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/events/event_te/o00_1040Mm_te.txt visualization/events/event_20240801_o00_1040Mm_Te.mp4
```

### `event_ensemble_c05_101Mm_closed.mp4`, `event_ensemble_o00_1040Mm_open.mp4`
Ensemble-member flare animations: closed loop (vs arc length) and open line (vs height).
```bash
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 \
  .venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/events/event_ensemble/c05_101Mm.txt  visualization/events/event_ensemble_c05_101Mm_closed.mp4
FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 \
  .venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/events/event_ensemble/o00_1040Mm.txt visualization/events/event_ensemble_o00_1040Mm_open.mp4
```

---

## 5. Synthetic Doppler and two-fluid diagnostics

**Naming rule.** `synth_doppler.py` and `compare_doppler_obs.py` take `<loop.dat> <state.txt> [out_stem]` and derive `<stem>` from `out_stem` (or the `.txt` basename). They then write fixed-suffix files into `visualization/`:
- `synth_doppler.py` → `<stem>_doppler.png`, `<stem>_dopplergram.png` (+ `outputs/<stem>_doppler.csv`).
- `compare_doppler_obs.py` → `<stem>_vs_obs.png`.

The shaded beam window is read from `FLARE_T_ON`/`FLARE_DUR` (defaults 2.0 s / 10.0 s); `_dur6` variants pass `FLARE_DUR=6`.

### `c05_101Mm_te_{doppler,dopplergram}.png`
Synthetic IRIS Doppler (closed 101 Mm loop, T_e≠T_i run).
```bash
.venv/bin/python util/synth_doppler.py scenarios/data/event_ensemble/c05_101Mm.dat outputs/events/event_te/c05_101Mm_te.txt c05_101Mm_te
```

### `o00_1040Mm_te_{doppler,dopplergram}.png`
Synthetic IRIS Doppler (open 1040 Mm line, T_e≠T_i run).
```bash
.venv/bin/python util/synth_doppler.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/events/event_te/o00_1040Mm_te.txt o00_1040Mm_te
```

### `event_20240801_AR13768_te_{doppler,dopplergram}.png`
Synthetic IRIS Doppler (2024-08-01 AR13768 line, T_e≠T_i run).
```bash
.venv/bin/python util/synth_doppler.py scenarios/data/event_20240801_AR13768.dat outputs/events/event_te/event_20240801_AR13768_te.txt event_20240801_AR13768_te
```

### `sol20140329_10Mm_1T_{doppler,dopplergram}.png`
SOL2014-03-29 10 Mm loop, single-temperature (T_e=T_i) run.
```bash
.venv/bin/python util/synth_doppler.py scenarios/data/sol20140329_10Mm.dat outputs/events/sol2014/sol20140329_10Mm_1T.txt sol20140329_10Mm_1T
```

### `sol20140329_10Mm_te_{doppler,dopplergram}.png`
SOL2014-03-29 10 Mm loop, T_e≠T_i run.
```bash
.venv/bin/python util/synth_doppler.py scenarios/data/sol20140329_10Mm.dat outputs/events/sol2014/sol20140329_10Mm_te.txt sol20140329_10Mm_te
```

### `sol20140329_10Mm_te_stab_{doppler,dopplergram}.png`
SOL2014-03-29 10 Mm loop, T_e≠T_i, numerically-stabilized variant.
```bash
.venv/bin/python util/synth_doppler.py scenarios/data/sol20140329_10Mm.dat outputs/events/sol2014/sol20140329_10Mm_te_stab.txt sol20140329_10Mm_te_stab
```

### `sol20140329_dur6_te_{doppler,dopplergram}.png`
SOL2014-03-29 10 Mm loop, T_e≠T_i, 6 s beam-duration variant.
```bash
FLARE_DUR=6 .venv/bin/python util/synth_doppler.py scenarios/data/sol20140329_10Mm.dat outputs/events/sol2014/sol20140329_dur6_te.txt sol20140329_dur6_te
```

### `sol20140329_10Mm_te_vs_obs.png`, `sol20140329_dur6_te_vs_obs.png`
Synthetic loop Doppler (Fe XXI / Si IV / C II / Mg II) overlaid on observed IRIS values for SOL2014-03-29 (Young+2015, Li+2015).
```bash
.venv/bin/python util/compare_doppler_obs.py scenarios/data/sol20140329_10Mm.dat outputs/events/sol2014/sol20140329_10Mm_te.txt sol20140329_10Mm_te
FLARE_DUR=6 .venv/bin/python util/compare_doppler_obs.py scenarios/data/sol20140329_10Mm.dat outputs/events/sol2014/sol20140329_dur6_te.txt sol20140329_dur6_te
```

### `two_fluid_diag_{c05_closed,o00_open,starter}.png`
Two-fluid decoupling diagnostic (slip-Mach |V−U|/c_s and |T_i−T_n|/T_i heatmaps). Args `<loop.dat> <state.txt> [tag]`; `tag` → `two_fluid_diag_<tag>.png`.
```bash
.venv/bin/python util/diag_two_fluid.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/two_fluid/c05_closed.txt c05_closed
.venv/bin/python util/diag_two_fluid.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/two_fluid/o00_open.txt   o00_open
.venv/bin/python util/diag_two_fluid.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/two_fluid/starter.txt   starter
```

### `two_fluid/compare_{c05_closed,o00_open}.png`
Neutrals-ON (two-fluid) vs neutrals-OFF (`SINGLE_FLUID=1`) comparison: ΔT(s,t), ΔV(s,t) heatmaps + peak-departure profiles. Args `<loop.dat> <on.txt> <off.txt> [tag]`; `tag` → `two_fluid/compare_<tag>.png`.
```bash
.venv/bin/python util/compare_two_fluid.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/two_fluid/c05_neutrals_on.txt outputs/two_fluid/c05_neutrals_off.txt c05_closed
.venv/bin/python util/compare_two_fluid.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/two_fluid/o00_neutrals_on.txt outputs/two_fluid/o00_neutrals_off.txt o00_open
```

### `two_fluid/{c05,o00}_neutrals_{on,off}.mp4`
Per-run movies of the neutrals-on / neutrals-off two-fluid runs (the same `outputs/two_fluid/` dumps used by the comparison PNGs above), animated with `animate_loop_flare.py` and written under `visualization/two_fluid/`.
```bash
.venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/two_fluid/c05_neutrals_on.txt  visualization/two_fluid/c05_neutrals_on.mp4
.venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/two_fluid/c05_neutrals_off.txt visualization/two_fluid/c05_neutrals_off.mp4
.venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/two_fluid/o00_neutrals_on.txt  visualization/two_fluid/o00_neutrals_on.mp4
.venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/two_fluid/o00_neutrals_off.txt visualization/two_fluid/o00_neutrals_off.mp4
```

---

## 6. Physics reference plots (scripts live in `visualization/reference/`)

### `c7_ioniz_recomb_rates.png`
Volumetric H ionization (photo / multilevel-CR / Voronov) and recombination (radiative / three-body / molecular) rates + Route-B photoionization closure, vs height for Model C7. No args; writes next to itself.
```bash
.venv/bin/python visualization/reference/plot_c7_ioniz_recomb_rates.py
```

### `heat-conductivies-c7-wide-range.png`, `heat-conductivities-c7-narrow-range.png`
Electron/neutral heat conductivities along the C7 column: full-range log panel (filename intentionally spelled `conductivies`) and the 2138–2153 km linear zoom. No flags → writes both; `--wide`/`--narrow` write one; `-o/--outfile` overrides (only with `--wide`/`--narrow`).
```bash
.venv/bin/python visualization/reference/plot_heat_conductivity_c7.py
```
> Note: a correctly-spelled `heat-conductivities-c7-wide-range.png` also exists in `visualization/reference/` but is **not** emitted by this script's default paths (the wide default is the misspelled `conductivies`). It is a renamed copy — `mv` the wide output if you want the corrected name.

---

## 7. Prerequisite input generators (write to `outputs/` & `scenarios/data/`, not `visualization/`)

### `util/make_loop_dat.py`
Builds a chromosphere→TR→corona loop `.dat` (C7 chromosphere + hydrostatic corona) — the `.dat` consumed by all loop plots.
```bash
.venv/bin/python util/make_loop_dat.py --output scenarios/data/loop_full_fine.dat --topology full --half-length-mm 25 --apex-T-MK 2.0 --fp-B-G 100
```

### `util/extract_field_line.py`
Extracts a single PFSS field line from a GONG synoptic magnetogram and writes the `pfss_field_line` `.dat`.
```bash
.venv/bin/python util/extract_field_line.py --magnetogram util/data/mrzqs190801t0014c2220_229.fits --output scenarios/data/pfss_qs_20190801.dat --ns 100 --top-height-km 986
```

### `util/extract_event_loop.py`
Extracts one flare-loop geometry (half-length, topology, footpoint |B|) from a real-event magnetogram via PFSS; emits a geometry JSON and prints the next `make_loop_dat.py` command.
```bash
.venv/bin/python util/extract_event_loop.py --magnetogram util/data/adapt40311_044012_202408010600_i00005600n1.fts.gz --ar-lon 74 --ar-lat -16 --obstime 2024-08-01T07:09:00 --out-json outputs/events/event_20240801_AR13768_geometry.json
```

### `util/extract_event_ensemble.py`
Traces an ensemble of PFSS lines around an AR, builds a loop `.dat` per line, and emits the `run_ensemble.py` manifest + the 3D lines index (`event_ensemble_lines.json`, consumed by `visualize_event_3d.py`).
```bash
.venv/bin/python util/extract_event_ensemble.py --magnetogram util/data/adapt40311_044012_202408010600_i00053600n1.fts.gz --ar-lon 74 --ar-lat -16 --obstime 2024-08-01T07:09:00 --n-closed 10 --n-open 4
```

### `util/run_ensemble.py`
Embarrassingly-parallel ensemble runner: one `chromo_main` per field line, writing each `<out>.txt` (+ `.run.log`).
```bash
# Mode A (manifest):
.venv/bin/python util/run_ensemble.py --manifest outputs/events/event_ensemble_jobs.json --workers 16
# Mode B (glob, shared env):
.venv/bin/python util/run_ensemble.py --dats "scenarios/data/event_ensemble/*.dat" --out-dir outputs/events/event_ensemble --env PFSS_FLARE=1 FLARE_BEAM_FLUX=3e7 FLARE_DELTA=4 FLARE_E_CUT=20
```

---

## 8. Presentation movie posters (linked, not embedded)

The "Evolution movies" slides in `docs/presentations/2026-07-11-group-meeting/presentation.tex` show a **first-frame poster** that is a
clickable link to the movie file: `\href{run:<movie>.mp4}{\includegraphics{poster}}`.
Clicking fires a PDF `/Launch` action → Adobe Acrobat/Reader opens the `.mp4` in the system default
video player (a one-time "open this file?" confirmation appears). **macOS Preview cannot** — it
sandbox-blocks launch actions; present in Acrobat Reader. The movies are **linked, not embedded**,
so the PDF stays ~1.2 MB.

**Same-folder requirement (Acrobat security).** Adobe restricts `/Launch` file paths to the PDF's own
directory and rejects `../` traversal — so the link is a **bare filename** and the two `.mp4`s must sit
**next to the PDF**. After building, copy them in (re-copy after any `latexmk -C`); to present
elsewhere, move `presentation.pdf` + both `.mp4`s together into one folder:
```bash
cp visualization/events/event_20240801_c05_101Mm_Te.mp4 visualization/events/event_20240801_o00_1040Mm_Te.mp4 docs/presentation-build/
```

The linked movies are regenerated *near-lossless* (the stock encode is ~400 kbps and visibly soft);
`animate_loop_flare.py` honours `ANIM_CRF` / `ANIM_PIX_FMT` / `ANIM_PRESET` (encode, see `util/_anim_parallel.py`),
`ANIM_FPS` / `ANIM_MAX_FRAMES` (playback rate / frame budget), `ANIM_UNIFORM_TIME=1`, and the per-phase
speed-ups `ANIM_SPEED_PRE` / `ANIM_SPEED_POST`.

**Current runs span 0–100 s, 60 fps, uniform-time, with a 2× faster cooling tail.** Simulations rerun with
`time_mult=2.0` (manifest `outputs/events/event_te_2x_jobs.json`); NB `time_mult` scales the *step cap* and dt grows
~10× post-flare, so the runs overshoot (c05→322 s, o00→371 s of physics). The movie is capped at **100 s** via
the 4th positional arg (`t_cap`) — both runs have data well past that. **`ANIM_UNIFORM_TIME=1` is essential
here:** the dump writes a frame every N steps, so with the post-flare dt blow-up the frames pile up in the
first ~30 s and the long cooling tail flashes past (the clip *looks* like it stops ~50 s). Uniform-time
resampling lays a uniform-in-sim-time target grid and snaps each point to the nearest dump frame, so the whole
span plays at a steady rate. `ANIM_SPEED_PRE=1 ANIM_SPEED_POST=2` then allocate target frames per phase
∝ Δt/speed, so the pre-beam + beam-on window (t ≤ T_OFF = 12 s) plays at 1× and the post-beam cooling tail
plays **2× faster** (215 frames pre/beam, 785 post-beam). `ANIM_MAX_FRAMES=1000` (the movie is *linked*, so
frame count doesn't grow the PDF); both clips run 16.7 s.

### `poster_c05_closed.jpg`, `poster_o00_open.jpg`  (+ regenerated `event_20240801_{c05_101Mm,o00_1040Mm}_Te.mp4`)
```bash
# (0) 2x-length simulations (overshoot, capped at render time)
.venv/bin/python util/run_ensemble.py --manifest outputs/events/event_te_2x_jobs.json --workers 2
# (1) render near-lossless, 0–100 s (t_cap), 60 fps, uniform-time, post-beam 2x (env inline — zsh won't word-split a var)
ENV="FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 ANIM_UNIFORM_TIME=1 ANIM_SPEED_PRE=1 ANIM_SPEED_POST=2 ANIM_FPS=60 ANIM_MAX_FRAMES=1000 ANIM_CRF=12 ANIM_PIX_FMT=yuv444p ANIM_PRESET=slow"
env ${=ENV} .venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/events/event_te/c05_101Mm_te.txt  visualization/events/event_20240801_c05_101Mm_Te.mp4 100
env ${=ENV} .venv/bin/python util/animate_loop_flare.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/events/event_te/o00_1040Mm_te.txt visualization/events/event_20240801_o00_1040Mm_Te.mp4 100
# (2) extract the first-frame poster shown on each slide
ffmpeg -y -i visualization/events/event_20240801_c05_101Mm_Te.mp4  -frames:v 1 -q:v 2 visualization/events/poster_c05_closed.jpg
ffmpeg -y -i visualization/events/event_20240801_o00_1040Mm_Te.mp4 -frames:v 1 -q:v 2 visualization/events/poster_o00_open.jpg
```

---

## 9. Gentle conduction-driven evaporation (docs/studies/evaporation/gentle_evaporation_plan.md)

> **Reconciliation note (2026-07-11, docs/design/scenario_reconciliation_plan.md):** the standalone
> `model_gentle` scenario was merged into **`model_column`**; `model_gentle` is now a backward-compat
> alias that applies the documented stable resolved-corona full-physics preset
> (`ISO_CORONA ISO_H_BASE=1003 ISO_HEAT_FLUX ISO_COOLING ISO_IONIZATION ISO_TWO_FLUID`). The former
> `GENTLE_*` knobs are renamed: `GENTLE_Q_*` → `ISO_QFLUX*`, `GENTLE_HEAT_*` → `ISO_CHEAT*`,
> `GENTLE_T_BOOST*` → `ISO_TBOOST*`. The `§9`/`pfss_field_line GENTLE` overlay commands below (which
> drive gentle evaporation on a closed loop) are unaffected — that overlay lives in pfss_field_line.cpp.

The `GENTLE=1` overlay (scenarios/pfss_field_line.cpp) adds the ambient coronal heating H(s)=E_H0·exp(−d/s_H) that makes a resolved chromosphere→TR→corona loop a true steady state, then a ×3 heating ramp drives gentle (subsonic, conduction-driven) evaporation. Run from the project root with `.venv` active.

**Improvement (2026-06-25): the well-balanced reconstruction (`Grid::well_balanced`, docs/studies/evaporation/gentle_evaporation_downflow.md) is now ON in the GENTLE overlay** (default on; `GENTLE_WELL_BALANCED=0` reverts). It undoes the φ_g MUSCL-slope inconsistency that gave every hydrostatic atmosphere a spurious (γ−1)g downforce, so the relaxed preflare baseline sits closer to V≈0 (ρ-weighted mean|V| 112 → **68 m/s**) and the gentle upflow is measured against a quieter zero. The 5958-test baseline is byte-for-byte unchanged (GENTLE is opt-in). The commands below are unchanged; the numbers were re-run with the fix on.

```bash
# (0) build the heating-consistent loop IC (prints the auto-calibrated E_H0, s_H)
python util/make_loop_dat.py --output scenarios/data/loop_closed_gentle.dat \
    --topology closed --half-length-mm 25 --apex-T-MK 2.0 --fp-B-G 100

# (1) Phase-2 relaxation to the steady preflare state (steady heating, enhance=1)
GENTLE=1 ./build/chromo_main outputs/model_column/gentle_relax.txt full ionization \
    pfss_field_line scenarios/data/loop_closed_gentle.dat 12 cooling
# → visualization/model_column/gentle_relax_convergence.png : max|V|(t), mean|V|(t), max|dlnT/dt|, T(t)
python util/check_relaxation.py outputs/model_column/gentle_relax.txt visualization/model_column/gentle_relax_convergence.png

# (2) Phase-3 combined run: relax (to t=1000 s) then ramp heating ×3 over 200 s
GENTLE=1 GENTLE_ENHANCE=3 GENTLE_T_ON=1000 GENTLE_RAMP=200 ./build/chromo_main \
    outputs/model_column/gentle_evap.txt full ionization pfss_field_line \
    scenarios/data/loop_closed_gentle.dat 16 cooling
# → visualization/model_column/gentle_evap_signature.png : coronal EM rise, TR upflow, T step, upflow Mach (≪1)
python util/gentle_evaporation_diag.py outputs/model_column/gentle_evap.txt visualization/model_column/gentle_evap_signature.png 1000 200
# → visualization/model_column/gentle_evap_profiles.png : split-axis preflare-vs-evaporation T, n_e, V profiles
python util/gentle_evap_profiles.py scenarios/data/loop_closed_gentle.dat \
    outputs/model_column/gentle_evap.txt visualization/model_column/gentle_evap_profiles.png 1000

# (3) evolution movies (split-axis T/V/n/ionization; parallel frames → MP4).
#     animate_gentle.py reuses animate_loop_flare's transform but relabels the
#     timeline RELAXING / HEATING RAMP / EVAPORATION (no beam) and overlays H(s,t).
#     Trailing args = t_on ramp enhance (match the run; t_on=1e9 ⇒ always RELAXING).
python util/animate_gentle.py scenarios/data/loop_closed_gentle.dat \
    outputs/model_column/gentle_relax.txt visualization/model_column/gentle_relax_evolution.mp4 1e9 1 1
python util/animate_gentle.py scenarios/data/loop_closed_gentle.dat \
    outputs/model_column/gentle_evap.txt  visualization/model_column/gentle_evap_evolution.mp4 1000 200 3
```

- `gentle_relax_convergence.png` — Phase-2 steady-state convergence monitor (well_balanced ON: ρ-weighted mean|V| → **68 m/s**, apex T → 2.24 MK, coronal⟨T⟩ → 1.41 MK; max|dlnT/dt| oscillates around the 3e-3 target — bounded TR-cell spikes, not a growing drift; the residual few-km/s coronal downflow is the physical energy-budget flow in low-density gas).
- `gentle_evap_signature.png` — Phase-3 gentle-evaporation fingerprint (well_balanced ON): coronal EM **×2.53**, apex T 2.38→3.12 MK, peak upflow **14.38 km/s at Mach 0.07** (subsonic ⇒ gentle, Antiochos & Sturrock 1978).
- `gentle_evap_profiles.png` — field-aligned T/n_e/V on the split x-axis, preflare (t=1000 s, apex 2.19 MK, max +V 2.9 km/s) vs evaporation (t=1200 s, apex 2.56 MK, max +V 14.8 km/s): TR burns down, corona fills, sustained subsonic upflow.
- `gentle_relax_evolution.mp4` — relaxation evolution movie (steady heating): the IC-settling transient decays toward the V≈0 steady preflare state.
- `gentle_evap_evolution.mp4` — full evolution movie: relaxation → ×3 heating ramp → gentle conduction-driven evaporation (subsonic upflow fills the corona, TR burns down), phase-labeled with the heating overlay H(s,t).

### §9b. Gentle evaporation v2 — corona-as-boundary at 2.628 Mm (NEGATIVE result; see §9c for the fix)

First corona-as-boundary attempt: single C7 column lifted only to h ≈ 2.628 Mm with real Avrett & Loeser Table 26 data, corona imposed as the `q(T)` boundary, driven by **ramping q(T)**. Kept as a documented negative result — the boundary sat at the top of the TR with no coronal volume, so no evaporation appeared (§9c lifts it into a real corona and fixes this).

```bash
# Phase B — relaxation to the steady preflare state (steady q, enhance=1 default)
./build/chromo_main outputs/model_column/gentle_v2_relax.txt full ionization model_gentle - 6 cooling
python util/check_relaxation.py outputs/model_column/gentle_v2_relax.txt visualization/model_column/gentle_v2_relax_convergence.png

# Phase B+C combined — relax (to t=1000 s) then ramp q(T) ×3 over 120 s
GENTLE_Q_ENHANCE=3 GENTLE_Q_T_ON=1000 GENTLE_Q_RAMP=120 \
    ./build/chromo_main outputs/model_column/gentle_v2_evap.txt full ionization model_gentle - 8 cooling
# evolution movie (6 height panels: T / total gas pressure / V / n_e,n_n /
#   ionization / field-aligned conductive flux q∥; phase-labeled by the q(T) ramp)
python util/animate_gentle_column.py outputs/model_column/gentle_v2_evap.txt \
    visualization/model_column/gentle_v2_evolution.mp4 1000 120 3
```

- `gentle_v2_relax_convergence.png` — Phase-B convergence: ρ-weighted mean|V| → 48 m/s, max|dlnT/dt| → 3e-5 /s, top settles at 0.45 MK (Mach-capped outflow; a free outflow collapses the top).
- `gentle_v2_evolution.mp4` — relaxation → q(T)-ramp ×3. The conductive-flux panel shows ~−500 W/m² being driven *down* into the TR, yet the velocity panel stays in ~−4 km/s downflow with **no upflow** and the column-integrated density is flat — i.e. **no evaporation**. The corona-as-boundary at 2.628 Mm sits right at the top of the TR, leaving no coronal volume for evaporated mass to fill (the EM-rise that *is* evaporation, Antiochos & Sturrock 1978). The pressure panel shows the gas pressure building toward a hotter static balance with no lower-pressure corona above to expand into. To see evaporation the corona must be a resolved *volume* (§9c/v3 below, or §9a/v1).

---

### §9c. Gentle evaporation v3 — boundary lifted HIGH into a resolved corona (model_gentle)

The fix for v2's null result (docs/studies/evaporation/gentle_evaporation_plan.md, option A): the extended C7 table (Avrett & Loeser 2008 Table 26) now reaches 68 Mm / 1.59 MK, and `model_gentle` lifts the upper `q(T)` boundary to `GENTLE_TOP_KM` (default 10 Mm, T ≈ 1.1 MK after relaxation) — a genuine resolved coronal **volume** for the evaporated mass to fill. `GENTLE_NS ≈ 300` keeps both the TR (TRAC-broadened) and the corona resolved.

```bash
# Phase B — relaxation: does the resolved C7 corona hold under a q(T) top BC?
GENTLE_NS=300 ./build/chromo_main outputs/model_column/gentle_v3_relax.txt full ionization model_gentle - 6 cooling

# Phase B+C combined — relax (to t=3000 s) then ramp q(T) ×3 over 300 s
GENTLE_NS=300 GENTLE_Q_ENHANCE=3 GENTLE_Q_T_ON=3000 GENTLE_Q_RAMP=300 \
    ./build/chromo_main outputs/model_column/gentle_v3_evap.txt full ionization model_gentle - 12 cooling
# evolution movie — Nature-styled, split x-axis (zoom TR / compress corona),
#   slow-motion in the evaporation-onset window. Default slow window = [t_on, t_on+1000].
python util/animate_gentle_column.py outputs/model_column/gentle_v3_evap.txt \
    visualization/model_column/gentle_v3_evolution.mp4 3000 300 3
# tune slow-mo: SLOW_T0 / SLOW_T1 [s] window, SLOW_FACTOR (default 4), ANIM_MAX_FRAMES (500)

# the imposed conductive-flux ramp q(t) itself (the gentle-evaporation driver)
python util/plot_q_ramp.py        # defaults reproduce gentle_v3: q0=230, x3, t_on=3000, ramp=300

# writeup summary figure (driver+EM response, and corona-filling profiles).
# NOTE: writes to docs/writeup-overleaf/ (a committed paper figure), not visualization/
python util/plot_gentle_summary.py
```

- `gentle_v3_evolution.mp4` — relaxation → q(T)-ramp ×3 in the resolved corona. **Evaporation now appears**: across the ramp the coronal (h > 2.3 Mm) emission measure ∫n_e² ds **× 2.0**, coronal mass ∫n_e ds ×1.42, T_top 1.14 → 1.50 MK, n_e(top) ×1.47, TR burns down 1.595 → 1.526 Mm, with a transient +1 km/s upflow during the ramp. Contrast §9b (same driver, boundary at 2.628 Mm): EM flat, no evaporation. **Confirms the v2 null was geometric** — the boundary needed coronal volume above the TR, not a different mechanism. Caveat: the relaxed baseline carries a steady ~8 km/s cap-throttled coronal downflow (no volumetric heating to pin V=0), so the *sustained* signature is the EM/density rise rather than a fast bulk upflow; the resolved-corona-with-heating run (§9a/v1) gives the cleaner upflow (15.6 km/s, EM ×2.6).
  - 7 panels (T / P / ion velocity V / bulk velocity (n_iV+n_nU)/(n_i+n_n) — both upflow-downflow shaded / n_e,n_n / ionization / conductive flux q∥) on a **split x-axis** (zoom sub-Mm chromosphere+TR, compress the corona; region bands chromo|TR|corona), with a live phase chip + time readout and a bottom timeline marking the **slow-motion ×4** window [3000,4000] s that stretches the evaporation onset. Nature-figure styling (restrained blue/red palette, despined axes, direct labels).
- `q_ramp.png` / `.pdf` — the imposed conductive-flux driver q(t) vs time: flat at q₀ ≈ 230 W m⁻² through relaxation, raised-cosine (Hann) ramp ×3 over [t_on, t_on+Δ] = [3000, 3300] s, then flat at 3q₀ ≈ 690 W m⁻². The ramp envelope is a code construct (smooth-start, not from a paper — see docs/studies/boundaries/q_of_T_boundary_condition.md); only q₀ and the slow-ramp requirement are physics.

### §9d. Heat-flux–driven (A&S78) gentle evaporation — hotter corona, NO ramp (`gentle_hotcorona_evolution.mp4`)

User-directed variant (2026-06-25): drive the evaporation by the **coronal conductive heat flux** (Antiochos & Sturrock 1978), exactly as the resolved-corona `iso_corona` did — **not** by ramping anything. `model_gentle` is already the right scenario: it is two-fluid + ionization + radiative cooling, and its only driver is the conductive flux `q(T)` (it has **no** volumetric heating term — that `H(s)` ramp belongs to the `pfss_field_line` `GENTLE` overlay). The new `GENTLE_T_BOOST` knob applies a **one-time coronal superheat** to the IC (smooth tanh selector above ~3×10⁵ K, so only the corona/upper-TR is heated; chromosphere untouched; density held fixed) so the corona starts HOTTER than its own conduction+radiation balance (~1.14 MK). Its own Spitzer flux (κ_e ∝ T^{5/2}) then conducts down the TR and drives evaporation from t=0, with `q(T)` held **steady** (`GENTLE_Q_ENHANCE=1`). Default `GENTLE_T_BOOST=1` ⇒ model_gentle unchanged.

**Two boundary/interior fixes were applied to `model_gentle` and A/B-tested against the downflow — both default ON, both leave the coronal downflow intact, jointly proving it is physical:**

1. **`well_balanced`** (`GENTLE_WELL_BALANCED`, default on) — the φ_g interior-reconstruction fix, previously missing here. A/B (off→on): halves the *chromospheric* mean downflow (−0.29 → −0.14 km/s) and removes the upper-chromosphere/TR-base downflow (h≈1.8 Mm: −2.83 → +0.08 km/s) — that part was the reconstruction artifact. *Coronal* downflow unchanged (mean −6.13 → −6.16 km/s).
2. **Hydrostatic outer BC** (`GENTLE_HYDROSTATIC_BC`, default on) — `model_gentle` now calls `model_c7_ic(extended=true, tr_jump_bc=true)`, so the coronal-top ghost imposes the iso-style **hydrostatic pressure** (p_{ns-2}−p_ghost)/2Δs = ρ_{ns-1} g + EOS density + continuous ghost T (heat still via the imposed q(T)), instead of the old Neumann (zero-gradient) RTV reservoir whose dp/ds≈0 left gravity unbalanced at the top. This makes V=0 a boundary fixed point — the **resolved corona is kept** (so evaporation is still possible). A/B (RTV→hydrostatic): coronal mean downflow −6.16 → **−6.13 km/s**, top cell −7.81 → −7.72, identical to ~0.3 km/s at every coronal height. So the boundary closure is **not** the cause either.

**Conclusion (three independent levers — well_balanced, hydrostatic BC, and the heat-flux driver — all leave it):** the big ~4–8 km/s coronal downflow is the **physical energy-budget downflow** — an under-heated corona (boundary q(T), no distributed heating) sheds its leftover radiative loss by advecting enthalpy down to the TR, the flow driven from the *interior* and exiting through the Mach-capped outflow (saturating ~0.05·c_s ≈ 8 km/s at the tenuous top). It is neither a reconstruction artifact nor a boundary artifact; pinning coronal V≈0 needs **distributed volumetric heating H(s)** (RTV static balance), which a boundary flux cannot provide. (This default-BC change also affects §9c's `model_gentle` runs; `GENTLE_HYDROSTATIC_BC=0` recovers the old RTV-reservoir behavior.)

```bash
GENTLE_NS=300 GENTLE_T_BOOST=1.5 GENTLE_Q_ENHANCE=1 \
    ./build/chromo_main outputs/model_column/gentle_hotcorona.txt full ionization model_gentle - 6 cooling
# steady run ⇒ animate with enhance=1 (phase chip reads RELAXING throughout); slow-mo the transient
SLOW_T0=0 SLOW_T1=200 SLOW_FACTOR=4 .venv/bin/python util/animate_gentle_column.py \
    outputs/model_column/gentle_hotcorona.txt visualization/model_column/gentle_hotcorona_evolution.mp4 1e9 1 1
```

- `gentle_hotcorona_evolution.mp4` — IC corona boosted ×1.5 → **1.29 MK** (above the ~1.14 MK balance), steady q(T). **Result: a genuine but TRANSIENT conduction-driven gentle evaporation.** For t ≈ 0–160 s the over-hot corona conducts its excess down the TR and drives evaporation upflow pulses peaking at **+13 km/s (Mach ≈ 0.09 — subsonic ⇒ *gentle*, A&S78)**, interleaved with condensation downflow fronts (−11…−15 km/s). With **no sustained heating**, the corona then cools 1.29 → 1.03 MK back to balance, the upflow dies, and it settles into the same mild steady drainage (~−3 to −8 km/s coronal downflow) as the unheated baseline; coronal mass/EM dip during the transient then recover and hold. So a hot corona's conductive flux *does* drive gentle evaporation (confirming A&S78), but a one-time energy increment gives only a transient — a *sustained* upflow needs sustained heating (the §9a/§9c drivers). Same 7-panel split-axis layout as §9c; phase chip reads RELAXING (steady q(T)); slow-motion ×4 over the [0,200] s evaporation window.

### §9e. Ambient volumetric heating H(s) added to model_gentle (`gentle_Hs_evolution.mp4`)

The RTV static-balance term `H(s) = E_H0·exp(−(s−s0)/s_H)` (Aschwanden & Schrijver 2002) was added to `model_gentle` (`GENTLE_HEAT=1`, default off; calibration knobs `GENTLE_HEAT_TMAX_MK`/`_SH_FRAC`/`_E0`/`_S0_MM`, ramp `GENTLE_HEAT_ENHANCE`/`_T_ON`/`_RAMP`). It reuses the existing `enable_coronal_heating` machinery (physics.hpp::`coronal_heating_rate`, integrators.cpp::`apply_coronal_heating_stage`). When on, the imposed boundary **q(T) is turned OFF** — the corona is now heated *inside* the domain (RTV way), so heating both volume and boundary would double-drive.

```bash
# uniform-ish heating (s_H = 3L), target apex 1.3 MK, steady (no ramp)
GENTLE_NS=300 GENTLE_HEAT=1 GENTLE_HEAT_SH_FRAC=3.0 GENTLE_HEAT_TMAX_MK=1.3 \
    ./build/chromo_main outputs/model_column/gentle_Hs_uniform.txt full ionization model_gentle - 10 cooling
SLOW_T0=0 SLOW_T1=300 SLOW_FACTOR=4 .venv/bin/python util/animate_gentle_column.py \
    outputs/model_column/gentle_Hs_uniform.txt visualization/model_column/gentle_Hs_evolution.mp4 1e9 1 1
```

- `gentle_Hs_evolution.mp4` — **H(s) works (heats the corona 0.86 → 1.53 MK, drives a strong evaporation transient: mean +22 km/s, peak +43 km/s at t ≈ 100 s) but does NOT pin V≈0.** After the transient it settles to a **steady ~−6.85 km/s coronal downflow** (stable t=2000–4000 s) — essentially the same as the q(T)-only run (−6.1 km/s), just with a hotter corona. **Both footpoint (`SH_FRAC=0.4`) and near-uniform (`SH_FRAC=3.0`) heating give the same residual downflow.** Root cause is **geometric, not the heating**: `model_gentle` is an OPEN column with a Mach-capped outflow top, which sustains a draining circulation (mass in at the top, down, out at the TR) regardless of how the corona is heated. v1's clean V≈0 baseline (§9a, mean|V| 112 m/s) came from the **CLOSED-loop confinement** (reflecting apex), not from H(s) per se. So pinning V≈0 needs the closed loop (§9a/pfss GENTLE), which already achieves it; the open column retains an intrinsic ~7 km/s downflow even with the RTV heating term present. (This completes the downflow investigation: it is neither the φ_g reconstruction, nor the boundary closure, nor missing volumetric heating — for the open column it is the open-top geometry.)

---

## 10. Unresolved / historical outputs

No current script in `util/` or `visualization/` writes these names (verified by grep). They were produced by an older or ad-hoc invocation; regenerate with the nearest matching script (likely `visualize_flare.py` / `visualize_flare_frames.py` / `animate_loop_flare.py`) using an explicit output path, then record the exact command here.

- `flare_explosive_evolution.mp4`
- `flare_explosive_evolution_more_mechanisms.mp4`
- `flare_spike_resolution.png`
- `footpoint_TR_burndown.png`


---

## 11. Isentropic chromosphere → conduction-driven evaporation (`model_isentropic`, docs/studies/evaporation/gentle_evaporation_downflow.md)

> **Reconciliation note (2026-07-11, docs/design/scenario_reconciliation_plan.md):** `model_isentropic`
> was renamed to **`model_column`** (the unified default scenario) with `model_isentropic` kept as a
> backward-compat alias, so every command below still runs. The "bestwb" numeric toggles
> (`ISO_LOG_RECON`, `ISO_MC3`/`ISO_MC3_BETA`, `ISO_EQ_WB`, `ISO_INNER_WB`, `ISO_INNER_WB_RHO`,
> `ISO_C7_IC`, `ISO_C7_HSE`) are now **always on** and ignored if set — a bare `model_column` run
> reproduces the "bestwb" scheme. The analytic-isentrope toy IC (`ISO_T_BASE`/`ISO_N_BASE`/`ISO_F_ION`),
> the base sponge (`ISO_SPONGE_*`), `ISO_TRAC_TC_FIXED`, and the `ISO_DIAG_*` diagnostics were removed;
> the `iso_stage1`/`iso_stage2*` analytic-isentrope runs below are retired (use the C7-IC runs).

First-principles restart of gentle evaporation: isentropic (adiabatic) chromosphere IC, hydrostatic-pressure + temperature-jump + EOS-density upper BC, well-balanced reconstruction (`Grid::well_balanced`, default-off). Staged by `ISO_HEAT_FLUX`. Default resolution **ns = 400** (Δs ≈ 750 m); the relaxation residual velocity is O(Δs) truncation error, so the commands below produce ns=400 without setting `ISO_NS`.

### `iso_stage1_relaxation.png`
Stage 1 (single-fluid hydro only; conduction/cooling/ionization off): the isentropic atmosphere relaxes to a steady state — T linear vs height (R²≈0.9999, flattening near the ghosts), V≈0 (interior ~2 m/s, max 7.3 m/s at the bottom boundary cell), pressure hydrostatic. Three panels: T(h) IC vs relaxed + linear fit; V(h); p(h).
```bash
ISO_HEAT_FLUX=0 ./build/chromo_main outputs/_archive/iso_stage1.txt full no-ionization model_isentropic - 20 no-cooling
.venv/bin/python util/iso_plot_stage1.py
```

### `iso_stage2_evaporation.png`
Stage 2 (conduction + the upper-BC temperature jump on): the downward conductive flux heats the top of the column; with no radiative sink the heat front sweeps the whole thin domain (no localized evaporation upflow). Panels: T(h), V(h), ρ(h) at several times. Default jump a=10, b=5 (env `ISO_TJUMP_A/_B`); gentle variant a=2, b=1.5.
```bash
ISO_HEAT_FLUX=1 ISO_TJUMP_A=10 ISO_TJUMP_B=5 ./build/chromo_main outputs/_archive/iso_stage2.txt full no-ionization model_isentropic - 20 no-cooling
.venv/bin/python util/iso_plot_stage2.py
```

### `iso_stage2b_evaporation.png`
Stage 2b (conduction + temperature jump + radiative cooling `ISO_COOLING=1` + TRAC): the conduction-vs-radiation balance localizes a steady transition region (cool ~5 kK base, sharp ~20 kK top) and a condensation downflow front (≈−1.4 km/s) propagates down during TR formation; the evaporation upflow stays weak (frozen ionization + unresolved corona). Panels: T(h), V(h), ρ(h) at several times.
```bash
ISO_HEAT_FLUX=1 ISO_COOLING=1 ISO_TJUMP_A=10 ISO_TJUMP_B=5 ./build/chromo_main outputs/_archive/iso_stage2b.txt full no-ionization model_isentropic - 20 no-cooling
.venv/bin/python util/iso_plot_stage2b.py
```

### `iso_stage1_evolution.mp4`, `iso_stage2_evolution.mp4`, `iso_stage2b_evolution.mp4`
Time-evolution animations of the three staged `model_isentropic` runs (single linear height axis — the domain is a thin sub-Mm slab, no split x-axis). 2×2 panels: T(h), V(h) in km/s, ρ(h) and p(h) (both log). `util/animate_isentropic.py` takes any iso output, subsamples uniformly to `ANIM_MAX_FRAMES` (default 400) frames, and renders them in parallel via `_anim_parallel.py`. Stage 1 shows the relaxation to V≈0; Stage 2 the heat front sweeping the whole column; Stage 2b the TR localizing with the condensation downflow front.
```bash
.venv/bin/python util/animate_isentropic.py outputs/_archive/iso_stage1.txt
.venv/bin/python util/animate_isentropic.py outputs/_archive/iso_stage2.txt
.venv/bin/python util/animate_isentropic.py outputs/_archive/iso_stage2b.txt   # default input if omitted
```

### `iso_stage2_ext_evolution.mp4`
Stage 2 (heat flux + jump, no cooling) on the **extended, photosphere-anchored** domain: `ISO_C7_IC=1` replaces the analytic isentrope with the **real Model C7 atmosphere** (Avrett & Loeser 2008, Table 26) so the lower boundary sits at **h = 0** and the domain spans the whole chromosphere up to the TR base (**h = 0 → 1303 km**, ns = 600, Δs ≈ 2.2 km). Same physics/BC as `iso_stage2` (single-fluid, well-balanced, conduction + ghost-T jump a=10/b=5, no cooling/ionization). The IC reproduces C7 (T minimum ≈ 4400 K at ~560 km, n_H 1.19×10²³ → 3.8×10¹⁸ m⁻³). With no radiative sink the conductive heat front **sweeps down the whole column** until the entire chromosphere is heated, then relaxes toward the ghost-driven equilibrium (T>20 kK boundary marches 745 km → photosphere ~7–9 km; transient mid-column peak ~140 kK at ~2900 s subsides to Tmax ≈ 66 kK = 10·T_top by ~6100 s; base ~10 kK), with the draining downflow growing monotonically to ~−4 km/s (an early ~17 km/s top transient as conduction switches on). The lower boundary uses the **Neumann (insulating, zero-conductive-flux) temperature BC** (`ISO_INNER_T_NEUMANN=1`); vs the default Dirichlet (fixed-photospheric-T) base the difference is small (base ~10 kK vs ~9.4 kK — the dense photosphere's heat capacity dominates, and the run is top-driven). Note the naive Neumann *hydro* ghost (floating ghost pressure) is unstable with no sink; the BC is imposed inside the Stage-D conduction solver instead. **The run begins with a 1500 s adiabatic relaxation** (`ISO_RELAX_TIME=1500`: conduction + jump off, V=0 reservoir top) so the C7 IC settles to a near-isentropic hydrostatic V≈0 equilibrium (max|V|→~20 m/s; the relaxed profile develops a ~16 kK bump near h≈1130 km as C7's outward-rising entropy rearranges with the top pinned) before Stage 2 switches on (the fixed jump re-anchored to the relaxed top T). Stage 2 from the relaxed state reproduces the no-relax result, confirming the evaporation signature isn't an IC artifact. Run with `time_mult=60` (step-capped at 6×10⁵ steps; relax to t≈1500 s + Stage 2 to t≈6900 s). Single linear height axis (chromospheric slab, no corona). The C7 rows below h=1003 km live in `scenarios/model_isentropic.cpp`; `chromo_main` labels heights from `Grid::out_base_km` (0 km here).
```bash
# Relax adiabatically for 1500 s (ISO_RELAX_TIME), then Stage 2; Neumann (insulating) base T.
# Drop ISO_RELAX_TIME for immediate Stage 2, or ISO_INNER_T_NEUMANN for the Dirichlet-sink base.
ISO_C7_IC=1 ISO_RELAX_TIME=1500 ISO_INNER_T_NEUMANN=1 ISO_HEAT_FLUX=1 ISO_TJUMP_A=10 ISO_TJUMP_B=5 ./build/chromo_main outputs/model_column/iso_stage2_ext.txt full no-ionization model_isentropic - 60 no-cooling
.venv/bin/python util/animate_isentropic.py outputs/model_column/iso_stage2_ext.txt
```

### `iso_t22k_evolution.mp4` — C7 chromosphere extended to the 22 kK TR base (h≈2153 km)
Photosphere-anchored real-C7 domain extended **up to the height where C7 reaches 22 kK** — the lower-TR base at **h ≈ 2152.6 km** (interpolated between the 2152 km/20.5 kK and 2153 km/23.1 kK table rows; top of the standard C7 column), `ISO_DH=2152.6`, `ISO_NS=1000`, Δs ≈ 2.15 km. Unlike the earlier 0→1300 km version (which truncated at 6.6 kK and imposed 22 kK as a hot reservoir above a cool top), this **resolves the lower TR rise (6.6 → ~19 kK) inside the domain**, with the boundary sitting at the natural 22 kK level (`ISO_T_TOP=22000` pins the outer ghost there). Density is rebuilt hydrostatically (`ISO_C7_HSE=1`, keeping C7's T(h) and ionization fraction): **required** here — the raw off-HSE C7 density (~20–26% off hydrostatic) launches a runaway flow at the steep TR gradient and NaNs in a few hundred steps; the HSE rebuild gives a clean V=0 start and runs stably the full t ≈ 13 460 s (no NaN, no TRAC/cooling needed at these ≤22 kK temperatures). Conduction on from t=0 (`ISO_HEAT_FLUX=1`), inner-BC well-balanced (`ISO_INNER_WB=1`), log-space reconstruction (`ISO_LOG_RECON=1`), default Dirichlet base, single fluid, no cooling/ionization. **Result:** T holds ~6.6 kK through the chromosphere (h ≲ 1900 km) then rises through the resolved TR to ~19 kK at the top cell; a **~+1 km/s upflow** develops in the TR (the gentle-evaporation–like conductive response, now resolved rather than imposed), with the base near V≈0 aside from the residual ~16 m/s conduction-driven drainage. Single linear height axis (2.15 Mm < 3 Mm split threshold; the TR occupies the top ~15% of the axis). `util/animate_isentropic.py` 2×3 panels: T, V [km/s], conductive heat flux q [kW/m²], ρ (log), p (log), mass flux ρV [kg/m²/s]; 400 frames @ 20 fps.
```bash
# Domain extended to where C7 T=22 kK (h≈2152.6 km, the lower-TR base); HSE-rebuilt density (ISO_C7_HSE)
# is REQUIRED — raw off-HSE C7 density NaNs at the steep TR.
# Current movie is ns=4000 (Δs≈0.54 km) to t≈2080 s — a refinement run. ns=4000's CFL timestep is 4×
# smaller, so reaching the full 13 500 s is ~16× the ns=1000 wall-time (impractical); 2080 s captures
# the rise→peak→decay-onset of the upflow. ns=1000/time_mult=60 reaches the full ~13 460 s.
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_NS=4000 ./build/chromo_main outputs/model_column/iso_t22k.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2080s
.venv/bin/python util/animate_isentropic.py outputs/model_column/iso_t22k.txt
.venv/bin/python util/plot_iso_frame.py outputs/model_column/iso_t22k.txt   # final-frame snapshot

# γ=1.05 (near-isothermal) variant at ns=4000. ISO_GAMMA must be passed to BOTH the run and the animator
# (the .txt stores no γ; default decode is 5/3). Manually stopped at t≈2140s.
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_NS=4000 ISO_GAMMA=1.05 ./build/chromo_main outputs/model_column/iso_t22k_gamma105.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2140s
ISO_GAMMA=1.05 .venv/bin/python util/animate_isentropic.py outputs/model_column/iso_t22k_gamma105.txt visualization/model_column/iso_t22k_evolution_gamma105.mp4   # -> iso_t22k_evolution_gamma105.mp4

# ns=2000 companion (same config, half resolution; for the grid-convergence comparison). Auto-stopped ~t=2208s.
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_NS=2000 ./build/chromo_main outputs/model_column/iso_t22k_ns2000.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2080s (overshot to 2208s)
.venv/bin/python util/animate_isentropic.py outputs/model_column/iso_t22k_ns2000.txt visualization/model_column/iso_t22k_evolution_ns2000.mp4   # -> iso_t22k_evolution_ns2000.mp4

# γ=1.05 (near-isothermal) variant of the SAME ns=2000 run. ISO_GAMMA sets grid.gamma_mono for the whole
# run (IC/BC/EOS). The plot decode MUST get the same γ (the .txt stores no γ): pass ISO_GAMMA=1.05 to the
# animator too. Manually stopped at t≈2300s (≈ baseline's t≈2208s) for a fair γ=5/3-vs-1.05 comparison.
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_NS=2000 ISO_GAMMA=1.05 ./build/chromo_main outputs/model_column/iso_t22k_ns2000_gamma105.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2300s
ISO_GAMMA=1.05 .venv/bin/python util/animate_isentropic.py outputs/model_column/iso_t22k_ns2000_gamma105.txt visualization/model_column/iso_t22k_evolution_ns2000_gamma105.mp4   # -> iso_t22k_evolution_ns2000_gamma105.mp4
```
`iso_t22k_evolution_ns2000_gamma105.mp4` repeats the ns=2000 conduction-only run with a **reduced adiabatic index γ=1.05** (vs 5/3) — a near-isothermal polytropic closure that stands in for the radiative sink this conduction-only run leaves out (c_v=k/(γ−1) is ~13× larger, so the gas absorbs conductive heat at nearly constant T). γ is now a single knob (`ISO_GAMMA`) threaded through the whole code: every energy↔pressure / heat-capacity factor derives from `grid.gamma_mono` via `Grid::gm1()`/`inv_gm1()`/`half_gm1()` (γ=5/3 reproduces the byte-baseline, 5958/5958). **Result vs the γ=5/3 baseline** (both at t≈2200–2300 s): the upper chromosphere (~1900 km) stays **pinned at ~6.6 kK** (init 6661 → final 6616 K) instead of conduction-heating to **~9.8 kK** (6669 → 9805 K) — the implicit thermostat — while a gentle evaporation upflow still develops (**+0.95 km/s** vs +1.05 km/s, Mach≈0.07). The top-cell density fills only ×1.17 (vs ×5.99). No NaNs (672 frames). (NB: the ion–neutral collisional energy-exchange coefficient 3 k_Bα/(m_i+m_n) is a kinetic-theory factor and is intentionally NOT tied to `gamma_mono`.)
`iso_t22k_evolution_ns2000.mp4` is the **ns=2000** companion to `iso_t22k_evolution.mp4` (ns=4000): identical production config, Δs ≈ 1.08 km (2× coarser), run to a comparable t≈2208 s. Use the pair to see the resolution dependence — the base downflow is ~2× larger at ns=2000 than ns=4000 (O(Δs) truncation), while the +1 km/s evaporation upflow is grid-converged.

**γ=1.05 near-isothermal variants** (`iso_t22k_evolution_gamma105.mp4` at ns=4000, `iso_t22k_evolution_ns2000_gamma105.mp4` at ns=2000). Reduced adiabatic index (`ISO_GAMMA=1.05`, vs 5/3) as a polytropic stand-in for the radiative sink this conduction-only run omits: c_v=k/(γ−1) is ~13× larger, so the gas absorbs the conductive heat flux at nearly constant T. γ is a single knob threaded through the whole code (`grid.gamma_mono` via `Grid::gm1()/inv_gm1()/half_gm1()`; γ=5/3 reproduces the 5958/5958 byte-baseline). **vs the γ=5/3 baseline** (compared at t≈2.1 ks): the upper chromosphere (~1900 km) stays **pinned ~6.6 kK** (6661→6613 K at ns=4000) instead of conduction-heating to **~9.7 kK** (6663→9727 K) — the implicit thermostat — while a gentle evaporation upflow still develops (**+0.95 km/s** vs +1.11; Mach≈0.07) and the top fills only ×1.19 (vs ×5.48). **Grid-converged:** ns=2000 vs ns=4000 γ=1.05 agree to ~0.5% on the upflow (+0.95 km/s both), T pinning (6616 vs 6613 K), and ρ_top (×1.17 vs ×1.19); only the (numerical, O(Δs)) base downflow halves with resolution (−8.1 → −4.4 m/s). The γ=5/3 baselines (`iso_t22k.txt`, `iso_t22k_ns2000.txt`) are preserved — γ=1.05 written to separate `*_gamma105.txt` files.

### `iso_t22k_evolution_ns2000_gamma105_mc3.mp4` / `..._mc3_wbrho.mp4` — base mass-flux fixes (MC3 limiter + inner-BC density well-balancing)
Two ns=2000 γ=1.05 variants probing the **base mass-flux truncation error** (the O(Δs) chromospheric downflow near the ρV zero crossing). Both extend the standard iso_t22k config; decode with `ISO_GAMMA=1.05`.
- **`..._mc3.mp4`** — MC3 / Koren limiter (`ISO_MC3=1 ISO_MC3_BETA=2`) instead of the default minmod. The asymmetric third-order (κ=1/3) limiter (BATSRUS `mc3`; `flux_lim_mc3_plus/minus` in `src/state.cpp`, gated by `grid.mc3_limiter`) is less diffusive, so it clips the under-resolved TR gradient less. Flat-ρ inner ghost (no density fix).
- **`..._mc3_wbrho.mp4`** — MC3 **plus** the inner-BC density well-balancing (`ISO_INNER_WB_RHO=1`), the current **best** config. The default inner reservoir pins the ghost density FLAT while the outer ghost is EOS-stratified; that flat ρ makes the inner-face reconstruction give ρ_L≠ρ_R → a Rusanov/LLF diffusive MASS flux at V=0 → a spurious base-cell downflow (see `docs/studies/boundaries/boundary_conditions_plan.md`). `ISO_INNER_WB_RHO` sets the ghost density from the EOS at the reservoir T₀ (ρ∝p, isothermal), mirroring the outer ghost, so ρ_L=ρ_R and the leak → O(Δs²).

**Base ρV at t≈2300 s (ns=2000, γ=1.05), by config** (min ρV [kg m⁻² s⁻¹] / base V[0]):

| config | min ρV | V[0] | TR upflow |
|---|---|---|---|
| flat-ρ, minmod (orig baseline) | −1.62×10⁻³ | −8175 m/s | +0.954 km/s |
| flat-ρ, **MC3** | −9.73×10⁻⁴ | −4914 m/s | +1.002 km/s |
| **fix** (WB_RHO), minmod | −4.17×10⁻⁴ | −2098 m/s | +0.947 km/s |
| **fix + MC3** | −4.19×10⁻⁴ | −1589 m/s | +0.994 km/s |

The inner-BC density fix is the dominant lever (removes the boundary spike, −74% on min ρV); MC3 stacks on top (further lowers V[0]); the physical TR upflow (~+0.95–1.0 km/s) is preserved throughout. `ISO_INNER_WB_RHO` (requires `ISO_INNER_WB`) and `ISO_MC3` are both default-off (test baseline 5972/5972 unchanged). Runs were manually stopped near t≈2300–2420 s.
```bash
# MC3 limiter only (flat-ρ inner ghost)
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_NS=2000 ISO_GAMMA=1.05 ISO_MC3=1 ISO_MC3_BETA=2 ./build/chromo_main outputs/model_column/iso_t22k_ns2000_gamma105_mc3.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2300s
ISO_GAMMA=1.05 .venv/bin/python util/animate_isentropic.py outputs/model_column/iso_t22k_ns2000_gamma105_mc3.txt visualization/model_column/iso_t22k_evolution_ns2000_gamma105_mc3.mp4

# fix (inner-BC density WB) + MC3 — best config
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_NS=2000 ISO_GAMMA=1.05 ./build/chromo_main outputs/model_column/iso_t22k_ns2000_gamma105_mc3_wbrho.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2300s
ISO_GAMMA=1.05 .venv/bin/python util/animate_isentropic.py outputs/model_column/iso_t22k_ns2000_gamma105_mc3_wbrho.txt visualization/model_column/iso_t22k_evolution_ns2000_gamma105_mc3_wbrho.mp4
```

### `iso_t22k_evolution_ns2000_gamma105_bestwb.mp4` + `..._bestwb_condcompare.mp4` — full well-balancing stack + conduction on/off mass-flux demo
The **best well-balanced** ns=2000 γ=1.05 config: `ISO_INNER_WB=1` (2nd-ghost pressure) + `ISO_INNER_WB_RHO=1` (ghost density) + `ISO_EQ_WB=1` (equilibrium-reference δ-form: subtract the frozen explicit-RHS residual at the IC so V=0 is an EXACT discrete fixed point for the non-isothermal C7 column) + `ISO_MC3=1` (Koren limiter) + `ISO_LOG_RECON=1`. Together these remove the O(Δs) numerical base-drainage entirely (root-caused in `docs/studies/boundaries/boundary_conditions_plan.md`).
- **`..._bestwb.mp4`** — the newest best result (conduction ON): the standard 2×3 `animate_isentropic.py` panels, decoded with `ISO_GAMMA=1.05`. The base ρV panel now sits ~flat while the TR evaporation upflow (+0.90 km/s) proceeds.
- **`..._bestwb_condcompare.mp4`** — the well-balancing demonstration: overlays ρV from two best-WB runs differing ONLY in conduction. `util/animate_iso_condcompare.py` (3 panels: ρV full domain, ρV zoomed h>1800 km, V zoomed; ρV/V decoded from RHO/MOM ⇒ NO γ needed). **conduction OFF stays flat at ~round-off everywhere (ρV∈[−5×10⁻⁶,+1×10⁻⁵], mean|V|base≈35 m/s, "upflow" ≈2 m/s) — the equilibrium is held to round-off with no driver — while conduction ON develops the physical TR evaporation (ρV up to +2.7×10⁻⁴, upflow +0.90 km/s).** So with the stack on, essentially ALL the mass flux is physical (conduction-driven); the former base "downflow" (−8×10³ m/s at cell 0 in the original flat-ρ minmod run) was numerical truncation and is gone.
```bash
# best well-balanced, conduction ON — the "newest best result" (item 1)
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_NS=2000 ISO_GAMMA=1.05 ./build/chromo_main outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2300s
ISO_GAMMA=1.05 Q_PHYS_YMAX_WM2=0.05 ANIM_MAX_FRAMES=400 MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/animate_isentropic.py outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt visualization/model_column/iso_t22k_evolution_ns2000_gamma105_bestwb.mp4 20

# same config, conduction OFF (ISO_HEAT_FLUX=0) — the equilibrium reference for the comparison
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=0 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_NS=2000 ISO_GAMMA=1.05 ./build/chromo_main outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condoff.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2300s

# mass-flux comparison movie (item 2): conduction ON vs OFF overlaid (defaults to the two files above)
.venv/bin/python util/animate_iso_condcompare.py   # -> visualization/model_column/iso_t22k_ns2000_gamma105_bestwb_condcompare.mp4
```

### `iso_t22k_evolution_ns2000_gamma105_tracon.mp4` — TRAC "try" experiment (forced broadening)
Same best-WB ns=2000 γ=1.05 config as `..._bestwb.mp4`, but with TRAC **force-enabled** via the new default-off knobs (`ISO_TRAC=1`, `ISO_TRAC_TCHROM=8000`, `ISO_TRAC_TC_FIXED=22000`) so the charged-fluid conductivity κ_e is broadened by ε=(T_c/T)^{5/2} throughout the 8–22 kK near-top rise. **Result: no measurable change.** Adaptive TRAC does not engage on this domain (T_c pins to the floor — the sub-22 kK profile is well-resolved at ns=2000, so the δ=½ criterion never fires, and Johnston's 0.2·T_peak cap forbids T_c from reaching the TR; the steep TR TRAC targets lives above the 22 kK top, outside the domain). Even the forced broadening changes the solution by <0.3% (global max|V| 1215→1218 m/s at the top). The standard 2×3 `animate_isentropic.py` panels are visually identical to `..._bestwb.mp4`: tiny chromospheric flow (V@250≈+6.6 m/s, peak|V|[60–700 km]≈9 m/s) and the physical +0.9 km/s A&S78 evaporation upflow at the 22 kK top. Confirms TRAC is neither triggered nor needed here (it belongs on the resolved-TR `ISO_CORONA` domain). New knobs: `ISO_TRAC` (force enable/disable, -1=auto=cooling‖corona), `ISO_TRAC_TCHROM` (T_b), `ISO_TRAC_TCMAXFRAC` (lift the 0.2·T_peak cap → `Grid::trac_Tc_max_frac`), `ISO_TRAC_TC_FIXED` (diagnostic fixed cutoff → `Grid::trac_fixed_cutoff_T`).
```bash
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_NS=2000 ISO_GAMMA=1.05 ISO_TRAC=1 ISO_TRAC_TCHROM=8000 ISO_TRAC_TC_FIXED=22000 ./build/chromo_main outputs/model_column/iso_t22k_ns2000_gamma105_tracon.txt full no-ionization model_isentropic - 4 no-cooling
ISO_GAMMA=1.05 .venv/bin/python util/animate_isentropic.py outputs/model_column/iso_t22k_ns2000_gamma105_tracon.txt visualization/model_column/iso_t22k_evolution_ns2000_gamma105_tracon.mp4
```

### `iso_t22k_250km_convergence.png` — the ~250 km mass-flux maximum is NUMERICAL (O(Δs))
Grid-convergence of the positive lower-chromosphere maximum that survives well-balancing (seen in the full-domain ρV panel of `..._bestwb_condcompare.mp4`, and present in BOTH bestwb and nowb). Three bestwb conduction-ON runs at **identical config**, ns=2000/3000/4000 (Δs=1.08/0.72/0.54 km), compared at t≈100 s. **Both V@250 and the peak ρV(h∈[40,700] km) scale ∝ Δs (first order) and extrapolate to ≈0 as Δs→0** (V@250: 10.7→7.1→5.4 m/s, Δs→0 intercept +0.07 m/s; ρV peak 8.5→5.8→4.3×10⁻⁴, intercept +0.16×10⁻⁴), and the ρV-peak location drifts toward the boundary (152→126 km). ⇒ the maximum is a **split conduction-to-hydro coupling artifact**, not a physical flow. It is conduction-gated (absent when conduction is off) but well-balancing can't touch it: WB balances the explicit hydrostatic/advective operator (a static, conduction-off property), whereas this is exposed when the next hydro stage differentiates the pressure changed by the separate conduction stage — hence identical in bestwb and nowb. Two panels: V@250 vs Δs and peak ρV vs Δs, each with a linear fit and the Δs→0 star. **The physical +0.9 km/s TR evaporation upflow is ~100× larger and grid-converged — unaffected.** Numerical remedies are grid/local refinement or improved conduction–hydro coupling; cooling is a physical sensitivity test, not a numerical fix.
```bash
# three matched bestwb cond-ON runs (delete the large .txt after plotting)
for NS in 3000 4000; do ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_NS=$NS ISO_GAMMA=1.05 ./build/chromo_main outputs/model_column/iso_bestwb_condon_ns$NS.txt full no-ionization model_isentropic - 2 no-cooling; done
.venv/bin/python util/plot_iso_250km_convergence.py 100   # ns=2000 uses the existing bestwb_condon.txt
```

### `iso_cond_hydro_convergence_t100.png` — matched-time split conduction-to-hydro diagnostic
Uses the Stage-D diagnostic CSVs from the best-WB ns=2000/3000/4000 sweep. It linearly interpolates each profile to t=100 s, then plots the time-normalized total conduction pressure change, total momentum-RHS difference, and total mass flux near 250 km. Each row includes a value-at-250-km versus Δs linear fit.
```bash
.venv/bin/python util/plot_cond_hydro_diagnostic.py --time 100 --out visualization/_archive/iso_cond_hydro_convergence_t100.png \
  outputs/_archive/iso_cond_hydro_ns2000.csv outputs/_archive/iso_cond_hydro_ns3000.csv outputs/_archive/iso_cond_hydro_ns4000.csv
```

### `iso_cond_hydro_local4_convergence_t100.png` — static-local-refinement split conduction-to-hydro diagnostic
The same best-WB conduction-ON diagnostic, but with **static local refinement** (`docs/studies/numerics/static_local_refinement.md`): each run refines the lower domain 0–700 km by 4× and retains its outer coarse spacing (transition 700–800 km, ratio ≤ 1.1). The coarse-equivalent counts `ISO_NS=2000/3000/4000` therefore share the SAME outer column while the local cell width at 250 km shrinks 4× (to the fine spacing 0.269/0.179/0.135 km), so the fit is against the **interpolated local `cell_width_m` at 250 km** — not the median domain spacing. **Result (t≈100 s):** the total mass flux ρ_iV+ρ_nU at 250 km is 6.37/4.50/3.24×10⁻⁴ kg m⁻² s⁻¹, a clean linear scaling in the local Δs (R²=0.994) with a Δs→0 intercept 2.5×10⁻⁵ (≈0). This isolates the ~250 km positive-ρV maximum as a **split conduction-to-hydro coupling artifact** (conduction-gated, well-balancing-invariant; scales with the local cell width and vanishes under local base refinement) rather than a physical flow or a pure conduction-operator artifact. The per-step dp_cond/dt at 250 km instead *rises* on finer grids (physical local-gradient resolution, R²≈0.74) and the total momentum-RHS difference sits at the ~10⁻⁷ noise floor (R²≈0.01). The physical +0.9 km/s TR evaporation upflow is ~100× larger and grid-converged. Cooling is a physical sensitivity test, not part of the numerical remedy (grid/local refinement or improved conduction–hydro coupling).
```bash
# three matched best-WB conduction-ON runs, refined 0–700 km by 4× (delete /tmp .txt after).
# time_mult is set so the 10000·time_mult step cap reaches ≳100 s at the fine-cell dt.
for pair in "2000 500 2.0" "3000 750 3.0" "4000 1000 4.0"; do set -- $pair; NS=$1 EVERY=$2 TM=$3
  ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 \
    ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_GAMMA=1.05 \
    ISO_REFINE_FACTOR=4 ISO_REFINE_S_LO_KM=0 ISO_REFINE_S_HI_KM=700 ISO_REFINE_TRANSITION_KM=100 \
    ISO_NS=$NS ISO_DIAG_COND_HYDRO=1 ISO_DIAG_COND_HYDRO_EVERY=$EVERY \
    ISO_DIAG_COND_HYDRO_OUT=outputs/_archive/iso_cond_hydro_local4_outer$NS.csv \
    ./build/chromo_main /tmp/_iso_local4_$NS.txt full no-ionization model_isentropic - $TM no-cooling; done
.venv/bin/python util/plot_cond_hydro_diagnostic.py --time 100 \
  --out visualization/_archive/iso_cond_hydro_local4_convergence_t100.png \
  outputs/_archive/iso_cond_hydro_local4_outer2000.csv outputs/_archive/iso_cond_hydro_local4_outer3000.csv outputs/_archive/iso_cond_hydro_local4_outer4000.csv
```

### `iso_dynamic_mix_hse_t100.png` — total-mixture dynamic-HSE diagnostic
Bottom-anchored **total-mixture** HSE reference (`p_mix`, `rho_mix`) for the best-WB ns=2000/3000/4000 sweep. The script linearly interpolates each diagnostic to t=100 s, plots `dp_cond_mix/dt`, `dpi_cond_mix/dt`, post-conduction mixture perturbation pressure, and total mass flux over 0–700 km, and prints both 250-km fits and window max/L1/peak metrics. The charged/neutral pressure fractions are partition diagnostics only, not separate HSE references.
```bash
.venv/bin/python util/plot_dynamic_mix_hse_diagnostic.py --time 100 --out visualization/_archive/iso_dynamic_mix_hse_t100.png \
  outputs/_archive/iso_dynamic_mix_hse_ns2000.csv outputs/_archive/iso_dynamic_mix_hse_ns3000.csv outputs/_archive/iso_dynamic_mix_hse_ns4000.csv
```

### `iso_t22k_ns2000_gamma105_innerwb_condcompare.mp4` — inner-BC WB only (eq-WB removed): isolates what eq-WB does
Same as `..._bestwb_condcompare.mp4` (ρV conduction ON vs OFF overlaid, 3 panels: full-domain ρV, TR ρV h>1800 km, TR V h>1800 km) but with **eq-WB removed** — keeps `ISO_INNER_WB=1`+`ISO_INNER_WB_RHO=1`+`ISO_MC3`+`ISO_LOG_RECON`, drops `ISO_EQ_WB`. **What it demonstrates:** with only the inner-BC well-balancing, the conduction-OFF trace is **no longer flat** — both cond-ON and cond-OFF droop to ρV≈−3.8×10⁻⁴ at the base (vs the round-off-flat blue line in the bestwb movie). That base drainage is the **O(Δs) non-isothermal reconstruction truncation** that eq-WB cancels; it is present in BOTH conduction states ⇒ static / conduction-independent. The hierarchy (cond-OFF base min ρV): no WB −2.3×10⁻³ → **inner-WB only −3.8×10⁻⁴** → full bestwb (eq-WB on) −4.6×10⁻⁶. The positive ~250 km bump appears only in cond-ON (the split conduction-to-hydro coupling artifact, unaffected by eq-WB), and the physical +1 km/s TR evaporation upflow shows in the cond-ON TR panels. Title set via `ANIM_LABEL`. Runs use `tm=120` and were stopped (`TaskStop`) at t≈2650 s to match the bestwb pair's duration.
```bash
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_NS=2000 ISO_GAMMA=1.05 ./build/chromo_main outputs/model_column/iso_t22k_ns2000_gamma105_innerwb_condon.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2650s
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=0 ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_NS=2000 ISO_GAMMA=1.05 ./build/chromo_main outputs/model_column/iso_t22k_ns2000_gamma105_innerwb_condoff.txt full no-ionization model_isentropic - 120 no-cooling   # stop ~t=2650s
ANIM_LABEL="inner-BC WB only (INNER_WB+RHO, MC3; NO eq-WB)" .venv/bin/python util/animate_iso_condcompare.py outputs/model_column/iso_t22k_ns2000_gamma105_innerwb_condon.txt outputs/model_column/iso_t22k_ns2000_gamma105_innerwb_condoff.txt visualization/model_column/iso_t22k_ns2000_gamma105_innerwb_condcompare.mp4
```

### `iso_t22k_gamma105_massflux_compare.mp4` — γ=1.05 mass-flux grid convergence (ns=1000/2000/4000)
Overlays the **three γ=1.05 resolutions** (ns=1000 `outputs/model_column/iso_t22k_ns1000_gamma105.txt`, ns=2000, ns=4000) to show grid convergence of the mass flux ρV. Three panels, all three resolutions overlaid: **(a)** ρV over the full domain (h>1800 km shaded), **(b)** ρV zoomed to **h>1800 km** (the TR), **(c)** V zoomed to h>1800 km. Synced on a common time axis [0, min(t_max)=2141 s]; nearest-frame matching across the runs' different cadences. ρV=(ρ_i+ρ_n)·V and V are pure RHO/MOM decodes (**no γ** — so this animator needs no ISO_GAMMA, unlike animate_isentropic.py). **What it shows:** in the full-domain panel the dense-base ρV (h≲500 km) **diverges** with resolution (ns=1000 largest excursions, ns=4000 smallest) — the non-converging O(Δs) numerical/relaxation artifact; in the h>1800 km panels the **TR evaporation** (ρV≈5–7×10⁻⁷, V→~0.95 km/s) sits nearly on top of each other across all three grids — grid-converged & physical. (Quantitatively, mean over t∈[1000,2000] s: TR max ρV 7.46→7.19→6.98×10⁻⁷ for ns=1000/2000/4000 — converging; global/base max ρV 1.97→1.72→1.15×10⁻⁴ — successive differences *grow*, i.e. NOT converging.) `util/animate_iso_gamma_compare.py [out.mp4]`; 300 frames @ 20 fps.
```bash
.venv/bin/python util/animate_iso_gamma_compare.py   # -> visualization/model_column/iso_t22k_gamma105_massflux_compare.mp4
```
**Resolution dependence of the base downflow (ns=1000 → 4000).** At a fixed time (t≈2000 s) the base downflow falls **≈−14 → ≈−3.3 m/s** — a ~4× reduction for 4× finer grid, i.e. **≈O(Δs)**. So the base downflow is *dominantly numerical*: the under-resolved steep TR produces a spurious conduction-driven compression that propagates down as a first-order truncation artifact, larger than the conduction-OFF residual alone implied. The **upflow** (gentle-evaporation signature, V_top≈+1 km/s, ρV peak ≈+2.9×10⁻⁴) is **grid-converged** (consistent with the earlier ns=600≈1200 check) and decays in time via self-limiting conduction (the TR warms toward the 22 kK boundary ⇒ ∂T/∂s flattens ⇒ flux and overpressure fall). Net: upflow = physical & converged; base downflow = mostly numerical & shrinks ~linearly with Δs.
**Inner-BC well-balancing** (`ISO_INNER_WB=1`, `Grid` inner reservoir): the second inner ghost was a flat *copy* of the first, so the ghost-side pressure slope was zero while the interior was hydrostatic — the inner-face MUSCL reconstruction went first-order (limiter r=0), the reconstructed face pressure jumped by ½ρgΔs, and the boundary cell felt a spurious force ⇒ a base downflow **even in HSE with conduction off**. `ISO_INNER_WB=1` continues the reservoir hydrostatically (p₀+2ρgΔs), so the ghost-side slope matches the interior, r→1, and V=0 is an exact discrete fixed point. **Clean proof on the analytic isentrope** (exact Euler steady state, conduction off): bottom-cell V 7.3→2.5 m/s, spike eliminated (base now = interior O(Δs) drift). In this conduction-on C7 run the base improves 20.9→15.0 m/s (mass-flux spike −4.1→−3.0×10⁻³); the residual is the C7-IC relaxation flow + genuine conduction drainage, not the inner BC. Default-off keeps the test baseline byte-for-byte (5958/5958).
**Log-space reconstruction** (`ISO_LOG_RECON=1`, `Grid::log_reconstruct`): the explicit MUSCL step limits log ρ / log p (the exponentially-stratified slots) instead of ρ / p; V, U stay linear. Default-off keeps the test baseline byte-for-byte (5958/5958). Effect here is **concentrated at the TR/top** (where the gradient is under-resolved and reverses — ρ differs ~7% at h≈1300 km, Tmax 17.9→18.2 kK), and **negligible in the chromosphere/base**: ds/H ≈ 1% there, so linear reconstruction is already accurate and minmod doesn't clip a monotonic exponential (r=e^(ds/H)>1 ⇒ full slope). The ~−20 m/s base downflow is therefore **not** ρ-limiter diffusion — it's the inner-BC drainage + the known O(Δs) bottom-boundary-cell truncation artifact (≈7 m/s even in Stage-1 with conduction off).

### `iso_t22k_final_frame.png` — final-frame snapshot (zoomed mass flux)
Static 6-panel snapshot of the **final frame** (t ≈ 13 460 s) of the `iso_t22k` run, same panels as the movie (T, V, q, ρ, p, ρV). The Temperature panel shows the flat ~6.6 kK chromosphere then the resolved TR rise to ~19 kK in the top ~350 km; the conductive-flux panel shows the downward flux concentrated at the TR. The mass-flux panel ρV is **zoomed** to its own 0.5–99.5 percentile range (the movie's shared scale is dominated by the early transient): downward **drainage** through the chromosphere (ρV down to ~−9×10⁻⁴ kg/m²/s at the dense base) easing toward ~0 in the TR, with only a thin near-top upflow. `util/plot_iso_frame.py [input.txt] [frame_index=-1] [output.png]` plots any single frame of any iso run as a static PNG (mass-flux panel always autoscaled robustly to that frame).
```bash
.venv/bin/python util/plot_iso_frame.py outputs/model_column/iso_t22k.txt
```

### `iso_corona_evolution.mp4`
**Resolved-corona evaporation testbed** (`ISO_CORONA=1`): the single-fluid model_isentropic domain is extended from the photosphere up into a **resolved ~1 MK corona** (h = 0 → 10 Mm) initialized from the real C7 profile with a **realistic frozen ionization profile** (n_i=n_e, n_n=n_HI via `c7_full_profile` — neutral chromosphere, ionized corona) and TRAC on. This is what was missing from the chromosphere-only `iso_stage2_ext` run: only an *ionized* corona conducts via Spitzer κ_e (a uniform frozen f leaves it neutral ⇒ κ_e ~10⁴× too small), and only a resolved coronal *volume* gives evaporated material somewhere to fill. Result: **conduction now drives a clear evaporation upflow** — V≈0 in the chromosphere rising through the TR to **+75 km/s in the corona** (T → ~1.2 MK), the Antiochos & Sturrock (1978) signature, opposite the chromosphere-only downflow. Conduction is on from t=0 (a conduction-off relaxation of the steep TR is unstable); HSE IC gives a clean V=0 start. No radiative cooling, no ionization dynamics, single fluid. **Split x-axis** (chromosphere+TR zoomed | corona compressed). The no-sink corona evaporates then drains and goes numerically unstable at the late drainage shock (~300 s); the animator stops at the last finite frame, so the movie covers the evaporation phase. **The ~300 s cap is physical, not a grid artifact:** at 2× resolution (`ISO_NS=1200`) the NaN moves only 302 → 319 s, and extra thermal diffusivity (`ISO_DIFF_MULT=12`) made it slightly *worse* (290 s) — it is the unmaintained 1 MK corona draining to vacuum at the top (a hydrodynamic, not thermal, instability). Running 3× longer would require the radiative sink (Stage 2b) so the corona reaches a quasi-steady state. This movie uses the longest stable window (the 1200-cell run, t→318.6 s). **2×3 panel layout:** T, V, **conductive heat flux** `q = -(κ_e+κ_n) dT/ds` (kW/m², >0 up; the un-broadened Spitzer flux computed in `animate_isentropic.py` from the same κ_e/κ_n as `physics.hpp`), ρ (log), p (log); the 6th cell is blank. The heat-flux panel shows the downward conductive flux spike at the TR — the driver of the A&S78 evaporation upflow — turning positive (upward) near the top where T falls toward the boundary.
```bash
ISO_CORONA=1 ISO_CORONA_TOP_KM=10000 ISO_C7_HSE=1 ISO_HEAT_FLUX=1 ISO_NS=1200 ./build/chromo_main outputs/model_column/iso_corona.txt full no-ionization model_isentropic - 30 no-cooling
.venv/bin/python util/animate_isentropic.py outputs/model_column/iso_corona.txt
```

### `iso_corona_full_evolution.mp4` — resolved corona with the FULL physics (radiation sink + two-fluid + ionization/recomb)
The resolved-corona run repeated with the three physics pieces the clean testbed left off — **radiative sink** (`ISO_COOLING=1`, Stage R optically-thick + thin), **two separate fluids** (`ISO_TWO_FLUID=1`, ion/neutral drag + T-equilibration instead of slaving), and the **ionization/recombination network** (`ISO_IONIZATION=1`, Stage E; the frozen C7 split is now just the START). When ionization is on, model_isentropic populates `grid.photoionization_rate_i` by the same **Route-B equilibrium inversion** as `model_c7_ic` (blended to the frozen Chae-2021 FAL-C rate above ~1500 km) so the chromosphere is a fixed point and doesn't drift. All three are env-gated and **default off** ⇒ the clean testbed above and the 5958-test baseline are byte-for-byte unchanged.

**The lower boundary is raised to h = 1003 km** (`ISO_H_BASE=1003`, model_c7's validated floor) rather than the photosphere: the Stage-E channel counts (`dt·n_i²·α`, `dt·n_i n_n·S`) and the radiative loss overflow **float32** at photospheric density (n ~ 10²³ ⇒ n² ~ 10⁴⁶ ≫ FLT_MAX ≈ 3.4×10³⁸), NaN-ing at step 0; at h ≥ 1003 km (n ≲ 10¹⁹) they are in range — the regime where model_c7's ionization/cooling are validated. The resolved corona above (1003 km → 10 Mm) is unchanged.

**Result — the sink fundamentally changes the corona's fate.** Where the no-sink testbed (above) ran a vigorous **+75 km/s** evaporation then crashed at ~300 s, the full-physics run is **stable for the entire 2445 s** (8× longer): the radiative sink regularizes the corona instead of letting it drain to vacuum. With **no sustained coronal heating**, though, the unmaintained corona sheds its energy by radiating and **slowly condenses/drains** — only a brief early evaporation transient (upflow building to ~3 km/s by t ≈ 200 s), then a persistent few-km/s coronal downflow (deepest, ~−15 km/s, during the conduction switch-on at t ≲ 120 s). Over the run the coronal (h > 2.3 Mm) mass ∫n_e ds falls **9.1 → 4.0 ×10²¹ m⁻² (−56%)** and the emission measure ∫n_e² ds **1.3×10³⁷ → 2.3×10³⁶ (−83%)**; the top stays pinned near its C7 value (~0.81 MK, continuous top ghost). This is the physically-correct picture A&S78 imply: radiation lets the corona lose energy, so an **unheated** corona cools and drains rather than evaporating — which is exactly why a sustained gentle evaporation needs ambient coronal heating (`model_gentle`'s H(s), §9). **Split x-axis** (chromosphere+TR zoomed | corona compressed); same 5-panel layout (T, V labeled up=evaporation/down=condensation, conductive heat flux, ρ log, p log).
```bash
# base raised to 1003 km (ISO_H_BASE) — ionization/cooling n² overflow float32 at the photosphere
ISO_CORONA=1 ISO_CORONA_TOP_KM=10000 ISO_H_BASE=1003 ISO_C7_HSE=1 ISO_HEAT_FLUX=1 \
    ISO_COOLING=1 ISO_IONIZATION=1 ISO_TWO_FLUID=1 ISO_NS=600 \
    ./build/chromo_main outputs/model_column/iso_corona_full.txt full ionization model_isentropic - 10 cooling
.venv/bin/python util/animate_isentropic.py outputs/model_column/iso_corona_full.txt
```

### `iso_combined_stage1_stage2.mp4`
Publication-styled (nature-figure) combined movie that concatenates **Stage 1 (relaxation)** then **Stage 2 (heat conduction enabled, = Stage 2b)** with a stage banner caption, persistently labelling it **single fluid**. 2×2 panel grid: **a** temperature (the TR forming), **b** velocity (downflow/condensation front, filled; minor gridlines + a +0.1 km/s tick), **c** density (TR cliff, log), **d** pressure (log). Common y-axes across both stages so the cool→hot contrast is honest; a q(T) arrow marks the conductive-flux boundary in Stage 2. Early times are sampled densely (power-law) to capture the transient/TR formation. Stage 2 = the `iso_stage2b` run (conduction + jump + radiative cooling). Single linear height axis (thin sub-Mm slab). 24 fps, ~10 s.
```bash
.venv/bin/python util/animate_iso_combined.py
# or: util/animate_iso_combined.py <stage1.txt> <stage2b.txt> <out.mp4> [fps]
```

### `docs/posters/2026-07-11-shine/figs/iso_t22k_snapshots.pdf` (+ `.png`) — SHINE poster Result 2 figure
Poster-styled (nature-figure) **initial-vs-final overlay** snapshot of the **γ=1.05** `iso_t22k_ns2000_gamma105` run (real Model-C7 atmosphere, photosphere → 22 kK lower-TR base, h = 0 → 2.15 Mm, HSE-rebuilt density, frozen ionization, conduction-driven from a hydrostatic V=0 start; **near-isothermal γ=1.05** as a polytropic stand-in for the omitted radiative sink). Four panels, each overlaying the **initial** (t=0, navy) and **final** (t≈2300 s, orange) frames on a single linear height axis (Mm): **(a)** Temperature — the upper chromosphere is **held near-isothermal at ~6.6 kK** (conductive heat absorbed at constant T; init/final overlap); **(b)** Velocity (hero) — a gentle **+0.95 km/s** subsonic upflow develops at the TR (Mach ≈ 0.07); **(c)** Density (log) — the TR/top fills modestly **~1.2×**; **(d)** Mass flux ρV — **whole-domain** view (the dense base carries a large **numerical O(Δs)** drainage spike that dominates the scale) **plus a zoomed inset for h > 1.8 Mm** (the TR), where the grid-converged evaporative upflow (~10⁻⁶) lives. The figure decode reads γ via `ISO_GAMMA` (default 1.05 here; `animate_isentropic.primitives` is reused, with `A.GM1` overridden). ρV and V are γ-independent; only T/p need the γ. UM navy/maize palette/fonts match `iso_corona_poster.py`. The upper-chromosphere+TR region is shaded. Written into `docs/poster-shine/figs/` (the poster's `\graphicspath`); this is the **Result 2 figure** in `poster.tex` (which was likewise updated to the γ=1.05 setup/numbers). `util/iso_t22k_poster.py [input.txt] [output.pdf]` (defaults to the γ=1.05 ns=2000 file).
```bash
.venv/bin/python util/iso_t22k_poster.py   # defaults to outputs/model_column/iso_t22k_ns2000_gamma105.txt, γ=1.05
# γ=5/3 version: ISO_GAMMA=1.6666667 .venv/bin/python util/iso_t22k_poster.py outputs/model_column/iso_t22k_ns2000.txt <out.pdf>
# rebuild the poster: cd docs/poster-shine && latexmk -lualatex -output-directory=build poster.tex
```

### `iso_local4_ns2000_bestwb_evolution.mp4` — validated short run, best-WB C7-HSE + static local refinement 0–700 km ×4
Validation re-run (2026-07-10) of the full best-well-balanced `model_isentropic` stack (`ISO_LOG_RECON`+`ISO_INNER_WB`+`ISO_INNER_WB_RHO`+`ISO_EQ_WB`+`ISO_MC3`, γ=1.05) combined with **static local refinement** (`docs/studies/numerics/static_local_refinement.md`) of the lower domain 0–700 km by 4× (transition 700–800 km, ratio ≤1.1), same config as the `iso_cond_hydro_local4_*` sweep. Build: `cmake --build build` (no changes, already up to date). Tests: `./build/chromo_tests` → **6351/6351 passed**. Sanity run (`ISO_NS=2000`, time_mult=1, step_cap=10000, t≈63 s, ~33 s wall) confirmed no NaN/Inf across all 1001 frames before committing to the longer run. Main run extends to t≈189 s (time_mult=3, step_cap=30000, ~101 s wall) with the split conduction-to-hydro diagnostic CSV enabled. **Result: numerically sane** — base held at the C7 photospheric 6578 K (well-balanced, no drift), upper chromosphere temperature-minimum (~4400 K @570 km) and plateau (~6600 K) intact, top relaxes toward the 20–22 kK target; velocity is ~round-off through most of the column with a physical **+0.86 km/s** evaporative upflow developing at the top (matches the documented bestwb signature); the well-known **positive ρV maximum near ~150–250 km** (the split conduction-to-hydro coupling artifact, `outputs/_archive/iso_cond_hydro_local4_ns2000_main.csv`, 604 804 rows, 0 NaN/Inf) is visibly present but small (peak ≈5.8–8.5×10⁻⁴ kg m⁻² s⁻¹, decaying over the run) — consistent with local refinement suppressing it relative to an unrefined ns=2000 run. All 3001 native frames finite; animation used 400 uniformly-sampled frames, none dropped.
```bash
# sanity check (t≈63 s, ~33 s wall) — run first, confirm no NaN before the longer run
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 \
    ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_GAMMA=1.05 \
    ISO_REFINE_FACTOR=4 ISO_REFINE_S_LO_KM=0 ISO_REFINE_S_HI_KM=700 ISO_REFINE_TRANSITION_KM=100 \
    ISO_NS=2000 \
    ./build/chromo_main outputs/_archive/iso_local4_sanity.txt full no-ionization model_isentropic - 1 no-cooling

# main run (t≈189 s, ~101 s wall) with the split conduction-to-hydro diagnostic CSV
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 \
    ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_GAMMA=1.05 \
    ISO_REFINE_FACTOR=4 ISO_REFINE_S_LO_KM=0 ISO_REFINE_S_HI_KM=700 ISO_REFINE_TRANSITION_KM=100 \
    ISO_NS=2000 ISO_DIAG_COND_HYDRO=1 ISO_DIAG_COND_HYDRO_EVERY=200 \
    ISO_DIAG_COND_HYDRO_OUT=outputs/_archive/iso_cond_hydro_local4_ns2000_main.csv \
    ./build/chromo_main outputs/_archive/iso_local4_ns2000_main.txt full no-ionization model_isentropic - 3 no-cooling

ISO_GAMMA=1.05 .venv/bin/python util/animate_isentropic.py \
    outputs/_archive/iso_local4_ns2000_main.txt visualization/_archive/iso_local4_ns2000_bestwb_evolution.mp4
```

### `iso_local4_ns2000_2000s_evolution.mp4` — same run extended to t≈2000 s
Same config as `iso_local4_ns2000_bestwb_evolution.mp4` above, extended to `time_mult=32` (step_cap=320000, dt≈0.0063 s) to reach **t≈2022 s** (~18 min wall time, run in background). Diagnostic CSV stride widened to `ISO_DIAG_COND_HYDRO_EVERY=1000` to keep the file size reasonable over the longer run (530 MB, 0 NaN/Inf). **Result: numerically sane over the full 2000 s** — all 5001 native frames finite (frame_stride=64), base pinned at the C7 photospheric ~6580 K the entire run (no drift), top settles near 18.3 kK (still slowly approaching the 22 kK target). A transient TR condensation dip (~−0.17 to −0.61 km/s, peaking near t≈280–350 s at heights 1300–1900 km) relaxes away by t≈2022 s into a smooth monotonic upflow topping out at +0.84 km/s. The **~250 km lower-chromosphere artifact continues decaying** as the run settles: peak ρV ≈5×10⁻⁴ kg m⁻²s⁻¹ at t≈353 s → ≈2.5×10⁻⁴ kg m⁻²s⁻¹ at t≈2022 s (half), and the value at 250 km itself falls from ≈4×10⁻⁸ (MOM_I, t≈189 s) to ≈2×10⁻⁸ (t≈2022 s) — consistent with it being a transient split conduction-to-hydro coupling artifact that keeps shrinking as the column relaxes toward a smoother quasi-steady balance, not a growing instability.
```bash
ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 \
    ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_GAMMA=1.05 \
    ISO_REFINE_FACTOR=4 ISO_REFINE_S_LO_KM=0 ISO_REFINE_S_HI_KM=700 ISO_REFINE_TRANSITION_KM=100 \
    ISO_NS=2000 ISO_DIAG_COND_HYDRO=1 ISO_DIAG_COND_HYDRO_EVERY=1000 \
    ISO_DIAG_COND_HYDRO_OUT=outputs/_archive/iso_cond_hydro_local4_ns2000_2000s.csv \
    ./build/chromo_main outputs/_archive/iso_local4_ns2000_2000s.txt full no-ionization model_isentropic - 32 no-cooling

ISO_GAMMA=1.05 .venv/bin/python util/animate_isentropic.py \
    outputs/_archive/iso_local4_ns2000_2000s.txt visualization/_archive/iso_local4_ns2000_2000s_evolution.mp4
```

---

## 12. Classical pure-H Saha / CRASH Γ₁ table (Stages 0–2)

`util/eos/tabulate_gamma.f90` builds the approved classical pure-hydrogen
cross-validation table with excitation, Fermi, and Coulomb corrections disabled.
The raw 11-column table retains `zAv`, energy γ_E, Γ₁, energy, and `Cv` for
validation; the tracked production table at
`data/eos/gamma1_hydrogen_v1.dat` contains only
`log10(T), log10(n_H), Γ₁`, with strict v1 metadata and a SHA-256 sidecar.
The grid is `T=3.2×10³–10⁸ K` (501 log points) by
`n_H=10¹²–10²⁶ m⁻³` (57 log points).

```bash
# Generate outputs/eos_gamma/gamma_hydrogen.dat and refresh the tracked
# data/eos/gamma1_hydrogen_v1.dat, verify its checksum, then run the strict
# whole-grid Saha/energy/Cv/Gamma1 validation.
bash util/eos/build_and_run.sh

# Generate all three figures below (project root as working directory).
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/plot_eos_gamma.py
```

### `eos_gamma/gamma_hydrogen.png` — γ_E, Γ₁, and ionization vs T

Energy γ_E (solid), sound-speed Γ₁ (dashed), and equilibrium ionization
for four densities spanning `10¹⁴–10²⁶ m⁻³`. Both indices recover to 5/3
at the neutral and fully ionized limits and soften through the Saha transition.

### `eos_gamma/gamma_hydrogen_2d.png` — production Γ₁(T,n_H) surface

Heatmap of the sound-speed index over the entire runtime grid with
`x={0.01,0.1,0.5,0.9,0.99}` ionization contours. This is the column loaded by
the Stage-2 C++ interpolator; CRASH pressure and energy are not runtime inputs.

### `eos_gamma/c7_eos_consistency.png` — C7 track and Stage-0 consistency

Walks the complete Model C7 atmosphere parsed from `scenarios/model_c7.cpp` and
shows C7 ionization against textbook Saha, CRASH γ_E and Γ₁, and the required
central-`ln T` `C_V^eff=(∂e_int/∂T)_rho` against CRASH `Cv`. The chromosphere/TR
is expanded and the extended corona compressed with a split height scale.

### `eos_gamma/c7_gamma_profile_classical_saha.png` — signed-off classical EOS along Model C7

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/plot_c7_gamma_profile.py --classical-saha
```

Preserves the earlier `c7_gamma_profile.png` and generates a separately named
height profile for the signed-off classical pure-H EOS. It overlays the
textbook-Saha energy index `γ_E=1+P/e` and the CRASH sound-speed index `Γ₁`
with excitation, Fermi, and Coulomb corrections disabled, alongside the same
Model C7 temperature and density context used by the original figure.

### `eos_gamma/c7_ionization_profile.png` — Model C7 ionization degree

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/plot_c7_ionization_profile.py
```

Model C7 hydrogen ionization degree `x=n_e/(n_HI+n_e)` vs height, alongside
the same temperature and mass-density context as `c7_gamma_profile.png`.

### `eos_gamma/c7_ionization_profile_classical_saha.png` — classical-Saha ionization degree

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/plot_c7_ionization_profile.py
```

Classical pure-H Saha ionization degree evaluated at the local Model C7
temperature and hydrogen-nuclei density, plotted separately from the C7 curve.

### `eos_gamma/c7_ionization_profile_compare.png` — C7 vs classical-Saha vs excitation-enabled CRASH

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/plot_c7_ionization_profile.py --mode compare
```

Direct overlay of the Model C7 ionization degree, analytic equilibrium
classical-Saha fraction, and the excitation-enabled CRASH code's tabulated
`zAv` evaluated along the same C7 atmosphere. An inset magnifies the
2119–2147 km transition-region interval where C7 and CRASH briefly diverge by
more than 20%; a small rectangle marks that interval on the main panel and an
arrow links it to the inset. The CRASH diagnostic table is generated separately,
preserving the signed-off excitation-off production table:

```bash
bash util/eos/build_and_run.sh excitation
```

The default plotting command writes all three ionization-profile figures.

### `eos_gamma/crash_excitation_full_table_c7.png` — full CRASH excitation comparison with C7 track

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/plot_crash_excitation_full_table.py
```

Three-panel view of the complete CRASH `(T,n_H)` grid: equilibrium hydrogen
ionization with bound excitation off, with excitation on, and the percentage
reduction caused by excitation. The complete Model C7 thermodynamic track is
overlaid on every panel and labeled at the photosphere, temperature minimum,
upper chromosphere, TR base, and upper TR. Contours identify the 1%, 10%, and
30% excitation-sensitive regions.

---

## 13. Pure-H Saha / CRASH-Γ₁ model-column evolution

### `model_column/model_column_saha_gamma_evolution.mp4`

Default 600-cell, well-balanced `model_column` evolved with the signed-off
classical pure-H EOS table. The supported gamma-table configuration uses the
full Euler integrator in single-fluid/common-temperature mode, with the
non-equilibrium ionization network and extra heating/cooling mechanisms off.
The EOS-aware `.gamma_diag` sidecar supplies physical temperature, pressure,
equilibrium composition, Γ₁, and conductivity to the six-panel animation.

```bash
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat ./build/chromo_main \
    outputs/model_column/model_column_saha_gamma.txt full no-ionization model_column - 1.0
Q_PHYS_YMAX_WM2=0.05 ANIM_MAX_FRAMES=400 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
    .venv/bin/python util/animate_isentropic.py \
    outputs/model_column/model_column_saha_gamma.txt \
    visualization/model_column/model_column_saha_gamma_evolution.mp4 20
```

### `model_column/model_column_saha_gamma_t22k_ns2000_evolution.mp4`

Extended-domain counterpart on the exact historical
`iso_t22k_evolution_ns2000_gamma105_bestwb` geometry: 0--2152.6 km, 2,000
uniform cells, a fixed 22 kK outer reservoir, and conduction enabled. The
well-balanced/log-MUSCL/MC3 stack is now intrinsic to `model_column`; the
finite-rate ionization network and radiative cooling remain off. The standard
10,000-step run reaches t=115.835 s and supplies 1,001 finite native frames;
the movie samples 400 at 20 fps.

```bash
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat ISO_H_BASE=0 ISO_DH=2152.6 \
    ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=2000 \
    ./build/chromo_main \
    outputs/model_column/model_column_saha_gamma_t22k_ns2000.txt \
    full no-ionization model_column - 1.0 no-cooling
ANIM_MAX_FRAMES=400 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
    .venv/bin/python util/animate_isentropic.py \
    outputs/model_column/model_column_saha_gamma_t22k_ns2000.txt \
    visualization/model_column/model_column_saha_gamma_t22k_ns2000_evolution.mp4 20
```

The visibly piecewise heat-flux profile in the shorter-domain movie is not an
EOS-table discontinuity. `c7_full_profile()` linearly interpolates temperature
between sparse Model C7 height knots, so `dT/dh` is piecewise constant and the
diagnostic `q=-kappa*dT/dh` changes slope abruptly at those knots. That earlier
run also had conduction disabled, making the panel a diagnostic estimate rather
than an active solver flux. This extended run enables conduction explicitly.

### `model_column/model_column_saha_gamma_t22k_ns2000_pchip_1000s_evolution.mp4`

Option-2 rerun of the extended column after replacing the gamma-table Model C7
temperature initialization with a shape-preserving monotone cubic Hermite
(PCHIP) profile. The historical fixed-gamma/linear-C7 path is unchanged. PCHIP
preserves every tabulated C7 temperature while making both temperature and its
first derivative continuous, so the initial conductive-flux profile is smooth
instead of inheriting steps from piecewise-constant `dT/dh`. The active-
conduction run reaches exactly 1,000 s in 85,892 steps and writes 431 native
frames; the movie samples 400 frames at 20 fps.

```bash
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=2000 \
    CHROMO_T_END=1000 CHROMO_FRAME_STRIDE=200 \
    ./build/chromo_main \
    outputs/model_column/model_column_saha_gamma_t22k_ns2000_pchip_1000s.txt \
    full no-ionization model_column - 9.0 no-cooling
ANIM_MAX_FRAMES=400 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
    .venv/bin/python util/animate_isentropic.py \
    outputs/model_column/model_column_saha_gamma_t22k_ns2000_pchip_1000s.txt \
    visualization/model_column/model_column_saha_gamma_t22k_ns2000_pchip_1000s_evolution.mp4 20
```

### `model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_evolution.mp4`

High-chromosphere/transition-region-only rerun of the same gamma-table,
22 kK conductive-reservoir experiment. The physical domain is exactly
1600--2153 km with 1,000 uniform cells (0.553 km spacing); the run reaches
exactly 1,000 s in 177,896 steps and stores all 357 native frames.

```bash
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=1000 \
    CHROMO_T_END=1000 CHROMO_FRAME_STRIDE=500 \
    ./build/chromo_main \
    outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s.txt \
    full no-ionization model_column - 30.0 no-cooling
ANIM_MAX_FRAMES=400 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
    .venv/bin/python util/animate_isentropic.py \
    outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s.txt \
    visualization/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_evolution.mp4 20
```

Generated with the **pre-fix CENTERED outer ghost** (`p_g0 = p[ns-2] − 2Δs·ρ_top·g`);
kept as the before-side of the upper-BC comparison below. The current code produces
the `_onesidedbc` run instead.

### `model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_onesidedbc_evolution.mp4`

Same run with the reworked outer ghost: a one-sided hydrostatic ladder
(`p_{k+1} = p_k − Δs·½(ρ_k+ρ_{k+1})·g`) anchored on a fixed reservoir
back-pressure captured from the IC top cell, instead of extrapolating across the
boundary cell from `p[ns-2]`. This removes the inverted top-face pressure
gradient, but the mass-flux ripple over 2130--2153 km is **unchanged**
(rel-std 0.101 → 0.105; V[last] −12.4 → −15.9 m/s), so the ripple is not caused
by the ghost pressure anchor. Identical command apart from the output paths.

```bash
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=1000 \
    CHROMO_T_END=1000 CHROMO_FRAME_STRIDE=500 \
    ./build/chromo_main \
    outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_onesidedbc.txt \
    full no-ionization model_column - 30.0 no-cooling
ANIM_MAX_FRAMES=400 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
    .venv/bin/python util/animate_isentropic.py \
    outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_onesidedbc.txt \
    visualization/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_onesidedbc_evolution.mp4 20
```

### Upper-BC hydro-T / conduction-T decoupling (`ISO_HYDRO_T_DECOUPLE`)

Controlled pair on the same 1600--2153 km, `ns=1000` transition-region column. The
**baseline** keeps the outer hydro ghost pinned at the 22 kK wall; the **decoupled**
run zero-gradient-extrapolates the hydro ghost temperature from the live top cell
(`T_hydro_g0 = T_hydro_g1 = T_top`) while Stage-D conduction keeps the identical
22 kK Dirichlet wall via `grid.outer_conduction_temperature`. Everything else --
domain, IC, EOS table, back-pressure `kOuterPRef = 0.0102845 Pa`, Mach cap, cooling,
TRAC, CFL, output cadence -- is identical. Result: the last-cell velocity reversal
disappears (`V_last` -31.6 -> +91.3 m/s) but the 2130--2150 km mass-flux ripple is
unchanged-to-slightly-worse. Full analysis in
`docs/studies/boundaries/upper_bc_hydro_temperature_decoupling_recap.md`.

`upper_bc_hydroT_decouple_500s_top90.png` -- top-90-cell `rho V`, `V` and `T`
profiles of the two runs at 500 s (the figure that shows the ripple is common to
both and only the final point differs).
`upper_bc_hydroT_decouple_500s_metrics.json` -- all twelve requested metrics at
t = 100/300/400/500 s plus the full `q_top(t)`, `E_cond(t)`, `M(t)` series.
`upper_bc_hydroT_decouple_1000s_top90.png` / `..._1000s_metrics.json` -- the same
figure and metrics for the 1000 s pair (t = 700/1000 s), whose baseline reproduces
the archived `_onesidedbc` run to 4--5 significant figures.
`model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_hydroT{base,decoupled}_evolution.mp4`
-- the two 1000 s evolution movies, kept separate from the archived runs.

```bash
# the run pair (ISO_HYDRO_T_DECOUPLE=0 baseline / =1 decoupled; 500 s and 1000 s)
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=1000 ISO_HYDRO_T_DECOUPLE=0 \
    CHROMO_T_END=500 CHROMO_FRAME_STRIDE=500 \
    ./build/chromo_main \
    outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_500s_hydroTbase.txt \
    full no-ionization model_column - 30.0 no-cooling \
    > outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_500s_hydroTbase.console.log 2>&1
# ... same with ISO_HYDRO_T_DECOUPLE=1 -> ..._500s_hydroTdecoupled.{txt,console.log}
# ... same with CHROMO_T_END=1000       -> ..._1000s_hydroTbase / ..._1000s_hydroTdecoupled

# metrics table + JSON + top-region comparison figure
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/upper_bc_decouple_diag.py \
    --run baseline=outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_500s_hydroTbase.txt \
    --run decoupled=outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_500s_hydroTdecoupled.txt \
    --decoupled decoupled --p-back baseline=0.0102845 --p-back decoupled=0.0102845 \
    --times 100 300 400 500 \
    --json outputs/model_column/upper_bc_hydroT_decouple_500s_metrics.json \
    --plot visualization/model_column/upper_bc_hydroT_decouple_500s_top90.png --plot-time 500
# ... same with the _1000s_ runs, --times 700 1000, and the _1000s_ json/png names

# evolution movies (1000 s pair)
ANIM_MAX_FRAMES=400 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
    .venv/bin/python util/animate_isentropic.py \
    outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_hydroTbase.txt \
    visualization/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_hydroTbase_evolution.mp4 20
ANIM_MAX_FRAMES=400 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
    .venv/bin/python util/animate_isentropic.py \
    outputs/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_hydroTdecoupled.txt \
    visualization/model_column/model_column_saha_gamma_h1600_2153_ns1000_pchip_1000s_hydroTdecoupled_evolution.mp4 20
```

### Top-region face mass-flux decomposition (`CHROMO_FACE_FLUX_DIAG`)

Answers whether the cell-centred `rho V` ripple below the top boundary also exists
in the finite-volume Rusanov face mass flux. It does not -- see
`docs/studies/numerics/top_ripple_face_flux_diagnosis.md`. All numbers come from a read-only
copy-out of the production `rhs_explicit_mixture` arrays (`.faceflux` sidecar),
never from a Python reimplementation of the reconstruction.

Metric definitions (also in the script docstring): index `i` is the UPPER face
`i+1/2` of cell `i`; `f_central = 0.5[(rho v)_L + (rho v)_R]`;
`f_diff = -0.5 a (rho_R - rho_L)`; `f_total = f_central + f_diff` (read from the
production flux, not recomputed); `f_ref` = `f_total` at t=0, i.e. the frozen
`eq_wb` reference flux; **`f_eff = f_total - f_ref` is the flux the update actually
sees, because `eq_wb` subtracts the reference residual from the continuity rows**
(HISTORICAL: true for the runs catalogued here, which predate the reference-free
release; `eq_wb` has since been retired, so for new runs `f_eff` is a diagnostic
baseline only -- see `docs/studies/numerics/reference_free_release_recap.md`);
`R = -(F[i]-F[i-1])/ds`; `roughness(x) = mean|D2 x| / mean|x|` (same definition as
`util/upper_bc_decouple_diag.py`).

`faceflux_500s_decomposition.png` / `faceflux_1000s_decomposition.png` -- four
panels over the top 90 cells: cell `rho V`, face total vs central flux, face
diffusive flux, and `f_eff`. The bottom panel is the result: `f_eff` is smooth and
monotone straight through the ripple region in both runs.
`faceflux_500s_limiter_compare.png` -- MC3(beta=2) vs minmod vs first order.
`faceflux_500s_cfl_compare.png` -- CFL 0.25 / 0.125 / 0.0625.
`faceflux_{500,1000}s_report.txt` -- the printed metric tables plus the per-peak
state dumps. `faceflux_*_metrics.json` -- the same as JSON.

```bash
# the four base runs (baseline/decoupled x 500/1000 s). CHROMO_FRAME_STRIDE is huge
# so only the first and last snapshots land on disk; the .faceflux sidecar captures
# step 0 (the eq_wb reference) and the final step.
for D in 0 1; do for T in 500 1000; do
  [ $D = 0 ] && L=base || L=decoupled
  O=outputs/model_column/faceflux_${T}s_hydroT${L}
  GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=1000 ISO_HYDRO_T_DECOUPLE=$D \
    CHROMO_T_END=$T CHROMO_FRAME_STRIDE=1000000 \
    CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=150 \
    ./build/chromo_main $O.txt full no-ionization model_column - 30.0 no-cooling \
    > $O.console.log 2>&1
done; done

# scheme-comparison legs: same as the baseline above plus one extra variable each.
# ISO_LIMITER is a DIAGNOSTIC-ONLY override; unset (production) = MC3 beta=2.
#   ... ISO_LIMITER=minmod   -> outputs/model_column/faceflux_500s_lim_minmod.txt
#   ... ISO_LIMITER=first    -> outputs/model_column/faceflux_500s_lim_first.txt
#   ... CHROMO_CFL=0.125     -> outputs/model_column/faceflux_500s_cfl0125.txt
#   ... CHROMO_CFL=0.0625    -> outputs/model_column/faceflux_500s_cfl00625.txt

# metrics + per-peak dumps + figure (repeat with _1000s_ names and --times equivalents)
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/face_flux_diag.py \
    --run baseline=outputs/model_column/faceflux_500s_hydroTbase.txt.faceflux \
    --run decoupled=outputs/model_column/faceflux_500s_hydroTdecoupled.txt.faceflux \
    --json outputs/model_column/faceflux_500s_metrics.json \
    --plot visualization/model_column/faceflux_500s_decomposition.png \
    > outputs/model_column/faceflux_500s_report.txt 2>&1

MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/face_flux_diag.py \
    --run mc3=outputs/model_column/faceflux_500s_hydroTbase.txt.faceflux \
    --run minmod=outputs/model_column/faceflux_500s_lim_minmod.txt.faceflux \
    --run first=outputs/model_column/faceflux_500s_lim_first.txt.faceflux \
    --json outputs/model_column/faceflux_500s_limiter_metrics.json \
    --plot visualization/model_column/faceflux_500s_limiter_compare.png

MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/face_flux_diag.py \
    --run cfl0.25=outputs/model_column/faceflux_500s_hydroTbase.txt.faceflux \
    --run cfl0.125=outputs/model_column/faceflux_500s_cfl0125.txt.faceflux \
    --run cfl0.0625=outputs/model_column/faceflux_500s_cfl00625.txt.faceflux \
    --json outputs/model_column/faceflux_500s_cfl_metrics.json \
    --plot visualization/model_column/faceflux_500s_cfl_compare.png
```

### Resolution convergence of the top-region ripple (`util/resolution_scan_diag.py`)

Follow-up to the section above: does the cell-centred `rho V` ripple in
2130-2150 km converge as `Delta h -> 0`, while `f_eff = f_total - f_ref` stays
smooth? See `docs/studies/validation/resolution_convergence_scan_recap.md`. Same physical domain
(1600-2153 km), same C7/pchip IC, same gamma/Saha table, same decoupled-hydro-T
upper BC, same 22 000 K conduction wall, cooling off, TRAC off, MC3 beta=2,
CFL 0.25, Mach cap 0.1, production `numerical_diffusivity = 2e3 * Delta h`
(deliberately left to scale with the mesh). `ns` is the only physics-affecting
variable; the diagnostic strides are set in STEPS so that all three runs capture
at the same ~6 s physical cadence and the same 83 km physical extent.

Answer: `A_M` (relstd of `rho V`) falls 0.150 -> 0.137 -> 0.067 and the
face-to-cell mismatch `A_mismatch` tracks it almost exactly, while `f_eff` stays
smooth at every resolution -- Result A. But `mean(f_eff)` falls ~linearly with
`Delta h` (3.69 -> 2.40 -> 1.84e-9 at 1000 s) with no sign of flattening --
Result D, concurrently. The mean evaporation rate is NOT grid-converged.

`resconv_profiles_{500,1000}s.png` -- four panels vs PHYSICAL height: cell
`rho V`, the same normalised by the window-mean `f_eff` (so the ripple
convergence is visible independently of the mean-level shift), `f_eff`, and
`delta F_diff`. `resconv_convergence_loglog.png` -- `A_M`/`Q_M`/`A_diff`/
`A_mismatch`/`A_eff` vs `Delta h` with `Delta h` and `Delta h^2` guides.
`resconv_mean_flux.png` -- mean mass flux and `R_eff` vs resolution.
`resconv_time_compare.png` -- 500 s vs 1000 s trends side by side.
`resconv_report.txt` / `resconv_metrics.json` -- full tables, Richardson
extrapolation of the mean flux, extremum positions and cross-resolution peak
matching.

```bash
# the three runs. time_mult=60.0 (argv[6]) ONLY raises step_cap (CHROMO_T_END
# overrides total_time); ns=2000 needs >300k steps. Verified byte-identical to
# the same runs at time_mult=30.0 for ns=500/1000.
#   ns=500  -> CHROMO_FRAME_STRIDE=2500  CHROMO_FACE_FLUX_TOP=75  CHROMO_FACE_FLUX_STRIDE=500
#   ns=1000 -> CHROMO_FRAME_STRIDE=5000  CHROMO_FACE_FLUX_TOP=150 CHROMO_FACE_FLUX_STRIDE=1000
#   ns=2000 -> CHROMO_FRAME_STRIDE=10000 CHROMO_FACE_FLUX_TOP=300 CHROMO_FACE_FLUX_STRIDE=2000
for NS in 500 1000 2000; do
  O=outputs/model_column/resconv_ns${NS}_dec_1000s
  GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=$NS ISO_HYDRO_T_DECOUPLE=1 \
    CHROMO_T_END=1000 CHROMO_FRAME_STRIDE=$((NS*5)) \
    CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=$((NS*150/1000)) \
    CHROMO_FACE_FLUX_STRIDE=$NS \
    ./build/chromo_main $O.txt full no-ionization model_column - 60.0 no-cooling \
    > $O.console.log 2>&1
done

# metrics, report, JSON and all five figures in one call
MPLCONFIGDIR=/tmp/chromosphere2026-mpl python util/resolution_scan_diag.py \
    --run 500=outputs/model_column/resconv_ns500_dec_1000s.txt \
    --run 1000=outputs/model_column/resconv_ns1000_dec_1000s.txt \
    --run 2000=outputs/model_column/resconv_ns2000_dec_1000s.txt \
    --times 500 1000 --window 2130 2150 \
    --json outputs/model_column/resconv_metrics.json \
    --report outputs/model_column/resconv_report.txt \
    --plot-dir visualization/model_column --plot-prefix resconv
```

### Is grid-scaled numerical conduction why the evaporation flux does not converge? (`util/numerical_conduction_diag.py`)

Follow-up to the section above, targeting Result D only. `numerical_diffusivity =
2000 * Delta h` enters the conduction solve as `kappa_num = numerical_diffusivity *
C_V`, so it scales with the mesh. The new default-off `.outercond` sidecar
(`CHROMO_OUTER_COND_DIAG=1`) captures the OUTER face of the FINAL CONVERGED
conduction Newton iteration, so `q_phys` / `q_num` / `q_total` are the production
face quantities, not a cell-centred proxy. The new diagnostic override
`ISO_NUMERICAL_DIFFUSIVITY_MULT=0` zeroes ONLY that coefficient (physical
conductivity, hydro, boundaries, EOS, limiter, CFL, well-balancing untouched);
unset = production. 2x2 = (ns 1000, 2000) x (production, `D_num=0`). See
`docs/studies/conduction/numerical_conduction_convergence_recap.md`.

Answer: `q_num/q_total` = 0.577 (ns=1000) and 0.385 (ns=2000), and with `D_num=0`
the outer flux collapses ~7x and `mean(F_eff)` ~6x — so most of the production
evaporation drive is the artificial conductivity. But `D_num=0` is still NOT
converged (`mean(F_eff)` +34 % / +16 % from ns=1000 to 2000, opposite sign to
production's -33 % / -25 %), and `D_num=0` at ns=2000 aborts at t = 921.5 s. Only
a partial explanation. The `dnum0_ns2000` figures/report are therefore quoted at
500 s and ~900 s, never 1000 s.

```bash
# the four runs. ISO_NUMERICAL_DIFFUSIVITY_MULT=1 is identical to unset (verified:
# prod ns=1000/2000 reproduce the resconv runs' mean(F_eff) to all printed digits).
# ncond_ns2000_dnum0 ends at t = 921.5 s (uncaught std::domain_error, see the recap).
for NS in 1000 2000; do for MULT in 1 0; do
  TAG=$([ $MULT = 1 ] && echo prod || echo dnum0)
  O=outputs/model_column/ncond_ns${NS}_${TAG}_1000s
  GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=$NS ISO_HYDRO_T_DECOUPLE=1 \
    ISO_NUMERICAL_DIFFUSIVITY_MULT=$MULT \
    CHROMO_T_END=1000 CHROMO_FRAME_STRIDE=$((NS*5)) \
    CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=$((NS*150/1000)) \
    CHROMO_FACE_FLUX_STRIDE=$NS \
    CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=$NS \
    ./build/chromo_main $O.txt full no-ionization model_column - 60.0 no-cooling \
    > $O.console.log 2>&1
done; done

# one table + one figure
MPLCONFIGDIR=/tmp/chromosphere2026-mpl python util/numerical_conduction_diag.py \
    --run prod_ns1000=outputs/model_column/ncond_ns1000_prod_1000s.txt \
    --run prod_ns2000=outputs/model_column/ncond_ns2000_prod_1000s.txt \
    --run dnum0_ns1000=outputs/model_column/ncond_ns1000_dnum0_1000s.txt \
    --run dnum0_ns2000=outputs/model_column/ncond_ns2000_dnum0_1000s.txt \
    --times 500 900 --window 2130 2150 \
    --report outputs/model_column/ncond_report.txt \
    --json outputs/model_column/ncond_metrics.json \
    --plot visualization/model_column/ncond_numerical_conduction.png
```

`ncond_numerical_conduction.png` -- left: `q_total` at the outer face vs `N` for
both branches (plus the `q_num` component of the production branch); right:
`mean(F_eff)` over 2130-2150 km vs `N`. Faint = 500 s, solid = ~900 s.

### Fixed gamma=1.05 vs gamma-table physical conductive flux at 1000 s

Matched-height, matched-time comparison using the shared cell-centred
`q_physical = -(kappa_e+kappa_n) dT/ds` diagnostic. The figure excludes TRAC,
numerical diffusivity, and imposed boundary flux; the JSON records below-500 km
maxima for the temperature gradient, densities, conductivity components, and flux.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/plot_physical_flux_compare.py --fixed outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt --gamma outputs/model_column/model_column_saha_gamma_t22k_ns2000_pchip_1000s.txt --time 1000 --fixed-gamma 1.05 --out visualization/model_column/iso_t22k_fixed105_vs_gamma_physical_flux_t1000.png --metrics-out outputs/model_column/iso_t22k_fixed105_vs_gamma_physical_flux_t1000_metrics.json
```

### Outer (upper-TR) static local refinement: does it reduce the top-region `rho V` ripple?

`docs/studies/numerics/static_local_refinement.md` (§ `outer` profile) + `docs/studies/numerics/outer_tr_refinement_recap.md`.
Same physical domain 1600-2153 km, same coarse-equivalent `ISO_NS=1000`, same
gamma-table / decoupled-hydro-T / conduction-wall / cooling-off configuration as the
resolution-convergence and numerical-conduction stages -- only the mesh changes. The
fine band is 2100-2153 km with a graded transition just below it (2080-2100 km).
`CHROMO_*_STRIDE` are scaled by the refinement factor so all runs capture at the same
PHYSICAL cadence (dt drops with the smallest local cell).

```bash
# U (uniform), R2 (outer refinement x2), R4 (outer refinement x4)
for TAG in U R2 R4; do
  case $TAG in U) F=1; S=1;; R2) F=2; S=2;; R4) F=4; S=4;; esac
  O=outputs/model_column/outref_${TAG}_500s
  env GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
    ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
    ISO_HEAT_FLUX=1 ISO_NS=1000 ISO_HYDRO_T_DECOUPLE=1 \
    ISO_NUMERICAL_DIFFUSIVITY_MULT=1 \
    $( [ $TAG != U ] && echo "ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=$F ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20" ) \
    CHROMO_T_END=500 CHROMO_FRAME_STRIDE=$((5000*S)) \
    CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=320 CHROMO_FACE_FLUX_STRIDE=$((1000*S)) \
    CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=$((1000*S)) \
    ./build/chromo_main $O.txt full no-ionization model_column - 60.0 no-cooling \
    > $O.console.log 2>&1
done

# the one reduced-numerical-diffusivity test, on R4 only
O=outputs/model_column/outref_R4_mult025_500s
env GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
  ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 \
  ISO_HEAT_FLUX=1 ISO_NS=1000 ISO_HYDRO_T_DECOUPLE=1 \
  ISO_NUMERICAL_DIFFUSIVITY_MULT=0.25 \
  ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20 \
  CHROMO_T_END=500 CHROMO_FRAME_STRIDE=20000 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=320 CHROMO_FACE_FLUX_STRIDE=4000 \
  CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=4000 \
  ./build/chromo_main $O.txt full no-ionization model_column - 60.0 no-cooling \
  > $O.console.log 2>&1

# one table + one figure (metrics imported from resolution_scan_diag / numerical_conduction_diag)
MPLCONFIGDIR=/tmp/chromosphere2026-mpl python util/outer_refine_diag.py \
    --run U=outputs/model_column/outref_U_500s.txt \
    --run R2=outputs/model_column/outref_R2_500s.txt \
    --run R4=outputs/model_column/outref_R4_500s.txt \
    --run R4_mult025=outputs/model_column/outref_R4_mult025_500s.txt \
    --time 500 --window 2130 2150 \
    --report outputs/model_column/outref_report.txt \
    --json outputs/model_column/outref_metrics.json \
    --plot visualization/model_column/outref_profiles_500s.png
```

`outref_profiles_500s.png` -- top: cell width vs physical height (log `Delta s`,
2000-2153 km) showing the coarse region, the graded transition and the fine band
reaching the outer face; bottom: cell-centred `rho V` at t = 500 s for uniform / x2 /
x4, with the 2130-2150 km metric window shaded.

`outref_R4_mult025_500s_evolution.mp4` -- 6-panel time evolution (T, V, physical
`q_par`, rho, p, rho*V vs height) of the latest outer-refined run: x4 upper-TR
refinement with `ISO_NUMERICAL_DIFFUSIVITY_MULT=0.25`. Only 19 frames, because that
run used `CHROMO_FRAME_STRIDE=20000`.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl python util/animate_isentropic.py \
    outputs/model_column/outref_R4_mult025_500s.txt \
    visualization/model_column/outref_R4_mult025_500s_evolution.mp4 6
```
### Face-local numerical-conduction validation (U and R4-local, 500 s)

No new visualization is required for this numerical-method change. The durable
artifacts are `outputs/model_column/outref_local_{U,R4}_500s.txt` and their
`.gamma_diag`, `.faceflux`, `.outercond`, and `.console.log` sidecars, plus
`outputs/model_column/outref_local_report.txt` and `outref_local_metrics.json`.

```bash
rtk env GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
  ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=1000 \
  ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_FACTOR=1 ISO_NUMERICAL_DIFFUSIVITY_MULT=1 \
  CHROMO_T_END=500 CHROMO_FRAME_STRIDE=5000 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=320 CHROMO_FACE_FLUX_STRIDE=1000 \
  CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=1000 \
  ./build/chromo_main outputs/model_column/outref_local_U_500s.txt \
  full no-ionization model_column - 60.0 no-cooling \
  > outputs/model_column/outref_local_U_500s.console.log 2>&1

rtk env GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
  ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=1000 \
  ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 \
  ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20 \
  ISO_NUMERICAL_DIFFUSIVITY_MULT=1 CHROMO_T_END=500 CHROMO_FRAME_STRIDE=20000 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=320 CHROMO_FACE_FLUX_STRIDE=4000 \
  CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=4000 \
  ./build/chromo_main outputs/model_column/outref_R4_local_500s.txt \
  full no-ionization model_column - 60.0 no-cooling \
  > outputs/model_column/outref_R4_local_500s.console.log 2>&1

rtk /opt/miniconda3/bin/python util/outer_refine_diag.py \
  --run U=outputs/model_column/outref_local_U_500s.txt \
  --run R4_local=outputs/model_column/outref_R4_local_500s.txt \
  --time 500 --window 2130 2150 \
  --report outputs/model_column/outref_local_report.txt \
  --json outputs/model_column/outref_local_metrics.json
```

The uniform result exactly reproduces the saved pre-change baseline. R4-local
uses `chi_num=2.7604e5 m2/s` at the 138 m outer face, gives
`q_num/q_total=0.2308`, `mean(F_eff)=1.057e-9 kg m-2 s-1`, `A_M=0.02954`, and
`Q_M=0.00684`, and remains stable through 500 s.

### Higher-resolution R4-local N=2000 convergence check (500 s)

No visualization was generated. Output:
`outputs/model_column/outref_R4_local_N2000_500s.txt` with `.gamma_diag`,
`.faceflux`, `.outercond`, and `.console.log` sidecars. The `120.0` positional
time multiplier only raises the driver step cap; `CHROMO_T_END=500` fixes the
physical stop time.

```bash
rtk env GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
  ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=2000 \
  ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 \
  ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20 \
  ISO_NUMERICAL_DIFFUSIVITY_MULT=1 CHROMO_T_END=500 CHROMO_FRAME_STRIDE=40000 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=640 CHROMO_FACE_FLUX_STRIDE=8000 \
  CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=8000 \
  ./build/chromo_main outputs/model_column/outref_R4_local_N2000_500s.txt \
  full no-ionization model_column - 120.0 no-cooling \
  > outputs/model_column/outref_R4_local_N2000_500s.console.log 2>&1

rtk /opt/miniconda3/bin/python util/outer_refine_diag.py \
  --run R4_local_N1000=outputs/model_column/outref_R4_local_500s.txt \
  --run R4_local_N2000=outputs/model_column/outref_R4_local_N2000_500s.txt \
  --time 500 --window 2130 2150 \
  --report outputs/model_column/outref_local_N1000_N2000_report.txt \
  --json outputs/model_column/outref_local_N1000_N2000_metrics.json
```

Result: N=2000 produces 2638 cells with 276.45 m coarse and 69.10 m
outer-fine spacing, completes 712,320 steps, and gives
`mean(F_eff)=6.1003e-10 kg m-2 s-1`, a 42.3% decrease from N=1000. This is a
stable but still under-resolved result, not a convergence claim.

### Third convergence point: R4-local N=4000 (OpenMP, 500 s)

No visualization was generated; the metrics answer the question. Output:
`outputs/model_column/outref_R4_local_N4000_500s.txt` with `.gamma_diag`,
`.faceflux`, `.outercond`, and `.console.log` sidecars. Run with the committed
OpenMP build at the workspace-default 12 threads (already established
byte-identical to serial, so no equivalence run was repeated). The `240.0`
positional multiplier only raises the driver step cap; `CHROMO_T_END=500` fixes
the physical stop time. Diagnostic stride and capture width are scaled so the
physical diagnostic coverage and sampling cadence match N=1000/N=2000.

```bash
cmake -S . -B build_omp -DCMAKE_BUILD_TYPE=Release \
  -DCHROMO_ENABLE_OPENMP=ON -DLIBOMP_ROOT="$(brew --prefix libomp)"
cmake --build build_omp -j4

env OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 OMP_NUM_THREADS=12 \
  ctest --test-dir build_omp --output-on-failure

env OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 OMP_NUM_THREADS=12 \
  GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
  ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=4000 \
  ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 \
  ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20 \
  ISO_NUMERICAL_DIFFUSIVITY_MULT=1 CHROMO_T_END=500 CHROMO_FRAME_STRIDE=80000 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=1280 CHROMO_FACE_FLUX_STRIDE=16000 \
  CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=16000 \
  ./scripts/run_chromo_omp.sh \
  outputs/model_column/outref_R4_local_N4000_500s.txt \
  full no-ionization model_column - 240.0 no-cooling \
  > outputs/model_column/outref_R4_local_N4000_500s.console.log 2>&1

/opt/miniconda3/bin/python util/outer_refine_diag.py \
  --run R4_N1000=outputs/model_column/outref_R4_local_500s.txt \
  --run R4_N2000=outputs/model_column/outref_R4_local_N2000_500s.txt \
  --run R4_N4000=outputs/model_column/outref_R4_local_N4000_500s.txt \
  --time 500 --window 2130 2150 \
  --report outputs/model_column/outref_local_N1000_N2000_N4000_report.txt \
  --json outputs/model_column/outref_local_N1000_N2000_N4000_metrics.json
```

`util/outer_refine_diag.py` gained an `outer_boundary_layer` helper and the
`ds_face`, `T_wall - T_top`, `G_wall` and boundary-layer report rows; every
pre-existing metric is unchanged and reproduces the saved N=1000/N=2000 values.

Result: N=4000 produces 5275 cells (138.2 m coarse, 34.55 m outer-fine),
completes 1,397,913 steps to 500 s, and gives `chi_num=6.9100e4 m2/s`,
`q_num/q_total=0.0685`, and `mean(F_eff)=6.2250e-10 kg m-2 s-1` -- a **+2.04%**
change from N=2000 versus **-42.27%** from N=1000 to N=2000, i.e. the successive
difference shrank 35.8x. Emerging convergence; the limiting evaporation rate
remains provisional and must not be quoted as a physical result.

### 1000 s R4-local N=2000 production column (latest settings)

Full 1000 s extension of the latest `model_column` configuration: 1600--2153 km
domain, gamma-table EOS, decoupled hydro-T upper BC with the 22 kK conduction
wall, x4 face-local outer refinement, production numerical diffusivity, cooling
and the finite-rate ionization network off. 2638 cells (276.45 m coarse,
69.10 m outer-fine), 1,395,699 steps to exactly t = 1000 s on the OpenMP build
at 12 threads (~28 min wall). 350 native frames, all used in the movie.

The `240.0` positional multiplier only raises the driver step cap;
`CHROMO_T_END=1000` fixes the physical stop time.

```bash
env OMP_DYNAMIC=FALSE OMP_MAX_ACTIVE_LEVELS=1 OMP_NUM_THREADS=12 \
  GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
  ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=2000 \
  ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 \
  ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20 \
  ISO_NUMERICAL_DIFFUSIVITY_MULT=1 CHROMO_T_END=1000 CHROMO_FRAME_STRIDE=4000 \
  CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_TOP=640 CHROMO_FACE_FLUX_STRIDE=8000 \
  CHROMO_OUTER_COND_DIAG=1 CHROMO_OUTER_COND_STRIDE=8000 \
  ./scripts/run_chromo_omp.sh \
  outputs/model_column/outref_R4_local_N2000_1000s.txt \
  full no-ionization model_column - 240.0 no-cooling \
  > outputs/model_column/outref_R4_local_N2000_1000s.console.log 2>&1

ANIM_MAX_FRAMES=350 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
  /opt/miniconda3/bin/python util/animate_isentropic.py \
  outputs/model_column/outref_R4_local_N2000_1000s.txt \
  visualization/model_column/outref_R4_local_N2000_1000s_evolution.mp4 20

/opt/miniconda3/bin/python util/outer_refine_diag.py \
  --run R4_N2000_1000s=outputs/model_column/outref_R4_local_N2000_1000s.txt \
  --time 1000 --window 2130 2150 \
  --report outputs/model_column/outref_R4_local_N2000_1000s_report.txt \
  --json outputs/model_column/outref_R4_local_N2000_1000s_metrics.json
```

`outref_R4_local_N2000_1000s_evolution.mp4` -- 6-panel time evolution (T, V,
physical `q_par`, rho, p, rho*V vs height), 350 frames at 20 fps (17.5 s).

`outref_R4_local_N2000_1000s_final_profiles.png` -- final saved 6-panel frame
used by the release paper. It is extracted from the verified movie above, so it
inherits the same run configuration and plotted quantities.

```bash
ffmpeg -y \
  -i visualization/model_column/outref_R4_local_N2000_1000s_evolution.mp4 \
  -vf "select=eq(n\\,349)" -frames:v 1 \
  visualization/model_column/outref_R4_local_N2000_1000s_final_profiles.png
```

Metrics at t = 1000 s over 2130--2150 km: `A_M = 0.0102`, `Q_M = 0.00184`,
`mean(F_eff) = 6.742e-10 kg m-2 s-1`, `q_num/q_total = 0.1289`,
`T_wall - T_top = 92.8 K`, boundary layer 9.67 km / 141 cells. The ripple is
~3x smaller than the same configuration at 500 s (`A_M = 0.0295`), and
`mean(F_eff)` is +10.5% over the 500 s value -- still a provisional number, not
a converged physical evaporation rate.

### Release N=500 mesh and initial temperature

`model_column_N500_mesh_temperature_initial.png` and `.pdf` -- exact cell width and
EOS-aware initial (`t=0`) temperature versus height for the **release**
coarse-equivalent N=500 column. The static x4 outer refinement above 500 km gives
661 cells over 1600--2153 km: coarse 1104.768 m grading to fine 276.042 m, with the
fine region starting at 2100 km. Initial temperature spans 6631.35--22674.60 K.

Provenance is the release smoke run `release_roe_lnp_N500_100s`, which was launched
with no environment overrides, so the mesh and initial state are the scenario's own
release defaults AT THE TIME (mixture Roe, `(ln rho,V,ln p)`, MC3 beta=2; the release flux is now `swmf-godunov` -- confirmed by the
console log's `flux=roe-local reconstruction=lnrho-v-lnp limiter=mc3 beta=2`). The
script independently rebuilds the faces from `scenarios/mesh.cpp` and cross-checks
them against the run header (max face error 4.995 m, within the float32
six-significant-digit header rounding tolerance).

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
  util/plot_model_column_mesh_temperature.py \
  --run outputs/model_column/release_roe_lnp_N500_100s.txt \
  --ns-coarse 500 \
  --output visualization/model_column/model_column_N500_mesh_temperature_initial
```

### Historical N=2000 previous-production mesh and initial temperature

`model_column_ns2000_mesh_temperature_initial.png` and `.pdf` -- exact cell width
and EOS-aware initial (`t=0`) temperature versus height for the historical
coarse-equivalent N=2000 previous-production column. The static x4 outer refinement adds
cells, giving 2638 cells over 1600--2153 km.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
  util/plot_model_column_mesh_temperature.py
```

### `model_column_2000s_evolution.mp4` -- historical 2000 s previous-production run

Full-length (t = 0 -> 2000 s) evolution of the **historical previous-production
`model_column` configuration** (2638 refined cells: `ISO_NS=2000` with the
`outer` R4 local refinement above 500 km, domain h = 1600 -> 2153 km, gamma
table EOS, hydro-T-decoupled upper BC, conduction on, no cooling). Launched
through `scripts/run_chromo_realtime.sh` so it inherits the validated
`CHROMO_CFL=0.50` and the 12-thread OpenMP runtime; `CHROMO_T_END=2000` gives a
true end-time termination and `CHROMO_FRAME_DT=4` keeps the snapshot count at
501 instead of tens of thousands. Standard 2x3 `animate_isentropic.py` panels
(T, V [km/s], physical `q_par`, rho, p, rho*V vs height), 500 frames @ 25 fps
(20 s). Run summary: `termination=end_time`, final step 1,395,700,
mean dt 1.4054e-3 s, wall time 1826.5 s (0.91x real time), no NaN / clamp /
floor activation.

```bash
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=2000 \
ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 \
ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20 ISO_NUMERICAL_DIFFUSIVITY_MULT=1 \
CHROMO_T_END=2000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=4 \
/usr/bin/time -p scripts/run_chromo_realtime.sh \
  outputs/model_column/realtime_2000s_cfl050.txt \
  full no-ionization model_column - 20.0 no-cooling \
  > outputs/model_column/realtime_2000s_cfl050.console.log 2>&1

ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/realtime_2000s_cfl050.txt \
  visualization/model_column/model_column_2000s_evolution.mp4 25
```

### `model_column_4000s_evolution.mp4` -- historical 4000 s previous-production run

Full-length (t = 0 -> 4000 s) continuation of the same historical
previous-production N=2000 configuration documented above. The x4 outer
refinement gives 2638 cells, and `CHROMO_FRAME_DT=8` records 501 synchronized
solution/EOS frames without excessive ASCII output. The run reached exactly
4000 s with `termination=end_time` at step 2,761,034 in 3463.52 s wall time;
the saved fields and EOS sidecar are finite, with no warning, clamp, or floor
activation. Standard 2x3 `animate_isentropic.py` panels (T, V [km/s], physical
`q_par`, rho, p, rho*V vs height), 500 frames @ 25 fps (20 s).

```bash
GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
ISO_H_BASE=1600 ISO_DH=553 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 ISO_NS=2000 \
ISO_HYDRO_T_DECOUPLE=1 ISO_REFINE_PROFILE=outer ISO_REFINE_FACTOR=4 \
ISO_REFINE_S_LO_KM=500 ISO_REFINE_TRANSITION_KM=20 ISO_NUMERICAL_DIFFUSIVITY_MULT=1 \
CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=8 \
/usr/bin/time -p scripts/run_chromo_realtime.sh \
  outputs/model_column/realtime_4000s_cfl050.txt \
  full no-ionization model_column - 20.0 no-cooling \
  > outputs/model_column/realtime_4000s_cfl050.console.log 2>&1

ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/realtime_4000s_cfl050.txt \
  visualization/model_column/model_column_4000s_evolution.mp4 25
```

### `model_column_release_N1000_4000s_evolution.mp4` -- focused long-duration resolution test

Directly comparable 0--4000 s animation for the physical-conduction-only release
model with only the coarse-equivalent resolution overridden from N=500 to N=1000.
The run uses R4 outer refinement (1320 actual cells), the fixed 22,000 K physical-face
conductive boundary, Gamma/Saha EOS, hydro-T decoupling, and CFL=0.50. The animation
uses the same standard 2x3 `animate_isentropic.py` panels and frame cap as the N500
release visualization.

```bash
ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/model_column_release_N1000_4000s.txt \
  visualization/model_column/model_column_release_N1000_4000s_evolution.mp4 25
```

### `acc_lnp_N500_rho_v_p.mp4` — (ρ, V, p) evolution, lnP reconstruction acceptance run

Three-field (density, velocity, total gas pressure) 0–4000 s animation of the
`(ln rho, V, ln p)` reconstruction acceptance run at N=500 (661 actual cells),
**Rusanov** reference solver (the release solver was Roe at the time; it is now the SWMF exact-Riemann Godunov flux), Gamma/Saha EOS, R4 outer refinement, 22,000 K physical-face
conductive boundary, physical conduction only, CFL=0.50. Fields are read from the
`.gamma_diag` sidecar, so they are the EOS-consistent `x_eq`-derived quantities.
Left column is the full 1600–2153 km column; right column zooms the refined
2100–2153 km evaporation region, with the 2130–2150 km release analysis window
shaded in both. Velocity is shaded blue for upflow and red for downflow — it stays
strictly positive for the whole run, corroborating the `n(V<0)=0` acceptance result.

```bash
.venv/bin/python util/animate_column_rvp.py \
  outputs/model_column/acc_lnp_N500.txt \
  visualization/model_column/acc_lnp_N500_rho_v_p.mp4 20
```

### `acc_lnp_rusanov_N500_4000s_evolution.mp4` — lnP reconstruction under the Rusanov reference solver

Standard 2×3 `animate_isentropic.py` panels (T, V, q∥, ρ, p, ρV) for the
`(ln rho, V, ln p)` reconstruction acceptance run at N=500 (661 actual cells),
**Rusanov** reference solver, Gamma/Saha EOS, R4 outer refinement, 22,000 K
physical-face conductive boundary, physical conduction only, CFL=0.50, 0–4000 s.
Deliberately the same script, panels, frame cap, and fps as
`lnp_roe_N500_4000s_evolution.mp4`, so the Rusanov-reference and Roe versions
of the *same* reconstruction are directly comparable.

```bash
ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/acc_lnp_N500.txt \
  visualization/model_column/acc_lnp_rusanov_N500_4000s_evolution.mp4 25
```

### Backfilled entries — Roe-local evolution movies

These three were produced during the Roe/reconstruction study but were never entered
here. Commands below are **reconstructed** from the run files and the standard
release-movie invocation, and reproduce the figures; they were not recorded at creation.
All were run with an explicit `ISO_RIEMANN=roe-local`, which was still a diagnostic
override at the time; Roe later became the release default and then, at the Godunov cutover, a reference override again; `ISO_RIEMANN=roe-local` reproduces these runs exactly.

`lnp_roe_N500_4000s_evolution.mp4` — `(ln rho,V,ln p)` reconstruction + Roe-local, N500, 4000 s.
`roe_N500_4000s_evolution.mp4` / `roe_N1000_4000s_evolution.mp4` — `(ln rho,V,ln T)` + Roe-local.

```bash
ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/lnp_roe_N500_4000s.txt \
  visualization/model_column/lnp_roe_N500_4000s_evolution.mp4 25
ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/roe_N500_4000s.txt \
  visualization/model_column/roe_N500_4000s_evolution.mp4 25
ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/roe_N1000_4000s.txt \
  visualization/model_column/roe_N1000_4000s_evolution.mp4 25
```

### `lnp_roe_N1000_4000s_evolution.mp4` — N=1000 companion to `lnp_roe_N500_4000s_evolution.mp4`

Same solver as `lnp_roe_N500_4000s_evolution.mp4` — mixture Roe flux (the release flux at the time, now a reference override),
`(ln rho, V, ln p)` MUSCL reconstruction with MC3/beta=2, Gamma/Saha EOS, R4 outer
refinement, 22,000 K physical-face conductive boundary, hydro-T decoupling, physical
conduction only, CFL=0.50 — with only the coarse-equivalent resolution raised from
N=500 to N=1000 (1320 actual cells). Roe and lnP were the release defaults when this was run, so the run
needs no `ISO_RIEMANN`/`ISO_RECONSTRUCTION` override. Run reached
`termination=end_time` at 4000 s, step 1,455,143, in 1310.1 s wall (12 threads).
Snapshot cadence 20 s gives 201 frames; standard 2x3 `animate_isentropic.py` panels
(T, V [km/s], physical `q_par`, rho, p, rho*V vs height) at 25 fps, matching the N500
version for direct comparison.

```bash
ISO_NS=1000 \
CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
/usr/bin/time -p scripts/run_chromo_realtime.sh \
  outputs/model_column/lnp_roe_N1000_4000s.txt \
  full no-ionization model_column - 20.0 no-cooling \
  > outputs/model_column/lnp_roe_N1000_4000s.console.log 2>&1

ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/lnp_roe_N1000_4000s.txt \
  visualization/model_column/lnp_roe_N1000_4000s_evolution.mp4 25
```

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
python util/plot_column_sloshing.py \
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

### `model_column/*_mass_flux.mp4` — cell-centred `rho V` vs the face mass flux

Read-only comparison movie for the release finite-volume mass transport:
cell-centred `rho V` against the production face mass flux and its
equilibrium-reference decomposition. Everything is read from the
`<run>.faceflux` sidecar (`CHROMO_FACE_FLUX_DIAG=1`, release solver only), so all
four series come from the SAME `mixture_rhs_explicit` evaluation at the same step
and the panels are synchronous by construction.

Definitions (index `i` = the UPPER face `i+1/2` of cell `i`; see also the
face-flux decomposition section above): `rhoV_cell = rho_cell * v_cell` plotted at
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
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/animate_face_mass_flux.py \
  outputs/model_column/lnp_roe_N500_4000s.txt \
  visualization/model_column/lnp_roe_N500_4000s_mass_flux.mp4 --fps 20

# current release configuration, 100 s smoke capture (11 records)
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/animate_face_mass_flux.py \
  outputs/model_column/release_roe_lnp_N500_100s.txt \
  visualization/model_column/release_roe_lnp_N500_100s_mass_flux.mp4 --fps 4
```

Options: `--zoom LO HI`, `--window LO HI`, `--stride N` / `--max-frames N`
(or `ANIM_STRIDE` / `ANIM_MAX_FRAMES`), `--faceflux PATH` for a sidecar that does
not sit next to the run, `--final-figure PATH` / `--no-final-figure`,
`--no-check-output`. A sidecar whose first record is not `t = 0` is rejected, since
`f_ref` is undefined without it.

---

## 12. Educational / tutorial animations (`visualization/godunov_tutorial/`, `visualization/roe_tutorial/`)

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

## 13. Release Godunov cutover — N500 4000 s evolution and base mass-flux diagnosis

### Mass-flux panel default changed

`util/animate_isentropic.py` now plots the **conservative numerical face mass flux
`f_total`** from the `<run>.faceflux` sidecar in its mass-flux panel, at the face
heights, whenever that sidecar exists. The cell-centred product `rho*V` is only a
reconstruction of the flux: it carries the full cell-centred ripple and is
systematically misleading in the boundary cells, where the ghost closure fixes the
face flux and not the cell-centred product. Set `ANIM_MASS_FLUX=cell` to force the
old panel; without a sidecar the script falls back to `rho*V` and labels the panel
accordingly. **Every movie in sections 1-12 above predates this change and shows
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
.venv/bin/python util/animate_isentropic.py \
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
  util/plot_base_mass_flux_diagnosis.py \
  "release Godunov (current code)=outputs/model_column/lnp_godunov_N500_4000s.txt" \
  "Roe, pre predictor-source fix=outputs/model_column/lnp_roe_N500_4000s.txt" \
  visualization/model_column/godunov_N500_4000s_base_mass_flux_diagnosis.png
```

---

## 14. Storage-precision comparison (HISTORICAL as written; superseded by section 14a)

> **Superseded 2026-08-17.** Double precision is now the RELEASE storage and `float32` is the
> diagnostic build; `CHROMO_STATE_FLOAT64` is retired and CMake hard-errors on it. Everything below
> is kept as the record of what was run at the time, but its "float64 control" leg was an
> incoherent hybrid — `Vec` was promoted and the pack, predictor, RHS, timestep and physical
> constants were not — so **two of its numbers are wrong** (the 2 % spread and the ~3 % evaporation
> shift). Use section 14a. Do not re-run the commands in this section; they will not configure.

### 14 (as originally written)

A matched control for the lower-chromosphere mass-flux gradient of section 13.
Everything is identical to the `lnp_godunov_N500_4000s` release run — physics,
grid, boundary conditions, reconstruction, Godunov flux, CFL, diagnostics, and
the launcher — except that `chromosphere::Vec` is promoted from `float` to
`double` by the default-off CMake option `CHROMO_STATE_FLOAT64`, built into a
separate tree. **Release defaults are untouched and this build must never be
used for a quoted release number.** Note the promotion also covers the static
mesh metric arrays, which share the `Vec` type; they are geometry, computed once,
so this only removes their representation error.

Both legs ran to `termination=end_time` at 4000 s in the **identical 712,314
steps** with 941,679,108 exact face solves and zero fallbacks; 502 s (float32)
against 507 s (float64) wall on 12 threads.

```bash
cmake -S . -B build_omp_f64 -DCMAKE_BUILD_TYPE=Release \
      -DCHROMO_ENABLE_OPENMP=ON -DCHROMO_STATE_FLOAT64=ON
cmake --build build_omp_f64 -j 12

CHROMO_BINARY=$PWD/build_omp_f64/chromo_main \
CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
CHROMO_FACE_FLUX_DIAG=1 CHROMO_FACE_FLUX_STRIDE=3300 \
/usr/bin/time -p scripts/run_chromo_realtime.sh \
  outputs/model_column/lnp_godunov_N500_4000s_f64.txt \
  full no-ionization model_column - 20.0 no-cooling \
  > outputs/model_column/lnp_godunov_N500_4000s_f64.console.log 2>&1
```

### `precision_control_N500_4000s.png`

Six panels: `f_total(h)` at 2000 s and at 4000 s for both legs; a base zoom
comparing `f_total` against cell-centred `rho*V`; the final velocity profile; the
column mass budget — the recorded density change over the run against what each
leg's own time-integrated `f_total` divergence implies it should have been; and
the per-step mass-density increment measured in ULPs of each leg's storage type.
Result:

| | float32 (release) | float64 (control) |
| --- | --- | --- |
| `f_total` base/top, t = 2000 s | 2.682 | **1.008** |
| `f_total` base/top, t = 4000 s | 2.811 | **0.996** |
| 1600–1850 km spread of `f_total`, t = 4000 s | 119 % | **2 %** |
| base velocity, t = 4000 s | 1.392 m/s | **0.506 m/s** |
| top velocity, t = 4000 s | 13.18 m/s | 13.58 m/s |
| cells with per-step mass increment below 0.5 ULP | 619 / 660 | **0 / 660** |
| `f_diff` at face 0, as a fraction of `f_total[0]` | 16.3 % | 15.9 % |
| cell 0 `rho*V` / `f_total[0]` | 1.280 | 1.275 |
| mass budget over 1605–1850 km, predicted / observed | **-215** | **+0.94** |

The broad gradient is **entirely a float32 artifact**: in double precision the
column reaches a genuinely constant mass flux from base to top, which is what
continuity demands of a quasi-steady state on a constant-area tube, and the
column mass budget closes to 6 % instead of failing by a factor of 200 with the
wrong sign. The first-interior-face artifact is **precision-independent** and
survives unchanged, so it is a real scheme/boundary-closure effect. The
evaporation observables at the top of the domain are only mildly affected (~3 %).
Full record: `docs/studies/numerics/float32_precision_control_experiment.md`.

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
  util/plot_precision_control.py \
  outputs/model_column/lnp_godunov_N500_4000s.txt \
  outputs/model_column/lnp_godunov_N500_4000s_f64.txt \
  visualization/model_column/precision_control_N500_4000s.png
```

### `lnp_godunov_N500_4000s_f64_evolution.mp4`

The control run animated with the same script, panels, frame cap and fps as
`lnp_godunov_N500_4000s_evolution.mp4`, so the two movies are directly
comparable. The mass-flux panel is the conservative face flux in both; in the
float64 movie it is flat from the first frame onward instead of decaying with
height.

```bash
ANIM_MAX_FRAMES=500 MPLCONFIGDIR=/tmp/chromosphere2026-mpl \
.venv/bin/python util/animate_isentropic.py \
  outputs/model_column/lnp_godunov_N500_4000s_f64.txt \
  visualization/model_column/lnp_godunov_N500_4000s_f64_evolution.mp4 25
```

---

## 14a. Storage-precision cutover comparison (CURRENT — double is the release)

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

Same six panels as section 14 (`f_total(h)` at 2000 s and 4000 s; base zoom of
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
closes, and no cell is frozen. **Two corrections to section 14, both against
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
  util/plot_precision_control.py \
  outputs/model_column/diag_float32_N500_4000s.txt \
  outputs/model_column/release_double_N500_4000s.txt \
  visualization/model_column/precision_cutover_N500_4000s.png
```

---

## 15. Presentation figure crops (2026-08-17 group meeting)

The deck at `docs/presentations/2026-08-17-group-meeting/` keeps its own figures so it
builds from a fresh clone — `visualization/` is gitignored, so a deck reaching into it
would not build for anyone else. Two provenances:

**Copied unchanged** from `visualization/reference/`:

```bash
cp visualization/reference/c7_ionization_profile_compare.png \
   docs/presentations/2026-08-17-group-meeting/figures/
```

Note that `util/plot_c7_ionization_profile.py --mode compare` writes to
`visualization/eos_gamma/`, and the `visualization/reference/` copy was placed by hand;
regenerating it also needs `outputs/eos_gamma/gamma_hydrogen_excitation.dat`, which is
not in the tree (`bash util/eos/build_and_run.sh excitation`).

**Cropped** from the six-panel comparison of section 14 — one panel per slide, because the
whole figure is unreadable at slide size:

```bash
.venv/bin/python - <<'PY'
from PIL import Image
src = Image.open("visualization/model_column/precision_control_N500_4000s.png")
W, H = src.size
top = int(0.045 * H); mid = top + (H - top) // 2; pad = 34   # 0.045 = the suptitle band
rows = [(top, mid + pad), (mid + pad, H)]
cols = [(0, W // 3), (W // 3, 2 * W // 3), (2 * W // 3, W)]
out = "docs/presentations/2026-08-17-group-meeting/figures/"
for name, (r, c) in {"precision_ftotal_4000s": (0, 1),   # (b) f_total at 4000 s
                     "precision_base_zoom":    (0, 2),   # (c) base zoom
                     "precision_mass_budget":  (1, 1),   # (e) column mass budget
                     "precision_ulp_mechanism":(1, 2)}.items():   # (f) the ULP mechanism
    src.crop((cols[c][0], rows[r][0], cols[c][1], rows[r][1])).save(out + name + ".png")
PY
```

Only `precision_ftotal_4000s.png` appears in the deck — the precision result is a single
aside slide — and the other three are backup-slide assets. See
`docs/presentations/2026-08-17-group-meeting/figures/README.md`.

### `presentation.pdf`

```bash
cd docs/presentations/2026-08-17-group-meeting
latexmk -pdf -output-directory=build presentation.tex
```
