# Visualization catalog — two-temperature column (`model_column_2t`, experimental)

The experimental `T_e != T_i` study of the release column, and the equilibration-time figure that motivates it. Outputs land in `visualization/model_column/` and `visualization/model_column_2t/`. See [README.md](README.md) for conventions.

---

## Electron–ion temperature-equilibration time (`util/tau_ei_equilibration.py`)

### `tau_ei_equilibration.png`

Measured (not estimated) `tau_eq^{ei}(h)` for the release column, against the timescales it has to beat for `T_e = T_i` to be a defensible closure. Post-processing only — no new run. The source is the canonical release baseline `outputs/model_column/twall_sweep/twall_22kK_cfl0.50.txt` (N=500/R4, 661 cells, CFL 0.50, double-precision build, physical conduction only, 22 kK outer conductive reservoir, 4000 s, `termination=end_time`, 0 Godunov fallbacks); `T`, `x`, `n_e` and `n_HI` are read straight out of its `.gamma_diag` sidecar, so the ionization state is the run's own Saha closure and nothing is re-derived.

Top panel: `tau_eq^{ei}` from the NRL Plasma Formulary rate `nu_eps = 3.2e-9 Z^2 n_i lnLambda / (mu T_e^{3/2})` with a **per-cell** Coulomb logarithm on the NRL branch that actually applies (branch 1, `23 - ln(n_e^{1/2} Z T_e^{-3/2})`, everywhere — the column never exceeds 10 Z^2 eV), plus the electron–neutral and ion–neutral energy-relaxation estimates and the three competing timescales: the solver `dt` parsed from the console log, the cell sound crossing `ds/c_s` at the release frozen sound speed, and the cell conduction time `ds^2/chi` from the run's own `kappa_physical`. `t = 0` dashed, `t = 4000 s` solid. Bottom panel: the actual `n_e(h)` against the 5e9–1e11 cm^-3 hand-estimate band, with `x(h)` on the right axis.

```bash
.venv/bin/python util/tau_ei_equilibration.py
```

Also writes `outputs/model_column/tau_ei_equilibration.json` (every tabulated scalar). Override the run with `--run`, the sampled heights with `--heights`, the cross-sections with `--sigma-en` / `--sigma-in`, and the labelled 1 MK **extrapolation** row with `--tr-ne` / `--tr-T`.

---

## EXPERIMENTAL two-temperature column (`util/animate_two_temp.py`, `visualization/model_column_2t/`)

One script, two modes (movie by default, `--snapshot` for a publication PNG + PDF), both drawing the **same twelve panels** off the `<run>.twotemp` sidecar written by the `model_column_2t` scenario. Panels: **a,b** `T_e` and `T_i` (plus the release single-temperature `T` when `--ref` is given), full column and zoomed on the R4-refined top; **c** the decoupling `T_e - T_i` [K] in the zoom against the sidecar's own ten-significant-digit print resolution; **d** `|T_e - T_i|/T_i` over the whole column, log, with that same floor shaded; **e** `|Q_ei|` against the exchange a print-floor temperature difference alone would produce; **f** signed `Q_ei` in the zoom, symlog with the print-resolution band as the linear band; **g** `tau_eq` against the cell acoustic crossing `ds/c_s` (the timescale the exchange has to beat); **h** the two conduction channels `kappa_e`, `kappa_i`; **i-l** hydro context `rho`, `V`, `x`/`p_e/p`, and the `p`/`p_i`/`p_e` split.

Panels **c, d, e, f, k** also carry a dashed marker at the single physical `T_e = T_i` sign change, taken as the zero crossing between the most negative and the most positive cell and suppressed entirely unless the profile is an order of magnitude above the print floor (so print-floor flicker never draws one). In the 4000 s runs that marker lands on the **ionization front at 2137.8 km**: below it the electron pool is paying for ionization and `T_e < T_i` by up to 1.1e-2 K, above it the electrons are conductively heated from the 22 kK wall and hand energy to the heavies, `T_e > T_i` by up to 7.1e-2 K. Drawing it on panel k as well is what makes the coincidence with `x -> 1` visible.

The domain is 1600-2153 km, so the height axis is **plain linear** — the split/compressed-axis convention applies only to runs that reach into the corona. The refined top is instead given its own zoom panels and is shaded on every full-column panel. `--ref` matches the release `.gamma_diag` frames to the 2T frames by **nearest time** (index-for-index at the shared 20 s cadence, worst offset 0.01 s; a warning is printed if it exceeds 1 s) and refuses to overlay if the two meshes differ. The `.gamma_diag` velocity column is named `v_cm` for historical reasons but holds SI m/s, so no conversion is applied. Frames are rendered in parallel through `util/_anim_parallel.py`. A trailing partially written frame is dropped, so both modes are safe to run against a file a live simulation is still appending to.

**The release reference is `outputs/model_column/release_N500_4000s_2tref.txt.gamma_diag`** — a matched run with the current binary at identical mesh, cadence, CFL and end time (top velocity 9.50 m/s at 4000 s). Do **not** use `outputs/model_column/lnp_godunov_N500_4000s_f64.*`: that is the retired, withdrawn incoherent `CHROMO_STATE_FLOAT64` hybrid build (top velocity 13.58 m/s) and overlaying it manufactures a ~30 % disagreement that does not exist.

### `twotemp_N500_4000s_Ti_neumann.mp4`, `twotemp_N500_4000s_Ti_dirichlet.mp4`

201-frame (20 s cadence, 0-4000 s, 20 fps) evolution of the two production runs, both `termination=end_time`, zero Godunov fallbacks, zero clamps — outer `T_i` Neumann (the solver's documented default: the fully ionized plasma above the domain has no heavy-particle conductive channel to supply) and the outer `T_i` Dirichlet control — each with the corrected release single-temperature column overlaid.

**The two movies are supposed to look identical.** After the outer heavy-channel wall fix, the two legs agree to ~1e-9 relative in every field, because at 22 kK the gas is fully ionized and `kappa_i` collapses to ~1.3e-8 of `kappa_e`, leaving no heavy conductive flux for the BC to control. That is the result, not a rendering error.

```bash
for BC in neumann dirichlet; do
  .venv/bin/python util/animate_two_temp.py \
    "outputs/model_column_2t/twotemp_N500_4000s_Ti_$BC.txt.twotemp" \
    --ref outputs/model_column/release_N500_4000s_2tref.txt.gamma_diag \
    --fps 20
done
```

Without `--out` the movie defaults to `visualization/model_column_2t/<run label>.mp4`, which is exactly the two filenames above.

### `twotemp_N500_4000s_Ti_neumann_t4000s.png` / `.pdf` (and the `_dirichlet` pair)

Publication figure: the final (t = 4000 s) frame of the same twelve panels at 300 dpi, plus a vector PDF. The header line carries the headline numbers — `max|T_e - T_i| = 0.0715 K` at 2138.4 km, `max|T_e - T_i|/T_i = 9.01e-6`, 645/661 cells above print resolution, and the `T_e = T_i` crossing at 2137.8 km.

```bash
for BC in neumann dirichlet; do
  .venv/bin/python util/animate_two_temp.py \
    "outputs/model_column_2t/twotemp_N500_4000s_Ti_$BC.txt.twotemp" \
    --ref outputs/model_column/release_N500_4000s_2tref.txt.gamma_diag \
    --snapshot
done
```

Without `--out` the snapshot defaults to `visualization/model_column_2t/<run label>_t<time>s.png` (plus the matching `.pdf`), i.e. the two filenames above. `--snapshot` takes the last frame by default; `--time <seconds>` picks the nearest frame to another time. Other options: `--fps`, `--stride N`, `--max-frames N`, `--zoom LO HI` (default `2080 2153`).

### `twotemp_transient_5s_t0p0222s.png` / `.pdf` — the study's headline frame

**The peak decoupling of the entire study is a startup transient, not a steady state.** From the dense 0.02 s-cadence run `outputs/model_column_2t/twotemp_transient_5s.txt.twotemp` (251 frames, 0-5 s, outer `T_i` Neumann), the frame at t = 0.0222 s has the top cell at `T_e = 22339.3 K`, **63.4 K BELOW** `T_i = 22402.8 K` (`max|T_e - T_i|/T_i = 2.83e-3`, ~4600x the 4000 s steady value). It decays roughly as 1/t: 43.1 K at 0.044 s, 18.0 K at 0.10 s, 1.7 K at 0.5 s. Panel c shows it as a single-cell downward spike at the very top; panel d (log) carries the decay.

```bash
.venv/bin/python util/animate_two_temp.py \
    outputs/model_column_2t/twotemp_transient_5s.txt.twotemp --snapshot --time 0.022
```

`--time` selects the nearest frame, and sub-second frames get a `t0p0222s`-style filename tag rather than being rounded to `t0s`. No `--ref` overlay: the release reference has no frames inside the first 20 s.

### `twotemp_transient_5s.mp4`

The same 251 frames as a 25 fps movie, so the 1/t collapse of the 63 K spike is watchable. Global y-limits are set by the peak frame, which is the point.

```bash
.venv/bin/python util/animate_two_temp.py \
    outputs/model_column_2t/twotemp_transient_5s.txt.twotemp --fps 25
```

### `smoke5s.mp4`, `smoke5s_t5s.png` / `.pdf`

Development/validation artifacts from the 5 s, 6-frame smoke run (`outputs/model_column_2t/smoke5s.txt.twotemp`, outer `T_i` Neumann). Kept because they are the cheap end-to-end check of both modes: 6 frames, 661 cells, `max|T_e - T_i| = 0.447 K` at 2150.5 km, `max|T_e - T_i|/T_i = 2.74e-5`, 260/661 cells above print resolution. No `--ref` overlay — the 5 s run shares no frame time with the 4000 s release reference.

```bash
.venv/bin/python util/animate_two_temp.py outputs/model_column_2t/smoke5s.txt.twotemp --fps 4
.venv/bin/python util/animate_two_temp.py outputs/model_column_2t/smoke5s.txt.twotemp --snapshot
```
