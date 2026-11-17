# Visualization catalog — presentation assets

Posters, crops and copies made for a talk. Each deck under `docs/presentations/` keeps its own figures so it builds from a fresh clone; these are the commands that put them there. See [README.md](README.md) for conventions.

---

## 2026-07-11 group meeting: movie posters (linked, not embedded)

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
`animate_loop_flare.py` honours `ANIM_CRF` / `ANIM_PIX_FMT` / `ANIM_PRESET` (encode, see `util/lib/_anim_parallel.py`),
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
.venv/bin/python util/run/run_ensemble.py --manifest outputs/events/event_te_2x_jobs.json --workers 2
# (1) render near-lossless, 0–100 s (t_cap), 60 fps, uniform-time, post-beam 2x (env inline — zsh won't word-split a var)
ENV="FLARE_DELTA=4 FLARE_E_CUT=20 FLARE_T_ON=2 FLARE_DUR=10 ANIM_UNIFORM_TIME=1 ANIM_SPEED_PRE=1 ANIM_SPEED_POST=2 ANIM_FPS=60 ANIM_MAX_FRAMES=1000 ANIM_CRF=12 ANIM_PIX_FMT=yuv444p ANIM_PRESET=slow"
env ${=ENV} .venv/bin/python util/legacy/animate_loop_flare.py scenarios/data/event_ensemble/c05_101Mm.dat  outputs/events/event_te/c05_101Mm_te.txt  visualization/events/event_20240801_c05_101Mm_Te.mp4 100
env ${=ENV} .venv/bin/python util/legacy/animate_loop_flare.py scenarios/data/event_ensemble/o00_1040Mm.dat outputs/events/event_te/o00_1040Mm_te.txt visualization/events/event_20240801_o00_1040Mm_Te.mp4 100
# (2) extract the first-frame poster shown on each slide
ffmpeg -y -i visualization/events/event_20240801_c05_101Mm_Te.mp4  -frames:v 1 -q:v 2 visualization/events/poster_c05_closed.jpg
ffmpeg -y -i visualization/events/event_20240801_o00_1040Mm_Te.mp4 -frames:v 1 -q:v 2 visualization/events/poster_o00_open.jpg
```

---

## 2026-08-17 group meeting: figure crops

The deck at `docs/presentations/2026-08-17-group-meeting/` keeps its own figures so it
builds from a fresh clone — `visualization/` is gitignored, so a deck reaching into it
would not build for anyone else. Two provenances:

**Copied unchanged** from `visualization/reference/`:

```bash
cp visualization/reference/c7_ionization_profile_compare.png \
   docs/presentations/2026-08-17-group-meeting/figures/
```

Note that `util/release/plot_c7_ionization_profile.py --mode compare` writes to
`visualization/eos_gamma/`, and the `visualization/reference/` copy was placed by hand;
regenerating it also needs `outputs/eos_gamma/gamma_hydrogen_excitation.dat`, which is
not in the tree (`bash util/eos/build_and_run.sh excitation`).

**Cropped** from the six-panel comparison of the original storage-precision section (now in [_historical.md](_historical.md)) — one panel per slide, because the
whole figure is unreadable at slide size:

> **WITHDRAWN** — contaminated hybrid float32/float64 build (`CHROMO_STATE_FLOAT64` promoted `Vec` only); its numbers are wrong. See `docs/studies/numerics/state_precision_release_cutover.md`; the source figure `precision_control_N500_4000s.png` is from that run; use the cutover comparison in [model_column.md](model_column.md) instead.

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
aside slide — and the other three are backup-slide assets.

**Cropped** from the wall-temperature sweep figures in [model_column.md](model_column.md), same reasoning:

```bash
.venv/bin/python - <<'PY'
from PIL import Image
out = "docs/presentations/2026-08-17-group-meeting/figures/"

# 2x3 summary. Two hand-tuned offsets: TRIM drops the sliver of the next panel's
# y-label that matplotlib leaves just inside the nominal W/3 boundary, and SPLIT
# is the row cut, placed just above the second row's panel titles so neither row
# clips or inherits the other's title.
TRIM, SPLIT = 52, 16
src = Image.open("visualization/model_column/twall_sweep_N500_4000s_summary.png")
W, H = src.size
top = int(0.055 * H); mid = top + (H - top) // 2
rows = [(top, mid + SPLIT), (mid + SPLIT, H)]
cols = [(0, W // 3 - TRIM), (W // 3 - TRIM, 2 * W // 3 - TRIM), (2 * W // 3 - TRIM, W)]
for name, (r, c) in {"twall_response":        (0, 0),   # V vs T_wall
                     "twall_local_slope":     (0, 2),   # log-log + local exponent
                     "twall_enthalpy_budget": (1, 2)}.items():
    src.crop((cols[c][0], rows[r][0], cols[c][1], rows[r][1])).save(out + name + ".png")

# 1x3 support: left panel only (the height-independent mass flux)
sup = Image.open("visualization/model_column/twall_sweep_N500_4000s_support.png")
W2, H2 = sup.size
sup.crop((0, int(0.055 * H2), W2 // 3 - 46, H2)).save(out + "twall_massflux_uniform.png")
PY
```

`twall_response.png` and `twall_enthalpy_budget.png` carry the two main upper-boundary
slides; `twall_massflux_uniform.png` is on the caveats backup slide (sized by `height`,
not `width` — the crop is near-square); `twall_local_slope.png` is an unused backup asset.
See `docs/presentations/2026-08-17-group-meeting/figures/README.md`.

### `presentation.pdf`

```bash
cd docs/presentations/2026-08-17-group-meeting
latexmk -pdf -output-directory=build presentation.tex
```

---

---

## 2026-09-29 tp-sweep deck

### Presentation crops (`docs/presentations/2026-09-29-tp-sweep/figures/`)

The deck keeps its own figures so it builds from a fresh clone. `analytic_J_integrand.png` is copied unchanged; the rest are single-panel crops (panel titles removed, captions are on the slides) of the section 20–21 figures:

```bash
cp visualization/model_column/analytic_J_integrand.png docs/presentations/2026-09-29-tp-sweep/figures/
.venv/bin/python - <<'PY'
from PIL import Image
V = 'visualization/model_column/'; D = 'docs/presentations/2026-09-29-tp-sweep/figures/'
crops = {
    'tp_sweep_sim_vs_analytic.png': {
        'sweep_J_vs_p.png': (0, 95, 745, 672), 'sweep_parity.png': (745, 95, 1495, 672),
        'sweep_hc_ratio.png': (1495, 95, 2240, 672), 'sweep_ER.png': (0, 745, 745, 1330),
        'sweep_FE_nom.png': (745, 745, 1495, 1330), 'sweep_FLJ_nom.png': (1495, 745, 2240, 1330)},
    'tp_sweep_time_series.png': {
        'series_ER.png': (0, 50, 2600, 496), 'series_unsteady.png': (0, 925, 2600, 1365)},
    'tp_corner_test.png': {
        'corner_std_profile.png': (180, 140, 1095, 610), 'corner_top_flux.png': (1095, 140, 2000, 610),
        'corner_spacetime.png': (180, 600, 2060, 1150)},
    'analytic_J_table.png': {
        'Jtable_T.png': (0, 110, 805, 660), 'Jtable_p.png': (805, 110, 1650, 660),
        'Jtable_Tc.png': (0, 670, 805, 1275), 'Jtable_exponent.png': (805, 670, 1650, 1275)},
}
for src, outs in crops.items():
    im = Image.open(V + src)
    for name, box in outs.items():
        im.crop(box).save(D + name)
PY
```

The pixel boxes assume the current figure sizes (2240×1330, 2600×1365, 2210×1235, 1650×1275); re-check them if a parent figure is redrawn at another size.

---

## 2026-10-11 hydrostatic-limiter deck (advisor follow-up)

### Figures (`docs/presentations/2026-10-11-hydrostatic-limiter/figures/`)

Unchanged copies of the three study figures catalogued in [model_column.md](model_column.md). `hse_recon_evaporation_bars.pdf` is Result 2; `hse_recon_evaporation.pdf` is its backup slide.

```bash
.venv/bin/python util/release/plot_hse_recon_comparison.py
cp visualization/model_column/hse_recon_static_residual.pdf \
   visualization/model_column/hse_recon_evaporation.pdf \
   visualization/model_column/hse_recon_evaporation_bars.pdf \
   docs/presentations/2026-10-11-hydrostatic-limiter/figures/
```

The two schematics (what the limiter compares; the hydrostatic version) are TikZ inside `presentation.tex`.

### `presentation.pdf`

```bash
cd docs/presentations/2026-10-11-hydrostatic-limiter
latexmk -pdf -output-directory=build presentation.tex && cp build/presentation.pdf .
```

---

## 2026-10-27 outer-BC and TR-pulse deck (advisor follow-up)

### Figures (`docs/presentations/2026-10-27-outer-bc-and-tr-pulses/figures/`)

Unchanged copies of two figures catalogued in [model_column.md](model_column.md) (outer pressure closure rerun and hot-wall TR pulses).

```bash
.venv/bin/python util/release/plot_outer_bc_and_tr_pulses.py      # on the mac-mini, where the runs live
cp visualization/model_column/outer_bc_rerun_timeseries.pdf \
   visualization/model_column/tr_pulses_resolution.pdf \
   docs/presentations/2026-10-27-outer-bc-and-tr-pulses/figures/
```

### Movies (`docs/presentations/2026-10-27-outer-bc-and-tr-pulses/movies/`)

Unchanged copies of four movies catalogued in [model_column.md](model_column.md) (same section).

```bash
cp visualization/model_column/outer_bc_rerun_22kK.mp4 \
   visualization/model_column/outer_bc_rerun_22kK_pressure.mp4 \
   visualization/model_column/outer_bc_rerun_100kK.mp4 \
   visualization/model_column/tr_pulses_N500_vs_N1000_100kK.mp4 \
   docs/presentations/2026-10-27-outer-bc-and-tr-pulses/movies/
```

### `presentation.pdf`

```bash
cd docs/presentations/2026-10-27-outer-bc-and-tr-pulses && latexmk presentation.tex
```
