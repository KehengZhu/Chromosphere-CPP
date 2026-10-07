# Visualization catalog — physics reference plots

Reference plots whose scripts live next to their output in `visualization/reference/`. See [README.md](README.md) for conventions.

---

## Physics reference plots (scripts live in `visualization/reference/`)

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
