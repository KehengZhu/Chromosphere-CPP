# Visualization catalog

The exact command that regenerates every file under `visualization/`. It was a single root file, `visualize_commands.md`, until the 2026-10 reorganization split it by scenario area.

**Rule: every file in `visualization/` has an entry here.** Whenever you generate a figure or movie, add or update its entry in the matching file below in the same change: the output filename, the exact command (run from the project root), and a one-line description. If you change a script's invocation, update its entries too.

## Files

| File | What it covers |
| --- | --- |
| [model_column.md](model_column.md) | The release `model_column` solver: release mesh and face-flux tools, the Godunov cutover, the storage-precision cutover, the conduction-wall (`T_wall`) sweep, the outer-face closure swap, the analytic `J` table, the `(T_wall, p_top)` sweep and the electron–ion equilibration time |
| [eos.md](eos.md) | The pure-hydrogen Saha closure and the CRASH `Gamma1` table |
| [reference.md](reference.md) | Physics reference plots in `visualization/reference/` |
| [tutorials.md](tutorials.md) | Educational Godunov and Roe face-flux animations |
| [presentations.md](presentations.md) | Posters, crops and copies made for the dated decks under `docs/presentations/` |
| [_historical.md](_historical.md) | Retired scenarios and configurations: two-fluid C7, open-BC / `rho_tot`, flare and loop, PFSS and event campaigns, synthetic Doppler, gentle evaporation, `model_isentropic`, the pre-release `model_column` development, the **withdrawn** hybrid float32/float64 comparison, and the two-temperature column `model_column_2t` |

Put a new entry in the file for its scenario area. When a configuration is retired, move its entries to `_historical.md` in the same change that marks its study `SUPERSEDED` or `HISTORICAL`.

## Conventions

- Run every command **from the project root**.
- The interpreter is the project venv: **`.venv/bin/python`** (see the root `README.md`).
- Plotting and animation scripts live in `util/`, grouped as `util/release/` (release `model_column` diagnostics), `util/legacy/` (retired scenarios, including `model_column_2t`), `util/run/` (sweep and ensemble drivers) and `util/lib/` (shared modules). Two physics-reference scripts live in `visualization/reference/` and write next to themselves. Outputs land in `visualization/` (gitignored).
- Commands recorded before the reorganization named the flat `util/<script>.py` path; they have been rewritten to the script's current location.
- **Both `outputs/` and `visualization/` are organized into per-scenario subdirectories** mirroring `scenarios/`: `model_c7/`, `model_column/` (the unified column — absorbs the old `iso_*`/`gentle_*` runs), `model_flare/`, `analytic_canopy/`, `loops/`, `pfss/`, `events/` (real-event campaigns: 2024-08-01 AR13768, SOL2014, ensembles, `event_te/`, `sol2014/`, …), and `two_fluid/`. `visualization/` also has `reference/` (physics-reference plots). Dumps and figures produced by now-deleted diagnostic code (`dynamic_hse_*`, `cond_hydro_*`) and retired analytic-isentrope stages live under `_archive/`.
- **Paths in this catalog are the paths a command *writes* to, which for an archived figure is no longer where the file sits.** Figures move to `visualization/<scenario>/_archive/<bucket>/` once the configuration they measured stops being current; the commands are left pointing at the top-level path so they still regenerate the figure in place. To locate an existing file, consult the `_archive/INDEX.md` of its scenario — `visualization/model_column/_archive/INDEX.md` lists which figures stayed release-facing and which bucket each of the rest went to (most recently on 2026-08-18, when the pre-Godunov Roe figures, the pre-cutover float32 figures, and the `q_ramp` two-fluid driver were archived).
- **Code-generation folders.** When the release numerics change, everything at the top level of `visualization/model_column/` that was made from the earlier code moves to a folder named `<date>_<last commit of that code>/`, with an `INDEX.md` saying what the folder contains. The top level then holds only outputs of the current release. The first such folder is `2026-10-12_7a0e1af/`, which holds everything made before the hydrostatic reconstruction became the release. Commands keep their top-level output paths, as with `_archive/`.
- Input simulation dumps live under `outputs/<scenario>/`; loop/field-line geometry `.dat` files live under `scenarios/data/`; magnetograms under `util/data/`. Some `outputs/` subdirectories are gitignored — adjust the input paths to your actual run dumps where noted.
- **`lnp_godunov_N500_4000s_f64` is WITHDRAWN.** It was built with the retired `CHROMO_STATE_FLOAT64` option, which promoted only `Vec` and left the rest of the update chain in `float` — an incoherent hybrid. Do not quote it or overlay it; see `docs/studies/numerics/state_precision_release_cutover.md`.
