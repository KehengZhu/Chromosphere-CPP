# Figures for the 2026-09-29 tp-sweep talk

Kept inside the talk's own folder so the deck stays buildable from a fresh clone.
`visualization/` is gitignored, so a deck that reached across into it would not build for anyone else.
The crop commands are recorded in `docs/viz_catalog/presentations.md` (2026-09-29 tp-sweep deck).

| File | Provenance |
| --- | --- |
| `analytic_J_integrand.png` | Copy of `visualization/model_column/analytic_J_integrand.png` (`util/release/plot_J_integrand.py`). |
| `Jtable_T.png`, `Jtable_p.png`, `Jtable_Tc.png`, `Jtable_exponent.png` | Panels (a), (b), (c), (d) of `visualization/model_column/analytic_J_table.png` (`util/release/analytic_evaporation_table.py`). |
| `sweep_J_vs_p.png`, `sweep_parity.png`, `sweep_hc_ratio.png` | Panels (a), (b), (c) of `visualization/model_column/tp_sweep_sim_vs_analytic.png` (`util/release/plot_tp_sweep_vs_analytic.py`). |
| `sweep_ER.png`, `sweep_FE_nom.png`, `sweep_FLJ_nom.png` | Panels (d), (e), (f) of the same figure. `sweep_FE_nom.png` is **not used**; kept as a backup asset. |
| `series_ER.png`, `series_unsteady.png` | Top and bottom rows of `visualization/model_column/tp_sweep_time_series.png`. |
| `corner_std_profile.png`, `corner_top_flux.png`, `corner_spacetime.png` | Panels (a), (b) and (c)–(f) of `visualization/model_column/tp_corner_test.png` (`util/release/tp_corner_test_compare.py`). |

Studies behind them: `docs/studies/boundaries/twall_mass_flux_relaxation_timescales.md` (derivation of J), `docs/studies/boundaries/analytic_evaporation_J_table.md` (the table), `docs/studies/boundaries/tp_sweep_evaporation_closure_validation.md` (the sweep).
