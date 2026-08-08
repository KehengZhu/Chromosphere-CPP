# Pressure-based primitive reconstruction for the Gamma/Saha `model_column` path

> **Superseded by the model_column release-numerics promotion.** Roe (`roe-local`) and
> `(ln rho, V, ln p)` are now the `model_column` release defaults and need no environment
> override; `ISO_RIEMANN=rusanov` and `ISO_RECONSTRUCTION=lnrho-v-lnt` are reference-only.
> Statements below describe the configuration in force **when this study was run** and are
> kept for chronology. Current release description: `docs/model_column_release_numerics_recap.md`.

**Status: GO** as a diagnostic option, with a recommendation to promote. Release defaults unchanged.

## Question

The Roe-local corrector (`ISO_RIEMANN=roe-local`, commit `5271ca2`) removed the short-time upper-TR velocity ripple, but an N500 run to 4000 s developed a large, stable, grid-scale velocity/`rho V` oscillation localized at the equilibrium ionization transition near 2141 km. The conservative effective mass flux stayed smooth throughout, and N1000 was much cleaner, so the mode was numerical.

The hypothesis under test: the active reconstruction limits `(ln rho, V, ln T)` component-wise and only afterwards obtains the face pressure through the nonlinear Saha map `p(rho,T)`. Two independently limited variables therefore set one mechanical quantity, and a mechanically smooth state acquires an artificial reconstructed pressure jump wherever `x_eq` turns over. Rusanov's acoustic-scale dissipation hides it; a low-dissipation Roe solver does not.

## Implementation

`Grid::pressure_reconstruct` (default false) selects the thermal MUSCL variable of the Gamma/Saha path:

| flag | reconstructed | face closure |
|---|---|---|
| false (release) | `(ln rho, V, ln T)` | `equilibrium_mixture_face_state_from_logs` |
| true | `(ln rho, V, ln p)` | `equilibrium_mixture_face_state_from_log_pressure` |

`p` is the authoritative total `p_i + p_n = (1+x) n_H k_B T` — the same total the momentum flux and the Roe pressure use, with the electron pressure already inside `p_i`.

The face temperature is recovered by `equilibrium_temperature_from_density_pressure(rho, p)`. The closure supplies its own exact bracket: with `T1 = p m_H/(rho k_B)`, `x in (0,1)` puts `T` in `(T1/2, T1]`, and `dp/dT|_rho = n_H k_B (1 + x + T dx/dT) > 0` makes `p(T)` strictly monotone there. The solver is a safeguarded Newton iteration that bisects on any step leaving that factor-of-two bracket, and throws rather than returning invalid thermodynamics. No table bracket, no floors, no clamps.

Both decode sites (`decode_mixture_field_into` and the predictor `decode_predicted_caloric_primitives_into`) write the same slot-2 quantity, so predictor and corrector are consistent; ghosts go through the same `put()`. Everything downstream — conservative update, equilibrium projection, sources, conduction, boundaries, `eq_wb` residual, Roe solver — is untouched.

Selector: `ISO_RECONSTRUCTION=lnrho-v-lnt` (default) | `lnrho-v-lnp`.

## Evidence

All runs: `model_column`, physical conduction only, N500 = 661 cells / N1000 = 1320 cells, R4, CFL 0.50, 22,000 K physical conductive face, `ISO_RIEMANN=roe-local`. Analysis: `util/pressure_reconstruction_diag.py`. `eps_p` is the reconstructed corrector face-pressure mismatch `|p_R-p_L| / mean` over 2139.5–2143.5 km; `Q_V`/`Q_M` are the established roughness metrics over 2130–2150 km.

### t = 1000 s

| | lnT N500 | lnP N500 | lnT N1000 | lnP N1000 |
|---|---|---|---|---|
| `eps_p` max | 0.6217 % | **0.0020 %** | 0.1085 % | 0.0004 % |
| `Q_V` | 0.1685 | **0.00285** | 0.0128 | 0.00080 |
| `Q_M` | 0.2262 | **0.00152** | 0.0161 | 0.00043 |
| min V [m/s] | −21.38 | **+3.68** | +0.40 | +3.86 |
| mean V [m/s] | 6.038 | 6.639 | 6.755 | 6.827 |
| `F_eff` [kg m^-2 s^-1] | 4.286e-10 | 4.276e-10 | 4.476e-10 | 4.464e-10 |
| `T_top` [K] | 21904.7 | 21904.7 | 21950.9 | 21950.9 |
| BL thickness [km] | 11.87 | 11.87 | 11.73 | 11.73 |
| `q_phys` [W m^-2] | 0.4562 | 0.4563 | 0.4701 | 0.4697 |
| cumulative [J m^-2] | 569.5 | 569.2 | 587.4 | 587.1 |

### N500 to t = 4000 s

| | lnT | lnP |
|---|---|---|
| `eps_p` max | 0.8017 % | **0.0023 %** |
| `Q_V` | 0.2330 | **0.00311** |
| `Q_M` | 0.3212 | **0.00185** |
| min / max V [m/s] | −30.15 / 19.23 | **+2.97 / 11.60** |
| mean V [m/s] | 4.894 | 5.693 |
| `F_eff` | 3.546e-10 | 3.571e-10 |
| `T_top` [K] | 21913.3 | 21912.8 |
| BL thickness [km] | 12.14 | 12.14 |
| cumulative conduction [J m^-2] | 1849.7 | 1853.8 |
| mass drift | −1.44e-4 | −1.63e-4 |

### Reading

1. The mismatch drops by a factor of ~350 and no longer grows with time (0.0014 % at 500 s to 0.0023 % at 4000 s, versus 0.27 % to 0.80 %). It is now at the limiter's own truncation level.
2. The late grid mode is removed, not damped: `Q_V` falls 75x, `Q_M` 174x, and there is **no negative-velocity cell at any captured time** where the old reconstruction had 3–4.
3. `Q_eff` (the conservative effective mass flux) is ~6e-4 in every run — it was already smooth and stays smooth. `F_eff` moves by 0.7 %.
4. Thermal and conductive observables are unchanged to 4 significant figures at fixed resolution. The change is confined to hydro reconstruction, which is the signature of removing a discretization defect rather than altering physics.
5. Resolution agreement improves sharply. N500 to N1000 in mean V: 6.038 to 6.755 (+11.9 %) for lnT, 6.639 to 6.827 (+2.8 %) for lnP. In min V: 21.8 m/s versus 0.19 m/s. **lnP at N500 is 4.5x smoother than lnT at N1000.**
6. Short-time behaviour is preserved: at 100 s lnP-Roe keeps `Q_V` 0.00133 / `Q_M` 0.00070 against Rusanov's 0.0138 / 0.0078, so the Roe advantage is retained and no Rusanov-scale diffusion was reintroduced.

### Equilibrium control (conduction off, `ISO_RELAX_TIME=1e9`, 1000 s)

Residual velocity saturates in both: lnT `max|V|` 0.024 m/s (rms 6.4e-3), lnP 0.097 m/s (rms 2.5e-2). Both are bounded with no growth after ~300 s, and both are two orders below the physical evaporation signal. lnP holds `T_top` better (−0.4 K versus −2.5 K over 1000 s). The larger lnP velocity residual is a truncation-amplitude difference in an artificial control configuration, not an instability, and it does not appear in the production configuration, where lnP is the smoother solution by two orders of magnitude.

## Root cause

Confirmed. Reconstructing `(ln rho, V, ln T)` through the Saha EOS was the dominant source of the late coarse-grid Roe ripple. The mismatch magnitude, its location (exactly where `x_eq` turns over), its growth history, and its resolution scaling all track the ripple; removing it removes the ripple while leaving every conservative and thermal observable in place. The earlier Property-U diagnostic (normalized residual ~5e-4 to 1e-3, not growing with the ripple) is consistent with this: the Roe average was never the problem.

## Reproduction

```bash
ISO_NS=500 ISO_RIEMANN=roe-local ISO_RECONSTRUCTION=lnrho-v-lnp \
CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FACE_FLUX_DIAG=1 \
CHROMO_OUTER_COND_DIAG=1 CHROMO_FRAME_DT=20 CHROMO_FACE_FLUX_STRIDE=2000 \
CHROMO_OUTER_COND_STRIDE=400 \
scripts/run_chromo_realtime.sh outputs/model_column/lnp_roe_N500_4000s.txt \
    full no-ionization model_column - 20.0 no-cooling

python3 util/pressure_reconstruction_diag.py --targets 500,1000,4000 \
  --run "lnT-N500=outputs/model_column/lnt_roe_N500_4000s.txt" \
  --run "lnP-N500=outputs/model_column/lnp_roe_N500_4000s.txt"
```

## Recommendation

Promote `pressure_reconstruct` to the `model_column` default in a separate, explicitly scoped release-policy change, after a Rusanov N500/N1000 acceptance pass (the release solver is still Rusanov, and the 100 s Rusanov comparison already shows lnP is neutral-to-better there). This task does not change any production default.
