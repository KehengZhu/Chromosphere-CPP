# Analytic steady-evaporation lookup table $J(T_\mathrm{top}, p_\mathrm{top}, T_c)$

**Status: CURRENT**

**Question.** Can the steady evaporation response of the release column to a boundary state $(T, p)$ be tabulated without a simulation per table point? This is the first step toward a steady lookup table for an AWSoM-R lower boundary. Physical placement of the interface is set aside here.

**Answer.**

1. **Yes, for $J = F L_\mathrm{TR}$, not for $F$.** The release TR integral (1) of `twall_mass_flux_relaxation_timescales.md` needs one closure for the chromospheric enthalpy, $h_c = h_\mathrm{Saha}(T_c, p_\mathrm{top})$. With it, $J$ follows from $(T_\mathrm{top}, p_\mathrm{top}, T_c)$ alone. $F$ still needs a length $L_\mathrm{TR}$, which TR physics does not set in the release.
2. **The closure is accurate; the table's error is the simulation's own $F L/J$ offset.** At the validated line ($p_\mathrm{top} = 0.0103$ Pa, 15–100 kK, 4000 s):
   - With the run's actual $T_c$, the table reproduces the $J$ formed from the measured $h_c$ to 0.1–3.3 %.
   - With the nominal $T_c = 6.7$ kK, the measured $F L_\mathrm{TR}$ is **0.87–0.94** of the table.
3. **$T_c$ is a hidden variable, and it matters at low pressure.** At 0.0103 Pa, $T_c = 6.7$ kK sits in the middle of Saha hydrogen ionization ($x = 0.64$), so $h_c$ is exponentially sensitive to it. At 22 kK, $T_c = 6.35 / 7.0$ kK changes $J$ by $\times 0.70 / 1.47$. At 1 Pa the same change is only $\times 0.94 / 1.06$. AWSoM does not supply $T_c$, so it must either be a table dimension or be declared a fixed model choice.
4. **$p_\mathrm{top}$ is a weak parameter where the table is well defined.** $\partial\ln J/\partial\ln p$ lies between −0.39 and +0.11 everywhere outside the marker-dependent region (point 5). Over four decades, 1e-3 to 10 Pa, $J$ changes by 0.78–1.19 at 95 kK and 0.70–1.38 at 50 kK.
5. **At high $p$ and low $T_\mathrm{top}$, $J$ depends on how the TR base is defined.** There the Saha ionization zone moves into the TR: $x(7.5\ \mathrm{kK}) = 0.32$ at 1 Pa and 0.11 at 10 Pa. Raising the marker $T_1$ from 7.5 to 10 kK then removes 42–92 % of $J$ at 15–22 kK and 1–10 Pa. The steep rise of $J$ with $p$ at 15 kK in panel (b) is this effect. In that region, $J$ and $L_\mathrm{TR}$ are meaningful only as a matched pair defined at the same isotherm.

Evidence:
- script: `util/release/analytic_evaporation_table.py`
- table: `outputs/model_column/analytic_J_table/J_table.{npz,csv}`, `J_table_summary.json`
- figure: `visualization/model_column/analytic_J_table.png` (`docs/viz_catalog/model_column.md`)

Post-processing only; no new run.

## Method

$$J(T_\mathrm{top}, p_\mathrm{top}, T_c)=\int_{T_1}^{T_\mathrm{top}}\frac{\kappa(T,p_\mathrm{top})}{h(T,p_\mathrm{top})-h_c}\,dT,\qquad h_c=h_\mathrm{Saha}(T_c,p_\mathrm{top}),\quad T_1=7.5\ \mathrm{kK}.$$

The ingredients are exactly those of the relaxation study:
- the pure-hydrogen Saha EOS and `j_theory` of `util/release/twall_relaxation_diag.py`;
- the Spitzer-plus-neutral $\kappa$ of `util/lib/conductive_flux.py`.

The integral assumes the same things as there: steady, subsonic, isobaric, thin, no radiation, no heating, $q = 0$ at the chromosphere top. Evaluating $h_c$ at $p_\mathrm{top}$ rather than at the chromosphere-top pressure is part of the isobaric-TR assumption.

The grid has 29 values of $T_\mathrm{top}$ from 15 to 200 kK (log-spaced), 33 values of $p_\mathrm{top}$ from $10^{-3}$ to 10 Pa (log-spaced), and $T_c \in \{6.0, 6.35, 6.7, 7.0\}$ kK. The table stores $J$, $h_c$, $x_c$ and the local exponents with respect to $\ln T$ and $\ln p$, plus $J$ with $T_1 = 10$ kK for the marker check.

## End-to-end check against the sweep (4000 s, $p_\mathrm{top} = 0.0103$ Pa)

| $T_\mathrm{wall}$ [kK] | run $T_c$ [K] | $h_{c,\mathrm{run}}/h_c(6.7\,\mathrm{kK})$ | $J_\mathrm{tab}(T_{c,\mathrm{run}})/J(h_{c,\mathrm{run}})$ | $F L_\mathrm{meas}/J_\mathrm{tab}(6.7\,\mathrm{kK})$ |
|---:|---:|---:|---:|---:|
| 15 | 6708 | 0.986 | 1.033 | 0.869 |
| 18 | 6713 | 0.998 | 1.021 | 0.916 |
| 22 | 6715 | 1.006 | 1.011 | 0.923 |
| 30 | 6690 | 0.983 | 1.006 | 0.920 |
| 45 | 6645 | 0.938 | 1.003 | 0.921 |
| 70 | 6579 | 0.868 | 1.002 | 0.929 |
| 100 | 6504 | 0.787 | 1.001 | 0.936 |

"Run $T_c$" is the temperature at the cell where the diagnostic reads $h_c$ ($|q| < 1\,\%\ q_\mathrm{wall}$). At the hot walls it falls below 6.7 kK by 4000 s through the adiabatic contraction. That fall is part of why the fixed-$T_c$ table drifts from the run, and in a coupled use $T_c$ is itself time dependent in the release.

The remaining 6–13 % comes from the measured $F L/J$ of the relaxation study: 0.89–1.01 at 4000 s and 0.95–0.98 at $t_\mathrm{ready}$. It combines two things: the energy cost of ionizing the chromospheric gas the lengthening TR absorbs, and, at 15–18 kK, an incompletely relaxed TR.

## Sensitivities (from `J_table_summary.json`)

**$T_c$**, as $J(T_c)/J(6.7\,\mathrm{kK})$ for $T_c = 6.0 / 6.35 / 7.0$ kK:

| $T_\mathrm{top}$ | 0.001 Pa | 0.0103 Pa | 0.1 Pa | 1 Pa |
|---|---|---|---|---|
| 22 kK | 0.41 / 0.66 / 1.26 | 0.57 / 0.70 / 1.47 | 0.82 / 0.88 / 1.22 | 0.90 / 0.94 / 1.12 |
| 50 kK | 0.64 / 0.83 / 1.09 | 0.71 / 0.81 / 1.21 | 0.88 / 0.92 / 1.13 | 0.95 / 0.97 / 1.05 |

This is not an artefact of $T_c$ approaching the 7.5 kK marker. With $T_1 = 10$ kK the 22 kK, 0.0103 Pa row is 0.58 / 0.71 / 1.43.

**Marker**, as $J(T_1 = 10\,\mathrm{kK})/J(7.5\,\mathrm{kK})$ at $T_c = 6.7$ kK. It is $\geq 0.86$ for $T_\mathrm{top} \geq 22$ kK and $p \leq 0.1$ Pa, and $\geq 0.94$ for $T_\mathrm{top} \geq 50$ kK and $p \leq 1$ Pa. It falls to 0.58 / 0.25 at 22 kK and 1 / 10 Pa, and to 0.26 / 0.08 at 15 kK. The dashed contour in panel (d) bounds this marker-dependent region.

**Exponents.** $\partial\ln J/\partial\ln T_\mathrm{top} = 2.3$–3.1 wherever the marker does not matter, tending to the Spitzer 5/2. Outside the marker-dependent region $\partial\ln J/\partial\ln p$ lies between −0.39 and +0.11, and it is negative at $p \leq 0.1$ Pa. Two effects compete: a higher $p$ lowers $x_c$, hence $h_c$, which widens $h - h_c$ and lowers $J$; it also lowers $x$ inside the TR, which narrows $h - h_c$ and raises $J$. The first wins at low $p$. Inside the marker-dependent region the exponent reaches $+0.5$.

## Limits

- **Validated on the sweep grid only.** The original simulation checks were at $p_\mathrm{top} = 0.0103$ Pa, $T_c \approx 6.7$ kK and 15–100 kK. `tp_sweep_evaporation_closure_validation.md` extends them to 18–70 kK × 0.003–0.08 Pa with consistently rebuilt ICs (`ISO_P_TOP`). There, where the column is quiet and relaxed, $F L/J$ with the nominal $T_c$ is 0.93–1.03 at $t_\mathrm{ready}$. By 5000 s it drifts to 0.84–0.92 at $p_\mathrm{top} \leq 0.0103$ Pa, because the run's $T_c$ falls (point 3 above, now measured). Everything outside that grid is still the integral alone: 100–200 kK, pressures below 0.003 or above 0.08 Pa, and any other $T_c$. In the corner $T_\mathrm{top} \geq 50$ kK with $p_\mathrm{top} \geq 0.02$–0.04 Pa, the N=500 flux carries grid-scale noise from the unresolved TR foot. It is valid there only as a ≳ 500 s mean; at N=1000 the ratio is 0.975. The A/B construction cross-check has not been done.
- **Release physics only.** No radiation and no heating. Saha LTE sets $h_c$, so the $T_c$ sensitivity is a property of pure-hydrogen Saha at the chromosphere top. A real, non-LTE chromosphere would respond differently.
- **Steady only.** Use it only when the boundary changes more slowly than $t_\mathrm{ready}$: 520–2300 s on the C7 line, and up to 4500 s at 0.04–0.08 Pa and 18–35 kK, where the TR is still relaxing at 5000 s.
- **Not a flux.** Converting $J$ to $F$ needs an $L_\mathrm{TR}$ defined at the same isotherm $T_1$. In the release that length is geometric.
- **Below 15 kK** (near the flow reversal) the table is not defined; the sweep reverses flow at ≈ 10.4 kK from IC drainage, which the steady integral does not describe.

## Reproduce

```bash
.venv/bin/python util/release/analytic_evaporation_table.py
```

This needs the sweep runs and `outputs/model_column/twall_sweep/twall_relaxation_metrics.json` from `util/release/twall_relaxation_diag.py`. Runtime is about 6 s.
