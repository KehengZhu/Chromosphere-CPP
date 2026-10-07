# Simulated against analytic evaporation flux over the $(T_\mathrm{wall}, p_\mathrm{top})$ sweep

**Status: CURRENT**

**Question.** The steady-TR closures of `twall_mass_flux_relaxation_timescales.md` and `analytic_evaporation_J_table.md` had been checked against the simulation on one pressure line only, the C7 column at $p_\mathrm{top} = 0.0103$ Pa. Do they still hold when the top pressure changes, over the 0.003–0.08 Pa that an AWSoM-R transition-region top can take?

The closures are:
- the energy closure $F_E = q_\mathrm{wall}/(h_\mathrm{top} - h_c)$, which borrows $q_\mathrm{wall}$ from the run;
- the TR closure $F_J = J(T_\mathrm{top}, p_\mathrm{top}, h_c)/L_\mathrm{TR}$, which borrows $L_\mathrm{TR}$ from the run.

Each is evaluated with the $h_c$ measured in the run and with the nominal $h_c = h_\mathrm{Saha}(6.7\ \mathrm{kK}, p_\mathrm{top})$ that a lookup table would use. Notation is that of `twall_mass_flux_relaxation_timescales.md`. In particular, $\mathrm{ER} = F(h_\mathrm{top} - h_c)/q_\mathrm{wall}$ is $F_\mathrm{sim}/F_E$ with the measured $h_c$.

**Answer.**

1. **The sweep is numerically clean.** All 31 runs reached `termination=end_time` at 5000 s with zero Godunov fallbacks.
2. **The closures hold off the C7 line wherever the column is quiet and its TR has relaxed.** That is 21 of the 30 N=500 points (class "relaxed" below). In those:
   - $\mathrm{ER}$ is **0.93–0.99 at $t_\mathrm{ready}$** and 0.90–0.98 at 5000 s;
   - $F L_\mathrm{TR}/J$ with the nominal $T_c$ is 0.93–1.03 at $t_\mathrm{ready}$.

   The measured-$h_c$ energy closure is never exceeded: the simulation carries 1–10 % less enthalpy flux than the conducted heat. The 0.93–0.98 found on the C7 line therefore extends unchanged across a factor of 27 in pressure.
3. **A fixed-$T_c$ table drifts low at low pressure as the run ages. The run's $T_c$ is what moves.** By 5000 s:
   - $F L_\mathrm{TR}/J$ with the nominal $T_c$ has fallen to 0.84–0.92 at $p_\mathrm{top} \leq 0.0103$ Pa, and it stays at 0.94–1.01 above that.
   - The measured $h_c$ has fallen below the nominal one, to $0.80\times$ at 70 kK and 0.003 Pa: the lifted chromosphere cools adiabatically, and the release has no chromospheric energy source to stop it.
   - The closure with the *measured* $h_c$ is unaffected: $\mathrm{ER} = 0.98$ at the same point.

   This measures in the simulation the hidden-variable sensitivity that `analytic_evaporation_J_table.md` point 3 predicted from the integral: at low pressure, $h_c$ is exponentially sensitive to $T_c$.
4. **At high pressure and low wall temperature, the TR relaxes more slowly than the run lasts.**
   - The initial C7-shaped TR carries excess energy ($\mathrm{ER} \approx 1.5$–2 early on).
   - It drains on an e-folding time that grows with pressure: at 18 kK, from about 370 s at 0.003 Pa to about 3500 s at 0.08 Pa.
   - Four runs are therefore still relaxing at 5000 s: 18 and 25 kK at 0.04 and 0.08 Pa. $\mathrm{ER}$ is still falling monotonically in all four.
   - Their ratios above 1 (up to $\mathrm{ER} = 1.17$ and $F L/J = 1.26$) are upper bounds on the steady value. They are not a closure failure.

   These runs pass the $t_\mathrm{ready}$ test, but only because it compares $\mathrm{ER}$ with its own end value. That test cannot detect a relaxation that outlasts the run.
5. **At high pressure and a hot wall, the N=500 upper chromosphere is unsteady, and the unsteadiness is a resolution artefact of the TR base.** This affects 50 kK at 0.04 and 0.08 Pa, and 70 kK at 0.02, 0.04 and 0.08 Pa.
   - **What the sweep shows.** The face flux fluctuates intermittently and strongly below the TR, and at 70 kK and 0.04 / 0.08 Pa the top-cell flux reverses in 21 / 37 of 500 frames. The base stays quiet.
   - **What the follow-up at 70 kK and 0.08 Pa shows** (four runs, 1500 s each):
     - **CFL 0.25** reproduces the bursts at the same times and amplitude, so the time step is not the cause.
     - **Moving the mesh grade** from 2079–2101 km to 1931–1949 km does not move the bursts either. The unsteadiness now fills the whole fine band down to 1950 km. The grade only marks where the coarse mesh starts to damp the grid-scale noise.
     - **N=1000** cuts the peak unsteadiness by 4× (6.7 → 1.7, in units of the column flux scale) and removes every top-face reversal.
     - **The bursts lock to the grid.** The burst onsets coincide, to within 10–30 s, with the 7 kK isotherm at the foot of the TR stepping into the next cell. That foot is steep: 7.0–7.5 kK spans 0.18 km at N=500, below one 0.276 km cell.
   - **The mean flux is unaffected.** Over 1000–1500 s, all four variants give the same mean $F$ to 1.5 %.
   - **At N=1000 the corner obeys the closure like every quiet run:** $\mathrm{ER} = 0.97$–0.975 and $F L/J = 0.97$–0.975. At N=500, 300 s windows scatter by ±0.15–0.2, so the sweep's corner values (0.91–1.06) are noisy but unbiased.
   - **Why only high pressure is not established.** At 0.08 Pa, 32 % of the enthalpy jump $h_\mathrm{top} - h_c$ is gained between 6.7 and 10 kK, inside the sub-cell foot, against 13–21 % at ≤ 0.0103 Pa. There the foot is even thinner (about 0.05 km), yet it stays quiet.
6. **Resolution is otherwise fine where it was checked.** At 50 kK and 0.0103 Pa, the N=1000 control agrees with N=500 to 0.2 % in $F$ and to 0.002 in $\mathrm{ER}$.

**What this means for an AWSoM-R lookup.** The release's steady closure is usable over 18–70 kK × 0.003–0.08 Pa, to 2–10 % low, **provided $h_c$ is right**. In the corner $T_\mathrm{wall} \geq 50$ kK with $p_\mathrm{top} \geq 0.02$–0.04 Pa, the N=500 flux is correct only as an average over ≳ 500 s. A run that needs the instantaneous flux there needs N=1000, and even N=1000 is not noise-free (point 5).

A table at fixed $T_c = 6.7$ kK inherits the drift of point 3 from the release. Whether a coupled model has that drift depends on its chromospheric energy balance, and the release has none. None of this tests the coupled evaporation. At the same $(T, p)$, the AWSoM-R radiative TR is estimated to carry 3–30 times more heat flux than this column (see Limits).

Evidence:
- runs: `outputs/model_column/tp_sweep/` on the mac-mini (6.4 GB; `driver.log`, per-run console logs);
- compact products, also copied locally: `tp_sweep_compare.json`, `tp_sweep_series.npz`, `compare.log`;
- scripts: `util/run_tp_sweep.sh`, `util/tp_sweep_compare.py`, `util/plot_tp_sweep_vs_analytic.py`;
- figures: `visualization/model_column/tp_sweep_sim_vs_analytic.png`, `tp_sweep_time_series.png`, `tp_sweep_flux_compare.png` (`docs/viz_catalog/model_column.md`);
- corner follow-up: runs in `outputs/model_column/tp_sweep_corner/` on the mac-mini (`launch.sh`); `corner_test.json` copied locally; `util/tp_corner_test_compare.py`; `visualization/model_column/tp_corner_test.png`.

## Setup

- **Model.** The canonical release: N=500/R4, 661 cells, single fixed outer temperature, physical conduction only, CFL 0.50, double precision. Two things change per case:
  - `ISO_T_TOP` $= T_\mathrm{wall}$ takes 18, 25, 35, 50 and 70 kK.
  - `ISO_P_TOP` takes 0.003, 0.005, 0.0103, 0.02, 0.04 and 0.08 Pa.
- **Pressure rebuild.** `ISO_P_TOP` rebuilds the initial column with the same C7 $T(h)$. A secant iteration on the base pressure puts the initial top cell at the target, and both reservoirs are captured from the rebuilt column. The base pressure is 0.0151 Pa for 0.003 Pa ($0.22\times$ C7) and 0.822 Pa for 0.08 Pa ($12.1\times$).
- **What the construction fixes.** The initial TR has the C7 shape at every pressure, so the initial $L_\mathrm{TR}$ and $T_c$ are fixed by construction. The initial TR energy content, however, scales with the density, which is why the relaxation of point 4 lengthens with $p$.
- **Control.** One N=1000 run at 50 kK and 0.0103 Pa.
- **Runs.** 5000 s each, with snapshots every 10 s and face-flux captures every 5.6–10 s.
- **Cost.** 4 jobs × 2 OpenMP threads on the mac-mini (M2 Pro), 00:14–04:24, about 4 h 10 min for all 31 runs.

**Reading.** Each run goes through `twall_relaxation_diag.analyse`. The time $t_\mathrm{ready}$ is the later of two times:
- $F_\mathrm{base}/F_\mathrm{top}$ has stayed within 5 %;
- $\mathrm{ER}$ has stayed within 0.05 of its end value.

The ratios are 300 s means read from $t_\mathrm{ready}$. A run is marked not ready (`*`) when $t_\mathrm{ready} + 300$ s exceeds 5000 s; it is then read over the last 300 s. Every quantity is also read over 4700–5000 s (the "end" columns), because the flux keeps drifting after $t_\mathrm{ready}$ and because of point 4.

**Classes.**
- **Unsteady:** $\max_z \mathrm{std}_t F/|\mathrm{mean}_t F| > 0.5$ after 1000 s. The gap is wide: the quiet runs reach at most 0.27, always at the base, where it is just the secular decay of $F$. The unsteady runs reach 0.96–12.6, always at 2089–2108 km.
- **Relaxing:** quiet, but with $\mathrm{ER} > 1.03$ at $t_\mathrm{ready}$ and still falling at 5000 s.
- **Relaxed:** everything else. This includes 50 kK at 0.02 Pa: it is quiet and its $\mathrm{ER}$ is flat at 0.96–0.97 from 500 s on, but its base/top ringing stays near 10 %, so it never passes the 5 % wave test.

## Results

"nom" means the nominal $T_c = 6.7$ kK. "end" means the 4700–5000 s window. `*` means not ready before 5000 s.

| $T_\mathrm{wall}$ [kK] | $p_\mathrm{top}$ [Pa] | class | $t_\mathrm{ready}$ [s] | $F_\mathrm{sim}$ [kg m⁻² s⁻¹] | $L_\mathrm{TR}$ [km] | ER | ER (end) | $F_\mathrm{sim}/F_E$ nom | nom (end) | $FL/J$ nom | nom (end) | $h_c/h_c^\mathrm{nom}$ (end) |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 18 | 0.003 | relaxed | 560 | 2.546e-10 | 14.1 | 0.928 | 0.899 | 0.951 | 0.864 | 0.967 | 0.837 | 0.979 |
| 18 | 0.005 | relaxed | 1150 | 2.107e-10 | 14.0 | 0.934 | 0.895 | 0.948 | 0.878 | 0.959 | 0.862 | 0.988 |
| 18 | 0.0103 | relaxed | 2330 | 1.727e-10 | 13.4 | 0.955 | 0.913 | 0.957 | 0.909 | 0.968 | 0.905 | 0.996 |
| 18 | 0.02 | relaxed | 3310 | 1.574e-10 | 12.6 | 0.988 | 0.948 | 0.984 | 0.942 | 1.001 | 0.947 | 0.993 |
| 18 | 0.04 | relaxing | 4140 | 1.547e-10 | 12.0 | 1.042 | 1.010 | 1.036 | 1.004 | 1.069 | 1.025 | 0.989 |
| 18 | 0.08 | relaxing | 4500 | 1.731e-10 | 11.9 | 1.174 | 1.154 | 1.169 | 1.149 | 1.261 | 1.231 | 0.987 |
| 25 | 0.003 | relaxed | 384 | 6.070e-10 | 13.8 | 0.933 | 0.935 | 0.951 | 0.885 | 0.953 | 0.853 | 0.958 |
| 25 | 0.005 | relaxed | 554 | 5.667e-10 | 13.2 | 0.935 | 0.931 | 0.970 | 0.900 | 0.984 | 0.880 | 0.970 |
| 25 | 0.0103 | relaxed | 1030 | 5.025e-10 | 12.7 | 0.966 | 0.930 | 1.006 | 0.925 | 1.032 | 0.917 | 0.993 |
| 25 | 0.02 | relaxed | 2260 | 4.303e-10 | 12.6 | 0.984 | 0.944 | 1.008 | 0.957 | 1.028 | 0.959 | 1.022 |
| 25 | 0.04 | relaxing | 3430 | 4.183e-10 | 11.8 | 1.031 | 0.992 | 1.045 | 1.004 | 1.074 | 1.018 | 1.028 |
| 25 | 0.08 | relaxing | 4160 | 4.583e-10 | 11.2 | 1.128 | 1.097 | 1.138 | 1.107 | 1.208 | 1.160 | 1.029 |
| 35 | 0.003 | relaxed | 370 | 1.258e-09 | 15.5 | 0.955 | 0.957 | 0.949 | 0.892 | 0.942 | 0.855 | 0.925 |
| 35 | 0.005 | relaxed | 537 | 1.211e-09 | 14.5 | 0.952 | 0.955 | 0.954 | 0.908 | 0.951 | 0.882 | 0.938 |
| 35 | 0.0103 | relaxed | 547 | 1.260e-09 | 12.6 | 0.950 | 0.952 | 0.991 | 0.931 | 1.007 | 0.920 | 0.964 |
| 35 | 0.02 | relaxed | 2216 | 9.850e-10 | 13.6 | 0.953 | 0.951 | 0.971 | 0.949 | 0.975 | 0.946 | 0.997 |
| 35 | 0.04 | relaxed | 3440 | 9.678e-10 | 12.7 | 0.961 | 0.957 | 0.989 | 0.976 | 0.998 | 0.980 | 1.056 |
| 35 | 0.08 | relaxed | 4430 | 1.009e-09 | 11.6 | 0.992 | 0.967 | 1.020 | 0.994 | 1.034 | 1.005 | 1.108 |
| 50 | 0.003 | relaxed | 519 | 2.208e-09 | 21.4 | 0.971 | 0.971 | 0.949 | 0.895 | 0.935 | 0.852 | 0.872 |
| 50 | 0.005 | relaxed | 524 | 2.464e-09 | 17.9 | 0.970 | 0.969 | 0.956 | 0.909 | 0.948 | 0.879 | 0.888 |
| 50 | 0.0103 | relaxed | 534 | 2.774e-09 | 14.4 | 0.968 | 0.968 | 0.970 | 0.932 | 0.970 | 0.916 | 0.917 |
| 50 | 0.02 | relaxed | 4971* | 1.753e-09 | 20.0 | 0.964 | 0.964 | 0.948 | 0.948 | 0.941 | 0.941 | 0.954 |
| 50 | 0.04 | unsteady | 4973* | 2.107e-09 | 15.5 | 0.966 | 0.966 | 0.966 | 0.966 | 0.964 | 0.964 | 0.997 |
| 50 | 0.08 | unsteady | 4774* | 2.340e-09 | 13.0 | 0.952 | 0.952 | 0.965 | 0.965 | 0.961 | 0.961 | 1.069 |
| 70 | 0.003 | relaxed | 605 | 3.615e-09 | 30.1 | 0.981 | 0.983 | 0.949 | 0.897 | 0.930 | 0.852 | 0.799 |
| 70 | 0.005 | relaxed | 526 | 4.404e-09 | 23.7 | 0.979 | 0.979 | 0.956 | 0.910 | 0.945 | 0.876 | 0.822 |
| 70 | 0.0103 | relaxed | 523 | 5.394e-09 | 17.9 | 0.977 | 0.978 | 0.966 | 0.932 | 0.961 | 0.911 | 0.856 |
| 70 | 0.02 | unsteady | 5000* | 3.139e-09 | 26.9 | 0.943 | 0.943 | 0.917 | 0.917 | 0.908 | 0.908 | 0.894 |
| 70 | 0.04 | unsteady | 5000* | 4.605e-09 | 19.6 | 1.061 | 1.061 | 1.050 | 1.050 | 1.047 | 1.047 | 0.948 |
| 70 | 0.08 | unsteady | 4990* | 5.528e-09 | 15.2 | 1.042 | 1.042 | 1.037 | 1.037 | 1.035 | 1.035 | 0.966 |
| 50 | 0.0103 | N=1000 control | 533 | 2.780e-09 | 14.5 | 0.966 | 0.968 | 0.969 | 0.933 | 0.968 | 0.916 | 0.919 |

Ranges by class, at $t_\mathrm{ready}$ and at the end:

| class (runs) | ER | $F_\mathrm{sim}/F_E$ nom | $FL/J$ measured $h_c$ | $FL/J$ nom |
|---|---|---|---|---|
| relaxed (21) | 0.928–0.992 → 0.895–0.983 | 0.948–1.020 → 0.864–0.994 | 0.925–1.006 → 0.884–0.982 | 0.930–1.034 → 0.837–1.005 |
| relaxing (4) | 1.031–1.174 → 0.992–1.154 | 1.036–1.169 → 1.004–1.149 | 1.057–1.269 → 1.003–1.238 | 1.069–1.261 → 1.018–1.231 |
| unsteady (5) | 0.943–1.061 (last window only) | 0.917–1.050 | 0.944–1.064 | 0.908–1.047 |

On the C7 line (0.0103 Pa), the relaxed ratios reproduce the relaxation study's 0.95–0.98 at $t_\mathrm{ready}$.

### Why the nominal closures drift at low pressure (point 3)

The relaxed runs at $p_\mathrm{top} \leq 0.0103$ Pa end with $F L/J$ (nom) at 0.837–0.920, against 0.941–1.005 above 0.0103 Pa. The measured-$h_c$ columns barely move. Three facts in the table and in panel (c) of `tp_sweep_sim_vs_analytic.png` explain this:

1. **$h_c$ falls with time only in hot, low-pressure runs.** At 5000 s, $h_c/h_c^\mathrm{nom}$ = 0.80–0.92 at 50–70 kK and $\leq 0.0103$ Pa, against 0.95–1.11 above. These runs evaporate fastest relative to their column mass, and $L_\mathrm{TR}$ doubles between $t_\mathrm{ready}$ and 5000 s (70 kK, 0.003 Pa: 30 → 62 km). The lifted chromosphere cools adiabatically. The chromospheric gas that the TR draws on is therefore colder and less ionized than the C7 plateau.
2. **At low pressure, $h_c$ is dominated by ionization energy and is exponentially sensitive to $T_c$** (`analytic_evaporation_J_table.md`, $T_c$ sensitivity table). A few hundred kelvin of adiabatic cooling is enough.
3. **The energy closure with the measured $h_c$ does not drift.** ER stays at 0.95–0.98 at 35–70 kK. So the TR balance still holds; only the table's assumed $T_c$ is wrong.

At 18–25 kK the measured-$h_c$ ER also drifts, from 0.93–0.99 to 0.90–0.95. This is the slow late decline of $F L/J$ that the relaxation study reported at the cool walls.

### Slow TR relaxation at high pressure (point 4)

The table gives the e-folding time of $\mathrm{ER}(t) = a + b\,e^{-t/\tau}$ fitted from 1000 s. It is quoted only where $b$ is well constrained, i.e. where $\mathrm{ER}$ visibly relaxes after 1000 s.

| $T_\mathrm{wall}$ | 0.003 Pa | 0.005 Pa | 0.0103 Pa | 0.02 Pa | 0.04 Pa | 0.08 Pa |
|---|---:|---:|---:|---:|---:|---:|
| 18 kK | 370 s | 650 s | 880 s | 1100 s | 1850 s | 3450 s |
| 25 kK | – | – | 440 s | 820 s | 1190 s | 2100 s |
| 35 kK | – | – | – | – | 630 s | 1270 s |

At 18 kK and 0.08 Pa, $\tau$ is comparable to the fit span. Extrapolating the fit to $t \to \infty$ gives $\mathrm{ER} \to 0.79$ but $F L/J \to 0.96$. The two should agree to about 0.02 in a steady TR, so neither asymptote is constrained, and none is quoted as a result. What can be stated: every relaxing run is monotonically falling toward the relaxed band, so its tabulated value is an upper bound.

The trend has a plausible reading, which was not checked quantitatively:
- The initial TR has the C7 shape at every pressure, so its mass and stored enthalpy scale with $p_\mathrm{top}$.
- The conductive drive $q_\mathrm{wall}$ that drains it depends mainly on $T_\mathrm{wall}$.

Reading off $\tau$ at those points therefore needs a run of order 10 ks, or an initial TR built in conduction–enthalpy balance.

### The unsteady hot, high-pressure corner (point 5)

**In the sweep** (`tp_sweep_time_series.png` and the face-flux captures):

- **Where it is.** Normalised by its local mean, the unsteadiness grows with height from a quiet base and peaks at 2089–2108 km. The TR base (7.5 kK) is at 2138–2143 km at 70 kK and 0.08 Pa. The peak location is partly an artefact of that normalisation (below, a single column-wide scale is used). The coincidence with the R4 grade (2079–2101 km), which first suggested the mesh, turned out not to be causal.
- **How strong it is.** At 70 kK and 0.08 Pa, local face fluxes reach $-1.6\times10^{-7}$ kg m⁻² s⁻¹, 25 times the mean and reversed. The top face reverses in 75 of 893 captures. The dominant periods of the top-face flux after 1000 s are about 23 s, 167 s and 364 s; the last is the acoustic-cutoff period of the column.
- **When it starts.** The first top-face reversal at 70 kK and 0.08 Pa is at 572 s. The disturbance is intermittent, with bursts of order 100 s separated by quiet stretches: 50 kK at 0.08 Pa is quiet from 4500 s to the end.
- **The TR balance survives it.** $\mathrm{ER}(t)$ is noisy but centred on 1, and the 300 s window means of the ratios lie at 0.91–1.06. The TR responds to the disturbance quasi-statically; the disturbance is not the TR failing.
- **The earlier CFL check could not see it.** The CFL 0.25 against 0.50 comparison at this corner (top mass flux within $2\times10^{-4}$) covered only 200 s, before the first burst.
- **A Schwarzschild check is inconclusive.** $\nabla > \nabla_\mathrm{ad}$ appears only in frames that are already disturbed, at 2070–2140 km, and never in the quiet frames or the quiet runs. In those disturbed frames, $d\ln p$ is itself perturbed, so any excess is an effect of the disturbance. The follow-up below makes a convective origin moot: a physical instability would not weaken 4× under refinement.

**Follow-up: resolution, CFL and mesh grade.** The case 70 kK, 0.08 Pa was repeated for 1500 s in three variants. The baseline is the first 1500 s of the sweep run. The bursts start by about 500 s, so 1500 s covers seven of them. The table gives statistics over 500–1500 s. $S = \mathrm{median}_z\,|\mathrm{mean}_t F(z)|$ is one flux scale for the whole column; a local mean can pass through zero and blow up a local ratio.

| case | cells | fine cell [km] | grade [km] | $\max_z \mathrm{std}_t F/S$ | at [km] | base | top-face $F \leq 0$ | $\min F/S$ | mean $F$, 1000–1500 s [kg m⁻² s⁻¹] | $\mathrm{std}_t F_\mathrm{top}/\bar F$ | ER, 1000–1500 s |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| base: N500, CFL 0.50 | 661 | 0.276 | 2081–2099 | 6.7 | 2138 | 0.15 | 9.0 % | −27 | 6.892e-9 | 0.78 | 0.985 ± 0.79 |
| CFL 0.25 | 661 | 0.276 | 2081–2099 | 7.2 | 2101 | 0.15 | 9.2 % | −25 | 6.909e-9 | 1.06 | 0.984 ± 1.02 |
| grade at 1950 km (`ISO_REFINE_S_LO_KM=350`) | 1068 | 0.276 | 1931–1949 | 7.2 | 1944 | 0.16 | 5.0 % | −33 | 6.971e-9 | 1.00 | 1.007 ± 1.01 |
| N1000/R4 | 1320 | 0.138 | 2081–2100 | 1.7 | 2132 | 0.20 | 0 | −8 | 6.883e-9 | 0.26 | 0.975 ± 0.25 |

All four runs end with `termination=end_time` and zero fallbacks. What the variants show:

- **The time step is not the cause.** CFL 0.25 matches the baseline burst for burst. The top-face bursts start at 496, 572, 676, 796, 950, 1130 and 1348 s, against 505, 572, 678, 802, 959, 1133 and 1340 s in the baseline.
- **The grade is not the cause, but it bounds the noise.** With the grade moved down by 150 km, the bursts keep the same times (510–1335 s), and the TR is untouched. The large-amplitude noise now fills the whole fine band: its unsteadiness plateau is about 6 from 1950 km up, where the baseline plateau only begins near 2080 km. In both N500 meshes it decays inside the coarse cells. The noise is therefore grid scale: made at the TR, carried by the fine mesh, and damped by the coarse mesh.
- **Resolution controls it.** Halving every cell cuts the peak unsteadiness by 4× and the top-face fluctuation by 3×, and removes the reversals. It does not remove the noise: 1.7 is still well above the ≤ 0.27 of the quiet runs.
- **The bursts are locked to the TR foot crossing cells.** Two facts show this:
  - **Timing.** The 7 kK isotherm, which moves *up* by 0.7 km over 300–1500 s while the 7.5 kK isotherm moves down by 2.7 km, enters a new cell at 800, 970, 1150 and 1370 s. The baseline bursts start at 802, 959, 1133 and 1340 s. The 7.5 kK crossings (540, 630, 750, 900, 1070, 1270 s) do not match. At N1000 the cells are half as wide, and the crossings and the bursts are about twice as frequent but much weaker.
  - **Sub-cell width.** The foot is under-resolved. By 1500 s, 7.0–7.5 kK spans 0.18 km at N500 and 0.12 km at N1000, against cells of 0.276 and 0.138 km. The steady-TR integral, $\Delta s = \int_{7.0\,\mathrm{kK}}^{7.5\,\mathrm{kK}} \kappa\,dT/[F(h - h_c)]$, predicts 0.13–0.14 km.
- **The mean flux is robust.** The 1000–1500 s mean $F$ agrees across all four runs to 1.5 %. At N1000 the closure ratios sit where the quiet 70 kK runs do ($\mathrm{ER} = 0.975$). The noise is zero-mean in a 500 s average, although 300 s windows at N500 scatter by ±0.15–0.2.

**What stays open.** First, it is not known why only high pressure is affected. By the same integral, every ≥ 50 kK run has a sub-cell foot, and the quiet low-pressure ones have the thinnest (about 0.05 km at 70 kK and 0.0103 Pa). What grows with pressure is the share of $h_\mathrm{top} - h_c$ gained between 6.7 and 10 kK, i.e. inside the unresolved foot (evaluated at 70 kK):
- 0.13 at 0.003 Pa;
- 0.21 at 0.0103 Pa;
- 0.25 at 0.02 Pa;
- 0.32 at 0.08 Pa.

That share is the ionization energy of a chromosphere-top gas that is only 29 % ionized at 0.08 Pa, against 84 % at 0.003 Pa. It is a plausible control parameter, but it has not been tested. Second, the noise has not been shown to vanish under further refinement: N2000 was not run. Third, a TR-broadening treatment such as TRAC would address exactly this kind of under-resolved front. TRAC is not release physics, however, and adding it is a design decision, not a fix.

## Limits

- **Release physics only.** There is no radiation, no chromospheric heating, and $h_c$ comes from Saha LTE. The $T_c$ drift of point 3 is a property of this energy balance.
- **Initial TR by construction.** At every pressure the initial $L_\mathrm{TR}$ and $T_c$ are the C7 values. The high-pressure relaxation of point 4 is partly an artefact of starting from a TR that is out of balance at that pressure.
- **Not a coupled test.** The AWSoM-R lower boundary is a static conduction–radiation TR with a prescribed length of about 1 Mm (Sokolov et al. 2021, Eqs. 19–22, `docs/supporting-papers/Sokolov_2021_ApJ_908_172.pdf`). At the same $(T, p)$, its heat flux was estimated at 3–30 times this column's when the sweep was designed. That estimate is not recorded elsewhere in the repository and was not re-derived here. The sweep tests the release's closure, not the coupled evaporation rate.
- **Resolution.** Resolution was checked at one mid-grid point and at the 70 kK / 0.08 Pa corner (1500 s). At the corner, N500 carries grid-scale noise at the unresolved TR foot but the same 500 s mean flux (point 5). The N500 long-duration evaporation mass flux is still not established as converged to better than 10 %.
- **CFL.** At the extreme corner, CFL 0.50 and 0.25 agree in the 1000–1500 s mean flux to 0.2 % and reproduce the same bursts.
- **No plateau.** Every $F$ here drifts 1–6 % per 500 s at the end, consistent with $F = J/L_\mathrm{TR}$. Report $J$ or the ratios together with their reading time, not $F$ alone.

## Reproduce

On the mac-mini (the runs are 6.4 GB and stay there):

```bash
SWEEP_JOBS=4 SWEEP_THREADS=2 util/run_tp_sweep.sh outputs/model_column/tp_sweep 5000 10
.venv/bin/python util/tp_sweep_compare.py --dir outputs/model_column/tp_sweep \
    --out-dir visualization/model_column --workers 8 > outputs/model_column/tp_sweep/compare.log
```

Anywhere `tp_sweep_compare.json` and `tp_sweep_series.npz` are (about 1.6 MB together):

```bash
.venv/bin/python util/plot_tp_sweep_vs_analytic.py --dir outputs/model_column/tp_sweep
```

The sweep takes about 4 h 10 min on an M2 Pro, the comparison about 12 s with 8 workers, and the plots about 3 s.

Corner follow-up, on the mac-mini. It uses 4 + 2 + 2 threads and takes 30 min; the N1000 run is the slowest.

```bash
outputs/model_column/tp_sweep_corner/launch.sh
.venv/bin/python util/tp_corner_test_compare.py
```

`launch.sh` runs `util/run_tp_sweep.sh` three times. Each call has `SWEEP_CASES` set to the 70 kK, 0.08 Pa case at N=1000 or N=500, plus `CHROMO_CFL=0.25` or `ISO_REFINE_S_LO_KM=350` where the variant needs it, with `T_END=1500` and `FRAME_DT=10`.
