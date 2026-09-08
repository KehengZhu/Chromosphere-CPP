# Upper-boundary conduction-wall temperature sensitivity of `model_column`

**Question.** With the C7 initial condition and every other piece of release physics and numerics held fixed, how does the field-aligned velocity respond to the imposed conduction-wall temperature $T_{\rm wall}$ at the 2153 km outer face?

**Answer.** Strongly nonlinear, with a flow reversal near $T_{\rm wall}\approx 10$ kK and **no saturation anywhere in the range tested**. On the upflow branch $V_{\rm top}$ follows an approximate power law $V\propto T_{\rm wall}^{3.1}$ whose *local* exponent falls monotonically from $\approx 12$ just above the reversal to $2.6$ at 100 kK. Over $10\to100$ kK the top-of-domain velocity spans $-0.10 \to +792$ m/s, i.e. **four decades of velocity for one decade of wall temperature**. The controlling quantity is the outer-face conductive flux $q_{\rm wall}$, and the response is quantitatively the Antiochos & Sturrock (1978) conduction-driven evaporation balance: the release has no radiative sink and no coronal heating, so a quasi-steady column must export $q_{\rm wall}$ as the enthalpy flux of the evaporating material.

Evidence: `outputs/model_column/twall_sweep/`, metrics in `twall_sweep_metrics.json`, figures `visualization/model_column/twall_sweep_N500_4000s_{V_of_t,V_of_h,summary,support}.png`, and `visualization/model_column/twall_sweep_N500_4000s_overlay.mp4` — all 13 accepted walls overlaid on the same six panels as the release evolution movie (`util/animate_twall_sweep.py`; catalogued in `visualize_commands.md` §16).

## What was varied, and what was not

`ISO_T_TOP` with `ISO_HYDRO_T_DECOUPLE=1` (both release presets) moves exactly one thing: `Grid::outer_conduction_temperature`, the fixed thermal datum the implicit conduction stage sees at the **physical** outer face. The hydro ghost temperature remains a zero-gradient extrapolation of the live top cell, so the outer face still imposes only its one legitimate external condition — the reservoir back-pressure `kOuterPRef`. That back-pressure is captured from the same C7 IC top cell in every case (`kOuterPRef = 0.0102856 Pa`, identical across the sweep, confirmed in every console log), so the sweep really is one-parameter.

Everything else is the canonical release configuration: N=500/R4 (661 cells), `h_base = 1600` km, `dh = 553` km, production `Gamma1` table, SWMF exact-Riemann Godunov flux on $(\ln\rho, V, \ln p)$ with the MC3 limiter, physical conduction only, `CHROMO_CFL=0.50`, 4000 s, **double-precision state** (`build_omp`, `CHROMO_STATE_FLOAT32=OFF`). Every accepted run reported `termination=end_time` and `godunov.fallbacks=0`.

Driver: `util/run_twall_sweep.sh`. Analysis and figures: `util/twall_sweep_diag.py`.

"Quasi-steady" below means the mean over the last 500 s of the 4000 s run. That is not a true steady state — see *Time convergence* — so it is reported with its drift.

## Result

| $T_{\rm wall}$ [kK] | $V_{\rm top}$ [m/s] | $V$@2100 km | $V$@1800 km | $(\rho V)_{\rm top}$ [kg m⁻² s⁻¹] | $q_{\rm wall}$ [W m⁻²] | $T_{\rm wall}-T_{\rm top}$ [K] | $\kappa_{\rm wall}$ | drift / 500 s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8   | **−0.435** | −0.212 | −0.0735 | −3.41e−11 | −0.0012 | −2.3 | 0.070 | −21.5 % |
| 10  | **−0.104** | −0.041 | −0.0145 | −6.52e−12 | +0.0117 | 17.6 | 0.092 | −62.1 % |
| 11  | +0.180 | 0.064 | 0.0214 | +1.02e−11 | 0.0211 | 24.9 | 0.117 | +25.1 % |
| 12  | +0.540 | 0.176 | 0.0598 | +2.81e−11 | 0.0330 | 31.3 | 0.145 | +4.2 % |
| 13  | +0.976 | 0.293 | 0.100 | +4.69e−11 | 0.0475 | 36.9 | 0.177 | −0.1 % |
| 14  | +1.495 | 0.417 | 0.143 | +6.68e−11 | 0.0649 | 42.0 | 0.214 | −1.7 % |
| 15  | +2.104 | 0.548 | 0.188 | +8.77e−11 | 0.0856 | 46.6 | 0.254 | −2.4 % |
| 18  | +4.561 | 0.991 | 0.339 | +1.58e−10 | 0.169 | 58.2 | 0.400 | −3.0 % |
| **22 (release)** | **+9.645** | 1.717 | 0.587 | +2.74e−10 | 0.338 | 70.5 | 0.661 | −3.4 % |
| 30  | +27.62 | 3.613 | 1.234 | +5.75e−10 | 0.908 | 87.4 | 1.435 | −4.3 % |
| 45  | +93.51 | 8.207 | 2.785 | +1.30e−09 | 2.893 | 101.0 | 3.954 | −5.7 % |
| 70  | +310.2 | 17.70 | 5.970 | +2.77e−09 | 9.209 | 106.5 | 11.93 | −6.4 % |
| 100 | +791.9 | 33.07 | 10.94 | +4.94e−09 | 23.01 | 109.1 | 29.11 | −6.3 % |

The 22 kK row reproduces the documented release result exactly: its final-instant top velocity is **9.50082 m/s** against the 9.50 m/s recorded for the double-precision cutover run. The sweep baseline is the release, not a lookalike.

### Shape of the response

- **Flow reversal.** $V_{\rm top}$ crosses zero at $T_{\rm wall}\approx 10.4$ kK (linear interpolation between the 10 and 11 kK points). Below it the whole column drains downward; above it the whole column evaporates upward. $q_{\rm wall}$ crosses zero earlier, at $\approx 8.2$ kK.
- **Not linear, not saturating.** The local logarithmic slope $n = d\ln V/d\ln T_{\rm wall}$ falls monotonically: 12.6 (11→12 kK), 5.75 (13→14), 4.24 (15→18), **3.73 (18→22, i.e. at the release point)**, 3.39 (22→30), 3.01 (30→45), 2.71 (45→70), 2.63 (70→100). A single power law fitted over $T_{\rm wall}\ge 15$ kK gives $V\propto T_{\rm wall}^{3.09}$. The softening is a gradual change of exponent, **not** an approach to a plateau: $V$ is still rising by a factor 2.6 per factor 1.43 in $T_{\rm wall}$ at the hot end.
- **Threshold-like, not a true threshold.** The apparent threshold at 10 kK is just where a smooth curve crosses zero; the enormous local exponent there is the artefact of taking $d\ln V$ through a sign change.
- **The response is a mass-flux response.** $\rho V$ is height-independent to within 14 % across the entire 553 km column in every case (top/1800 km ratio 1.05–1.14). $V(h)$ rises steeply with height only because $\rho(h)$ falls; the physical response variable is the mass flux, and $V$ at any given height is that mass flux divided by the local density.
- **Only the top ~40 km moves.** The 6.7 kK C7 chromospheric plateau is unchanged by $T_{\rm wall}$ at every value tested; the entire temperature response is confined to the resolved TR above ~2115 km (`..._support.png`, middle panel).

## What controls the velocity — physical interpretation

Three links, each measured:

**1. The wall thermostats the top cell.** $T_{\rm top}$ tracks $T_{\rm wall}$ to within 0.1–0.3 % at every wall temperature (109 K at $T_{\rm wall}=100$ kK). The discrete Dirichlet wall sits at the physical face, half a refined cell (138 m) from the top cell centre, where $\kappa$ is large, so the top cell simply equilibrates to the wall. The temperature drop $T_{\rm wall}-T_{\rm top}$ therefore grows only as $\approx T_{\rm wall}^{0.8}$ — this is the self-limiting behaviour the boundary was designed to have.

**2. That fixes the conductive flux, and Spitzer supplies the nonlinearity.** $q_{\rm wall} = \kappa_{\rm wall}(T_{\rm wall}-T_{\rm top})/\Delta s_{\rm face}$ with $\kappa_{\rm wall}\propto T_{\rm wall}^{2.50}$ (measured exponent 2.50 over $10\to100$ kK, i.e. exactly Spitzer). Combined with the sub-linear $\Delta T$, $q_{\rm wall}\propto T_{\rm wall}^{3.3}$ — this is where essentially all of the nonlinearity comes from. Nothing in the hydrodynamics is nonlinear here: the flow is deeply subsonic (quasi-steady Mach number $\le 0.015$ everywhere, even at 100 kK).

**3. A quasi-steady column exports that flux as an enthalpy flux.** The release has no radiative sink and no volumetric heating, so in steady state the conducted energy has nowhere to go but out of the top with the evaporating material:

$$q_{\rm wall}\;\simeq\;\rho_{\rm top} V_{\rm top}\left(\frac{\gamma}{\gamma-1}\frac{p_{\rm top}}{\rho_{\rm top}} + x\,\frac{\chi_H}{m_H}\right).$$

Measured against this one-line estimate, $V_{\rm measured}/V_{\rm predicted}$ = 1.97, 1.92, **1.79**, 1.61, 1.42, 1.26, 1.16 at 15, 18, 22, 30, 45, 70, 100 kK — agreement to better than a factor 2 across four decades of velocity, tightening toward unity as the column approaches steady state at the hot end. The residual excess is the part of the mass flux still being supplied from the column's own stored mass and internal energy rather than from $q_{\rm wall}$ (see *Time convergence*). The two enthalpy terms are comparable at the release point ($9.05\times10^8$ thermal versus $1.30\times10^9$ J kg⁻¹ ionization at 22 kK) and the thermal term dominates by 100 kK ($4.12\times10^9$); the hydrogen ionization reservoir is a first-order part of the evaporation energy budget in the chromospheric range, not a correction.

**Why $V$ is steeper than $\rho V$.** The fixed reservoir back-pressure pins $p_{\rm top}$, so $\rho_{\rm top}\propto 1/T_{\rm wall}$ (measured: $6.24\times10^{-11}\to6.24\times10^{-12}$ kg m⁻³ over $10\to100$ kK, exactly a factor 10) — the constant-pressure transition-region relation. The mass flux scales as $T_{\rm wall}^{2.8}$ and $V_{\rm top}=(\rho V)/\rho_{\rm top}$ picks up one extra power, giving $T_{\rm wall}^{3.8}$ over the same interval.

**The reversal.** Below $T_{\rm wall}\approx 10$ kK the imposed wall is no longer hotter than the material immediately beneath the top cell (the C7 IC is 10.2 kK at 2146 km and 7.5 kK at 2140 km), so the boundary stops driving heat into the TR, the top of the column cools and contracts, and the fixed coronal back-pressure above pushes material back down: $\rho V<0$ at every height. The reversal is a **change of sign of the boundary's thermal forcing relative to the underlying chromosphere**, not a dynamical instability.

## Limits of validity and caveats

- **Time convergence.** The runs are quasi-steady, not steady. Over the last 500 s, $V_{\rm top}$ is still falling by 2–6 % per 500 s for every $T_{\rm wall}\ge13$ kK. The quoted magnitudes are therefore mild **upper bounds**; the *trend* with $T_{\rm wall}$ is robust because the drift is nearly the same fraction for every hot case. Near the reversal (8–11 kK) the fractional drift is large only because $V\to0$; the absolute drift there is $\lesssim 0.1$ m/s. **The reversal temperature is bracketed by the data as 8–11 kK; the 10.4 kK interpolation is not converged in time.**
- **Cold-wall limit: $T_{\rm wall}=6$ kK is outside the model.** That run aborts at $t\approx380$ s with `Gamma1 table query out of bounds: T=3199.86 K` — the top of the column expands and cools below the EOS table floor of 3200 K. It is excluded from all figures and metrics. The coldest usable wall is 8 kK.
- **Hot-wall numerics: 100 kK required `CHROMO_CFL=0.25`.** At CFL 0.50 the implicit conduction Newton iteration fails to converge on the very first step (`mixture conduction Newton did not converge`); the threshold is between 0.30 (converges) and 0.35 (fails). This is a solver-conditioning limit, not a physical one. A **matched CFL control at 70 kK** joins the two groups: $V_{\rm top}$ = 310.212 m/s at CFL 0.50 versus 310.204 m/s at CFL 0.25, a relative difference of $2.6\times10^{-5}$. The 100 kK point is therefore directly comparable with the CFL 0.50 group. **This does not change the release default, which remains CFL 0.50 for the validated N=500/R4 configuration.**
- **The outflow Mach cap never binds.** `ISO_VCAP = 0.1` caps the outer *ghost* velocity at $0.1 c_s$. The quasi-steady top-cell $|V|$ reaches at most 15 % of that cap (100 kK); the transient peak reaches 63 % during the first ~100 s of the 100 kK run and never touches it. The high-$T_{\rm wall}$ softening of the exponent is physics, not the boundary clipping.
- **Every `rho V` number above is the cell-centred product, not the conservative face flux.** The sweep was later re-run with the read-only face-flux diagnostic (`SWEEP_FACE_FLUX=1`) to add `.faceflux` sidecars for the overlay movie; all 13 `.gamma_diag` files came back **byte-identical**, so nothing in this recap changed, and the re-run stands as an exact reproducibility check of the whole sweep. Comparing the two measures at t = 4000 s: they agree at the top face to 0.1 % at every wall (so the headline top-of-domain numbers are unaffected), but `rho*V[0]/f_total[0]` = **1.268, 1.275, 1.282, 1.292** at 8, 22, 70, 100 kK — the known first-interior-face artifact, reproducing the 1.28 of the release run and **essentially independent of $T_{\rm wall}$**, which is further evidence that it is a ghost/reconstruction closure issue rather than something the boundary drives. In the interior the discrepancy grows with the wall, from 0.96 % at 8 kK to 13.7 % at 100 kK, so the "height-independent to within 14 %" statement above should be read as a property of the cell-centred product; prefer `f_total` for any quantitative interior mass-flux claim at high $T_{\rm wall}$.
- **Grid convergence was not tested here.** The known `model_column` limitation — N=500 long-duration evaporation mass flux is not established as <10 % grid-converged — applies unchanged to every number in this study. The sweep is a controlled one-parameter comparison at fixed resolution; the *relative* response is far better determined than the absolute magnitudes.
- **Single-parameter caveat.** Real transition-region physics would not change $T_{\rm wall}$ at fixed back-pressure: raising the coronal temperature at fixed coronal density raises the pressure too. This study deliberately holds `kOuterPRef` at its C7 value so that only the thermal boundary condition moves.

## Reproduce

```bash
# broad log sweep + refinement + the CFL control (about 3 h wall on 16 cores)
util/run_twall_sweep.sh                                     # 10,15,18,22,30,45,70,100 kK
SWEEP_CASES="6000:0.50 8000:0.50 11000:0.50 12000:0.50 13000:0.50 14000:0.50" \
    util/run_twall_sweep.sh                                  # reversal refinement

# figures + metrics
.venv/bin/python util/twall_sweep_diag.py --cfl 0.50 \
    --run 100000=outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt \
    --heights 1800 2000 2100 --avg-window 500 \
    --out-dir visualization/model_column --prefix twall_sweep_N500_4000s \
    --json outputs/model_column/twall_sweep/twall_sweep_metrics.json
```
