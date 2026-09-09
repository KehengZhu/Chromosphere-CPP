# Two-temperature (T_e != T_i) test study on the release column

**Status: complete. EXPERIMENTAL, non-release.** This records a test study, not a release
capability. The release single-fluid solver is untouched and byte-identical with or without
`src/two_temp/` compiled in.

**Question.** The release is a single-fluid equilibrium mixture with one temperature. How far do
the electron and heavy-particle temperatures actually drift apart in the modelled 1600–2153 km
column, and is the single-temperature assumption safe for the evaporation observables the model
exists to measure?

**Answer.** They barely drift at all. The worst decoupling anywhere, at any step, is
`max |T_e - T_i| / T_i = 3.28e-3`, and it lasts tens of milliseconds during the boundary
switch-on. In quasi-steady state the two temperatures agree to about one part in `1e5`. The
two-temperature run reproduces the release top-of-domain velocity and mass flux to **0.18 %**.
For this configuration the release's single-temperature assumption is validated.

---

## 1. What was built

A third solver directory, `src/two_temp/`, with conserved state `U = (rho, rho u, E, E_e)`,
selected only by the new scenario `model_column_2t`. Full formulation: `src/two_temp/two_temp.hpp`.
Narrative documentation: @ref two_temp_physics (physics), @ref numerics-two-temp (numerics and the
boundary-condition derivation), @ref scenario-model-column-2t (the scenario).

`E` keeps its release meaning — the total energy — so mass, momentum and total energy stay in
conservation form and the release numerical flux applies to those three rows unchanged. The heavy
internal energy is derived, which makes `T_i` purely algebraic. The electron pool owns the
ionization reservoir and Saha is evaluated at `T_e`, because hydrogen ionization here is an
electron-impact process.

**The comparison is controlled by construction.** `model_column_2t` takes exactly the release
model and grid defaults and shares one IC/BC implementation with `model_column`, so the temperature
split is the only variable. It presets the identical `ISO_*` list, runs at the same CFL 0.50, and
rejects `ISO_GAMMA`, `ISO_RIEMANN` and `ISO_RECONSTRUCTION`.

## 2. Boundary conditions, from the characteristics

This was the design question. The hyperbolic and parabolic stages must be counted separately.

**Hyperbolic.** In primitive variables `(rho, u, p, p_e)` the electron energy equation reduces to
the adiabatic electron-pressure law `d_t p_e + u d_s p_e + (5/3) p_e d_s u = 0`, and the flux
Jacobian has characteristic polynomial `mu^2 (mu^2 - c^2)` with `mu = u - lambda`,
`c = sqrt(5/3 p/rho)`. Eigenvalues `{u - c, u, u, u + c}`; the doubled `u` is **not** defective —
its eigenspace is spanned by the entropy/density mode and by a pure electron-heavy partition mode.

At the outer face the flow is subsonic outflow, so three characteristics leave and only `u - c`
enters: **exactly one** condition may be imposed, and the release already spends it on the fixed
reservoir back-pressure. **Both hydrodynamic temperatures must therefore free-float**, and are
zero-gradient extrapolated. At the 1600 km base three characteristics are incoming, so the
Dirichlet reservoir in density and both temperatures with `v = 0` is within budget.

**Parabolic.** Conduction is a separate operator needing one condition per channel per end:

| Face | Electron channel | Heavy channel |
| --- | --- | --- |
| Outer 2153 km | Dirichlet `T_e = 22,000 K` | Neumann (zero flux) |
| Inner 1600 km | Dirichlet `T_0` | Dirichlet `T_0` |

The outer electron Dirichlet is the physically meaningful thing the split buys: conduction down
from the TR/corona is electron-conducted, and the release necessarily applies that reservoir to one
lumped temperature. The heavy channel is insulating because the modelled heavy conductivity is
*neutral-hydrogen* conduction and there are no neutrals above the domain. The base is Dirichlet on
both because `tau_eq ~ 0.09 ms` there — a property of the reservoir, not an added assumption.

`TT_OUTER_TI=dirichlet` runs the controlled alternative (the 22 kK wall on both channels).

## 3. Runs

All at N=500/R4, 661 cells, CFL 0.50, `double` storage, 4000 s, production `Gamma1` table.
All terminated `termination=end_time` at exactly t = 4000 s.

| Run | Output |
| --- | --- |
| Two-temperature, outer `T_i` Neumann (the physical case) | `outputs/model_column_2t/twotemp_N500_4000s_Ti_neumann.txt{,.twotemp}` |
| Two-temperature, outer `T_i` Dirichlet (controlled comparison) | `outputs/model_column_2t/twotemp_N500_4000s_Ti_dirichlet.txt{,.twotemp}` |
| Matched release reference, current binary | `outputs/model_column/release_N500_4000s_2tref.txt{,.gamma_diag}` |
| Dense-cadence startup transient, 0.02 s frames over 5 s | `outputs/model_column_2t/twotemp_transient_5s.txt{,.twotemp}` |

```bash
CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 CHROMO_T_END=4000 \
scripts/run_chromo_realtime.sh outputs/model_column_2t/twotemp_N500_4000s_Ti_neumann.txt \
    full no-ionization model_column_2t - 20.0 no-cooling
# TT_OUTER_TI=dirichlet for the comparison leg; scenario model_column for the reference.
```

**Run integrity, both two-temperature legs.** `godunov_fallbacks = 0` of `9.413e8` face solves,
`clamps = 0` (face_pe 0, decode_Te 0, decode_Ti 0), `conduction_iterations_max = 4`.

**Measured CFL equivalence.** 2T took 712031 steps, the release 712029 — a difference of 2 steps in
712k, confirming the prediction that the frozen signal speed `sqrt(5/3 p/rho)` depends only on the
total pressure, so splitting the temperature does not change `dt`.

## 4. Results

### 4.1 Exact degeneracy at `T_e = T_i` (control)

At t = 0 the two-temperature state and the release state are **bit-identical** in all three release
rows: maximum absolute difference exactly `0.0` in `rho`, `rho u` and `E` across all 661 cells.
`E_e` spans 0.0435–0.367 J m^-3. This is the strongest available check that the closure reduces to
the release closure exactly, and it is what makes the rest of the comparison meaningful.

### 4.2 The decoupling is small and transient

- **Global maximum over every cell and all 712031 steps: `max |T_e - T_i|/T_i = 3.28e-3`.** The
  identical value appears in a 5 s smoke run and in the 4000 s run, which proves the global maximum
  occurs within the first 5 s.
- The dense-cadence run localizes it: the peak is at **t ~ 0.022 s in the topmost cell**, where
  `T_e = 22339.3 K` sits **63.4 K below** `T_i = 22402.8 K`. The sign is diagnostic — the IC top
  cell is at 22674.8 K and the imposed electron wall is *colder* at 22,000 K, so the electron
  channel alone conducts heat out of the top cell while the heavy channel is insulated. It decays
  roughly as `1/t`: 43.1 K at 0.044 s, 18.0 K at 0.10 s, 5.0 K at 0.26 s, 1.7 K at 0.5 s.
- **Quasi-steady:** `2.05e-5` at 20 s falling to `9.01e-6` at 4000 s, i.e.
  `max |T_e - T_i| = 0.0715 K` at 2138.4 km.

**Final-frame structure.** `T_e - T_i` is smooth and monotone with one physically meaningful sign
change, at 2137.8–2138.4 km, exactly where `x -> 1`. Below it, through the partial-ionization zone
(`x` rising 0.62 to 0.99), `T_e < T_i` by up to `1.1e-2 K` — the electron pool is paying for
ionization. Above it, fully ionized, `T_e > T_i` by up to `7.1e-2 K` — electrons are conductively
heated from the boundary and hand energy to the heavies.

### 4.3 Why: timescale separation

`tau_eq` runs from `9.25e-5 s` at the 1600 km base to `3.97e-3 s` in the top cell, against
`dt = 5.6e-3 s` and a conduction/evaporation timescale of order `1e3 s`. Five to seven orders of
magnitude separate the collisional coupling from the physics of interest, so collisions enforce a
common temperature. That is the whole mechanism.

Note the top cell is the interesting one: `tau_eq/dt = 0.70` there, so the exchange is only
*marginally* resolved rather than stiff, which is precisely where a real decoupling is physically
plausible — and it is where the transient peak appears.

### 4.4 Two-temperature versus the release

| Quantity at 4000 s | Two-temperature | Release | Difference |
| --- | --- | --- | --- |
| top velocity | 9.4841 m/s | 9.5008 m/s | **-0.176 %** |
| top mass flux `rho v` | 2.69498e-10 | 2.69974e-10 kg m^-2 s^-1 | **-0.176 %** |
| top density | 2.84157e-11 | 2.84159e-11 kg m^-3 | -0.0006 % |
| total pressure, whole column | — | — | max rel diff **3.11e-5** |
| temperature | — | — | max abs diff 72.91 K at 2138.7 km |

The 72.91 K (0.80 %) temperature difference is **not** a bulk disagreement: it is localized at the
ionization front / TR onset where the vertical temperature gradient is extremely steep, so a
sub-cell displacement of the front produces a large local difference. The column-wide pressure
agreement of `3e-5` shows the mechanical state is essentially identical. The same front
displacement explains the only other outlier, `kappa_e + kappa_i` versus the release `kappa`
(max rel diff `8.4e-2`, at 2138.1 km): at the formula level the identity is exact, because
`two_temp_kappa_e`/`two_temp_kappa_i` are the release `physical_kappa_e`/`physical_kappa_n`
expressions verbatim (`n_i = n_e` for pure hydrogen).

### 4.5 The outer heavy-channel BC is immaterial — and that is the finding

Neumann versus Dirichlet on the outer `T_i` agree **to round-off**: max relative difference
`4.7e-10` in `T_e`, `8.7e-10` in `T_i`, `3.5e-10` in `rho`, `1.0e-9` in `v`, and exactly `0` in
`p_total`. Top-cell velocity is 9.4841 m/s in both.

The reason is clean: at the 22,000 K wall the gas is fully ionized, `n_HI -> 0`, and the
neutral-hydrogen conductivity collapses — `kappa_e/kappa_i = 7.6e7` in the top cell, so the heavy
channel carries about `1.3e-8` of the electron channel. There is no heavy conductive flux for the
boundary condition to control. The Neumann choice is therefore not merely defensible but
demonstrably immaterial, and the electron channel is the only one that matters at the top. This is
a positive robustness result for the electron-conducted TR closure.

## 5. Two bugs found and fixed during the study

**1. The `E_e -> T_e` inversion could fail to converge and abort a run.** Observed at
`rho = 8.26e-11`, `x = 0.999965`, which killed a 4000 s run at t ~ 56 s. Near full ionization the
caloric curve has a very sharp knee: at `x -> 1` the electron capacity collapses to the
translational value (`1.03e-6 J m^-3 K^-1`) while just below the knee `dx/dT_e` makes it orders of
magnitude larger. A Newton step from the flat side overshoots the knee by ~5500 K, and the original
safeguard — which bisected only when the candidate *left* the bracket — never fired, because each
iterate landed strictly inside by ~`1e-10 K`. The bracket shrank by ~`1e-10 K` per pass and the
100-iteration limit was exhausted. Replaced with safeguarded Newton with **guaranteed bisection
progress** (Numerical Recipes `rtsafe`), which also bisects whenever Newton fails to halve the
bracket. Verified behaviour-preserving where the old path worked (identical 5 s smoke results).

**2. `TT_OUTER_TI=dirichlet` did not do what it documented.** `outer_conduction_wall` substituted
the imposed wall temperature into the electron slot only and left the heavy slot at the
*extrapolated* ghost `T_i`, so the "controlled comparison" was a near-copy of the Neumann case and
a null difference would have proved nothing. Fixed the code to match the documented intent (the
wall now enters both slots when the heavy channel is Dirichlet). The Neumann leg is bit-identical
under the fix; the Dirichlet leg was re-run. §4.5 is the result from the corrected comparison, and
it is a null for a demonstrated physical reason rather than by construction.

Also corrected: the stiffness figures in the `conduction.cpp` header (it claimed "tens of
milliseconds" and "stiff by three to four orders of magnitude everywhere"; measured `dt` is
5.53–5.62 ms and `tau_eq/dt` is 0.016 / 0.029 / 0.70 at base / mid / top).

## 6. Conclusion and limits

For this model — a 1600–2153 km chromospheric column, physical conduction only, 22 kK
electron-conducted outer reservoir — **the release's single-temperature equilibrium-mixture
assumption is validated to about 0.2 % in the evaporation observables**, with a worst-case
transient decoupling of 0.33 % lasting tens of milliseconds and a quasi-steady decoupling of `1e-5`.

**This does not license the assumption elsewhere.** It measured one configuration. In a flare,
beam-heated, much hotter or much more rarefied regime `tau_eq` rises and the electron heat flux is
far larger, and the conclusion would have to be re-measured.

**Reporting limitation.** The `.twotemp` sidecar prints 10 significant digits, so below about
2080 km the printed `T_e - T_i` sits at the print floor (~`1e-6 K` at 6.6 kK). The interior
decoupling there is an **upper bound**, not a measurement. The `3.28e-3` global maximum is
unaffected, being accumulated in-solver in full double precision over every step.

**Do not quote an interior `Q_ei` from this sidecar.** `Q_ei = g_ei (T_i - T_e)` is proportional to
the very difference that is at the print floor — in the interior the median `|T_e - T_i|` is about
4 print units — so the last printed digit of `T` swings `Q_ei` by tens of per cent. Between the two
BC legs, which agree to `1e-9` in every temperature, the printed interior `Q_ei` differs by up to
35 %. That is quantization, not physics. Fixing it properly means printing a dedicated
`T_e - T_i` column or raising the sidecar precision; the precision was deliberately left at 10 for
this study so that both legs are directly comparable.

**Not re-measured here:** the release CFL acceptance matrix, grid convergence of the decoupling
(only N=500 was run), and any regime other than the one above.
