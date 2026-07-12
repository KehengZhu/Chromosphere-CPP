# Gentle (conduction-driven) chromospheric evaporation — plan

**Status:** IMPLEMENTED (all four phases). **Date:** 2026-06-22. **Purpose:** "add heating, then relax" workflow for a steady chromosphere→TR→corona initial condition, then drive gentle evaporation. **Owner:** Keheng Zhu (advisor-directed).

---

## v2 (CURRENT DIRECTION, 2026-06-22): corona as a BOUNDARY, not resolved

**Decision (advisor-directed):** abandon the resolved-corona loop. Stay in the single C7 column, **lift the outer boundary from h ≈ 2.153 Mm → 2.6 Mm**, and represent the corona entirely through the existing `model_c7` upper boundary condition — the RTV downward conductive flux `q(T)` plus a reservoir `(ρ, T_top)`. **Drive gentle evaporation by ramping `q(T)` up.** This is conduction-driven by construction (Antiochos & Sturrock 1978: a coronal conductive flux drives the upflow).

**Why (over v1's resolved corona):**
- **No TNE drainage.** There is no long resolved corona to drain — the Klimchuk & Luna (2018) apex/footpoint criterion is moot. The drainage that plagued v1 cannot occur.
- **Reuses tested machinery.** `impose_outer_heat_flux`, the fixed TR-base reservoir, TRAC, the grid-scaled numerical diffusivity, and the vacuum floor are all already validated for `model_c7`.
- **The corona is needed only as a BC** (see below) — `q(T)` both (1) holds the steady TR/upper-chromosphere up against radiative+ionization cooling (drop it → condensation runaway, §2) and (2) is the evaporation driver when ramped.
- The volumetric `H(s)` term and the `GENTLE` overlay from v1 stay in the codebase (off); they are not the driver here.

**Settled with the user:** extend only to **2.6 Mm**; driver = **ramp the imposed `q(T)` flux** (not raise T_top).

**Key correction (mined from the paper):** Avrett & Loeser (2008) **Table 26 is the COMPLETE atmosphere**, h = −100 km → **68,084 km (68 Mm)**, T = 9.4 kK → **1.586 MK** — the full chromosphere + TR + corona. The code's `MODEL_C7` array merely truncated it at row 68 (h = 2152.9 km, 23.1 kK). So the 2.153 → 2.6 Mm extension uses **REAL C7 data, not a prescribed ramp**. Reference rows (h[km], T[K], n_e[cm⁻³], n_HI[cm⁻³]): row 68 = 2152.9, 2.31e4, 1.24e10, 2.77e8 (current top); row 47 = 2175.3, 7.63e4, 4.05e9, 1.61e5; row 32 = 2269.7, 1.59e5, 2.05e9, 9.48e3; **row 20 = 2628.3, 2.93e5, 1.14e9, 1.55e3 (≈ 2.6 Mm)**; row 5 = 21133, 1.08e6 (1 MK at 21 Mm); row 1 = 68084, 1.59e6 (model top). **At 2.6 Mm T ≈ 0.29 MK — upper TR, NOT corona; the real corona (≥1 MK) is at ≥21 Mm.**

### Phase A — extend the C7 domain to 2.6 Mm (real Table 26 data)
- **Append Avrett & Loeser Table 26 rows 67 → ~20 (h = 2153 → 2628 km)** to the `MODEL_C7` array in `model_c7_ic` — measured T, n_e, n_HI; no prescribing. New top reservoir = the real row-20 state (2.6 Mm, 0.29 MK).
- Boundary flux `q₀`: either the existing RTV estimate, or — better, now that we have C7's resolved TR gradient — **computed directly from C7's own `q = κ₀ T^{5/2} dT/ds` at 2.6 Mm (≈ 100 W/m² downward)**, self-consistent with the model.
- `enable_coronal_heating` stays OFF — pure boundary-driven; no other physics changes.
- (Aside: since Table 26 has the full corona to 68 Mm, a *resolved*-corona variant from real C7 data is also possible later — but the agreed direction is corona-as-boundary at 2.6 Mm.)

### Phase B — relax to the steady preflare IC
- Run extended-`model_c7` with steady `q₀` + conduction + cooling, beam off; relax to V≈0 (`check_relaxation.py`, now reporting the upper-domain max|V|, not just the mass-weighted mean). Characterize the residual TR flow (the `q(T)`↔cooling pointwise-imbalance drift) — it is the baseline; evaporation is the departure from it.

### Phase C — drive gentle evaporation by ramping `q(T)`
- From the relaxed IC, ramp `q(T) = q₀·enhance(t)` up (×2–5, slowly, ≫ the TR sound-crossing time); keep sub-threshold so it stays gentle (v ≪ c_s). Free-outflow top BC so the upflow exits at 2.6 Mm.
- Diagnostics: height profiles of V, T, n (before vs during); TR upflow speed + Mach; **mass flux through the top** (what would fill a real corona); Doppler blueshift via `util/synth_doppler.py`.

### Code changes
1. **Extend the C7 IC to 2.6 Mm** in `model_c7_ic` (prescribed TR ramp + reservoir above 2.153 Mm; recompute the top ghost).
2. **Time-ramp on `q(T)`** — Grid fields `outer_heat_flux_t_on/_ramp/_enhance`, applied as `q(t)=q₀·enhance(t)` in `model_c7_update_bc` (mirrors the coronal-heat ramp already in the tree).
3. **Driver overlay + env knobs** (`GENTLE_Q_ENHANCE/_T_ON/_RAMP`, free-outflow) to relax-then-ramp in one run.

### Honest tradeoff
The corona is outside the box, so we **do not see the coronal EM rise** inside the domain (v1's headline signature). We instead measure the TR upflow speed/Mach, the TR burndown, the upper-domain density rise, and the **mass flux out the top**. The upflow + Doppler are the primary observables anyway.

### IMPLEMENTATION + RESULTS (2026-06-22)
**Built & tested (5950 pass).** Phase A: `model_c7_ic(extended=true)` appends the real Table-26 rows to h ≈ 2.628 Mm and computes the C7-gradient `q₀ ≈ 130 W/m²`. The `q(T)` ramp (`outer_heat_flux_base/_t_on/_ramp/_enhance`, applied in `model_c7_update_bc`) + the new `model_gentle` scenario (`GENTLE_Q_ENHANCE/_T_ON/_RAMP`, env-driven) are the driver. `outer_mach_cap` Grid field added (default 0.05; `model_c7` tests pin 0.05, unaffected). Commands in `visualize_commands.md` §10.

**Phase B (relaxation) CONVERGED:** with the Mach-capped outflow (cap 0.05), the column relaxes to a quiet steady preflare state — ρ-weighted mean|V| = **48 m/s**, max|dlnT/dt| = **3×10⁻⁵/s**, top settled at a stable **0.45 MK** (`gentle_v2_relax_convergence.png`).

**KEY BC FINDING:** a free/Neumann outflow at the top is **ill-posed once the residual is a downflow** — it feeds mass in and the upper TR catastrophically collapses (dense 60 kK slug raining at 190 km/s). The Mach-capped outflow is required. A looser cap (0.5) lets more throughflow → top runs hotter (0.76 MK) and doesn't settle; the tight 0.05 cap is cleanest.

**Phase C (q(T) ramp ×3) — stable, conduction-driven response:** the TR **burned down** (the 0.04 MK isotherm at 1.58 Mm heated to 0.27 MK; the TR shifted ~0.1 Mm lower), the upper chromosphere/TR **densified** (n_e ×1.4–1.9), and the top stepped **0.45 → 0.62 MK**. Movie: `gentle_v2_evolution.mp4`.

**Geometry caveat confirmed:** in corona-as-boundary AT 2.628 Mm, the q(T)-ramp response is **conductive TR-burndown + in-place densification (the evaporate accumulates)**, NOT a fast Doppler upflow — and the column-integrated density stays **flat** (no evaporation). Diagnosis: the boundary at 2.628 Mm sits right at the top of the TR with essentially **no coronal volume above it inside the domain**, so the evaporated mass has nowhere to accumulate. Evaporation *is* the coronal EM rise (Antiochos & Sturrock 1978) — pin the corona at the boundary and the phenomenon is removed by construction. The driver (downward q(T)) was correct; the geometry forbade the response. → fixed in option A below.

### Option A — lift the boundary HIGH into a resolved corona (2026-06-22, WORKS)
The minimal fix that keeps the corona-as-boundary philosophy: don't put the boundary at the TR — put it high in a real corona, with a genuine coronal **volume** between the TR and the top for the evaporated mass to fill.

**Implemented:** extended the C7 table (`MODEL_C7`, N_C7 51→70) with the real Table-26 **coronal** rows up to h ≈ 68 Mm / 1.59 MK. New Grid field `corona_top_km`; `model_c7_ic(extended)` truncates the table at it. `model_gentle` sets it from `GENTLE_TOP_KM` (default **10 Mm**, T ≈ 1.1 MK relaxed). `q(T)` at the top is still C7's own conductive gradient there (≈230 W/m² base). Run at `GENTLE_NS ≈ 300` (TR via TRAC + corona both resolved). 5950 tests still pass (model_c7 `extended=false` unchanged).

**Phase B (relaxation) — the corona HOLDS, does not drain.** Over ~3000 s the resolved C7 corona under the q(T) top BC settles to a stable T_top ≈ 1.14 MK, flat coronal EM, TR at 1.595 Mm. (It carries a steady ~8 km/s cap-throttled coronal downflow — a mild draining tendency the boundary throttles to a stable EM, since there's no volumetric heating to pin V=0. Steady, not drifting → usable preflare baseline.)

**Phase C (q(T) ramp ×3 at t=3000 s) — EVAPORATION APPEARS.** Coronal (h>2.3 Mm) emission measure ∫n_e² ds **×2.0**, coronal mass ∫n_e ds ×1.42, T_top 1.14→1.50 MK, n_e(top) ×1.47, TR burns down 1.595→1.526 Mm, transient +1 km/s upflow during the ramp. The corona gains mass — the A&S 1978 signature, **absent in v2**. Comparable to v1 (EM ×2.6). Movie `gentle_v3_evolution.mp4`; commands `visualize_commands.md` §9c.

**Remaining limitation:** the steady baseline downflow means the *sustained* signature is the EM/density rise, not a fast bulk upflow; v1 (resolved corona + volumetric heating) gives the cleaner sustained upflow (15.6 km/s). Pinning V≈0 in the corona-as-boundary baseline would need a small volumetric coronal heating term (→ v1's H(s)), which is the deliberate tradeoff this branch avoids.

---

## v1 (SUPERSEDED — resolved-corona + volumetric H(s); code retained, off)

All four phases are built and verified; commands catalogued in `visualize_commands.md` §9.

- **Phase 0 — heating term.** New `Grid` fields `enable_coronal_heating`, `coronal_heat_E0/_sH/_s0/_two_sided` and ramp fields `coronal_heat_t_on/_ramp/_enhance` ([chromosphere.hpp](../chromosphere.hpp)). `coronal_heating_rate(grid)` → W/m³ in [physics.hpp](../physics.hpp): `H(s)=E_H0·exp(−d(s)/s_H)·enhance(t)`, with `d(s)` the distance to the nearest footpoint (one-sided for open/half-loop, two-sided for a full loop). `apply_coronal_heating_stage()` ([src/integrators.cpp](../src/integrators.cpp)) mirrors the beam (electron pool when `enable_Te`, else heat-capacity-shared; TRAC-broadened; runs before cooling) plus a CFL guard. Default off ⇒ baseline byte-for-byte; **5950 tests pass** (added `test_coronal_heating_rate_profile`, `test_coronal_heating_partitions_by_heat_capacity`).
- **Phase 1 — calibration.** `GENTLE=1` env overlay in [pfss_field_line.cpp](../scenarios/pfss_field_line.cpp) (alternative to `PFSS_FLARE`) inverts the RTV (1978) + Aschwanden & Schrijver (2002) scaling: `p₀=(T_max/1400)³/L_cm`, `E_unif=9.8e4·p₀^{7/6}·L_cm^{−5/6}`, `E_H0=E_unif·exp(½L/s_H)`, ×0.1→SI. Knobs `GENTLE_TMAX_MK` (2.0), `GENTLE_SH_FRAC` (0.4 — keep ≳ 0.3, the Serio/Martens static-existence limit), `GENTLE_E0/_S0_MM/_ENHANCE/_T_ON/_RAMP`. `data_file_parser` now reads `loop_half_length_m`. For L=25 Mm, T_max=2 MK ⇒ **E_H0≈6.0×10⁻⁴ W/m³, s_H=10 Mm**. `make_loop_dat.py` prints this auto-calibration.
- **Phase 2 — relaxation.** New `util/check_relaxation.py` tracks `max|V|(t)`, ρ-weighted `mean|V|(t)`, `max|∂lnT/∂t|`, apex/coronal T. The closed L=25 Mm loop **converges**: ρ-weighted mean |V| → 112 m/s, `max|dlnT/dt|` → 1.7×10⁻³ /s, apex → 2.26 MK (`gentle_relax_convergence.png`). Residual ~10 km/s is a bounded TR cell, not a growing drift.
- **Phase 3 — gentle evaporation.** Ramp ×3 over 200 s from the relaxed state (one combined run, ramp gated to `t_on=1000 s`). Result is textbook gentle/conduction-driven: coronal **EM ×2.6**, apex T **2.4→3.1 MK**, peak upflow **15.6 km/s at Mach 0.07** (subsonic ⇒ gentle, A&S 1978), TR burns down, corona fills (`gentle_evap_signature.png`, `gentle_evap_profiles.png`; diagnostics `util/gentle_evaporation_diag.py`, `util/gentle_evap_profiles.py`).

**Settled decisions:** (1) s_H = 0.4 L (not 0.2 L — static-existence margin). (2) E_H0 from the RTV/AS inversion above; T_max is the relaxed outcome (2.26 MK vs 2.0 target — acceptable). (3) Ramp ×3 over 200 s stays gentle (M≈0.07 ≪ 1).

---

**Original plan (pre-implementation handover) follows.**

**Purpose:** handover document so the next session can implement the "add heating, then relax" workflow for a steady chromosphere→TR→corona initial condition and then drive gentle evaporation. **Owner:** Keheng Zhu (advisor-directed).

---

## Background

The project evolves a 1.5D field-aligned two-fluid (plus optional separate electron temperature) model of the solar atmosphere. The `model_c7` scenario ([scenarios/model_c7.cpp](../scenarios/model_c7.cpp)) initializes from the semi-empirical Model C7 atmosphere (Avrett & Loeser 2008, Table 26), covering the chromosphere up into the lower transition region (h ≈ 1003 → 2153 km, T ≈ 6225 → 23100 K). The `pfss_field_line` scenario ([scenarios/pfss_field_line.cpp](../scenarios/pfss_field_line.cpp)) runs the same physics on a PFSS-traced loop geometry and, via `util/make_loop_dat.py`, can extend the domain through a transition region into a resolved corona out to the loop apex.

**Science goal (advisor's directive):** study **gentle, conduction-driven chromospheric evaporation** — the generic mechanism where a hot corona drives a downward conductive heat flux into the transition region and upper chromosphere, heating that material until it expands upward into the corona. This is explicitly **not** the electron-beam-induced *explosive* evaporation already implemented in `model_flare` (Fisher, Canfield & McClymont 1985; see [docs/explosive_evaporation_plan.md](explosive_evaporation_plan.md)). The relevant physics is the gentle/conductive regime of Antiochos & Sturrock (1978).

**Conventions** (for whoever reads the output): arc length / height `s` increases upward, so **negative velocity = downflow**. Output columns are `[RHO_I, RHO_N, MOM_I, MOM_N, E_I, E_N]` (+`E_E` when `ENABLE_TE=1`); `v` = ion velocity, `u` = neutral velocity (see `util/plot_output.py::primitives`).

---

## Problem

### 1. A flare/evaporation study needs a steady-state initial condition

Evaporation observables (upflow speeds, Doppler blueshifts, coronal emission-measure rise) are measured as **departures from the preflare state**. If the IC is not a steady state, the atmosphere is already moving on its own, and that baseline drift cannot be separated from the driven evaporation signal. This is standard flare-loop radiative-hydrodynamics methodology: Fisher, Canfield & McClymont (1985) start from a preflare loop atmosphere in hydrostatic equilibrium, and RADYN models (Allred et al. 2005) relax to a radiative-equilibrium preflare atmosphere before applying the driver. **The relaxation must use the same physics as the production run** — relaxing under a different physics set yields an IC that immediately drifts when the production physics is switched on.

### 2. The current `model_c7` run does not relax to V = 0

The `model_c7_evolution` movie (`outputs/model_c7_movie.txt`, t → 588 s) reaches only a *driven quasi-steady* state, not a static one. The chromosphere proper settles to a small smooth flow (density-weighted mean ≈ −33 m/s), while the TR top holds a persistent ≈ −1.37 km/s downflow that sits right at the boundary velocity cap. Three drivers, none of which a different IC would remove:

- **C7 is ~20–26% off discrete hydrostatic equilibrium under the code's EOS** (it was built from NLTE radiative balance, not ideal-gas HSE; see comment at [scenarios/model_c7.cpp](../scenarios/model_c7.cpp) `model_c7_ic`). Only the inner ghost is engineered so V=0 is an exact fixed point; the interior carries C7's tabulated (n,T), which is not in discrete HSE.
- **q(T) ↔ radiative-cooling imbalance.** The imposed coronal conductive flux `q ≈ 29 W/m²` at the top and the volumetric radiative cooling do not balance *pointwise*; the residual local imbalance is continuously converted into bulk motion (advection closes the energy budget that static conduction+radiation cannot). This is genuine thermal non-equilibrium.
- **Open Mach-0.05 outflow boundary.** The outer face is transparent (zero-gradient ρ, T, V) with the velocity magnitude clamped to `0.05 c_s ≈ 1.26 km/s`; it never imposes V=0 and sustains a throughflow rather than damping it.

So "the IC doesn't matter because it reaches steady state" is only partly true: the steady state is set by the BCs + source terms, and the residual V ≠ 0 comes from the driving + open BC, not from the IC choice. Disabling q(T) is **not** the fix — q(T) is real coronal conductive heating (the leading-order heat input wherever the domain reaches the TR), and removing it either leaves the corona/TR thermally unsupported (q(T) off + cooling off → artificial collapse of the TR temperature structure) or causes catastrophic cooling (cooling on + q(T) off → condensation runaway). It also wouldn't address the open BC.

### 3. A resolved corona still needs an IC — and the code can't currently make a steady one

The advisor's gentle-evaporation mechanism (hot corona → conductive flux → evaporation) needs the corona *inside* the domain (so evaporated plasma has somewhere to go, and so coronal heating is where it physically belongs, rather than imposed as a boundary `q(T)` crutch). `util/make_loop_dat.py` already builds such a domain: C7 chromosphere + cosine-ramp TR + **hydrostatic** corona at a prescribed apex temperature (default `--apex-T-MK 2.0`), with `util/check_loop_ic.py` checking the HSE residual.

**Key finding:** that construction is *mechanical* equilibrium only — the coronal temperature is **prescribed**, not the result of an energy balance. The code has **no ambient coronal volumetric heating term** — only the flare beam (`apply_beam_heating_stage` in [src/integrators.cpp](../src/integrators.cpp)) and the boundary `q(T)` (a truncated-domain trick). A corona is not a static equilibrium on its own: it exists only because of continuous heating. With the corona resolved and no heating term, it radiates and conducts its energy away and drains/cools. **Therefore a steady chromosphere→TR→corona IC is impossible without adding a heating term.** The static loop equilibrium is the Rosner–Tucker–Vaiana (1978) balance:

- hydrostatic:  `dp/ds = −ρ g_∥(s)`
- energy:    `d/ds(κ T^{5/2} dT/ds) + H(s) − n² Λ(T) = 0`  (conduction + heating − radiation)

`make_loop_dat.py` gives the first line; the second needs the heating term `H(s)`, which is missing.

### 4. The evaporation paradox

"If the IC is relaxed with heating on, won't re-running show no evaporation because it's already steady?" — Yes, for an *unchanged* heating rate; that is the definition of steady state. **Gentle evaporation is the response to a *change* (an increase) in the coronal heating, not to its presence.** Relax with the quiet/preflare heating to get the V≈0 baseline, then raise the heating: the enhanced downward conductive flux exceeds what the upper chromosphere can radiate, heating it toward coronal temperatures and driving the upflow. The amount of evaporation is set by the *increment*, so a steady IC is the required zero baseline, not an obstacle.

---

## Solution — "add heating, then relax," then ramp the heating

Four phases. Phase 0's structure does not depend on the calibration numbers (Phases 1–3), so it can be implemented before the reference papers are read.

### Phase 0 — Add an ambient coronal heating term `H(s)`
- **Parametrization** (Aschwanden & Schrijver 2002): `H(s) = E_H0 · exp(−(s − s₀)/s_H)`, footpoint-anchored, with a uniform-heating limit (`s_H ≫ L`) as a control. Martens (2010) shows `T(s)` is only weakly sensitive to the heating shape, so a simple exponential is defensible; Serio et al. (1981) give the gravity-stratified scaling tying `E_H0, s_H, L → T_max`.
- **Code:**
  - [chromosphere.hpp](../chromosphere.hpp): new `Grid` flags/fields — `enable_coronal_heating`, `coronal_heat_E0` (W/m³), `coronal_heat_sH` (m), `coronal_heat_s0` (m), plus time-ramp fields for Phase 3 (`coronal_heat_t_on`, `coronal_heat_ramp`, `coronal_heat_enhance`).
  - [physics.hpp](../physics.hpp): `coronal_heating_rate(grid, s)` → W/m³.
  - [src/integrators.cpp](../src/integrators.cpp): new `apply_coronal_heating_stage()` modeled on `apply_beam_heating_stage` (~L754) — deposit into `p_i` (and the electron pool when `enable_Te`); reuse the beam CFL guard pattern (~L27–51). Call it in `advance_Euler_state`, always-on for steady runs (unlike the gated beam).
- **Validation:** `enable_coronal_heating = false` must reproduce the existing baseline byte-for-byte (gate like the beam/floor stages), so the test suite is unchanged.

### Phase 1 — Calibrate `E_H0`, `s_H` to a target corona
Pick a target apex temperature (quiet/AR `T_max ≈ 1.5–2 MK`) and loop half-length `L` from the PFSS trace, then invert the Aschwanden & Schrijver (2002) scaling `E_H0(L, s_H, T_max)` to set the magnitude; cross-check the absolute value against Klimchuk (2006). Start with `s_H ≈ 0.2 L` (footpoint-concentrated; Aschwanden et al. 2001 measured `s_H ≈ 12 ± 5 Mm`, `s_H/L ≈ 0.2 ± 0.1`) and a uniform-heating control. Literature magnitude anchor: `E_0 ≈ 10⁻² erg cm⁻³ s⁻¹ ≈ 10⁻³ W m⁻³` (Ishigami et al. 2024).

### Phase 2 — Relax to the steady preflare IC
- **Initial guess:** `util/make_loop_dat.py` (C7 chromosphere + cosine TR + hydrostatic corona at `apex_T = target`).
- **Run:** `pfss_field_line` scenario, **beam off**, with `enable_coronal_heating` + conduction + radiative cooling on; integrate until steady.
- **Convergence criterion:** extend `util/check_loop_ic.py` to track `max|V|(t)` and `max|∂T/∂t|`; declare the fixed point when `max|V|` is far below the expected evaporation speed (target ≲ a few hundred m/s) and the profiles are stationary. Verify global energy balance `∫ H ds ≈ ∫ n²Λ ds + conductive losses`. Save the relaxed dump as the preflare IC.

### Phase 3 — Drive gentle, conduction-driven evaporation
- From the relaxed dump, **ramp `H(s)` up** (raise `E_H0`, or add a localized coronal heating enhancement) over a few seconds. Keep it **sub-threshold** so it stays gentle: per Antiochos & Sturrock (1978), gentle evaporation is the regime where conduction dominates radiation and `v ≪ c_s`; a large/fast jump tips into the Fisher (1985) explosive regime — avoid.
- Use an outflow-permitting / resolved-corona configuration so the modest (~tens of km/s) upflow is not throttled by the Mach-0.05 cap.
- **Diagnostics:** arc-length profiles with the split x-axis (see `CLAUDE.md` and `util/animate_loop_flare.py`), `V`, `T`, `n`, and coronal emission-measure rise (the classic gentle-evaporation signature); optionally feed `util/synth_doppler.py` for the blueshift.

### Decisions to settle once the papers are read
1. **`s_H`** — footpoint (`~0.2 L`) vs. uniform (Martens 2010: `T(s)` barely cares, but density/stability do).
2. **`E_H0`** absolute magnitude and the `T_max` it yields (Serio 1981 / Aschwanden & Schrijver 2002 scaling vs. Klimchuk 2006 observed range).
3. **Phase-3 ramp amplitude/duration** to stay in the Antiochos & Sturrock (1978) gentle regime.

---

## Code map (key locations)

- [scenarios/model_c7.cpp](../scenarios/model_c7.cpp) — C7 IC, boundary closures, `q(T)` setup (`impose_outer_heat_flux = enable_radiative_cooling`), photoionization inversion, vacuum floor.
- [scenarios/pfss_field_line.cpp](../scenarios/pfss_field_line.cpp) — loads the loop `.dat` (`pfss_ic`), outer-face reflecting/open logic, optional `PFSS_FLARE` beam overlay.
- [scenarios/scenario.cpp](../scenarios/scenario.cpp) — `apply_open_bcs` (zero-gradient outer ghost; reflecting variant).
- [src/integrators.cpp](../src/integrators.cpp) — stage sequence in `advance_Euler_state`; `apply_beam_heating_stage` (template for the new heating stage); `apply_conduction_stage` (Stage D, where `q(T)` enters); `apply_radiative_cooling_stage` (Stage R).
- [src/rhs.cpp](../src/rhs.cpp) — explicit MUSCL/Rusanov flux + gravity/pressure source; implicit drag/conduction RHS.
- `util/make_loop_dat.py` — builds the chromosphere→TR→corona loop `.dat` (hydrostatic, prescribed-T). **Add an option to set the relaxation target / heating-consistent apex_T here.**
- `util/check_loop_ic.py` — HSE-residual check. **Extend with `max|V|(t)` and `dT/dt` convergence monitors.**
- `util/animate_loop_flare.py` — split-axis loop visualization (auto-detects half/full topology).

---

## References

**Papers** — target folder `docs/supporting-papers/`. Already present: `RosnerTuckerVaiana1978.pdf`, `FisherCanfieldMcClymont1985.pdf`. **Still to download** (Consensus links below; save with the listed filenames): Serio 1981, Aschwanden & Schrijver 2002, Martens 2010, Klimchuk 2006, Antiochos & Sturrock 1978.

- Rosner, Tucker & Vaiana 1978, ApJ 220, 643 — `RosnerTuckerVaiana1978.pdf`. Foundational static-loop scaling laws `T_max ≈ 1.4×10³ (pL)^{1/3}`; loops are stable only with `T_max` at the apex and heating scale length ≳ loop size.
- Serio, Peres, Vaiana, Golub & Rosner 1981, ApJ 243, 288 — `Serio1981.pdf`. Gravity-stratified generalization of RTV; scaling laws relating base pressure/heating to `s_H`, `s_p`, `L`. [Consensus](https://consensus.app/papers/details/5c88f75c0dde584ab19dad20f856cd24/)
- Aschwanden & Schrijver 2002, ApJS 142, 269 — `AschwandenSchrijver2002.pdf`. The `E_H0 · exp(−s/s_H)` heating parametrization; analytic `T(s)`, `n(s)`, `p(s)`; scaling laws `p₀(L, s_H, T_max)` and `E_H0(L, s_H, T_max)`. Primary methods-and-parameters reference. [Consensus](https://consensus.app/papers/details/5a8a27bfb83757b3a374f108bbe7ac99/)
- Martens 2010, ApJ 714, 1290 — `Martens2010.pdf`. Analytic temperature profiles/scaling for non-uniform heating; `T(s)` weakly depends on heating distribution; no quasi-static solution for excessive footpoint concentration. [Consensus](https://consensus.app/papers/details/d52618683eaa54338949d5fcf54fcd71/)
- Klimchuk 2006, Sol. Phys. 234, 41 — `Klimchuk2006.pdf`. Coronal-heating review; heating magnitudes/parametrization, impulsive-vs-steady perspective. [Consensus](https://consensus.app/papers/details/5ffca13c5d7450b797b4c5fb35bcd6fb/)
- Antiochos & Sturrock 1978, ApJ 220, 1137 — `AntiochosSturrock1978.pdf`. Gentle/conduction-driven evaporation: conductive losses dominate radiation, evaporative velocities `≪ c_s`; the target regime for this study. [Consensus](https://consensus.app/papers/details/d66d7b383e225b16bad2f464696b1943/)
- Fisher, Canfield & McClymont 1985, ApJ 289, 425 — `FisherCanfieldMcClymont1985.pdf` (already present). Explosive-vs-gentle threshold; preflare hydrostatic atmosphere methodology. Context / contrast (this is the beam-driven path we are *not* using).
- Allred, Hawley, Abbett & Carlsson 2005, ApJ 630, 573 — RADYN preflare-relaxation methodology. (Secondary; fetch if needed.)

**Parameter anchors from the literature search:** Aschwanden et al. (2001, ApJ 550, 1036) — `s_H ≈ 12 ± 5 Mm`, `s_H/L ≈ 0.2 ± 0.1` from TRACE loops. Ishigami et al. (2024) — `E_0 ≈ 10⁻² erg cm⁻³ s⁻¹`, `F_H ≈ 10⁷ erg cm⁻² s⁻¹`, `s_H ≈ 10 Mm`.

**Related project docs / memory:** [docs/explosive_evaporation_plan.md](explosive_evaporation_plan.md) (the beam path), [docs/boundary_conditions_plan.md](boundary_conditions_plan.md) (the `q(T)` upper-BC closure), and memory note `gentle-evaporation-add-heating-relax`.
