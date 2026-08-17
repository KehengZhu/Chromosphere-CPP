# Explosive chromospheric evaporation — scenario plan

## Goal

Test whether a physically-plausible flare energy-deposition driver, added to the existing C7 chromosphere column, can reproduce **explosive chromospheric evaporation** in the sense of Fisher, Canfield & McClymont (1985, ApJ 289, 425): a short, intense burst of energy deposited in the upper chromosphere that the gas cannot radiate away, so it is heated to coronal temperatures and blasts upward into the corona at ~100s of km/s, with a matching downward condensation recoil.

## The physics we must reproduce (Fisher et al. 1985)

1. **Energy-flux threshold (their §II).** Explosive evaporation occurs when the heating rate exceeds what the upper chromosphere can radiate near the 10⁵ K radiative-loss peak. Their threshold for a VAL/F preflare atmosphere is `F_crit ≈ 7×10⁹ erg cm⁻² s⁻¹ = 7×10⁶ W m⁻²`. Below it → *gentle* evaporation (slow, conduction-driven, weak upflow); above it → *explosive*. **This means the radiative sink must be active** — the threshold is defined by the competition heating-vs-radiation.
2. **Timescale condition (their eq. 1):** `3kT/Q ≲ L₀/c_s` — the heating must be faster than the hydrodynamic expansion time. Deposit over ~10 s.
3. **Velocity ceiling (their eq. 11–12):** upflow `v ≤ 2.35 c_s` of the evaporated (coronal-T) material; their `F = 5×10¹⁰ erg cm⁻² s⁻¹` run peaks at ~500 km s⁻¹, T ~ 10⁶–10⁷ K.
4. **Driver:** a thick-target nonthermal electron beam, heating rate per particle `Q(N) ≈ F_e/N*`, deposited where the beam stops in the upper chromosphere.

## Design decisions (made autonomously)

- **Driver = volumetric beam-heating source term, not a hot-wall BC.** Explosive evaporation is fundamentally driven by *energy deposited inside* the upper chromosphere (a beam), not by a boundary temperature. A Dirichlet hot wall over-conducts on this coarse grid (we learned this building the C7 upper BC). So the "boundary condition" the user asked for is realized as a **transient interior heating term** representing the precipitating flare beam — the standard RADYN/Fisher driver — plus an **opened outer boundary** so the evaporated plasma can actually leave the top (the existing C7 Mach-0.05 cap would choke a 100 km/s upflow).
- **Deposition profile.** Thick-target stopping (~10²⁰ cm⁻² for tens of keV) is *sub-grid* on the 11 km C7 mesh (a single top cell already holds ~10²⁶ m⁻²), so we deposit the flux over a finite upper-chromosphere layer of thickness `Δh_beam` (default 80 km below the top) weighted ∝ local density `n_tot` (densest reachable layer absorbs most — the signature of beam stopping), normalized so the column integral equals `F_e`. Honest about the coarse-grid limitation; documented as a known idealization.
- **Heating goes into the electron/ion thermal pool** (`p_i = 2 n_i k T_e`), as a beam heats ambient electrons. Added as a separate operator-split stage `apply_beam_heating_stage` (explicit additive source, exact since it is T-independent) placed **immediately before radiative cooling** so the heating-vs-radiation threshold competition is captured.
- **Everything else physical stays on:** ionization (Stage E), optically-thick + optically-thin radiative cooling (Stage R, defines the threshold), TRAC, the C7 inner discrete-HSE base, the ambient RTV q(T) (negligible vs the beam).
- **Two runs:** an **explosive** case `F_e = 5×10¹⁰ erg cm⁻² s⁻¹ = 5×10⁷ W m⁻²` (well above threshold, Fisher's middle explosive run) and a **gentle control** `F_e = 1×10⁹ erg cm⁻² s⁻¹ = 1×10⁶ W m⁻²` (below threshold). Same everything else → isolates the threshold.
- **Timing:** beam on at `t_on = 2 s` (let the IC settle one acoustic time), duration `τ = 10 s`, 1 s cosine ramps. Run ~80–120 s total to capture onset → upflow → quasi-steady.

## Implementation

1. **`chromosphere.hpp` Grid fields:** `enable_beam_heating`, `beam_flux` [W/m²], `beam_t_on`, `beam_duration`, `beam_ramp` [s], `beam_deposition_height` [km], `outer_free_outflow` [bool], `sim_time` [s]. Declare `apply_beam_heating_stage`.
2. **`physics.hpp`:** `beam_heating_rate(grid, n_i, n_n)` → Vec [W/m³], column-depth-from-top deposition over `Δh_beam`, ∝ n_tot, normalized to `beam_flux`, gated by the temporal window `g(sim_time)`.
3. **`src/integrators.cpp`:** `apply_beam_heating_stage` (Δp_i = (2/3)·dt·Q_beam); call it in `advance_Euler_state` before cooling. TRAC-broaden consistently (divide by ε like the thin loss).
4. **`scenarios/model_flare.cpp/.hpp`:** `ic` = `model_c7_ic` + set beam params + `outer_free_outflow=true` + force cooling/ionization on; `update_bc` = `model_c7_update_bc` (which honors `outer_free_outflow` → Neumann V instead of Mach cap).
5. **`model_c7_update_bc`:** honor `outer_free_outflow` (skip the Mach cap → free supersonic outflow).
6. **`chromo_main.cpp`:** set `grid.sim_time = time` each step.
7. **Register** `model_flare` in `make_scenario`; add to `CMakeLists.txt`; add tests (beam profile normalization + gentle/explosive monotonicity).
8. **Run** both cases; **visualize** the height–time evolution of T, V, n (the classic evaporation diagram) + an 8-panel evolution movie, into `visualization/`.

## Plausibility checks (post-run)

- Explosive upflow reaches ~10²–10³ km s⁻¹ but `≤ 2.35 c_s(T_peak)`.
- Top heated to ~10⁶–10⁷ K; gentle control stays cool with only a weak (≲30 km/s) upflow.
- Beam fluxes are in the observed flare range (10⁹–10¹¹ erg cm⁻² s⁻¹).
- A condensation downflow appears below the evaporation front.
