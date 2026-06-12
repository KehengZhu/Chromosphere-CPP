# Separate Electron Temperature (T_e ≠ T_i) — Implementation Plan

**Status:** ✅ implemented (runtime toggle `ENABLE_TE=1`; default off reproduces the single-temperature baseline bit-for-bit). See **§7 Implementation notes** below for the as-built design, which differs from the original §2 sketch in one structural choice (total-energy + electron-energy split rather than proton-only `E_I`). The ν_ei coefficient is the NRL/Spitzer energy-equilibration rate (`physics.hpp::nu_ei`).

## 1. Motivation

The code is currently a **two-fluid** model (ions + neutrals) but with a single charged-fluid temperature: electrons are lumped with the protons at one temperature. This is encoded by the factor of 2 in `T_i = p_i / (2 n_i k_B)` ([integrators.cpp:129](../src/integrators.cpp#L129), [:173](../src/integrators.cpp#L173)) — `p_i = (n_i + n_e) k_B T = 2 n_i k_B T` under quasi-neutrality (n_e = n_i) — and by the charged heat capacity `C_i = 3 n_i k_B` ([integrators.cpp:177](../src/integrators.cpp#L177)) being that of ions+electrons sharing one temperature T_e = T_i.

Unlike the ion–neutral split (which a 1-D field-aligned model largely cannot exercise, since the dominant two-fluid effect — ambipolar diffusion — is cross-field; see `SINGLE_FLUID` experiment in `outputs/two_fluid_test/`), separating T_e has a **qualitative hook that the chromospheric-evaporation physics actually exercises**: the thick-target electron beam deposits its energy into *ambient electrons* first. With T_e = T_i this is shared instantly with protons and neutrals; with a free T_e, the beam can drive T_e ≫ T_i at the low-density loop-top / impulsive onset before Coulomb collisions equilibrate it. This is exactly the regime where the literature reports decoupling (Manchester 2012: electrons reach only ~10% of the proton temperature in low-density shocks; Russell 2025: T_ion can reach ≥60 MK, several × T_e, in flare onset / above-the-loop-top).

In the dense, cool chromosphere the electron–ion Coulomb equilibration time is short (ν_ei ∝ n_e lnΛ / T_e^{3/2}), so T_e ≈ T_i there — like the neutral coupling, the decoupling lives at the TR / coronal / onset end. But the beam→electron channel makes it worth resolving.

## 2. Equations

Add a **7th conserved variable** `e_e` (electron internal energy). Electrons keep the charged-fluid velocity V (electron inertia negligible) and n_e = n_i by quasi-neutrality — so **no new momentum equation and no new continuity equation**, only an energy equation. The result is a three-temperature model: T_e, T_i (proton-only), T_n.

### Electron internal energy
ε_e = (3/2) n_e k_B T_e = (3/2) p_e, advected at V:

```
∂ε_e/∂t + ∂_s(ε_e V) + p_e ∂_s V =
      ∂_s(κ_e^Spitzer ∂_s T_e)              [electron heat conduction]
    + Q_beam                                 [thick-target beam → electrons]
    − n_e n_H Λ(T_e)                         [optically-thin radiative loss]
    − χ_H (S n_e n_HI − α n_e²)              [ionization cost / recombination release]
    + ν_ei C_e (T_i − T_e)                   [Coulomb e–i equilibration]
    + ν_en C_e (T_n − T_e)                   [e–n collisional equilibration]
```

### Proton internal energy (now proton-only; `p_i = n_i k_B T_i`, factor 1 not 2)
ε_i = (3/2) n_i k_B T_i:

```
∂ε_i/∂t + ∂_s(ε_i V) + p_i ∂_s V =
      ∂_s(κ_i ∂_s T_i)                       [ion conduction ≈ κ_e·√(m_e/m_i) ≈ κ_e/43; optional]
    + ν_ei C_e (T_e − T_i)                   [Coulomb e–i back-reaction]
    + ion–neutral drag / frictional heating  [already in Stage B, integrators.cpp:119]
```

### Charged-fluid momentum (carries the SUM of partial pressures)
```
∂_t(ρ_i V) + ∂_s(ρ_i V² + p_e + p_i) = (p_e + p_i) B ∂_s(1/B) − ρ_i ∂_s φ_g
```
i.e. wherever `flux.cpp` currently uses `p_i` ([flux.cpp:17](../src/flux.cpp#L17), [:22](../src/flux.cpp#L22), [:85](../src/flux.cpp#L85)), it becomes `p_e + p_i`.

### The four physics changes that matter (vs current single-charged-T)
1. **Beam → electrons.** Q_beam goes to ε_e, not split by heat capacity into charged+neutral. (Current: [integrators.cpp:607](../src/integrators.cpp#L607) `apply_beam_heating_stage`.)
2. **κ_e on T_e.** Spitzer κ_e ∝ T_e^{5/2} acts on the electron temperature ([physics.hpp:38](../physics.hpp#L38)). Ion conduction κ_i ≈ κ_e/43 — drop or carry cheaply.
3. **Radiation & ionization → electrons.** n_e n_H Λ(T_e) and the χ_H (13.6 eV) bookkeeping are electron-pool terms (Stage R [integrators.cpp:463](../src/integrators.cpp#L463); Stage E [integrators.cpp:632](../src/integrators.cpp#L632)).
4. **e–i equilibration ν_ei** — the load-bearing parameter. Spitzer (1962) / Braginskii (1965) / NRL formulary rate, ν_ei ∝ n_e lnΛ / T_e^{3/2}. Fast (T_e≈T_i) in the dense chromosphere; slow (decoupling) at the hot tenuous loop-top/onset. **Lock the exact coefficient from Bradshaw 2006 / Manchester 2012.**

The existing Stage C (`apply_temperature_stage`, [integrators.cpp:163](../src/integrators.cpp#L163)) generalizes from a **2-way** T_i↔T_n point-implicit relaxation to a **3-way** (T_e, T_i, T_n) relaxation with pairwise ν_ei / ν_en / ν_in. The current point-implicit solver is the exact template — a 3×3 linear relaxation per cell instead of 2×2.

## 3. Implementation steps

| change | effort | reuses |
|---|---|---|
| `num_of_eq` 6→7; add `cons::E_E` / `prim::P_E`; update `state.cpp` cons↔prim | mechanical | — |
| `flux.cpp`: momentum uses `p_e + p_i`; add electron advective flux `(ε_e + p_e) V` | small | the `p_i` pattern |
| conduction stage: route κ_e onto the new T_e row | small | `ghost_extended_T` per-species solve ([integrators.cpp:229](../src/integrators.cpp#L229)) |
| Stage C: 2-temperature → 3-temperature relaxation (add ν_ei, ν_en) | small | the point-implicit T_i↔T_n solver |
| re-apportion beam / radiation / ionization source terms to ε_e | medium (the physics) | — |
| output 6→7 columns → update `util/animate_loop_flare.py`, `util/diag_two_fluid.py`, `util/compare_two_fluid.py` parsers | mechanical | — |

Estimated ~1 focused day. Highest-risk piece is the source re-apportionment (3), not the plumbing.

## 4. Toggle design — `ENABLE_TE`

Use the **same pattern as `SINGLE_FLUID`** ([chromo_main.cpp:50](../chromo_main.cpp#L50), [integrators.cpp:139](../src/integrators.cpp#L139)/[:185](../src/integrators.cpp#L185)). Always carry the 7th variable; when `ENABLE_TE=0`, force T_e = T_i every step (the instantaneous-Coulomb-equilibration limit — slave electrons to ions exactly as `SINGLE_FLUID` slaves neutrals to ions). This gives a clean runtime on/off **without a variable-size state vector**, and `ENABLE_TE=0` must recover the current single-charged-temperature physics (verify to round-off, as was done for `SINGLE_FLUID`: `n_n` median 0.000% in `outputs/two_fluid_test/`).

Architecturally symmetric — the final code is a three-temperature model whose two decoupling channels are each independently toggleable:
- `SINGLE_FLUID=1` → slave neutrals to ions (drop ion–neutral decoupling)
- `ENABLE_TE=0` → slave electrons to ions (drop electron–ion decoupling)

Heat-capacity bookkeeping for exact reproduction: split `C_e = (3/2) n_e k_B`, `C_i = (3/2) n_i k_B` (sum = `3 n_i k_B` = current charged `C_i`). With T_e = T_i forced, electrons+ions act as the current single 3 n_i k_B pool. Caveat: exact reproduction of the *beam* path requires care — currently the beam is shared charged↔neutral by heat capacity instantaneously, whereas the new scheme routes beam→electrons then equilibrates outward at finite ν_en/ν_in. Expect SINGLE_FLUID-level (negligible) differences, not bit-identical; verify.

## 5. Verification plan

1. Build with `ENABLE_TE=0`, confirm it reproduces a current baseline run to round-off (closed `c05`, open `o00`).
2. Turn `ENABLE_TE=1` and run the same lines; reuse `util/compare_two_fluid.py` (generalized for the T_e column) to map where/when T_e departs from T_i — expect the largest split at the impulsive onset (beam→electrons) and the low-density loop-top, small in the dense chromosphere.
3. Energy-budget check: the beam energy must still close (now routed through electrons → equilibration → ions/neutrals), and total energy conserved across the 3-temperature relaxation stage.

## 6. Where it pays off / why build it

Beam-driven onset and the tenuous loop-top, where electron–ion equilibration is slow — precisely the chromospheric-evaporation onset physics this project studies. Tractable in the 1-D field-aligned geometry (Bradshaw 2006 does exactly this). Strictly a generalization: `ENABLE_TE=0` recovers today's model, so it is pure upside modulo the (modest) numerical machinery.

## References

In `docs/supporting-papers/` (download via the find-paper manifest):
- **Bradshaw2006.pdf** — Bradshaw, Cargill et al. 2006, A&A 458, 987, "Explosive heating of low-density coronal plasma" — 1-D field-aligned, separate e/ion temperatures + evaporation. **Primary methods template.**
- **Manchester2012.pdf** — Manchester et al. 2012, ApJ 756, 81 — Michigan two-temperature (electron + proton) model; e–p equilibration treatment; ~10% coupling in low-density shocks.
- **Russell2025.pdf** — Russell et al. 2025, ApJL, "Solar Flare Ion Temperatures" — physical motivation (T_i ≫ T_e in onset / above-the-loop-top).
- **MacNeice1984.pdf** — MacNeice et al. 1984, Sol. Phys. — historical analog: separate electron vs ion/neutral temperatures *with* time-dependent ionization/recombination in beam-driven evaporation.
- **Liu2009.pdf** (optional) — Liu et al. 2009, ApJ 702, 1553 — Fokker–Planck + hydro (separate e/ion temps), evaporation.

Already in library:
- **Sokolov_2021_ApJ_908_172.pdf** — AWSoM (SWMF), separate electron/proton temperatures.
- Al Shidi et al. 2019 — SWMF two-fluid MHD.

Foundational transport coefficients (ν_ei, κ_e, κ_i): Spitzer 1962; Braginskii 1965; NRL Plasma Formulary.

## 7. Implementation notes (as built)

The model is implemented behind `Grid::enable_Te` (env `ENABLE_TE=1`), carrying the
7th variable always and slaving electrons to ions when off. One structural choice
differs from the §2 sketch and is worth recording:

**Total-energy + electron-energy split (not proton-only `E_I`).** Rather than
redefining `cons::E_I` as proton-only energy with momentum carrying `p_e + p_i`,
`E_I` keeps its original meaning — the **total charged-fluid energy** (protons +
electrons + bulk KE + gravity) — and the new `cons::E_E` carries the electron
internal energy `ε_e = (3/2) p_e`. The proton temperature is *derived*,
`T_i = (p_total − p_e)/(n_i k_B)`, with `p_total = ⅔(E_I − KE − φ)` exactly as
before. This is the standard two-temperature MHD formulation (AWSoM/BATSRUS,
Sokolov 2021) and has two decisive advantages over the proton-only form:

1. **Exact baseline reproduction.** The conservative MUSCL/Rusanov flux, the
   pressure-area + gravity source, the sound speed, and `cal_dt_i` all read
   `p_total` (= the old `p_i`), so they are byte-for-byte unchanged. `E_E` never
   feeds back into `E_I`/momentum/ρ, so `ENABLE_TE=0` reproduces the
   single-temperature run to **0.000e+00** relative error (verified over a full
   flare run, not merely round-off).
2. **Exact total-energy conservation.** Every source stage that moves energy
   between pools updates `E_E` and mirrors the same Δ into `E_I`, so the protons
   (`p_total − p_e`) absorb exactly the complement — no separate proton-energy
   bookkeeping that could drift.

As-built mapping of the four "physics changes that matter" (§2.46):

- **Electron energy equation** — `E_E` advects at `V` (flux `ε_e V`, `flux.cpp`)
  with explicit compression work `−p_e ∇·V` (`cal_source_state`, gated on
  `enable_Te`); the complementary `−p_proton ∇·V` is automatic via `E_I`.
- **Beam → electrons** — `apply_beam_heating_stage` routes all `Q_beam` into the
  electron pool (true path); `cal_dt_i` caps Δt on `ε_e` so the onset spike is
  bounded. Single-T path keeps the heat-capacity charged/neutral split.
- **κ_e on T_e** — `apply_conduction_stage` is parametrized: it conducts
  `T_charged` with `C = 3 n_i k_B` (single-T) or `T_e` with `C_e = 1.5 n_i k_B`
  (three-T, κ_i ≈ κ_e/43 dropped); TRAC keys off the conducted temperature.
- **Radiation & χ_H → electrons** — `apply_radiative_cooling_stage` drains the
  electron pool; `apply_ionization_stage` attributes the χ_H ionization
  cost/return to electrons (thermal/KE inheritance stays with the protons).
- **ν_ei equilibration** — Stage C generalizes the 2-way `T_charged ↔ T_n`
  point-implicit relaxation to a **symmetric 3×3** `(T_e, T_i, T_n)` solve with
  pairwise conductances `g_ei = C_e ν_ei`, `g_en = C_e ν_en`, `g_in = β` (the
  original ion–neutral β, now proton↔neutral). `physics.hpp::nu_ei` is the
  Spitzer/NRL energy rate `≈ 2.03×10⁻⁴³ n_e lnΛ /(k_B T_e)^{3/2}` (fast in the
  dense chromosphere, slow at the tenuous loop-top → decoupling).

**Verification.** Tests stay green (5926). `ENABLE_TE=0` ≡ baseline (Δ = 0).
`ENABLE_TE=1` on the 2024-08-01 PFSS event lines is stable and shows T_e running
up to ~6× T_i in the beam-heated footpoints/loop-top during the impulsive onset,
relaxing toward T_e = T_i as the loop fills and ν_ei speeds up — the predicted
beam-onset / low-density decoupling. Movies in `util/visualization/event_20240801_*_Te.mp4`.
