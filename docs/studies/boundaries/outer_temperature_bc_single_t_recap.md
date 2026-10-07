# Outer temperature BC: one fixed ghost-cell temperature (release as of 2026-09-22)

**Status: CURRENT**

**Decision.** The `model_column` release now imposes a single fixed outer temperature, `T_c = ISO_T_TOP = 22 000 K`, in the upper hydro ghost cells, and conduction reads the same ghost cell as its Dirichlet datum (`ISO_HYDRO_T_DECOUPLE=0`). This replaces the decoupled closure (`=1`), in which the hydro ghost temperature was extrapolated from the top cell and the conduction datum sat at the physical outer face. The outer back-pressure `kOuterPRef`, the Mach-capped velocity and the EOS-closed ghost density are unchanged. The only code change is the preset default in `scenarios/scenario.cpp`. `model_column_2t` keeps `=1`, because its closed two-temperature study was run that way.

**Why.** The single-T boundary is one sentence in the paper: the ghost cells hold a coronal reservoir at `(p_c, T_c)`, and T is treated like the other variables. On the release mesh it gives the same solution as the decoupled closure to within 1 %. The decoupled closure needed a paragraph about two temperatures that contributed nothing to the physics the paper is about.

## 1. What the two closures impose

| | decoupled (previous release) | single T (release) |
|---|---|---|
| hydro ghost T | `T_top` (zero gradient) | `T_c` (fixed) |
| hydro ghost ρ | `ρ_EOS(p_g, T_top)` | `ρ_EOS(p_g, T_c)` |
| conduction Dirichlet point | physical face, `d = Δs/2 = 138 m` | ghost centre, `d = Δs = 276 m` |
| conduction face κ | κ at the reconstructed face state | harmonic mean of κ_top and κ_ghost |
| outer pressure | `kOuterPRef` (same) | `kOuterPRef` (same) |

Formally the single-T face imposes two conditions (p and T) on a subsonic outflow, which has one incoming characteristic. Sections 2 and 3 show why that makes no measurable difference here.

## 2. Why fixed and free hydro ghost temperatures are indistinguishable

**(a) Conduction resets the top-cell temperature far faster than the flow renews it.** Treat the top cell as a heat capacity `ρ c_v Δs` connected to the datum through a conductance `κ_face/d`:

    τ_th = ρ c_v Δs d / κ_face,     c_v = (3/2)(1+x) k_B/m_H = 2.47e4 J kg^-1 K^-1  (x = 1 at the top)

Measured top-cell state (N=500/R4, 4000 s runs):

| | τ_th | τ_th / Δt | top-cell acoustic time Δs/c_s | top-cell advection time Δs/V | advection / τ_th |
|---|---|---|---|---|---|
| decoupled, 4000 s | 0.041 s | 7 | 0.011 s | 29 s | 720 |
| single T, 20 s | 0.084 s | 15 | 0.011 s | 7.1 s | 84 |
| single T, 4000 s | 0.082 s | 15 | 0.011 s | 29 s | 360 |

(`Δt = 5.63e-3 s`, `ρ_top ≈ 2.85e-11 kg m^-3`, `κ_face ≈ 0.655 W m^-1 K^-1`, `Δs = 276 m`.) This estimate is an upper bound, because conduction to the cell below speeds the relaxation further. Gas entering the top cell has its temperature reset about 100–700 times before the outflow carries it through the cell. The top-cell temperature is therefore set by the conduction balance, `T_top = T_c − q d/κ`, whatever the hydro ghost carries.

**(b) So the two ghost temperatures differ only by the small conductive deficit.** Fixed minus free ghost T is exactly `T_c − T_top = q d/κ`: 70 K (0.32 %) with the face datum and 139 K (0.64 %) with the ghost-centre datum at 4000 s. At the fully ionized top, fixed p gives `ρ ∝ 1/T`, so the ghost density and entropy differ by the same fraction and the acoustic impedance `ρc` by half of it.

**(c) At an outflow face the ghost entropy barely enters the flux anyway.** With `u > 0` the contact moves out of the domain, so the exact Riemann solution sampled at `x/t = 0` is on the interior side of the contact and carries the interior entropy. The ghost state enters only through `p*` and `u*`, via its impedance, and the flow Mach number there is 4e-4.

**Measured.** Pinning the hydro ghost T while keeping the face datum (`ISO_OUTER_T_HYDRO_WALL=1`, `outputs/model_column/outer_bc_swap/tfix_pref_4000s`) moves V_top from 9.501 to 9.515 m/s (+0.15 %), with a largest column-wide |ΔV| of 0.014 m/s at 4000 s. The top-cell difference decays in time: 0.075, 0.096, 0.028 and 0.014 m/s at 20, 200, 1000 and 4000 s.

## 3. What does differ: where the Dirichlet datum sits

The single-T closure also moves the conduction datum from the physical face to the ghost centre, so `T_c` is imposed half a cell (138 m) outside the domain. In quasi-steady state the column's enthalpy export fixes `q`, and `T_c − T_top = q d/κ` scales with the stencil distance:

| 4000 s | decoupled | single T |
|---|---|---|
| d [m] | 138.0 | 276.0 |
| q_face [W m^-2] | 0.3335 | 0.3310 (−0.76 %) |
| T_c − T_top [K] | 69.7 | 139.4 (×2.00) |
| V_top [m/s] | 9.501 | 9.419 (−0.86 %) |
| max column-wide \|ΔV\| | — | 0.097 m/s |
| max ΔT/T, Δρ/ρ | — | 0.78 %, 0.80 % |
| M/M(0) | 0.999937 | 0.999938 |

This is a first-order-in-Δs offset in the boundary location. It is permanent in time and halves when the top cell is halved, so it is an accuracy cost that converges away, not an inconsistency. It is far inside the standing N=500 long-duration mass-flux convergence limitation. Split into its two parts at 4000 s: datum location gives −1.0 % in V_top (9.515 → 9.419) and hydro ghost T gives +0.15 % (9.501 → 9.515).

Over 0–1500 s (`visualization/model_column/outer_bc_singleT_vs_release_1500s.mp4`), the V_top gap is −4.0 % at 20 s, during the initial adjustment, and has settled to −0.8 to −0.9 % by 400 s. T differences above 0.1 % are confined to the top ~13 km.

## 4. Where the equivalence fails: no conduction, so no thermostat

Everything in section 2 rests on τ_th being short. With conduction off, τ_th is infinite. The release ghost then holds 22 000 K against the C7 initial top cell (22 674.8 K at N=500, 22 880.3 K at N=1000), and nothing relaxes that 3–4 % entropy jump. The conduction-off hydrostatic check (20 s, CFL 0.50, `outputs/model_column/single_t_bc/`) shows this directly:

| conduction off, 20 s | decoupled | single T, T_c = 22 kK | single T, T_c = T_top(0) |
|---|---|---|---|
| N=500 t=0 max\|F^ρ\| [kg m^-2 s^-1] | 1.7782e-13 @ 2149.1 km | 2.7288e-13 @ 2152.7 km (outer face) | 1.7782e-13 @ 2149.1 km |
| N=500 max\|V\| [m/s] | 5.727e-3 | 1.380e-2 | 5.727e-3 |
| N=1000 t=0 max\|F^ρ\| | 3.9112e-14 | 1.2223e-13 | 3.9110e-14 |
| N=1000 max\|V\| | 1.529e-3 | 6.748e-3 | 1.529e-3 |
| order, max\|F^ρ\| / max\|V\| | 2.19 / 1.90 | 1.16 / 1.03 | 2.19 / 1.90 |
| base-face \|F^ρ\| (N=500) | 1.2457e-13 | 1.2457e-13 | 1.2457e-13 |
| 50 s undriven, N=1000, max\|V\| | 2.655e-3 | 7.328e-3 | 2.655e-3 |

With the ghost temperature matched to the initial top cell, the single-T closure reproduces the decoupled numbers to four digits. The degradation is therefore entirely the mismatch between the initial state and the boundary temperature, not a defect of the interior scheme, and the base face is untouched. **Consequence for the documentation:** the conduction-off hydrostatic convergence evidence (orders 2.19/1.95/1.90 and the six-rung ladder) measures the interior discretization with a boundary temperature consistent with the initial state. It must not be described as the conduction-off behaviour of the release boundary at 22 kK. The release always runs with conduction on, where section 2 applies.

**The historical N=1000 uniform reversal does not reproduce in current numerics.** The original decoupling study (`upper_bc_hydro_temperature_decoupling_recap.md`; N=1000 uniform, 553 m cells, float32, pre-Godunov numerics) found a last-cell velocity reversal under the single-T closure: V_last = −31.6 m/s against +91.3 m/s decoupled at 500 s. The same pair rerun in current release numerics (double, Godunov, `(ln ρ, u, ln p)`; N=1000, `ISO_REFINE_FACTOR=1`, CFL 0.25, 500 s; `outputs/model_column/single_t_bc/uniform_n1000_500s_*`) gives:

| N=1000 uniform, 500 s | single T | decoupled |
|---|---|---|
| V_last [m/s] | 19.206 | 19.557 |
| V_last / V_{N−2} | 1.025 | 1.019 |
| (ρV)_last / (ρV)_mid | 1.076 | 1.067 |
| T_top [K] | 21 577.8 | 21 789.8 |
| cells with V < 0 | 0 | 0 |

There is no reversal in either closure. The single-T gap in V_top is −8.0 %, −4.7 % and −1.8 % at 20, 100 and 500 s, about twice the R4 gap at matching times, as a first-order datum offset with a top cell twice as wide should be. The reversal was therefore a property of the old numerics, not of the single-T closure. The equivalence still depends on the mesh through the size of that offset (∝ Δs_top), so re-measure it before relying on it on a much coarser top cell.

## 5. What was not re-validated

- The wall-temperature sweep (8–100 kK, `upper_wall_temperature_sensitivity_recap.md`) and the outer-closure swap study were run with the decoupled closure. Their conclusions (V_top ∝ T_wall^3.1, the pressure closure decides well-posedness) are closure-independent in mechanism but were not re-run. The 100 kK point is exactly where a pinned ghost starts 4.4x less dense than the top cell.
- The CFL 0.50 sweep was not repeated. The step is limited by the top cell's sound speed and is unchanged (Δt = 5.63e-3 s in both closures), and the 4000 s single-T run completed with `termination=end_time` and zero Godunov fallbacks.

## Release validation

`scripts/release_validation.sh` passes under the new default (`outputs/model_column/release_validation_singleT.log`). The default and explicit paths are bit-identical. The N500 100 s smoke has `termination=end_time` and 2.34e7 exact face solves with 0 fallbacks, with `eps_p` max 0.0002 %, `Q_V` 0.00131, `Q_M` 0.00073, `Q_eff` 7.23e-4, mean V 13.733 m/s, `n(V<0)` 0, `F_eff` 1.0175e-9, `T_top` 21677.2 K, BL 12.150 km, `q_phys` 0.7585 W m^-2 and mass drift −3.407e-4. The previous-closure values were 14.070 m/s, 1.0414e-9, 21838.6 K and 11.870 km. The 4000 s window diagnostic (`util/release/pressure_reconstruction_diag.py --targets 4000`) moves `F_eff` 2.4880e-10 → 2.4695e-10 and mean V 4.741 → 4.690 m/s, and leaves `Q_V`, `Q_M`, BL thickness and `n(V<0)=0` unchanged. The unit suite passes 11647/11647. The 5-step conduction-off regression now runs with `ISO_HYDRO_T_DECOUPLE=1` (§4) and reads max|V| 4.9e-3 m/s.

## Runs

- `outputs/model_column/release_singleT_N500_4000s.*` — canonical release run under the new default (face-flux and outer-conduction sidecars, 20 s frames).
- `outputs/model_column/outer_bc_swap/singleT_ghost_22kK_4000s.*` — the same configuration run explicitly with `ISO_HYDRO_T_DECOUPLE=0` before the default changed.
- `outputs/model_column/twall_sweep/twall_22kK_cfl0.50.*` — the decoupled (previous release) reference.
- `outputs/model_column/outer_bc_swap/tfix_pref_4000s.*` — face datum with pinned hydro ghost T.
- `outputs/model_column/single_t_bc/hse_*`, `undriven50_singleT_n1000.*` — conduction-off checks of section 4.
