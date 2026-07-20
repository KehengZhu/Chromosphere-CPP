# CRASH γ-only closure for Chromosphere2026 — implementation plan

Status: Stage 0 (offline diagnostic), Stage 1 (stripped pure-H CRASH table + `zAv` check + raw-table validation), the standalone Saha utility, and the Stage-2 γ-table loader (data structure / parser / interpolation tests) are approved to implement; this document is approved. Six advisor review rounds are incorporated (see "Approval gating" for the resolved lists); the plan is judged implementation-ready. The hydro-solver stages (Stage 3+ touching `state.cpp` / `rhs.cpp` / flux reconstruction) are gated on Stage 0 confirming the Saha/γ consistency, on the equilibrium-mixture MUSCL path (Stage 5.5), and on the review-4 fixes below — all now in the plan.

## Goal and non-goals

Goal (第一版): Chromosphere2026 keeps the fluid equations, conduction, heating, cooling, and boundaries, and uses its **own classical-Saha caloric closure** (`x`, and `q_E`/`γ_E` computed analytically, cross-validated against CRASH) — while **CRASH supplies the equilibrium sound-speed index `Γ₁`(ρ,T)** at runtime. The title "γ-only closure" still fits: the one genuinely CRASH-derived runtime quantity is `Γ₁`. A later variant may take `γ_E` from the table too (Stage 2), but then this Goal line and the Runtime-inputs line must both say so.

Non-goals (explicitly out of scope — this is **not** a CRASH EOS port): no runtime read of CRASH pressure or internal energy; no inversion of an `e_table(ρ,T)=e` CRASH energy table; no general `EosState{T,p,e,x,gamma1,cv}`. The ionization fraction `x` is computed in-code by textbook Saha, not read from CRASH.

## Governing closure (total mixture, γ-only)

The closure acts on the **total mixture** `(ρ_tot, e_int,tot, p_tot)`, then splits back to the state rows by `x`.

- `x = x_Saha(ρ_tot, T)` — Chromosphere2026's own textbook pure-hydrogen Saha.
- `p_tot = (1 + x)(ρ_tot/m_H) k_B T` — ideal partial-ionization pressure.
- `e_int,tot = (3/2) p_tot + x_eq n_H χ_H`. Equivalently `q_E = e_int/p` and `γ_E = 1 + 1/q_E`. In 第一版 this caloric closure is computed **analytically** from the in-code Saha solution and cross-validated against the CRASH `Gamma` column; CRASH supplies `Γ₁` at runtime. **This is the total internal energy, NOT a thermal energy** — it already **includes the ionization energy** `x_eq n_H χ_H`, and is named `e_int` (never `e_th`) so no source stage re-adds `x_eq n_H χ_H` and double-counts.
- Temperature inversion: solve `F(T) = (3/2)(1 + x_eq(ρ,T)) n_H k_B T + x_eq(ρ,T) n_H χ_H − e_int,target = 0` (root-find; the analytic form avoids interpolating the CRASH `Gamma`), with `x_eq = x_Saha`, `n_H = ρ/m_H`; seed with the previous-step T (a speedup only, Stage 3).
- Sound speed: one mixture value `c_s² = Γ₁(ρ,T) (p_i + p_n)/(ρ_i + ρ_n)`.

Runtime inputs from CRASH: **`Γ₁(ρ,T)`** (the equilibrium sound-speed index). In 第一版 `q_E=1/(γ_E−1)` (hence `e_int`) is computed analytically from the in-code Saha `x`, and the CRASH `γ_E` column is used only as a startup/test cross-check — see Stage 2.

## Energy-row mapping (approved, explicit)

The state always carries the 7-variable two-fluid layout (`cons::{RHO_I,RHO_N,MOM_I,MOM_N,E_I,E_N,E_E}`), where `E_I` is the *total charged energy* (protons+electrons+KE+φ) and `E_E` a separately-stored electron internal energy that does **not** enter the total energy. The total internal energy maps to the rows as:

```
e_I,int = (3/2) p_i + x_eq n_H χ_H   # charged row: proton+electron translational + ionization energy
e_N,int = (3/2) p_n                  # neutral row: neutral-H translational only
e_E     = (3/2) p_e                  # electron internal energy (diagnostic row)
```
with the **physics** pressures from `x_eq` — `p_i = 2 x_eq n_H k_B T`, `p_n = (1−x_eq) n_H k_B T`, `p_e = x_eq n_H k_B T` — and the kinetic/gravity terms from the **storage** rows `ρ_i = x_row ρ`, `ρ_n = (1−x_row) ρ` (the `x_eq`/`x_row` split, Stage 5.6):
```
E_I = e_I,int + ½ ρ_i v² + ρ_i Φ
E_N = e_N,int + ½ ρ_n v² + ρ_n Φ
```
This partition — neutral-H translational in `E_N`, proton+electron translational **plus the whole ionization buffer** in `E_I` — is the one consistent with the existing `E_I = total charged energy` definition. Note `p_i` uses `x_eq` while `ρ_i` uses `x_row`, so in the trace-floor limit `E_I` is **not** the physical energy of the stored charged mass (carrier-semantics note, Stage 5.6).

**`E_E` must never use CRASH γ_E.** It is free-electron translational energy and stays at 5/3: `e_E = (3/2) p_e`. The end-of-step slaving `p_e = ½ p_i` (`integrators.cpp:1219-1228`) remains valid at a common temperature for pure hydrogen, but the `prim2cons` conversion for `E_E` must keep the fixed `(3/2) p_e` and must **not** call the variable `inv_gm1()`.

## Stages

### Stage 0 — Offline diagnostic (approved to run; no C++)

Extend `util/plot_eos_gamma.py` to walk the C7 column and output `x_C7(h)` vs `x_Saha(h)`, `γ_E(h)`, `Γ₁(h)`, and `C_V,γ^eff(h)` vs `C_V,CRASH(h)`.

The `C_V` panel is an **implementation-consistency verification, not a design fork.** With (i) the CRASH switches off, (ii) the C++ Saha identical to CRASH's, and (iii) the same energy zero-point, `e_γ(ρ,T) = p_Saha/(γ_E−1)` equals the CRASH internal energy *mathematically*, so `∂e_γ/∂T` must equal `C_V,CRASH`. A significant disagreement means a bug: inconsistent Saha, statistical weights, an uncopied CRASH low-T cutoff/floor, wrong units, a misused γ column, or an interpolation/numerical-derivative error. Compute `C_V^eff` by central difference in `ln T` and compare the normalized `C_V/(n_H k_B)` (the table's `Cv` is k_B per heavy particle; the conduction stage uses a volumetric heat capacity).

The official physical conduction capacity is `C_V^eff = (∂e_int/∂T)_ρ`. The prototype form `C_V ≈ (p/T)/(γ_E−1)` is `e/T`, not `de/dT`; it drops `∂x/∂T` and `∂γ_E/∂T`, i.e. exactly the ionization buffering the γ model is meant to carry. It is retained only as a debug/sensitivity branch and is **not** a candidate for physical runs. Stage 0 does not choose between them.

**Stage 0 gating (whole table, not just the C7 track):** require `C_V^eff > 0`, `∂e_int/∂T > 0`, and `Γ₁ > 0` everywhere — these are the monotonic-inversion and positive-sound-speed preconditions the runtime solver depends on. A failure here is a table/Saha bug to fix before any solver work.

### Stage 1 — Consistent table model + Saha↔zAv verification (table approved)

第一版 (the recommended first version): in `util/eos/tabulate_gamma.f90` set
```fortran
UseExcitation        = .false.
UseFermiGas          = .false.
UseCoulombCorrection = .false.
```
so CRASH emits classical pure-hydrogen Saha `γ_E(T,n_H)` and `Γ₁(T,n_H)`, matching the in-code textbook Saha used for `x`. Alternative (not 第一版): keep full CRASH physics, but then the C++ side must table CRASH `zAv` and drop its own Saha; never mix `x_textbook-Saha` with `γ_excitation+Fermi+Coulomb`.

Two additions to the table build:
- **Extend the range** beyond the current `2×10³–5×10⁶ K` to cover all expected states — flare/evaporation reaches `~10⁷ K` — e.g. `10³–10⁸ K` with a density margin, so the high-T end reaches the fully-ionized ideal limit inside the table rather than at a clamp.
- **Still output `zAv`** in 第一版, for offline validation only (the runtime loader does not read it).

New C++ component: `saha_ionization_fraction(ρ,T)` — pure-H textbook Saha, double precision, log-domain to avoid exponential underflow (none exists today; `apply_ionization_stage` is finite-rate and `c7_route_b_photoionization` is the Route-B equilibrium, not bare Saha).

**Verification test (required):** `x_C++,Saha(T,n_H)` vs `zAv_CRASH(T,n_H)` over the **whole grid**, not just the C7 trajectory. Turning the three switches off does not by itself guarantee agreement — CRASH may still differ in ground-state statistical weight handling, the electron spin factor, `UseGroundStatWeight`, a `set_zero_ionization` branch at very low T, the `LogGeMinBoltzmann` floor, or population tolerances. Tolerance requirement: bounded absolute/relative error through the ionization-transition zone; compare `log x` where `x` is tiny so relative error stays meaningful. Only after this passes is the correct claim: *both sides are configured to represent the same classical Saha model and are verified numerically against each other* (not "share one model by construction").

### Stage 2 — γ-table loader

Runtime table = `log T, log n_H, Γ₁` — the **only** quantity interpolated at runtime (the caloric `q_E`/`γ_E` is analytic; see the governing closure and Stage 3). Requirements: double-precision interpolation (cast to `float` at the boundary); bilinear in `(log T, log n_H)`; **no `Zbar≈0/1 → 5/3` override** (with corrections on, neither the neutral nor the ionized edge is exactly 5/3). Build the table at least one margin wider than the expected runtime range (Stage 1); any runtime `Γ₁` out-of-bounds is a **hard error** (below). `Grid::eos_gamma_table` default off ⇒ the original fixed-γ path runs untouched. Reuse idioms: `physics.hpp:270-284 lin_interp` (1-D → new 2-D), the `data_file_parser.cpp` ifstream loader, and an `env_f` hook (`GAMMA_TABLE=path`). The table's `Gamma`/`zAv`/`Cv` columns are loaded only for validation (Stage 0/1), never for production interpolation.

**第一版: analytic `q_E`/`γ_E` from the code's own Saha is authoritative; the table supplies `Γ₁`.** For the stripped classical pure-H model `q_E = e_int/p = 3/2 + xχ_H/((1+x)k_BT)` and `γ_E = 1 + 1/q_E` are computable directly from the in-code Saha `x` (already validated against `zAv`). Using the analytic `q_E` removes `γ_E` interpolation error, makes `e_int` automatically consistent with the C++ Saha, stabilizes the temperature inversion, and simplifies `C_V^eff` — so in 第一版 the runtime authoritative values are `q_E` (analytic) and `Γ₁` (CRASH table, the equilibrium sound-speed index that genuinely needs the table). The table's `γ_E` column is retained only for a **strict startup/test cross-check** against the analytic value. (If instead you want to claim "CRASH supplies both gammas," interpolate `q_E=1/(γ_E−1)` from the table — never `γ_E` itself, since the `1/(γ_E−1)` amplification blows up near `γ_E→1` — and keep the analytic↔table mismatch as a hard startup assertion.) Either way require, at table-build time, `e_int(ρ,T)` strictly increasing in T on every density column (unique inversion root) and `Γ₁>0`.

**Out-of-table policy (decided): hard error on `Γ₁`.** Because `q_E`/`γ_E` is analytic (valid at any positive `(ρ,T)`, so there is no caloric OOB), `Γ₁` is the only quantity that can leave the table. Any runtime `Γ₁(T,n_H)` OOB aborts with the cell index, time, and `ρ,T,e_int`; the IC is verified fully in-table at initialization and the min/max `(T,n_H)` are logged each step. An edge clamp exists only in an explicitly-enabled debug mode and may not be used in production. `Γ₁=5/3` is **not** a general fallback — correct only in the fully-neutral/fully-ionized limits, never for a density overflow while T is still in the ionization transition.

**Mode-compatibility guard (decided): `GAMMA_TABLE` and `ISO_GAMMA` are mutually exclusive.** Gamma-table mode asserts `single_fluid==true`, `enable_Te==false`, `enable_ionization==false`, and aborts (never silently runs) on `ISO_TWO_FLUID=1`, `ENABLE_TE=1`, a finite-rate Stage E, or a simultaneous `ISO_GAMMA`.

**Timing + lifecycle (requires restructuring `model_column_ic`).** Today `model_column_ic` builds the grid and C7/HSE profile and packs the IC (`~159-205`) *before* it sets `single_fluid`/`enable_ionization`/`enable_Te` (`~285-303`) — so "finish IC, then load the table, then build HSE" is impossible as written, and the column would be built with the old C7 composition before gamma mode is even known. Restructure the IC into an explicit lifecycle so the gamma mode and toggles are finalized **before** the C7/HSE construction: **A.** parse all env vars / mode selection; **B.** set Grid physics toggles; **C.** load the γ table; **D.** validate mode compatibility; **E.** build the Saha-HSE profile; **F.** pack the IC; **G.** build inner/outer ghost states; **H.** freeze `eq_state`/`eq_residual`.

### Stage 3 — Narrow mixture closure + read-only decode vs conservative projection

Not a full `EosState`. Two clearly separated operations, so the equilibrium projection is never hidden inside an ordinary `cons2prim`:

```cpp
struct GammaState { double x_eq, gamma_energy, gamma_sound; };
GammaState gamma_state(double rho_total, double T);                       // x_Saha + analytic γ_E + table Γ₁
double temperature_from_rho_eint(double rho_total, double e_int_total,    // root-find on F(T)=0
                                 double T_guess = NAN);                    // guess is a speedup only

// The single physical entry point. Every gamma-table physics / transport / source path consumes this
// (or its explicit fields) and must NEVER re-derive composition from the stored ρ_i, ρ_n:
//   n_e = x_eq n_H, n_HI = (1−x_eq) n_H are physical; ρ_i/m_H = x_row n_H is only a storage carrier.
struct MixtureThermo {
    double rho, T, x_eq, x_row, n_H, n_e, n_HI, p_i, p_n, p_e, gamma1;
};
// read-only: reads totals, inverts T, fills MixtureThermo — does NOT mutate state.
// Takes the gravitational potential the state ACTUALLY carries (see well-balanced note in Stage 7);
// used by flux, RHS, spectral radius, transport, cooling, diagnostics, boundary calculations.
MixtureThermo decode_equilibrium_mixture(rho_i, rho_n, mom_i, mom_n, E_i, E_n, phi_of_this_state);

// conservative: keeps ρ_tot, ρ_tot·v, and E_I+E_N invariant; rebuilds ρ_i=x_row ρ, ρ_n=(1−x_row)ρ,
// V=U=v, and the E_I/E_N/E_E split (p from x_eq, mass from x_row). The only place rows are rewritten.
void project_equilibrium_single_fluid(...);
```
`cons2prim` must not silently mutate `ρ_i, ρ_n` — it is called by many reconstruction and ghost-state paths, and hidden mutation would make well-balancing and the RK stages impossible to reason about.

**Temperature-inversion lifecycle (correctness must not depend on a previous-step T).** The conserved state carries no temperature, and `decode_equilibrium_mixture` runs on shifted-neighbour, MUSCL-reconstructed, ghost, RK-intermediate, and boundary-trial states that have no associated previous-step T. So `T_guess` is only an accelerator, never a precondition. `temperature_from_rho_eint` must be a self-sufficient bracketed solver: take `[T_min, T_max]` from the table range, check `e_int(ρ,T_min) ≤ e_target ≤ e_int(ρ,T_max)`, run a safeguarded Newton that falls back to bisection whenever a step leaves the bracket, and use the previous-step T (if any) only as the initial iterate. Per-real-cell T may be cached, but reconstructed/shifted states must invert independently without a cache.

### Stage 4 — State reconstruction + energy mapping

In `single_fluid=true`, gamma-table mode:
```
ρ_tot = ρ_i + ρ_n ;  v_cm = (ρ_i V + ρ_n U)/ρ_tot ;  μ = ρ_i ρ_n / ρ_tot ;  n_H = ρ_tot/m_H
e_int,proj = E_I + E_N − ½ ρ_tot v_cm² − ρ_tot Φ        # = e_int,old + ½ μ (V−U)²  (see below)
T = temperature_from_rho_eint(ρ_tot, e_int,proj)
x_eq  = x_Saha(ρ_tot, T)                                # physics
x_row = clamp(x_eq, f_min, 1−f_min)                     # storage only
ρ_i = x_row ρ_tot ; ρ_n = (1−x_row) ρ_tot ; V = U = v_cm            # storage rows (numerical carriers)
p_i = 2 x_eq n_H k_B T ; p_n = (1−x_eq) n_H k_B T ; p_e = x_eq n_H k_B T   # physics pressures
```
then the row energies per the approved mapping above. No donor-velocity bookkeeping (V=U=v makes the reconstruction trivial). Touch points: `state.cpp:110-141` (a gamma-table branch), `model_column.cpp:196-205` energy packing.

**The projection must subtract center-of-mass KE, not the two separate KEs.** With `½ρ_iV² + ½ρ_nU² = ½ρ_tot v_cm² + ½μ(V−U)²`, setting `V=U=v_cm` deletes the relative-drift term `½μ(V−U)²`; a conservative projection must **thermalize** it, i.e. invert T from `e_int,proj = E_I+E_N − ½ρ_tot v_cm² − ρ_tot Φ = e_int,old + ½μ(V−U)²`. This is the key distinction from the **read-only** `decode_equilibrium_mixture`, which — when it is merely interpreting a genuine `V≠U` non-equilibrium state — subtracts the two original kinetic energies. Only the conservative projection folds the relative KE into internal energy.

### Stage 5 — Replace the old Stage B/C with one conservative projection (gamma-table mode)

The current single-fluid path lets the explicit hydro step separate `V,U,T_i,T_n`, then rejoins them with infinite-drag Stage B (`V=U`) and a fixed-heat-capacity Stage C (`T_i=T_n`, using `C_i = 2 n_i k_B/(γ−1)`, `C_n = n_n k_B/(γ−1)` at `integrators.cpp:232-233`). Under a variable-γ equilibrium closure those two constant component capacities lose their meaning — the total heat capacity includes the equilibrium ionization response and cannot be split into two constants. The gamma-table branch instead:

1. explicit hydro → predicted state;
2. keep total mass, total momentum, and total energy;
3. one `project_equilibrium_single_fluid()` enforcing `V=U=v`, `T_i=T_n=T`, `x=x_Saha`;
4. skip the old Stage B/C formulas.

`old Stage B + Stage C → conservative single-fluid equilibrium projection`. This stays entirely Chromosphere2026's own fluid handling and avoids the ill-defined component capacities.

**Conservation guard.** The projection must satisfy `(E_I+E_N)_after = (E_I+E_N)_before` to float tolerance (together with `ρ_tot` and `ρ_tot v_cm` invariant) — assert it in tests.

**Drag heating is subsumed by the projection.** Since gamma-table mode skips Stage B/C, the thermalization of the relative-drift KE `½μ(V−U)²` is done once by the conservative projection (Stage 4). Do **not** additionally list drag-generated heating as a separate source in Stage 7 — that would thermalize `½μ(V−U)²` twice.

**Integrator support scope (第一版).** The projection placement is part of the contract, so the mode's integrator support is restricted up front: *gamma-table equilibrium mode initially supports only `advance_Euler_state`, with one conservative projection after the explicit predictor and after each active energy-source stage.* `advance_RK4` and `advance_Euler_explicit_state` must either be rejected with an explicit error when `eos_gamma_table` is on, or run an equilibrium projection before every RK substage — never silently execute the undefined combination (an RK intermediate state that violates `V=U`, `x=x_Saha`, and the equilibrium pressure partition would feed a wrong RHS).

### Stage 5.5 — Equilibrium-mixture MUSCL reconstruction (blocking for the solver)

`rhs_explicit_state` does not build fluxes from the cell-centered state directly: it independently reconstructs `ρ_i, ρ_n, V, U, p_i, p_n, p_e` and calls `prim2cons()` per face. Even when both neighbour cells satisfy `V=U`, `x=ρ_i/(ρ_i+ρ_n)=x_Saha(ρ,T)`, `p_i=2x n_H k_B T`, `p_n=(1−x)n_H k_B T`, MUSCL-limiting each of those slots **separately** pushes the face state OFF the equilibrium manifold: `x_face = ρ_i,face/(ρ_i,face+ρ_n,face) ≠ x_Saha(ρ_face,T_face)`, and the independently reconstructed `p_i,p_n,ρ_i,ρ_n` generally imply different temperatures. The end-of-step projection cannot repair this — the inconsistent face states have already fed the Rusanov flux, the MUSCL half-step predictor, the second-order corrector, and the spectral radius. And RHS has its OWN internal predictor–corrector (`face primitive → prim2cons → flux → half-step conserved predictor → cons2prim → second reconstruction`) that is not part of `advance_RK4`, so the Stage-5 outer-integrator restriction is not sufficient.

Fix: in gamma-table single-fluid mode, reconstruct **equilibrium mixture variables**, not two species primitives. Recommended slot set `{log ρ_tot, v, log T}` (Saha and the γ table both take `(ρ,T)` directly; `{log ρ_tot, v, log p_tot}` is a valid alternative). At each face obtain `ρ_face, v_face, T_face`, then build the equilibrium face state algebraically:
```
x_eq,face  = x_Saha(ρ_face, T_face)                          # physics
x_row,face = clamp(x_eq,face, f_min, 1−f_min)                # storage only
n_H,face   = ρ_face / m_H ;  V_face = U_face = v_face
ρ_i,face = x_row,face ρ_face ;  ρ_n,face = (1−x_row,face) ρ_face
p_i,face = 2 x_eq,face n_H,face k_B T_face ;  p_n,face = (1−x_eq,face) n_H,face k_B T_face ;  p_e,face = x_eq,face n_H,face k_B T_face
E_I,face, E_N,face, E_E,face  per the approved row mapping (chemical + p from x_eq, mass from x_row)
```
so every face state lies on the equilibrium manifold by construction. **Both** the initial face reconstruction and the internal MUSCL predictor–corrector must use this path. It composes with `well_balanced` / `log_reconstruct` / `mc3_limiter` (already primitive-slot limiters) — the limited slots simply become `{log ρ_tot, v, log T}`, and `log ρ_tot` / `log T` reuse the existing log-reconstruction machinery. Touch points: `rhs.cpp` reconstruction + its internal predictor; `state.cpp` `prim2cons` (the face path); the sound speed (Stage 6) then reads a manifold-consistent `p_i+p_n`, `ρ_i+ρ_n`. **This is the one item that must be in the plan before `state.cpp`, `rhs.cpp`, or the flux reconstruction is touched.**

Explicit predictor–corrector path (the internal predictor must re-decode through the mixture decoder, **never** the legacy `cons2prim`):
```
conserved cell state
  ↓ decode_equilibrium_mixture(state, φ_state)
mixture primitive {log ρ_tot, v, log T}
  ↓ MUSCL reconstruction (limits these slots)
on-manifold face state
  ↓ mixture_to_conserved(face, φ_face)
face flux → half-step conservative predictor → predicted conserved state
  ↓ decode_equilibrium_mixture(predicted, φ_cell)
predicted mixture primitive
  ↓ second reconstruction
on-manifold corrected face state
```
Use a dedicated `struct MixturePrimitive { double rho, velocity, temperature; }` (or three `Vec`s), not the seven-slot primitive, so the legacy `LOG_PRIM_VARS={ρ_i,ρ_n,p_i,p_n,p_e}` cannot leak into the gamma-table path.

**Positivity floors move to mixture space.** The legacy face/predictor floor clips `ρ_i,ρ_n,p_i,p_n` independently, which would again push the state off `p_i=2x n_H k_B T`, `p_n=(1−x)n_H k_B T`. In gamma-table mode all face and predictor positivity safeguards act on **total density and mixture internal energy/temperature** (`ρ_tot≥ρ_min`, `T_min≤T≤T_max`, or `e_int≥e_int(ρ,T_min)` for the conserved predictor); the species rows are then rebuilt by the Saha closure. The legacy independent `(ρ_i,ρ_n,p_i,p_n)` floor is bypassed in this mode.

**Common face potential.** At each interface, **both** reconstructed left and right conserved states are packed with the **same** face potential — `phi_g_iph` at (i+½), `phi_g_imh` at (i−½): `E^{L/R}_face = e^{L/R}_int + ½ρ^{L/R} v² + ρ^{L/R} Φ_face`. Never pack the left state with the left cell's `Φ_i` and the right with `Φ_{i+1}`, or the Rusanov dissipation `−½a(U_R−U_L)` would diffuse `ρ_R Φ_{i+1} − ρ_L Φ_i` as a spurious gravitational-energy jump even for two thermodynamically identical states. Cell-centered predicted states use the cell-centered potential; boundary ghost states use the corresponding boundary-face potential. (This is the face-reconstruction landing of the Stage-7 "decode with the state's own Φ" rule.)

**Carry a predecoded face object (do not re-invert T in the flux).** Reconstruction already knows `(ρ, v, T, x_eq, p_*, Γ₁)` on-manifold, so do **not** pack to a seven-row conserved state and then have `cal_flux_state` / `cal_spectral_radius_state` root-find T again (twice per face, with round-trip error). Keep a `struct MixtureFaceState { double rho, velocity, temperature, x_eq, x_row, p_total, p_i, p_n, gamma1; /* + packed conserved rows for the Rusanov jump */ };` and compute `mixture_flux(face)` and `mixture_spectral_radius(face)` from it directly; the Rusanov jump still uses the packed conserved vector, but the physical flux and wave speed skip the decode. (Cell-centered real cells may cache T for performance, never for correctness — Stage 3.) Not a physics blocker, but it removes most of Stage 5.5's cost and complexity.

### Stage 5.6 — Trace-species floor (split `x_eq` from `x_row`)

Saha gives `x→0` at low T and `1−x→0` at high T, while the code still forms `V=ρ_iV/ρ_i`, `U=ρ_nU/ρ_n` and uses `n_i,n_n` in transport — so a minority row can go ill-conditioned in float even when the mathematical `x` is nonzero. A single clamped `x` used for *everything* would break Saha↔CRASH consistency: the table's `q_E=1/(γ_E−1)` and `Γ₁` come from the **un-clamped** equilibrium, so `e_int = p(x_clamped)·q_E(x_Saha)` would no longer be the CRASH internal energy and `C_V^eff, Γ₁` would belong to a different EOS. Split the two roles:

- **`x_eq = x_Saha(ρ,T)`** (un-clamped) drives ALL physics: `p_tot`, `e_int`, `γ_E`, `Γ₁`, `n_e=x_eq n_H`, `n_HI=(1−x_eq)n_H`, `κ_e`, `κ_n`, radiative cooling, HSE, diagnostics.
- **`x_row = clamp(x_eq, f_min, 1−f_min)`** is ONLY the numerical storage split `ρ_i=x_row ρ`, `ρ_n=(1−x_row)ρ`, keeping both rows nonzero (no divide-by-zero / float conditioning). The minority velocity is `v_cm`.

In gamma-table mode every physics function reads `x_eq, n_e, n_HI, T` from the mixture decoder — it must **not** re-derive composition from the stored `ρ_i, ρ_n`. The stored rows are numerical carriers in the extreme-ionization limits, not exact physical species fractions; this is consistent with the mode already forbidding two-fluid dynamics. (Insisting on one clamped `x` for everything instead defines a *clipped* EOS and would require regenerating `γ_E, Γ₁` consistent with the clip — the raw CRASH table could no longer be used.) Bound the storage error by `f_min` and verify the thermodynamics is independent of it (Stage 9 storage-floor-independence test).

**Carrier semantics.** When `x_row ≠ x_eq` (trace-floor limit) `ρ_i = x_row ρ` but `p_i = 2 x_eq n_H k_B T` and `e_I,int = (3/2)p_i + x_eq n_H χ_H`, so `E_I` is *charged translational energy + chemical reservoir stored in the numerical charged row*, not the physical energy of `ρ_i`. Only the row **sums** and the decoded `MixtureThermo` carry direct physical meaning — never compute `E_I/ρ_i` as a charged specific energy.

### Stage 6 — Common mixture sound speed + operator-splitting convergence

`flux.cpp:60-61` (the only raw `gamma_mono` site) → `c_s² = Γ₁ (p_i+p_n)/(ρ_i+ρ_n)`, the same `c_s` for both Rusanov rows in single-fluid mode. (BC Mach caps `model_column.cpp:505`, `model_c7.cpp:693,731` and tests `chromo_tests.cpp:248-249,1400` use the mixture form.)

Note the residual operator-splitting artifact: the explicit hydro still applies the `p_i` and `p_n` gradients to the two momentum rows separately. Their sum is the correct total pressure gradient, but within a substep this produces a small relative velocity that the projection removes. Add a **Δt→0 convergence test**: the pre-projection relative kinetic energy and the post-projection extra thermalization must → 0.

### Stage 7 — Total-energy source updates (conduction, heating, cooling, drag, floors)

Every active gamma-dependent source currently uses `Δp = (γ−1)Δe`, which under variable γ with ionization buffering is only a frozen-γ approximation. In gamma-table mode convert them all to total-energy updates:

1. add/subtract `Δe` to the total internal energy;
2. keep `ρ`;
3. re-solve `e_int(ρ, T_new) = e_int,new`;
4. re-derive `x, p_i, p_n`.

Applies to beam heating, coronal heating, radiative cooling, and numerical/floor energy corrections. (Drag-generated heating is **not** in this list — it is thermalized by the conservative projection, Stage 5.) Project after each source stage, or merge several local energy sources and project once (Lie splitting — document which).

**Conduction is a single mixture solve, not two per-species rows.** In equilibrium single-fluid mode there is one common temperature, so charged and neutral conduction cannot each carry their own `C_V^eff` and then be projected (that duplicates or arbitrarily splits the heat capacity). The single-fluid conduction equation is
```
∂ e_int(ρ,T)/∂t = ∇·[ κ_tot(T) ∇T ] ,   κ_tot = κ_e + κ_n
```
where `n_e = x n_H`, `n_HI = (1−x) n_H` still set `κ_e` and `κ_n` separately (Spitzer κ∝T^{5/2} stays γ-free) but only their **sum** conducts on the single T. Solve it as the nonlinear implicit residual
```
R(T^{n+1}) = e_int(ρ, T^{n+1}) − e_int(ρ, T^n) − Δt ∇·[ κ_tot(T^{n+1}) ∇T^{n+1} ] = 0
```
with `C_V^eff = ∂e_int/∂T` as the Jacobian of the Newton/Picard linearization — **not** a table value dropped into the existing tridiagonal `cfac`. This keeps the conductive energy change consistent with the EOS internal energy. (`integrators.cpp:420-641` is the current per-species backward-Euler-on-T tridiagonal being replaced in this mode.)

**Well-balanced correction must not use a local γ (`rhs.cpp:21`).** The current `d=(γ−1)(φ_src−φ_cell); p−=ρd` relies on the linear `δp=(γ−1)δe` and is invalid for the nonlinear variable-γ closure. The fix is not to substitute a local `γ_E`, but to have the decoder use the gravitational potential the conserved state **actually carries**: a shifted state's energy holds the shifted cell's `ρφ_src`, so pass `φ_src` into `decode_equilibrium_mixture` (the `phi_of_this_state` argument in Stage 3) rather than decoding with the local `φ_cell` and patching via a linear pressure correction. This deletes the `wb_correct_shifted()` algebra in the gamma-table branch and is cleaner in both fixed-γ and variable-γ modes.

Classification of every γ call site (`gm1()/inv_gm1()/half_gm1()/gamma_mono`):

| Category | Sites |
|---|---|
| (a) total-mixture EOS (→ closure) | `state.cpp:110-112,139-141`; `flux.cpp:20-21,48-49,86-87,107`; `rhs.cpp:383-432` |
| (b) conduction C_V → `C_V^eff` | `integrators.cpp:437-441,540,590,513-515` |
| (c) heating/cooling → total-energy update | `integrators.cpp:189-190,695,853,861-862,901,907-908` |
| (d) drag / T-relax → replaced by projection | `integrators.cpp:232-233,274-276` |
| (e) ionization stage (off in this mode) | `integrators.cpp:1072-1130` |
| (f) boundary Mach cap → mixture c_s | `model_column.cpp:505`; `model_c7.cpp:693,731` |
| (g) well-balanced correction | `rhs.cpp:21` |
| (h) IC/BC/ghost energy packing | `model_column.cpp:63-65,203-205,240-242,532-543` — `E_I/E_N` → new row mapping; `E_E` → fixed `3p_e/2`; ghost composition `x=x_Saha`; ghost energy uses the ghost's own `Φ` |
| (i) other scenario BC audit | grep **all** `inv_gm1()` under `scenarios/`; explicitly mark scenarios incompatible with gamma-table mode as unsupported |
| **γ-free — never touch** | Spitzer κ∝T^{5/2}; collisional `3k_B` β (`integrators.cpp:236,279`); `E_E` electron translational (always 5/3) |

Without (h)/(i) the cell solver could be correct while the boundaries keep injecting fixed-γ energy.

**All transport/source physics must consume `MixtureThermo` (`x_eq`), never `ρ_i/m_H`.** In gamma-table mode `ρ_i/m_H = x_row n_H` is a storage carrier, not `n_e`; the physical densities are `n_e = x_eq n_H`, `n_HI = (1−x_eq) n_H`. Every path that currently derives `n_i,n_n` from the stored rows must instead take the decoded `MixtureThermo` (or its explicit `n_e, n_HI, T, p_*`). Audit and convert at least: `κ_e`, `κ_n` (`physics.hpp`), radiative cooling (Stage R), the `n_e`-dependent parts of beam/coronal heating, the collision frequencies `ν_ei/ν_en/ν_in`, TRAC (`compute_trac_cutoff_T`), the outer heat-flux `q(T)`, boundary diagnostics, the Mach caps, and the ionization fraction written to output/visualization. Enforce via the single decode entry point so no path re-derives composition from `ρ_i, ρ_n`.

### Stage 8 — Saha-consistent C7 HSE

Keep `T(h)=T_C7(h)`; replace the frozen `x=n_i/(n_i+n_n)` (`model_column.cpp:174`) with `x=x_Saha`; re-integrate `dp/dh=−ρg` with `p=(1+x_Saha)ρk_BT/m_H` (the existing trapezoid on `ln p`, `:159-185`, gains a per-cell inner Saha iteration since `x` now depends on `ρ`). Re-freeze the `eq_wb` reference (`:396-397`). No CRASH pressure table involved.

**Fix the base anchor** or "re-integrate HSE" is under-determined. Keep C7's base total hydrogen density `n_H,base^C7`, then set `x_base = x_Saha(n_H,base, T_base)` and `p_base = (1+x_base) n_H,base k_B T_base`, and integrate upward from there. (The alternative — keep C7 base pressure and back-solve a Saha-consistent density — is equally valid but must be chosen explicitly; density-anchored is recommended for continuity with the C7 table.)

### Stage 9 — Validation

1. feature off ⇒ original regression suite unchanged (within tolerance);
2. neutral / high-T limits `γ_E, Γ₁ → 5/3`;
3. table-interpolation continuity;
4. `e_int → T → e_int` round-trip;
5. static Saha-HSE column holds (V≈0);
6. small-amplitude acoustic-wave speed matches `Γ₁`;
7. conduction-pulse energy conservation;
8. evaporation run vs fixed-5/3: compare mass flux, T, p;
9. `x_C++,Saha` vs `zAv_CRASH` full-grid (Stage 1);
10. Δt→0 operator-splitting convergence (Stage 6);
11. conservative projection conserves `E_I+E_N`, `ρ_tot`, `ρ_tot v_cm`, and thermalizes `½μ(V−U)²` for a `V≠U` predictor (Stage 4/5);
12. temperature inversion succeeds on reconstructed / ghost states with **no** T cache (self-sufficient bracketed solver, Stage 3);
13. `advance_RK4` / `advance_Euler_explicit_state` are rejected (or projected per substage) in gamma-table mode, never silently run (Stage 5);
14. equilibrium-manifold reconstruction: a face reconstructed from two on-manifold cells is itself on-manifold (`x_face = x_Saha(ρ_face,T_face)`) to tolerance (Stage 5.5);
15. trace-species floor-convergence (`f_min → 0`) — mass/pressure/n_e error vanishes (Stage 5.6);
16. table-resolution / between-node accuracy: compare `γ_E, Γ₁, C_V^eff, T(ρ,e)` from the loaded grid against a 2×-refined `(T,n_H)` table and against direct CRASH calls at randomly sampled off-node points (runtime error lives between nodes, not at them);
17. boundary consistency: every inner/outer ghost satisfies `x=x_Saha(ρ,T)` and `E_I+E_N = e_int + ½ρv² + ρΦ` with the ghost's own `Φ`;
18. face total-flux identity: for each equilibrium face `F_ρi+F_ρn = ρv` and `F_EI+F_EN = (E_tot+p_tot)v` (the split rows preserve single-fluid total conservation);
19. storage-floor independence: at `f_min = 1e-6, 1e-8, 1e-10` the total mass / energy / evaporation mass flux converge and **all `x_eq`-derived quantities** — `p, T, n_e, n_HI, κ_e, κ_n, Q_rad` — are strictly unchanged by the storage floor (the `x_eq`/`x_row` split, Stages 5.6 / 7);
20. analytic `Γ₁` cross-check (Stage 0/1): `Γ₁,analytic-Saha` (from the classical pure-H `p(ρ,T), e(ρ,T)` and their derivatives) vs the CRASH `GammaS` column over the grid — confirms the `GammaS` column semantics, density/temperature normalization, no column-index error, and statistical-weight assumptions consistent across all thermodynamic derivatives;
21. common-face-potential cancellation: two identical mixture face states reconstructed from cells carrying different cell-centered potentials pack to identical conserved face energies under the shared face `Φ`, and the Rusanov gravitational-energy diffusion term vanishes (Stage 5.5).

## Baseline promise

Default-off ⇒ regression tests unchanged **within existing tolerances**. "Byte-identical" is only promised where the original code expressions are left in place; routing an expression through a shared closure can change float operation order, Armadillo expression fusion, temporaries/rounding, and OpenMP details even when the return value is 5/3.

## Approval gating

Approved to start now: Stage 0 diagnostic; the Stage-1 stripped pure-H CRASH table + `zAv` output; the C++-Saha-vs-CRASH-`zAv` full-grid test; the standalone Saha utility (`saha_ionization_fraction`); the Stage-2 γ-table loader; this document.

Resolved in review 1: `e_int` naming; the explicit `E_I/E_N/E_E` mapping (incl. `E_E` fixed at 5/3); physical `C_V^eff` as the official conduction capacity (not chosen by Stage 0); the read-only-decode vs conservative-projection split; Stage B/C → conservative projection; all active heating/cooling/conduction sources → total-internal-energy updates.

Resolved in review 2 (the five solver items + three minor fixes): conservative-projection energy uses **center-of-mass** KE so `½μ(V−U)²` thermalizes (`E_I+E_N` invariant); conduction is a **single nonlinear mixture residual** `R(T^{n+1})=0` with `C_V^eff` as the Jacobian, `κ_tot=κ_e+κ_n`; the well-balanced correction is replaced by decoding with the **state-carried φ** (deletes `wb_correct_shifted` in this mode); the temperature root-find is a **self-sufficient bracketed solver** (previous-T is a speedup only); the mode is **restricted to `advance_Euler_state`** with defined projection points (RK/explicit paths error out or project per substage). Minor: out-of-table high-T uses the analytic Saha asymptote (not strict 5/3); HSE anchored on C7 base `n_H`; drag heating folded into the projection (not a separate source).

Resolved in review 3 (the architectural gap + numerical protections): **Stage 5.5 equilibrium-mixture MUSCL reconstruction** — reconstruct `{log ρ_tot, v, log T}` and rebuild the species/energy rows on-manifold, for both the initial face reconstruction and the internal RHS predictor–corrector (the true blocker before touching the solver); **Stage 5.6 trace-species floor** (review 3 introduced a single-clamped-composition floor; **review 4 superseded it** with the final `x_eq`/`x_row` split — see the review-4/5 entries) with a floor-convergence test; loader interpolates `q_E=1/(γ_E−1)` (not `γ_E`) with a table-build monotonicity/`Γ₁>0` check and an analytic-`γ_E` cross-check/fallback; a **defined out-of-table `Γ₁` fallback** (no blanket 5/3); an init-time **mode-compatibility guard** (`single_fluid`, `!enable_Te`, `!enable_ionization`; `GAMMA_TABLE` vs `ISO_GAMMA` exclusivity); whole-table Stage-0 gating (`C_V^eff>0`, `∂e_int/∂T>0`, `Γ₁>0`); and a between-node table-resolution validation.

Resolved in review 4 (the five blocking fixes + recommendations): the trace floor now splits **`x_eq` (all physics) from `x_row` (storage only)** so the clamp cannot desync the Saha/CRASH EOS; the out-of-table policy is **decided as a hard error** (edge-clamp is debug-only, never `Γ₁=5/3`); `GAMMA_TABLE`↔`ISO_GAMMA` are **decided mutually exclusive**, with the mode guard run *after* the scenario sets its toggles; the **IC/BC/ghost energy packing** (`model_column.cpp:63-65,203-205,240-242,532-543`) and a `scenarios/`-wide `inv_gm1()` audit are added to the edit surface; Stage 5.5 gains the **explicit predictor path** (internal predictor re-decodes via the mixture decoder, a dedicated `MixturePrimitive` type, no legacy `cons2prim`) and moves **positivity floors into mixture space**; and 第一版 makes **analytic `q_E`/`γ_E` authoritative** with the table supplying `Γ₁` (the CRASH `γ_E` column becomes a strict cross-check). Validation items 17–19 (boundary consistency, face total-flux identity, storage-floor independence) added.

Resolved in review 5 (three blocking fixes + clarifications): the `x_eq`/`x_row` split is now applied **consistently** in Stage 4, Stage 5.5, and the energy-row mapping (`p`/chemical from `x_eq`, mass rows from `x_row`), with `GammaState.x_eq` and a `MixtureThermo{...,x_eq,x_row,n_e,n_HI,...}` carrying both; the `model_column_ic` **initialization lifecycle is reordered** (A parse → B toggles → C load table → D validate → E Saha-HSE → F IC → G ghosts → H freeze eq_state) so gamma mode is fixed before HSE construction; a **unified `MixtureThermo` API** is mandated for all transport/source paths (`κ_e,κ_n`, cooling, heating `n_e` terms, collisions, TRAC, outer `q(T)`, Mach caps, output ionization) so none re-derive composition from `ρ_i,ρ_n`. Clarifications: Goal/Runtime-inputs aligned (analytic `γ_E` authoritative, CRASH supplies `Γ₁`); carrier-semantics note for `E_I` in the trace-floor limit; validation 19 extended to all `x_eq`-derived transport coefficients; new validation 20 = analytic `Γ₁` vs CRASH `GammaS`.

Resolved in review 6 (one convention + two cleanups): the **common-face-potential rule** (both left/right reconstructed states packed with the shared `phi_g_iph`/`phi_g_imh`, cell-centered predicted states with the cell potential, ghosts with the boundary-face potential — no spurious Rusanov gravitational-energy diffusion; validation 21) plus the recommended **predecoded `MixtureFaceState`** so the flux/wave-speed skip re-inverting T; the governing closure now states `e_int = (3/2)p + x_eq n_H χ_H` with analytic `q_E`/`γ_E` (the CRASH `Gamma` column is validation-only); Stage 2 runtime table is **`Γ₁`-only** with a hard error on `Γ₁` OOB (no production caloric fallback). The stale review-3 clamp wording is corrected.

Cleared to implement now: Stage 0; Stage 1 (table generation + raw-table validation); the `saha_ionization_fraction` utility; the Stage-2 loader (data structure, parser, interpolation unit tests). The hydro-solver stages (Stage 3–9) are approved once Stage 0 confirms the Saha/γ consistency; the plan needs no further extension — implement and review per stage.
