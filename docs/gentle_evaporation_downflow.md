# Why the relaxed corona settles to a steady downflow ($V<0$)

In the corona-as-boundary gentle-evaporation runs (`model_gentle`, see `gentle_evaporation_plan.md`, option A) the resolved corona relaxes to a *stationary* profile — $T_i$, $n_e$, and the emission measure stop changing — yet it carries a persistent downflow of a few km/s ($V<0$, since the field-aligned coordinate $s$ increases upward). This note explains why, using only the symbols of the six-equation system in `paper.tex`.

## The balance that governs the corona

Take the charged-fluid energy row of the field-aligned system (`paper.tex`, the six-equation system and its source vector). In the corona the gas is fully ionized, so $U=V$, $T_n=T_i$, the neutral terms and the ion–neutral drag $\alpha$ vanish, and the geometry is a single straight line ($B=1$). The only energy source in $Q^{e_i}$ is the optically-thin radiative loss $n_e n_H\,\Lambda(T)$ (there is **no** ambient volumetric heating in this branch — coronal energy enters *only* as the conductive flux imposed at the top boundary). In steady state ($\partial_t=0$) the row reduces to

$$\partial_s\!\big[(e_i + 2 n_i k_B T_i)\,V\big] \;=\; \partial_s\!\big[(\kappa_e+\kappa_i)\,\partial_s T_i\big]\;-\;n_e n_H\,\Lambda(T).$$

The left side is the divergence of the **advected energy (enthalpy) flux** $(e_i + 2 n_i k_B T_i)\,V$; the first term on the right is the divergence of the **conductive flux** $(\kappa_e+\kappa_i)\,\partial_s T_i$; the last term is the radiative loss. With $e_i = \tfrac32\,2 n_i k_B T_i + \tfrac12\rho_i V^2 + \rho_i\phi_g$, the advected flux is $V$ times $\big[\tfrac52\,2 n_i k_B T_i + \tfrac12\rho_i V^2 + \rho_i\phi_g\big]$ — the $\tfrac52\,2 n_i k_B T_i$ piece is the enthalpy the flow carries.

## A static corona ($V=0$) is not available here

If $V=0$ the advected flux vanishes and the balance becomes $\partial_s[(\kappa_e+\kappa_i)\,\partial_s T_i] = n_e n_H\,\Lambda(T)$: the divergence of the conductive flux must equal the radiative loss at every point. In general that requires a volumetric heating term to make up the difference (the standard "a static corona must be heated" result). Here the corona has no such heating — the conductive flux is fixed only at the top via $q(T)$. For the Model C7 density/temperature profile the flux that $q(T)$ delivers does **not** exactly equal the integrated radiative loss down the column; it falls a little short. So $V=0$ cannot satisfy the balance — the corona is mildly under-heated.

## The downflow is how the budget closes

When conduction under-supplies the radiation, the only remaining term that can carry the difference is the advected (enthalpy) flux $(e_i + 2 n_i k_B T_i)\,V$. Because the imbalance is a net *loss*, the gas must move energy **down** toward the dense transition region, where $n_e n_H\,\Lambda(T)$ is large and the energy is radiated away. Moving enthalpy downward (toward smaller $s$) means $V<0$ — a downflow. Equivalently, in the momentum row the under-heated gas loses pressure support, so the gravitational pull $\rho_i\,\partial_s\phi_g$ slightly wins and the plasma settles. This is **steady, not transient**: it is continuously driven by the energy deficit, and the outer boundary tops the corona back up from its reservoir at the capped rate, so $T_i$ and $n_e$ (hence the emission measure $\int n_e n_H\,\Lambda\,\mathrm{d}s$, dominated by $\int n_e^2\,\mathrm{d}s$) stay fixed while $V<0$ persists.

## Why $|V|$ grows with height

The downward enthalpy flux that must reach the transition region is, to leading order, the same at every height (energy is conserved as it is carried down). That flux is $\approx \tfrac52\,2 n_i k_B T_i\,V$. The charged-fluid pressure $2 n_i k_B T_i$ falls with height, so to keep the product roughly constant, $|V|$ must **grow** toward the tenuous top — matching the relaxed profile, which runs from about $-0.2$ km/s at the base to about $-8$ km/s at the top.

## Role of the two outer-boundary settings

$q(T)$ fixes the conductive flux $(\kappa_e+\kappa_i)\,\partial_s T_i$ entering at the top — the energy input. The Mach-capped outflow limits the velocity at the top to $|V|\le 0.05\,c_{s,i}$ (with $c_{s,i}=\sqrt{2\gamma k_B n_i T_i/\rho_i}$, $\gamma=5/3$), which throttles how fast mass can be exchanged with the reservoir. The steady downflow that emerges is just whatever enthalpy flux is needed to carry off the radiation that conduction alone cannot supply, throttled to the cap. The flow pins near the cap at the top, where the pressure is lowest and the required $|V|$ is highest.

## One-line summary

A corona with no volumetric heating and only a boundary conductive flux $q(T)$ is slightly under-heated; it cannot stand still, so it sheds the leftover radiative loss by carrying enthalpy downward to the transition region — a steady $V<0$ whose magnitude grows with height as the pressure $2 n_i k_B T_i$ drops. Adding a volumetric heating term to balance the radiation is what would pin $V\approx 0$ (the resolved-corona-with-heating variant, v1).


## New explanation:
goal: we want to see chromospheric evaporation is related to heat flux from upper BC.

Problem: pressure gradient does not balance gravity at boundary.

upper BC:
1. pressure: pressure gradient should balance gravity otherwise we will see a downflow from upper boundary. use this at boundary: (p_{i-1} - p_{i+1})/(2dx) = rho_i g.

2. temperature: move the upper boundary downward to the bottom of TR (do not include corona). use a jump temperature at ghost cell: for example $T_{ghost1} = a*Tn, T_{ghost2} = b*T_{ghost1}$. a and b are chosen such that the resulting heat flux q(T) is physically reasonable.

3. at the boundary, because pressure and temperature are imposed, we need to adjust density accordingly.

simulation:
first we disable ionization and recombination and radiative loss and heat flux. only keep single fluid.
we need to start from an isentropic chromopshere. refer to Landau Lifshitz Vol6 Hydrodynamics. T decreases linearly with height. find a reasonable physical IC. we will relax it to steady state and see temperature decreases linearly with height except for near the ghost cells. Ideally V should be near zero in this case (allow tens of m/s).

after this simulation, make a table of what physics is added as summary.

Then we enable heat flux. we shall see T increases. and there's an upflow of velocity, and a downflow in the front. (typical gentle chromosphric evaporation).

only then we plan to add ionization back. don't add back until I approve.

---

## Implementation & results (2026-06-23)

Implemented as the new scenario **`model_isentropic`** ([scenarios/model_isentropic.cpp](../scenarios/model_isentropic.cpp)), env-staged by `ISO_HEAT_FLUX` (0 = Stage 1, 1 = Stage 2). The old corona-as-boundary `model_gentle` is left intact. Run `chromo_main … model_isentropic … no-cooling` (the IC forces single-fluid + all source physics off regardless of CLI args). Knobs: `ISO_NS` (grid cells, **default 400**), `ISO_T_BASE, ISO_N_BASE, ISO_DH, ISO_F_ION, ISO_HEAT_FLUX, ISO_TJUMP_A, ISO_TJUMP_B, ISO_VCAP`.

## The load-bearing discovery — the explicit scheme was not well-balanced

This is the most important result of this work, so it gets a full treatment. **The persistent chromospheric downflow that motivated this whole note — "the pressure gradient does not balance gravity at the boundary" — is not principally a boundary effect. It is a gravity (well-balancing) bug in the explicit MUSCL reconstruction that gives _every_ hydrostatic atmosphere a spurious downward force, everywhere in the domain.**

### Symptom

Stage 1 was set up to be the cleanest possible test: a single-fluid atmosphere whose initial condition is an *analytically exact* hydrostatic, isentropic profile ($T$ linear, $\rho \propto T^{3/2}$, $p \propto T^{5/2}$; see below), with every source term (conduction, cooling, ionization) switched off. Such an atmosphere is an exact stationary solution of the Euler equations with gravity, so it should simply sit still ($V \equiv 0$). Instead the whole column accelerated **downward** as a coherent body, $V$ growing roughly linearly in time at a uniform $\sim\!-180\ \mathrm{m\,s^{-1}}$ per second across the interior, until the top evacuated and the run NaN'd. The atmosphere was in free-fall, not relaxation.

### Diagnosis

Instrumenting the right-hand side directly (`rhs_explicit_state`) at $t = 0$, cell by cell, decomposed the acceleration into the gravity source and the pressure-flux divergence:

- The **gravity source** is correct: $-g$ exactly.
- The **pressure-flux divergence** supplies only **$\tfrac13$** of the support it should. For the (mass-dominant) neutral fluid the flux force was **$0.33g$** instead of $1.0g$, leaving a **net $-0.67g$** downward. For the ion fluid (which carries the electron pressure, so its hydrostatic gradient is $2\rho g$) the flux gave $1.33g$ instead of $2.0g$ — again exactly **$0.67g$ short**.

Three properties pinned it down as a *structural* error, not truncation:

1. **Resolution-independent.** The deficit was $-0.667g$ at ds = 6000 m and still $-0.667g$ at ds = 375 m (a 16× refinement). A consistent scheme's truncation error must vanish as $\Delta s \to 0$; this does not.
2. **Time-step-independent.** Reducing the CFL factor 1000× left the $t = 0$ flux force at $1.3313g$ (ion) unchanged — so it is not the predictor half-step / operator splitting.
3. **The magnitude is exactly $(\gamma-1)g$.** For $\gamma = 5/3$ the deficit is $\tfrac23 g = (\gamma-1)g$ — a clean rational factor, the fingerprint of a specific algebraic error, not numerical noise.

Printing the *reconstructed face pressures* at full precision showed the cell-centred pressures were correct (their centred gradient $= -\rho g$ to round-off), but the neighbour pressures fed into the MUSCL slopes were biased: the recovered $p$ at index $i$ of the shifted array `ip1(xn)` was **higher** than the true cell-$(i{+}1)$ pressure by exactly **$\tfrac23\rho g \Delta s$**, and the `im1` neighbour was **lower** by the same amount.

### Root cause

The TVD-MUSCL reconstruction limits **primitive** variables. It obtains the neighbour primitives by calling `cons2prim` on the **spatially shifted** conserved states (`ip1(xn)`, `im1(xn)`, `ip2`, `im2` in `rhs_explicit_state`). The trouble: this code stores the **total** energy, which includes the gravitational term,

$$E = \tfrac32 p + \tfrac12 \rho v^2 + \rho \phi_g,$$

and `cons2prim` recovers the thermal pressure by subtracting it back off,

$$p = \tfrac23 E - \tfrac13 \rho v^2 - \tfrac23 \rho\,\phi_g(\text{local index}).$$

But `cons2prim` always uses the **cell-local** potential $\phi_g(i) = \tfrac12\big(\phi_{i-1/2}(i)+\phi_{i+1/2}(i)\big)$. When it is handed a *shifted* array — whose value at index $i$ is actually cell $(i{+}1)$'s conserved data, carrying cell $(i{+}1)$'s $\rho\,\phi_g(i{+}1)$ inside its energy — it subtracts the **wrong** potential, $\phi_g(i)$ instead of $\phi_g(i{+}1)$. The recovered pressure is therefore off by

$$\Delta p = \tfrac23 \rho\,\big(\phi_g(i{+}1) - \phi_g(i)\big) = \tfrac23 \rho g \Delta s$$

(using $\phi_g(i{+}1) - \phi_g(i) = g\,\Delta s$ for the uniform grid). So the upper neighbour reads $\tfrac23\rho g \Delta s$ too **high**, the lower neighbour $\tfrac23\rho g \Delta s$ too **low**.

Now form the reconstructed hydrostatic slope. The true slope across a cell is $\mathrm{d}p/\mathrm{d}s\cdot\Delta s = -\rho g \Delta s$. The MUSCL slope uses $(p_{i+1} - p_i)$ and $(p_i - p_{i-1})$; with the biased neighbours,

$$p_{i+1}^{\text{recovered}} - p_i = (p_{i+1} + \tfrac23\rho g \Delta s) - p_i = (-\rho g \Delta s) + \tfrac23\rho g \Delta s = -\tfrac13\rho g \Delta s,$$

and identically for the lower face. **The reconstructed pressure gradient is exactly $\tfrac13$ of the true gradient**, so the face-pressure flux delivers only $\tfrac13\rho g$ of support, leaving the $(\gamma-1)g = \tfrac23 g$ spurious downforce. Because the bias ($\tfrac23\rho g \Delta s$) and the true gradient ($\rho g \Delta s$) both scale with $\Delta s$, their ratio is fixed — hence the error is **resolution-independent**, exactly as observed.

The same $\phi_g$ mismatch also corrupts the **predictor-shifted** states (`prim_xt_state_ip1/im1`, lines that re-shift the half-step state), which feed the 2nd-order corrector — that is a second, independent $\tfrac13$-contribution; both must be corrected to fully balance.

Tellingly, **Stage-D heat conduction does _not_ have this bug**: `rhs_implicit_state` explicitly builds the shifted potentials `phi_g_ip1`/`phi_g_im1` (pinning the boundary ghosts to the face potentials) before forming the neighbour temperatures. Only the *explicit* hydro reconstruction omitted that shift.

### Why this was hidden until now

`model_c7` and the other scenarios run with conduction, numerical diffusivity, optically-thin/thick cooling, TRAC, and a Mach-capped reservoir boundary all active. Those terms damp the spurious $\tfrac23 g$ forcing into a **bounded, driven quasi-steady downflow** rather than a runaway — which is precisely the "persistent $\sim 1\ \mathrm{km\,s^{-1}}$ transition-region downflow" `model_c7` has always shown and which earlier notes attributed to the open boundary and C7's ~20–26 % off-HSE initial condition. Stage 1, by deliberately removing **all** of that damping, exposed the underlying scheme error with nothing left to mask it.

### Fix and verification

A new flag `Grid::well_balanced` (default **false**) gates the correction in `rhs_explicit_state`: when set, each shifted primitive's $p_i$ and $p_n$ are corrected with the $\phi_g$ actually baked into that data's energy (`wb_correct_shifted`), using shifted potential arrays built exactly as `rhs_implicit_state` does (interior → neighbour cell potential; the outer/inner ghosts → the face potentials the boundary used when packing the ghost energies). $p_e$ carries no gravitational term and is shift-invariant, so it is untouched. With the flag on, the neutral net force is **$-0.001g$** (balanced to 0.1 %) at every resolution from ds = 6000 m to 375 m, and Stage 1 relaxes to $V \approx 0$ as it should.

Default-off means **every existing scenario and the full 5950-test regression baseline are byte-for-byte unchanged**; only `model_isentropic` turns it on. It is, however, a strong **candidate to enable globally** — it is the physically correct reconstruction, and globalizing it would remove the spurious downflow from `model_c7`/`model_flare`/`pfss_field_line` (at the cost of re-baselining the tests and the previously-catalogued numbers). That is a deliberate, separate decision (kept off for now).

### Physics added — staged summary

| Physics | Stage 1 (relax) | Stage 2 (heat flux) | Stage 3 (later) |
|---|---|---|---|
| Single-fluid hydro (mass, momentum, energy + gravity) | ✅ on | ✅ on | ✅ on |
| Well-balanced reconstruction (`well_balanced`) | ✅ on | ✅ on | ✅ on |
| Hydrostatic-pressure + EOS-density outer ghost | ✅ on | ✅ on | ✅ on |
| Outer top velocity | V=0 reservoir | Mach-capped outflow | tbd |
| Ion–neutral drag / T-equilibration (Stage B/C) | slaved (single-fluid) | slaved | slaved |
| Heat conduction (Stage D) + numerical diffusivity | ❌ **off** | ✅ **on** | ✅ on |
| Upper-BC temperature jump $T_\text{ghost} = a\cdot T_\text{ref},\ b\,a\cdot T_\text{ref}$ | $a=b=1$ (none) | $a,b>1$ (downward $q(T)$) | $a,b>1$ |
| Radiative cooling (Stage R) + TRAC + vacuum floor | ❌ off | `ISO_COOLING`: off (Stage 2) / on (Stage 2b) | ✅ on |
| Ionization / recombination (Stage E) | ❌ off | ❌ off | ⏸ **awaiting approval** |

### Stage 1 — RESULT: clean steady state (as predicted)
Isentropic IC $T(z)=T_\text{base}(1 - z/H_\text{ad})$ (L&L Fluid Mechanics §4), $T_\text{base}=6500$ K, domain 300 km ($< H_\text{ad}\approx 490$ km so $T>0$). At the **default resolution ns = 400** ($\Delta s \approx 750$ m) it relaxes to a steady state with **$\max|V| = 7.3$ m/s** (a single spike at the bottom boundary cell; the interior is $\approx -2$ m/s), **$T$ linear** ($R^2 = 0.9999$, lapse $-13.2$ K/km $=$ the adiabatic $-g/c_p$) flattening only **near the ghost cells** — exactly as anticipated. Pressure hydrostatic and stable. Plot: `visualization/iso_stage1_relaxation.png`; evolution movie: `visualization/iso_stage1_evolution.mp4`.

**The residual velocity is O($\Delta s$) truncation error — it converges away.** After the well-balancing fix removed the resolution-*independent* $(\gamma-1)g$ force, what remains is the ordinary finite-grid representation of the *curved* isentrope ($p \propto T^{5/2}$): a linear MUSCL reconstruction reproduces the curvature only to first order, leaving a tiny hydrostatic imbalance $\propto\Delta s$ that the Rusanov dissipation holds at a steady (non-accelerating) weak drift. A 32× refinement sweep confirms it halves with $\Delta s$ (no float32 floor):

| ns | $\Delta s$ | interior $\max|V|$ | interior mean $V$ | boundary cell-0 $V$ |
|---:|---:|---:|---:|---:|
| 50 | 6000 m | 16.0 m/s | $-15.2$ | $-55.6$ |
| 100 | 3000 m | 8.07 m/s | $-7.64$ | $-28.2$ |
| 400 (default) | 750 m | 2.14 m/s | $-1.97$ | $-7.27$ |
| 800 | 375 m | 1.15 m/s | $-0.96$ | $-3.74$ |
| 1600 | 187.5 m | 0.73 m/s | $-0.49$ | $-2.07$ |

Convergence order $\log(16.0/0.73)/\log(32) \approx 0.89$ (≈ first-order; the slight sub-linearity is the fixed-ghost boundary residual, which converges marginally slower than the interior). **ns = 400 was adopted as the default** — the interior is already below ~2 m/s, well under any physical evaporation signal (the Stage-2b downflow front is ~1.4 km/s, ~$10^3\times$ larger), while staying cheap enough for rapid staging.

### Stage 2 — RESULT: top heats, but no clean upflow without a sink
The fixed-reference temperature jump ($a\cdot T_\text{ref}$, self-limiting — a jump tracking the *live* heating top $T$ runs away) drives a downward conductive flux that **heats the top strongly** ($T_\text{top}$ 2574 → 25 kK for $a=10, b=5$; → 5.1 kK for the gentle $a=2, b=1.5$). But the response is a downward conduction front sweeping the thin column + a draining downflow — **not** the textbook localized evaporation upflow. Diagnosis: with **no radiative sink** (Stage 2 is heat-flux-only by design) and no resolved coronal volume above, the conductive heat simply heats the whole 300-km column toward the ghost temperature; gentle evaporation is the **conduction-vs-radiation competition** (Antiochos & Sturrock 1978), so the localized TR + upflow needs the radiative sink. Plot: `visualization/iso_stage2_evaporation.png`.

### Stage 2b — RESULT: adding the radiative sink localizes a TR + a downflow front
`ISO_COOLING=1` turns on Stage R (CL2012 optically-thick + optically-thin radiative loss) and TRAC, paired as in `model_c7`. Now the conductive heating and the radiative sink balance: the **base stays cool (~5–6.5 kK)** while the top forms a **steady, sharp transition region (~20 kK over the last ~50 km)** — the heat no longer sweeps the whole column. During the formation a **condensation downflow front** ($-1.4$ km/s at $t\approx 10$ s) propagates down into the chromosphere — the user's "downflow in the front". The evaporation **upflow into the corona is still weak** ($\lesssim$ tens of m/s), limited by (i) the **frozen ionization** ($f = 10^{-3}$, so the heated TR gas cannot ionize toward coronal values and the optically-thin loss $n_e n_H \Lambda$ is suppressed by $\sim 10^3$) and (ii) the **unresolved corona** above (the upflow has nowhere to fill). So: $T$-rise ✅ and downflow front ✅, but the sustained upflow needs the next physics (ionization, and/or a resolved coronal volume). Plot: `visualization/iso_stage2b_evaporation.png`. **→ this is the natural place to add ionization back (Stage 3) — awaiting approval.**

## Carrying the well-balanced fix to the corona loop scenario — GENTLE (2026-06-25)

The well-balanced reconstruction (`Grid::well_balanced`) was applied to the **corona loop scenario** — the resolved chromosphere→TR→corona loop run via `pfss_field_line` with the `GENTLE=1` overlay (conduction-driven gentle evaporation, **no beam**; docs/gentle_evaporation_plan.md §9a/v1). The overlay now sets `grid.well_balanced = true` ([scenarios/pfss_field_line.cpp](../scenarios/pfss_field_line.cpp), env `GENTLE_WELL_BALANCED`, default on). It is scoped to the opt-in GENTLE branch, so the plain `pfss_field_line` scenario, the `PFSS_FLARE` (explosive-beam) overlay, and the 5958-test baseline are byte-for-byte unchanged. Unlike `model_c7`, no boundary closure is changed — the loop heats the corona *inside* the domain with `H(s)` and keeps its open/reflecting BC; only the φ_g reconstruction fix is added.

**Result — the relaxed preflare baseline is quieter; the gentle signature is preserved.** Re-running the L = 25 Mm closed half-loop (`loop_closed_gentle.dat`, ns = 241) with the fix on:

| | original (well_balanced off) | well_balanced ON |
|---|---|---|
| relaxed ρ-weighted mean\|V\| | 112 m/s | **68 m/s** |
| relaxed apex T | 2.26 MK | 2.24 MK |
| evap coronal EM rise | ×2.6 | ×2.53 |
| evap peak upflow | 15.6 km/s (Mach 0.07) | 14.38 km/s (Mach 0.07) |
| evap apex T step | 2.4 → 3.1 MK | 2.38 → 3.12 MK |

The φ_g fix removes the spurious (γ−1)g downforce from the dense chromosphere, halving the mass-weighted baseline downflow (112 → 68 m/s) so the gentle (Antiochos & Sturrock 1978) upflow is measured against a quieter zero. The bounded ~8 km/s TR-cell spike and the few-km/s coronal downflow survive — the latter is the physical energy-budget downflow (the corona is kept up only by `H(s)` + conduction), not the reconstruction artifact. The evaporation fingerprint (subsonic EM rise + apex heating) is essentially unchanged, confirming the fix cleans the baseline without distorting the physics. Movies re-animated: `visualization/gentle_relax_evolution.mp4`, `visualization/gentle_evap_evolution.mp4`; plots `gentle_relax_convergence.png`, `gentle_evap_signature.png`, `gentle_evap_profiles.png`. Commands: `visualize_commands.md` §9a.

## Carrying the upper BC back to the full `model_c7` (2026-06-24)

Having validated `well_balanced` + the hydrostatic/EOS upper BC on the clean single-fluid `model_isentropic` testbed, the same closure was applied to the **full** `model_c7` scenario (real C7 atmosphere, two-fluid, ionization, cooling). It is gated behind a new `model_c7_ic(…, tr_jump_bc=true)` argument and a `Grid::c7_tr_jump_bc` flag, set **only** by the bare `model_c7` scenario — `model_flare` / `analytic_canopy` / `model_gentle` (which share `model_c7_ic` / `model_c7_update_bc`) and the 5958-test baseline are byte-for-byte unchanged.

**The heat input stays the imposed flux, NOT a Dirichlet T-jump.** A literal port of the New-explanation temperature jump was tried first and **fails on the coarse C7 grid**: the 100-cell grid puts the top cell at chromospheric ~6.6 kK (the thin TR rise to 23 kK is unresolved), so any ghost jump hot enough to deliver a physical ~30 W/m² is ~10× hotter than the top, over-conducts ($\kappa_e\propto T^{5/2}$), and rockets the top to 67–140 kK with spurious ±km/s flows — exactly the failure mode model_c7's imposed-flux BC was built to avoid. So the **hybrid** keeps the parts that fix the boundary downflow (well-balanced reconstruction + hydrostatic ghost pressure + EOS ghost density, continuous ghost temperature) and keeps the coronal heat coming from the imposed Neumann flux $q(T)$.

**Result — the spurious downflow is roughly halved, confirming it was partly the well-balancing bug.** Side-by-side (single column, `time_mult=2`, cooling on):

| | original BC | hybrid (`tr_jump_bc`) |
|---|---|---|
| settled interior $V$ mean | $-510$ m/s | $-122$ m/s |
| settled interior $V$ min | $-2000$ m/s | $-1030$ m/s |
| top $T$ | ~200 kK (deep corona) | ~10–20 kK (TR base) |
| stability | hit step cap, $dt\to0.036$ | ran to completion |

The residual ~1 km/s downflow that survives is the **physical** energy-budget downflow this very note opens with (a corona heated only by a boundary conductive flux is slightly under-heated and sheds the leftover radiative loss by carrying enthalpy downward) — *not* the reconstruction bug, which `well_balanced` removes. The top also relocates from a ~200 kK corona to a ~10–20 kK TR base — a more physically-located "bottom of TR" boundary, consistent with the New-explanation intent. Movie re-animated to `visualization/model_c7_evolution.mp4` (rest → relaxation upflow → mild residual downflow). Jump factors are env-tunable (`C7_TJUMP_A`/`C7_TJUMP_B`, default 1 ⇒ continuous), but they only scale the EOS ghost density — they do not drive conduction.

## Extended (photosphere-anchored) Stage 2 — real C7 IC (2026-06-24)

The Stage-2 run was repeated with the **lower boundary moved from h = 1003 km down to the photosphere (h = 0)**. The analytic isentrope cannot span that tall a column — anchored at the C7 photosphere (T ≈ 6583 K) the adiabatic lapse drives T → 0 by ~470 km — so the IC is instead sampled from the **real Model C7 atmosphere** (Avrett & Loeser 2008, Table 26) over h = 0 → 1303 km. The Table-26 rows below h = 1003 km (absent from the embedded table in `model_c7.cpp`) are carried in `model_isentropic.cpp` and enabled with `ISO_C7_IC=1`; the IC reproduces C7's structure (temperature minimum ≈ 4400 K near 560 km rising to ≈ 6.6 kK at the top; n_H from 1.19×10²³ at the photosphere to 3.8×10¹⁸ m⁻³ at the top). All Stage-2 physics and both boundary closures are unchanged (single-fluid, `well_balanced`, conduction + the fixed-reference ghost-T jump a = 10/b = 5, no cooling/ionization); the analytic-isentrope path and the 5958-test baseline are byte-for-byte unchanged. `chromo_main` now labels output heights from a new `Grid::out_base_km` (1003 km for every other scenario, 0 km here) instead of the previously hardcoded +1003 km offset.

**Result — same no-sink signature, now over the full chromosphere.** Exactly as the thin-slab Stage 2, with **no radiative sink** the downward conductive flux from the ghost jump has nowhere to localize, so the heat front **sweeps down the whole column** until the entire chromosphere is heated and then relaxes toward the ghost-driven equilibrium. The T > 20 kK boundary marches all the way from ~745 km down to the photosphere (~9 km by ~3800 s, then the last few km by ~6100 s); a transient mid-column peak reaches ~140 kK (~2900 s) before subsiding, and by t ≈ 6100 s the whole column has settled toward the imposed ghost temperature (Tmax ≈ 66 kK = 10·T_top, base ~9 kK), with the draining downflow growing monotonically to ~−4 km/s (after an early ~17 km/s top transient as conduction switches on — the rarefied top adjusting, decaying within ~100 s). The localized TR + evaporation upflow again require the radiative sink (Stage 2b) — adding it on this extended domain is the natural next step. A clean V≈0 hydrostatic start (C7 temperature, density re-integrated hydrostatically) is available via `ISO_C7_HSE=1` if the early transient needs to be removed. Run with `time_mult=60` (step-capped at 6×10⁵ steps, t ≈ 6100 s). Movie: `visualization/iso_stage2_ext_evolution.mp4`.

### Lower-boundary temperature BC — Dirichlet sink vs Neumann insulating (2026-06-24)

The default lower BC holds the base at the photospheric reservoir temperature (Dirichlet T), so the innermost conduction face conducts heat *out* of the domain once the base warms above it — a heat sink. Switching it to **Neumann** (`ISO_INNER_T_NEUMANN=1`, zero conductive flux ⇒ insulating base, dT/ds = 0) is a **conduction** boundary condition and must be imposed *inside the Stage-D solver* (new `Grid::inner_conduction_neumann`, which drops the innermost face coupling from the tridiagonal). The first, naive attempt instead made the *hydro* inner ghost follow a Neumann temperature (ghost pressure floating with the live cell-0 T); that **reflects the conduction-driven pressure/thermal waves**, and with no radiative sink to damp them they grow without bound — `max|V|` bounces 16 → 4 → 14 → … km/s through the column and the run NaNs at t ≈ 375 s, whereas the Dirichlet base damps the same switch-on transient monotonically to ~1 km/s. Imposing the zero-flux BC in the conduction operator (leaving the hydro reservoir a stable fixed wall) runs cleanly to t ≈ 6100 s.

**Result — the lower-T BC barely matters here.** The run is top-driven (the heat is injected at the ghost jump) and the photospheric base has enormous heat capacity (n_H ≈ 10²³ m⁻³), so insulating the base only raises its final temperature from ~9.4 kK (Dirichlet sink) to ~10.0 kK (Neumann); the equilibrium Tmax (66 kK = 10·T_top), the T > 20 kK front depth (~7 km), and the draining downflow (−4 km/s) are unchanged.

### Relax-then-Stage-2 — adiabatic relaxation before the heat flux (2026-06-24)

The raw C7 profile is not in single-fluid gas-pressure hydrostatic balance (C7's support includes turbulent pressure), so switching on Stage-2 conduction immediately mixes the heat-flux response with an ~3–17 km/s settling transient. To separate them, `ISO_RELAX_TIME=1500` runs a **Stage-1 adiabatic relaxation first** (conduction + jump off, V=0 reservoir top — single-fluid, no cooling/ionization, i.e. isentropic per parcel) for the first 1500 s, then switches on Stage 2 (the fixed ghost-T jump re-anchored to the *relaxed* top temperature). It is time-gated off `grid.sim_time` in the BC update, the same mechanism as model_flare's beam.

**Relaxation:** the C7 IC settles to a clean V≈0 equilibrium — max|V| damps 3.6 km/s (t≈300 s) → ~20 m/s by t≈1200 s and holds. The relaxed profile is *not* C7: adiabatic rearrangement of C7's outward-rising entropy (T varies modestly while ρ falls ~4 decades, so s = ln T − (γ−1) ln n_H rises steeply with height) plus the top reservoir pinning the boundary at the IC top T produces a smooth profile with a **~16 kK temperature bump near h ≈ 1130 km** (and the cool 4.4 kK C7 temperature-minimum dip washed out to ~5.9 kK). This is the genuine near-isentropic hydrostatic equilibrium reachable from C7 under adiabatic dynamics, not a transient.

**Stage 2 from the relaxed state reproduces the no-relax run** — heat front sweeps 1127 → 7 km, interior peaks ~134 kK then equilibrates to ~67 kK (= 10·T_top), draining downflow grows to −4 km/s — confirming the evaporation signature is not an artifact of starting from the off-equilibrium C7 IC. This relax→Stage-2 (+ Neumann base) is the configuration now in `outputs/iso_stage2_ext.txt` / the movie. The relaxation's top-T BC (the IC-pinned reservoir) is the lever if a smoother relaxed top is wanted.

## Why the chromosphere-only run shows no evaporation, and resolving the corona fixes it (2026-06-24)

The chromosphere-only runs above show only a **downflow** (the heat front sweeps the column, then it drains), never the textbook evaporation **upflow**. Antiochos & Sturrock (1978, *Evaporative cooling of flare plasma*) makes clear why: their conduction-driven evaporation is a **coronal** process. Their model base sits at *T > 10⁵ K*, the energy/mass exchange *with the chromosphere is explicitly negligible*, and "cooling is due principally to the evaporation (redistribution) of plasma which is already at coronal temperatures." The large Spitzer flux that drives it (κ_e ∝ T^{5/2}) comes from a resolved 10⁶–10⁷ K corona. The model_isentropic chromosphere-only domain (≤ ~66 kK, corona = a boundary ghost) is simply not in that regime, so a downward conductive flux only heats and drains the chromosphere.

**Resolving the corona produces the evaporation.** `ISO_CORONA=1` extends the single-fluid domain up into a resolved ~1 MK corona (h = 0 → 10 Mm) from the real C7 profile (`scenarios/model_c7.cpp::c7_full_profile`, Avrett & Loeser Table 26 photosphere → 68 Mm). Two additions are forced — neither a dynamical mechanism:

1. **Realistic frozen ionization profile** (n_i = n_e, n_n = n_HI per cell). A 1 MK corona must be ionized for Spitzer conduction; with the uniform frozen f = 10⁻³ the corona is (wrongly) neutral and `kappa_e` collapses to its e–n-limited value (~10⁴× too small at 1 MK), so it cannot conduct. The profile is just the IC — no ionization/recombination network.
2. **TRAC** broadens the steep C7 TR so Spitzer κ_e does not over-conduct and blow up on the coarse grid.

Conduction is on **from t = 0** — a conduction-off relaxation of the steep, stiff chromosphere↔corona TR is unstable (it NaNs at ~250 s with or without an HSE IC); conduction + TRAC are what stabilize it. The HSE IC (`ISO_C7_HSE=1`) gives a clean V = 0 start.

**Result — the upflow appears.** With the ionized corona resolved, the conductive flux drives a clear evaporation upflow: V ≈ 0 in the chromosphere, rising through the TR to **+60–75 km/s in the corona** (T → ~1.2 MK), with a condensation downflow at the TR — the A&S78 evaporation+condensation pattern, the opposite of the chromosphere-only downflow. So the earlier "no evaporation" was a missing-corona artifact, not a contradiction of A&S78. The evaporation here is **vigorous (near-transonic, M ~ 0.6), not "gentle"**: with no radiative sink the 1 MK corona drives hard, then (unmaintained) cools and drains, and the drainage shock goes numerically unstable at ~300 s on the coarse grid — consistent with A&S78's gentle regime requiring M ≪ 1 / the conduction-radiation balance. The movie (`visualization/iso_corona_evolution.mp4`, split x-axis) covers the evaporation phase. A genuinely gentle, sustained, fully-stable version needs the radiative sink (Stage 2b) and sustained coronal heating — i.e. the full `model_gentle` physics.

**The ~300 s cap is physical, not a numerical knob.** An attempt to run 3× longer pinned this down. The death is the *unmaintained* 1 MK corona draining to vacuum at the domain top — a **hydrodynamic** instability (the NaN'd cell has ρ_n ≈ 7×10⁻¹⁶ kg/m³ and V ≈ −13 km/s, i.e. a tenuous, fast downflow whose neutral energy goes negative), not a thermal one. Two levers confirm it: (i) **resolution** — doubling the grid (`ISO_NS=1200`) moves the NaN only from 302 s to 319 s (≈ 5 % for 2× cells: it is tracking a physical drainage timescale, not Δs); (ii) **thermal diffusivity** — tripling it (`ISO_DIFF_MULT=12`) made it slightly *worse* (290 s), since smoothing the temperature does nothing for a velocity/density drainage shock. So the testbed cannot be pushed to ~900 s without giving the corona a way to reach a quasi-steady state — i.e. the radiative sink (and ideally sustained coronal heating). The current `visualization/iso_corona_evolution.mp4` uses the longest stable window obtained (the 1200-cell run, t → 318.6 s); it remains the no-sink testbed by design.

### Adding the full physics — radiation sink + two-fluid + ionization/recomb (2026-06-25, user-directed)

The clean testbed was rerun with the three pieces it deliberately omitted, behind new default-off env gates (the 5958-test baseline is byte-for-byte unchanged): **radiative sink** `ISO_COOLING=1` (Stage R), **two separate fluids** `ISO_TWO_FLUID=1` (ion/neutral drag + T-equilibration via Stage B/C, `single_fluid=false`), and the **ionization/recombination network** `ISO_IONIZATION=1` (Stage E — the frozen C7 split is now only the IC). When the network is on, `model_isentropic_ic` populates `grid.photoionization_rate_i` by the **same Route-B equilibrium inversion as `model_c7_ic`** (blended to the frozen Chae-2021 FAL-C rate above the 1500 km NEQ marker), so the C7 chromosphere is a fixed point and does not drift; the corona's near-zero n_HI is floored in the (n_e/n_n) factor and discarded by the w_eq → 0 coronal weight.

**Float32 forces the base up to h = 1003 km.** The Stage-E channel counts (`dt n_i² α`, `dt n_i n_n S`) and the radiative loss carry bare n² products; at the photosphere (n ~ 10²³ ⇒ n² ~ 10⁴⁶) these overflow float32 (FLT_MAX ≈ 3.4×10³⁸) and NaN at step 0. So the full-physics run anchors at `ISO_H_BASE=1003` (model_c7's validated floor, n ≲ 10¹⁹), keeping the resolved corona 1003 km → 10 Mm. The single-fluid/no-sink testbed is exempt (no n² terms) and keeps its photosphere anchor.

**Result — the sink changes the corona's fate from "evaporate then crash" to "cool and drain, stably."** The full-physics run is **stable for the entire 2445 s** — 8× past the no-sink ~300 s cap — confirming the prediction above that the radiative sink is what lets the corona reach a quasi-steady state instead of draining to a NaN. But with **no sustained coronal heating**, the unmaintained corona now sheds its energy by *radiating* (it could not before), so instead of the vigorous +75 km/s evaporation it **slowly condenses and drains**: a brief early evaporation transient (upflow → ~3 km/s by t ≈ 200 s, after a ~−15 km/s conduction-switch-on downflow at t ≲ 120 s), then a persistent few-km/s coronal downflow; the coronal (h > 2.3 Mm) mass ∫n_e ds falls **−56 %** (9.1 → 4.0 ×10²¹ m⁻²) and EM ∫n_e² ds **−83 %** (1.3×10³⁷ → 2.3×10³⁶) over the run, while the top stays pinned near its C7 value (~0.81 MK, continuous ghost). This is exactly the A&S78 logic: radiation lets a corona lose energy, so an **unheated** corona cools/drains rather than evaporating — the missing ingredient for a *sustained, gentle* upflow is ambient coronal heating (`model_gentle`'s H(s); see gentle_evaporation_plan.md §9a/v1). Movie: `visualization/iso_corona_full_evolution.mp4` (split x-axis); command in `visualize_commands.md`.

### Why the corona must be ionized for it to conduct (the frozen ionization profile, in plain terms)

Heat conduction here is **electron** thermal conduction — the electrons carry the heat. The code uses the Spitzer form ([physics.hpp:45](../physics.hpp#L45))

$$\kappa_e = \frac{9.2048\times10^{-12}\,n_e\,T^{5/2}}{n_e \;+\; 2.836\times10^{-11}\,n_n\,T^{2}}.$$

The denominator counts what slows the heat-carrying electrons down: the $n_e$ term is electron–ion (Coulomb) collisions, the $2.836\times10^{-11}\,n_n T^2$ term is electron–**neutral** collisions.

- **Fully ionized** (real corona): $n_n \to 0$, the neutral term vanishes, and $\kappa_e \to 9.2\times10^{-12}\,T^{5/2}$ — the classic Spitzer conductivity, which is *enormous* at 1 MK.
- **Mostly neutral** (chromosphere, or a corona we *pretend* is neutral): the $n_n$ term dominates the denominator and $\kappa_e \approx 0.32\,(n_e/n_n)\,T^{1/2}$ — far smaller, both because $n_e/n_n \ll 1$ **and** because the temperature scaling collapses from $T^{5/2}$ to $T^{1/2}$.

So whether the corona conducts depends entirely on its **ionization state**. But this is a clean single-fluid testbed — the ionization/recombination network is **off** ([model_isentropic.cpp:262](../scenarios/model_isentropic.cpp#L262)), so the code does not compute $n_e$ and $n_n$; we must *hand it* the split as part of the initial condition. The earlier runs used **one uniform fraction** $f = 10^{-3}$ everywhere (right for the weakly-ionized chromosphere). Applied to a 1 MK corona that makes it 99.9 % **neutral** ($n_n \approx 1000\,n_e$), and the e–n term swamps the denominator: at $T = 10^6$ K, $2.836\times10^{-11}\,n_n T^2 / n_e \approx 2.836\times10^{-11}\cdot1000\cdot10^{12} \approx 3\times10^{4}$ — so $\kappa_e$ is **~10⁴× too small** and the "corona" cannot conduct. No Spitzer flux, no evaporation. It is purely an artifact of pretending the corona is neutral.

The fix is the **frozen ionization profile**: at each cell take the *real* C7 split — set the ion (= electron) density $n_i = n_e^{\rm C7}(h)$ and the neutral density $n_n = n_{\rm HI}^{\rm C7}(h)$ from `c7_full_profile`. C7 is nearly neutral low down ($n_{\rm HI} \gg n_e$) and essentially fully ionized in the corona ($n_e \gg n_{\rm HI}\to 0$), so now the corona conducts via full Spitzer and the chromosphere stays weakly conducting — both physically correct. **"Frozen"** = set once as the IC and never evolved (no ionization dynamics); **"profile"** = varies with height (per cell), not one uniform number.

### Physics included in the resolved-corona run

Run: `ISO_CORONA=1 ISO_CORONA_TOP_KM=10000 ISO_C7_HSE=1 ISO_HEAT_FLUX=1 ISO_NS=600 ./build/chromo_main outputs/iso_corona.txt full no-ionization model_isentropic - 10 no-cooling`

| Mechanism | Status | Notes |
|---|---|---|
| Single-fluid hydro (mass, momentum, energy + gravity) | ✅ ON | `single_fluid=true`; one velocity & temperature |
| Well-balanced reconstruction (`well_balanced`) | ✅ ON | φ_g fix — a hydrostatic atmosphere is a discrete fixed point (V≈0) |
| **Heat conduction (Stage D, Spitzer κ_e)** | ✅ **ON, the entire run** | from t=0 (no relax phase); **this is what drives the evaporation** |
| Numerical diffusivity | ✅ ON | folded into Stage D; `ISO_DIFF_MULT=4` (corona default) |
| TRAC (Johnston 2020) | ✅ ON | broadens the unresolved C7 TR so κ_e doesn't over-conduct on the coarse grid |
| Vacuum floor | ✅ ON | keeps the hot tenuous corona well-conditioned |
| Frozen ionization profile (n_i=n_e, n_n=n_HI per cell) | ✅ ON (as IC only) | from `c7_full_profile`; **required** for coronal Spitzer conduction — not a dynamical network |
| Radiative cooling (Stage R) | ❌ OFF | `no-cooling` — the "no sink"; why the run is vigorous and eventually drains |
| Ionization / recombination (Stage E) | ❌ OFF | clean testbed — awaiting approval |
| Ion–neutral drag / two-temperature (Stage B/C) | ❌ OFF (slaved) | single-fluid |
| Beam heating | ❌ OFF | — |
| Ambient (volumetric) coronal heating | ❌ OFF | corona is unmaintained → cools and drains once evaporation runs |
| Electron temperature ($T_e \neq T_i$) | ❌ OFF | — |
| Upper-BC temperature jump (a, b) | a=b=1 (none) | conductive flux is the **resolved corona's own internal Spitzer flux** down the TR, not an imposed jump |

### Boundary conditions in the resolved-corona run

| Boundary | Quantity | Condition |
|---|---|---|
| **Lower** (photosphere, h = 0) | density ρ_i, ρ_n | Dirichlet — fixed C7 reservoir |
| | pressure / energy | Dirichlet — fixed reservoir (p_i, p_n) |
| | velocity | V = 0 reservoir wall (ghost momentum = 0) |
| | temperature (conduction) | Dirichlet — base held at reservoir T (`ISO_INNER_T_NEUMANN=0`) |
| **Upper** (corona top, h = 10 Mm) | pressure | Hydrostatic — ghost p set so the centered gradient balances gravity, $(p_{ns-2}-p_{\rm ghost})/2\Delta s = \rho g$ ⇒ V=0 fixed point |
| | temperature | Continuous (jump a=b=1): ghost T = top-cell T |
| | density | EOS — $\rho_{\rm ghost} = p\,m/((1+x)\,k_B T)$ from the imposed p, T |
| | velocity | Mach-capped outflow, $\lvert V\rvert \le 0.1\,c_s$ |

