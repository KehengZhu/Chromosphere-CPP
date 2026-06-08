# Plan: physically-grounded boundary conditions (photospheric base + lower-TR top)

## Scope and geometry

The field-aligned 1-D chromosphere code needs two boundaries replaced with physically-grounded closures, replacing the current ad-hoc ghost-cell prescriptions in `scenarios/model_c7.cpp` (`model_c7_update_bc`, line 234).

- **Lower boundary — photosphere.** A $V=0$ steady state that is *physically grounded*: a discrete hydrostatic-equilibrium reservoir fixed to the Model C7 base state, wave-transparent, with numerical (isotropic) diffusion permitted for stability.
- **Upper boundary — lower transition region (TR).** Located at $h \approx 2153$ km, $T \approx 22000$ K — just above the $\kappa_e/\kappa_n$ conductivity crossover ($h = 2149$ km in the writeup, `paper.tex` §Closures), i.e. where electron heat conduction takes over from neutral conduction. Here we impose the density $\rho$, field-aligned velocity $V$, and conductive heat flux $q$ **as functions of temperature $T$**, representing the corona feeding the TR from above.

These choices are consistent with the literature: Rosner, Tucker & Vaiana (1978) independently place the loop/TR base boundary at exactly $T_0 = 2\times10^4$ K with $v\to0$, and the TRAC method (Johnston et al. 2019, 2020) uses a chromospheric floor temperature $T_\mathrm{chrom} = 2\times10^4$ K. So $T\approx22000$ K is the natural top of a chromosphere-only computational column.

## Current state of the code (what we are replacing)

`model_c7_update_bc` (`scenarios/model_c7.cpp:234`) currently does:

- **Inner (photospheric) ghost:** Dirichlet-pinned to the IC snapshot — $\rho_i,\rho_n$ from C7, $V=U=0$, energies from the C7 $(n,T)$ and gravitational potential (`model_c7.cpp:188–203`). This is a reasonable $V=0$ start but is *not* in hydrostatic equilibrium under the code's own EOS: the C7 profile is off discrete-HSE by a nearly uniform 20–26% ($|dp/dh + \rho g|/|\rho g|\approx0.25$, see the comment at `model_c7.cpp:83–91`), which drives a bulk drift on an acoustic-crossing time of order 100 s (writeup §Limitations, `paper.tex:969`).
- **Outer (TR) ghost:** an ad-hoc "$(T, 2T)$ temperature cascade + $V/2, V/4$ velocity halving" (`model_c7.cpp:241–274`). The writeup itself describes this as "mimicking a TR-side outflow" (`paper.tex:880`). This is the prescription we replace with the $\rho(T), V(T), q(T)$ closures below.

The C7 table hardcoded at `model_c7.cpp:57` has 15 points and **tops at $h = 1989$ km, $T = 6674$ K** — *below* the new upper boundary. See the data gap in the Implementation section.

## Supporting papers and the equations they contribute

PDFs in `docs/supporting-papers/`. Equation/page references are from the papers as read (Hansteen 1993 via OCR text `Hansteen1993.txt`, the rest via their text layers / page images).

### Rosner, Tucker & Vaiana 1978 (RTV), ApJ 220, 643 — upper BC backbone

Static loop energy balance (their eq. 3.4, $v=0$, $p=$ const): $E_H + E_R - \nabla\!\cdot\!F_c = 0$, with Spitzer flux $F_c = -\kappa_0 T^{5/2}\,dT/ds$, $\kappa_0\sim10^{-6}$ erg cm$^{-1}$ s$^{-1}$ K$^{-7/2}$ (cgs).

- **Boundary conditions (eqs. 3.5–3.6):** $T(s{=}0) = T_0 = 2\times10^4$ K at the TR base; $F_c(s_\mathrm{max}) = 0$ at the apex; $v\to0$ at both boundaries. (This directly supports our upper-boundary $T\approx22000$ K with near-zero base flux contribution.)
- **Constant pressure across the TR** (the TR is thin vs. the pressure scale height $H = 5\times10^3\,T$ cm): with $p = nkT = $ const, $n \propto p/T$, i.e. **$\rho \propto 1/T$**. This is our $\rho(T)$ closure.
- **Conductive flux as a function of $T$** (integrated energy eq. 3.12, approximation 3.15): $F_c^2(T) \approx f_R(T) - f_H(T)$, where $f_R(T) = (\kappa p^2/2k^2)\int_{T_0}^{T} T'^{1/2} P(T')\,dT'$ (radiation) and $f_H(T) = 2\kappa\int_{T_0}^{T} T'^{5/2} E_H(T')\,dT'$ (heating). The base value $F_c^2(T_0)$ is *negligible* ($<10^4$ vs. $>10^{10}$ in the corona): the TR base does not set the flux — it **absorbs the downward coronal conductive flux and re-radiates it**. Physically this licenses imposing $q$ at the boundary as a coronal-determined value rather than a free one.
- **Radiative loss** $P(T)$ piecewise power law (their Appendix A, eq. A1); compact form used for $f_R$: $\propto T^3$ below $\sim10^{5.1}$ K, $\propto T$ above.
- **Scaling laws:** $T_\mathrm{max} \approx 1.4\times10^3 (pL)^{1/3}$ (eq. 4.3); $E_H \approx 9.8\times10^4\,p^{7/6}L^{-5/6}$ (eq. 4.4).

### Johnston et al. 2019 (ApJL 873, L22) + 2020 (A&A 635, A168) — TRAC, the coarse-grid closure

Field-aligned HD equations (2020 eqs. 2–5) with $F_c = -\kappa_\parallel T^{5/2}\,dT/ds$, $\kappa_0 = 10^{-11}$ W m$^{-1}$ K$^{-7/2}$ (SI), EOS $P = 2k_B nT$, optically-thin radiation $n^2\Lambda(T)$. The explicit-conduction stability limit $\Delta t < \min[k_B n (\Delta s)^2/\kappa_\parallel]$ is what makes a resolved TR unaffordable — directly relevant to our 100-cell grid.

- **Energy jump condition at the TR top (2020 eq. 27)** — the single relation that ties our three boundary quantities together:
$$\frac{\gamma}{\gamma-1}P_c v_c \;+\; \tfrac12\rho_c v_c^3 \;+\; F_{c,c} \;=\; Q_\mathrm{trac} - R_\mathrm{trac},$$
i.e. (enthalpy flux + kinetic-energy flux + downward conductive flux) evaluated at $T_c$ equals (integrated heating $-$ integrated radiation) across the unresolved TR. Subscript $c$ = top of the TR/TRAC region.
- **Cutoff temperature $T_c$:** the highest $T$ where the grid is unresolved, $L_R/L_T > \delta = \tfrac12$ with $L_T = T/|dT/ds|$, $L_R = \Delta s$; bounded $T_\mathrm{chrom} \le T_c \le 0.2\,T_\mathrm{peak}$, $T_\mathrm{chrom} = 2\times10^4$ K. Re-evaluated each step (a fixed $T_c$, as in Lionello 2009 / Mikić 2013, fails for strong transients).
- **Broadening below $T_c$** (2020 eqs. 21–26): $\kappa_\parallel'(T) = \kappa_0 T_c^{5/2}$, $\Lambda'(T) = \Lambda(T)(T/T_c)^{5/2}$, $Q'(T) = Q(T)(T/T_c)^{5/2}$. This conserves $\int\kappa_\parallel\Lambda$ and $\int\kappa_\parallel Q$ (proved analytically, 2020 §2.2), so the integrated radiative losses $R_\mathrm{trac}$ and heating $Q_\mathrm{trac}$ — hence the corona's density response — are preserved even on a coarse grid. The heating modification $Q'$ is essential for footpoint-heated cases.
- The $T_c$-location algorithm is inherited from the Johnston et al. 2017 "jump condition" method (where the TR was treated as a discontinuity and eq. 27 imposed directly as a BC).

### Bradshaw & Emslie 2020, ApJ 904, 141 — the $V(T)$ closure

(Title "Scaling Laws for Dynamic Solar Loops"; correct authors are **Bradshaw & Emslie**, not "Viall" — PDF renamed to `BradshawEmslie2020.pdf`.)

- **Enthalpy flux (eq. 13):** $F_E = \big[\tfrac{\gamma}{\gamma-1}P + \tfrac12\rho V^2\big]V$. For subsonic flow the kinetic term ($\propto M^3$) is negligible vs. the thermal enthalpy term ($\propto M$), so $F_E \approx \tfrac{\gamma}{\gamma-1}PV$.
- **Flow-modified RTV** keeps the static functional form with a Mach-number prefactor: $T_M \propto (1+20M)^{1/3}(P_M L)^{1/3}$ (Appendix A, eq. A2). **For $M<0.05$ the correction is $<10\%$** — a slow, near-static outflow boundary is well justified.
- **Boundary velocity:** $V = M\,c_s(T)$ with $c_s = \sqrt{2\gamma k_B T/m_p}$ (apex sound speed). There is no closed-form local $V(T)$ through the TR; $M$ is a single loop parameter. So at the boundary: pick $M$ (slow), and combine with mass-flux continuity (below).
- Conduction $\kappa_0 T^{5/2}$ ($\kappa_0\sim1.7\times10^{-6}$ cgs); radiation $\Lambda = \chi T^{\alpha}$, $\alpha=-1/2$ over $10^5$–$10^7$ K.

### Hansteen 1993, ApJ 402, 741 — time-dependent TR flows (magnitude check + lower-BC template)

(Scanned/image-only; OCR'd to `docs/supporting-papers/Hansteen1993.txt`.)

- Full HD set (eqs. 2–4) with **un-limited Spitzer** $F_c = -\kappa_0 T^{5/2}\,dT/dz$, $\kappa_0 = 1.1\times10^{-6}$ cgs, and an explicit caveat that Spitzer "may have to be revised for the **lower** transition region" — exactly where our upper boundary sits (flux limiter / non-local conduction is a known caveat there).
- **Lower (footpoint) BC template:** fixed footpoint pressure $p_E = 0.7$ dyn cm$^{-2}$ and $T = 1.5\times10^4$ K, with the boundary **transparent to outgoing waves** (characteristic-variable extrapolation, Korevaar & Van Leer 1988; reflection only a few %). A chromospheric heating term $S = 3\times10^{-23}n_H$ holds the footpoint temperature.
- **TR flow magnitudes:** net downflow (redshift) $\sim1$ km s$^{-1}$ in lower-TR lines (C IV, O IV, O VI; Table 1); the velocity is set by the local balance of conductive-flux divergence vs. radiative losses, not a simple $V(T)$. This is the magnitude/sign sanity check for our $V(T)$.

### Pandey et al. 2024, A&A — lower BC grounding + numerical diffusion

- **Discrete hydrostatic stratification:** integrate $dp/dz = -\rho g$ (ideal gas $p = (k_B/\mu m_p)\rho T$) **on the actual numerical grid**, because "the analytical solution of hydrostatic equilibrium is not identical to the numerical solution" — the analytic integral leaves a residual that launches waves. The $V=0$ steady state is then a true fixed point of the *discretized* scheme.
- **Isotropic numerical diffusion is stability-critical:** Spitzer conduction alone is numerically unstable at the TR; an isotropic $\chi\,\nabla^2$-type term is required. Counterintuitive finding: **higher resolution needs *larger* $\chi$** (finer grids resolve steeper TR gradients). Benchmark: $\chi \gtrsim 30\times10^8$ m$^2$ s$^{-1}$ at $\Delta z = 40$ km. This is the literature backing for "allow for numerical diffusion."
- Coronal heating function = balance of the radiative + conduction losses evaluated on the initial grid (not an analytic scale height); the model is run with velocity damping until it "settles to a numerical equilibrium." So $V=0$ is achieved by *settling*, not a hard clamp.

## Refined boundary-condition closures

### Lower boundary — photosphere, physically-grounded $V=0$ steady state

1. **Discrete-HSE reservoir.** Replace the "pin to IC" inner ghost (`model_c7.cpp:188–203`) with a state in *discrete* hydrostatic equilibrium: integrate $dp/ds = -\rho g$ on the face grid using the same finite-difference stencil the solver uses, seeded by C7's base $T(s)$, so that $V=U=0$ is an exact fixed point. This removes the ~25% HSE residual and the associated 100 s drift. (Pandey 2024 §3.1.)
2. **Wave-transparent + diffusive.** Hold $V=U=0$ at the base but make the ghost wave-transparent (characteristic extrapolation, Hansteen 1993 §2.2) and permit isotropic numerical diffusion tuned to $\Delta s$ (Pandey 2024 §3.3; Gudiksen et al. 2011 hyperdiffusion, already cited) so residual imbalances damp rather than ring. The base $(n,T)$ stay fixed to C7 / VAL–FAL (Avrett & Loeser 2008), the observationally-constrained static photosphere.

### Upper boundary — lower TR, impose $\rho(T)$, $V(T)$, $q(T)$ from the corona

1. **$\rho(T)$ — constant-pressure TR (RTV).**
$$\rho(T) = \frac{\mu m_H}{k_B}\,\frac{P_\mathrm{TR}}{T},$$
with $P_\mathrm{TR}$ the TR pressure taken from C7 at $h\approx2153$ km (Avrett & Loeser 2008 Table 26 — see data gap). $\rho \propto 1/T$.
2. **$q(T)$ — downward coronal conductive flux.** $F_c = -\kappa_0 T^{5/2}\,dT/ds$, with the value set either by
   - (a) the RTV integral $F_c^2(T) = f_R(T) - f_H(T)$ evaluated at $T_\mathrm{top}$, or
   - (b) — preferred for the coarse grid — the **TRAC jump condition** (Johnston 2020 eq. 27), which makes $q$ consistent with the radiation/heating of the unresolved TR below the boundary. The writeup already frames the outer BC as "an imposed heat flux" (`paper.tex:623`); this gives that flux a defensible value.
3. **$V(T)$ — mass-flux continuity + slow Mach cap.** From $\rho V B = $ const along the flux tube (here $B$ uniform) with $\rho(T)$ from step 1, capped at a slow Mach number $V = M c_s(T)$, $M \lesssim 0.05$ (Bradshaw & Emslie 2020; near-static regime), magnitude $\sim1$ km s$^{-1}$ consistent with Hansteen 1993. Replaces the $V/2,V/4$ cascade.
4. **(Recommended) TRAC broadening.** Below the adaptive cutoff $T_c$, rescale $\kappa_\parallel, \Lambda, Q$ by the factors above so the 100-cell grid recovers the correct coronal density / enthalpy exchange without resolving the TR.

## Recommended method (ρ_n→0 robustness + the strong-coupling top)

The upper boundary sits in the lower TR (T≈22000 K) where hydrogen is nearly fully ionized, so the neutral density ρ_n becomes very small. A naive two-fluid treatment blows up there: the neutral velocity reconstruction U_n = (ρ_n U_n)/ρ_n divides a near-zero momentum by a near-zero density, producing spurious large neutral velocities, and the collisional momentum/energy exchange becomes stiff as the coupling strengthens. The resolution is to recognize that **the upper boundary is physically single-fluid**: Gómez Míguez et al. 2024 show hydrogen is almost fully ionized for T>10⁴ K and "the plasma behaves as a single-fluid," with charge–neutral drift terms only relevant where the neutral fraction is already a trace. Every closure in the upper-BC plan above (RTV ρ∝1/T, the TRAC q jump condition, Bradshaw & Emslie V=M c_s) is itself a single-fluid construction, so the boundary closure and the gas physics agree: apply the closures to the charge fluid and *slave* the neutral fluid to it.

**(1) Single-fluid closure on the charge fluid; slave the neutral ghost.** Impose ρ(T), V(T), q(T) on the charged fluid (which carries essentially all the mass and energy at T≥10⁴ K). For the neutral ghost, do not impose an independent neutral momentum BC and never reconstruct U_n = ρU_n/ρ_n at the boundary; instead set ρ_n = (1−f_eq(T)) ρ_tot from the local Route-B ionization equilibrium (tiny but strictly positive — no underflow), U_n = V_c (collisionally locked, the strong-coupling limit), and T_n = T_c. This is the two-fluid-consistent way to apply a single-fluid TR closure and retires the ad-hoc V/2,V/4 cascade.

**(2) Solver-level robustness (implemented).** Stage E (`apply_ionization_stage`) clamps f to [F_FLOOR, 1−F_FLOOR] (F_FLOOR=1e-8) so neither ρ_i nor ρ_n underflows, and guards the velocity reconstruction so a species below a few×F_FLOOR of ρ_tot takes the dominant fluid's velocity (collisional locking) instead of the ill-conditioned momentum/density ratio. The drag stage (`apply_drag_stage`) is **already** robust: it works in center-of-mass + drift variables (V_cm, w) with a backward-Euler relaxation w_new = w/(1+Δt λ), and because the collision coefficient α∝n_n cancels the ρ_n in λ_drag, the update locks U_n→V smoothly as ρ_n→0 with no small-Δt restriction. Optional accuracy upgrade: replace the backward-Euler drift relaxation with the exact exponential collisional integrator of Braileanu & Keppens 2022 (their Eq. 25–26, the Hillier et al. 2016 analytic form), w_new = w·exp(−λ Δt) — L-stable for any Δt, exact velocity locking in the strong-coupling limit. Not required for stability (the CM/drift backward-Euler already locks correctly), so deferred.

**(3) The q(T) coronal flux is required, not optional — it balances ionization cooling.** A run of the full Route-B model exposed a distinct instability that is *not* the ρ_n→0 blowup: with the current ad-hoc upper BC, the top-cell pressure collapses by ×80 over ~90 s (2.6e-2 → 3.2e-4 Pa) with accelerating downflow, going singular near t≈104 s. The cause is energy balance: Route B's collisional ionization (S_i + multilevel S_CR) drains χ_H from the electron thermal pool at the hot top, and three-body super-elastic return is negligible at the low top density, so there is a net thermal sink. The current ad-hoc BC supplies no coronal heat input to balance it, so the top cools, loses pressure support, and drains. Diagnostic confirmation: the no-ionization run is stable to 474 s (no χ_H sink), while the ionization run collapses. This makes the upper-BC **q(T) downward coronal conductive flux the physical heating that holds the top against ionization cooling** — so step "Upper BC" q(T) below is load-bearing, not cosmetic, and (with the constant-pressure ρ(T)) is the actual fix for this collapse. Until then, Route-B model_c7 runs are stable only for t ≲ 90 s.

## Implementation status (2026-06, full plan minus TRAC)

The full plan below is **implemented**, with the user's "RTV now, TRAC after" choice for `q(T)`. Summary of what landed and the empirical outcome:

- **Data gap closed.** `MODEL_C7` (`scenarios/model_c7.cpp`) extended from 15 rows (top 1989 km) to 32 rows up to **2153 km / 2.31×10⁴ K**, from Avrett & Loeser 2008 Table 26 (`n_e`, `n_HI`→`n_n`, `T`; cm⁻³→m⁻³). Spline output clamped to positive floors to suppress overshoot across the steep TR.
- **Lower BC.** Discrete-HSE inner ghost: Dirichlet ρ at C7 base, `V=U=0`, ghost pressures `p_ghost = p_0 + ρ_0 g Δs`. **Result: the ~25%/100 s base drift is gone — ρ_i[0] drifts −0.1% over 328 s.**
- **Upper BC.** Retired the ad-hoc (T,2T)/(V/2,V/4) cascade. New: ρ,T Neumann (constant-pressure TR), Mach-capped (M≤0.05) mass-flux outflow, neutral slaved (U=V). `q(T)` is imposed as a **Neumann conductive-flux** outer BC on the ion conduction row (new `Grid::impose_outer_heat_flux`, `outer_heat_flux`), NOT a Dirichlet hot wall — a hot wall over-conducts via κ_e∝T^{5/2} on the coarse grid and drove a spurious supersonic downflow (verified). `q = (2/7)κ₀T_cor^{7/2}/L ≈ 90 W/m²` (RTV quiet-Sun).
- **Isotropic numerical diffusion** (`Grid::numerical_diffusivity`, Pandey 2024) folded into the Stage-D face conductivities as `K_num = χ·C`; model_c7 sets χ≈5.5×10⁷ m²/s (grid-scaled).
- **q(T) ↔ cooling are paired.** `impose_outer_heat_flux = enable_radiative_cooling`: an imposed coronal flux has no steady sink without radiation. So the **physically-complete model_c7 runs with cooling ON**; no-cooling is a diagnostic baseline.
- **Validation.** *cooling + q(T):* stable to the step cap (t=328 s), top supported (P_top ×1.5), base steady, |V|≲3 km/s — the top forms a self-consistent TR/low-corona (~3×10⁵ K) that re-conducts the flux down. *no-cooling (q off):* now stable to 553 s (the old collapse at ~104 s is gone; top drains only ×6, never singular). Movie: `util/visualization/model_c7_evolution.mp4` (cooling+q run, 1001 frames). All 5834 C++ tests pass.
- **TR optically-thin radiation — ADDED (2026-06).** The "known limitation" below (no TR-temperature radiative sink) is now resolved: `physics.hpp::radiative_loss_thin` adds `Q_thin = n_e n_H Λ(T)` (coronal-abundance Λ, peak ≈4×10⁻³⁵ W m³ at log T≈5.1; CHIANTI-class, Klimchuk 2008), stitched to CL2012 by a ~2×10⁴ K tanh switch, summed into Stage R. With it the column self-consistently forms **chromosphere → thin radiating TR → low corona**: the chromosphere proper stays 5–9×10³ K (it is no longer the chromosphere that heats — the >10⁵ K gas is corona *above* the TR). The TR-temperature radiation (column-integrated ~30–60 W/m²) re-radiates the imposed q; the enthalpy flux (the other dominant TR term) was already in the hydro. `q(T)` was retuned from 90→**~29 W/m²** (L_cor 1.5e7→4.5e7), the sweet spot of the stable window q∈[~30,90]: it puts the TR base as high as possible (~1.77×10³ km) with the corona capped ~2×10⁵ K. Below ~10 W/m² the run hits the radiative-loss-peak thermal/condensation instability and fails. Stable to 411 s; 5834 tests pass; movie regenerated. Refs added to bib: Klimchuk2008, BradshawCargill2010, Landi2012 (NEQ Λ underestimate ≲3×). Writeup updated (paper.tex §closures, main.tex cooling appendix).
- **TRAC κ/Λ/Q broadening + eq.27 jump condition — ADDED (2026-06).** Johnston et al. 2019/2020 TRAC implemented: `integrators.cpp::compute_trac_cutoff_T` finds the adaptive cutoff T_c each step (max T where L_R/L_T > δ=1/2, L_T=T/|dT/ds|, clamped [trac_T_chrom=2e4, 0.2 T_peak]); `physics.hpp::trac_broadening_factor` returns ε=(T_c/T)^{5/2} in [T_b,T_c); Stage D multiplies κ_e by ε, Stage R divides the optically-thin loss by ε. This keeps κΛ and κQ invariant → TR-integrated radiation/heating (coronal density response) conserved on the coarse grid = **the eq.27 jump condition enforced analytically** by the broadening (not as a separate discontinuity BC — that's the 2019 approach, incompatible with our finite-volume conduction). Gated `enable_trac = enable_radiative_cooling`. **Cutoff limiter (Johnston 2020 App. A.2) was REQUIRED:** without it a sudden T_c jump (2e4→3.4e4) caused a conduction shock and NaN at ~105 s; rate-limiting |ΔT_c| to ≤3%/step (both directions) makes engagement smooth and stable. With the limiter, the ad-hoc numerical diffusivity was reduced 5.5e7→2.2e7 (TRAC, which conserves the integrals, now does the principled TR broadening; diffusion is just a stability floor). Stable to step cap (406 s), TRAC engages smoothly (T_c lifts to ~2.1e4 when the TR sharpens, then relaxes), chromosphere 5.1–6.2 kK, TR base ~1.74e3 km, corona ~2e5 K, drift −1.5%. 5852 tests pass (+2 TRAC tests: broadening conservation, cutoff detection+limiter). Writeup: paper.tex §bc + main.tex Stage-D both have TRAC paragraphs; Johnston2019/2020 already in bib. Movie regenerated. Added T_c to chromo_main step diagnostic.
- **Finding:** TRAC is near-latent in the steady model_c7 (the TR is already marginally resolved by the grid+diffusion), correctly acting only when the TR sharpens. It does NOT move the TR base to 2153 km — TR location is set by the q-vs-radiation energy balance, not resolution; TRAC makes that balance grid-converged. Pushing the TR to the boundary would need a different q / coronal density, or a longer column. Secondary omitted sinks remain: ambipolar H-ionization energy flux in the lower TR (Fontenla 1990), He ionization (Golding) — minor.
- **Writeup + bib done.** `paper.tex` §4.5 (`sec:bc`) + §6.1 setup and `main.tex` Stage-D BC + Numerical-Results setup rewritten; `reference.bib` gains RosnerTuckerVaiana1978, Hansteen1993, BradshawEmslie2020, Johnston2019/2020, Pandey2024, GomezMiguez2024. Both documents compile clean.

## Implementation plan

1. **Data gap (blocking for $\rho(T)$):** extend the hardcoded C7 table (`model_c7.cpp:57`, currently topping at 1989 km / 6674 K) up to $\ge2153$ km using Avrett & Loeser 2008 Table 26 (the full $(h,T,n_H,n_{H1},n_e)$ table, p. 20–21). Extract $T\approx22000$ K, $n_e$, $n_H$, and hence $P_\mathrm{TR}$ at 2153 km.
2. **Lower BC:** add a discrete-HSE integrator to build the inner-ghost reservoir; make the inner ghost wave-transparent; expose the isotropic-diffusion coefficient as a grid-scaled parameter.
3. **Upper BC:** rewrite the outer-face block of `model_c7_update_bc` (`model_c7.cpp:241–274`) to set ghost $\rho$ from $P_\mathrm{TR}/(2k_BT)$, ghost $V$ from mass-flux continuity (Mach-capped), and the conduction tridiagonal's outer Dirichlet/flux from the RTV/TRAC $q(T)$.
4. **Optional:** implement adaptive $T_c$ + TRAC broadening of $\kappa_\parallel,\Lambda,Q$ in the conduction/source stages.
5. **Validation:** confirm the discrete-HSE base no longer drifts (the 100 s bulk drift should vanish); check TR-base outflow is $\sim1$ km s$^{-1}$; check the conductive flux magnitude against the RTV integral.

## Writeup changes (both `main.tex` and `paper.tex`, per CLAUDE.md)

- Rewrite §4.5 "Boundary and initial conditions" and the §6.1 setup paragraph (`paper.tex:618`, `paper.tex:870–883`), replacing "mimicking a TR-side outflow" / "ghost temperature set to twice the outermost value" with the $\rho(T), V(T), q(T)$ closures and citations.
- Note the consistency of the $h=2153$ km boundary with the $\kappa_e/\kappa_n$ crossover at 2149 km.

## `reference.bib` entries to add

`RosnerTuckerVaiana1978` (ApJ 220, 643), `Hansteen1993` (ApJ 402, 741), `Johnston2019` (ApJL 873, L22), `Johnston2020` (A&A 635, A168), `BradshawEmslie2020` (ApJ 904, 141), `Pandey2024` (A&A). Already present and reusable: `Spitzer1962`, `Gudiksen2011`, `Fontenla1993`, `Gabriel1976`, `Song2023`, `AvrettLoeser2008`, `CarlssonStein2002`.

## Open questions / next actions

- (a) Pull the C7 Table-26 value at 2153 km to fix $P_\mathrm{TR}$.
- (b) Draft the `reference.bib` entries + the §4.5 writeup rewrite.
- (c) Implement the BC code changes in `model_c7_update_bc`.
- Decide between RTV-integral $q(T)$ vs. full TRAC (jump condition + broadening). TRAC is more robust on the coarse grid but is a larger code change.

## References (papers in `docs/supporting-papers/`)

- Rosner, Tucker & Vaiana 1978, ApJ 220, 643 — `RosnerTuckerVaiana1978.pdf`
- Hansteen 1993, ApJ 402, 741 — `Hansteen1993.pdf` (+ OCR `Hansteen1993.txt`)
- Johnston, Bradshaw, et al. 2019, ApJL 873, L22 — `Johnston2019.pdf`
- Johnston, Cargill, et al. 2020, A&A 635, A168 — `Johnston2020.pdf`
- Bradshaw & Emslie 2020, ApJ 904, 141 — `BradshawEmslie2020.pdf`
- Pandey et al. 2024, A&A — `Pandey2024.pdf`
