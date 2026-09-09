# Numerical method {#numerics}

This page describes the discretization. For the equations being discretized see @ref physics_model.

## The release scheme

| Element | Release setting |
| --- | --- |
| Model | 1D field-aligned hydrodynamics on a straight field line, gravity on, `B = 1` |
| State | single-fluid equilibrium mixture, `U = (rho, rho u, E)`; carrier quantities derived, never stored |
| Storage precision | `double` throughout — `chromosphere::Real` is the one knob and every stored physical scalar is declared in terms of it |
| Update path | `U^n -> MUSCL-Hancock/Godunov hydro -> U* -> implicit physical conduction -> U^{n+1}` |
| Splitting | first-order Lie (not Strang) |
| Thermodynamics | Gamma/Saha equilibrium closure with ionization energy, production `Gamma1` table |
| Reconstruction | MUSCL on `(ln rho, u, ln p)`, MC3/Koren limiter, `beta = 2` |
| Face temperature | exact inversion of the same Saha closure, `T(rho, p)` |
| Numerical flux | SWMF-style **exact-Riemann Godunov** flux at frozen composition (`gamma = 5/3`), non-ideal energy carried in the passive offset `E0`, face-local Rusanov fallback |
| Riemann signal speed | frozen-composition `sqrt(5/3 p/rho)`; the CFL uses the same speed |
| Conduction | physical only, structurally: the operator has no TRAC factor and no artificial diffusivity term at all |
| Upper boundary | 22,000 K external conductive reservoir at the physical face |
| Mesh | coarse-equivalent `N = 500` with R4 outer refinement, giving 661 actual cells |
| CFL | 0.50 (validated for this configuration; the hard-coded solver default stays 0.25) |

None of this is selected by an environment variable. The `model_column` scenario sets it in `model_column_ic`, so the code is the single source of truth, and the `Grid` struct defaults stay false so the non-release two-fluid scenarios are unaffected.

## Finite-volume discretization

The mesh is a static non-uniform finite-volume grid of cell widths `Grid::ds_i`, with the cell centre at the geometric midpoint of its two faces. "Static local refinement" means exactly that: the mesh is generated once at initial-condition time and never regridded, remapped or adapted at runtime. This is not AMR. See `chromosphere::build_static_mesh_faces`.

Grid::broadcast() derives every metric from `ds_i` alone: face arc lengths `s_face`, cell centres `s_i`, and the centre-to-centre distances `ds_iph_i` and `ds_imh_i` with the boundary ghost widths mirrored. On a uniform mesh those all collapse to `ds_i`, `Grid::uniform_mesh` is true, and every spacing-sensitive operator keeps its uniform code path, so a uniform run is byte-for-byte what it was before non-uniform support existed. The graded transition between coarse and fine bands keeps the adjacent cell-width ratio at or below `chromosphere::MESH_MAX_RATIO` (1.1).

Fluxes are differenced with the flux-tube area factor, so the scheme is conservative for variable `B(s)`. The MUSCL reconstruction weights and slope-ratio metrics for a non-uniform mesh are functions of the static mesh alone, so they are built once per mesh generation rather than once per timestep and are invalidated by `Grid::metrics_generation()`.

## Reconstruction

The reconstructed variables are `(ln rho, u, ln p)`.

**Why logarithms.** In a gravitationally stratified atmosphere `rho` and `p` are exponential in height, so a linear TVD slope misrepresents them and a limiter clips the steep gradient down to first order, injecting numerical diffusion that is worst at the dense lower boundary. Limiting the logarithms makes an isothermal hydrostatic column piecewise linear, the slope ratio goes to one, the limiter passes the full slope, and the reconstructed face jump collapses. The signed velocity stays linear — a logarithm is undefined for it and it has no large dynamic range. Faces are exponentiated back before the flux, so conservation is unaffected.

**Why pressure rather than temperature.** The alternative set `(ln rho, u, ln T)` obtains the face pressure afterwards from the nonlinear Saha mapping `p(rho, T)`, so two independently limited variables determine one mechanical quantity, and a mechanically smooth constant-pressure state is *not* reproduced across the partial-ionization transition. Reconstructing pressure and recovering the face temperature by inverting the same authoritative closure removes that manufactured face-pressure mismatch. `ISO_RECONSTRUCTION=lnrho-v-lnt` is retained as a reference configuration for controlled comparison, never for production.

**The limiter.** MC3 / Koren (BATSRUS `ModFaceValue` 'mc3') is the asymmetric third-order (`kappa = 1/3`) monotonized-central limiter. Because it is asymmetric, the two faces of a cell take *different* limited slopes — `chromosphere::flux_lim_mc3_plus` and `chromosphere::flux_lim_mc3_minus`. It is markedly less diffusive than minmod, so it clips an under-resolved transition-region gradient far less. `beta = 2` is the classic Koren value.

## Predictor and numerical flux

A **MUSCL-Hancock** half step advances the reconstructed states to the time-centred level before the flux. Crucially, the predictor applies the **same momentum source** as the corrector — the flux-tube pressure term plus gravity. Applying the source in the corrector only is what left the dominant discrete hydrostatic inconsistency in earlier versions of this scheme.

Both sides of every face are then built through the authoritative equilibrium face builder, so a face state is on the Saha manifold by construction rather than by interpolation.

### The release flux: SWMF-style exact-Riemann Godunov at frozen composition

The release numerical flux follows SWMF `util/CRASH/src/test_godunov.f90::get_godunov_flux`, together with the exact Riemann solver `share/Library/src/ModExactRS.f90` (Toro's exact two-rarefaction/two-shock pressure iteration), ported line for line to `src/single_fluid/exact_rs.hpp`. Per face:

1. the **ideal-gas** Riemann problem is solved exactly at a fixed `gamma = 5/3`;
2. every non-ideal part of the internal energy — here the Saha ionization energy, plus this model's gravitational potential — is carried in a passive **specific energy offset** `E0 = (E - 1/2 rho u^2 - p/(gamma - 1))/rho`, taken from whichever side the contact says is upwind and re-added when the star state's total energy is assembled;
3. the self-similar solution is sampled at `x/t = 0` and the flux assembled from the sampled state;
4. any face whose exact solve fails (vacuum guard, non-positive star pressure, iteration budget) or whose sample is inadmissible falls back to the face-local Rusanov flux, and is counted in `chromosphere::GodunovFluxStats`. The driver prints the tally as `godunov.*` lines at the end of every release run and warns on stderr if it is nonzero.

**Why a fixed `gamma = 5/3` is the right index here, and not an approximation being tolerated.** The equilibrium Saha internal energy of this mixture splits exactly as `e_int = p/(5/3 - 1) + e_ion`: a monatomic translational part plus an ionization reservoir. Solving the face Riemann problem at `gamma = 5/3` with `e_ion` riding inertly in `E0` is therefore the **frozen-composition** limit — the ionization state is held fixed while the waves cross the face. That is the appropriate limit when the ionization/recombination relaxation time is long compared with the dynamical time of the disturbance, which is the chromospheric case for hydrogen: the relaxation time is dominated by slow leakage out of the ground state and runs to `10^3`–`10^5 s` in the mid-chromosphere against a chromospheric hydrodynamic timescale of order a minute (Carlsson & Stein 2002; Leenaarts et al. 2007; Leenaarts 2020).

The competing closure — Saha equilibrium re-established *instantaneously within the wave* — instead gives the equilibrium acoustic index `Gamma1`, which the production CRASH table puts at a median `1.090` in this column. **Neither index is universally "wrong"**: `Gamma1` is the equilibrium acoustic response and `5/3` the frozen one, and the release chooses the frozen one on the physical argument above and on the SWMF methodology. The consequence is that the release flux propagates waves about `1.24x` faster than the equilibrium system does, which is why the timestep is sized from the frozen speed as well (see **Timestep control** below).

Measured behaviour on the release model: zero fallbacks over `2.8 x 10^7` face solves, and the star-pressure Newton solve converging in a single iteration at every face. Evidence: `docs/studies/numerics/swmf_godunov_flux_experiment.md`.

### Reference face solvers

Both are retained for regression and controlled numerical comparison. Neither is used in production, and both print an override line naming the release default when selected.

**Roe (`ISO_RIEMANN=roe-local`)** — the previous release flux, and the reference solver against which the Godunov cutover was quantified. It is a 3x3 local characteristic linearization of `U = (rho, rho u, E)` about the **equilibrium** system: its energy characteristic uses the general-EOS response `(dp/de_int)_rho` and its acoustic eigenvalues carry `c^2 = Gamma1 p/rho`. It is a general-EOS Roe-*type* linearization built from tabulated EOS derivatives, **not** a Roe average constructed to satisfy the exact jump condition, so the classical Roe properties (Property U, exact discrete shock capture) are not claimed for it. Whenever the average is not admissible it falls back to local Rusanov on that face only. Selecting it also reverts the timestep to the equilibrium signal speed, so the reference mode stays self-consistent.

**Rusanov (`ISO_RIEMANN=rusanov`)** — local Lax-Friedrichs, and also the automatic per-face fallback of both other solvers. Not for production: its acoustic-scale dissipation is far too large at the very low Mach numbers of this problem and leaves a persistent transition-region velocity ripple.

## Reference-free hydrostatic balance

**The release interior hydrodynamic operator needs no frozen global hydrostatic reference state.** `Grid::eq_wb` and `ISO_EQ_WB` — the equilibrium-reference "delta-form" device that subtracted a frozen explicit-RHS residual every step — have been **retired from `model_column`**. They still exist in `chromosphere.hpp` for the two-fluid research solver, and must never be set on a release Grid.

What replaced them: the predictor source term described above, plus a second-order trapezoidal, EOS-closed lower ghost ladder. Together these cut the raw discrete hydrostatic face mass-flux defect by roughly a factor of 560.

At the release double-precision storage the remaining defect is `1.78e-13 kg m^-2 s^-1` at `N = 500` (20 s, conduction off, release Godunov flux), and its maximum sits in the upper transition region at 2149 km rather than on the base face. Two things follow, and both are stronger statements than the release could make before the storage cutover:

- **The residual has no representation component.** It is about `1e8` times the double round-off floor of the stored state, so it is entirely discretization and boundary closure. Under the previous `float` storage the same quantity was `1.08e-12` and sat at `1.07x` its own round-off floor — that is, it *was* the representation limit.
- **It is second-order convergent.** Refining `N = 500 -> N = 1000` reduces it 4.55x (observed order 2.19), with order 1.95 at the base face and 1.90 in the 20 s maximum residual velocity. In `float` the same refinement made the defect *worse* (apparent order `-1.0`), the signature of per-face ULP noise. The older statement that no clean truncation-error trend survived between the two resolutions was therefore a property of `float` storage, not of the scheme.

What survives is a genuine defect, not noise: in `double` the base face is still dominated by the dissipative part of the flux (`f_diff = 1.29e-13` of an `f_total = 1.25e-13`), i.e. Godunov dissipation acting on a genuinely nonzero reconstructed jump. That is the same localized first-interior-face artifact listed under Known limitations in @ref validation, and it is open. Measurements: `docs/studies/numerics/reference_free_release_recap.md`.


Be precise about what this claims. The release is **not** a well-balanced scheme in the constructive sense. It removes the dominant discretization inconsistencies and leaves a small, resolution-dependent residual. The boundary closures still capture fixed reservoir data once from the initial condition, as any truncated-domain problem must. Evidence: `docs/studies/numerics/reference_free_release_recap.md`.

## Implicit conduction

The conduction stage is a nonlinear backward-Euler solve on the mixture energy row alone. Newton iterates on the cell temperatures with a tridiagonal Jacobian; the effective heat capacity, the conductivity inputs and the residual energy all come from **one fused Saha evaluation per cell per pass** rather than three independent solves.

Mass and momentum are untouched, and the accepted energy is the flux-updated internal energy plus the unchanged kinetic and gravitational parts, so the stage is conservative by construction. Because the solve is unconditionally stable it imposes no timestep constraint — the acoustic CFL condition is the only constraint the release has, since it carries no volumetric source term.

Face conductivities on a refined mesh use the width-weighted series-resistance form `chromosphere::face_conductivity_series` rather than an arithmetic average, because the flux crosses two half-cells in series. On a uniform mesh callers keep the legacy arithmetic average so results are unchanged.

## State storage precision

The conserved state is stored, and the update `U^{n+1} = U^n + dt R` accumulated, in **`double`**. This is part of the numerical method, not a build detail: the release integrates a stratified layer whose mass-loading timescale exceeds the CFL step by `~3.5e7`, there is no compensated summation and no deviation (perturbation-variable) formulation, so the per-step continuity increment has to be representable against the stored `rho`. In `float` it is not — it falls below half a ULP of `rho` through most of the lower chromosphere, and the density there cannot evolve at all. `double` leaves about eight orders of magnitude of margin, and matches SWMF/BATS-R-US, whose shipped builds all set `PRECISION = ${DOUBLEPREC}` and whose `State_VGB` update is a float64 accumulation.

The EOS and conduction solves were always in double regardless of storage. The single knob is `chromosphere::Real`; the diagnostic-only `CHROMO_STATE_FLOAT32` build reverts the storage to `float` and must not be used for production or any quoted number. See @ref architecture for the build flag and @ref validation for what the cutover fixed and what it did not. Evidence: `docs/studies/numerics/state_precision_release_cutover.md`.

## EOS inversion performance

Every decode inverts the caloric closure for `T`, so the inversion is the hot path. Three things keep it cheap:

- **A safeguarded bracket that cannot be left.** With `T1 = p m_H/(rho k_B)` the fully-neutral temperature, `x` in `(0,1)` puts `T` in `(T1/2, T1]` for every physical state, and `dp/dT` at fixed density is strictly positive there. Newton therefore runs inside an exact factor-of-two bracket and bisects on any rejected step.
- **Per-cell temperature hints.** `Grid::eos_T_hint` records the last temperature decoded in each cell by any stage. It is purely an accelerator: it seeds Newton, never changes the converged root or the tolerance, and a stale or absent hint costs only the ordinary safeguarded solve. Measured on an `h1600`, `ns = 1000` column, a seeded call converges in 1.13 iterations with zero bisections against 8.68 iterations and 4.05 bisections unseeded.
- **State-bound decode caching.** Each immutable conserved state is decoded once into a two-ghost-layer chromosphere::MixtureField, and the CFL, source and reconstruction stages share that field. Predicted and post-source states get separate cache generations. The conduction workspace and the MUSCL arrays are owned by the Grid and reused between steps, so there is no function-static state and no repeated hot-loop allocation.

## Timestep control

`chromosphere::mixture_timestep` returns the CFL-limited step from `|u| + c`, where **`c` is the signal speed of the numerical flux in force**, not unconditionally the equilibrium one:

- release Godunov flux — the frozen-composition speed `sqrt(5/3 p/rho)`;
- Roe and Rusanov reference fluxes — the equilibrium Saha speed `sqrt(Gamma1 p/rho)`.

This is what keeps `CHROMO_CFL` meaning the CFL of the scheme actually running. Implemented as `max(Gamma1, 5/3)` on the release path, so the bound stays conservative if the table ever returns an index above the monatomic value.

**Measured effect on the canonical release configuration: none.** The CFL-limiting cell of the `N = 500`/R4 column is the topmost cell (2153 km, ~22 kK, hydrogen fully ionized), where the tabulated `Gamma1` is already `5/3`; the 1.24x speed gap lives in the partial-ionization interior, which is nowhere near the limit. The rule therefore leaves `dt`, the step count (17781 over 100 s) and the end time unchanged on this model, and the residual solution difference was at the storage round-off level of the build it was measured in — the pre-cutover float32 state. **This also corrects an earlier claim** that running the Godunov flux at `CHROMO_CFL=0.50` gave an effective CFL of about 0.62: measured against the frozen wave speeds, the equilibrium-sized step gave an effective CFL of exactly 0.500. The rule matters as a guarantee for other mesh or model shapes, where a partial-ionization cell could become limiting and the step would then be up to 1.24x shorter than the equilibrium rule would give.

`CHROMO_CFL=0.50` is the validated production value for the reduced release configuration; the conservative hard-coded default remains 0.25 and stays the comparison reference. Rerun the sweep before trusting 0.50 on different hardware or a materially different model shape. The original sweep was run with the Roe flux and the equilibrium signal speed, and it carries over to the release flux only because the step is unchanged on this configuration, as measured above. Evidence: `docs/studies/conduction/coarse_model_column_physical_conduction_recap.md`, `docs/studies/numerics/swmf_godunov_flux_experiment.md`.

Two run-control variables matter for long runs. `CHROMO_T_END` sets an absolute stop time in physical seconds and disables the legacy `time_mult`-derived step cap, so a long run cannot be silently truncated. `CHROMO_FRAME_DT` sets the snapshot cadence in physical seconds rather than steps, and is required whenever output is enabled on a long run — without it the wall time is dominated by ASCII formatting. Every run prints `termination=end_time|step_cap|other`; always check that a long run ended with `end_time`. See @ref configuration.

## The experimental two-temperature scheme {#numerics-two-temp}

**Not release numerics.** `src/two_temp/` advances the four-row state `U = (rho, rho u, E, E_e)` with the release's two-stage Lie split, the single scalar conduction solve replaced by a coupled 2x2-block one:

@verbatim
    U^n -> explicit MUSCL-Hancock/Godunov hydro -> U* -> implicit coupled
           (T_e, T_i) conduction + collisional exchange -> U^{n+1}
@endverbatim

**The acoustics are unchanged, and that is a result rather than an approximation.** Both electrons and heavy particles are monatomic, so at frozen composition the total internal energy still splits exactly as

@verbatim
    e_int = e_h + e_e = (3/2) p_h + (3/2) p_e + e_ion = p/(5/3 - 1) + e_ion
@endverbatim

*independently of how the thermal energy is shared between the two pools*. The release SWMF-style exact-Riemann Godunov flux at `gamma = 5/3`, carrying the ionization energy in SWMF's passive specific-energy offset `E0`, is therefore exactly as valid here as for the release; the frozen signal speed is still `sqrt(5/3 p/rho)`; and `two_temp_timestep` is the release CFL rule verbatim. The stiff collisional exchange imposes no step limit because it is solved implicitly.

The extra physical content rides on a **second characteristic at `lambda = u`**. In primitive variables `(rho, u, p, p_e)` the electron energy equation reduces to the adiabatic electron-pressure law `d_t p_e + u d_s p_e + (5/3) p_e d_s u = 0`, and the flux Jacobian has characteristic polynomial `mu^2 (mu^2 - c^2)` with `mu = u - lambda`. The eigenvalues are `{u - c, u, u, u + c}`, and the doubled `u` is not defective: its two-dimensional eigenspace is spanned by the usual entropy/density mode and by a **pure electron-heavy partition mode**. The three rows the release flux acts on are exactly the rows whose characteristic structure is unchanged; the fourth row is advected on the contact, upwinded with the exact mass flux the Riemann solve already produced.

Reconstruction adds a fourth slot: `(ln rho, u, ln p, ln T_e)` with the same MC3/Koren limiter at `beta = 2`. The **total** pressure remains the directly reconstructed mechanical variable, so the electron/heavy split is the only thing the new slot changes. Face states are built algebraically from `(rho, u, p, T_e)` with no temperature inversion at all: `x = x_Saha(n_H, T_e)`, `p_e = x n_H k_B T_e`, `p_i = p - p_e`, `T_i = p_i/(n_H k_B)`. If a reconstruction leaves `p_e >= p` the **electron** share is capped and the event is counted, because the total pressure is authoritative and must be preserved exactly.

**Conduction and collisional exchange** are solved together by Newton on `(T_e, T_i)` with a 2x2-block tridiagonal Jacobian: the two channels never couple across cells, so the sub- and super-diagonal blocks are diagonal and only the 2x2 diagonal block carries the exchange. The accepted energies are the backward-Euler targets rather than an EOS round-trip of the converged temperatures, so the discrete energy balance holds exactly and the antisymmetric exchange cancels between the pools. Mass and momentum are untouched and both conductive fluxes are re-summed into the total-energy row, so the stage conserves total energy by construction.

The electron caloric pair is taken as exact differences of the *public release EOS functions* — `E_e = e_int(rho,T_e) - (3/2) n_H k_B T_e` and likewise for the capacity — so the two solvers can never drift apart thermodynamically, and at `T_e = T_i` the pools re-sum to `e_int`. Inverting `E_e` for `T_e` uses a **safeguarded Newton with guaranteed bisection progress** (`rtsafe`). The guarantee is not decorative: near full ionization the caloric curve has a very sharp knee, because at `x -> 1` the electron capacity collapses to the translational value while just below the knee `dx/dT_e` makes it orders of magnitude larger. A Newton step taken from the flat side overshoots the knee by thousands of kelvin, and an inversion safeguarded *only* against leaving the bracket oscillates between the two bracket ends indefinitely — each iterate landing strictly inside by ~1e-10 K, so the out-of-bracket test never fires and the bracket shrinks by ~1e-10 K per pass. That was an observed run-ending failure (at `rho = 8.26e-11`, `x = 0.999965`), not a hypothetical.

Every guard the solver can trip is counted in `TwoTempStats` and printed as a `twotemp.*` tally at the end of the run: face `p_e` caps, both decode clamps, Godunov fallbacks, the maximum conduction Newton pass count, and the headline `max_relative_decoupling`. A nonzero clamp or fallback count invalidates the run.

### Two-temperature boundary conditions {#numerics-two-temp-bc}

The hyperbolic and parabolic stages need separate accounting, and conflating them is the usual way to get this wrong.

**Hyperbolic stage — count the characteristics.** At the outer face the flow is subsonic outflow (`0 < u < c`), so of `{u - c, u, u, u + c}` three characteristics *leave* the domain and only `u - c` enters. Exactly **one** condition may be imposed from outside, and the release already spends it on the fixed reservoir back-pressure `kOuterPRef`. Both hydrodynamic temperatures — equivalently the entropy mode and the electron/heavy partition mode, which are precisely the two modes riding the doubled `lambda = u` — must therefore **free-float**, and are zero-gradient extrapolated from the live top cell. Imposing a temperature there would over-specify the problem. At the 1600 km base the flow enters the domain subsonically, so three characteristics are incoming, three conditions are admissible, and the scenario's Dirichlet reservoir in density and *both* temperatures with `v = 0` is within budget.

**Parabolic stage — one condition per channel per end.** Conduction is a separate operator and needs its own four conditions:

| Face | Electron channel | Heavy channel | Why |
| --- | --- | --- | --- |
| Outer (2153 km) | **Dirichlet** `T_e = 22,000 K` | **Neumann** (zero flux) | Conduction down from the transition region and corona is **electron**-conducted: the plasma above is fully ionized and Spitzer `kappa_e ~ T^{5/2}` carries it. The heavy channel modelled here is *neutral-hydrogen* conduction, and there are no neutrals above the domain to supply a flux, so zero heavy flux is the consistent statement. `TT_OUTER_TI=dirichlet` runs the controlled comparison. |
| Inner (1600 km) | **Dirichlet** `T_0` | **Dirichlet** `T_0` | The dense base is collisionally equilibrated — the local equilibration time is about 0.09 ms there — so `T_e = T_i = T_0` is a property of the reservoir, not an extra assumption. `ISO_INNER_T_NEUMANN=1` makes both channels insulating instead. |

This is also the point of the study: assigning the 22,000 K conductive reservoir to the electron channel *specifically* is the physically meaningful thing the temperature split buys, where the release necessarily applies it to one lumped temperature.

**Separate conductivities.** `two_temp_kappa_e` is the release Spitzer expression evaluated at `T_e`; `two_temp_kappa_i` is the release neutral-hydrogen expression evaluated at `T_i`. Spitzer *proton* conduction is deliberately **not** added: it is about 2.3 % of `kappa_e` for hydrogen, and omitting it makes `kappa_e + kappa_n` reproduce the release total conductivity *exactly* at `T_e = T_i`, so the experiment changes only which temperature each existing channel acts on. The exchange conductance `g_ei` is the translational electron capacity `(3/2) n_e k_B` times `nu_ei + nu_en`; both channels are kept because the weakly ionized lower chromosphere has `n_HI >> n_e`, where electron-neutral collisions are the faster route to a common temperature.

## The two-fluid scheme

For contrast: the historical solver uses TVD-MUSCL with the symmetric minmod limiter, the Rusanov / local Lax-Friedrichs flux, and a semi-implicit operator-split driver with the stage sequence listed in `src/two_fluid/integrators.cpp`. Its implicit branch is currently disabled at the call site (see the known issue in @ref physics_model). Explicit RK4 is also available for the explicit-only path. Those are two-fluid entry points; the release timestep is `chromosphere::mixture_advance`.

## Rejected, and not to be reintroduced without new evidence

- The full SWMF `ModTransitionRegion` momentum pressure/gravity split. Measurably no better than the predictor source fix alone on this constant-area column, and not generally conservative for variable-area flux-tube geometry.
- A gravity-aware characteristic (non-reflecting) lower boundary. It changed the hydrostatic residual not at all and the evaporation velocity by at most 0.4 percent. The 1600 km face is not an important acoustic reflector for this model.

## See also

@ref validation for what is verified and what is not, @ref release_solver for the API, `docs/studies/numerics/model_column_release_numerics_recap.md` for the full decision record.
