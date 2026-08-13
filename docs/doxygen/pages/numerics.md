# Numerical method {#numerics}

This page describes the discretization. For the equations being discretized see @ref physics_model.

## The release scheme

| Element | Release setting |
| --- | --- |
| Model | 1D field-aligned hydrodynamics on a straight field line, gravity on, `B = 1` |
| State | single-fluid equilibrium mixture, `U = (rho, rho u, E)`; carrier quantities derived, never stored |
| Update path | `U^n -> MUSCL-Hancock/Roe hydro -> U* -> implicit physical conduction -> U^{n+1}` |
| Splitting | first-order Lie (not Strang) |
| Thermodynamics | Gamma/Saha equilibrium closure with ionization energy, production `Gamma1` table |
| Reconstruction | MUSCL on `(ln rho, u, ln p)`, MC3/Koren limiter, `beta = 2` |
| Face temperature | exact inversion of the same Saha closure, `T(rho, p)` |
| Numerical flux | mixture Roe characteristic flux (3x3 local linearization), face-local Rusanov fallback |
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

Both sides of every face are then built through the authoritative equilibrium face builder, so a face state is on the Saha manifold by construction rather than by interpolation. The numerical flux is the **mixture Roe characteristic flux**: a 3x3 local linearization of `U = (rho, rho u, E)` whose energy characteristic uses the general-EOS response `(dp/de_int)_rho`, not `Gamma1`. Whenever a Roe average is not admissible the scheme falls back to local Rusanov **on that face only**.

Rusanov (`ISO_RIEMANN=rusanov`) is kept for regression and controlled comparison and must not be used in production: its acoustic-scale dissipation is far too large at the very low Mach numbers of this problem and leaves a persistent transition-region velocity ripple.

## Reference-free hydrostatic balance

**The release interior hydrodynamic operator needs no frozen global hydrostatic reference state.** `Grid::eq_wb` and `ISO_EQ_WB` — the equilibrium-reference "delta-form" device that subtracted a frozen explicit-RHS residual every step — have been **retired from `model_column`**. They still exist in `chromosphere.hpp` for the two-fluid research solver, and must never be set on a release Grid.

What replaced them: the predictor source term described above, plus a second-order trapezoidal, EOS-closed lower ghost ladder. Together these cut the raw discrete hydrostatic face mass-flux defect by roughly a factor of 560, down to a level comparable to the estimated float32 round-off scale of the stored state at release resolution — no clean truncation-error trend survives between `N = 500` and `N = 1000`.

Be precise about what this claims. The release is **not** a well-balanced scheme in the constructive sense. It removes the dominant discretization inconsistencies and leaves a small, resolution-dependent residual. The boundary closures still capture fixed reservoir data once from the initial condition, as any truncated-domain problem must. Evidence: `docs/reference_free_release_recap.md`.

## Implicit conduction

The conduction stage is a nonlinear backward-Euler solve on the mixture energy row alone. Newton iterates on the cell temperatures with a tridiagonal Jacobian; the effective heat capacity, the conductivity inputs and the residual energy all come from **one fused Saha evaluation per cell per pass** rather than three independent solves.

Mass and momentum are untouched, and the accepted energy is the flux-updated internal energy plus the unchanged kinetic and gravitational parts, so the stage is conservative by construction. Because the solve is unconditionally stable it imposes no timestep constraint — the CFL condition comes from the acoustic signal speed `|u| + c_s` alone, which is the only constraint the release has, since it carries no volumetric source term.

Face conductivities on a refined mesh use the width-weighted series-resistance form `chromosphere::face_conductivity_series` rather than an arithmetic average, because the flux crosses two half-cells in series. On a uniform mesh callers keep the legacy arithmetic average so results are unchanged.

## EOS inversion performance

Every decode inverts the caloric closure for `T`, so the inversion is the hot path. Three things keep it cheap:

- **A safeguarded bracket that cannot be left.** With `T1 = p m_H/(rho k_B)` the fully-neutral temperature, `x` in `(0,1)` puts `T` in `(T1/2, T1]` for every physical state, and `dp/dT` at fixed density is strictly positive there. Newton therefore runs inside an exact factor-of-two bracket and bisects on any rejected step.
- **Per-cell temperature hints.** `Grid::eos_T_hint` records the last temperature decoded in each cell by any stage. It is purely an accelerator: it seeds Newton, never changes the converged root or the tolerance, and a stale or absent hint costs only the ordinary safeguarded solve. Measured on an `h1600`, `ns = 1000` column, a seeded call converges in 1.13 iterations with zero bisections against 8.68 iterations and 4.05 bisections unseeded.
- **State-bound decode caching.** Each immutable conserved state is decoded once into a two-ghost-layer chromosphere::MixtureField, and the CFL, source and reconstruction stages share that field. Predicted and post-source states get separate cache generations. The conduction workspace and the MUSCL arrays are owned by the Grid and reused between steps, so there is no function-static state and no repeated hot-loop allocation.

## Timestep control

`chromosphere::mixture_timestep` returns the CFL-limited step from `|u| + c_s`. `CHROMO_CFL=0.50` is the validated production value for the reduced release configuration; the conservative hard-coded default remains 0.25 and stays the comparison reference. Rerun the sweep before trusting 0.50 on different hardware or a materially different model shape. Evidence: `docs/coarse_model_column_physical_conduction_recap.md`.

Two run-control variables matter for long runs. `CHROMO_T_END` sets an absolute stop time in physical seconds and disables the legacy `time_mult`-derived step cap, so a long run cannot be silently truncated. `CHROMO_FRAME_DT` sets the snapshot cadence in physical seconds rather than steps, and is required whenever output is enabled on a long run — without it the wall time is dominated by ASCII formatting. Every run prints `termination=end_time|step_cap|other`; always check that a long run ended with `end_time`. See @ref configuration.

## The two-fluid scheme

For contrast: the historical solver uses TVD-MUSCL with the symmetric minmod limiter, the Rusanov / local Lax-Friedrichs flux, and a semi-implicit operator-split driver with the stage sequence listed in `src/two_fluid/integrators.cpp`. Its implicit branch is currently disabled at the call site (see the known issue in @ref physics_model). Explicit RK4 is also available for the explicit-only path. Those are two-fluid entry points; the release timestep is `chromosphere::mixture_advance`.

## Rejected, and not to be reintroduced without new evidence

- The full SWMF `ModTransitionRegion` momentum pressure/gravity split. Measurably no better than the predictor source fix alone on this constant-area column, and not generally conservative for variable-area flux-tube geometry.
- A gravity-aware characteristic (non-reflecting) lower boundary. It changed the hydrostatic residual not at all and the evaporation velocity by at most 0.4 percent. The 1600 km face is not an important acoustic reflector for this model.

## See also

@ref validation for what is verified and what is not, @ref release_solver for the API, `docs/model_column_release_numerics_recap.md` for the full decision record.
