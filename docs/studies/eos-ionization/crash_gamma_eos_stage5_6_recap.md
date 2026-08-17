# CRASH gamma-only EOS: Stage 5–6 implementation recap

## Outcome

Stages 5, 5.5, 5.6, and 6 are implemented for the isolated gamma-table hydro
operator and projected Euler integrator. The fixed-gamma path remains the
default and continues through its original code.

Production scenario activation is deliberately still fail-closed. Stages 7–8
must first convert active energy sources, Saha-HSE initial states, and boundary
ghost packing to the same EOS. This avoids the unsafe intermediate state where
the new hydro decoder would consume legacy fixed-gamma IC or ghost energies.

## Stage 5: projected Euler split

`advance_Euler_state()` now has a gamma-table branch which:

1. validates the single-fluid compatibility contract before evaluating the RHS;
2. advances the equilibrium-mixture explicit hydro operator;
3. applies one conservative `project_equilibrium_single_fluid()` to the
   predictor; and
4. returns without entering the legacy drag and fixed-capacity temperature
   stages.

The branch conserves the Stage-4 total mass, total momentum, and `E_I+E_N`
invariants. Relative carrier drift generated during the explicit split is
thermalized only by the projection, so drag heating is not counted twice.

The initial implementation rejects gamma-table mode when conduction, TRAC,
legacy floors, beam/coronal heating, or radiative cooling are active. Those
updates still contain fixed-γ pressure increments and belong to Stage 7.
`advance_Euler_explicit_state()` and `advance_RK4()` also reject the mode rather
than evolving an off-manifold intermediate state.

## Stage 5.5: equilibrium-mixture MUSCL reconstruction

The gamma-table RHS is a separate early branch in `rhs_explicit_state()`.
It reconstructs only:

```text
{ log(rho_total), velocity, log(T) }
```

The original seven primitive variables are never passed to the gamma branch.
Both spatial reconstruction passes use double-precision mixture vectors,
double-precision minmod or MC3 limiter algebra, and double-precision
non-uniform-mesh weights on those three mixture slots. Conversion to the
project's packed float `Vec` occurs only when conserved face rows, physical
fluxes, and spectral radii are written back to the solver boundary.

Every reconstructed face is then built algebraically with
`equilibrium_mixture_face_state()`:

- `x_eq` comes from the unclamped Saha closure;
- `x_row` is the trace-floor-clamped carrier split;
- both velocities equal the reconstructed mixture velocity;
- pressures and chemical energy use `x_eq`;
- carrier masses, momenta, kinetic energy, and gravitational energy use
  `x_row`; and
- `E_E = 3 p_e / 2`.

The internal MUSCL predictor follows the required path:

```text
on-manifold face fluxes
  -> conservative half-step state
  -> decode_equilibrium_mixture(state, phi_state)
  -> second mixture reconstruction
  -> on-manifold corrected faces
```

It never calls legacy `cons2prim()`. Density and temperature remain positive by
log-space reconstruction; an invalid conserved predictor fails through the
strict mixture decoder rather than being independently clipped by species.

Shifted states are decoded with the gravitational potential actually carried by
the shifted cell or ghost. At each interface, both left and right states are
packed with the same face value: `phi_g_iph` at the right face and `phi_g_imh`
at the left face. This removes a gravitational-energy jump from the Rusanov
dissipation.

`MixtureFaceState` carries the already known `(rho,v,T,x_eq,x_row,p_*,Gamma1)`,
packed conserved rows, total pressure, and sound speed. Physical flux and
spectral-radius evaluation therefore do not repeat the caloric inversion.

## Stage 5.6: trace-carrier semantics

The face builder preserves the Stage-4 split everywhere:

- `x_eq` drives `p_i`, `p_n`, `p_e`, chemical energy, `n_e`, `n_HI`, caloric
  energy, and CRASH `Gamma1`;
- `x_row` drives only the stored carrier densities, momenta, and their
  kinetic/gravitational energy terms.

Tests vary the storage floor and confirm that face thermodynamics and the common
sound speed are unchanged while only the trace carrier row changes.

## Stage 6: common mixture wave speed

Gamma-table fluxes now use one wave speed for every row:

```text
c_s^2 = Gamma1 * (p_i + p_n) / (rho_i + rho_n)
a     = |v| + c_s
```

The scalar face API, packed `cal_spectral_radius_state()` branch, and Rusanov
face bundles all use this same definition. The tests also verify the total flux
identities

```text
F_rho_i + F_rho_n = rho v
F_mom_i + F_mom_n = rho v^2 + p_total
F_E_i + F_E_n     = (E_I + E_N + p_total) v
```

The operator-splitting convergence test evaluates three successively halved
timesteps. The pre-projection relative kinetic energy is positive and decreases
quadratically: both measured ratios are within the test tolerance of `0.25`.
The scalar projection's internal-energy increment matches that relative kinetic
energy, confirming that projection heating vanishes as `dt -> 0`.

## Runtime lifecycle and guards

`chromo_main` now recognizes `GAMMA_TABLE=path`, loads the strict production
table, and rejects an explicitly supplied `ISO_GAMMA` as mutually exclusive.
After successful loading it currently stops with an explicit Stage-7/8 message
before scenario IC construction. This is intentional: the controlled solver
tests exercise the completed Stage-5/6 operator, while production scenarios
cannot accidentally activate it with old IC or boundary packing.

## Review follow-up: table-edge precision

The review identified that the first implementation stored `log(rho)`, `v`, and
`log(T)` in the project's float `Vec`. A valid table-edge state could therefore
move outside the strict double-precision EOS domain after
`double -> float(log) -> double(exp)`, even without a physical or limiter
overshoot.

The reconstruction container is now `arma::Col<double>`. Decode output, logs,
limiter ratios, MC3/minmod evaluation, metric weighting, face extrapolation,
and predictor time averaging all remain double. The production hard-error
policy and its double-precision table tolerance were not relaxed.

A production-table regression covers exact mathematical `T_min`, `T_max`,
`n_H,min`, and `n_H,max` through the double face API. For the packed solver path,
it constructs the closest float-representable state on the legal side of each
endpoint, then runs both `rhs_explicit_state()` and `advance_Euler_state()`.
This distinguishes reconstruction drift from the unavoidable quantization of
two float carrier-energy or density rows.

## Verification

The following checks passed on 2026-07-19:

- build completed for `chromosphere`, `chromo_main`, both EOS validators, and
  `chromo_tests`;
- unit suite: **6,556 checks passed, 0 failed**;
- CTest: **3/3 passed**;
- full CRASH validation: **28,557 rows**;
- C++ Saha vs CRASH transition relative error: `0.00173614`;
- runtime Gamma1 node maximum absolute error: `1.39888e-14`;
- production-table checksum:
  `758dfe315cc907097a15368c1cfe242d279251c1f9c2f72ac1df7544d1635908`;
- production validator passed the Saha, analytic closure, heat-capacity, and
  positive-Gamma1 gates;
- Python syntax check, shell syntax check, and `git diff --check` passed.

No commit was created.

## Remaining work

Stage 5–6 does not claim production scenario readiness. The next activation
blockers are the planned Stage-7 total-energy source/conduction conversions and
Stage-8 Saha-HSE IC plus EOS-consistent ghost/boundary packing. Once those are
implemented, the fail-closed `chromo_main` gate can be replaced by the finalized
mode lifecycle and per-step EOS-domain diagnostics.
