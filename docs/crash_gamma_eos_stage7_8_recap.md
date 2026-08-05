# CRASH gamma-only EOS: Stage 7-8 implementation recap

Revised: 2026-07-20

## Result

Stages 7 and 8 are implemented for equilibrium gamma-table operation through
`advance_Euler_state` and the `model_column` scenario. Fixed-gamma execution
continues through the original source, conduction, IC, and boundary paths.

Gamma-table production activation remains deliberately narrow:

- supported scenario: `model_column` (`model_isentropic` compatibility alias);
- supported integrator: the semi-implicit Euler driver;
- required physics mode: single fluid, common temperature, equilibrium
  ionization (`ENABLE_TE=0`, explicit `no-ionization` argument);
- unsupported scenarios, explicit-only stepping, RK4, `ISO_GAMMA`, and
  two-fluid/non-equilibrium combinations fail before evolution.

The review follow-up tightened this contract before scenario initialization:
explicit `ENABLE_TE=1`, explicit `SINGLE_FLUID=0`, and a finite-rate
`ionization` command-line request now terminate instead of being silently
overwritten by `model_column_ic()`.

## Stage 7: variable-EOS source coupling

### Total-energy source updates

The gamma-mode Euler split now uses this order:

1. explicit hydro predictor;
2. conservative equilibrium projection;
3. TRAC cutoff refresh, when enabled;
4. nonlinear mixture conduction;
5. beam heating;
6. coronal heating;
7. radiative cooling.

Beam and coronal heating add `dt Q` to the mixture total-energy row sum, then
immediately call the equilibrium projection. Radiative cooling computes the
physical optically thick plus optically thin loss, applies TRAC only to the
thin term, removes total internal energy with a positivity-preserving backward
relaxation, bounds it by the EOS energy at the table minimum temperature, and
re-projects. No gamma-dependent pressure increment is used in this path.

The legacy flare density/pressure floor is not applied in gamma mode. Trace-row
conditioning remains owned by the equilibrium projection's `x_row` storage
floor, leaving physical `x_eq` unchanged.

### Float-safe EOS-boundary packing

All gamma energy operators now use one internal-energy target helper. It:

1. computes the requested internal energy in double precision;
2. rejects non-finite or physically out-of-table targets with cell, time,
   density, target-energy, valid-energy, and temperature-range context;
3. packs the conservative `E_N` remainder into the float state;
4. moves that remainder with `nextafter` only when a physically legal double
   target is rounded outside the strict interval by float storage;
   outside the strict EOS interval;
5. projects, repeats the legality check on the projection's final float rows,
   and performs a strict decoder verification.

Cooling explicitly owns its approved Tmin floor before calling the helper.
Heating and conduction are never silently truncated: a physical target beyond
Tmax/Tmin remains a production hard error. After the final remainder adjustment,
`E_E` is re-slaved to `3 p_e/2` from the strict decoded state.

### Physical composition for source and transport coefficients

Gamma source paths decode each cell through `decode_equilibrium_mixture()` and
use its physical `T`, `n_e`, `n_HI`, and total internal energy. Beam deposition,
radiative loss, TRAC detection/broadening, electron conductivity, and neutral
conductivity therefore do not reinterpret the trace carrier rows as physical
species densities.

The gamma-mode heating timestep guards likewise use decoded total equilibrium
internal energy and physical electron/neutral number densities.

### One nonlinear mixture conduction solve

The gamma conduction stage solves one common-temperature backward-Euler
residual,

```text
e_int(rho,T_new) - e_int(rho,T_old)
  - dt div[(kappa_e + kappa_n) grad(T_new)] = 0.
```

The quasi-Newton/Picard tridiagonal system uses the analytic
`equilibrium_heat_capacity(rho,T) = d e_int/dT` on its diagonal. Conductivities
and Saha composition are refreshed at every nonlinear iteration. TRAC broadens
only `kappa_e`; `kappa_n` is unchanged. Nonuniform meshes retain the existing
series-resistance face conductivity, numerical temperature diffusivity is
converted with the same effective heat capacity, and both Dirichlet/Neumann
inner and imposed-flux outer boundary choices are supported. The converged
internal-energy change is applied to total energy and projected once.

Convergence is declared only by the full nonlinear residual. A small update
forces one more full residual evaluation at the updated temperature: it is
accepted if that residual passes and otherwise reports stagnation. The final
iteration is also evaluation-only, so its preceding update cannot be rejected
without checking the new residual. Exhausting the iteration limit reports
non-convergence; no path returns an edge-damped state as a solution. A public
residual diagnostic recomputes the complete
backward-Euler residual from the before/after packed states for regression use.

## Stage 8: Saha-consistent C7 HSE and boundaries

### Initialization lifecycle and HSE

`model_column_ic()` now establishes gamma mode before constructing the
atmosphere. It loads `GAMMA_TABLE` when needed, rejects `ISO_GAMMA`, parses and
validates the single-fluid/ionization toggles, then builds the HSE profile.

The C7 temperature profile is retained. The base anchor is the C7 total
hydrogen density,

```text
n_H,base = n_i,C7 + n_HI,C7,
x_base   = x_Saha(n_H,base,T_base),
p_base   = (1+x_base)n_H,base k_B T_base.
```

Upward integration uses the trapezoidal logarithmic-pressure equation. Each
cell performs an inner pressure-to-density Saha solve so the local scale height
and composition are mutually consistent. The iteration records convergence and
throws a cell/height/temperature/pressure diagnostic if its 60-iteration limit
is reached. The resulting IC is packed with
`equilibrium_mixture_face_state()` at the cell-centered gravitational
potential. The optional coronal temperature boost uses the same mapper rather
than rewriting fixed-gamma energy rows.

### Ghost states and sound-speed cap

All four column ghosts are rebuilt from total density, temperature, velocity,
and their own boundary-face gravitational potential through the equilibrium
row mapper. Outer ghost density is obtained by solving

```text
p = (1+x_Saha(rho,T)) rho k_B T / m_H
```

within the production table density domain. Inner ghosts extend the base
isothermal HSE reservoir using total mixture pressure and the same Saha density
solve. Interior boundary diagnostics decode the equilibrium mixture, and the
outer velocity cap uses `sqrt(Gamma1 p/rho)`.

Ghost construction precedes freezing `grid.eq_state`, so equilibrium-reference
well balancing sees the completed gamma-consistent IC and boundary lifecycle.

### Production diagnostics

Gamma runs retain the seven-row conserved snapshot for restart/raw use and mark
it with `EOS_MODE=gamma_table`. A companion `<output>.gamma_diag` frame contains
decoded `rho_total`, `v_cm`, `T`, `x_eq`, `n_e`, `n_HI`, `p_total`, `Gamma1`,
and cell-centered solver-effective conductivity, including TRAC broadening and
the configured numerical diffusivity. Its header records the actual table path,
runtime SHA-256, and gravitational-potential convention. The canonical column
animation prefers a current sidecar and therefore never interprets `x_row` as
physical ionization.

## Files changed in Stages 7-8

- `eos.hpp`, `src/eos.cpp`: public analytic equilibrium heat-capacity API.
- `src/integrators.cpp`: decoded gamma thermodynamic fields, nonlinear mixture
  conduction, TRAC cutoff, total-energy heating/cooling, and Euler split activation.
- `scenarios/model_column.cpp`: Saha-HSE integration, EOS IC/boost/ghost packing,
  decoded boundary state, and mixture Mach cap.
- `chromo_main.cpp`: restricted production activation plus pre-IC and post-IC
  compatibility checks, portable runtime table SHA-256, and decoded gamma
  diagnostic sidecar output. The raw stream also declares `EOS_MODE=gamma_table`.
- `CMakeLists.txt`: process-level rejection tests for incompatible gamma requests.
- `tests/chromo_tests.cpp`: Stage 7 source/conduction and Stage 8 HSE/ghost/Euler tests.
- `util/animate_gentle_column.py`: automatic gamma-sidecar ingestion and
  physical-composition, TRAC-aware solver-effective conductivity diagnostics.
- `README.md`: activation contract and current regression coverage.

## Regression coverage added

The new tests check:

- gamma conduction is accepted while `ENABLE_TE=1`, `SINGLE_FLUID=0`, and
  finite-rate ionization process requests fail before IC;
- coronal heating changes total energy by the expected packed-float `dt sum(H)`
  while total cell mass is unchanged;
- radiative cooling decreases total energy and remains finite;
- strong cooling at four density decades reaches Tmin without crossing the
  strict caloric domain, while strong heating beyond Tmax throws instead of
  losing source energy; a non-finite heating target also throws;
- a genuinely zero-flux, two-ended conduction problem on a nonuniform,
  variable-area mesh conserves the physical weighted domain energy
  `sum((E_I+E_N) ds/B)`;
- the packed conduction result satisfies a separately recomputed nonlinear
  backward-Euler residual;
- a strong nonlinear `10^4`--`10^6 K` TRAC profile converges, remains in the
  table, conserves closed-domain weighted energy, and passes the packed residual;
- all Model C7 gamma cells decode on the equilibrium manifold;
- the integrated column satisfies the trapezoidal Saha-HSE residual;
- all four ghosts satisfy Saha composition and the total-energy identity using
  their own face potential;
- five `cal_dt_i()` gamma `model_column` steps keep maximum velocity and total
  momentum near zero and preserve the equilibrium-reference state.
- a zero-time production smoke run writes the decoded gamma sidecar and records
  the actual runtime table SHA-256 and potential convention.

## Validation performed

```text
cmake -S . -B build
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

Result: 7/7 CTests passed. The main C++ suite reports 67 cases and 6,708 checks,
all passing.

```text
bash util/eos/build_and_run.sh
```

Result: 28,557 CRASH rows validated; C++ Saha and runtime Gamma1 node checks
passed; whole-grid thermodynamic gates passed. The tracked runtime table SHA-256
is:

```text
758dfe315cc907097a15368c1cfe242d279251c1f9c2f72ac1df7544d1635908
```

OpenMP was not found by this local CMake configuration, so the validation used
the serial build. This is a configuration note, not a test failure.
