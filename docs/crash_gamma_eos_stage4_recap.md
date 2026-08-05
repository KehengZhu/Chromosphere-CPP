# CRASH gamma-only EOS: Stage 4 implementation recap

Date: 2026-07-19

Scope implemented: the conservative equilibrium projection and explicit
seven-row energy mapping specified by Stage 4 of
`docs/crash_gamma_eos_plan.md`.

## Projection APIs

The double-precision EOS layer now exposes:

```cpp
struct ProjectedMixture {
    double rho_i, rho_n, momentum_i, momentum_n;
    double energy_i, energy_n, energy_e;
    MixtureThermo thermo;
};

ProjectedMixture project_equilibrium_single_fluid(...);
```

`state.cpp` provides a packed-float overload:

```cpp
Vec project_equilibrium_single_fluid(const Grid&, const Vec&);
```

The packed overload applies the scalar projection independently to every cell
using that cell's centered gravitational potential. It is available for Stage
5 but is not called by an integrator in Stage 4.

## Conservative center-of-mass projection

For each input state the projection computes:

```text
rho = rho_i + rho_n
momentum = momentum_i + momentum_n
v_cm = momentum/rho
e_int,proj = E_I + E_N - 1/2 rho v_cm^2 - rho Phi
```

It then solves the Stage-3 caloric inversion for `T`. Because the projection
subtracts center-of-mass kinetic energy rather than the two original row
kinetic energies,

```text
e_int,proj = e_int,old + 1/2 mu (V-U)^2,
mu = rho_i rho_n/rho.
```

Thus the removed relative-drift kinetic energy is thermalized exactly once.
The returned rows share `V=U=v_cm`. The original total mass, momentum, and
`E_I+E_N` are conserved; the neutral energy is stored as the conservative
remainder after constructing the charged row from the analytic mapping.

## Composition and row mapping

The projection keeps the Stage-3 distinction:

```text
x_eq  = x_Saha(rho,T)                         physical composition
x_row = clamp(x_eq,f_min,1-f_min)             storage composition only
```

The default carrier floor is `f_min=1e-8` (`Grid::eos_trace_fraction_floor`).
Only `rho_i`, `rho_n`, and their common-velocity momentum/gravity terms use
`x_row`. All thermodynamics use unclamped `x_eq`:

```text
p_e = x_eq n_H k_B T
p_i = 2 p_e
p_n = (1-x_eq) n_H k_B T

E_I = 3/2 p_i + x_eq n_H chi_H + 1/2 rho_i v_cm^2 + rho_i Phi
E_N = 3/2 p_n                 + 1/2 rho_n v_cm^2 + rho_n Phi
E_E = 3/2 p_e
```

`E_E` is overwritten from the fixed monatomic electron relation and remains a
diagnostic row outside the conserved `E_I+E_N` total. It never uses CRASH
energy gamma.

## Validation

The Stage-4 unit test starts from a genuine unequal-velocity, off-equilibrium
carrier split and verifies:

- conservation of total mass, momentum, and `E_I+E_N`;
- common post-projection velocity equal to `v_cm`;
- the internal-energy increase equals `1/2 mu(V-U)^2`;
- inversion of that energy gives a hotter equilibrium state;
- charged, neutral, chemical, kinetic, gravity, and electron row mappings;
- decoding the projected rows recovers the projected thermodynamics;
- applying the projection twice is idempotent to inversion tolerance;
- invalid trace floors are rejected;
- floors `1e-6` and `1e-10` alter carrier mass only, while `x_eq`, `T`, `p_i`,
  `p_n`, and `n_e` are unchanged;
- the packed float wrapper conserves row sums to float tolerance, enforces a
  common velocity, and replaces a deliberately invalid old `E_E` with
  `3p_e/2`.

## Review follow-up and formal sign-off coverage

The conditional-review test blocker is closed. Tests now explicitly verify:

- caloric inversion rejects density above the table maximum as well as below
  its minimum;
- scalar projection rejects otherwise-valid states below and above the
  production density range before constructing output rows;
- a two-cell packed projection throws when cell 1 is density-OOB after cell 0
  has been processed into the function-local result, while the caller's input
  `Vec` remains exactly unchanged;
- the nearly fully ionized upper clamp gives `x_row=1-f_min`, changes only the
  neutral carrier rows between `f_min=1e-6` and `1e-10`, and leaves `x_eq`, `T`,
  `p_i`, `p_n`, `p_e`, `n_e`, `n_HI`, and `E_I+E_N` unchanged; the projected
  `p_n→0` state remains decodable;
- a nonzero relative drift with exactly zero total momentum projects to
  `M_I=M_N=0` and converts all original kinetic energy into internal energy.

These tests establish density-OOB failure before caller-visible mutation and
cover both lower and upper trace-species clamps. Stage 4 is therefore no longer
conditional on additional test coverage.

The real production-table CTest also runs a non-equilibrium conservative
projection and checks mass, momentum, energy, drift thermalization, and
electron-energy mapping.

Final verification:

- `chromo_tests`: 6,490 checks passed, 0 failed;
- CTest: 3/3 passed;
- the full 28,557-row CRASH/Saha/Gamma1 validation passed unchanged;
- production checksum remains
  `758dfe315cc907097a15368c1cfe242d279251c1f9c2f72ac1df7544d1635908`;
- Python syntax, shell syntax, and `git diff --check` passed.

During Stage 4, `Grid::init` was also aligned with the shared
`eos_constants` source, eliminating its remaining reintroduction of the older
electron-mass literal after default construction.

## Deliberately not activated

The projection is not yet called by `advance_Euler_state`; the old Stage B/C
split remains unchanged. `cons2prim`, `prim2cons`, fluxes, spectral radii, RHS,
sources, boundaries, C7 HSE, and scenario lifecycle are also unchanged.

Activation must wait for Stage 5 and the Stage-5.5 equilibrium-manifold MUSCL
path. Calling the projection at the end of a step before face reconstruction is
made manifold-consistent would still allow off-EOS face states to feed the
Rusanov flux and internal predictor.
