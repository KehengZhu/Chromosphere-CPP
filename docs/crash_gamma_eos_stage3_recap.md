# CRASH gamma-only EOS: Stage 3 implementation recap

Date: 2026-07-19

Scope implemented: the narrow analytic pure-hydrogen mixture closure,
self-contained caloric temperature inversion, and read-only conserved-row
decoder specified by Stage 3 of `docs/crash_gamma_eos_plan.md`.

## Public closure API

`eos.hpp` / `src/eos.cpp` now provide:

- `GammaState {x_eq, gamma_energy, gamma_sound}`;
- `MixtureThermo {rho, T, x_eq, x_row, n_H, n_e, n_HI, p_i,
  p_n, p_e, gamma1, internal_energy}`;
- `equilibrium_internal_energy(rho,T)`;
- `gamma_state(table,rho,T)`;
- `temperature_from_rho_eint(table,rho,e_int,T_guess)`;
- `decode_equilibrium_mixture(table,rho_i,rho_n,mom_i,mom_n,E_i,E_n,phi)`.

The caloric closure is analytic and uses the same double-precision Saha
solution and shared constants as Stages 1–2:

```text
n_H   = rho / m_H
p     = (1 + x_eq) n_H k_B T
e_int = (3/2) p + x_eq n_H chi_H
gamma_E = 1 + p/e_int
```

CRASH remains authoritative only for `Gamma1(T,n_H)`, obtained through the
strict tracked-table interpolator.

## Safeguarded temperature inversion

The inverse solves `e_int(rho,T)=e_target` using the production table's
temperature interval as a hard bracket. It:

- validates positive, finite density and target energy;
- verifies the target lies between the analytic endpoint energies;
- uses analytic `C_V = partial e_int / partial T` for Newton steps;
- falls back to bisection whenever Newton leaves the bracket;
- converges without a cached or previous temperature;
- treats a finite in-bracket `T_guess` only as an accelerator;
- rejects caloric targets below or above the supported temperature range.
- rejects density outside the production Γ₁ table before starting the root-find;
  an explicit debug-only bypass is available for controlled diagnostics.

The residual tolerance is scaled by the requested energy, not the energy at
the table's hot endpoint. This matters for accurate cool-state inversions over
the table's more than four-decade temperature range.

## Read-only mixture decode

The decoder computes

```text
e_int = E_i + E_n
        - mom_i^2/(2 rho_i) - mom_n^2/(2 rho_n)
        - (rho_i + rho_n) phi_of_this_state
```

and then inverts for temperature. This deliberately subtracts the two original
kinetic energies. It does not use center-of-mass kinetic energy, thermalize
relative drift, modify composition, or rewrite the input rows; those operations
belong to the Stage-4 conservative projection.

`x_row = rho_i/(rho_i+rho_n)` records the numerical storage split, while the
unclamped Saha value `x_eq` drives every returned physical quantity:

```text
n_e  = x_eq n_H
n_HI = (1-x_eq) n_H
p_e  = x_eq n_H k_B T
p_i  = 2 p_e
p_n  = (1-x_eq) n_H k_B T
```

The supplied gravitational potential is the potential actually carried by the
decoded state; no fixed-gamma well-balanced pressure correction is applied.

## Verification

The Stage-3 unit test checks:

- energy-to-temperature round trips at both bracket boundaries and through the
  ionization regime, with no guess and with a deliberately poor guess;
- analytic `x_eq` and energy gamma plus interpolated CRASH Gamma1;
- a non-equilibrium two-row state with unequal velocities and `x_row != x_eq`;
- recovery of `T`, `e_int`, `n_e`, `n_HI`, `p_i`, `p_n`, and `p_e`;
- rejection of below/above-bracket energy, Gamma1 density OOB, and invalid rows.

The production-table CTest additionally performs 24 caloric round trips over
`n_H={1e12,1e16,1e20,1e26} m^-3` and
`T={3.2e3,5e3,1e4,1e5,1e7,1e8} K`, including all table edges.

The post-review hardening tests also cover:

- density-OOB rejection directly at `temperature_from_rho_eint`, plus the
  explicit debug bypass;
- `NaN`, infinite, below-bracket, and above-bracket guesses producing the same
  result as a no-cache inversion;
- energies `e_min(1+1e-10)` and `e_max(1-1e-10)` returning interior
  temperatures rather than being mistaken for exact endpoints;
- `EosGammaTable::contains` at exact boundaries and outside both axes.

The table exposes `contains(T,n_H)` and
`require_n_h_in_bounds(n_H,debug_clamp)`. Density errors report `rho`, `n_H`,
and the permitted table range, allowing future Stage-4 callers to add cell/time
context at their boundary.

## Build portability hardening

The checksum CTest and `util/eos/build_and_run.sh` now use
`cmake/check_eos_checksum.cmake`, which computes the hash with CMake's built-in
`file(SHA256)`. Project configuration no longer fails on Linux/HPC systems that
lack the external macOS-style `shasum` utility.

Final checks:

- `chromo_tests`: 6,434 checks passed, 0 failed;
- CTest: 3/3 passed;
- the full 28,557-row CRASH validation remains unchanged;
- the tracked runtime-table SHA-256 remains
  `758dfe315cc907097a15368c1cfe242d279251c1f9c2f72ac1df7544d1635908`.

## Deliberately not implemented

Stage 3 does not activate `GAMMA_TABLE`, branch `cons2prim`/`prim2cons`, project
carrier rows, change fluxes or spectral radii, reconstruct equilibrium face
states, alter sources/conduction/boundaries, or restructure the C7 HSE. Those
are Stage 4–8 tasks and require the plan's conservative projection and
equilibrium-manifold reconstruction rather than a partial Gamma1 substitution.
