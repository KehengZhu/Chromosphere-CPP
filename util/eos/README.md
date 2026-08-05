# Classical pure-hydrogen Saha / Γ₁ tables from SWMF/CRASH

This tool builds the validation and runtime inputs for stages 0–9 of
`docs/crash_gamma_eos_plan.md`. CRASH supplies the equilibrium sound-speed
index `Γ₁(T,n_H)`; Chromosphere2026 computes the pure-H Saha ionization and
caloric closure analytically.

## EOS configuration

CRASH's statistical-sum EOS can include excitation, electron degeneracy, and
Coulomb corrections. For this project's first implementation all three are
explicitly disabled:

- `UseExcitation = .false.`
- `UseFermiGas = .false.`
- `UseCoulombCorrection = .false.`

Hydrogen ground-state statistical weights are explicitly enabled. The neutral
ground degeneracy 2 cancels CRASH's explicit electron-spin factor 2, matching
the textbook Saha function already used in `physics.hpp`:
`(2πm_e k_B T/h²)^(3/2) exp(-χ_H/k_B T)`. The output is pure hydrogen; the
former H/He comparison is not part of this closure.

CRASH returns two useful indices:

- energy `γ_E = 1 + P/e`, retained only to validate the analytic caloric
  closure;
- adiabatic `Γ₁ = (∂lnP/∂lnρ)_S`, the runtime sound-speed table.

`zAv`, pressure, energy, and `Cv` are also retained in the raw file for strict
whole-grid checks. Radiation, conduction, cooling, photoionization, and
finite-rate ionization are not part of this LTE table. Gamma-table equilibrium
mode must keep the existing finite-rate ionization stage off to avoid counting
the hydrogen ionization energy twice.

## Files and range

- `tabulate_gamma.f90` sweeps `T=3.2e3..1e8 K` and
  `n_H=1e12..1e26 m^-3`. The 3200 K edge is just above CRASH's built-in
  `0.02 chi_H` forced-neutral cutoff and below the coolest expected C7 state.
- `build_and_run.sh` compiles the driver against `libCRASH.a`, `libSHARE`, and
  `libTIMING`, generates the production table plus a disposable 2x-refined
  Stage-9 grid, then runs strict validation. It invokes
  `eos_validate_saha` so every CRASH `zAv` row is compared with the actual C++
  log-domain Saha function and every raw `GammaS` row is compared with an
  independent analytic C++ `Gamma1` derived from Saha pressure and caloric
  derivatives. It verifies the tracked table checksum and then runs the broader
  Python thermodynamic checks. `eos_validate_refined_table`
  exhaustively checks all 28,000 production-cell midpoint nodes from the
  refined direct-CRASH grid, including coarse-grid Gamma1 interpolation,
  analytic gamma_E and Cv, and temperature inversion.

The generated files are:

- `outputs/eos_gamma/gamma_hydrogen.dat`: ignored 11-column validation table;
- `data/eos/gamma1_hydrogen_v1.dat`: tracked three-column production table;
- `data/eos/gamma1_hydrogen_v1.dat.sha256`: reviewed production checksum.
- `outputs/eos_gamma/gamma_hydrogen_refined.dat`: ignored 1001x113 direct-CRASH
  Stage-9 reference (every production interval split in two);
- `outputs/eos_gamma/gamma1_hydrogen_refined.dat`: ignored refined runtime table.

## How to run

```bash
# one-time: build the CRASH EOS library (gfortran + Open MPI required)
cd SWMF/util/CRASH/src && make LIB && cd -

# compile, generate, and validate
bash util/eos/build_and_run.sh

# validate again and generate Stage-0 figures
python util/plot_eos_gamma.py
```

The raw validation columns are:

```
1 T[K]  2 n_H[m^-3]  3 Rho[kg/m^3]  4 zAv  5 P[Pa]  6 Edens[J/m^3]
7 GammaE  8 Gamma1  9 Gammae  10 GammaSe  11 Cv[k_B/heavy-particle]
```

The runtime table contains exactly:

```
1 log10(T[K])  2 log10(n_H[m^-3])  3 Gamma1
```

Its v1 key/value header is part of the file contract. The C++ loader rejects
missing, duplicate, or unknown keys, unsupported physics/configuration values,
wrong constants, invalid declared dimensions, malformed rows, nonascending
axes, and nonpositive/nonfinite Gamma1. See `data/eos/README.md` for provenance
and update policy.

The BATSRUS `test_eosgodunov` target is not required. That demonstration needs
the optional `srcUserExtra/ModUserEos.f90`, which is absent from this checkout;
the standalone tabulator calls the already-built CRASH EOS library directly.

## Stage-3 consumer

`eos.hpp` and `src/eos.cpp` use this table in the narrow read-only mixture
closure. The analytic Saha model supplies `x_eq`, total internal energy, and
energy gamma; this file supplies only Γ₁. The self-contained temperature
inversion uses the table temperature bounds as a bracket and safeguarded
Newton/bisection, with any supplied previous temperature treated only as an
optional initial guess.

`decode_equilibrium_mixture` subtracts the two stored rows' original kinetic
energies and the gravitational potential carried by that state. It does not
project or mutate the rows. Conservative center-of-mass projection is Stage 4;
Stages 5--8 connect the closure to reconstruction, fluxes, source/conduction
splitting, and Model C7 initial/boundary states.

Stage 4 adds `project_equilibrium_single_fluid`. Unlike the read-only decoder,
the projection subtracts center-of-mass kinetic energy, so relative drift is
thermalized while total mass, momentum, and `E_I+E_N` remain conserved. It
rebuilds both carrier rows with the approved pressure/chemical-energy mapping
and keeps `E_E=3p_e/2`. The packed-float wrapper is used by the restricted
gamma-table Euler path and remains separate from read-only decode.

The production density range is enforced before caloric root-finding, even
though the analytic energy formula itself can be evaluated outside the table.
This prevents a future projection from completing where no production Γ₁
exists. An explicit debug flag bypasses the density gate and pairs with the
table's debug-only Γ₁ edge clamp; production code must leave it disabled.

Checksum verification uses `cmake/check_eos_checksum.cmake` and CMake's built-in
`file(SHA256)`, so configuring or testing the fixed-γ solver does not require an
external `shasum` or `sha256sum` executable.
