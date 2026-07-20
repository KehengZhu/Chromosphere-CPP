# Effective-γ from the SWMF/CRASH equation of state

This tool extracts two effective thermodynamic indices of a partially ionized
plasma from the SWMF/CRASH statistical-sum EOS. They can be compared with the
constant `ISO_GAMMA=1.05` currently used in the chromosphere model, but a table
lookup of γ alone is **not** a thermodynamically complete replacement for the
EOS.

## What the CRASH EOS is

`SWMF/util/CRASH/src/` is the equation-of-state package of the CRASH
radiation-hydrodynamics code (originally for xenon/beryllium ICF targets, but
the underlying physics is general). For a given element/mixture at electron
temperature `Te` and heavy-particle number density `Na` it solves **Saha
ionization equilibrium** and returns thermodynamics that include:

- ideal translation of ions and electrons,
- **ionization energy** (the dominant effect — this is what softens γ),
- bound-state **excitation** (H, He tabulated),
- electron **Fermi degeneracy**,
- **Coulomb / Debye** correction.

This is a local-thermodynamic-equilibrium (LTE) material EOS. Radiation energy,
radiative loss, thermal conduction, photoionization, and finite-rate
ionization/recombination are *not* folded into γ. Those are separate radiation,
transport, or kinetic source terms.

## The two gammas (`get_gamma` in `ModStatSum.f90`)

- **energy γ** `= 1 + P/e` — reproduces the EOS internal energy at one state via
  `e = P/(γ−1)`. This is the closest match to the model's present algebraic
  energy closure. `GammaOut` in the code.
- **adiabatic Γ₁** `= (∂lnP/∂lnρ)_S` — governs the sound speed / wave dynamics.
  `GammaSOut` in the code.

Both dive from 5/3 toward ~1.1 through the ionization zone and recover to 5/3
once neutral or fully ionized; γ is a function of `(T, n, element)`, with the
ionization degree `Z = Z(T, n)` derived internally.

These two quantities generally differ in the ionization zone. Substituting the
energy γ into every place where the current code uses one constant γ would make
the stored energy more EOS-like, but would not automatically give the correct
sound speed, Jacobian, or conservative flux. A full implementation should
tabulate at least `P(ρ,e)`, temperature, and thermodynamic derivatives (including
Γ₁), rather than treating γ as the only state variable.

## Interaction with the existing ionization stage

`src/integrators.cpp::apply_ionization_stage` already evolves a finite-rate
hydrogen network and explicitly accounts for the 13.6 eV ionization potential.
The CRASH LTE EOS also includes ionization energy in `e(T,ρ)`. Therefore:

- using the CRASH EOS is appropriate for an LTE-equilibrium closure in which
  the EOS determines the ionization state; or
- keeping the present non-equilibrium ionization stage requires an EOS whose
  independent composition/ionization fraction is supplied explicitly.

Using the current CRASH equilibrium γ table together with the existing
ionization-energy source unchanged would double-count latent ionization energy.

## Files

- `tabulate_gamma.f90` — standalone driver: sweeps a `(T, Na)` grid for pure
  hydrogen and an H0.9/He0.1 mix, writes `outputs/eos_gamma/*.dat`.
- `build_and_run.sh` — compiles the driver against the prebuilt `libCRASH.a`
  (+ `libSHARE`, `libTIMING`) and runs it.

## How to run

```bash
# one-time: build the CRASH EOS library (gfortran + Open MPI required)
cd SWMF/util/CRASH/src && make LIB && cd -

# compile the tabulator and generate the tables
bash util/eos/build_and_run.sh

# plot
python util/plot_eos_gamma.py     # -> visualization/eos_gamma/*.png
```

## Note on `make test_eosgodunov`

The BATSRUS recipe `make test_eosgodunov` (GM/BATSRUS) cannot be built from this
SWMF checkout: it needs the user module `ModUserEos.f90`, which lives in the
optional/restricted `srcUserExtra` repository that is not present here. The
standalone `make GODUNOV` in `util/CRASH/src` has the same problem — its exact
Riemann reference routines (`pu_star`, `sample`) are absent from the checkout.
Neither is needed for tabulating γ: the shock tube is only a solver
demonstration, whereas the γ values come straight from the CRASH EOS library,
which builds and runs cleanly (`make LIB`).

## Table columns

```
1 T[K]  2 Na[m^-3]  3 Rho[kg/m^3]  4 Zbar  5 P[Pa]  6 Edens[J/m^3]
7 Gamma  8 GammaS  9 Gammae  10 GammaSe  11 Cv[k_B/atom]
```
