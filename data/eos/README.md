# Versioned EOS runtime data

`gamma1_hydrogen_v1.dat` is the Stage-2 production input for the classical,
equilibrium, pure-hydrogen CRASH equation of state. It contains only
`log10(T [K])`, `log10(n_H [m^-3])`, and the adiabatic sound-speed index
`Gamma1`. Its strict machine-readable header records the physics switches,
Saha convention, constants, dimensions, generator, and provenance.

The full 11-column validation table remains generated and ignored at
`outputs/eos_gamma/gamma_hydrogen.dat`. Regenerate both files with:

```bash
bash util/eos/build_and_run.sh
```

The tracked table is accepted only when its metadata matches the v1 loader
contract. Its committed SHA-256 is stored in
`gamma1_hydrogen_v1.dat.sha256`; the build-and-run workflow and CTest verify
that checksum through CMake's built-in cross-platform SHA-256 implementation
(`cmake/check_eos_checksum.cmake`), without an external checksum program. If an
intentional generator change alters the bytes, review
the physics and validation results before updating the sidecar.

This table is loaded and tested in Stage 2 but is not connected to hydro state
conversion, fluxes, sources, or boundaries until Stage 3.
