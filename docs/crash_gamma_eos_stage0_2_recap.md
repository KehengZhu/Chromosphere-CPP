# CRASH gamma-only EOS: Stage 0-2 implementation recap

Date: 2026-07-19

Scope implemented: the approved offline diagnostic (Stage 0), stripped
classical pure-H CRASH table and whole-grid Saha validation (Stage 1), the
standalone double-precision C++ Saha closure, and the Stage-2 Gamma1 table data
structure/parser/interpolator with tests.

## Stage 0: diagnostics and gates

`util/plot_eos_gamma.py` now fails closed before plotting. Over the entire raw
table it verifies:

- complete, finite, density-major rectangular axes;
- CRASH `zAv` against the analytic classical-Saha fraction, including `log(x)`
  where ionization is tiny;
- CRASH internal energy and energy gamma against the analytic caloric closure;
- `C_V^eff=(partial e_int/partial T)_rho`, evaluated by a central difference in
  `ln(T)`, against CRASH `Cv/(n_H k_B)`;
- analytic classical-Saha Gamma1 against CRASH `GammaS` (column-semantics check);
- strictly increasing `e_int(T)` on each density column, positive `C_V`, and
  positive `Gamma1` everywhere;
- exact raw-to-runtime extraction of `log10(T), log10(n_H), Gamma1`.

The generated C7 diagnostic walks the complete C7 profile from
`scenarios/model_c7.cpp`. It plots C7 vs Saha ionization, energy gamma, Gamma1,
and both heat capacities. The chromosphere/TR is expanded and the extended
corona compressed on a split height axis.

Generated figures (gitignored, cataloged in `visualize_commands.md`):

- `visualization/eos_gamma/gamma_hydrogen.png`
- `visualization/eos_gamma/gamma_hydrogen_2d.png`
- `visualization/eos_gamma/c7_eos_consistency.png`

## Stage 1: stripped CRASH table

`util/eos/tabulate_gamma.f90` now generates pure hydrogen only, on a
`57 x 501` grid:

- `T = 3.2e3 .. 1e8 K`;
- `n_H = 1e12 .. 1e26 m^-3`.

The lower temperature is just above CRASH's built-in forced-neutral threshold
(`0.02 chi_H/k_B`, about 3157 K) and below the coolest expected C7 state. The
upper and density bounds include margins around the expected runtime domain.

Configuration is explicit and cannot inherit future CRASH defaults silently:

- excitation off;
- Fermi gas off;
- Coulomb correction off;
- Boltzmann statistical-weight floor disabled for the requested classical gas;
- neutral-H ground degeneracy applied explicitly, cancelling CRASH's electron
  spin factor and matching the textbook coefficient-one Saha function;
- hydrogen ionization energy and the exact SI Boltzmann constant aligned with
  Chromosphere2026;
- population tolerance tightened and output precision increased to 16 digits;
- every CRASH EOS error aborts table generation.

Outputs:

- `outputs/eos_gamma/gamma_hydrogen.dat`: ignored full 11-column validation table;
- `data/eos/gamma1_hydrogen_v1.dat`: versioned stripped three-column runtime table;
- `data/eos/gamma1_hydrogen_v1.dat.sha256`: reviewed SHA-256 sidecar.

The runtime table is available in a clean clone. Its v1 header records the
configuration, statistical-weight convention, constants, dimensions,
generator, generation date, and source revision. `data/eos/README.md` documents
the regeneration and review policy.

`util/eos/validate_saha.cpp` reads all 28,557 raw rows and calls the actual C++
Saha function. It also loads the production table through the actual Stage-2
parser and verifies every Gamma1 grid node against the raw CRASH column.

## Stage 2: C++ Saha and Gamma1 loader

`eos.hpp` / `src/eos.cpp` provide:

- `saha_ionization_fraction(rho,T)` and an `n_H` overload, in double precision;
- one `eos_constants` namespace for `k_b`, `m_e`, `m_h`, `h`, and `chi_h`,
  also used to initialize the corresponding float `Grid` defaults;
- a log-domain, overflow/underflow-safe positive Saha root;
- strict positive/finite input validation;
- `EosGammaTable::load`, restricted by an exact machine-readable v1 schema to
  the stripped classical pure-H three-column format;
- rejection of missing, duplicate, unknown, malformed, or unsupported metadata
  and declared dimensions that disagree with the data;
- strict finite/rectangular/ascending/positive-Gamma1 validation;
- double-precision bilinear interpolation in `(log T, log n_H)`;
- a hard out-of-bounds error by default, with an explicit debug-only clamp;
- an empty/default table value in `Grid`, so the fixed-gamma solver remains
  untouched.

The loader is deliberately not wired to a `GAMMA_TABLE` scenario environment
hook yet. Enabling the table in hydro before the Stage-3 caloric inversion and
Stage-5.5 equilibrium-manifold reconstruction would create an invalid partial
EOS mode. Runtime activation, mode-compatibility checks, and the restructured
`model_column_ic` lifecycle therefore remain part of the gated Stage-3+ solver
work.

## Verification results

The revised generated table has SHA-256
`758dfe315cc907097a15368c1cfe242d279251c1f9c2f72ac1df7544d1635908` and passed:

| Check | Result |
|---|---:|
| C++ Saha rows checked | 28,557 |
| max absolute `zAv` difference | `5.2564e-7` |
| max transition relative difference | `1.73614e-3` |
| max tiny-ionization `abs(delta log x)` | `1.73464e-3` |
| runtime-loader Gamma1 node max difference | `1.39888e-14` |
| analytic-energy max relative difference | `1.140e-6` |
| energy-gamma max absolute difference | `2.100e-5` |
| central-difference `C_V` max relative difference | `8.255e-4` |
| analytic-vs-CRASH Gamma1 max absolute difference | `5.177e-4` |
| minimum `C_V/(n_H k_B)` | `1.5` |
| minimum Gamma1 | `1.05872` |

The small residual CRASH/Saha differences are bounded by CRASH's internal
ten-iteration population solve. The checks reject the materially wrong cases
encountered during implementation (factor-two statistical weight, legacy
Boltzmann floor, low-temperature forced-zero region, or mismatched constants).

Build and tests:

```bash
bash util/eos/build_and_run.sh
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

- `chromo_tests`: 6,393 checks passed; shared-constant, metadata rejection,
  malformed-dimension/axis/value, exact-boundary, density-OOB, and
  two-dimensional debug-clamp coverage added;
- CTest: 3/3 passed (`chromo_tests`, production-table validation, and checksum
  validation);
- Python syntax, shell syntax, and `git diff --check`: passed.

No Stage-3+ hydro, flux, RHS, source, boundary, HSE, or integrator behavior was
changed in this implementation.
