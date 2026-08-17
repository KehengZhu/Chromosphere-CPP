# Physical conductive-flux diagnostic recap

## Outcome

Fixed-gamma and gamma-table plots now use the same cell-centred diagnostic,

\[
q_{\rm phys}=-(\kappa_e+\kappa_n)\,\frac{\partial T}{\partial s},
\]

with the exact electron and neutral conductivity formulas from `physics.hpp`
and one shared `numpy.gradient(T, h_km*1e3)` implementation. The plotted value
excludes TRAC broadening, numerical diffusivity, imposed boundary flux, and all
other solver-only terms.

## Implementation

- `physics.hpp` exposes scalar physical electron, neutral, and total
  conductivities plus a separately named solver-effective diagnostic.
- `chromo_main.cpp` extends gamma sidecars from 9 to 10 columns. The stable
  physical-field prefix is followed by both `kappa_physical` and
  `kappa_solver`; the header names both explicitly.
- `util/conductive_flux.py` is the single plotting implementation of the
  conductivity components and cell-centred physical heat flux.
- `util/animate_isentropic.py` and `util/animate_gentle_column.py` route both
  closure modes through that helper and explicitly label the panel as physical
  flux with TRAC and numerical diffusion excluded. Legacy nine-column gamma
  sidecars remain plottable because physical conductivity is recomputed from
  `T`, `n_e`, and `n_HI`, rather than interpreting their old last column as
  physical.
- `util/plot_physical_flux_compare.py` produces a matched-time, matched-grid
  comparison and machine-readable metrics.
- `cmake/check_gamma_output.cmake`, `tests/chromo_tests.cpp`, and
  `tests/test_conductive_flux.py` cover the new schema and solver-term
  independence.

## Matched comparison at 1000 s

The gamma-table snapshot exists at exactly 1000 s. The fixed-gamma conserved
state is linearly interpolated between its adjacent stored frames to exactly
1000 s. Both use the identical 2000-cell height array spanning 1.0763 to
2152.6 km and identical plot limits.

Below 500 km:

| maximum | fixed gamma=1.05 | gamma table | gamma vs fixed |
|---|---:|---:|---:|
| `abs(dT/ds)` [K/m] | 1.37072e-2 | 1.84893e-2 | +34.89% |
| `n_e` [m^-3] | 8.27432e19 | 7.40315e19 | -10.53% |
| `n_HI` [m^-3] | 1.18532e23 | 1.18666e23 | +0.11% |
| `kappa_e` [W m^-1 K^-1] | 1.16995e-2 | 1.08629e-2 | -7.15% |
| `kappa_n` [W m^-1 K^-1] | 1.62429 | 1.62360 | -0.04% |
| `abs(q_physical)` [W/m^2] | 2.24248e-2 | 3.02201e-2 | +34.76% |

The lower-atmosphere difference therefore remains after removing the mixed
diagnostic, but it is not a conductivity/TRAC/numerical-diffusion effect.
Neutral conductivity dominates and differs by only 0.04%; the 34.76% flux
maximum difference tracks the 34.89% steeper maximum temperature gradient in
the gamma-table run. Its lower electron density also reduces `kappa_e`, but
`kappa_e` is subdominant there. These are maxima over the region and need not
occur at the same cell.

Artifacts:

- `visualization/model_column/iso_t22k_fixed105_vs_gamma_physical_flux_t1000.png`
- `outputs/model_column/iso_t22k_fixed105_vs_gamma_physical_flux_t1000_metrics.json`

## Validation

- `cmake --build build -j4`: passed.
- `./build/chromo_tests`: 7,211 checks passed, 0 failed.
- `ctest --test-dir build --output-on-failure`: 8/8 passed, including the gamma
  sidecar smoke test and the shared Python diagnostic test.
- Two-frame gamma-table animation smoke render: passed.
- `git diff --check`: passed.
- The new PNG was visually inspected; labels, curves, matched axes, and the
  below-500 km panels render correctly.
