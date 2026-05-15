# Chromosphere2026

A 1.5D field-aligned two-fluid (ion + neutral) MHD solver for the solar
chromosphere. C++14 with the vendored header-only Armadillo in `include/`.

The physics is the field-aligned reduction of the 2D two-fluid MHD model of
[Al Shidi et al. 2019] (Qusai's CoMFi). Equations and numerical scheme are
documented in [docs/Keheng_s_chromosphere.pdf](docs/Keheng_s_chromosphere.pdf).

## Physics

Assuming flow strictly parallel to **B** (`V = V·b̂`, `U = U·b̂`), the 3D
divergence collapses into a 1D arc-length derivative with a B-weighted
Jacobian (writeup eq 61):

```
∂/∂t  (ρ_i, ρ_n, ρ_i V, ρ_n U, e_i, e_n)^T
  + B · ∂/∂s [ flux vector / B ]
  = S
```

Electrons are folded into the ion equations via quasi-neutrality (`n_e = n_i`,
`T_e = T_i`), so ion pressure carries an extra factor of 2: `p_i = 2 n_i k_b T_i`
(writeup eq 38).

**Numerical scheme:** TVD-MUSCL with minmod limiter, Rusanov / local
Lax–Friedrichs flux, 2nd-order time accuracy via a predictor-corrector
reconstruction (writeup §2.6–2.8). Semi-implicit Euler driver
([chromosphere.cpp:369](chromosphere.cpp#L369)); the implicit branch is
currently disabled at [chromosphere.cpp:388](chromosphere.cpp#L388). RK4 is
also available ([advance_RK4](chromosphere.cpp#L398)) for the explicit-only path.

The 6 conserved variables per cell (writeup eq 61):

| Index | Symbol | Meaning |
| --- | --- | --- |
| `CNI` | `ρ_i` | ion mass density |
| `CNN` | `ρ_n` | neutral mass density |
| `CNV` | `ρ_i V` | ion momentum density |
| `CNU` | `ρ_n U` | neutral momentum density |
| `CEI` | `e_i` | ion total energy (thermal + kinetic + gravitational, including electron contribution) |
| `CEN` | `e_n` | neutral total energy |

## File layout

```
chromosphere.{hpp,cpp}   Solver core (flux, MUSCL, Rusanov, time integrators)
chromo_init.cpp          Constants, Model C7 IC, outer-BC update
chromo_main.cpp          Main entry point
util/chromo_util.cpp     Helpers (cons2prim, ip1/im1, get_scalar, ...)
include/armadillo        Vendored header-only Armadillo
fortran/                 Legacy Fortran main + SWMF couplers (not built)
tests/chromo_tests.cpp   Test suite (physics + numerics)
docs/                    Al Shidi 2019 paper + Keheng's writeup
CMakeLists.txt           Build config
```

## Build

### With CMake (preferred)

```sh
brew install cmake          # or: sudo apt-get install cmake — not yet on this machine
mkdir build && cd build
cmake ..                    # configure (Release by default)
make -j                     # build chromo_main and chromo_tests
./chromo_main               # run the simulation
ctest --output-on-failure   # run the test suite (or ./chromo_tests directly)
```

To switch to a Debug build: `cmake -DCMAKE_BUILD_TYPE=Debug ..`.

To regenerate the build after editing `CMakeLists.txt`: re-run `cmake ..` from
`build/` (no need to delete the directory unless you change the generator).

### Without CMake (direct clang++)

If `cmake` is not installed, the build is two short commands:

```sh
clang++ -std=c++14 -O2 -Iinclude -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_BLAS \
    chromosphere.cpp chromo_init.cpp util/chromo_util.cpp chromo_main.cpp \
    -o chromo_main

clang++ -std=c++14 -O2 -Iinclude -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_BLAS \
    chromosphere.cpp chromo_init.cpp util/chromo_util.cpp tests/chromo_tests.cpp \
    -o chromo_tests
```

The `ARMA_DONT_USE_LAPACK`/`ARMA_DONT_USE_BLAS` defines keep Armadillo
header-only — the solver only uses element-wise ops, `sqrt`, `abs`, `min`/`max`,
and `sub2ind`, none of which require LAPACK/BLAS.

## Run

`chromo_main` defaults: `ns = 100`, `CFL = 0.25`, Model C7 initial condition
([docs/Keheng_s_chromosphere.pdf](docs/Keheng_s_chromosphere.pdf) Table 1).
Runs until `t = 10·L/Cs ≈ 474 s` or 10000 steps, whichever comes first.

Outputs:

- `output.log` — total simulation time, header info
- `output.txt` — line 1: cumulative cell heights (km, offset by 1003);
  remaining `ns` lines: the 6 conserved variables per cell, space-separated.

Plotting helpers (untouched): `util/plotting.ipynb`, `util/plotting.m`.

## Tests

```sh
cd build && ./chromo_tests    # or: ctest
```

The suite ([tests/chromo_tests.cpp](tests/chromo_tests.cpp)) covers three
layers:

**Helper invariants**
- `scalar_to` / `get_scalar` are inverses (pack/unpack the flat
  `ns × num_of_eq` layout)
- `ip1` / `im1` shift values one cell over (interior; ghost slots come from
  the boundary buffers, NaN handling included)
- `flux_lim` matches the minmod formula `max(0, min(1, r))` (writeup eq 15;
  the docstring claims `ospre` but the code is minmod)
- `cons2prim ∘ prim2cons = identity`

**Physics formulas vs. writeup**
- Pressure relation `p_i = 2 n_i k_b T_i`, `p_n = n_n k_b T_n` (writeup eq 38)
- Spectral radius at rest = `max(c_s,i, c_s,n)` with
  `c_s,i = sqrt(2γ k_b T_i / m_i)`, `c_s,n = sqrt(γ k_b T_n / m_n)`
  (writeup eq 64/65)
- Ion-neutral collision frequency `ν_in = (2 a₀)² n_n √(8π k_b (T_i+T_n) / m_i)`
  (writeup eq 57, target-density form)
- Electron heat conductivity (writeup eq 53) + sanity at chromosphere conditions
- Neutral heat conductivity (writeup eq 59) + neutral dominates electron at
  `T ≈ 6500 K` (writeup §3.2)

**Numerical correctness**
- `cal_dt_i` produces `dt = CFL · ds / max_v` and the CFL bound is satisfied
- `get_max_v_i` matches `find_spectral_radius_ii`
- A uniform, motionless, gravity-free, `B`-uniform state is an *exact* fixed
  point of `advance_Euler_ii` and `advance_RK4`
- Total ion and neutral mass are preserved across one explicit step on that
  fixed-point state

Currently 15 test cases / 294 individual checks; all passing.

## Known issues (worth attending to before more physics is added)

- [chromosphere.cpp:340-342, 349-350](chromosphere.cpp#L340) — debug `cout`s
  in `rhs_implicit_ii` will spam stdout once the implicit branch is re-enabled.
- [chromosphere.cpp:388](chromosphere.cpp#L388) — implicit RHS is zeroed inside
  the fix-point loop, so `advance_Euler_ii` is effectively explicit-only.
  Re-enable by switching to `RI_np1_kp1_ii = rhs_implicit_ii(xnp1_ii);`.
- `cal_S_ii` (writeup pressure-expansion source) and `rhs_implicit_ii`
  *both* add a `p · B · ∂(1/B)/∂s` term to `CNV`/`CNU`, with different
  coefficients. They are likely double-counting once the implicit branch is on.

## References

- Al Shidi, Q., et al. (2019). *Time-Dependent Two-Fluid Magnetohydrodynamic
  Model of the Solar Chromosphere.* [docs/Al Shidi et al. - 2019.pdf](docs/Al%20Shidi%20et%20al.%20-%202019%20-%20Time-Dependent%20Two-Fluid%20Magnetohydrodynamic%20Model.pdf)
- Keheng's chromosphere writeup (field-aligned reduction, numerical scheme,
  Model C7 setup): [docs/Keheng_s_chromosphere.pdf](docs/Keheng_s_chromosphere.pdf)
- Upstream CoMFi (MIT license): https://github.com/qalshidi/phd-dissertation
