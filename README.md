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
([src/integrators.cpp:38](src/integrators.cpp#L38)); the implicit branch is
currently disabled at [src/integrators.cpp:57](src/integrators.cpp#L57). RK4 is
also available ([advance_RK4](src/integrators.cpp#L67)) for the explicit-only path.

The 6 conserved variables per cell (writeup eq 61):

| Index | Symbol | Meaning |
| --- | --- | --- |
| `cons::RHO_I` | `ρ_i` | ion mass density |
| `cons::RHO_N` | `ρ_n` | neutral mass density |
| `cons::MOM_I` | `ρ_i V` | ion momentum density |
| `cons::MOM_N` | `ρ_n U` | neutral momentum density |
| `cons::E_I`   | `e_i` | ion total energy (thermal + kinetic + gravitational, including electron contribution) |
| `cons::E_N`   | `e_n` | neutral total energy |

## File layout

```
chromosphere.hpp         Public API: Grid struct, cons::/prim:: indices, function decls
physics.hpp              Inline collision frequency (nu_in) and heat conductivities (kappa_e, kappa_n)
chromo_main.cpp          Main entry point
src/
  grid.cpp               Grid::init, Grid::broadcast
  state.cpp              cons2prim, prim2cons, get_scalar, scalar_to, ip1/im1/ip2/im2, flux_lim
  flux.cpp               cal_flux_state, cal_spectral_radius_state, cal_source_state
  rhs.cpp                rhs_explicit_state (MUSCL+Rusanov), rhs_implicit_state
  integrators.cpp        advance_Euler_state, advance_RK4, cal_dt_i, cal_max_v_i, print_xn
scenarios/
  model_c7.{hpp,cpp}     Model C7 IC, BC update, cubic-spline interpolation, MODEL_C7 table
tests/chromo_tests.cpp   Test suite (physics + numerics)
include/armadillo        Vendored header-only Armadillo
util/                    Plotting helpers (plotting.ipynb, plotting.m)
fortran/                 Legacy Fortran main + SWMF couplers (not built)
docs/                    Al Shidi 2019 paper + Keheng's writeup
CMakeLists.txt           Build config
```

All solver state lives in `chromosphere::Grid` — physical constants, cell-centered
arrays, ghost-buffer storage, and the broadcast `_state` caches. Solver functions
take `const Grid&` (or `Grid&` for time integrators, which write the `dt_state`
scratch buffer). Tests construct a `Grid` per test; no global mutable state.

## Naming conventions

Goal: a C++ symbol should be readable straight off the writeup's symbol. One axis (grid location) on the suffix, one axis (physical quantity) on the body.

### Grid-location suffix

| Suffix | Length | Meaning |
| --- | --- | --- |
| `_i`     | `ns` | one scalar per cell, cell-centered |
| `_state` | `ns·num_of_eq` | packed state — all 6 equations stacked per cell |
| `_iph`   | `ns` | scalar at the i+1/2 (right) face |
| `_imh`   | `ns` | scalar at the i-1/2 (left) face |
| `_state_iph` / `_state_imh` | `ns·num_of_eq` | packed face state |
| `_ip1` / `_im1` / `_ip2` / `_im2` | matches input | shifted neighbour (output of `ip1`/`im1`/…) |

The suffix carries grid location only — packed-vs-scalar is spelled out as `_state` instead of doubling the `i`, so `xn_state` and `xn_i` are no longer one keystroke apart.

### Variable bodies — match the writeup symbol

| Concept (writeup) | C++ |
| --- | --- |
| ion velocity `V` | `V` |
| neutral velocity `U` | `U` |
| ion / neutral density `ρ_i, ρ_n` | `rho_i`, `rho_n` |
| ion / neutral momentum `ρ_i V, ρ_n U` | `rhoV_i`, `rhoU_n` |
| ion / neutral energy `e_i, e_n` | `e_i`, `e_n` |
| ion / neutral pressure `p_i, p_n` | `p_i`, `p_n` |
| ion / neutral temperature `T_i, T_n` | `T_i`, `T_n` |
| ion / neutral number density `n_i, n_n` | `n_i`, `n_n` |
| gravitational potential `φ_g` | `phi_g` |
| ion-neutral collision freq `ν_in` | `nu_in` |
| α = ρ_i ν_in | `alpha` |
| κ_e+κ_i at face | `Kei_iph` |
| ratio of specific heats | `gamma_mono` |
| elementary charge | `q_e` |
| state-vector length `ns·num_of_eq` | `n_state` |

The species marker stays on the right so it pairs cleanly with the grid suffix: `rho_i_iph` reads as "ρ_i at i+1/2". π is not a global — use `arma::datum::pi` directly to avoid shadowing local `p_i`.

### Index constants — namespaced

Conserved- and primitive-variable indices live in separate namespaces, so the same short name (`V`, `P_I`, …) can't silently collide at index 0.

| Variable | Conserved | Primitive |
| --- | --- | --- |
| `ρ_i` | `cons::RHO_I` | `prim::RHO_I` |
| `ρ_n` | `cons::RHO_N` | `prim::RHO_N` |
| ion mom / vel | `cons::MOM_I` (= `ρ_i V`) | `prim::V` |
| neutral mom / vel | `cons::MOM_N` (= `ρ_n U`) | `prim::U` |
| ion energy / pressure | `cons::E_I` | `prim::P_I` |
| neutral energy / pressure | `cons::E_N` | `prim::P_N` |

### Function names — one verb per role

| Role | Verb | Examples |
| --- | --- | --- |
| compute a derived physical or numerical quantity | `cal_*` | `cal_flux_state`, `cal_source_state`, `cal_spectral_radius_state`, `cal_max_v_i`, `cal_dt_i` |
| RHS of a conservation law | `rhs_*` | `rhs_explicit_state`, `rhs_implicit_state` |
| time integrator | `advance_*` | `advance_Euler_state`, `advance_RK4` |
| pure accessor / packer | `get_*` / `scalar_to` | `get_scalar`, `scalar_to` |
| index shift | bare name | `ip1`, `im1`, `ip2`, `im2` |

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
SOLVER_SRCS="src/grid.cpp src/state.cpp src/flux.cpp src/rhs.cpp src/integrators.cpp scenarios/model_c7.cpp"

clang++ -std=c++14 -O2 -I. -Iinclude -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_BLAS \
    $SOLVER_SRCS chromo_main.cpp -o chromo_main

clang++ -std=c++14 -O2 -I. -Iinclude -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_BLAS \
    $SOLVER_SRCS tests/chromo_tests.cpp -o chromo_tests
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
- `cal_max_v_i` matches `cal_spectral_radius_state`
- A uniform, motionless, gravity-free, `B`-uniform state is an *exact* fixed
  point of `advance_Euler_state` and `advance_RK4`
- Total ion and neutral mass are preserved across one explicit step on that
  fixed-point state

Currently 15 test cases / 294 individual checks; all passing.

## Known issues (worth attending to before more physics is added)

- The writeup's "semi-implicit" Picard iteration (writeup §3.7) is unstable
  for the implicit terms actually present here. At Model C7 conditions the
  drag is stiff (`α·dt ≈ 60` on the first step), and Picard iteration of the
  form `u^{n+1,k+1} = u^n + dt·(R_E + R_I(u^{n+1,k}))` only converges when
  `dt·|R_I'| < 1`. The body of `rhs_implicit_state` is correct (drag,
  collisional + frictional heating, and conservative field-aligned heat
  conduction, writeup §4.3–4.4), but the call site in
  [src/integrators.cpp:57](src/integrators.cpp#L57) keeps the residual zeroed —
  the step is pure explicit Euler until a real implicit solve (Newton on the
  `R_I` Jacobian) replaces the Picard iteration.
- `dinvB_ds_i` is allocated and used (writeup pressure-area term in
  `cal_source_state`), but never populated for a non-uniform `B`. `model_c7_ic`
  sets `B ≡ 1` and `dinvB_ds_i ≡ 0`, so the geometric expansion source is
  latent. Populating `dinvB_ds_i` from the actual `B(s)` is the next step
  before stratified or expanding flux-tube runs.

## References

- Al Shidi, Q., et al. (2019). *Time-Dependent Two-Fluid Magnetohydrodynamic
  Model of the Solar Chromosphere.* [docs/Al Shidi et al. - 2019.pdf](docs/Al%20Shidi%20et%20al.%20-%202019%20-%20Time-Dependent%20Two-Fluid%20Magnetohydrodynamic%20Model.pdf)
- Keheng's chromosphere writeup (field-aligned reduction, numerical scheme,
  Model C7 setup): [docs/Keheng_s_chromosphere.pdf](docs/Keheng_s_chromosphere.pdf)
- Upstream CoMFi (MIT license): https://github.com/qalshidi/phd-dissertation
