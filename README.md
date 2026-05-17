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

## Naming conventions

### Current (what the code uses today)

**Grid-location suffix** — encodes how long the vector is and where on the grid it lives. This is the part of the convention that already works well; keep it.

| Suffix | Length | Meaning |
| --- | --- | --- |
| `_i`   | `ns` | one scalar per cell, cell-centered |
| `_ii`  | `ns·num_of_eq` | packed state — all 6 equations stacked per cell |
| `_iph` | `ns` | scalar at the i+1/2 (right) face |
| `_imh` | `ns` | scalar at the i-1/2 (left) face |
| `_iiph` / `_iimh` | `ns·num_of_eq` | packed face state |
| `_ip1` / `_im1` / `_ip2` / `_im2` | matches input | shifted neighbour (output of `ip1`/`im1`/…) |

**Conserved-variable indices** (used with `get_scalar` / `scalar_to`): `CNI`=ρ_i, `CNN`=ρ_n, `CNV`=ρ_i V, `CNU`=ρ_n U, `CEI`=e_i, `CEN`=e_n. **Primitive** counterparts: `PNI`, `PNN`, `PV`, `PU`, `PPI`, `PPN`.

**Local variables inside solver bodies:** `rhoi`, `rhon` (densities); `rhov`, `rhou` (momenta — *not* densities, despite the `rho` prefix); `vv`, `uu` (writeup `V`, `U`); `ei`, `en` (energies); `pi`, `pn` (pressures); `Ti`, `Tn`; `ni`, `nn`; `phig` (φ_g); `nuin` (ν_in); `alpha` (α = ρ_i ν_in); `Ke`, `Ki`, `Kn` (κ_e, κ_i, κ_n); `Kt_iph` (κ_e+κ_i at right face); `RE_*` / `RI_*` (R_E, R_I residuals from writeup).

**Function prefixes:** `cal_*` (compute a derived quantity), `rhs_*` (right-hand side of conservation law), `advance_*` (time integrator), `find_*` / `get_*` (accessor).

**Known clashes / friction:**
- The global float `pi` (π) shadows local `pi` (ion pressure) in every function that uses both. Currently survives only because `arma::datum::pi` is used where π is actually needed.
- `CNV` / `CNU` start with `N` ("number density") but hold *momenta*. Carried over from the [Al Shidi 2019] state-vector layout.
- `PV` / `PU` break the `P + N/P + species` pattern of `PNI` / `PNN` / `PPI` / `PPN`.
- `vv`, `uu` double letters; `e_` (elementary charge) trailing underscore; `gammamono` runs two words together; `num_of_elem` is mixed snake/abbrev.

### Proposed (not yet applied — keep until the user signs off)

Goal: a C++ symbol should be readable straight off the writeup's symbol. One axis (grid location) on the suffix, one axis (physical quantity) on the body.

| Concept (writeup) | Current | Proposed |
| --- | --- | --- |
| ion velocity `V` | `vv` | `V` |
| neutral velocity `U` | `uu` | `U` |
| ion pressure `p_i` | `pi` (shadows π) | `p_i` (and rename the global π constant or drop it for `arma::datum::pi`) |
| neutral pressure `p_n` | `pn` | `p_n` |
| ion momentum `ρ_i V` | `rhov` | `rhoV_i` |
| neutral momentum `ρ_n U` | `rhou` | `rhoU_n` |
| ion / neutral density `ρ_i, ρ_n` | `rhoi`, `rhon` | `rho_i`, `rho_n` |
| ion / neutral energy `e_i, e_n` | `ei`, `en` | `e_i`, `e_n` |
| ion / neutral temperature | `Ti`, `Tn` | `T_i`, `T_n` |
| ion / neutral number density | `ni`, `nn` | `n_i`, `n_n` |
| index for ρ_i (cons) | `CNI` | `IRHO_I` |
| index for ρ_i V (cons) | `CNV` (misleading `N`) | `IMOM_I` |
| index for e_i (cons) | `CEI` | `IE_I` |
| index for V (prim) | `PV` | `IV` (prim-only — namespaced or in separate header) |
| index for p_i (prim) | `PPI` | `IP_I` |
| κ_e+κ_i at face | `Kt_iph` | `Kei_iph` |
| ratio of specific heats | `gammamono` | `gamma_mono` |
| ion-neutral collision freq | `nuin` | `nu_in` (already used in function names) |
| elementary charge | `e_` | `q_e` |
| state-vector length | `num_of_elem` | `n_state` |
| gravitational potential | `phig` | `phi_g` |

**Principles:**
1. Suffix carries grid location only (`_i`, `_ii`, `_iph`, `_imh`). Don't repurpose it.
2. Body matches writeup symbol with `_i` / `_n` standing in for ion/neutral subscripts — `ρ_i V` → `rhoV_i`, `p_i` → `p_i`, `e_n` → `e_n`. The species marker stays on the right, so it pairs cleanly with the grid suffix: `rho_i_iph` (ion density at the i+1/2 face) reads as "ρ_i at iph."
3. No double-letter names (`vv`, `uu`).
4. Index constants get a consistent `I*` prefix (uppercase macro-style) and drop the `C`/`P` muddle: `IRHO_I`, `IMOM_I`, `IE_I`. Keep conserved and primitive index sets in separate headers/namespaces so the same `IV` can't mean two things.
5. No silent shadows of global constants. If `pi` is the float π, the local pressure cannot also be `pi`.

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

- The writeup's "semi-implicit" Picard iteration (writeup §3.7) is unstable
  for the implicit terms actually present here. At Model C7 conditions the
  drag is stiff (`α·dt ≈ 60` on the first step), and Picard iteration of the
  form `u^{n+1,k+1} = u^n + dt·(R_E + R_I(u^{n+1,k}))` only converges when
  `dt·|R_I'| < 1`. The body of `rhs_implicit_ii` is correct (drag,
  collisional + frictional heating, and conservative field-aligned heat
  conduction, writeup §4.3–4.4), but the call site in
  [chromosphere.cpp:393](chromosphere.cpp#L393) keeps the residual zeroed —
  the step is pure explicit Euler until a real implicit solve (Newton on the
  `R_I` Jacobian) replaces the Picard iteration.
- `dinvB_ds_i` is allocated and used (writeup pressure-area term in
  `cal_S_ii`), but never populated for a non-uniform `B`. `model_c7_ic`
  sets `B ≡ 1` and `dinvB_ds_i ≡ 0`, so the geometric expansion source is
  latent. Populating `dinvB_ds_i` from the actual `B(s)` is the next step
  before stratified or expanding flux-tube runs.

## References

- Al Shidi, Q., et al. (2019). *Time-Dependent Two-Fluid Magnetohydrodynamic
  Model of the Solar Chromosphere.* [docs/Al Shidi et al. - 2019.pdf](docs/Al%20Shidi%20et%20al.%20-%202019%20-%20Time-Dependent%20Two-Fluid%20Magnetohydrodynamic%20Model.pdf)
- Keheng's chromosphere writeup (field-aligned reduction, numerical scheme,
  Model C7 setup): [docs/Keheng_s_chromosphere.pdf](docs/Keheng_s_chromosphere.pdf)
- Upstream CoMFi (MIT license): https://github.com/qalshidi/phd-dissertation
