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
physics.hpp              Inline collision frequency (nu_in), conductivities (kappa_e, kappa_n),
                         Stage E ionization/recombination rates (Voronov 1997, Hummer 1994)
chromo_main.cpp          Main entry point — CLI dispatcher
src/
  grid.cpp               Grid::init, Grid::broadcast
  state.cpp              cons2prim, prim2cons, get_scalar, scalar_to, ip1/im1/ip2/im2, flux_lim
  flux.cpp               cal_flux_state, cal_spectral_radius_state, cal_source_state
  rhs.cpp                rhs_explicit_state (MUSCL+Rusanov), rhs_implicit_state
  integrators.cpp        advance_Euler_state, advance_RK4, cal_dt_i, cal_max_v_i, Stage E driver
scenarios/
  scenario.{hpp,cpp}     Scenario dispatcher (make_scenario) + shared apply_open_bcs helper
  model_c7.{hpp,cpp}     Model C7 atmosphere IC + photospheric inner BC
  analytic_canopy.{hpp,cpp}  Exponential-canopy B(z) overlay on the C7 thermodynamic profile
  pfss_field_line.{hpp,cpp}  Tabulated PFSS field-line IC (reads scenarios/data/*.dat)
  data_file_parser.{hpp,cpp} ASCII [META]/[CELLS]/[GHOSTS] parser used by pfss_field_line
  data/                  Generated field-line .dat files (gitignored except .gitkeep)
tests/chromo_tests.cpp   Test suite (physics + numerics + scenarios)
include/armadillo        Vendored header-only Armadillo
util/                    Python visualizations + PFSS extraction pipeline
fortran/                 Legacy Fortran main + SWMF couplers (not built)
docs/                    Al Shidi 2019 paper + Keheng's writeup
outputs/                 Default destination for chromo_main snapshot files (.txt, .log)
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

### CLI

```
build/chromo_main [output_path] [mode] [ionization] [scenario] [data_path] [time_mult]
```

Each argument is positional and optional; defaults are shown in parentheses.

| Arg | Position | Values | Default | Meaning |
| --- | --- | --- | --- | --- |
| `output_path` | 1 | any path | `outputs/output.txt` | snapshot file (line 1 = `ns num_of_eq`; line 2 = cumulative heights in km; then `# t = T step = S` markers each followed by `ns` rows of conserved variables) |
| `mode`        | 2 | `full` \| `explicit` | `full` | `full` = semi-implicit driver (Stage A explicit, B drag, C T-equil, D conduction). `explicit` zeroes `R_I` and runs pure explicit Euler |
| `ionization`  | 3 | `ionization` \| `no-ionization` | `ionization` | toggles Stage E (Voronov 1997 ionization + Hummer 1994 recombination) |
| `scenario`    | 4 | `model_c7` \| `analytic_canopy` \| `pfss_field_line` | `model_c7` | which IC + BC pair to dispatch (see [Scenarios](#scenarios)) |
| `data_path`   | 5 | path to `.dat` \| `""` | `""` | required for tabulated scenarios (`pfss_field_line`); ignored otherwise |
| `time_mult`   | 6 | float | `1.0` | multiplier on the default total simulation time `10·L/Cs`. Step cap scales with this so longer runs aren't truncated |

A companion `outputs/output.log` is appended with the configured `total_time`.

### Examples

```sh
# Default: Model C7, semi-implicit, ionization ON
build/chromo_main

# Ionization OFF baseline for comparison
build/chromo_main outputs/out_c7_off.txt full no-ionization model_c7

# Run 10× longer (useful for relaxation studies)
build/chromo_main outputs/out_c7_10x.txt full ionization model_c7 "" 10.0

# Pure-explicit run (drops the implicit drag / conduction / T-equilibration stages)
build/chromo_main outputs/out_c7_explicit.txt explicit ionization model_c7

# Analytic exponential-canopy B(z) overlay on the C7 thermodynamic profile
build/chromo_main outputs/out_canopy.txt full ionization analytic_canopy

# PFSS field line generated from a GONG synoptic magnetogram
build/chromo_main outputs/out_pfss.txt full ionization pfss_field_line \
    scenarios/data/pfss_qs_20190801.dat
```

To compare ionization on/off side-by-side after a run, point the 8-panel
animation at the two snapshot files:

```sh
.venv/bin/python util/animate_ionization_compare.py \
    outputs/out_c7_10x_on.txt outputs/out_c7_10x_off.txt \
    util/c7_openbc_8panel_compare.mp4 24
```

Other ready-made viz scripts in `util/`: `visualize_rho_total.py`,
`visualize_rho_total_compare.py`, `visualize_canopy.py`, `visualize_pfss_3d.py`,
`visualize_pfss_lines_evolution.py`, `plot_output.py` (last-frame 8-panel
snapshot).

## Scenarios

A *scenario* fixes the field-line geometry, initial conditions, and per-step
boundary update rule. Scenarios share a single `Grid` (cells, B-field profile,
φ_g, ghost buffers) and a single dispatcher:

```cpp
auto sc = make_scenario(name, data_path);   // "model_c7" | "analytic_canopy" | "pfss_field_line"
Grid grid;
grid.init(sc.peek_ns(), cfl);
grid.enable_ionization = true;
Vec xn = sc.ic(grid);                       // initial conserved state
// per step:
sc.update_bc(grid, xn);
xn = advance_Euler_state(grid, xn, dt);
```

The shared `apply_open_bcs` helper in [scenarios/scenario.cpp](scenarios/scenario.cpp)
provides the default open-hyperbolic BC (inner reflecting wall + outer pure-Neumann
outflow); individual scenarios call it then overwrite ghosts that need scenario-specific
treatment.

### `model_c7`

Source: [scenarios/model_c7.cpp](scenarios/model_c7.cpp).

Quiet-sun atmosphere from Avrett & Loeser (2008, *ApJS* 175, 229), tabulated
at 15 heights from 1003 km to 1989 km above τ₅₀₀₀ = 1.

- **IC** — cubic-spline interpolation of (n_e, n_n, T) onto a uniform 100-cell
  grid; V = U = 0.
- **B-field** — `B ≡ 1 T` uniform, `∂(1/B)/∂s ≡ 0` (geometry overlays live in
  the other two scenarios).
- **Gravity** — φ_g(s) = g·(h(s) − h_base) with g = 274 m/s². This breaks
  ideal-gas HSE of C7 by ~25% (C7 was built for NLTE radiative balance), so
  the atmosphere relaxes; the resulting transient is the "quiet-sun
  relaxation" experiment.
- **Outer BC** — pure Neumann on ρ (via `apply_open_bcs`) with two overlaid
  cascades. *Velocity*: V_g0 = V/2, V_g1 = V/4 (same for U) — the jump across
  face ns−½ keeps the Rusanov flux dissipative. *Temperature*: T_g0 = T,
  T_g1 = 2·T for both species — the sharp T-jump at face ns+½ mimics the
  chromosphere→corona transition-region temperature rise.
- **Inner BC** — photospheric reservoir: ρ_n and T_n Dirichlet at the C7 base
  snapshot; n_i, T_i, V, U Neumann from cell 0 (open lower boundary).

### `analytic_canopy`

Source: [scenarios/analytic_canopy.cpp](scenarios/analytic_canopy.cpp).

Same thermodynamic IC as `model_c7`, but overlays the consensus
exponential-canopy magnetic-field recipe (Bellot Rubio & Orozco Suárez 2019;
Martínez-Sykora et al. 2019):

```
B(z) = B_∞ + (B_0 − B_∞) · exp(−(z − z_base) / H_B)
```

Defaults: `B_0 = 100 G` (footpoint), `B_∞ = 15 G` (canopy-merged), `H_B = 300 km`.

`∂(1/B)/∂s` is now non-zero, so the flux-tube expansion source term in
[src/flux.cpp](src/flux.cpp) becomes active — this is the **only** physics
difference from `model_c7`. By flux conservation the tube cross-section grows
as `A(z) ∝ 1/B(z)`, expanding by ~6.7× across the chromospheric domain.

Visualization: [util/visualize_canopy.py](util/visualize_canopy.py) renders
the 3D "trumpet" flux tube + B(z) + A(z) + 4-panel time evolution.

### `pfss_field_line`

Source: [scenarios/pfss_field_line.cpp](scenarios/pfss_field_line.cpp),
[scenarios/data_file_parser.cpp](scenarios/data_file_parser.cpp).

Reads a tabulated `[META]`/`[CELLS]`/`[GHOSTS]` ASCII file produced by
[util/extract_field_line.py](util/extract_field_line.py) from a GONG synoptic
magnetogram via PFSS extrapolation. One-time setup:

```sh
bash util/setup_venv.sh
bash util/fetch_magnetogram.sh
.venv/bin/python util/extract_field_line.py \
    --magnetogram util/data/mrzqs190801t0014c2220_229.fits \
    --output     scenarios/data/pfss_qs_20190801.dat \
    --ns 100 --top-height-km 986
```

**Honest caveat**: PFSS resolves r ∈ [1, 2.5] R⊙; the ~1 Mm chromospheric
domain is < 0.3% of that and falls inside the first PFSS radial cell. So
`pfss_field_line` gives a physically-motivated **footpoint |B|** and the
field-line **topology** (open / closed, arc length, inclination), but it does
**not** resolve B(s) inside the chromosphere — for canopy expansion below
~1 Mm use `analytic_canopy`. Closed loops are simulated only on one footpoint
side; the apex and far footpoint are absorbed into the outer BC.

Visualizations: [util/visualize_pfss_3d.py](util/visualize_pfss_3d.py) (global
sphere + local box renders), [util/visualize_pfss_lines_evolution.py](util/visualize_pfss_lines_evolution.py)
(side-by-side comparison of multiple field lines).

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

Currently 34 test cases / 5552 individual checks; all passing.

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
