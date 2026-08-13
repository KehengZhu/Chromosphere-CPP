# Naming conventions {#naming}

The goal is that a C++ symbol should be readable straight off the write-up's symbol. One axis — grid location — lives on the suffix; the other — physical quantity — lives on the body.

## Grid-location suffix

| Suffix | Length | Meaning |
| --- | --- | --- |
| `_i` | `ns` | one scalar per cell, cell-centered |
| `_state` | `ns * num_of_eq` | packed state, all equations stacked per cell |
| `_iph` | `ns` | scalar at the `i+1/2` (right) face |
| `_imh` | `ns` | scalar at the `i-1/2` (left) face |
| `_state_iph` / `_state_imh` | `ns * num_of_eq` | packed face state |
| `_ip1` / `_im1` / `_ip2` / `_im2` | matches input | shifted neighbour (output of `ip1`, `im1`, ...) |

The suffix carries grid location only. Packed-versus-scalar is spelled out as `_state` rather than doubling the `i`, so `xn_state` and `xn_i` are not one keystroke apart.

## Variable bodies

| Concept (write-up) | C++ |
| --- | --- |
| ion velocity `V` | `V` |
| neutral velocity `U` | `U` |
| ion / neutral density `rho_i`, `rho_n` | `rho_i`, `rho_n` |
| ion / neutral momentum `rho_i V`, `rho_n U` | `rhoV_i`, `rhoU_n` |
| ion / neutral energy `e_i`, `e_n` | `e_i`, `e_n` |
| ion / neutral pressure `p_i`, `p_n` | `p_i`, `p_n` |
| ion / neutral temperature `T_i`, `T_n` | `T_i`, `T_n` |
| ion / neutral number density `n_i`, `n_n` | `n_i`, `n_n` |
| gravitational potential `phi_g` | `phi_g` |
| ion-neutral collision frequency `nu_in` | `nu_in` |
| `alpha = rho_i nu_in` | `alpha` |
| `kappa_e + kappa_i` at a face | `Kei_iph` |
| ratio of specific heats | `gamma_mono` |
| elementary charge | `q_e` |
| state-vector length `ns * num_of_eq` | `n_state` |

The species marker stays on the right so it pairs cleanly with the grid suffix: `rho_i_iph` reads as "rho_i at i+1/2".

Pi is deliberately not a global. Use `arma::datum::pi` directly, so nothing shadows a local `p_i`.

## Index constants

Conserved and primitive indices live in separate namespaces, so the same short name cannot silently collide at index 0.

| Variable | Conserved | Primitive |
| --- | --- | --- |
| `rho_i` | `chromosphere::cons::RHO_I` | `chromosphere::prim::RHO_I` |
| `rho_n` | `chromosphere::cons::RHO_N` | `chromosphere::prim::RHO_N` |
| ion momentum / velocity | `cons::MOM_I` (= `rho_i V`) | `prim::V` |
| neutral momentum / velocity | `cons::MOM_N` (= `rho_n U`) | `prim::U` |
| ion energy / pressure | `cons::E_I` | `prim::P_I` |
| neutral energy / pressure | `cons::E_N` | `prim::P_N` |
| electron energy / pressure | `cons::E_E` | `prim::P_E` |

The release solver has its own, much smaller set in `chromosphere::mix`: `RHO`, `MOM`, `ENERGY`. The two index sets are never mixed — they belong to solvers with different state widths.

Note that `cons::E_I` is the **total** charged-fluid energy and `prim::P_I` the **total** charged pressure (protons plus electrons); the proton-only quantities are derived by subtracting the electron parts. See @ref physics_model.

## Function names

One verb per role:

| Role | Verb | Examples |
| --- | --- | --- |
| compute a derived physical or numerical quantity | `cal_*` | `cal_flux_state`, `cal_source_state`, `cal_spectral_radius_state`, `cal_max_v_i`, `cal_dt_i` |
| right-hand side of a conservation law | `rhs_*` | `rhs_explicit_state`, `rhs_implicit_state` |
| time integrator | `advance_*` | `advance_Euler_state`, `advance_RK4` |
| pure accessor / packer | `get_*` / `scalar_to` | `get_scalar`, `scalar_to` |
| index shift | bare name | `ip1`, `im1`, `ip2`, `im2` |

The release solver prefixes its entry points with `mixture_` instead — `mixture_decode`, `mixture_rhs_explicit`, `mixture_timestep`, `mixture_apply_conduction`, `mixture_advance` — which is what makes it visible at a call site which solver a line belongs to.

## Units

SI throughout, with two conventional exceptions that are always spelled out in the name or the comment: heights are quoted in kilometres in scenario parameters and output files (`_km` suffix, `Grid::out_base_km`), and a handful of tabulated rate coefficients are digitized in the cgs/eV units of their source paper and converted inline at the point of use, with the conversion factor commented.
