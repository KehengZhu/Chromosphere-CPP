# Output file formats {#io_formats}

Everything `chromo_main` writes is line-oriented ASCII. This page is the on-disk layout of each file. The writer is `chromo_main.cpp`; it is the authority for every column name quoted here. See @ref driver for the driver's control flow and @ref architecture for where the pieces live.

A run writes at most five files: the main snapshot at the path given as `argv[1]` (default `outputs/output.txt`), three sidecars derived from that path by suffix, and one appended log.

| File | Enabled by | Default |
| --- | --- | --- |
| `<out>` (main snapshot) | `CHROMO_OUTPUT` | on |
| `<out>.gamma_diag` | `CHROMO_GAMMA_DIAG`, release solver only | on |
| `<out>.faceflux` | `CHROMO_FACE_FLUX_DIAG`, release solver only | off |
| `<out>.outercond` | `CHROMO_OUTER_COND_DIAG`, release solver only | off |
| `outputs/output.log` | always | always |

`CHROMO_FACE_FLUX_DIAG` and `CHROMO_OUTER_COND_DIAG` throw at startup if the run is not in release mode, because their capture points live in `mixture_rhs_explicit` and `mixture_apply_conduction`.

## Main snapshot file

@verbatim
661 3
  1601.1  1602.21  1603.31   ...   (ns values)
# EOS_MODE=gamma_table diagnostics=outputs/model_column/run.txt.gamma_diag
# t = 0 step = 0
  9.70534e-10  0  0.447124
  9.65854e-10  0  0.446071
  ...                                (ns rows)
# t = 0.498707 step = 89
  ...                                (ns rows)
@endverbatim

**Line 1** is `ns state_rows`, separated by one space: the number of physical cells and the width of the conserved state.

**Line 2** is `ns` cumulative cell heights in kilometres, each preceded by two spaces. The value for cell `i` is

@verbatim
    h_i = out_base_km + sum_{j <= i} ds_j * 1e-3
@endverbatim

so it is the **upper face** height of cell `i`, not the cell centre — `ds_i` is in metres, hence the `1e-3`. `Grid::out_base_km` is the absolute height of the domain base and exists only to label these values. It defaults to `1003.0` (the Model C7 base); `model_column` sets it to its own actual base. It never enters the solver.

**The optional `# EOS_MODE=` line** is written next, and only when the run is in release mode *and* the `.gamma_diag` sidecar is enabled. Its exact form is `# EOS_MODE=gamma_table diagnostics=<out>.gamma_diag`. A reader that assumes the first `#` line after the heights is a frame marker must skip it.

**Then, repeated per frame:** a marker line `# t = <t> step = <step>` followed by exactly `ns` rows. Each row is one cell, in increasing cell index, holding `state_rows` values each preceded by two spaces, written in row-index order `k = 0 .. state_rows-1`. Values are the stored conserved-state numbers (`chromosphere::Real`, `double` in the release) printed at the default stream precision, i.e. **six significant digits** — so the on-disk format is a fixed ASCII text format and is unaffected by the storage precision the solver was built with. There is no binary state I/O anywhere in this code. The final state is always written even when it does not land on the output cadence, so the last frame may not be a multiple of the stride.

`state_rows` takes exactly two values.

### `state_rows == 3` — release snapshot

The single-fluid equilibrium-mixture state `U = (rho, rho u, E)`, indices `mix::`.

| Column | Symbol | Meaning | Units |
| --- | --- | --- | --- |
| 0 | `mix::RHO` | total mass density | kg m<sup>-3</sup> |
| 1 | `mix::MOM` | field-aligned momentum density `rho u` | kg m<sup>-2</sup> s<sup>-1</sup> |
| 2 | `mix::ENERGY` | total energy density `E = e_int(rho,T) + 1/2 rho u^2 + rho phi_g` | J m<sup>-3</sup> |

`e_int` includes the hydrogen ionization energy. The ionization fraction `x`, `n_e`, `n_HI`, `p_e` and `T` are **not** in this file — they are derived EOS diagnostics, and the `.gamma_diag` sidecar is where the driver writes them out.

### `state_rows == 7` — two-fluid snapshot

The historical carrier state, indices `cons::`.

| Column | Symbol | Meaning | Units |
| --- | --- | --- | --- |
| 0 | `cons::RHO_I` | ion mass density `rho_i` | kg m<sup>-3</sup> |
| 1 | `cons::RHO_N` | neutral mass density `rho_n` | kg m<sup>-3</sup> |
| 2 | `cons::MOM_I` | ion momentum density `rho_i V` | kg m<sup>-2</sup> s<sup>-1</sup> |
| 3 | `cons::MOM_N` | neutral momentum density `rho_n U` | kg m<sup>-2</sup> s<sup>-1</sup> |
| 4 | `cons::E_I` | **total** charged energy (protons + electrons + bulk KE + gravity) | J m<sup>-3</sup> |
| 5 | `cons::E_N` | neutral total energy | J m<sup>-3</sup> |
| 6 | `cons::E_E` | electron internal energy `(3/2) p_e` | J m<sup>-3</sup> |

`E_I` is the total charged energy, not the proton energy: the proton temperature is derived as `T_i = (p_total - p_e)/(n_i k_B)`, and `E_E` never feeds back into `E_I`, momentum or density, so a run with `ENABLE_TE=0` reproduces the single-temperature baseline exactly.

## `<out>.gamma_diag`

The EOS-decoded cell-centred state of a release run. Written independently of the main snapshot: `CHROMO_OUTPUT=0` with `CHROMO_GAMMA_DIAG=1` produces this file and no main file.

@verbatim
661 10
  1601.1  1602.21  1603.31   ...   (ns values, same convention as the snapshot)
# EOS_MODE=gamma_table
# GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat
# GAMMA_TABLE_SHA256=758dfe315cc907097a15368c1cfe242d279251c1f9c2f72ac1df7544d1635908
# potential=cell_center_mean_of_phi_g_imh_phi_g_iph
# release conduction is PHYSICAL only, so kappa_solver is
# identically kappa_physical; the column is kept for format
# stability with existing analysis scripts.
# columns=rho_total v_cm T x_eq n_e n_HI p_total Gamma1 kappa_physical kappa_solver
# t = 0 step = 0
  9.70534e-10  0  6631.35  0.273169  1.58505e+17  4.21742e+17  0.067637  1.09166  1.2196  1.2196
  ...                                (ns rows)
@endverbatim

Line 1 is `ns 10`; line 2 repeats the snapshot's height line verbatim (same `out_base_km` offset, same upper-face convention). Then eight `#` provenance lines:

- `# EOS_MODE=gamma_table` — the closure in force.
- `# GAMMA_TABLE=<path>` — the table path exactly as `GAMMA_TABLE` supplied it.
- `# GAMMA_TABLE_SHA256=<64 hex>` — SHA-256 of the table's bytes, computed by the driver at run time (its own implementation, no external tool). This pins which table produced the run.
- `# potential=cell_center_mean_of_phi_g_imh_phi_g_iph` — the gravitational-potential convention used to decode the state. `mixture_cell_phi(grid,i) = 0.5*(phi_g_imh(i) + phi_g_iph(i))`; `E` carries `rho phi_g`, so the decode is only meaningful with the matching gauge.
- Three lines recording that release conduction is physical only, so `kappa_solver` is identically `kappa_physical`. The column is retained solely for format stability with existing analysis scripts and carries no independent information.
- `# columns=...` — the ten column names.

Then, per frame, a `# t = <t> step = <step>` marker and `ns` rows of 10 values, two spaces before each. Every row is a fresh `decode_equilibrium_mixture` of the corresponding snapshot cell with **no** temperature hint (a NaN seed), so it is a clean re-derivation rather than a solver-internal cache dump.

| Column | Name | Meaning | Units |
| --- | --- | --- | --- |
| 0 | `rho_total` | total mass density | kg m<sup>-3</sup> |
| 1 | `v_cm` | field-aligned velocity `rho u / rho`. The name is historical; the value is SI, positive along increasing `s` | m s<sup>-1</sup> |
| 2 | `T` | common temperature | K |
| 3 | `x_eq` | Saha equilibrium ionization fraction `n_e / n_H` | — |
| 4 | `n_e` | electron number density `x n_H` | m<sup>-3</sup> |
| 5 | `n_HI` | neutral hydrogen number density `(1-x) n_H` | m<sup>-3</sup> |
| 6 | `p_total` | total pressure `(1+x) n_H k_B T` | Pa |
| 7 | `Gamma1` | CRASH adiabatic (acoustic) index | — |
| 8 | `kappa_physical` | `physical_conductivity(n_e, n_HI, T)` = `physical_kappa_e + physical_kappa_n` | W m<sup>-1</sup> K<sup>-1</sup> |
| 9 | `kappa_solver` | the identical call; see the header note | W m<sup>-1</sup> K<sup>-1</sup> |

## `<out>.faceflux`

A read-only copy-out of `Grid::face_flux_capture` (`MixtureFaceFluxCapture`), filled by `mixture_rhs_explicit` only when `Grid::capture_face_flux` is true and **never read back by any solver stage**. Enabling it therefore changes no numerical result bit-for-bit: the capture arms before the step's single RHS call, the record is written after the step, and the solver's arithmetic is untouched.

What it captures is one explicit RHS evaluation: the decoded cell-centred inputs the reconstruction saw, the one-sided reconstructed face states, the limiter inputs and outputs actually used for that face, and the production total-mass numerical flux with its central/dissipative split. Because every column comes from the same RHS call, one record is self-consistent and needs no cross-matching against the snapshot file.

**Face indexing.** Row `i` describes the **upper face `i+1/2` of cell `i`** — exactly the face the continuity row differences. `L` is the state extrapolated from cell `i`, `R` the state extrapolated from cell `i+1`. The last row is the outer boundary face. `cell_km` is the cell **centre** height; `face_km` is that upper-face height, which equals entry `i` of the snapshot's line 2.

**Cadence and extent.** `CHROMO_FACE_FLUX_STRIDE` is the capture cadence in steps (default: the snapshot stride). `CHROMO_FACE_FLUX_TOP=N` restricts each record to the top `N` cells, setting `first_cell = ns - N`; `0` (the default) means all cells. The final step of a run is always captured. The `t`/`step` on a record label the **start** of the step — the state the RHS was evaluated at. Numbers are written with `precision(10)`.

@verbatim
# face-flux diagnostic (read-only capture of mixture_rhs_explicit)
# ns=661 first_cell=0 uniform_mesh=0 mc3=1 beta=2 roe=0 swmf_godunov=1 pressure_reconstruct=1 ds_km=1.10477
# mass flux = the mix::RHO continuity row (total mass density; the release state has no carrier rows)
# columns=cell_km face_km rho_cell v_cell T_cell rho_L rho_R v_L v_R T_L T_R cs_L cs_R a_face f_central f_diff f_total r_rho phi_plus_rho r_ip1_rho phi_minus_rho r_v phi_plus_v r_T phi_plus_T p_cell p_L p_R
# t = 0 step = 0
1600.552368 1601.104736 9.705339865e-10 0 6631.354019 ...
@endverbatim

Header keys on line 2: `ns` (cell count), `first_cell` (index of the first row in each record), `uniform_mesh`, `mc3` (MC3/Koren limiter on), `beta` (limiter beta), `roe` (Roe characteristic reference flux on), `swmf_godunov` (release SWMF exact-Riemann Godunov flux on — `1` on a default release run), `pressure_reconstruct` (`(ln rho, u, ln p)` reconstruction on), `ds_km` (width of cell 0, km).

The 28 columns, in order:

| # | Column | Meaning |
| --- | --- | --- |
| 0 | `cell_km` | centre height of cell `i` [km] |
| 1 | `face_km` | upper-face height of cell `i` [km] |
| 2 | `rho_cell` | decoded cell-centred density [kg m<sup>-3</sup>] |
| 3 | `v_cell` | decoded cell-centred velocity [m s<sup>-1</sup>] |
| 4 | `T_cell` | decoded cell-centred temperature [K] |
| 5 | `rho_L` | reconstructed density at `i+1/2` from cell `i` |
| 6 | `rho_R` | reconstructed density at `i+1/2` from cell `i+1` |
| 7 | `v_L` | reconstructed velocity, left state |
| 8 | `v_R` | reconstructed velocity, right state |
| 9 | `T_L` | reconstructed temperature, left state |
| 10 | `T_R` | reconstructed temperature, right state |
| 11 | `cs_L` | **equilibrium** sound speed `sqrt(Gamma1 p/rho)` of the left face state [m s<sup>-1</sup>]. This is the closure derivative, not the frozen `sqrt(5/3 p/rho)` the release Riemann solve and the CFL use. |
| 12 | `cs_R` | same, right face state |
| 13 | `a_face` | equilibrium acoustic spectral radius `max(|v_L|+cs_L, |v_R|+cs_R)`. It is the Rusanov coefficient when Rusanov is selected, and the per-face fallback coefficient of both other solvers; under the release Godunov flux and under Roe it is otherwise diagnostic only. |
| 14 | `f_central` | central (averaged-physical-flux) part of the mass flux |
| 15 | `f_diff` | dissipative part of the mass flux. The release `swmf-godunov` flux and the Roe reference flux are only defined as a whole, so this is whatever they add to `f_central`; in Rusanov mode it is instead the closed form `-1/2 a (U_R - U_L)`. |
| 16 | `f_total` | the production numerical mass flux at this face, `f_central + f_diff` [kg m<sup>-2</sup> s<sup>-1</sup>] |
| 17 | `r_rho` | limiter slope ratio for the density slot, building `L` from cell `i` |
| 18 | `phi_plus_rho` | limiter value applied for that `L` state |
| 19 | `r_ip1_rho` | limiter slope ratio for the density slot, building `R` from cell `i+1` |
| 20 | `phi_minus_rho` | limiter value applied for that `R` state |
| 21 | `r_v` | limiter slope ratio, velocity slot (`L` side) |
| 22 | `phi_plus_v` | limiter value, velocity slot (`L` side) |
| 23 | `r_T` | limiter slope ratio, temperature/pressure slot (`L` side) |
| 24 | `phi_plus_T` | limiter value, that slot (`L` side) |
| 25 | `p_cell` | decoded cell-centred total pressure [Pa] |
| 26 | `p_L` | total pressure of the left corrector face state, rebuilt through the production face builder |
| 27 | `p_R` | total pressure of the right corrector face state |

`p_L` versus `p_R` is the reconstruction-induced pressure mismatch that the `(ln rho, u, ln p)` reconstruction targets; `p_cell`, `p_L` and `p_R` are rebuilt through the production face builder on capture steps only.

## `<out>.outercond`

A read-only copy-out of `Grid::outer_conduction_capture` (`OuterConductionCapture`): the outer-face quantities of the **final converged Newton iteration** of that step's implicit conduction solve — the ones that actually built the accepted `g_right[ns-1]`. Like the face-flux sidecar it is never read back by a solver stage, so enabling it changes no numerical result.

Unlike the other files this one is a **flat table with no per-frame markers**: one row per captured step, with `t` and `step` as the first two columns. `t`/`step` label the start of the step, matching the `.faceflux` convention. Cadence is `CHROMO_OUTER_COND_STRIDE` in steps (default: the face-flux stride); the final step is always captured. A row is emitted only when the capture is valid. Numbers are written with `precision(10)`.

@verbatim
# outer-face conduction diagnostic (read-only capture of the final converged mixture_apply_conduction iteration)
# ns=661 ds_km=0.276042 uniform_mesh=0 impose_outer_heat_flux=0
# release conduction is PHYSICAL only (kappa_e + kappa_n; no TRAC, no artificial diffusivity)
# ds_face = half the top cell when an external face temperature is imposed, else the centre-to-ghost-centre span
# q_face = kappa_face*(T_wall - T_top)/ds_face  [W m^-2, positive = into the top cell]
# columns=t step T_top T_wall kappa_face ds_face area_ratio q_face
0 0 22551.30289 22000 0.6608015671 138.0208282 1 -2.639469825
0.4993866384 89 21827.25415 22000 0.6608015671 138.0208282 1 0.8270543763
@endverbatim

Header keys on line 2: `ns`, `ds_km` (width of the **top** cell `ns-1`, km), `uniform_mesh`, `impose_outer_heat_flux`.

| # | Column | Meaning | Units |
| --- | --- | --- | --- |
| 0 | `t` | physical time at the start of the step | s |
| 1 | `step` | step index at the start of the step | — |
| 2 | `T_top` | top-cell centre temperature | K |
| 3 | `T_wall` | the thermal datum at the outer face | K |
| 4 | `kappa_face` | physical face conductivity `kappa_e + kappa_n` | W m<sup>-1</sup> K<sup>-1</sup> |
| 5 | `ds_face` | distance from the top-cell centre to the thermal datum | m |
| 6 | `area_ratio` | flux-tube factor `B_i(ns-1) / B_iph(ns-1)`, converting the face flux into the cell's volumetric divergence | — |
| 7 | `q_face` | outer-face conductive flux | W m<sup>-2</sup> |

### Sign and geometry conventions

@verbatim
    q_face = kappa_face * (T_wall - T_top) / ds_face
@endverbatim

**`q_face > 0` means heat flowing into the top cell** (the reservoir is hotter than the top cell); `q_face < 0` means the column is losing heat through the outer face. In the sample above the very first row is negative because the initial top cell sits above the 22,000 K reservoir and relaxes down onto it.

`ds_face` depends on which outer closure the scenario chose:

- When the scenario imposes an external conductive reservoir temperature at the **physical** outer face, `T_wall` is that imposed value (22,000 K in the release configuration), `kappa_face` is evaluated at that face state, and `ds_face = 0.5 * ds_i(ns-1)` — the top-cell centre-to-face distance.
- Without that override, the datum is the hydro ghost **centre** and `ds_face = ds_iph_i(ns-1)`, the centre-to-ghost-centre span.

The `imposed_neumann` flag of the capture struct is not written to the file; `impose_outer_heat_flux` in the header records the same configuration choice.

## `outputs/output.log`

Opened in **append** mode at the hard-coded relative path `outputs/output.log`, so the directory must exist relative to the working directory. Every run appends exactly one line:

@verbatim
[outputs/model_column/run.txt mode=full] Total time = 100
@endverbatim

That is the output path as given on the command line, the `mode` argument, and the resolved `total_time` in physical seconds (either `time_mult * 10 * sum(ds_i) / 2e4` or the `CHROMO_T_END` override). Nothing else is written to this file — the per-step progress lines, the profile report and the `termination=` block go to stdout, and the OpenMP banner and warnings to stderr.

## Godunov fallback tally on stdout

Every run using the release Godunov flux — i.e. every default `model_column` run — prints four extra stdout lines just before the profile report, so a fallback to Rusanov can never pass unnoticed:

@verbatim
godunov.faces=23506482
godunov.exact=23506482
godunov.fallbacks=0 (bad_input=0 vacuum=0 negative_p=0 no_converge=0 bad_sample=0)
godunov.iterations_max=1 mean=1
@endverbatim

`faces` counts every face solve attempted over the whole run (two per cell per RHS evaluation), `exact` those the exact Riemann solve supplied, and the parenthesised breakdown attributes each fallback to its cause. `iterations_max` / `mean` are the Newton-iteration counts of the star-pressure solve. A nonzero `fallbacks` additionally emits a stderr `WARNING`, and `scripts/release_validation.sh` fails the release smoke on it. These lines are absent only from a run that overrode the flux with `ISO_RIEMANN=roe-local` or `ISO_RIEMANN=rusanov`.

## Reading these in Python

There is no shared I/O package in `util/`: scripts share code by `sys.path`-inserting the `util/` directory and importing each other directly, so the canonical reader for each format is whichever module the others import from.

| Format | Canonical reader |
| --- | --- |
| main snapshot | `util/plot_output.py` → `read_frames(path)` → `(heights, [(t, step, xn[ns, state_rows])])`. Imported by `animate_output.py`, `animate_compare.py`, `animate_cooling_compare.py`, `animate_ionization_compare.py`, `visualize_canopy.py`, `visualize_pfss_lines_evolution.py`, `visualize_rho_total.py`, `visualize_rho_total_compare.py`. Tolerates a truncated trailing frame. |
| main snapshot, loop / field-line runs | `util/animate_loop_flare.py` → `read_frames(path)`, with `parse_dat(path)` and `read_meta(path)` for the `[META]`/`[CELLS]`/`[GHOSTS]` mesh file. Imported by `animate_gentle.py`, `gentle_evap_profiles.py`, `visualize_event_3d.py`. `util/synth_doppler.py` carries a near-duplicate. |
| `.gamma_diag` | `util/upper_bc_decouple_diag.py` → `load_gamma_diag(path)` → `(heights_km, [(t, arr[ns, 10])])`, plus `COLUMNS` and `pick(frames, t)`. Appends the suffix if absent and stops at the first non-finite or short frame. Imported by `resolution_scan_diag.py` and `outer_refine_diag.py`. |
| `.gamma_diag`, space-time form | `util/column_spacetime.py` → `load(path)` → `(h_km[ns], t[nt], d[nt, ns, 10])` with an `.npz` cache next to the source, plus `col()`, `at_height()`, `band()`, `xcorr_lag()`. Imported by `plot_column_sloshing.py`. |
| `.faceflux` | `util/face_flux_diag.py` → `read_faceflux(path)` → `(meta, [(t, step, {column: array})])`. Takes the column names from the sidecar's own `# columns=` line and falls back to its module-level `COLUMNS` for pre-header files. Imported by `resolution_scan_diag.py`, `numerical_conduction_diag.py`, `outer_refine_diag.py`, `animate_face_mass_flux.py`. `numerical_conduction_diag.read_faceflux_tolerant` wraps it for sidecars truncated by a crashed run. |
| `.outercond` | `util/numerical_conduction_diag.py` → `read_outercond(path)` → dict of column arrays plus `_meta`, with `pick_outer(oc, t)`. |

Two caveats hold for the current writer:

- `plot_output.py::read_frames` expects the first `#` line after the height line to be a frame marker, so it does not skip the optional `# EOS_MODE=` line a release run writes. Readers of release snapshots must skip `#` lines that are not `# t = ...`.
- `numerical_conduction_diag.py::read_outercond` accepts only rows of 12 (`OUTER_COLS`) or 11 (`LEGACY_OUTER_COLS`) values, both of which predate the release: they still carry `chi_num_face` / `kappa_num_face` / `q_phys` / `q_num`. The current release sidecar has the 8 columns documented above, and this reader will find no data rows in it.

Many other `util/` scripts re-implement snapshot parsing inline (`plot_gentle_summary.py`, `visualize_flare.py`, `verify_pfss_flare.py`, `compare_two_fluid.py`, and others). Those are per-script copies, not the shared format definition.
