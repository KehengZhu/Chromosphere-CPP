# Scenario Reconciliation & Toggle Cleanup Plan

Status: **EXECUTED (2026-07-11).** All 7 phases implemented and verified: `chromo_tests` passes 6351/6351 (unchanged from baseline); a bare `model_column` run with the bestwb physical knobs is **bit-identical** to the old `model_isentropic` bestwb command; `model_isentropic`/`model_gentle` aliases resolve (the latter reproduces the stable resolved-corona full-physics preset with no NaN); the default `./chromo_main` now runs `model_column`. ~460 lines of diagnostic/dead code removed from `src/integrators.cpp`; `model_gentle.{cpp,hpp}` and `model_isentropic.{cpp,hpp}` deleted; `model_column.{cpp,hpp}` added.

## Goal

The code accumulated three near-duplicate "single field-line, gravity-on, C7-derived chromosphere→corona column" scenarios (`model_c7`, `model_gentle`, `model_isentropic`) as the gentle-evaporation study evolved, plus a large surface of environment-variable toggles — many of which are diagnostics or select between an old and a new numerical method. The result is confusing and hard to reason about.

This plan:

1. **Merges the three column scenarios into one** — a new `model_column` (the current `model_isentropic`, renamed and generalized), which becomes the **default scenario** in `chromo_main`.
2. **Promotes `model_isentropic`'s newest settings to always-on defaults** (the "bestwb" configuration), deleting the env toggles that used to select the stale alternative numerical methods.
3. **Keeps every toggle that selects a distinct *physical* setting**, absorbing `model_gentle`'s physical drivers (the q(T)-flux ramp and the ambient coronal H(s) heating) as options.
4. **Deletes diagnostic-only toggles and their apparatus** (~420 lines of diagnostic code in `src/integrators.cpp`, plus `ISO_TRAC_TC_FIXED`, the base sponge, and the analytic-isentrope toy IC).

### Decisions taken (2026-07-11)

- **Unified scenario name:** `model_column` (with a `model_isentropic` compatibility alias in `make_scenario`, and a `model_gentle` alias mapping to the gentle preset, so existing saved commands keep working).
- **Pruning aggressiveness:** keep alternative *physical* drivers (q(T)-flux ramp, ambient coronal H(s) heating, IC T-boost) as options; drop the analytic-isentrope toy IC and all numeric-only / stale toggles.
- **C7 atmosphere:** kept as a shared library. `model_c7` remains (it is the quiet-Sun baseline AND the source of the `MODEL_C7` table + `c7_full_profile` + Route-B inversion that `model_flare`, `analytic_canopy`, and the unified scenario all depend on). Not merged away.

## Current state (verified)

### Scenario dependency web

- `model_flare.cpp` → calls `model_c7_ic` / `model_c7_update_bc` (flare = C7 + beam).
- `analytic_canopy.cpp` → calls `model_c7_ic` / `model_c7_update_bc` (overrides B-field only).
- `model_isentropic.cpp` → depends on `c7_full_profile` (defined in `model_c7.cpp`).
- `model_gentle.cpp` → thin wrapper over `model_c7_ic(extended=true)` + `model_c7_update_bc`; adds q(T) ramp, coronal H(s) heating, IC T-boost.
- `pfss_field_line.cpp` → self-contained (uses `apply_open_bcs`; only references C7 conventions in comments).

**Consequence:** the C7 atmosphere (table, spline, `c7_full_profile`, Route-B inversion) is load-bearing shared code and must survive. `model_flare` / `analytic_canopy` must keep compiling against `model_c7_*`.

### What is actually run (from README / visualize_commands.md / docs)

- `model_c7` — quiet-Sun baseline & default; actively run.
- `model_isentropic` — the main testbed; dozens of runs.
- `pfss_field_line` — loops, flares (`PFSS_FLARE`/`FLARE_*`), and gentle-on-loops (`GENTLE_*`); actively run.
- `analytic_canopy` — lightly run (mostly a static-figure overlay).
- `model_gentle` — semi-dormant (a few runs; partly a documented negative result).
- `model_flare` — **dormant as a direct scenario**: no `chromo_main … model_flare` command exists anywhere. Flare runs go through `pfss_field_line` + `PFSS_FLARE`. (Kept because tests/plot scripts reference the physics and it is cheap to keep.)

### Test baseline constraints (`tests/chromo_tests.cpp`)

- Directly calls `model_c7_ic`, `model_isentropic_ic`, `pfss_ic`, `analytic_canopy_ic` (+ their `update_bc`/`peek_ns`). Renaming/removing any of those four breaks compilation. `model_gentle`/`model_flare` are **not** called by tests.
- Pins current defaults: `enable_ionization` off (line ~1001), `photoionization_rate == 1e-4` (427), the C7 tr-jump IC leaves `well_balanced==true` / `c7_tr_jump_bc==true` / `impose_outer_heat_flux==enable_radiative_cooling` (1425–1427).
- Only bitwise check is a **run-to-run determinism** test (`test_uniform_regression_refinement_unset`, ~2016) — not a stored golden file. It sets `ISO_C7_IC` + `ISO_LOG_RECON`.
- iso tests read env: `ISO_C7_IC`, `ISO_EQ_WB`, `ISO_LOG_RECON`, `ISO_REFINE_*`, `GRID_REFINE_*`.

**Key safety principle:** the reconstruction/well-balancing switches (`well_balanced`, `log_reconstruct`, `mc3_limiter`, `eq_wb`, inner-WB) stay as **Grid fields defaulting to false**. Only the *unified scenario's IC* forces them true. `model_c7` / `pfss` / `analytic_canopy` / `model_flare` and the regression baseline keep the legacy path byte-for-byte. We do **not** flip global Grid defaults and do **not** delete the legacy solver branches (they are still exercised by the other scenarios and the tests).

## Target architecture

```
scenarios/
  model_c7.{hpp,cpp}        # SHARED C7 atmosphere library + quiet-Sun baseline scenario
                            #   - MODEL_C7 / MODEL_C7_PHOTO tables, interp1_spline
                            #   - c7_full_profile (consumed by model_column)
                            #   - Route-B P_phot inversion  -> extract to a shared helper
                            #   - model_c7_ic / model_c7_update_bc (used by flare + canopy)
  model_column.{hpp,cpp}    # THE unified column scenario (was model_isentropic), DEFAULT.
                            #   real-C7 IC + hydrostatic-ghost/T-jump BC, newest numerics
                            #   hardcoded, absorbs gentle's physical drivers.
  model_flare.{hpp,cpp}     # unchanged (C7 + beam)
  analytic_canopy.{hpp,cpp} # unchanged (C7 + canopy B)
  pfss_field_line.{hpp,cpp} # unchanged
  mesh.{hpp,cpp}            # unchanged (static local refinement)
  data_file_parser.{hpp,cpp}# unchanged
  scenario.cpp              # dispatch: default model_column; aliases model_isentropic,
                            #   model_gentle -> model_column(+gentle preset)
  # DELETED: model_gentle.{hpp,cpp}
```

`chromo_main` default scenario: `model_c7` → **`model_column`**.

## Toggle taxonomy

### A. PROMOTE to always-on (delete the env toggle; unified IC always sets it)

These are `model_isentropic`'s newest numerical methods — the "bestwb" set. In `model_column` they are unconditional; the Grid fields remain (other scenarios keep using them off).

| Old env toggle | Grid field | New behavior |
|---|---|---|
| `ISO_LOG_RECON` | `log_reconstruct` | always true |
| `ISO_MC3` / `ISO_MC3_BETA` | `mc3_limiter` / `limiter_beta` | always true, β = 2 |
| `ISO_EQ_WB` | `eq_wb` | always true |
| `ISO_INNER_WB` | (kInnerWB) | always true |
| `ISO_INNER_WB_RHO` | (kInnerWBRho) | always true |
| `ISO_C7_IC` | — | real-C7 IC is the only IC (see D) |
| `ISO_C7_HSE` | — | hydrostatic density re-integration always on |
| (`well_balanced`) | `well_balanced` | already always true in the scenario |

### B. KEEP — distinct physical settings (env toggle retained, `ISO_` prefix)

| Toggle | Controls |
|---|---|
| `ISO_HEAT_FLUX` | Stage-1 (relax, conduction off) vs Stage-2 (conduction + T-jump) |
| `ISO_COOLING` | radiative sink (Stage R) |
| `ISO_TWO_FLUID` | separate ion/neutral fluids vs single-fluid |
| `ISO_IONIZATION` | Stage-E ionization/recombination network |
| `ISO_TJUMP_A` / `ISO_TJUMP_B` | top ghost-T jump factors (conduction driver strength) |
| `ISO_T_TOP` | absolute imposed top-boundary temperature |
| `ISO_H_BASE` / `ISO_DH` | domain base height / height span |
| `ISO_CORONA` / `ISO_CORONA_TOP_KM` | extend into a resolved corona |
| `ISO_GAMMA` | polytropic index (radiative-thermostat stand-in) |
| `ISO_RELAX_TIME` | length of the Stage-1 relaxation phase |
| `ISO_VCAP` | outer outflow Mach cap |
| `ISO_INNER_T_NEUMANN` | conduction lower BC: Dirichlet vs Neumann |
| `ISO_NS` | grid resolution |
| `ISO_REFINE_*` | static local mesh refinement (see F for alias cleanup) |
| `ISO_TRAC` / `ISO_TRAC_TCHROM` / `ISO_TRAC_TCMAXFRAC` | TRAC TR-broadening (needed for the resolved corona) |

### C. ABSORB from `model_gentle` — distinct physical drivers (new `ISO_`-prefixed names)

These preserve the gentle-evaporation mechanisms (including the one the SHINE poster used). The Grid already has all the backing fields (`impose_outer_heat_flux`, `outer_heat_flux_*`, `enable_coronal_heating`, `coronal_heat_*`).

| New toggle | Was (gentle) | Controls |
|---|---|---|
| `ISO_QFLUX_ENHANCE` | `GENTLE_Q_ENHANCE` | q(T)-flux ramp peak multiplier (alt conduction driver) |
| `ISO_QFLUX_T_ON` | `GENTLE_Q_T_ON` | q-ramp onset |
| `ISO_QFLUX_RAMP` | `GENTLE_Q_RAMP` | q-ramp half-width |
| `ISO_QFLUX0` | `GENTLE_Q0` | steady base q(T) override |
| `ISO_CHEAT_E0` | `GENTLE_HEAT_E0` | ambient coronal H(s) amplitude |
| `ISO_CHEAT_SH_FRAC` | `GENTLE_HEAT_SH_FRAC` | H(s) scale-height fraction |
| `ISO_CHEAT_ENHANCE` | `GENTLE_HEAT_ENHANCE` | H(s) ramp multiplier |
| `ISO_CHEAT_T_ON` / `ISO_CHEAT_RAMP` | `GENTLE_HEAT_T_ON/_RAMP` | H(s) ramp timing |
| `ISO_TBOOST` / `ISO_TBOOST_TC` / `ISO_TBOOST_W` | `GENTLE_T_BOOST*` | one-time IC coronal superheat |

A `GENTLE_* → ISO_*` mapping table goes in the README migration note. `make_scenario("model_gentle")` aliases to `model_column` with the corona-extended + cooling-on preset so a bare `model_gentle` run reproduces the old default configuration.

### D. ABANDON — analytic-isentrope toy IC (delete the branch + its toggles)

The analytic isentrope was the Stage-1 development toy, fully superseded by the real-C7 IC. Delete the `T_of_z` / `n_H ∝ T^{3/2}` branch, the `DH_max` adiabatic cap, and:

- `ISO_T_BASE`, `ISO_N_BASE`, `ISO_F_ION` (only used by the toy IC).

This substantially simplifies `model_column_ic` (one profile path instead of three).

### E. ABANDON — diagnostic-only toggles + apparatus

| Toggle(s) | Apparatus to delete |
|---|---|
| `ISO_DIAG_COND_HYDRO` / `_EVERY` / `_OUT` | `CondHydroDiagnostic` struct + `cond_hydro_diagnostic()` (`src/integrators.cpp` ~434–579) |
| `ISO_DIAG_DYNAMIC_HSE` / `_EVERY` / `_OUT` | `dynamic_hse_*` helpers + `DynamicHseDiagnostic` (`src/integrators.cpp` ~581–869) |
| `ISO_TRAC_TC_FIXED` | `Grid::trac_fixed_cutoff_T` + its use (labeled "Diagnostic only" in the header) |
| `ISO_SPONGE_CELLS` / `ISO_SPONGE_RATE` | `Grid::sponge_base_cells/_rate` + the sponge stage in the integrator (superseded by `eq_wb`; the old sponge figures were a units bug) |

Also remove the diagnostic-recording hooks woven into `apply_conduction_stage` (`record_diagnostic` / `record_dynamic_hse`, ~lines 871–880 and ~1115–1125). Net: **~420 lines removed from `src/integrators.cpp`** (~24% of the file) with no change to physics results.

### F. ABANDON / fold — numeric crutches

- `ISO_DIFF_MULT` → fold into an internal automatic value (corona ⇒ 4, else ⇒ 1); remove the env override. It is a numerical diffusivity crutch, not a physical setting.
- `GRID_REFINE_*` alias → unify on `ISO_REFINE_*` only (update `refine_params_from_env` and `test_refine_params_env_aliases`). Low priority; can defer.

### G. Untouched (outside merge scope)

`ENABLE_TE`, `SINGLE_FLUID`, `CHROMO_CFL` (chromo_main-level); `C7_TJUMP_A/B` (model_c7 library); `FLARE_*` / `PFSS_FLARE` (flare physics); `GRID_REFINE_*`/`ISO_REFINE_*` mesh (kept, see F).

## Execution phases

Each phase ends green (`chromo_main` + `chromo_tests` build; test suite passes) before the next.

### Phase 0 — Safety net
- Record current `chromo_tests` pass count.
- Capture reference outputs for the current best `model_isentropic` "bestwb" run and a `model_gentle` run, to diff physics before/after the merge (expect bit-identical for bestwb since it becomes the default; expect equivalent for gentle via the alias preset).

### Phase 1 — Extract the shared C7 library (no behavior change)
- Extract the Route-B P_phot inversion (currently duplicated inline in `model_c7.cpp` and `model_isentropic.cpp`) into one shared helper in `model_c7.{hpp,cpp}` (e.g. `c7_route_b_photoionization(grid, T, n_i, n_n, h_km) -> Vec`). Have both callers use it.
- Confirm `c7_full_profile` is the single source of the full atmosphere. No functional change.

### Phase 2 — Delete the diagnostic apparatus (`src/integrators.cpp`)
- Remove `CondHydroDiagnostic`, `DynamicHseDiagnostic`, the `dynamic_hse_*` helpers, and the recording hooks in `apply_conduction_stage`.
- Remove `ISO_DIAG_*` reads. Verify tests unaffected (no test references them).

### Phase 3 — Remove the base sponge + `trac_fixed_cutoff_T`
- Delete `Grid::sponge_base_cells/_rate` and the sponge stage; delete `Grid::trac_fixed_cutoff_T` and its use. Remove `ISO_SPONGE_*` / `ISO_TRAC_TC_FIXED`. Verify tests.

### Phase 4 — Create `model_column` from `model_isentropic`
- Rename `model_isentropic.{hpp,cpp}` → `model_column.{hpp,cpp}`; rename `model_isentropic_ic/_update_bc` → `model_column_ic/_update_bc`.
- **Hardcode the newest numerics** (§A): always set `log_reconstruct`, `mc3_limiter` (β=2), `eq_wb`, inner-WB + inner-WB-rho true; C7 IC + HSE always on. Delete the corresponding env reads and the analytic-isentrope branch (§D).
- Simplify `peek_ns` (no analytic-isentrope default; C7 = 600, corona = 1000/1200).
- Rewrite the (very long) header doc comment for the merged scenario: one clear description of modes and the retained `ISO_` knobs (§B). Delete the historical staging narrative.
- Absorb gentle's drivers (§C): add the q(T)-flux ramp and coronal H(s) heating configuration + the new `ISO_QFLUX_*` / `ISO_CHEAT_*` / `ISO_TBOOST_*` knobs.
- `make_scenario`: dispatch `model_column`; add `model_isentropic` and `model_gentle` aliases (the latter sets the gentle preset). Delete the `model_gentle` case's body and delete `model_gentle.{hpp,cpp}` + its `#include`.
- Update `CMakeLists.txt` (drop `model_gentle.cpp`, add `model_column.cpp`).
- Fold `ISO_DIFF_MULT` to auto (§F).

### Phase 5 — Flip the default + downstream docs
- `chromo_main`: default scenario `model_c7` → `model_column`; update the usage comment.
- Update `tests/chromo_tests.cpp`: rename `model_isentropic_*` calls → `model_column_*`; drop now-defaulted `setenv("ISO_C7_IC"/"ISO_EQ_WB"/"ISO_LOG_RECON")` where they are now unconditional (keep the determinism test meaningful). Re-run; fix any assertions.
- Update `README.md` (CLI scenario list, examples), `visualize_commands.md` (migrate `model_isentropic`/`model_gentle` commands to `model_column`; the "bestwb" env soup collapses to a bare `model_column` run), and add the `GENTLE_*→ISO_*` migration table. Leave historical `docs/*_plan.md` as-is (they are a research record).
- Update `util/` parser defaults/labels that hardcode "model_c7" naming (`plot_output.py`, `animate_output.py`) to be scenario-agnostic or default to `model_column`. Note (do not necessarily fix here) the `H_BASE_KM=1003` hardcodes in `verify_pfss_flare.py` / `animate_pfss_flare.py` / `visualize_canopy.py` — `model_column` sets `out_base_km=0`, so those parsers must read the height column from line 2 (as `plot_output.py` already does) rather than assume 1003 km.

### Phase 6 — Verify
- `chromo_tests` passes with the same (or intentionally updated) count.
- The default `model_column` run is bit-identical to the old "bestwb" `model_isentropic` command.
- The `model_gentle` alias reproduces the old gentle configuration (equivalent physics).
- `/verify`-style end-to-end: run `./chromo_main` (default), a corona run, and a gentle-alias run; confirm expected evaporation signatures.

## Risks & mitigations

- **Silently changing actively-used runs.** Mitigated by keeping the reconstruction switches as Grid fields (defaults unchanged) so only `model_column` uses the new numerics; `model_c7`/`pfss`/`canopy`/`flare` and the regression baseline are byte-for-byte unchanged.
- **Broken saved commands.** Mitigated by `model_isentropic`/`model_gentle` aliases in `make_scenario` and the `GENTLE_*→ISO_*` migration table. (Old `ISO_*` numeric toggles like `ISO_LOG_RECON` become no-ops rather than errors — a removed `getenv` is simply not read; document this.)
- **Parser height mislabeling.** `model_column`'s `out_base_km=0` differs from C7's 1003 km; call out the 1003-hardcoded parsers (Phase 5) so plots of `model_column` output are not shifted.
- **Test churn from the rename.** Contained to the 3 iso tests + the alias test; enumerated in Phase 5.

## Out of scope (explicitly not doing now)

- Merging `model_flare` / `analytic_canopy` / `pfss_field_line`.
- Flipping global Grid reconstruction defaults or deleting the legacy minmod/flat-ghost solver branches.
- Renaming the `ISO_` env prefix to `COL_` (pure churn; deferred).
- Reworking the `ENABLE_TE` 7th-column handling in the generic plotters.
