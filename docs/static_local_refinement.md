# Static local-refinement support

A metric-aware, **static** non-uniform finite-volume mesh. This is static local
refinement, **not** AMR: the mesh is generated once at IC time and never
regridded, conservatively remapped, or adapted at runtime. Uniform meshes remain
the exact default — every legacy scenario and the test baseline are byte-for-byte
unchanged when refinement is disabled.

## What it does

Refinement coordinates are field-line distance from the domain **inner face**
(`s = 0` at the inner boundary), so the refined window is anchored the same way as
the photosphere-anchored `model_isentropic` C7 IC. The documented diagnostic
profile refines the lower domain **0–700 km by 4×**, retains the original outer
coarse spacing above 800 km, and transitions smoothly over **700–800 km** with an
adjacent cell-width ratio capped at **1.1**.

The face grid is built by `scenarios/mesh.{hpp,cpp}` from a coarse-equivalent cell
count `ns_coarse` (= `ISO_NS` for `model_isentropic`) and the refine parameters:

| segment | span | spacing |
|---|---|---|
| pre-fine coarse | `[0, s_lo]` | `coarse_ds = L/ns_coarse` (only if `s_lo > 0`) |
| fine | `[s_lo, s_hi]` | `≈ coarse_ds / factor` (uniform) |
| transition | `[s_hi, s_hi + τ]` | exponential (geometric) grade `fine → coarse`, ratio ≤ 1.1 |
| outer coarse | `[s_hi + τ, L]` | `≈ coarse_ds` (the retained outer spacing) |

The segment boundaries `0`, `s_lo`, `s_hi`, `s_hi + τ` and the domain top `L` are
exact faces. Refinement only **adds lower-domain cells**; `ISO_NS` stays the
coarse-grid-equivalent count and the rest of the column is never coarsened.

## Mesh controls (environment)

Generic (all static scenarios):

- `GRID_REFINE_FACTOR` — fine/coarse resolution ratio; **`1` disables** (default).
- `GRID_REFINE_S_LO_KM` — inner edge of the refined region (default `0`).
- `GRID_REFINE_S_HI_KM` — outer edge of the refined region (default `700`).
- `GRID_REFINE_TRANSITION_KM` — transition width (default `100`).

`model_isentropic` additionally accepts `ISO_REFINE_FACTOR`, `ISO_REFINE_S_LO_KM`,
`ISO_REFINE_S_HI_KM`, `ISO_REFINE_TRANSITION_KM`, which **override** the generic
`GRID_REFINE_*`. The documented diagnostic profile is `factor 4`, `s_lo 0`,
`s_hi 700`, `transition 100`.

**Scenario coverage.** The mesh builder, the Grid metric caches, and every
non-uniform-correct operator below are scenario-agnostic — any static scenario runs
correctly on a refined mesh. The refinement is currently *generated* by
`model_isentropic` (the diagnostic scenario), which builds its face grid through
the shared builder and `Grid::resize()`. PFSS field lines carry their own
externally-specified mesh and are never resampled. The other built-in scenarios
(`model_c7` / `model_flare` / `model_gentle` / `analytic_canopy`) build their fixed
uniform table-interpolated meshes as before; wiring them to consume `GRID_REFINE_*`
is a mechanical follow-up (route their interior face construction through
`build_static_mesh_faces` + `Grid::resize`, same pattern as `model_isentropic`).

## Correctness on non-uniform cells

Every spacing-sensitive operator has a non-uniform-correct path guarded so the
uniform mesh keeps its exact legacy code:

- **MUSCL reconstruction** limits center-to-center *gradients* and reconstructs
  each face with the owning cell's half-width (uniform ⇒ the legacy 0.5·Δ form).
- **Electron compression** `−p_e ∇·V` uses a metric-weighted face velocity.
- **Implicit conduction** (both paths — `rhs_implicit_state` and the Stage-D
  tridiagonal) uses width-weighted **series-resistance** face conductivities on a
  refined mesh; uniform meshes keep the arithmetic average.
- **TRAC** gradients, the radiative column integrals, the beam/coronal heating
  locations, diagnostics, and output coordinates all use the actual center
  separations / cell widths.
- The static **HSE ghost** constructions use the actual center-to-ghost distance
  (`ds_i(0)` at the fine base cell), so the discrete `V = 0` fixed point holds on
  the refined mesh (with `ISO_EQ_WB` it is exact to round-off; with the
  well-balanced reconstruction alone the residual base flow is m/s-scale).

The boundary model is otherwise **unchanged**: no free-outflow replacement and no
dynamic-HSE residual correction were introduced.

## Refined best-WB conduction-to-hydro diagnostic

The best well-balanced `model_isentropic` config (`ISO_INNER_WB` + `ISO_INNER_WB_RHO`
+ `ISO_EQ_WB` + `ISO_MC3` + `ISO_LOG_RECON`, γ = 1.05) holds `V = 0` to round-off,
so the only mass flux is physical (conduction-driven). A small **positive** ρV
maximum near ~250 km survives well-balancing in the conduction-ON runs. Refining
0–700 km by 4× lets us measure how that maximum scales with the **local** cell
width at 250 km.

Run the sweep (coarse-equivalent `ISO_NS = 2000, 3000, 4000`; writes one CSV per
run with the new `cell_width_m` column):

```bash
for NS in 2000 3000 4000; do
  ISO_C7_IC=1 ISO_C7_HSE=1 ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 ISO_HEAT_FLUX=1 \
    ISO_LOG_RECON=1 ISO_INNER_WB=1 ISO_INNER_WB_RHO=1 ISO_EQ_WB=1 ISO_MC3=1 ISO_MC3_BETA=2 ISO_GAMMA=1.05 \
    ISO_REFINE_FACTOR=4 ISO_REFINE_S_LO_KM=0 ISO_REFINE_S_HI_KM=700 ISO_REFINE_TRANSITION_KM=100 \
    ISO_NS=$NS ISO_DIAG_COND_HYDRO=1 ISO_DIAG_COND_HYDRO_EVERY=$((NS/4)) \
    ISO_DIAG_COND_HYDRO_OUT=outputs/iso_cond_hydro_local4_outer$NS.csv \
    ./build/chromo_main /tmp/_iso_local4_$NS.txt full no-ionization model_isentropic - 0.18 no-cooling
done
```

Verify each CSV has full snapshots bracketing ~100 s, then plot (the fit is against
the **interpolated local cell width at 250 km**, not the median domain spacing):

```bash
.venv/bin/python util/plot_cond_hydro_diagnostic.py --time 100 \
  --out visualization/iso_cond_hydro_local4_convergence_t100.png \
  outputs/iso_cond_hydro_local4_outer2000.csv \
  outputs/iso_cond_hydro_local4_outer3000.csv \
  outputs/iso_cond_hydro_local4_outer4000.csv
```

## Split conduction-to-hydro coupling artifact vs pure conduction-operator artifact

Two different numerical effects must not be conflated:

- **Split conduction-to-hydro coupling artifact.** The ~250 km positive ρV maximum
  is *conduction-gated* (absent when conduction is off) yet **invariant under
  well-balancing** — well-balancing corrects the explicit hydrostatic/advective
  operator (a static, conduction-off property), whereas this artifact is exposed
  when the *next* explicit hydro stage differentiates the pressure field that the
  *separate*, operator-split conduction stage just modified. It is a property of
  the **coupling between two split operators**, and it scales with the local cell
  width, so it shrinks under local refinement of the base.

- **Pure conduction-operator artifact.** An error internal to the conduction
  operator alone (e.g. the face-conductivity averaging or the tridiagonal solve on
  an under-resolved gradient) would appear in the conduction stage in isolation and
  is not what the split coupling above measures.

The **numerical remedy** for the split-coupling artifact is grid / local refinement
of the lower chromosphere (this feature) or an improved conduction–hydro coupling.
Radiative cooling is a **physical** sensitivity test, not a numerical fix, and is
deliberately left out of the numerical-remedy discussion. The physical A&S78 TR
evaporation upflow (~+0.9 km/s) is ~100× larger and grid-converged — unaffected by
either artifact.

### Analyzed scaling (`visualization/iso_cond_hydro_local4_convergence_t100.png`)

Refining 0–700 km by 4× drops the **local** cell width at 250 km to 0.269 / 0.179 /
0.135 km for coarse-equivalent `ISO_NS` = 2000 / 3000 / 4000 (the fine spacing,
`coarse_ds/4`), while the outer column is identical across the three runs. At
t ≈ 100 s the total mass flux ρ_iV + ρ_nU at 250 km is
6.37 / 4.50 / 3.24 ×10⁻⁴ kg m⁻² s⁻¹ — a clean linear scaling in the local Δs
(R² = 0.994) with a Δs → 0 intercept of 2.5×10⁻⁵ (≈ 0, ~4 % of the coarsest value).
So the ~250 km mass flux is a **local-Δs artifact that local base refinement
suppresses**, confirming the split conduction-to-hydro coupling reading. (The
per-step conduction pressure change dp_cond/dt at 250 km rises on finer grids as it
resolves the local conduction gradient better — a physical resolution effect, R²≈0.74
— while the total momentum-RHS difference stays at the ~10⁻⁷ noise floor, R²≈0.01.)
