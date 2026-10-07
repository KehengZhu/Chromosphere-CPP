# Visualization catalog — Saha closure and the CRASH `Gamma1` table

The pure-hydrogen Saha closure and the production `Gamma1` table the release reads. Outputs land in `visualization/eos_gamma/`. See [README.md](README.md) for conventions.

---

## Classical pure-H Saha / CRASH Γ₁ table (Stages 0–2)

`util/eos/tabulate_gamma.f90` builds the approved classical pure-hydrogen
cross-validation table with excitation, Fermi, and Coulomb corrections disabled.
The raw 11-column table retains `zAv`, energy γ_E, Γ₁, energy, and `Cv` for
validation; the tracked production table at
`data/eos/gamma1_hydrogen_v1.dat` contains only
`log10(T), log10(n_H), Γ₁`, with strict v1 metadata and a SHA-256 sidecar.
The grid is `T=3.2×10³–10⁸ K` (501 log points) by
`n_H=10¹²–10²⁶ m⁻³` (57 log points).

```bash
# Generate outputs/eos_gamma/gamma_hydrogen.dat and refresh the tracked
# data/eos/gamma1_hydrogen_v1.dat, verify its checksum, then run the strict
# whole-grid Saha/energy/Cv/Gamma1 validation.
bash util/eos/build_and_run.sh

# Generate all three figures below (project root as working directory).
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/plot_eos_gamma.py
```

### `eos_gamma/gamma_hydrogen.png` — γ_E, Γ₁, and ionization vs T

Energy γ_E (solid), sound-speed Γ₁ (dashed), and equilibrium ionization
for four densities spanning `10¹⁴–10²⁶ m⁻³`. Both indices recover to 5/3
at the neutral and fully ionized limits and soften through the Saha transition.

### `eos_gamma/gamma_hydrogen_2d.png` — production Γ₁(T,n_H) surface

Heatmap of the sound-speed index over the entire runtime grid with
`x={0.01,0.1,0.5,0.9,0.99}` ionization contours. This is the column loaded by
the Stage-2 C++ interpolator; CRASH pressure and energy are not runtime inputs.

### `eos_gamma/c7_eos_consistency.png` — C7 track and Stage-0 consistency

Walks the complete Model C7 atmosphere parsed from `scenarios/model_c7.cpp` and
shows C7 ionization against textbook Saha, CRASH γ_E and Γ₁, and the required
central-`ln T` `C_V^eff=(∂e_int/∂T)_rho` against CRASH `Cv`. The chromosphere/TR
is expanded and the extended corona compressed with a split height scale.

### `eos_gamma/c7_gamma_profile_classical_saha.png` — signed-off classical EOS along Model C7

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/plot_c7_gamma_profile.py --classical-saha
```

Preserves the earlier `c7_gamma_profile.png` and generates a separately named
height profile for the signed-off classical pure-H EOS. It overlays the
textbook-Saha energy index `γ_E=1+P/e` and the CRASH sound-speed index `Γ₁`
with excitation, Fermi, and Coulomb corrections disabled, alongside the same
Model C7 temperature and density context used by the original figure.

### `eos_gamma/c7_ionization_profile.png` — Model C7 ionization degree

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/plot_c7_ionization_profile.py
```

Model C7 hydrogen ionization degree `x=n_e/(n_HI+n_e)` vs height, alongside
the same temperature and mass-density context as `c7_gamma_profile.png`.

### `eos_gamma/c7_ionization_profile_classical_saha.png` — classical-Saha ionization degree

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/plot_c7_ionization_profile.py
```

Classical pure-H Saha ionization degree evaluated at the local Model C7
temperature and hydrogen-nuclei density, plotted separately from the C7 curve.

### `eos_gamma/c7_ionization_profile_compare.png` — C7 vs classical-Saha vs excitation-enabled CRASH

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/plot_c7_ionization_profile.py --mode compare
```

Direct overlay of the Model C7 ionization degree, analytic equilibrium
classical-Saha fraction, and the excitation-enabled CRASH code's tabulated
`zAv` evaluated along the same C7 atmosphere. An inset magnifies the
2119–2147 km transition-region interval where C7 and CRASH briefly diverge by
more than 20%; a small rectangle marks that interval on the main panel and an
arrow links it to the inset. The CRASH diagnostic table is generated separately,
preserving the signed-off excitation-off production table:

```bash
bash util/eos/build_and_run.sh excitation
```

The default plotting command writes all three ionization-profile figures.

### `eos_gamma/crash_excitation_full_table_c7.png` — full CRASH excitation comparison with C7 track

```bash
MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python util/release/plot_crash_excitation_full_table.py
```

Three-panel view of the complete CRASH `(T,n_H)` grid: equilibrium hydrogen
ionization with bound excitation off, with excitation on, and the percentage
reduction caused by excitation. The complete Model C7 thermodynamic track is
overlaid on every panel and labeled at the photosphere, temperature minimum,
upper chromosphere, TR base, and upper TR. Contours identify the 1%, 10%, and
30% excitation-sensitive regions.

---
