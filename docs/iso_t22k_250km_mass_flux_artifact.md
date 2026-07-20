# Numerical mass-flux artifact near 200–250 km in the `iso_t22k` run

## Purpose

This note summarizes the small positive mass-flux feature near 200–250 km in the run used to generate

`visualization/model_column/iso_t22k_evolution_ns2000_gamma105_bestwb.mp4`

from

`outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt`.

The feature is **not interpreted as a physical lower-chromospheric upflow**. Resolution studies show that it is a finite-grid numerical artifact associated with the operator-split coupling between thermal conduction and the following hydrodynamic update.

---

## 1. Run configuration

The original run command was

```bash
ISO_C7_IC=1 \
ISO_C7_HSE=1 \
ISO_H_BASE=0 \
ISO_DH=2152.6 \
ISO_T_TOP=22000 \
ISO_HEAT_FLUX=1 \
ISO_LOG_RECON=1 \
ISO_INNER_WB=1 \
ISO_INNER_WB_RHO=1 \
ISO_EQ_WB=1 \
ISO_MC3=1 \
ISO_MC3_BETA=2 \
ISO_NS=2000 \
ISO_GAMMA=1.05 \
./build/chromo_main \
  outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt \
  full no-ionization model_isentropic - 120 no-cooling
```

The animation was decoded using

```bash
ISO_GAMMA=1.05 \
python util/animate_isentropic.py \
  outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt \
  visualization/model_column/iso_t22k_evolution_ns2000_gamma105_bestwb.mp4
```

### Main numerical and physical settings

- One-dimensional, field-aligned vertical column.
- Domain: approximately 0–2152.6 km above the photospheric base.
- Uniform grid with 2000 cells.
- Cell width: approximately 1.076 km.
- Constant solar surface gravity.
- Nearly isothermal equation of state with `gamma = 1.05`.
- Single-fluid evolution with a frozen height-dependent ionization fraction.
- Thermal conduction enabled.
- Time-dependent ionization disabled.
- Radiative cooling disabled.
- No beam heating or volumetric coronal heating.
- Log-space reconstruction of density and pressure.
- MC3/Koren-type limiter with beta = 2.
- Full equilibrium-reference well-balancing enabled.

The old environment-variable names are retained above because they identify the exact historical run. The present `model_column` implementation has since absorbed this setup and enables most of the validated well-balanced stack by default.

---

## 2. Initial condition

The initial atmosphere was based on Model C7, but the C7 tabulated density was not inserted directly.

### Quantities retained from C7

At every height, the initialization sampled

- the C7 temperature profile, `T(h)`;
- the C7 local ionization fraction, obtained from the electron and neutral-hydrogen number densities.

### Hydrostatic reconstruction of density and pressure

The total pressure was re-integrated upward from the photospheric base using hydrostatic equilibrium,

```text
dp/ds = -rho g.
```

For hydrogen with local ionization fraction `x`, the total pressure is

```text
p = (2 n_i + n_n) k_B T
  = (1 + x) n_total k_B T.
```

The code therefore integrated `ln(p)` using the local temperature and ionization-dependent scale height, then reconstructed the total number density from the equation of state. The ion and neutral densities were finally obtained by repartitioning the total density according to the C7 ionization fraction.

The initial velocities were zero everywhere. Thus the intended initial state was a stationary, non-isothermal hydrostatic atmosphere that retained the C7 thermal and ionization structure while using a self-consistent hydrostatic density profile.

### Well-balanced reference state

The full initial conservative state was stored as an equilibrium reference. On the first explicit step, the code evaluated the discrete hydro RHS of this state and cached the residual. Subsequent explicit hydro updates used

```text
R_corrected(U) = R_hydro(U) - R_hydro(U_equilibrium).
```

Consequently, with conduction disabled, the initial atmosphere is an exact fixed point of the discrete explicit hydro operator up to floating-point and implementation-level residuals. This correction does not make a state modified later by conduction hydrostatic.

---

## 3. Boundary conditions used in this run

## 3.1 Lower boundary: photospheric reservoir

The lower hydrodynamic boundary was a fixed, stationary photospheric reservoir.

- Velocity in both lower ghost cells was set to zero.
- Ghost-cell pressures continued the hydrostatic column below the first interior cell.
- The first and second ghost pressures were approximately

```text
p_g1 = p_0 + rho_0 g Delta_s,
p_g2 = p_0 + 2 rho_0 g Delta_s.
```

- At the fixed base temperature, the ghost densities were scaled with pressure, so that `rho_g proportional to p_g`.
- This two-ghost construction made the lower MUSCL stencil consistent with a stationary hydrostatic reservoir and removed the much larger numerical base drainage seen with simpler boundary treatments.

For thermal conduction, `ISO_INNER_T_NEUMANN` was not specified. The run therefore used the default **Dirichlet-type lower temperature boundary**, in which the conduction solve couples the first interior cell to a fixed lower ghost temperature. It was not an insulating zero-flux boundary.

## 3.2 Upper boundary: hot reservoir at the top of the modeled column

The top of the domain represented the lower edge of an unresolved hotter region.

- `ISO_T_TOP=22000` imposed a fixed reference top temperature of 22,000 K for the outer thermal boundary.
- The default temperature-jump factors were unity, so the outer ghost temperatures were effectively held at 22,000 K.
- Ghost pressure was hydrostatically extrapolated from the upper interior cells.
- Ghost density was then reconstructed from the imposed ghost pressure, ghost temperature, and the local top-cell ionization fraction.
- Because conduction was enabled, the upper velocity boundary allowed an outflow copied from the top interior velocity, limited to a small fraction of the sound speed. This allowed the upper evaporation flow to leave the domain while suppressing unstable inflow.

This fixed-temperature upper reservoir drove a downward conductive heat flux into the modeled transition-region/chromospheric column.

---

## 4. Observed numerical feature

The conduction-on run develops a small positive mass-flux maximum in the lower chromosphere, broadly around 200–250 km. The exact peak location changes with resolution and can lie somewhat below 250 km.

This feature persists after applying the complete hydro well-balancing stack:

- corrected gravitational-potential handling during shifted-state reconstruction;
- log-space density and pressure reconstruction;
- MC3 limiter;
- hydrostatic two-ghost lower boundary;
- equilibrium-reference residual subtraction.

However, the corresponding conduction-off runs do not show the same feature at comparable amplitude. The feature is therefore conduction-gated.

---

## 5. Cause of the artifact

The artifact is caused by the **operator-split coupling between the implicit conduction stage and the next explicit hydro stage**.

The timestep ordering was approximately

```text
1. explicit MUSCL/Rusanov hydrodynamics,
2. ion-neutral drag,
3. temperature equilibration,
4. implicit thermal conduction,
5. optional ionization/heating/cooling stages.
```

The important point is that the conduction stage changes temperature and thermal pressure while, during that stage, leaving density and velocity fixed. Immediately after conduction, the state generally no longer satisfies

```text
dp/ds + rho g = 0.
```

On the following timestep, the explicit hydro operator differentiates this post-conduction pressure field. The finite-volume pressure gradient and gravity source do not cancel exactly for the conduction-modified, curved, strongly stratified profile. The resulting small residual momentum source produces the positive mass-flux feature near 200–250 km.

This is why the feature is not removed by static well-balancing. The equilibrium-reference correction balances the explicit hydro operator around the original stationary reference state. It does not rebalance the new pressure field created after every separate conduction solve.

Schematically,

```text
initial HSE state
    -> hydro update is well balanced
    -> conduction changes p at fixed rho and V
    -> state is no longer discretely hydrostatic
    -> next hydro update differentiates the modified p profile
    -> finite-grid residual acceleration
    -> small positive mass flux near 200–250 km.
```

The evidence indicates that this should be described as a **split conduction-to-hydro coupling artifact**, rather than as an error generated by the conduction operator in isolation. A pure conduction-operator error would already be identifiable within the conduction stage alone. Here the spurious flow appears when the following hydro stage responds to the pressure field produced by conduction.

---

## 6. Grid-convergence evidence

A matched set of conduction-on runs was performed with the same physical and numerical configuration but with 2000, 3000, and 4000 uniform cells.

The corresponding grid spacings were approximately

| Number of cells | Delta s near 250 km |
|---:|---:|
| 2000 | 1.08 km |
| 3000 | 0.72 km |
| 4000 | 0.54 km |

At approximately 100 s, the velocity interpolated near 250 km was

| Number of cells | V near 250 km |
|---:|---:|
| 2000 | 10.7 m/s |
| 3000 | 7.1 m/s |
| 4000 | 5.4 m/s |

A linear fit against `Delta s` extrapolated to approximately 0.07 m/s as `Delta s -> 0`, consistent with zero. The positive lower-chromospheric mass-flux peak likewise decreased approximately linearly with resolution:

| Number of cells | peak rho V in the lower chromosphere |
|---:|---:|
| 2000 | 8.5 x 10^-4 kg m^-2 s^-1 |
| 3000 | 5.8 x 10^-4 kg m^-2 s^-1 |
| 4000 | 4.3 x 10^-4 kg m^-2 s^-1 |

The extrapolated zero-grid-spacing intercept was small, and the peak location moved toward the lower boundary as the grid was refined. These are signatures of a finite-grid artifact rather than a converged physical flow.

A second test refined only the lower 0–700 km region by a factor of four while leaving the outer column effectively unchanged. At 250 km, the local cell widths became approximately 0.269, 0.179, and 0.135 km. The mass flux at about 100 s was

| Local Delta s at 250 km | rho V at 250 km |
|---:|---:|
| 0.269 km | 6.37 x 10^-4 kg m^-2 s^-1 |
| 0.179 km | 4.50 x 10^-4 kg m^-2 s^-1 |
| 0.135 km | 3.24 x 10^-4 kg m^-2 s^-1 |

The scaling with the **local** cell width was nearly linear, with `R^2 = 0.994`, and the extrapolated intercept was approximately `2.5 x 10^-5 kg m^-2 s^-1`, close to zero relative to the coarsest value. This local-refinement test isolates the lower-chromospheric spatial resolution as the controlling numerical parameter.

---

## 7. Distinction from the physical evaporation flow

The small lower-chromospheric feature should not be confused with the transition-region evaporation flow near the top of the domain.

- The 200–250 km feature is conduction-gated but decreases approximately as `O(Delta s)` and extrapolates toward zero.
- The upper transition-region evaporation reaches roughly `+0.9 km/s`, is about two orders of magnitude larger in velocity, and is grid-converged in the tested resolutions.

The first is a numerical coupling artifact; the second is the intended physical response to the imposed downward conductive heat flux.

---

## 8. Practical interpretation and remedies

The present evidence supports the following interpretation:

> The positive mass flux near 200–250 km is generated when the explicit hydro step acts on a pressure profile modified by the preceding, separately solved conduction step. The residual scales with the local grid spacing and vanishes under refinement, so it is a finite-grid operator-splitting artifact rather than a physical lower-chromospheric upflow.

The available numerical remedies are

1. refine the lower chromosphere, either globally or with static local refinement;
2. improve the conduction-hydro coupling, for example through a more tightly coupled or well-balanced treatment of the post-conduction pressure update;
3. consider higher-order or symmetric operator splitting as a diagnostic or future improvement, while verifying that the modified method preserves the intended energy transfer.

Changing radiative cooling is a physical sensitivity test, not a numerical correction for this artifact. Likewise, ordinary hydro well-balancing alone cannot remove it because the source is the interaction between two separately applied operators.

---

## 9. Concise summary for discussion

- The run begins from a zero-velocity C7-based atmosphere whose pressure and total density were re-integrated into hydrostatic equilibrium.
- The lower boundary is a fixed, zero-velocity hydrostatic photospheric reservoir; the lower conduction boundary is fixed-temperature Dirichlet.
- The upper boundary imposes a 22,000 K thermal reservoir with hydrostatic ghost pressure and a Mach-capped outflow.
- The explicit hydro solver is strongly well balanced and holds the conduction-off reference atmosphere close to stationary.
- Conduction changes pressure at fixed density and velocity.
- The following explicit hydro step differentiates this newly modified pressure field, leaving a finite-grid pressure-gravity residual near 200–250 km.
- The resulting positive mass flux scales nearly linearly with the local cell width and extrapolates to approximately zero.
- It is therefore a numerical split conduction-to-hydro coupling artifact, not a physical flow.
