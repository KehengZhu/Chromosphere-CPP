# Diagnosis: the broad low-chromosphere oscillation in the Roe + `(ln ρ, V, ln p)` N1000 run

> **Partly superseded (see `docs/reference_free_release_recap.md`).** The release hydrodynamic operator has since gained the MUSCL-Hancock predictor momentum source and a trapezoidal EOS-closed lower ghost ladder, and **equilibrium-reference well balancing (`eq_wb` / `ISO_EQ_WB`) has been retired from the release** — it survives only in the two-fluid research solver. Statements below that describe `eq_wb` as active release method, or that quantify the discrete hydrostatic defect, are historical. Everything else in this document still stands.


**Scope.** Bounded read-only diagnosis of the large-scale velocity / cell-centred mass-flux oscillation seen in `outputs/model_column/lnp_roe_N1000_4000s.txt`. No production code, defaults, or release decision were changed. One new bounded control run and two diagnostic scripts were added.

**Verdict up front.** The feature is a **conduction-excited, boundary-trapped standing oscillation of the sub-TR column**, present with the same amplitude under Rusanov as under Roe, and damped only numerically. It is bounded, smooth, resolved, and does **not** contaminate the evaporation flux. **Release impact: NON-BLOCKER.**

---

## Evidence base

| Run | Configuration | Role |
| --- | --- | --- |
| `lnp_roe_N1000_4000s` | Roe + lnP, N1000/R4, 1320 cells, 4000 s | primary |
| `lnp_roe_N500_4000s_opt` | Roe + lnP, N500/R4, 661 cells, 4000 s | resolution control |
| `lnp_roe_N1000_partial2080s` `.faceflux` | same as primary, interrupted at 2080 s | conservative flux, 5.5 s cadence |
| `lnp_roe_N1000_1000s` `.faceflux` | same, 1000 s | early transient, 2.8 s cadence |
| `acc_lnp_N1000` / `acc_lnp_N500` | **Rusanov** + lnP, 4000 s | solver control |
| `acc_lnt_N1000` | Rusanov + lnT, 4000 s | reconstruction control |
| `_archive/.../eqctl_lnp_N500_1000s` | Roe + lnP, `ISO_HEAT_FLUX=0` | conduction-off control (existing) |
| `eqctl_lnp_roe_N1000_1000s` | Roe + lnP, N1000, `ISO_HEAT_FLUX=0` | **new**, conduction-off control at N1000 |

New control command (only new simulation run for this task):

```bash
ISO_NS=1000 ISO_HEAT_FLUX=0 CHROMO_T_END=1000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FRAME_DT=20 \
  scripts/run_chromo_realtime.sh outputs/model_column/eqctl_lnp_roe_N1000_1000s.txt \
  full no-ionization model_column - 20.0 no-cooling
```

Console confirms `conduction.calls=0` and otherwise release defaults (`flux=roe-local reconstruction=lnrho-v-lnp limiter=mc3 beta=2`).

New tooling: `util/column_spacetime.py` (cached `(n_t, n_cell)` loader + height/phase primitives for `.gamma_diag`), `util/plot_column_sloshing.py` (space-time / amplitude / evaporation figure). Figure: `visualization/model_column/roe_n1000_lower_chromosphere_sloshing.png`.

Domain reminder: `model_column` spans **1600.55–2152.99 km** (552 km). The "low chromosphere" region 1739–1847 km is therefore ~140–250 km above the inner boundary, in the middle of the domain.

---

## 1. What the animation is actually showing

`util/animate_isentropic.py` plots `mflux = rho * V` from the `.gamma_diag` sidecar — a **cell-centred** diagnostic, not the finite-volume face flux that updates mass. Two separate effects inflate its visual impact:

- **The mass-flux panel is the only panel without a fixed y-limit.** `ylims["mflux"]` is computed (line 249) but the `set_ylim` call is commented out (line 170), while `T`, `V`, `q`, `ρ`, `p` all use fixed global limits. The panel therefore rescales every frame.
- A zero line is drawn (`axhline(0)`). N500's domain-wide `ρV` minimum stays positive at all times (+2.6e-10 … +4.5e-10 kg m⁻² s⁻¹), so its curve never touches that line; N1000 dips to −1.6e-10 at t ≈ 1000 s and −1.5e-10 at t ≈ 3000 s and visibly crosses it. The domain-wide `ρV` **span** is comparable in both runs (0.9–1.6e-9), so the difference the eye reads as "much larger variation" is largely the sign crossing plus per-frame rescaling.

**But the sign is real.** Comparing the cell-centred quantity against the production conservative face flux `f_eff = f_total − f_ref` (what the `eq_wb` δ-form update actually differences), Roe + lnP N1000, window 1700–1900 km:

| t [s] | ⟨ρV⟩ cell | min ρV | ⟨f_eff⟩ face | min f_eff |
| --- | --- | --- | --- | --- |
| 393 | +1.318e-09 | +9.69e-10 | +1.322e-09 | +9.73e-10 |
| 988 | +1.833e-11 | **−1.236e-10** | +2.175e-11 | **−1.202e-10** |
| 1993 | +1.118e-09 | +7.61e-10 | +1.122e-09 | +7.64e-10 |

The two agree to a few percent in mean, minimum and sign. **Unlike the previously diagnosed upper-TR ripple, here the cell-centred diagnostic is a faithful proxy**: the negative excursion is genuine (very small) conservative downward transport, not a decomposition artifact. Roughness (normalised second difference) is 4.0e-3 for cell `ρV` and 4.8e-3 for `f_eff` — both smooth.

---

## 2. Spatial character — resolved, not a grid ripple

At t = 1000 s (N1000) the negative-V region spans **1738.8–1846.6 km: 108 km, 196 cells**; at t = 3000 s, 1733.2–1849.9 km (212 cells).

Decomposing V in 1620–2050 km into components above/below 30 km:

| t [s] | broad (>30 km) p-t-p [m/s] | fine (<30 km) rms [m/s] | fine/broad |
| --- | --- | --- | --- |
| 100 | 9.51 | 0.0157 | 0.0016 |
| 1000 | 3.05 | 0.0181 | 0.0059 |
| 2000 | 1.48 | 0.0095 | 0.0064 |
| 4000 | 1.45 | 0.0177 | 0.0123 |

After the initial transient the cell-scale content is **0.5–1.2 % of the broad-scale amplitude**. Normalised curvature of `ln p` in 1620–2100 km is 9e-6 (N1000) and 2e-5 (N500) — the thermodynamic fields are smooth across the feature; there is no cell-to-cell sign alternation of significant amplitude.

Thermodynamic perturbations in the oscillating region are negligible: over the whole 4000 s, |Δp/p|, |Δρ/ρ|, |ΔT/T| at 1650/1800/1900/2000/2100 km all stay **below 1.2e-3**, and most below 5e-4. The mode is essentially a velocity-only, Mach ~1e-4 disturbance riding on a thermodynamically static column. All the real restructuring is in the top ~15 km (at 2140 km: Δρ/ρ = +21 %, ΔT/T = −8.5 %).

**Conclusion: a resolved O(100 km) / O(200 cell) structure, categorically different from a cell-scale numerical ripple.**

---

## 3. Temporal character — standing, non-acoustic period, acoustic launch

**Period.** Spectra of the 5.5 s-cadence series (t ≥ 150 s, 0–2080 s) peak at ~1.9e3 s at every height in 1650–2000 km, with no significant power at the acoustic cavity period. Zero-crossing estimates of the polynomial-detrended V(1800 km, t) over 4000 s give P ≈ 1.65e3 s (Roe N1000), 1.28e3 s (Roe N500), 1.34e3 s (Rusanov N1000), 0.93e3 s (Rusanov N500). Only ~2–3 cycles are available, so **the period is O(10³ s) but is not resolution- or scheme-converged.**

**Standing, not propagating.** High-pass-filtered V (t ≥ 300 s, N1000) correlated against V(1800 km):

| h [km] | 1620 | 1650 | 1700 | 1750 | 1800 | 1850 | 1900 | 1950 | 2000 | 2050 | 2100 | 2140 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rms [m/s] | 0.14 | 0.30 | 0.54 | 0.69 | 0.70 | 0.72 | 0.66 | 0.56 | 0.51 | 0.43 | 0.46 | 0.52 |
| corr @ lag 0 | +0.82 | +0.91 | +0.97 | +0.98 | +1.00 | +0.98 | +0.96 | +0.88 | +0.76 | +0.57 | +0.50 | +0.51 |
| best lag [s] | +40 | +20 | +20 | 0 | 0 | 0 | −20 | −60 | −60 | −80 | −80 | −80 |

Lags are ≤ 80 s against a ~1.7e3 s period (< 5 % of a cycle, at/near the 20 s sampling): the column oscillates **in phase**. Amplitude has a **near-node at the inner wall** (0.07 m/s at 1600 km) and peaks at 1830 km. Pressure leads velocity by 500–640 s ≈ a quarter period (p@2000 km vs V@1800 km: corr −0.50 at zero lag, +0.41 at −640 s) — the p–V quadrature of a genuine oscillator, not a forced in-phase response.

**The period is not acoustic.** Using the t = 0 state (c_s = 8.7–25 km/s, H_p ≈ 280 km):

- one-way base→top acoustic time ∫ds/c_s = **58.7 s**; round trip **117 s**;
- acoustic cutoff period 4πH/c_s ≈ **380–400 s** — the observed ~1.7e3 s oscillation is far below the cutoff frequency, i.e. non-propagating, consistent with the in-phase standing pattern;
- conductive diffusion across the column: ~1e5–1e6 s; across the hot top layer (10.5 km, χ = 4.6e5 m²/s): ~240 s.

None of these matches ~1.7e3 s. **The restoring mechanism that sets the period is not identified by this bounded diagnosis** (it necessarily involves the coupling between the conductive top boundary layer and the column, since the mode does not exist with conduction off — §5).

**The excitation is unambiguously acoustic.** The 2.8 s-cadence `.faceflux` capture resolves the conduction turn-on directly (figure, top-right panel):

| height [km] | 2100 | 2000 | 1900 | 1800 | 1700 | 1650 | 1605 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| front arrival [s] | ~8 | ~16 | ~27 | ~37 | ~52 | ~53 | ~58 |

Base arrival at ~58 s matches ∫ds/c_s = 58.7 s to better than 2 %.

---

## 4. N500 vs N1000 — same mode, different damping

Amplitude of high-passed V at 1800 km (t ≥ 300 s):

| run | rms [m/s] | mean V (t > 2000 s) |
| --- | --- | --- |
| Roe + lnP N1000 | 0.699 | +0.90 |
| Roe + lnP N500 | 0.305 | +1.08 |
| Rusanov + lnP N1000 | 0.729 | +1.17 |
| Rusanov + lnP N500 | 0.395 | +1.00 |

The **time-mean** transport agrees across all four runs to ~20 %; only the oscillation amplitude differs, by ~2× per halving of Δs. The amplitude-vs-height profiles (figure, bottom-middle) have identical shape and collapse onto one another above ~2050 km.

RMS amplitude in successive windows at 1800 km:

| run | 300–1550 s | 1550–2800 s | 2800–4000 s |
| --- | --- | --- | --- |
| Roe + lnP N1000 | 1.36 | 0.52 | 0.56 |
| Roe + lnP N500 | 0.71 | 0.16 | 0.08 |
| Rusanov + lnP N1000 | 1.40 | 0.91 | 0.49 |
| Rusanov + lnP N500 | 0.76 | 0.46 | 0.33 |

N500 decays with an e-folding time of order 1.1e3 s; N1000 is **nearly undamped after the first cycle** and is bounded, not growing.

The **excitation itself is resolution-independent**: velocity in the top cell at t = 20 s is 40.75 m/s (Roe N1000) vs 40.79 m/s (Roe N500), and 54.21 vs 53.08 (Rusanov). So the amplitude difference at 1800 km is **damping, not excitation**.

This model configuration contains **no physical damping for this mode**: radiative cooling is off (`ISO_COOLING=0`), there is no viscosity, and conduction acts on 1e5 s. The only dissipation is numerical. The N1000 result is therefore the more faithful solution *of the model as posed*; the N500 decay is a numerical artifact.

---

## 5. Excitation source — conduction, not the hydrostatic residual

New N1000 conduction-off control (`ISO_HEAT_FLUX=0`, all other release defaults, `conduction.calls=0`), after 1000 s:

| quantity | conduction OFF (N1000) | conduction ON (N1000) |
| --- | --- | --- |
| V at 1800 km | −2.7e-3 m/s | −0.14 m/s (and ±3 m/s over a cycle) |
| max &#124;V&#124; in 1620–2100 km | **0.098 m/s** | ~3 m/s |

The archived N500 control (`eqctl_lnp_N500_1000s`) agrees: max |V| ~ 5e-2 m/s, V(1800 km) = +7e-3 m/s.

**The equilibrium / well-balanced residual accounts for ≲ 3 % of the mode amplitude at N1000 and ≲ 0.3 % at 1800 km.** Hypothesis C (well-balanced defect) is ruled out at both resolutions. The `eq_wb` reference is verified frozen: the captured `eq_residual_mass` is bit-identical across all output times and matches −Δf_ref/Δs to 0.7 %.

---

## 6. Is Roe uniquely responsible? No.

| run | rms of high-passed V @ 1800 km |
| --- | --- |
| Roe + lnP N1000 | 0.699 |
| Rusanov + lnP N1000 | 0.729 |
| Rusanov + lnT N1000 | 0.783 |

Rusanov shows the mode at **slightly larger** amplitude than Roe at the same resolution, with the same ~2× resolution scaling and the same spatial profile. Reconstruction (lnP vs lnT) is likewise irrelevant. Phase and period differ by ~15 % between schemes, as expected when numerical dispersion is comparable, but the mode is present and of the same character in all of them. **Hypothesis D (Roe-specific low-dissipation numerical mode) is ruled out.** No solver-selection question is reopened.

---

## 7. Lower-boundary role — reflecting, and it traps the mode

`scenarios/model_column.cpp:840–866`: both inner ghosts are **fully frozen Dirichlet** — density, pressure and temperature captured from the IC (`kInnerRhoIGh`, `kInnerPiGh`, `kInnerTRef`, and the second-ghost equivalents), and `MOM_I = MOM_N = 0`. No interior state enters them. There is no characteristic decomposition, no outgoing-wave condition: **the boundary is perfectly reflecting for outgoing acoustic characteristics.**

Two independent lines of evidence that it reflects the mode:

1. **Direct capture of the reflection.** In the 2.8 s-cadence record, V at 1605 km goes to −0.045 m/s (t = 56 s), −0.35 m/s (t = 62 s), then rebounds to +0.23 (t = 67 s) and +0.46 (t = 79 s); at 1650 km, −0.25 → +0.22 → +0.66 over t = 56 → 67 s. The incident compression arrives, inverts at the wall, and travels back up.
2. **Node structure.** The mode amplitude profile has a near-node at the wall (0.07 m/s at 1600 km against 0.72 m/s at 1830 km) — the signature of a standing mode against a velocity-node boundary.

Because the model contains no physical damping for this mode and the floor is perfectly reflecting, an excited oscillation rings for the whole run. The boundary is not the *source*; it is what **traps** the conduction-launched disturbance.

*Characterisation note (not a defect):* the inner ghosts hold V = 0 but the inner-face Riemann problem still admits mass flux, so the inner boundary acts as a **mass reservoir**. Over 4000 s the column loses only 0.018 % of its mass while ~0.9 % of its mass leaves through the top — most of the evaporated mass is resupplied from the photospheric reservoir. This is the documented intent of the discrete-HSE reservoir BC.

---

## 8. Effect on evaporation — none material

Window 2130–2150 km, cell-centred ⟨ρV⟩ over t ≥ 1000 s:

| run | mean [kg m⁻² s⁻¹] | min | max | high-pass rms / mean | min V [m/s] |
| --- | --- | --- | --- | --- | --- |
| Roe + lnP N1000 | +3.855e-10 | +3.611e-10 | +4.494e-10 | **2.07 %** | +2.99 |
| Roe + lnP N500 | +3.761e-10 | +3.551e-10 | +4.271e-10 | **1.98 %** | +2.96 |
| Rusanov + lnP N1000 | +3.733e-10 | +3.451e-10 | +4.456e-10 | 2.01 % | +2.84 |
| Rusanov + lnP N500 | +3.464e-10 | +3.161e-10 | +4.146e-10 | 1.91 % | +1.65 |

- The mean evaporation flux is **resolution-consistent to 2.5 %** (Roe N1000 vs N500).
- Velocity and mass flux in the evaporation window are **positive at all times** in every run; there is no reversal.
- The oscillatory modulation of the evaporation flux is **2.07 % at N1000 vs 1.98 % at N500** — i.e. *the resolution difference that doubles the low-layer mode does not measurably change the evaporation-flux modulation.*

Conservative face flux at 2140 km (Roe + lnP N1000): +5.78e-10 (t = 393 s), +4.27e-10 (988 s), +3.83e-10 (1993 s) — uniformly positive, monotone, roughness ~5e-4, with min f_eff/max f_eff = 0.98 across 2130–2150 km at every sampled time.

The mode is **confined below ~2000 km**: above 2050 km the amplitude-vs-height curves of all four runs coincide within a few percent (figure, bottom-middle). The low-layer mode is dynamically decoupled from the evaporation layer at the accuracy of these diagnostics.

**Mass conservation is not degrading.** Column mass over 4000 s: −0.0183 % (N1000), −0.0172 % (N500), and the drift is *decreasing* in magnitude with time, not accumulating.

---

## 9. Physical vs numerical classification

**Boundary-coupled physical mode.** Confidence: **high** for the excitation and trapping mechanism; **medium** for the restoring physics that sets the ~1.7e3 s period.

Supporting the classification:
- vanishes (by 3 orders of magnitude) when conduction is off → physically excited, not a discretisation residual;
- present at equal or larger amplitude under Rusanov and under both reconstructions → not scheme-specific;
- spatially resolved over ~200 cells with smooth thermodynamics → not a grid instability;
- p–V quadrature, standing structure, node at the wall → a genuine oscillatory eigenmode of the domain-plus-boundary system;
- the initial disturbance propagates at exactly c_s and its reflection off the frozen V = 0 wall is directly captured;
- damping is purely numerical (the model has no physical damping channel for it), which is why N1000 rings for the whole run while N500 decays.

Against a purely "physical transient" reading (option A): the *trapping* is a boundary artifact. A non-reflecting lower boundary would let the conduction-launched pulse leave the domain instead of ringing. So this sits at **option B** in the mission's taxonomy, not A, C, or D.

Not supported: A (the ringing is not free physical evolution — a reflecting wall sustains it), C (conduction-off control), D (Rusanov control).

---

## 10. Release impact

**NON-BLOCKER** for the Roe + `(ln ρ, V, ln p)` release.

Checked against every blocker criterion:

| Criterion | Result |
| --- | --- |
| Conservative evaporation flux oscillatory or reversed | **No** — positive and monotone at all times; 2.1 % modulation at N1000 vs 2.0 % at N500 |
| TR solution materially modulated by the low-layer mode | **No** — amplitude profiles of all four runs coincide above 2050 km |
| Mode amplitude grows without bound | **No** — bounded, non-growing, mildly decaying |
| Mass conservation degrades | **No** — column mass −0.018 % over 4000 s, drift decreasing |
| Generated by an erroneous boundary or source treatment | **No** — the source (conduction) is the intended physics; the boundary is reflecting *by design* and its reflection does not reach the science observable |
| N1000 fails to converge toward a physical solution | **No** — mean evaporation flux converges to 2.5 % between N500 and N1000 |

The existing documented limitation stands unchanged and is unrelated to this mode: **N500 long-duration evaporation mass flux is not established as < 10 % grid-converged.**

---

## 11. Next action

**Smallest scientifically necessary step: none required for release.** Record the mode as a known, characterised feature of the current lower boundary condition.

If it is later worth removing (e.g. before quoting time-resolved chromospheric velocities, as opposed to the TR evaporation flux), the single bounded follow-up is:

> Re-run 1000–2000 s at N1000 with the inner boundary moved ~1–2 pressure scale heights deeper (`ISO_H_BASE` lowered, `ISO_DH` raised to keep the top fixed) and check whether the mode amplitude at 1800 km drops. If it does, the mode is boundary-trapped as diagnosed and a non-reflecting inner condition becomes a scoped future task; if it does not, the restoring mechanism is internal to the conduction–column coupling.

Do **not** add artificial diffusion, change the reconstruction, tune Roe, or implement NSCBC on the strength of this diagnosis.

---

## Unrelated issue found (recorded, not actioned)

The `.faceflux` sidecar's **flux values** agree with cell-centred `ρV` (§1) and with the frozen `eq_wb` reference (§5), but its **divergence is inconsistent with the density evolution recorded in the same file**. Taking two adjacent records 5.6 s apart in `lnp_roe_N1000_partial2080s.txt.faceflux`, `−Δf_total/Δs − eq_residual_mass` exceeds the measured `Δρ_cell/Δt` by ~60× at 1700–1800 km, ~23× at 2000 km, ~7× at 2100 km, and ~0.6× at 2140 km — right sign, wrong magnitude, with a systematic height trend. The capture is taken from the corrector-stage `flux_iph` in `capture_gamma_face_flux` (`src/rhs.cpp:579`) and armed once per step in `chromo_main.cpp:502`, so it is *intended* to be the flux the update differences.

This does not affect the present diagnosis — density and velocity were measured directly from the state, and the flux comparison in §1 uses flux *values*, not their divergence. But **`.faceflux` divergences (`R_total`, `R_eff` in `util/face_flux_diag.py`) should not be trusted quantitatively** until the capture's stage consistency is checked. Suggested separate bounded task: verify whether `advance_Euler_state` evaluates the RHS more than once per step, or whether the captured face index convention at the domain ends omits the boundary faces.
