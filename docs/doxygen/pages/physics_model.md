# Physical model {#physics_model}

Two physical models live in this repository. The @ref release_solver "release solver" advances a single-fluid equilibrium mixture; the @ref two_fluid_solver "historical two-fluid solver" advances separate ion, neutral and electron carriers. They share the field-aligned reduction and the mesh, and nothing else.

## Field-aligned reduction

Both models assume the flow is strictly parallel to **B**, so the 3D divergence collapses onto a one-dimensional arc-length derivative with a `B`-weighted Jacobian. Flux conservation makes the tube cross-section `A(s)` satisfy `A(s) B(s) = const`, so `1/A` and `B` are interchangeable as the geometric weight, and a variable `B(s)` produces a pressure-area source term rather than a transverse flux. The magnetic field is *prescribed*: it is a fixed background geometry, not an evolved variable, and there is no Lorentz force in the momentum equation beyond the geometric term.

Gravity enters as a prescribed potential `phi_g(s)` sampled at cell faces (`Grid::phi_g_imh`, `Grid::phi_g_iph`). For a straight vertical column `phi_g = g (h - h_base)`; for a traced field line it is whatever the geometry gives, so a closed loop has its potential maximum at the apex.

---

# The release model {#release_physics}

## Governing equations

One fluid, one velocity `u`, one temperature `T`, on the prescribed flux tube:

```
d_t(rho)   + (1/A) d_s(A rho u)         = 0
d_t(rho u) + (1/A) d_s(A (rho u^2 + p)) = p d_s(ln A) - rho d_s(phi_g)
d_t(E)     + (1/A) d_s(A (E + p) u)     = -(1/A) d_s(A q)
```

The conserved state is exactly the three rows of `chromosphere::mix` — `RHO`, `MOM`, `ENERGY`. There is no fourth row, no carrier storage, and no equilibrium-projection stage.

## Closure

```
E     = e_int(rho,T) + 1/2 rho u^2 + rho phi_g
e_int = 3/2 p + x n_H chi_H
p     = (1 + x) n_H k_B T
n_H   = rho / m_H
x     = x_Saha(n_H, T)
q     = -kappa(n_e, n_HI, T) d_s T
c_s   = sqrt(Gamma1(T, n_H) p / rho)
```

Three points about this closure carry most of the design:

**The ionization energy is inside the total energy.** `e_int` includes `x n_H chi_H`, so the ionization/recombination enthalpy is transported and conserved along with everything else, and no separate energy bookkeeping is needed when the gas ionizes. The identity `e_int - 3p/2 = x n_H chi_H` is a regression test.

**Every carrier quantity is derived, never advanced.** The ionization fraction `x`, the carrier densities `x rho` and `(1-x) rho`, the electron density `n_e = x n_H`, the neutral density `n_HI = (1-x) n_H` and the electron partial pressure `p_e = x n_H k_B T` are all functions of `(rho, e_int)` produced by the equilibrium closure at decode time. They exist as diagnostics — in the `.gamma_diag` sidecar and in chromosphere::MixtureThermo — and nowhere else.

**Three distinct adiabatic indices, and which one the release uses where.** In a partially ionized gas these are genuinely different numbers, and conflating them is a real modelling error.

- The *caloric* index relating `e_int` to `p` is analytic from the Saha closure.
- The *equilibrium acoustic* index `Gamma1 = (d ln p / d ln rho)_s` sets the isentropic sound speed **when the ionization balance is free to re-equilibrate inside the wave**. It is read from a tabulated CRASH `Gamma1(log T, log n_H)` grid (`data/eos/gamma1_hydrogen_v1.dat`, checksum-verified by a CTest), and through the hydrogen ionization zone it dips well below 5/3 — median `1.090` in the release column, below 1.5 over 94.9 % of cells. Using the caloric index in its place would misplace this speed.
- The *frozen* index is `5/3`, the monatomic translational value, and it is the correct acoustic index **when the ionization state cannot change on the timescale of the disturbance**. It follows directly from the exact split `e_int = p/(5/3 - 1) + e_ion`: with `e_ion` held fixed, only the translational term responds.

Neither acoustic index is universally correct; they are the two limits of the same closure. Chromospheric hydrogen ionization/recombination relaxes on `10^3`–`10^5 s` in the mid-chromosphere, dominated by slow leakage out of the ground state (Carlsson & Stein 2002), against a chromospheric hydrodynamic timescale of order a minute (Leenaarts 2020), so an individual compressive disturbance sees an essentially frozen composition. **The release therefore uses the frozen index `5/3` in the face Riemann problem and in the CFL**, and keeps `Gamma1` as the equilibrium derivative of the closure — reported in the `.gamma_diag` sidecar, used by the outer-boundary Mach cap, and linearized by the Roe reference flux. See @ref numerics. A fourth quantity, `(dp/de_int)_rho`, is what the Roe characteristic vectors need and is carried separately in chromosphere::MixtureFaceState.

**Open modelling question.** The Riemann solve freezes the composition, but the EOS decode re-imposes instantaneous Saha equilibrium at the end of every timestep. Those are separate modelling choices, and the second is not justified by the first. See the limitations discussion in the writeup (`docs/chromosphere_writeup.pdf`).

## Conduction

Field-aligned heat conduction is the only non-hydrodynamic physics in the release. It is **physical only**: Spitzer electron conduction plus neutral-hydrogen conduction, with no TRAC broadening factor and no artificial or numerical diffusivity anywhere in the operator. See `chromosphere::physical_kappa_e`, `chromosphere::physical_kappa_n` and `chromosphere::physical_conductivity`.

The conductive driving comes from an external reservoir temperature imposed at the **physical outer face**, not at a ghost centre — so the flux distance is half the top cell width, and `kappa` is evaluated at the face state reconstructed from the imposed external pressure and that face temperature. Keeping the thermal datum separate from the hydrodynamic ghost state is deliberate: `Grid::outer_conduction_temperature` and the hydro ghost buffers are never conflated, which lets the hydro ghost temperature free-float while the conductive wall stays fixed.

## Boundaries

The lower boundary of the release column sits at 1600 km. **It is a mid-chromospheric truncation of the modeled domain, not the photosphere.** It marks the height below which the Saha/LTE closure is not intended to be trusted; the atmosphere below it is outside the model. Its physical content is a quasi-static stratified reservoir at the truncation temperature with `v = 0`, and the two inner ghosts are a numerical closure for the MUSCL stencil — not a claim that LTE extends downward.

The boundary closures capture fixed reservoir data once from the initial condition (`kOuterPRef`, `kInnerTRef`, `kInnerPRes`), as any truncated-domain problem must. This is the only frozen reference data in the release; the interior operator itself is reference-free (see @ref numerics).

## What the release deliberately does not contain

TRAC, nonthermal beam heating, volumetric coronal heating, radiative cooling, and artificial or numerical conduction are **not release physics**. They are not present in the release timestep at all — not even behind default-off flags. They exist only in the two-fluid solver below.

---

# The experimental two-temperature model {#two_temp_physics}

**Not release physics.** This is a test study, run only by `model_column_2t` on the solver in `src/two_temp/`, asking one question the release cannot answer from inside its own assumptions: **how far do the electron and heavy-particle temperatures drift apart in the modelled 1600–2153 km column?** The release path is unaffected by its presence.

## Governing equations

One fluid, one velocity, **two** temperatures, on the same prescribed flux tube. Electrons are massless and comoving — quasi-neutrality plus zero field-aligned current — so there is still exactly one momentum equation. Protons and neutral hydrogen are lumped into a single **heavy** population at `T_i`: they exchange energy through charge exchange and elastic collisions far faster than anything resolved here, and the release model has no ion/neutral drift either.

```
d_t(rho)   + (1/A) d_s(A rho u)         = 0
d_t(rho u) + (1/A) d_s(A (rho u^2 + p)) = p d_s(ln A) - rho d_s(phi_g)
d_t(E)     + (1/A) d_s(A (E + p) u)     = -(1/A) d_s(A (q_e + q_i))
d_t(E_e)   + (1/A) d_s(A E_e u)         = -p_e (1/A) d_s(A u)
                                          -(1/A) d_s(A q_e) + Q_ei
```

`E` keeps its release meaning — the **total** energy — so mass, momentum and total energy stay in exact conservation form and the release numerical flux applies to those three rows unchanged. The heavy internal energy is then **derived**, `e_i = E - E_e - 1/2 rho u^2 - rho phi_g`, which makes `T_i = (2/3) e_i / (n_H k_B)` purely algebraic with no inversion at all. This "total energy plus electron energy" split is the standard two-temperature MHD formulation (AWSoM/BATS-R-US; van der Holst et al. 2010, Sokolov et al. 2021), and it is the same split the historical two-fluid `cons::E_E` row already uses.

## Closure

```
n_H  = rho / m_H                       (protons + neutrals)
x    = x_Saha(n_H, T_e)                (ionization is an ELECTRON process)
n_e  = x n_H,   n_HI = (1 - x) n_H
p_e  = n_e k_B T_e
p_i  = n_H k_B T_i
p    = p_e + p_i = n_H k_B (T_i + x T_e)
E_e  = 3/2 p_e + x n_H chi_H           (the electron pool OWNS the ionization reservoir)
E    = 3/2 p_i + E_e + 1/2 rho u^2 + rho phi_g
q_e  = -kappa_e(n_e, n_HI, T_e) d_s T_e
q_i  = -kappa_n(n_e, n_HI, T_i) d_s T_i
Q_ei = g_ei (T_i - T_e),  g_ei = 3/2 n_e k_B (nu_ei + nu_en)
```

**At `T_e = T_i = T` this degenerates exactly to the release closure**: `p = (1+x) n_H k_B T` and `E_e + (3/2) p_i = (3/2) p + x n_H chi_H = e_int`. That exact degeneracy is the experiment's control, and it is enforced in code by taking the electron caloric pair as literal differences of the public release EOS functions rather than re-deriving them.

**Why the ionization energy sits in the electron pool, and why Saha is evaluated at `T_e`.** Hydrogen ionization and recombination in this regime are electron-impact processes: the energy `chi_H` per ionization is taken from, and returned to, the electron thermal reservoir. Putting it anywhere else would mean the heavy population paid for a transition it does not mediate.

## Conduction and collisional coupling

The two conductive channels are the release's own two, each now evaluated at the temperature that actually drives it: Spitzer `kappa_e` at `T_e`, neutral-hydrogen `kappa_n` at `T_i`. Spitzer **proton** conduction is deliberately omitted — about 2.3 % of `kappa_e` for hydrogen — precisely so that `kappa_e + kappa_n` reproduces the release total conductivity *exactly* at `T_e = T_i`, leaving the temperature split as the only changed variable.

The equilibration conductance `g_ei` uses the translational electron capacity `(3/2) n_e k_B` times `nu_ei + nu_en`, with the Spitzer/NRL electron-ion rate and the Vranjes & Krstic (2013) electron-neutral cross-section. Both channels matter: the lower chromosphere is weakly ionized, so `n_HI >> n_e` there and electron-neutral collisions carry a large share of the electron thermalization. The resulting local equilibration time across this column is short — of order 0.09 ms at the base rising to a few ms at the top — which is why the initial condition sets `T_e = T_i` and lets the solver, not the IC, generate any split.

## Boundaries

The conductive reservoir at the top is assigned to the **electron** channel (Dirichlet `T_e = 22,000 K`), with the heavy channel insulating, because conduction down from the transition region and corona is electron-conducted and there are no neutrals above the domain to supply a heavy flux. The base is Dirichlet on both temperatures, which is a property of a collisionally equilibrated dense reservoir rather than an added assumption. Both *hydrodynamic* temperatures free-float at the outflow face, because the subsonic-outflow characteristic budget admits exactly one imposed condition and the reservoir back-pressure already spends it. The full derivation is in @ref numerics-two-temp-bc.

## What this model deliberately does not contain

No finite-rate ionization, no radiative cooling, no TRAC, no beam, no volumetric heating, no ion/neutral drift, no anisotropic ion temperature. **The release's open modelling question carries over unchanged and is if anything sharper here**: the Riemann solve freezes composition while the EOS decode re-imposes instantaneous Saha equilibrium after every step — now at `T_e`. The algorithm is still `Saha-equilibrated state -> frozen-composition hydro step -> Saha-equilibrated state`, an operator-split instantaneous-relaxation approximation that is not justified when the physical relaxation time is long.

---

# The historical two-fluid model {#two_fluid_physics}

Not the release path. Kept because it is the only path with finite-rate ionization, separate species temperatures and the full source-term catalogue, and because several scenarios still run on it.

## Conserved state

Seven rows (`chromosphere::cons`): `RHO_I`, `RHO_N`, `MOM_I`, `MOM_N`, `E_I`, `E_N`, `E_E`. Electrons are folded into the ion equations by quasi-neutrality (`n_e = n_i`), so the ion pressure carries an extra factor of two: `p_i = 2 n_i k_B T_i`, while `p_n = n_n k_B T_n`.

The three-temperature extension keeps `E_I` meaning the **total** charged-fluid energy (protons plus electrons plus bulk kinetic plus gravitational) and adds `E_E` as the electron internal energy. The proton temperature is then derived, `T_i = (p_total - p_e)/(n_i k_B)`, and `T_e = p_e/(n_i k_B)` with `p_e = (2/3) E_E`. This total-energy-plus-electron-energy split is the standard two-temperature MHD formulation and it guarantees that `ENABLE_TE=0` reproduces the single-temperature baseline exactly, because `E_E` never feeds back into `E_I`, momentum or density.

## Source physics

| Term | Where | Reference |
| --- | --- | --- |
| Ion-neutral drag and frictional heating | Stage B | collision frequency `chromosphere::nu_in` |
| Temperature equilibration `T_e`, `T_i`, `T_n` | Stage C | `chromosphere::nu_ei`, `chromosphere::nu_en`, `chromosphere::nu_in` |
| Field-aligned conduction per species | Stage D | `chromosphere::kappa_e`, `chromosphere::kappa_n` |
| Hydrogen ionization / recombination | Stage E | see below |
| Radiative cooling | Stage R | see below |
| Nonthermal beam heating | beam stage | `chromosphere::beam_heating_rate` |
| Ambient coronal heating `H(s)` | heating stage | `chromosphere::coronal_heating_rate` |
| TRAC conductivity broadening | Stages D and R | `chromosphere::trac_broadening_factor` |

### Ionization network

The complete network is the default, because no single channel dominates across the range this code spans — a cool weakly ionized 6 kK chromosphere, a dense flare condensation, and a fully ionized 10^7 K evaporated plasma. Each channel leads in some regime:

- **Ionization**: photoionization `P_phot` (a frozen-radiation-field rate, Chae 2021 FAL-C Table 1, via `chromosphere::photoionization_rate_chae`), multilevel collisional `S_CR` (`chromosphere::ionization_rate_S_CR` — the Rydberg-ladder channel, dominant in the cool chromosphere), and direct ground-state electron impact `S_i` (`chromosphere::ionization_rate_S`, Voronov 1997 — subdominant below 10^5 K but the physically robust channel above it, where the cool-plasma `S_CR` fit is an extrapolation).
- **Recombination**: case-B radiative `alpha_r` (`chromosphere::recombination_rate_alpha`, Hummer 1994) and three-body `alpha_c = kappa_c n_e` (`chromosphere::recombination_rate_alpha_c`, Hinnov & Hirschberg 1962 / Stevefelt 1975). The three-body rate scales as `n_e^3` and is the dominant sink wherever the gas is compressed — the dense base and, critically, the shock-compressed flare condensation, where dropping it destabilizes the run. Including it makes the Stage E solve cubic rather than quadratic.

`S_CR` and `alpha_c` are detailed-balance inverses through the Saha function `chromosphere::saha_phi`, which is why they share a coefficient.

### Radiative losses

Two recipes stitched with a smooth switch near 2x10^4 K so neither double-counts the other:

- **Optically thick chromospheric** cooling (`chromosphere::radiative_loss_thick`): the Carlsson & Leenaarts (2012) recipe summed over H I, Ca II and Mg II, with the escape probabilities driven by column integrals accumulated downward from the top of the column. Tables in `chromosphere::cl2012`.
- **Optically thin TR/coronal** cooling (`chromosphere::radiative_loss_thin`): `Q = n_e n_H Lambda(T)` with a CHIANTI-class coronal-abundance loss function tabulated in `chromosphere::optthin`. This is the sink the lower transition region needs in order to re-radiate an imposed downward coronal conductive flux; without it the imposed flux runs the top of the column up to coronal temperatures. Note that equilibrium `Lambda(T)` underestimates lower-TR losses by up to a factor of three under non-equilibrium ionization — a known conservative bias.

### Evaporation drivers

Two regimes, distinguished by how fast energy arrives:

- **Explosive** (Fisher, Canfield & McClymont 1985): a transient nonthermal electron beam deposits above the threshold `F_crit ~ 7x10^6 W/m^2` in the upper chromosphere, which then cannot radiate the deposited energy, is heated to coronal temperatures, and drives a supersonic upflow with a condensation downflow beneath it. The threshold is a competition between deposited heating and the upper-chromosphere radiative loss, so the radiative sink must be active for the threshold to mean anything. Two deposition models are available: a fixed density-weighted height window, or a self-consistent Emslie (1978) thick-target column integral that lets the stopping depth migrate up the loop as evaporation fills it.
- **Gentle** (Antiochos & Sturrock 1978): the ambient footpoint-anchored coronal heating `H(s) = E_H0 exp(-d(s)/s_H)` is ramped up slowly, staying sub-threshold so the flow remains far below the sound speed. `H(s)` is the missing term in the Rosner-Tucker-Vaiana static loop balance `d/ds(kappa_e T^{5/2} dT/ds) + H - n^2 Lambda = 0`; without it a resolved corona simply drains.

## Known issue

The write-up's semi-implicit Picard iteration is unstable for the implicit terms actually present here. At Model C7 conditions the drag is stiff, and Picard iteration converges only when `dt |R_I'| < 1`. The body of `chromosphere::rhs_implicit_state` is correct, but the call site in `src/two_fluid/integrators.cpp` keeps the residual zeroed, so the two-fluid step is effectively explicit Euler plus the point-implicit and tridiagonal stages until a real Newton solve on the `R_I` Jacobian replaces the Picard iteration.

## See also

@ref numerics for the discretization, @ref scenario_reference for which scenario runs which model, @ref eos for the closure API.
