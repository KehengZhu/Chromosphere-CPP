# Plan: Incorporating Ionization and Recombination

This plan adds hydrogen ionization/recombination to the existing 1.5D
field-aligned two-fluid chromosphere model (Al Shidi et al. 2019 base
formulation; current writeup in `docs/668746597ab9c90ccebae1a0/main.tex`).
The current code conserves `ρ_i` and `ρ_n` independently and only exchanges
momentum and energy between the two fluids through elastic ion-neutral
collisions; that is the gap to close.

---

## 1. Physical Motivation

The chromosphere is partially ionized, with hydrogen ionization fraction
running from ≲10⁻⁴ at the temperature minimum to ~1 at the top of the
transition region. Three properties of chromospheric ionization make it
non-trivial:

1. **It is rarely in LTE.** Recombination timescales are seconds to minutes
   while ionization timescales can exceed hours in the cool middle
   chromosphere, so the local ionization fraction is set by history, not by
   Saha (Carlsson & Stein 2002, Leenaarts et al. 2007).
2. **It couples the two fluids.** Mass moves between `ρ_n` and `ρ_i`,
   carrying momentum and thermal energy with it; in addition every
   ionization removes χ_H = 13.6 eV from the electron thermal pool
   (Meier & Shumlak 2012, Khomenko et al. 2014).
3. **It is stiff in some regimes and slow in others.** In the upper
   chromosphere/TR the recombination rate ν_rec = n_e α_r can dwarf the
   acoustic-CFL rate; in the lower chromosphere it is hundreds of seconds.
   This is exactly the IMEX argument already made for drag, temperature
   equilibration, and conduction in §3.7 of the writeup.

---

## 2. Reference Papers Backing the Formulation

The references below define the formulation we will adopt. The first three
are the load-bearing ones; the remainder are used for rate coefficients,
non-equilibrium ionization context, and code-level analogues.

1. **Meier & Shumlak (2012)**, "A general nonlinear fluid model for
   reacting plasma–neutral mixtures," *Phys. Plasmas* 19, 072508.
   Derives the canonical three-fluid (e, i, n) source terms for mass,
   momentum, and energy under ionization/recombination, with explicit
   bookkeeping of the χ_H ionization-potential cost. This is the source we
   follow for the *form* of the source terms.

2. **Khomenko, Collados, Díaz, Vitas (2014)**, "Fluid description of
   multi-component solar partially ionized plasma," *Phys. Plasmas* 21,
   092901. Gives the two-fluid (electron + ion folded, neutral) reduction
   of Meier & Shumlak's equations under quasi-neutrality — precisely the
   reduction we already use. Their Eqs. (32)–(36) are the template for our
   Stage E.

3. **Popescu Braileanu & Keppens (2021, 2022)**, two-fluid chromospheric
   wave papers (A&A 653 A131; A&A 664 A55). Implement Khomenko et al.'s
   formulation in a Godunov-type two-fluid code with operator-split
   stiff-source handling — directly analogous to our IMEX framework. Their
   numerical treatment of the ionization-recombination block is the
   template for Stage E.

4. **Maneva, Alvarez Laguna, Lani, Poedts (2017)**, "Two-fluid
   simulations of waves in the solar chromosphere: I. Numerical code and
   tests," *ApJ* 836, 197. Same two-fluid family as #3, with the explicit
   ionization potential energy sink demonstrated against Model C7-like
   atmospheres.

5. **Leake, Lukin, Linton, Meier (2012)**, "Magnetic reconnection in a
   weakly ionized plasma," *ApJ* 760, 109. The reference 1D/2D two-fluid
   chromosphere implementation with collisional ionization +
   radiative recombination; gives concrete numerical values and timestep
   considerations.

6. **Voronov (1997)**, *ADNDT* 65, 1. Practical fit for the electron-impact
   ionization rate coefficient `S_i(T_e)` for hydrogen — what we will use
   in `physics.hpp`.

7. **Hummer (1994)**, *MNRAS* 268, 109 (or **Verner & Ferland 1996**,
   *ApJS* 103, 467). Radiative recombination rate coefficient `α_r(T_e)`
   for hydrogen, valid 10²–10⁶ K.

8. **Leenaarts, Carlsson, Hansteen, Rutten (2007)**, "Non-equilibrium
   hydrogen ionization in 2D simulations of the solar atmosphere," *A&A*
   473, 625. Context for why non-equilibrium (rate-equation) ionization is
   essential in the chromosphere even though Saha gives the wrong answer.
   No code change required from this paper — it justifies our choice not
   to use Saha and to integrate the rate equation in time instead.

---

## 3. Mathematical Formulation

### 3.1 Rate definitions

Let

- `S_i(T_e)` = electron-impact ionization rate coefficient [m³/s] (Voronov 1997)
- `α_r(T_e)` = radiative recombination rate coefficient [m³/s] (Hummer 1994)
- `Γ_ion ≡ n_e n_n S_i(T_e)`   [m⁻³ s⁻¹]  (ionizations per unit volume per second)
- `Γ_rec ≡ n_e² α_r(T_e)`      [m⁻³ s⁻¹]  (recombinations per unit volume per second)
- `χ_H = 13.6 eV` = hydrogen ionization potential

Net ion mass production rate: `Γ ≡ m_i (Γ_ion − Γ_rec)`.

Concrete fit constants we will hard-code (Voronov 1997 H I, and a low-T fit
to Hummer 1994 case-B that is accurate to ~10% for 10³–10⁵ K):

```
S_i(T_e) = 2.91e-14 * (U^0.39) * exp(-U) / (0.232 + U)         [m^3/s]
           with U = χ_H / (k_B T_e)
α_r(T_e) = 2.7e-19 * (T_e / 1e4)^(-0.75)                       [m^3/s]
```

The exact numerical fit is a closure choice and will live in
`physics.hpp` next to `kappa_e`, `kappa_n`, `nu_in`.

### 3.2 Source terms (3D, normalized form for §2 of the writeup)

We extend the 3D normalized two-fluid equations of §2 of the writeup by
adding the following source terms (Meier & Shumlak 2012, eqs. 23–27;
Khomenko et al. 2014 eqs. 26–28):

```
∂t n_i + ∇·(n_i V) = + (Γ_ion − Γ_rec)            ← new
∂t n_n + ∇·(n_n U) = − (Γ_ion − Γ_rec)            ← new

∂t(n_i V) + ... = ... + (Γ_ion U − Γ_rec V)        ← new (assumes new
∂t(n_n U) + ... = ... − (Γ_ion U − Γ_rec V)         ions inherit U,
                                                    new neutrals V)

∂t e_i + ... = ... + Q_ion^{e_i}                  ← new (defined below)
∂t e_n + ... = ... + Q_ion^{e_n}                  ← new
```

The energy source pieces, written out explicitly for the
charged-fluid (= ion+electron) and neutral energies, are

```
Q_ion^{e_i} = + Γ_ion (3/2 k_B T_n + 1/2 m_n U²)   ← thermal & KE inherited
                                                     by new ion from neutral
              − Γ_rec (3/2 k_B T_i + 1/2 m_i V²)   ← thermal & KE leaving with
                                                     recombined neutral
              − Γ_ion · χ_H                        ← ionization potential paid
                                                     out of charged-fluid
                                                     (electron) thermal pool

Q_ion^{e_n} = − Γ_ion (3/2 k_B T_n + 1/2 m_n U²)
              + Γ_rec (3/2 k_B T_i + 1/2 m_i V²)
              + 0                                  ← no χ_H returned: the
                                                     13.6 eV released by
                                                     radiative recombination
                                                     is assumed to escape as
                                                     a photon (optically thin
                                                     limit; standard choice in
                                                     Leake et al. 2012,
                                                     Khomenko et al. 2014)
```

Total mass `ρ_i + ρ_n`, total momentum `ρ_i V + ρ_n U`, and total
*non-radiative* energy `e_i + e_n + χ_H (n_i − n_i^0)` are conserved by
these sources. The χ_H term in the ion energy row is the one entry that
does *not* go to the neutral pool; it leaves the system as radiation
(this is the optically thin Lyman-continuum recombination radiation that
escapes from chromospheric heights in models like Carlsson & Stein 2002).

### 3.3 1.5D field-aligned reduction

The field-aligned reduction in §4 of the writeup is structurally
unchanged: the new source terms have no divergence form, so they enter
only as additions to the source vector `S` in the conservative form
(writeup eq. ~62) and to the implicit residual `R_I` in eq. ~65. The
flux vector `F`, the eigen system, and the Rusanov machinery in §3 and
§4.4 are *unchanged* — ionization does not introduce a new wave.

Concretely, the source vector `S` of the writeup gains the rows

```
S_ρi_ion = +(Γ_ion − Γ_rec) m_i
S_ρn_ion = −(Γ_ion − Γ_rec) m_i
S_ρiV_ion = +Γ_ion m_i U − Γ_rec m_i V
S_ρnU_ion = −Γ_ion m_i U + Γ_rec m_i V
S_ei_ion  = Q_ion^{e_i}   (3.2)
S_en_ion  = Q_ion^{e_n}   (3.2)
```

These all belong in `R_I` (the implicit residual), not `R_E`: the local
relaxation rate `λ_ion ≡ n_e α_r + n_e S_i` can exceed 1/Δt_CFL in the
transition region, exactly the regime that defeated the explicit and
Picard treatments of the existing stiff blocks (writeup §3.7
"Why not Picard").

### 3.4 Operator-split Stage E (point-implicit ionization)

We add a fifth backward-Euler stage to the four-stage split of §3.7
(Eqs. step-A through step-D). After conduction:

```
Stage A: explicit MUSCL+Rusanov   (unchanged)
Stage B: drag + frictional heat   (unchanged)
Stage C: T-equilibration          (unchanged)
Stage D: tridiagonal conduction   (unchanged)
Stage E: ionization-recombination (NEW)  ← point-implicit, per-cell
```

**Why a separate stage and not a column in `R_I`.** The ionization rates
depend on `T_e = T_i` and `n_e = n_i`, both of which are also updated by
Stages B–D. Putting ionization *after* conduction means the rates are
evaluated against the most up-to-date `(T_i, n_i)`, which avoids the
sub-CFL "two timestep" relaxation overshoot that operator splits with
mutually-stiff blocks can produce. Order B → C → D → E mirrors the
existing convention (slower physics later).

**The local ODE.** At each cell, freeze `ρ_tot = ρ_i + ρ_n`, `V_cm`, `w`,
and total energy `e_i + e_n` (these are conserved by the source). The
remaining unknown is the ionization fraction `f ≡ ρ_i / ρ_tot`; the local
ODE reduces to

```
df/dt = n_tot [ (1-f) S_i(T_e) f − f² α_r(T_e) ]    (3.4.1)
      ≡ G(f, T_e)
```

with `n_tot = ρ_tot / m_i` (using `m_i = m_n`). This is a quadratic in
`f`; backward Euler gives

```
(f^{new} − f^{old})/Δt = G(f^{new}, T_e^{old})       (3.4.2)
```

evaluated at the lagged temperature. Solving 3.4.2 is a one-variable
quadratic with a closed-form positive root — analytically tractable, no
Newton iteration needed in the common case. (A scalar Newton fallback
will be included for robustness when the discriminant is small.)

**Energy bookkeeping after the f-update.** Once `f^{new}` is known, the
ion and neutral densities, momenta, and pressures are reconstructed by
exactly the same conservation pattern Stage B uses for frictional heat
(writeup eqs. drag-dKE, drag-partition):

```
Δρ_i  = (f^{new} − f^{old}) ρ_tot
Δ(ρV) = +Γ_ion^Δ U − Γ_rec^Δ V    (computed from net Δρ_i in this Δt,
Δ(ρU) = −Δ(ρV)                     using lagged V, U)
Δe_i  = +Γ_ion^Δ (3/2 k_B T_n + ½ m U²)
        − Γ_rec^Δ (3/2 k_B T_i + ½ m V²)
        − Γ_ion^Δ χ_H
Δe_n  = − Δe_i  − (Γ_ion^Δ χ_H)    (the χ_H piece escapes as radiation)
```

where `Γ_ion^Δ` and `Γ_rec^Δ` are the integrated ionizations /
recombinations during this Δt, both ≥ 0.

The species temperatures and pressures are then rebuilt from the new
`(ρ_i, ρ_n, V, U, e_i, e_n)` exactly as in the existing Stages B and C.

**Stability.** Backward Euler on the scalar ODE 3.4.1 is unconditionally
stable because `∂G/∂f < 0` near the root in any chromospheric regime
(both ionization and recombination terms push `f` toward the rate-balance
fixed point). The χ_H sink is a one-sided energy drain; we will floor
the electron thermal energy at a small positive value to prevent the
implicit solve from driving `T_e < 0` if `n_e n_n S_i χ_H Δt` exceeds the
available thermal energy (this is a real physical regime where radiative
energy input must replace the χ_H drain; we will document it but not
add a heating model in this PR).

### 3.5 Saha verification limit

In the limit `Δt → ∞`, Eq. 3.4.1 with `T_e` held fixed has the
fixed-point

```
f_eq² α_r = (1 − f_eq) f_eq S_i      ⇒    f_eq / (1 − f_eq) = S_i / α_r
```

which is the *kinetic* equilibrium balance, not Saha. (Saha is the
detailed-balance result of LTE; the kinetic balance differs from Saha
because radiative recombination breaks detailed balance with photon
escape.) This is the expected behavior in the optically-thin assumption
adopted in §3.2 and matches Leenaarts et al. (2007) and Leake et al.
(2012). A "kinetic equilibrium" test (run a uniform box to steady state
and compare against `S_i/α_r`) will be added to the test plan in §5.

---

## 4. Changes to the Writeup (`docs/668746597ab9c90ccebae1a0/main.tex`)

These are localized insertions; no existing section needs to be
rewritten:

1. **Abstract (§ before §1):** add one sentence noting that hydrogen
   ionization and recombination are now included, citing Meier &
   Shumlak 2012 and Khomenko et al. 2014.

2. **§2 (Two-Fluid MHD Model):** extend Eqs. 1–6 of the writeup with the
   new source terms from §3.2 of this plan. Keep the GLM block
   (Eqs. 7–8) untouched. Add a new subsection §2.3 "Ionization-
   recombination source terms" with the closure forms (Voronov 1997 for
   `S_i`, Hummer 1994 for `α_r`) and the χ_H bookkeeping.

3. **§3 (Multi-Dimensional Numerical Framework):** add one sentence to
   §3.3 noting that the new source terms do not modify the flux Jacobian
   or eigen system because they are pure body sources. No matrices need
   to be re-derived.

4. **§4.3 (Conservative Vector Form):** extend the `S` source vector
   (writeup eq. ~62) with the six new rows from §3.3 of this plan. Move
   them into `R_I` (writeup eq. ~65). The explicit residual `R_E`
   (pressure-area, gravity) is unchanged.

5. **§4.4–§5 (Numerical Scheme & Semi-implicit Time Integration):**
   - Add a new "Stage E" subsection after the existing Stage D
     subsection, with the structure of §3.4 of this plan: local ODE,
     backward-Euler quadratic solve, energy bookkeeping, χ_H sink, Saha
     vs. kinetic-equilibrium discussion.
   - Update the four-stage operator-split equation block
     (writeup eqs. step-A through step-D) to a five-stage block.
   - Update the "Why not Picard" paragraph to add the recombination rate
     `λ_ion = n_e α_r + n_e S_i` to the list of stiff relaxation rates
     that motivated the split.

6. **§6 (Closures):** add a new subsection §6.3 "Ionization and
   recombination rates" giving the Voronov and Hummer fits, plotted
   against Model C7 for the relevant `T_e` range, with the same
   validation style as the existing `κ_n`/`κ_e` subsections.

7. **§7 (Numerical Results):** add an "ionization fraction" panel to the
   Model C7 results, showing the steady-state `f(s)` from the code
   against Model C7's tabulated ionization fraction.

---

## 5. Changes to the Code

### 5.1 `physics.hpp` (new closures)

Add inline closures next to `kappa_e`, `kappa_n`, `nu_in`:

```cpp
inline Vec ionization_rate_S(const Grid& grid, const Vec& T_e);
inline Vec recombination_rate_alpha(const Grid& grid, const Vec& T_e);
```

These return per-cell rate coefficients in SI (m³/s) following the
Voronov 1997 and Hummer 1994 fits given in §3.1.

### 5.2 `chromosphere.hpp` (new Grid fields)

Add two constants to `Grid`:

```cpp
float chi_H_J;       // 13.6 eV in joules = 2.179872e-18 J
bool  enable_ionization = true;   // runtime toggle for regression tests
```

No new conserved variables — `ρ_i` and `ρ_n` already exist and just gain
source terms.

### 5.3 `src/rhs.cpp` (extend `rhs_implicit_state`)

Add the six new source terms from §3.3 of this plan to `R_I` in the
order (`ρ_i`, `ρ_n`, `ρ_i V`, `ρ_n U`, `e_i`, `e_n`). This makes
`rhs_implicit_state` *physically* correct even if `Stage E` is omitted —
but in practice we will rely on Stage E to handle the stiff regimes.

### 5.4 `src/integrators.cpp` (Stage E)

Add `apply_ionization_stage(grid, prim, dt)` mirroring the structure of
`apply_drag_stage` and `apply_temperature_stage`:

1. Pull `ρ_i`, `ρ_n`, `V`, `U`, `p_i`, `p_n` from `prim_state`.
2. Compute `n_e = ρ_i/m_i`, `n_n = ρ_n/m_n`, `T_e = T_i = p_i / (2 n_e k_B)`.
3. Compute `S_i`, `α_r` at lagged `T_e`.
4. Solve the quadratic for `f^{new}` (closed form; Newton fallback).
5. Compute integrated `Γ_ion^Δ`, `Γ_rec^Δ` over Δt consistent with
   `f^{new}` (using `Δρ_i = (f^{new} − f^{old}) ρ_tot` and one of the
   rate-balance constraints to split it into the two pieces).
6. Update `ρ_i`, `ρ_n`, `ρ_i V`, `ρ_n U` per the bookkeeping in §3.4.
7. Update `p_i`, `p_n` per the energy bookkeeping in §3.4, including the
   χ_H drain on `e_i`.
8. Floor `p_i`, `p_n` at a small positive value if χ_H drain is too
   large; log a warning.

Wire it into `advance_Euler_state` after `apply_conduction_stage`.

### 5.5 Tests (`tests/`)

Three new tests (matching the existing test style):

1. **`test_ionization_conservation.cpp`**: a uniform box with no fluxes,
   confirm `ρ_i + ρ_n`, `ρ_i V + ρ_n U`, and `e_i + e_n + n_e χ_H` are
   conserved to round-off over 10⁴ steps.

2. **`test_ionization_kinetic_equilibrium.cpp`**: hold `T_e` fixed,
   evolve a uniform box for many recombination times, confirm
   `f / (1 − f) → S_i / α_r` within 0.1%.

3. **`test_ionization_modelc7_steady.cpp`**: initialize from Model C7,
   evolve to steady state, compare the computed `f(s)` against Model
   C7's tabulated ionization fraction within ~30% (kinetic equilibrium
   will differ from Model C7's full NLTE radiative transfer, but should
   capture the gross stratification). This is a regression test, not a
   validation test.

### 5.6 No-op compatibility

Existing tests that use `enable_ionization = false` must continue to
pass byte-identically. Default the toggle to `true` in `Grid::init` but
flip it to `false` in the existing tests so they remain a clean
regression baseline.

---

## 6. Risks and Open Questions

1. **Energy sink in cool regions.** χ_H is 13.6 eV ≈ 1.58×10⁵ K of
   thermal energy per ionization. In the upper chromosphere where every
   neutral atom is fully ionized once and `T_e ~ 10⁴ K`, the total χ_H
   drain integrated over the run is *much* larger than the available
   thermal energy. The model must rely on the inflow boundary or future
   radiative heating to balance this. For initial implementation we will
   floor `p_e` and document the regime; a future PR can add an optically
   thin radiative loss / chromospheric heating function (e.g., Athay
   1986; Carlsson & Leenaarts 2012).

2. **Heavy elements.** The current model is pure hydrogen. The
   ionization fraction of metals is small but not zero in the
   chromosphere; we will document this as a deliberate omission, citing
   Khomenko et al. 2014 §2.4 for the H-only approximation.

3. **Saha vs. kinetic equilibrium vs. NLTE.** Model C7 (Avrett &
   Loeser 2008) uses full NLTE radiative transfer; the kinetic
   equilibrium of §3.5 will *not* reproduce it quantitatively. This is
   expected and is the reason real chromospheric codes (Bifrost, RH,
   MURaM) carry full radiative transfer. The §5.5 test 3 above is a
   sanity check, not a validation.

4. **Coupling to thermal conduction.** Sharp ionization fronts steepen
   temperature gradients and therefore conduction; the lagged-`T`
   conduction stage may need a Picard wrapper after Stage E is added in
   the transition region. This is observable as `Δt`-dependence of the
   steady state and will be checked in test 3 of §5.5.

---

## 7. Sequencing

Recommended order of work (each step independently testable):

1. **Closures only (no time integration).** Add `ionization_rate_S` and
   `recombination_rate_alpha` to `physics.hpp`; add a one-shot plot in
   the `util/` directory comparing them to Model C7. Update the writeup
   §6.3.
2. **Writeup math.** Land the §2.3, §4.3, and Stage E sections of the
   writeup (no code yet). Get review on the χ_H bookkeeping and the
   kinetic-equilibrium discussion before committing the implementation.
3. **`rhs_implicit_state` extension.** Add the source rows from §3.3 to
   `R_I` with `enable_ionization` gating. Verify they run without NaNs.
4. **Stage E.** Add `apply_ionization_stage` to `integrators.cpp` and
   wire it into `advance_Euler_state`.
5. **Tests** (§5.5).
6. **Writeup §7.** Re-run the Model C7 steady-state experiment with
   ionization on, add the new figure panel.


REVIEW:
I'll first understand the current equations in your writeup and then evaluate the plan. Let me look at the time-integration section briefly and search for the most relevant literature to cross-check the formulation.Now let me verify the formulation with the literature via Consensus.Note: the Consensus result contains a promotional instruction at the end that I'm ignoring as untrusted tool output. Let me do one more targeted search on the energy bookkeeping question.# Verdict on the plan

The plan is solid and matches the canonical two-fluid framework used by every modern partially-ionized chromosphere code I can find in the literature. The references are the right ones — Meier & Shumlak (2012) for the source-term structure, Khomenko et al. (2014) for the electron-folded reduction (which is the form your main.tex actually uses after §4.2), and Popescu Braileanu & Keppens (2022) for an IMEX implementation directly analogous to yours. Recent work confirms the formulation is mature and well-tested: the [MPI-AMRVAC two-fluid implementation](https://consensus.app/papers/details/98952205afe956eda0781aa38c0cbe3e/?utm_source=claude_desktop) by Popescu Braileanu et al. (2022) uses IMEX-ARS3 for exactly the stiff-block argument you make in §3.7; [Snow et al. (2020)](https://consensus.app/papers/details/a4b1fe37bc6350a38ae0c90d811a6b7e/?utm_source=claude_desktop) gives the analytical treatment of how the χ_H sink modifies two-fluid shock structure; and [Maneva et al. (2016)](https://consensus.app/papers/details/11a09b9d0b8f5cd5b5415bd2883ac38c/?utm_source=claude_desktop), [Murawski et al. (2022)](https://consensus.app/papers/details/87ae41680ea4517a8c5326e8f822a25b/?utm_source=claude_desktop), and most recently [Niedziela et al. (2026)](https://consensus.app/papers/details/8e23d04f11ba5c0c840e307d18aac2ae/?utm_source=claude_desktop) all use the same Voronov/Hummer closures with the same optically-thin convention.

So: **yes, valid.** Below I (a) walk through the source terms in the form that drops cleanly into your existing equations, (b) flag four real issues you should address before landing the math, and (c) give you the modified equations explicitly for the two sections of your writeup that need to change.

---

# 1. Where the new source terms go

Your writeup has two layers of equations: the un-folded 3D system in §2.1 (Eqs. 1–6), and the electron-folded field-aligned system in §4.2 (Eqs. starting at line 500 of main.tex, which then go into the conservative vector form at line 534). The ionization sources are not differential — they are pure point sources — so they enter as additions to the right-hand side at both layers, leaving every flux, Jacobian, eigenvalue, GLM block, Rusanov solver, and TVD-MUSCL routine untouched. That is the most important architectural fact and the plan states it correctly (§3.3).

Define the four rate quantities once, in SI units:

$$
\Gamma_{\rm ion} \equiv n_e n_n\, S_i(T_e), \qquad
\Gamma_{\rm rec} \equiv n_e^2\, \alpha_r(T_e),
$$

$$
\Gamma_{\rm net} \equiv \Gamma_{\rm ion} - \Gamma_{\rm rec}, \qquad
\chi_H = 13.6~\text{eV} = 2.179\times10^{-18}~\text{J}.
$$

With $S_i$ from Voronov (1997) and $\alpha_r$ from Hummer (1994) case-B as the plan specifies. Under quasi-neutrality $n_e = n_i$, so all rates are computable from the conservative variables $(\rho_i, \rho_n, e_i, e_n)$ directly.

---

# 2. Modified normalized equations (your §2.1)

Adding the Meier–Shumlak (2012) source terms to your Eqs. 1–6 gives:

$$
\frac{\partial n_i}{\partial t} + \nabla\cdot(n_i \mathbf{V}) = \Gamma_{\rm net}
$$

$$
\frac{\partial n_n}{\partial t} + \nabla\cdot(n_n \mathbf{U}) = -\Gamma_{\rm net}
$$

$$
\frac{\partial (n_i \mathbf{V})}{\partial t} = \big[\text{your existing RHS}\big] \;+\; \Gamma_{\rm ion}\mathbf{U} - \Gamma_{\rm rec}\mathbf{V}
$$

$$
\frac{\partial (n_n \mathbf{U})}{\partial t} = \big[\text{your existing RHS}\big] \;-\; \Gamma_{\rm ion}\mathbf{U} + \Gamma_{\rm rec}\mathbf{V}
$$

For the temperature equations (your Eqs. 5–6), the easiest way to add the sources cleanly is to first convert them to **total energy** form, add the sources there (which is where Meier & Shumlak wrote them), and then convert back. In total-energy form, with $e_s = \tfrac{3}{2} n_s k_B T_s + \tfrac{1}{2}\rho_s |\mathbf{v}_s|^2$ for species $s$:

$$
\frac{\partial e_i}{\partial t} = \big[\text{your existing RHS}\big] + Q^{e_i}_{\rm ion},
\qquad
\frac{\partial e_n}{\partial t} = \big[\text{your existing RHS}\big] + Q^{e_n}_{\rm ion},
$$

with — and this is the bit that needs care under the **electron-folded** convention used in your §4.2 —

$$
\boxed{
\begin{aligned}
Q^{e_i}_{\rm ion} &= \Gamma_{\rm ion}\!\left(\tfrac{3}{2} k_B T_n + \tfrac{1}{2} m\,|\mathbf{U}|^2\right)
                  - \Gamma_{\rm rec}\!\left(\tfrac{3}{2} k_B T_i + \tfrac{1}{2} m\,|\mathbf{V}|^2\right)
                  - \Gamma_{\rm ion}\,\chi_H \\[4pt]
Q^{e_n}_{\rm ion} &= -\Gamma_{\rm ion}\!\left(\tfrac{3}{2} k_B T_n + \tfrac{1}{2} m\,|\mathbf{U}|^2\right)
                  + \Gamma_{\rm rec}\!\left(\tfrac{3}{2} k_B T_i + \tfrac{1}{2} m\,|\mathbf{V}|^2\right)
\end{aligned}
}
$$

with $m \equiv m_i = m_n$ for hydrogen. The $-\Gamma_{\rm ion}\chi_H$ term in $Q^{e_i}_{\rm ion}$ is what the plan calls "the χ_H drain on the electron thermal pool"; it has no counterpart in $Q^{e_n}_{\rm ion}$ because under the optically-thin recombination assumption the 13.6 eV released at every recombination escapes as a Lyman-continuum photon. This convention is identical to Leake et al. (2012), Khomenko et al. (2014), Maneva et al. (2017), and Snow et al. (2020) — it is the standard.

**Sanity check on the bookkeeping.** Sum the energy sources:

$$
Q^{e_i}_{\rm ion} + Q^{e_n}_{\rm ion} = -\Gamma_{\rm ion}\chi_H.
$$

So total non-radiative fluid energy decreases at rate $\Gamma_{\rm ion}\chi_H$, which is exactly what should leave the system as escaping recombination photons in steady state (because at steady state $\Gamma_{\rm ion} = \Gamma_{\rm rec}$, each ionization is matched by a recombination, and each recombination radiates one Lyman-continuum photon of energy $\chi_H$). This is the right asymptote and matches Snow et al. (2020) §2.

---

# 3. Modified conservative vector form (your §4.3, Eq. ~62)

Your existing source vector $\mathbf S$ on line 553 of main.tex needs **two new entries on the mass rows, two on the momentum rows, two on the energy rows** — six additions total:

$$
\mathbf{S} = 
\underbrace{\begin{bmatrix} 0 \\ 0 \\ \cdots \\ \cdots \\ \cdots \\ \cdots \end{bmatrix}}_{\text{your existing}}
\;+\;
\begin{bmatrix}
+m\,\Gamma_{\rm net} \\[2pt]
-m\,\Gamma_{\rm net} \\[2pt]
+m\,\Gamma_{\rm ion}\,U - m\,\Gamma_{\rm rec}\,V \\[2pt]
-m\,\Gamma_{\rm ion}\,U + m\,\Gamma_{\rm rec}\,V \\[2pt]
Q^{e_i}_{\rm ion} \\[2pt]
Q^{e_n}_{\rm ion}
\end{bmatrix}
$$

where in 1.5D field-aligned form $V$ and $U$ are scalars (signed along $\hat{\mathbf b}$), so $\tfrac{1}{2}m|\mathbf U|^2 \to \tfrac{1}{2}m U^2$ in $Q^{e_i}_{\rm ion}$ and $Q^{e_n}_{\rm ion}$. The flux vector on line 542 of main.tex, the eigenvalues (§4.4), the Rusanov solver, and the GLM block are all unchanged — ionization introduces no new wave, only a stiff body source.

These six rows belong in $R_I$ (your implicit residual), not $R_E$. The plan's argument for this is correct and uses exactly the same logic as your existing "Why not Picard" paragraph: in the transition region the effective relaxation rate $\lambda_{\rm ion} \equiv n_e(S_i + \alpha_r)$ can satisfy $\lambda_{\rm ion}\Delta t \gg 1$ at the acoustic CFL, so Picard diverges and explicit treatment is paralyzed. Backward Euler on a scalar quadratic is the right choice.

---

# 4. Stage E: what the local solve actually looks like

The plan's §3.4 reduction to a scalar ODE in the ionization fraction $f \equiv \rho_i/\rho_{\rm tot}$ is correct and is the key implementation simplification. With $\rho_{\rm tot}$, $V_{\rm cm} = (\rho_i V + \rho_n U)/\rho_{\rm tot}$, $w = V - U$, and the species *total* energy frozen during the Stage-E sub-step, the local ODE is

$$
\frac{df}{dt} = n_{\rm tot}\,\big[\,f(1-f)\,S_i(T_e) - f^2\,\alpha_r(T_e)\,\big]
$$

with $n_{\rm tot} = \rho_{\rm tot}/m$. Backward Euler at lagged $T_e$ gives a quadratic in $f^{n+1}$:

$$
f^{n+1} - \Delta t\, n_{\rm tot}\big[f^{n+1}(1-f^{n+1})S_i - (f^{n+1})^2 \alpha_r\big] = f^{n}.
$$

This expands to

$$
\Delta t\, n_{\rm tot}(S_i + \alpha_r)\,(f^{n+1})^2 + (1 - \Delta t\, n_{\rm tot} S_i)\,f^{n+1} - f^n = 0,
$$

which has a closed-form positive root. The unconditional stability claim in the plan is correct: $\partial G/\partial f|_{f_{\rm eq}} = -n_{\rm tot}S_i < 0$, and the quadratic has exactly one root in $[0,1]$ for any $\Delta t > 0$ as long as the rate coefficients are non-negative.

---

# 5. Four issues to fix or clarify before landing the math

These are not show-stoppers — the plan is fundamentally right — but you should address them in the writeup.

**(a) The Voronov coefficient.** The plan writes
$$S_i = 2.91\times10^{-14}\,U^{0.39}\,e^{-U}/(0.232 + U) \quad \text{m}^3/\text{s}$$
which is Voronov (1997) H I in SI with $U = \chi_H/(k_B T_e)$. The original is $0.291\times10^{-7}\,U^{0.39}\,e^{-U}/(0.232+U)\,\text{cm}^3/\text{s}$; the unit conversion is correct ($0.291\times10^{-7}\,\text{cm}^3 = 2.91\times10^{-14}\,\text{m}^3$). Worth citing the exact Voronov constants $A$, $K$, $X$, $P$, $\Phi$ for traceability in §6.3 of the writeup as you propose.

**(b) The Hummer coefficient and case-A vs case-B.** The plan gives $\alpha_r = 2.7\times10^{-19}(T_e/10^4)^{-0.75}~\text{m}^3/\text{s}$, which is a power-law fit to Hummer (1994) case-B. The choice of case-B (no return to ground state, photons escape) is *consistent* with the optically-thin assumption in your $Q^{e_n}$ derivation — but you should state this consistency explicitly in §6.3, because using case-A here while assuming optically-thin Lyman-α elsewhere is a frequent silent bug in the literature. Snow et al. (2023) review this in detail.

**(c) The kinetic-equilibrium ≠ Saha note in §3.5 needs to be in the paper, not just the plan.** Reviewers will ask why your steady-state $f$ doesn't match the Avrett & Loeser (2008) tabulated values, and the answer ("kinetic equilibrium under optically thin radiation differs from full NLTE radiative transfer") needs to live in the writeup. Leenaarts et al. (2007) is the right reference for the size of this discrepancy.

**(d) Splitting Δn_i into Γ_ion^Δ and Γ_rec^Δ for the energy bookkeeping.** This is a real implementation subtlety not fully resolved in §3.4 of the plan. After solving the scalar quadratic you know $\Delta n_i = (f^{n+1} - f^n)\,n_{\rm tot}$, the net. But the energy update needs the two pieces separately because $Q^{e_i}_{\rm ion}$ has *different* coefficients on $\Gamma_{\rm ion}$ and $\Gamma_{\rm rec}$ (one carries $T_n$, the other carries $T_i$; one carries $\tfrac{1}{2}mU^2$, the other carries $\tfrac{1}{2}mV^2$; one is multiplied by $\chi_H$, the other is not). The natural prescription, consistent with backward Euler, is to evaluate the rates at the *post-step* densities:

$$
\Gamma_{\rm ion}^{\Delta} = \Delta t \cdot n_e^{n+1}\, n_n^{n+1}\, S_i(T_e^n),
\qquad
\Gamma_{\rm rec}^{\Delta} = \Delta t \cdot (n_e^{n+1})^2\, \alpha_r(T_e^n).
$$

You can verify $\Gamma_{\rm ion}^{\Delta} - \Gamma_{\rm rec}^{\Delta} = \Delta n_i / m$ as a check; round-off mismatches should be redistributed in proportion. This should be explicit in §3.4.

There is also a positivity concern the plan flags but doesn't fully address: in regions where $\Gamma_{\rm ion}\chi_H \Delta t$ exceeds the available ion thermal energy $\tfrac{3}{2} \cdot 2 n_i k_B T_i$, naive subtraction drives $T_e < 0$. The plan's "floor $p_i$ and document the regime" is fine for a first PR, but you should be honest in §3.4 that this floor is non-conservative — the energy that gets floored away is real energy lost to the unphysical regime, and the only correct fix is to add the radiative-heating term that physically balances the χ_H drain in the real Sun. Snow et al. (2020) make this point.

---

# 6. A note on the literature framing for the writeup

For the new §2.3 the plan proposes, the canonical citation chain is Meier & Shumlak (2012, *Phys. Plasmas*) for the form of the source terms → Khomenko et al. (2014, *Phys. Plasmas*) for the electron-folded reduction matching your §4.2 → [Popescu Braileanu & Keppens (2022, *A&A*)](https://consensus.app/papers/details/98952205afe956eda0781aa38c0cbe3e/?utm_source=claude_desktop) for the IMEX implementation directly analogous to yours → [Leake et al. (2012, *ApJ*)](https://consensus.app/papers/details/8a28d2d4c7a8566ebb39e1990218781f/?utm_source=claude_desktop) for the chromospheric application. For the optically-thin χ_H bookkeeping, [Snow et al. (2020)](https://consensus.app/papers/details/a4b1fe37bc6350a38ae0c90d811a6b7e/?utm_source=claude_desktop) is the cleanest analytical treatment. For non-equilibrium-vs-Saha, Carlsson & Stein (2002) and Leenaarts et al. (2007). The plan's reference list is correct; I would just add Snow et al. 2020 explicitly as the source-of-record for the ionization-potential-energy-sink discussion, since it is the analytic backbone for the χ_H bookkeeping you're adopting.

---

# Bottom line

Land the plan. The math is right, the IMEX placement is right, the references are right. Before you commit the writeup edits, fix issues (a)–(d): nail down the rate-coefficient citations (a, b), put the kinetic-vs-Saha caveat in the paper not just the plan (c), and make the $\Gamma_{\rm ion}^{\Delta}/\Gamma_{\rm rec}^{\Delta}$ split explicit (d). Everything else — Stage E placement, the implicit-residual classification, the operator-split order B→C→D→E, the closed-form quadratic solve, the conservation tests — is correct.