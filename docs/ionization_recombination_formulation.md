# Hydrogen Ionization & Recombination — Formulation

**Scope.** The physically complete formulation of the hydrogen ionization/recombination block for the 1.5D field-aligned two-fluid chromosphere model. It supersedes §3.1–§3.5 of [`ionization_plan.md`](ionization_plan.md), which used an optically-thin / coronal approximation (collisional ionization + radiative recombination only).

**Bottom line.** The one structural addition is **three-body (collisional) recombination** $e+e+\mathrm{H^+}\to\mathrm H+e$, scaling as $n_e^2 n_i$. It makes the equilibrium interpolate between the coronal/kinetic balance at low density and **Saha/LTE at the dense base** — a limit a radiative-only model can never reach.

**Conventions.** SI throughout. Pure hydrogen, quasi-neutral ($n_e=n_i$), single mass $m\equiv m_i=m_n$. The folded charged fluid carries one temperature $T_e=T_i$.

## Contents

1. [The physical picture: two regimes](#1-the-physical-picture-two-regimes)
2. [Microphysics: processes & rate coefficients](#2-microphysics-processes--rate-coefficients)
3. [Equilibrium limits](#3-equilibrium-limits)
4. [Where each channel dominates](#4-where-each-channel-dominates-in-the-chromosphere)
5. [Two-fluid source terms & energy bookkeeping](#5-two-fluid-source-terms-with-full-energy-bookkeeping)
6. [Local evolution & Stage E solve](#6-local-evolution-and-point-implicit-stage-e-solve)
7. [Validation tests](#7-validation-tests)
8. [Changes from `ionization_plan.md`](#8-changes-from-ionization_planmd)
9. [References](#9-references)

## 0. Notation

| Symbol | Meaning |
|---|---|
| $f \equiv \rho_i/\rho_\text{tot}$ | ionization fraction |
| $n_\text{tot}=\rho_\text{tot}/m$ | total number density |
| $n_e=n_i=f\,n_\text{tot}$ | electron / ion density |
| $n_n=(1-f)\,n_\text{tot}$ | neutral density |
| $T_e=T_i$ | charged-fluid temperature |
| $\chi_H=13.6$ eV | hydrogen ionization energy |
| $U\equiv\chi_H/k_BT_e$ | ionization energy in units of $k_BT_e$ |
| $\mathbf V,\ \mathbf U$ | ion, neutral bulk velocity (§5 only) |
| $S_i(T_e)$ | collisional (electron-impact) ionization coefficient |
| $P_\text{phot}$ | effective photoionization rate $[\mathrm{s^{-1}}]$ — frozen-field $R_{ik}(z)$ (§2.2) |
| $\alpha_r(T_e)$ | radiative recombination coefficient |
| $\alpha_c(T_e,n_e)$ | three-body / collisional-radiative recombination coefficient |
| $\alpha_\text{tot}=\alpha_r+\alpha_c$ | total effective recombination coefficient |
| $\Phi(T_e)$ | Saha function (§3.2) |

---

## 1. The physical picture: two regimes

The chromosphere is partially ionized, with $f$ running from $\lesssim10^{-4}$ at the temperature minimum to $\sim1$ at the top of the transition region. One formulation must bridge two regimes:

- **Dense lower chromosphere — collision-dominated, near LTE.** Ionization is set by *Saha equilibrium*. Reproducing Saha requires collisional ionization to be balanced by its exact inverse, **three-body recombination** $e+e+\mathrm{H^+}\to\mathrm H+e$ (scaling as $n_e^2 n_i$) — the term that was **absent** from the previous formulation.
- **Rarefied middle/upper chromosphere — out of equilibrium.** The ionization/recombination timescale exceeds the dynamical timescale, so $f$ is set by history, not Saha (Carlsson & Stein 2002; Leenaarts et al. 2007). Radiative recombination dominates and the coronal-type balance of `ionization_plan.md` is appropriate.

The previous formulation relaxed to the kinetic balance $f/(1-f)=S_i/\alpha_r$ at **all** heights — density-independent, so it can never equal Saha. Adding three-body recombination makes the equilibrium interpolate: coronal at low density, **Saha at high density** (§3.3 explains why a single-level model reaches Saha only approximately).

### 1.1 "Dense" means high $n_e$, driven by high $n_\text{tot}$

The three-body rate is controlled by the **electron/ion density** $n_e=n_i$ (two electrons and an ion must meet), not $n_\text{tot}$ directly. But $n_e=f\,n_\text{tot}$, and the largest $n_e$ sits at the *weakly ionized* dense base because $n_\text{tot}$ is so large there:

| Layer | $T$ (K) | $n_\text{tot}$ (m⁻³) | $f \approx n_e/n_\text{tot}$ | $n_e$ (m⁻³) |
|---|---|---|---|---|
| Lower chromosphere | 5180 | ~$10^{21}$–$10^{22}$ | ~$10^{-3}$ | $6.5\times10^{18}$ |
| Upper chromosphere | 6220 | ~$10^{18}$–$10^{19}$ | ~$2\times10^{-2}$ | $7.5\times10^{16}$ |

So the base is more neutral in fraction yet has ~90× the electron density. Separately, $n_\text{tot}$ enters the Saha equilibrium value $f^2/(1-f)=\Phi/n_\text{tot}$ (higher $n_\text{tot}\to$ lower $f$), which is why the dense base is weakly ionized in the first place. Net: **$n_e$ controls the rate competition; $n_\text{tot}$ sets the equilibrium.**

### 1.2 Corrections in this revision

Two errors in earlier drafts were caught and fixed here, both verified against the primary papers:

1. **The three-body coefficient is a collisional-radiative (CR) fit, not the ground-state $S_i/\Phi$.** Detailed balance holds *level by level*, and the process proceeds overwhelmingly through capture into high-$n$ Rydberg levels — so the ground-state-only estimate undershoots by ~$5\times10^4$. The correct single-level coefficient sums the capture-cascade over all levels (Hinnov & Hirschberg 1962; Stevefelt, Boulmer & Delpech 1975), which is why a single-level model recovers Saha only approximately (§2.4, §3.3).
2. **Neutral third bodies are not negligible at the dense base.** When $n_n\gg n_e$ (here $n_n/n_e\sim10^3$), neutral-stabilized recombination into excited levels can *exceed* the electron-stabilized process (Drawin 1968), and the molecular channel comes within a factor of a few of radiative (Cretu 2022; Nóbrega-Siverio 2019). These are now carried as potentially-comparable terms, not omitted (§2.5).

---

## 2. Microphysics: processes & rate coefficients

Ionization channels first (§2.1–§2.2), then recombination channels (§2.3–§2.5). Each subsection closes with the coefficient this model uses.

### 2.1 Collisional (electron-impact) ionization — ground-state `S_i` and multilevel `S_CR`

Two distinct collisional ionization rates matter, and conflating them is the ionization-side counterpart of the three-body error (§1.2): the ground-state rate is negligible, but the *level-summed* rate is the process that pulls the dense base to Saha.

**(a) Direct ground-state ionization** $\mathrm{H}(1) + e \to \mathrm{H^+} + 2e$, with $\Gamma_\mathrm{ion}^{\,\mathrm{col,1}} = n_e\, n_n\, S_i(T_e)$ and the Voronov (1997) fit currently in [`physics.hpp`](physics.hpp):

$$
S_i(T_e) = 2.91\times10^{-14}\;\frac{U^{0.39}\,e^{-U}}{0.232 + U}\ \ [\mathrm{m^3\,s^{-1}}],
\qquad U \equiv \frac{\chi_H}{k_B T_e},\quad \chi_H = 13.6\ \mathrm{eV}.
$$

As a *coefficient* it is density-independent and valid wherever electrons are thermal. But at chromospheric $T_e$ (5000–6000 K), $U\sim26$–31, so $S_i\propto e^{-U}$ is exponentially tiny: direct ground-state ionization is **negligible** in the chromosphere proper and only takes over in the transition region. (The actual chromospheric ground-state ionizer is Ly$\alpha$ excitation + Balmer-continuum *photo*ionization — §2.2; Carlsson & Stein 2002.)

**(b) Multilevel / collisional-radiative ionization** $\mathrm{H}(n) + e \to \mathrm{H^+} + 2e$ summed over the excited-level ladder, with $\Gamma_\mathrm{ion}^{\,\mathrm{col,CR}} = n_e\, n_n\, S_\mathrm{CR}(T_e,n_e)$. This is the **exact inverse of three-body recombination** (§2.4): just as capture proceeds through high-$n$ Rydberg levels, collisional ionization proceeds by *stepwise* excitation up the ladder and over the lowered continuum, where the per-step thresholds $\chi_n=\chi_H/n^2\to0$ carry no Boltzmann penalty. The level-summed coefficient is therefore ~$10^4$–$10^5\times$ the ground-state $S_i$ and, like $\alpha_c$, **grows with $n_e$** (more electrons drive more stepwise excitation). It is the collisional ionizer that balances three-body recombination to give Saha at the dense base — the process the §4.0 summary table labels "multilevel / Rydberg-ladder."

*Coefficient.* The level-resolved rates are given by **Johnson (1972)** (his eq. 39, the Maxwellian-averaged $S_e(n,c)$ for H, used by Leenaarts et al. 2007) and **Vriens & Smeets (1980)** (closed-form analytic rate coefficients for ionization from any level $n$, plus the matching three-body recombination); the level-summed CR coefficient $S_\mathrm{CR}(T_e,n_e)$ was first tabulated by **Bates, Kingston & McWhirter (1962)** and rests on the high-$n$ ladder rates of **Mansbach & Keck (1969)** — the same rates underlying the Stevefelt $\alpha_c$ already adopted in §2.4. In a single-level model the cleanest implementable form follows by detailed balance from the three-body coefficient already in hand (electrons are Maxwellian, so the pair balances *level by level*, §3.3):

$$
\boxed{\;S_\mathrm{CR}(T_e,n_e) = \alpha_c(T_e,n_e)\,\frac{\Phi(T_e)}{n_e}\;}
\qquad\Longrightarrow\qquad
\text{the pair }S_\mathrm{CR}\!\leftrightarrow\!\alpha_{3b}\text{ has Saha (§3.2) as its exact equilibrium.}
$$

This requires **no new data** — it reuses $\alpha_c$ (Hinnov/Stevefelt, §2.4) and the Saha function $\Phi$ (§3.2), and is exactly the detailed-balance closure offered as "option 2" in §3.3.

**Use.** Carry the total collisional ionization $S_i^\mathrm{eff} = S_i + S_\mathrm{CR}$: ground-state Voronov $S_i$ always (negligible except in the TR) **plus** the multilevel $S_\mathrm{CR}=\alpha_c\Phi/n_e$, which makes the high-density equilibrium exactly Saha (§3.3 option 2). Both are *collisional* ionization and so both drain $\chi_H$ from the electron pool (§5); in LTE the $S_\mathrm{CR}$ drain exactly cancels the three-body superelastic heating, as detailed balance requires.

### 2.2 Photoionization — `P_phot`

$$
\mathrm{H} + \gamma \;\longrightarrow\; \mathrm{H^+} + e,
\qquad
\Gamma_\mathrm{ion}^{\,\mathrm{phot}} = n_n\, P_\mathrm{phot}.
$$

Effective ground-state rate $[\mathrm{s^{-1}}]$ for the **dominant chromospheric ionization channel** — Ly$\alpha$ excitation + Balmer-continuum ($n=2\to$ continuum) photoionization, *not* direct Lyman-continuum absorption (Carlsson & Stein 2002). Its energy comes from the radiation field, so it carries **no** $\chi_H$ drain on the gas.

**The photon-flux closure.** This is where the photon flux enters the model — we do not solve radiative transfer. Carlsson & Stein (2002) show the chromospheric photoionizing field is nearly time-independent: its radiation temperature varies only ~1–2% (40–100 K out of ~6000 K) even through shocks, because the Balmer-continuum photons come from the slowly-varying deep layers. The radiative ionization rate can therefore be **frozen** into a function of height alone — the "fixed radiative-rate" recipe (Sollum 1999; Leenaarts et al. 2007; Leenaarts 2020). We adopt the FAL-C $R_{ik}(z)$ tabulation of **Chae (2021, Table 1)**, who inverts the FAL-C radiation field (Fontenla et al. 1993) under exactly this assumption, and set $P_\mathrm{phot}(s)=R_{ik}(z(s))$ by log-interpolating their 43-point table ($z=-100\to2110$ km) onto the grid.

The profile is strongly height-structured: $R_{ik}\approx7\times10^{-8}\ \mathrm{s^{-1}}$ at the temperature minimum (weakest radiation field), rising to $\sim7\times10^{-3}\ \mathrm{s^{-1}}$ at the top of the chromosphere — and with it photoionization tracks radiative recombination to within a factor of a few across the chromosphere (the radiation-driven balance, visualized in [`visualization/plot_c7_ioniz_recomb_rates.py`](../visualization/plot_c7_ioniz_recomb_rates.py)). Two caveats: (1) $R_{ik}$ is a FAL-C rate applied to the sibling C7 stratification and calibrated to a 1D *mean* atmosphere, so it is not self-consistent in strongly magnetic or shocked gas; (2) the frozen field omits the downward *coronal* EUV flux — but this is **benign for hydrogen**: H photoionization is Balmer-continuum dominated and the Lyman continuum (the only coronal-EUV route to H ionization, since the corona emits nothing at the Balmer edge) "does not play any significant role" (Carlsson & Stein 2002). Coronal EUV instead matters as upper-chromosphere/TR *heating* (handled by thermal conduction from the corona + the optional Carlsson & Leenaarts 2012 recipe) and as a *helium* ionizer (out of scope for this H-only treatment; Golding et al. 2016 for the self-consistent binned-EUV method).

**Use.** For `model_c7`, $P_\mathrm{phot}$ is the **C7-context inversion**: solve the local ionization balance of the full network ($S_i$, $S_\mathrm{CR}$, $\alpha_r$, $\alpha_c$) for the photo channel on C7's tabulated $(n_p, n_e, n_n, T)$, so C7 is an equilibrium fixed point of our case-B network in the dense lower/mid chromosphere, $P_\mathrm{phot}=\frac{n_p}{n_n}(\alpha_r+\alpha_c)n_e-(S_i+S_\mathrm{CR})n_e$ with $n_p=n_\mathrm{H}-n_\mathrm{H1}$ (AvrettLoeser 2008, Table 26). Blend above $z\approx1500$ km to the frozen FAL-C $R_{ik}(z)$ of Chae (2021, Table 1; log-interpolated, geometrically diluted above the chromosphere top per their Eq. 26). The inverted rate runs ~2–17× above Chae's *ground-state* $R_{ik}$ — the expected case-B-vs-ground-state recombination + C7-vs-FAL-C atmosphere difference, **not** metal contamination (metals are only ~5–9% of $n_e$, from $n_p=n_\mathrm{H}-n_\mathrm{H1}$; the earlier "30–50%" figure is not supported by Table 26). Chae's $R_{ik}$ is the NEQ-region frozen rate + an independent sanity comparison, not the target. Fallback: a uniform $10^{-4}\ \mathrm{s^{-1}}$ for scenarios without a tabulated atmosphere.

### 2.3 Radiative recombination — `α_r(T_e)`

$$
\mathrm{H^+} + e \;\longrightarrow\; \mathrm{H} + \gamma,
\qquad
\Gamma_\mathrm{rec}^{\,\mathrm{rad}} = n_e^2\, \alpha_r(T_e).
$$

Hummer (1994) Case-B fit, the closure in [`physics.hpp`](physics.hpp) — a *total* (level-summed) coefficient accurate ~10% over $10^3$–$10^5$ K:

$$
\alpha_r(T_e) = 2.7\times10^{-19}\left(\frac{T_e}{10^4\ \mathrm{K}}\right)^{-0.75}\ [\mathrm{m^3\,s^{-1}}].
$$

The binding energy leaves as an escaping Lyman-continuum photon (optically thin). Cross-check: Hinnov & Hirschberg (1962) give $\alpha_R \approx 2.7\times10^{-13}(k_BT_e[\mathrm{eV}])^{-3/4}\ \mathrm{cm^3\,s^{-1}} = 3.0\times10^{-19}\ \mathrm{m^3\,s^{-1}}$ at $10^4$ K — same law, within 10%.

**Use.** Hummer $\alpha_r$ as given.

### 2.4 Three-body (collisional-radiative) recombination — `α_c(T_e, n_e)` *(corrected)*

$$
\mathrm{H^+} + e + e \;\longrightarrow\; \mathrm{H}(n) + e \;\xrightarrow{\text{cascade}}\; \mathrm{H}(1) + e,
\qquad
\Gamma_\mathrm{rec}^{\,\mathrm{col}} = n_e^2\, n_i\, \alpha_{3b} = n_e\, n_i\,\big[\alpha_{3b}\,n_e\big].
$$

A second electron carries off the excess energy. The decisive fact: capture is into **high-$n$ Rydberg levels**, where the threshold $\chi_n = \chi_H/n^2 \to 0$ gives large cross-sections and no Boltzmann penalty, followed by collisional/radiative cascade to the ground state (Bates, Kingston & McWhirter 1962; McWhirter & Hearn 1963). The *total* coefficient is therefore dominated by the excited-level ladder and is ~$5\times10^4$ larger than the ground-state-only contribution.

**Why the naive $\alpha_{3b}=S_i/\Phi$ (ground state) fails.** Detailed balance holds *level by level*: $\alpha_{3b}(n) = S_i(n)/\Phi_n$, with $\Phi_n$ the Saha factor for level $n$ (ionization energy $\chi_n=\chi_H/n^2$). The total $\sum_n \alpha_{3b}(n)$ is dominated by high $n$. Using ground-state Voronov $S_i$ with the full Saha function $\Phi$ ($\chi_1=13.6$ eV) keeps only the $n{=}1$ term — the single *smallest* contribution. At 5180 K, $n_e=6.5\times10^{18}$: $S_i/\Phi$ gives $\alpha_c\approx2.6\times10^{-23}\,\mathrm{m^3/s}$ (ratio $6\times10^{-5}$ to $\alpha_r$), whereas the true CR coefficient gives $\alpha_c\approx1.4\times10^{-18}\,\mathrm{m^3/s}$ (ratio $\approx2.8$). The discrepancy is the missing Rydberg ladder.

**Coefficient to use** — a closed-form CR fit that already sums the capture-cascade. Two options, both giving $\alpha_c\propto n_e$:

- **Hinnov & Hirschberg (1962)**, valid $0.25<k_BT_e<1$ eV (≈2900–11600 K — covers the chromosphere):
  $$
  \alpha_c(T_e,n_e) \approx 5.6\times10^{-27}\,\big(k_BT_e[\mathrm{eV}]\big)^{-9/2}\,n_e[\mathrm{cm^{-3}}]\ \ [\mathrm{cm^3\,s^{-1}}].
  $$
- **Stevefelt, Boulmer & Delpech (1975)**, full CR coefficient (radiative + coupling + collisional), derived for $T_e<4000$ K so *extrapolated* at chromospheric $T$:
  $$
  \alpha_\mathrm{CR} = \underbrace{1.55\times10^{-10}T^{-0.63}}_{\text{radiative}}+ \underbrace{6.0\times10^{-9}T^{-2.18}n_e^{0.37}}_{\text{CR coupling}}+ \underbrace{3.8\times10^{-9}T^{-4.5}n_e}_{\text{pure 3-body}}\ [\mathrm{cm^3\,s^{-1}}],
  $$ . 
  $T$ in K, $n_e$ in cm⁻³; the last term matches Hinnov's $\alpha_c$ within a factor ~3.

The **total effective recombination coefficient** is then

$$
\boxed{\;\alpha_\mathrm{tot}(T_e, n_e) = \alpha_r(T_e) + \alpha_c(T_e, n_e)\;}
$$

— the CR coefficient of Bates, Kingston & McWhirter (1962), growing linearly (or, via the Stevefelt coupling term, sub-linearly) with $n_e$.

**Use.** Hinnov $\alpha_c$ (baseline) or the Stevefelt three-term fit; add to $\alpha_r$.

### 2.5 Neutral & molecular channels — *not negligible at the dense base*

When a heavy neutral, rather than a second electron, removes the excess energy, there are two channels. An earlier draft dismissed both as subdominant; the primary papers show that is wrong for the weakly-ionized chromosphere.

**Channel 1 — neutral-stabilized electron-ion recombination** $e + \mathrm{H^+} + \mathrm{H} \to \mathrm{H} + \mathrm{H}$, rate $\propto n_e\,n_i\,n_n$ (Bates 1965; Drawin 1968). Drawin's coefficient follows by detailed balance from the Thomson atom-atom ionization cross-section [his Eqs. (10)–(11)]:

$$
Q_{a\text{-}a} = \frac{g_i}{2g_+}\left(\frac{m_A+m_e}{m_e m_A}\right)^{3/2}
\frac{h^3}{(2\pi k_B T)^{3/2}}\,\langle\sigma v\rangle_{a\text{-}a}\,e^{+E_i/k_B T},
$$

and — exactly as in the electron case (§2.4) — the $e^{-E_i/k_BT}$ inside $\langle\sigma v\rangle_{a\text{-}a}$ cancels, so $Q_{a\text{-}a}$ is a *smooth* function of $T$. Drawin notes (p. 416) that when $n_a\gg n_e$, neutral stabilization can govern the populations of the highly-excited levels.

**Quantified (this work).** For *any* third body $X$, $Q_X(i)=\langle\sigma v\rangle_{\rm ion,X}(i)/\Phi_i$ (the Saha factor $\Phi_i$ is independent of $X$), so the neutral/electron ratio is just the ratio of *ionization* rate coefficients:

$$
\frac{Q_{a\text{-}a}}{Q_{a\text{-}e}} = \frac{\langle\sigma v\rangle_{a\text{-}a}}{\langle\sigma v\rangle_{a\text{-}e}}.
$$

Evaluating at the base (5180 K, ground state, $\omega_i=30.5$): the electron rate is Voronov $S_i=2.1\times10^{-28}\,\mathrm{m^3/s}$; the atom-atom rate from Drawin Eq. (8) — with $\xi_i=1$, $m_e/m_H=5.4\times10^{-4}$, $(k_BT/\pi m_H)^{1/2}=3.69\times10^3$ m/s, $32\pi a_0^2=2.8\times10^{-19}\,\mathrm{m^2}$, and $\Psi_{m_H}(\omega_i)=(1+2/\omega_i)e^{-\omega_i}$ (verified against Drawin's Table, p. 415) — is $\langle\sigma v\rangle_{a\text{-}a}=3.5\times10^{-32}\,\mathrm{m^3/s}$. Hence $Q_{a\text{-}a}/Q_{a\text{-}e}\approx1.7\times10^{-4}$, and with $n_n/n_e\approx10^3$ at the base,

$$
\frac{\nu_{\rm neutral}}{\nu_{ee}}=\frac{n_n}{n_e}\,\frac{Q_{a\text{-}a}}{Q_{a\text{-}e}}\approx 0.17,
\qquad
\alpha_n \equiv Q_{a\text{-}a}\,n_n \approx 2.3\times10^{-19}\,\mathrm{m^3/s}\ (\approx\tfrac12\alpha_r).
$$

So neutral-stabilized recombination is **~15–20% of the electron three-body rate** (≈ half radiative) at the dense base — a real correction, not negligible, but not dominant. This is likely a **lower bound**: it is computed from the *ground-state* ionization rate, whereas the geometric $\propto n^4$ Rydberg–neutral cross-section makes the neutral channel relatively more effective for the high-$n$ levels that dominate the cascade; a firm value needs a multilevel CR calculation with both electron *and* neutral collision rates. *(Consistent with Gousset et al. 1977.)*

**Channel 2 — molecular (ion-atom-atom) route** $\mathrm{H^+}+\mathrm{H}+\mathrm{H}\to\mathrm{H_2^+}+\mathrm{H}$ (rate $\propto n_i\,n_n^2$; Pérez-Ríos 2018; Cretu 2022), or radiative association $\mathrm{H^+}+\mathrm{H}\to\mathrm{H_2^+}+\gamma$, then fast dissociative recombination $\mathrm{H_2^+}+e\to2\mathrm{H}$. Cretu et al. give $k_3(T)\sim10^{-31}\,\mathrm{cm^6\,s^{-1}}$ at 3000 K, falling to $\sim10^{-32}$ at $10^4$ K, with a *smoother* $T$-dependence than the electron process ($\propto T^{-9/2}$) — so it gains relative importance as $T$ rises. Scaled to the dense base ($T\approx5200$ K, $k_3\approx4\times10^{-32}\,\mathrm{cm^6\,s^{-1}}$, $n_n\approx6.5\times10^{15}\,\mathrm{cm^{-3}}$), the per-ion rate $\nu_\mathrm{mol}=k_3 n_n^2\approx\mathbf{1.8\ \mathrm{s^{-1}}}$ — comparable to radiative ($\nu_\mathrm{rad}=\alpha_r n_e\approx2.9\,\mathrm{s^{-1}}$) and within ~5× of electron three-body ($\nu_{ee}\approx8.9\,\mathrm{s^{-1}}$). *(Order-of-magnitude; $k_3$ extrapolated slightly below Cretu's tabulated range, and the H$_2^+$ fate adds uncertainty.)*

**Direct H₂ formation.** Nóbrega-Siverio et al. (2019) include non-equilibrium $\mathrm{H_2}$ formation/dissociation (6-level NEQ hydrogen) in chromospheric Bifrost flux-emergence runs. H₂ forms in the coolest dense pockets ($T\approx3\times10^3$ K), releasing **4.48 eV per molecule** — exothermic heating that counteracts cooling — and the LTE assumption *vastly overestimates* the H₂ fraction. So the molecular channel is a genuine thermodynamic actor in cool, dense gas.

**Use / decision.** (a) Implement electron three-body now (§2.4). (b) Carry the neutral-stabilized term as an explicit, quantified *uncertainty* (Drawin Eq. 11). (c) Defer the molecular route to an H₂-chemistry module (à la Nóbrega-Siverio 2019) for cool-pocket fidelity. This is "include the electron term, flag the neutral terms as potentially comparable," **not** "omit and document."

---

## 3. Equilibrium limits

Setting $df/dt=0$ (§6) and dropping photoionization for clarity, the full rate balance is

$$
(1-f)\,S_i = f\,\big[\alpha_r + \alpha_{3b}\,f\,n_\mathrm{tot}\big].
$$

### 3.1 Low-density (coronal/kinetic) limit $\big(\alpha_c \ll \alpha_r\big)$

$$
\frac{f}{1-f} = \frac{S_i}{\alpha_r}\qquad\text{— density-independent (the old `ionization\_plan.md` balance).}
$$

### 3.2 High-density (Saha/LTE) limit $\big(\alpha_c \gg \alpha_r\big)$

If — and only if — the collisional recombination coefficient is the exact detailed-balance inverse of collisional ionization, $\alpha_{3b}=S_i/\Phi$, then

$$
(1-f)\,S_i = f^2\,\frac{S_i}{\Phi}\,n_\mathrm{tot}
\;\Longrightarrow\;
\boxed{\;\frac{f^2}{1-f} = \frac{\Phi(T_e)}{n_\mathrm{tot}}\;}
\;\Longleftrightarrow\;
\frac{n_e n_i}{n_n} = \Phi(T_e),
$$

the **Saha equation**, with the Saha function

$$
\Phi(T_e) = \left(\frac{2\pi m_e k_B T_e}{h^2}\right)^{3/2} e^{-\chi_H/k_B T_e}
\approx 2.42\times10^{21}\Big(\tfrac{T_e}{\mathrm K}\Big)^{3/2} e^{-1.578\times10^5\,\mathrm{K}/T_e}\ \ [\mathrm{m^{-3}}].
$$

The collisional rate $S_i$ cancels — equilibrium becomes pure thermodynamics, as LTE demands.

### 3.3 Why a single-level model recovers Saha only approximately

The collisional-ionization / three-body pair is in **mutual detailed balance at any density** (the electrons are Maxwellian), so $\alpha_{3b}=S_i/\Phi$ is exact *for whatever $S_i$ one uses*. The subtlety is which $S_i$:

- Exact Saha (§3.2) requires the **total, level-summed** collisional ionization to balance the **total** three-body recombination. Both proceed mainly through excited Rydberg levels, where rates are ~$10^4$–$10^5\times$ the ground-state values. Those large collisional rates are what overpower the (non-detailed-balanced) radiative pair $P_\mathrm{phot}\!\leftrightarrow\!\alpha_r$ at high density and pull the gas to Saha.
- A **single-level (ground + continuum) model omits that ladder.** Keeping ground-state Voronov $S_i$ and a CR *fit* for $\alpha_c$ (§2.4), the two are no longer mutual detailed-balance partners, so the high-density equilibrium is *near* Saha but not exact.

Three ways to handle it, in increasing fidelity:

1. **CR-fit recombination + Voronov ionization (this model's baseline).** Realistic *total* recombination (so the §4 crossover is right), approximate Saha; document the residual.
2. **Detailed-balance-enforced single level.** Add the multilevel collisional ionization $S_\mathrm{CR}=\alpha_c\,\Phi/n_e$ of §2.1b as the partner of the CR recombination $\alpha_c$, so the pair satisfies Saha *by construction*. Cheap; guarantees the LTE limit; the $S_\mathrm{CR}$ term *is* the multilevel ladder, expressed through detailed balance rather than computed level-by-level.
3. **Multilevel atom (production-grade; Johnson 1972; Leenaarts et al. 2007; Snow et al. 2023).** 5–6 levels with collisional + radiative rates and level-by-level detailed balance. Correct rates *and* exact Saha. The right target if the single-level LTE error proves significant in test §7.2.

**Recommendation:** implement option 1 now; add the option-2 toggle so the LTE limit can be enforced and the error bounded; reserve option 3 for later.

---

## 4. Where each channel dominates in the chromosphere

### 4.0 Summary: dominant ionization & recombination by layer

| Layer | $T_e$ (K) | $n_e$ (cm⁻³) | Regime | **Dominant ionization** | **Dominant recombination** | Brief reason |
|---|---|---|---|---|---|---|
| **Lower chromosphere** (dense base) | 5180 | $6.5\times10^{12}$ | Collision-dominated, near **LTE/Saha** | **Collisional (multilevel / Rydberg-ladder)** $S_\mathrm{CR}$ (§2.1b) — balances 3-body to set Saha | **Three-body (collisional-radiative)** $\alpha_c$ ($\alpha_c/\alpha_r\approx2.8$) | High $n_e$ (driven by huge $n_\text{tot}$): 3-body $\propto n_e^2 n_i$ wins; large collisional rates through excited levels overpower the radiative pair and pull gas to Saha. Neutral + molecular channels are comparable here (~20% of total) — none ignorable. |
| **Mid chromosphere** | 5030 | $7.5\times10^{10}$ | Rarefied, **non-equilibrium** | **Photoionization** (Ly$\alpha$ excitation + Balmer-continuum) | **Radiative** $\alpha_r$ ($\alpha_c/\alpha_r\approx0.04$) | $n_e$ ~90× lower → 3-body ($\propto n_e$) collapses. Direct collisional ionization exponentially suppressed ($U\sim26$–31, $S_i\propto e^{-U}$), so the radiation field drives ionization; $\chi_H$ escapes as a Lyman-continuum photon. |
| **Upper chromosphere** | 6220 | $7.5\times10^{10}$ | Rarefied, **non-equilibrium** | **Photoionization** (Ly$\alpha$ + Balmer-continuum) | **Radiative** $\alpha_r$ ($\alpha_c/\alpha_r\approx0.02$) | Same as mid: low $n_e$ kills 3-body; higher $T$ raises $n_e^\mathrm{cross}$, lowering the ratio further. Ionization set by radiative history, not Saha. |
| **Transition region** (top) | $\gtrsim10^4$, rising | rising, $f\to1$ | Coronal/kinetic | **Direct electron-impact** (Voronov $S_i$) | Radiative (shrinking as $f\to1$) | As $T$ climbs, $U=\chi_H/k_BT_e$ drops and the $e^{-U}$ penalty lifts, so thermal collisional ionization finally takes over; recombination is weak in the nearly-fully-ionized gas. |

**Organizing principle.** $n_e$ controls the *rate* competition (three-body $\propto n_e^2 n_i$ vs. radiative $\propto n_e^2$, so their ratio $\propto n_e$), while $T_e$ sets both the Saha *equilibrium* and the $e^{-U}$ suppression of direct collisional ionization. The crossover is $n_e^\mathrm{cross}\approx4.8\times10^{13}(k_BT_e[\mathrm{eV}])^{15/4}$ cm⁻³ — only the dense base sits above it, which is why three-body recombination (absent from the old radiative-only model) matters there and nowhere else.

### 4.1 Crossover detail (recombination)

Radiative/three-body crossover ($\alpha_c=\alpha_r$, Hinnov closed form): $n_e^\mathrm{cross} \approx 4.8\times10^{13}\,(k_BT_e[\mathrm{eV}])^{15/4}\ \mathrm{cm^{-3}}$. Against the VAL IIIC atmosphere (Vernazza, Avrett & Loeser 1981; values as in Snow 2023, Table 2):

| Layer | $T_e$ (K) | $n_e$ (cm⁻³) | $n_e^\mathrm{cross}$ (cm⁻³) | $\alpha_c/\alpha_r$ (Hinnov) | dominant recombination |
|---|---|---|---|---|---|
| **Lower chromosphere** | 5180 | $6.5\times10^{12}$ | $2.3\times10^{12}$ | **≈ 2.8** | **three-body** |
| Mid chromosphere | 5030 | $7.5\times10^{10}$ | $2.1\times10^{12}$ | ≈ 0.04 | radiative |
| Upper chromosphere | 6220 | $7.5\times10^{10}$ | $4.6\times10^{12}$ | ≈ 0.02 | radiative |

**Conclusion.** In the dense lower chromosphere three-body recombination is ~3× radiative (Hinnov, in its validity range); the radiative-only model underestimates total recombination by ~4× there and cannot reach Saha. Higher up it is a few-percent correction by the pure-three-body estimate. (The full Stevefelt coefficient, with its CR-coupling term, gives a *larger* base enhancement, ~5×, and hints at a tens-of-percent effect even mid-chromosphere — but Stevefelt is extrapolated above its $T<4000$ K range, so treat that as an upper bound. The robust statement: **three-body dominant at the base, modest to negligible above**; pinning down the mid-chromosphere value requires the multilevel calculation.) The fix matters precisely at the dense base — the inflow boundary, where wrong $f$ and $n_e$ contaminate pressure, the EOS, and every collisional coupling.

### 4.2 Full recombination budget at the dense base

Combining all channels (per-ion recombination frequency $\nu = $ rate per ion, at the VAL IIIC lower-chromosphere point: $T=5180$ K, $n_e=6.5\times10^{18}\,\mathrm{m^{-3}}$, $n_n\approx6.5\times10^{21}\,\mathrm{m^{-3}}$, $f\approx10^{-3}$):

| Channel | coefficient | $\nu$ (s⁻¹) | in current model? | source |
|---|---|---|---|---|
| Electron three-body | $\alpha_c=1.4\times10^{-18}\,\mathrm{m^3/s}$ | **8.9** | ✗ | Hinnov 1962 |
| Radiative | $\alpha_r=4.4\times10^{-19}\,\mathrm{m^3/s}$ | 2.9 | ✓ | Hummer 1994 |
| Molecular (H₂⁺) | $k_3 n_n^2$, $k_3\approx4\times10^{-32}\,\mathrm{cm^6/s}$ | 1.8 | ✗ | Cretu 2022 |
| Neutral-stabilized | $\alpha_n=2.3\times10^{-19}\,\mathrm{m^3/s}$ | 1.5 | ✗ | Drawin 1968 (§2.5) |
| **Total** | | **≈ 15** | **only ~19%** | |

**Headline.** The current radiative-only model captures only ~19% of the total recombination at the dense base. Electron three-body alone is the single largest channel (~60% of total); the two neutral channels together add another ~20%. All four are within a factor of ~6 of one another — none is ignorable there. (Each $\nu$ is order-of-magnitude: Hinnov is factor-2 accurate, $\nu_n$ is a likely lower bound, $\nu_\mathrm{mol}$ has $k_3$/H₂⁺-fate uncertainty.) Higher in the chromosphere ($n_e\sim10^{17}\,\mathrm{m^{-3}}$) all three collisional/neutral channels fall steeply and radiative again dominates.

---

## 5. Two-fluid source terms with full energy bookkeeping

Pure point sources (no flux divergence): they add to the RHS and leave every flux, Jacobian, eigensystem, and Riemann solver untouched. Integrated volumetric rates (all $\ge0$), split by *energy fate*:

$$
\Gamma_\mathrm{ion}  = \underbrace{n_e n_n S_i}_{\text{collisional}} + \underbrace{n_n P_\mathrm{phot}}_{\text{photo}},
\qquad
\Gamma_\mathrm{rec}  = \underbrace{n_e^2 \alpha_r}_{\text{radiative}} + \underbrace{n_e^2\,\alpha_c}_{\text{collisional (3-body)}}
                     \equiv \Gamma_\mathrm{rec}^{\,\mathrm{rad}} + \Gamma_\mathrm{rec}^{\,\mathrm{col}}.
$$

Keeping $\Gamma_\mathrm{rec}^{\,\mathrm{rad}}$ and $\Gamma_\mathrm{rec}^{\,\mathrm{col}}$ separate is essential — they dispose of $\chi_H$ differently.

### Mass and momentum

$$
\partial_t n_i + \nabla\!\cdot(n_i \mathbf V) = \Gamma_\mathrm{ion} - \Gamma_\mathrm{rec},
\qquad
\partial_t n_n + \nabla\!\cdot(n_n \mathbf U) = -(\Gamma_\mathrm{ion} - \Gamma_\mathrm{rec}),
$$
$$
\partial_t(n_i\mathbf V)\big|_\mathrm{src} = \Gamma_\mathrm{ion}\mathbf U - \Gamma_\mathrm{rec}\mathbf V,
\qquad
\partial_t(n_n\mathbf U)\big|_\mathrm{src} = -\Gamma_\mathrm{ion}\mathbf U + \Gamma_\mathrm{rec}\mathbf V.
$$

### Energy — with three-body superelastic heating

$$
\begin{aligned}
Q^{e_i}_\mathrm{ion} &= \Gamma_\mathrm{ion}\!\left(\tfrac32 k_B T_n + \tfrac12 m U^2\right)
   - \Gamma_\mathrm{rec}\!\left(\tfrac32 k_B T_i + \tfrac12 m V^2\right)
   - \Gamma_\mathrm{ion}^{\,\mathrm{col}}\,\chi_H
   \;\boxed{+\;\Gamma_\mathrm{rec}^{\,\mathrm{col}}\,\chi_H}\\[4pt]
Q^{e_n}_\mathrm{ion} &= -\Gamma_\mathrm{ion}\!\left(\tfrac32 k_B T_n + \tfrac12 m U^2\right)
   + \Gamma_\mathrm{rec}\!\left(\tfrac32 k_B T_i + \tfrac12 m V^2\right)
\end{aligned}
$$

The two $\chi_H$ terms:

- $-\Gamma_\mathrm{ion}^{\,\mathrm{col}}\chi_H$: **collisional** ionization pays $\chi_H$ from the electron pool. Photoionization does *not* appear (its energy came from radiation), so the drain uses $\Gamma_\mathrm{ion}^{\,\mathrm{col}}$ only.
- $+\Gamma_\mathrm{rec}^{\,\mathrm{col}}\chi_H$: **collisional (three-body)** recombination is *superelastic* — $\chi_H$ is handed to the third electron and stays as heat. Radiative recombination has no such term ($\chi_H$ escapes as a photon). *(Subtlety: a Rydberg capture releases $\chi_n$ promptly, the rest during cascade; in the collision-dominated dense limit the cascade is collisional and essentially all $\chi_H$ is recycled, so $+\Gamma_\mathrm{rec}^{\,\mathrm{col}}\chi_H$ is the correct high-density limit. At low density the cascade radiates and $\Gamma_\mathrm{rec}^{\,\mathrm{col}}\to0$ anyway.)*

**Net non-radiative budget:** $Q^{e_i}_\mathrm{ion} + Q^{e_n}_\mathrm{ion} = -\Gamma_\mathrm{ion}^{\,\mathrm{col}}\chi_H + \Gamma_\mathrm{rec}^{\,\mathrm{col}}\chi_H$.

- **Dense/LTE limit:** collisional ionization balances three-body recombination, net $\to0$ — *no spurious radiative sink in equilibrium*, as detailed balance requires. (The old "radiate $\chi_H$ on every recombination" convention violated this and would have over-cooled the deep chromosphere.)
- **Rarefied limit:** $\Gamma_\mathrm{rec}^{\,\mathrm{col}}\to0$, recovering the old budget — each radiative recombination radiates one $\chi_H$ photon.

---

## 6. Local evolution and point-implicit Stage E solve

Freezing $\rho_\mathrm{tot}$, center-of-mass velocity, drift, and total energy across the substep:

$$
\frac{df}{dt} = n_\mathrm{tot}\Big[f(1-f)\,S_i - f^2\big(\alpha_r + \alpha_c\big)\Big] + (1-f)\,P_\mathrm{phot},
\qquad \alpha_c = \alpha_c(T_e, n_e),\ n_e=f n_\mathrm{tot}.
$$

Because $\alpha_c$ depends on $n_e$ (linearly for Hinnov, or as $n_e^{0.37}$ for the Stevefelt coupling term), this is **nonlinear** in $f$ beyond the previous quadratic. Recommended solve:

**(A) Lagged-coefficient quadratic + Picard (primary).** Evaluate $\alpha_\mathrm{eff} = \alpha_r + \alpha_c(T_e, n_e^{\mathrm{old}})$, then the backward-Euler step is the **existing closed-form quadratic**:

$$
A f^{2} + B f - C = 0,\quad
A = \Delta t\,n_\mathrm{tot}(S_i + \alpha_\mathrm{eff}),\;
B = 1 - \Delta t\,n_\mathrm{tot} S_i + \Delta t\,P_\mathrm{phot},\;
C = f^{\mathrm{old}} + \Delta t\,P_\mathrm{phot}.
$$

Update $n_e\leftarrow f^{\mathrm{new}}n_\mathrm{tot}$, recompute $\alpha_\mathrm{eff}$, repeat 1–2 Picard iterations. The $n_e$-dependence is mild over one step, so this converges fast and reuses the current solver.

**(B) Scalar Newton on the full residual (fallback for stiff dense cells).** Solve $R(f)=f - f^{\mathrm{old}} - \Delta t\,[\,\dots\,]=0$ directly, including the exact $\alpha_c(T_e, f n_\mathrm{tot})$ dependence.

Either way the step is unconditionally stable ($\partial(df/dt)/\partial f<0$ at the root; all loss processes push $f$ toward the fixed point; exactly one root in $[0,1]$ for $\Delta t>0$).

**If enforcing the LTE limit (option 2 of §3.3):** also add the multilevel term on the collisional-ionization side, $S_i^\mathrm{eff} = S_i^\mathrm{Voronov} + S_\mathrm{CR}$ with $S_\mathrm{CR}=\alpha_c\,\Phi/n_e$ (§2.1b), so the collisional pair is detailed-balanced and the high-density root is exactly Saha. The matching $\chi_H$ drain in §5 then uses $\Gamma_\mathrm{ion}^{\,\mathrm{col}}=n_e n_n(S_i+S_\mathrm{CR})$, and its LTE value cancels the three-body heating $\Gamma_\mathrm{rec}^{\,\mathrm{col}}\chi_H$ to zero net.

### Energy update after the $f$-solve

Split the net $\Delta n_i$ into channels at post-step densities:

$$
\Gamma_\mathrm{ion}^{\,\mathrm{col},\Delta} = \Delta t\, n_e^{\mathrm{new}} n_n^{\mathrm{new}} S_i,\quad
\Gamma_\mathrm{ion}^{\,\mathrm{phot},\Delta} = \Delta t\, n_n^{\mathrm{new}} P_\mathrm{phot},\quad
\Gamma_\mathrm{rec}^{\,\mathrm{rad},\Delta} = \Delta t\,(n_e^{\mathrm{new}})^2 \alpha_r,\quad
\Gamma_\mathrm{rec}^{\,\mathrm{col},\Delta} = \Delta t\,(n_e^{\mathrm{new}})^2 \alpha_c.
$$

(Check $\sum\Gamma_\mathrm{ion}^\Delta - \sum\Gamma_\mathrm{rec}^\Delta = \Delta n_i$; redistribute round-off proportionally.) Apply $Q^{e_i},Q^{e_n}$ of §5, including the $-\Gamma_\mathrm{ion}^{\,\mathrm{col},\Delta}\chi_H$ drain and the $+\Gamma_\mathrm{rec}^{\,\mathrm{col},\Delta}\chi_H$ heating. Floor $p_i,p_n$ at a small positive value if the net drain exceeds available thermal energy; log it (non-conservative).

---

## 7. Validation tests

1. **Conservation.** Uniform box, no fluxes: $\rho_i+\rho_n$, total momentum, and total non-radiative energy $e_i+e_n+\chi_H(n_i-n_i^0)-\chi_H\!\int\!\Gamma_\mathrm{rec}^{\,\mathrm{col}}dt$ conserved to round-off (three-body heating recycles $\chi_H$; only radiative recombinations remove it).
2. **Saha limit (the key new test).** Hold $T_e$ fixed at high density ($n_e\gtrsim10^{19}\,\mathrm{m^{-3}}$), evolve to steady state, confirm $f^2/(1-f)\to\Phi(T_e)/n_\mathrm{tot}$. With option 1 of §3.3 expect *near* Saha (quantify the residual); with option 2 expect exact Saha. This is the test the previous formulation could not pass.
3. **Kinetic limit.** Low density ($n_e\lesssim10^{16}\,\mathrm{m^{-3}}$): confirm $f/(1-f)\to S_i/\alpha_r$ (three-body negligible; recovers old behavior).
4. **VAL IIIC stratification.** Initialize from VAL IIIC / Model C7, evolve to steady state; confirm $f(s)$ tracks the tabulated ionization fraction better in the lower chromosphere than the radiative-only model (which over-ionizes there).

---

## 8. Changes from `ionization_plan.md`

| Aspect | Old (`ionization_plan.md`) | New (this document) |
|---|---|---|
| Recombination | radiative only, $\alpha_r$ | $\alpha_\mathrm{tot}=\alpha_r+\alpha_c(T_e,n_e)$ |
| Three-body coeff. | — | **CR fit** (Hinnov 1962 / Stevefelt 1975), *not* ground-state $S_i/\Phi$ |
| Collisional ionization | ground-state Voronov $S_i$ only (negligible) | $+$ multilevel $S_\mathrm{CR}=\alpha_c\Phi/n_e$ (§2.1b; Johnson 1972, Bates et al. 1962, Vriens & Smeets 1980, Mansbach & Keck 1969) — the ionizer that reaches Saha |
| Equilibrium | kinetic everywhere | Saha at high $n_e$ (approx. in single level), kinetic at low $n_e$ |
| LTE recovery | impossible | option 1 approximate / option 2 exact / option 3 multilevel |
| $\chi_H$ energy | drained on ioniz., lost on every recomb. | drained on collisional ioniz.; **recycled** by collisional recomb. |
| Net loss in LTE | spurious $\Gamma\chi_H$ sink | $\to0$ (detailed balance) |
| Neutral third body | — | *potentially comparable* at dense base (Drawin 1968; Cretu 2022); carried as quantified uncertainty + future H₂ module, not omitted |
| Stage E solve | quadratic | lagged-$\alpha_\mathrm{eff}$ quadratic + Picard, or scalar Newton |
| Lower-chromosphere recomb. | underestimated ~4× | correct (three-body dominant) |

New physical inputs: the Saha function $\Phi(T_e)$ and a CR recombination fit ($\alpha_c$). The three-body coefficient is **not** simply the ground-state $S_i/\Phi$.

---

## 9. References

**Three-body recombination & multilevel collisional-radiative ionization (electron third body)**
- **Bates, Kingston & McWhirter (1962)**, *Proc. R. Soc. A* 267, 297 — founding CR theory; tabulates both the CR recombination coefficient $\alpha_\mathrm{CR}(T_e,n_e)$ and the CR **ionization** coefficient $S_\mathrm{CR}(T_e,n_e)$ (§2.1b), each growing with $n_e$. [`supporting-papers/BatesKingstonMcWhirter1962.pdf`]
- **McWhirter & Hearn (1963)**, *Proc. Phys. Soc.* 82, 641 — high-density Saha-Boltzmann → low-density radiative transition. [`supporting-papers/McWhirterHearn1963.pdf`]
- **Hinnov & Hirschberg (1962)**, *Phys. Rev.* 125, 795 — closed-form e-e-ion coefficient $\alpha_c\propto n_e T^{-9/2}$. [`supporting-papers/HinnovHirschberg1962.pdf`]
- **Stevefelt, Boulmer & Delpech (1975)**, *Phys. Rev. A* 12, 1246 — three-term CR recombination fit (radiative + coupling + collisional), built on the Mansbach & Keck ladder rates. [`supporting-papers/Stevefelt1975.pdf`]
- **Mansbach & Keck (1969)**, *Phys. Rev.* 181, 275 — classical Monte-Carlo trajectory rates for the high-$n$ ladder of hydrogenlike atoms; steady-state collisional recombination $\alpha=2.0\times10^{-27}n_e\,(k_BT_e[\mathrm{eV}])^{-9/2}\,\mathrm{cm^3/s}$ (§IV). This is the **radiationless / high-density limit** (their §IV assumes radiative transitions negligible; the authors flag radiative cascade and neutral collisions as significant otherwise) — hence used not bare but via the radiation-corrected Stevefelt (1975) fit it underlies. [`supporting-papers/MansbachKeck1969.pdf`]
- **Vriens & Smeets (1980)**, *Phys. Rev. A* 22, 940 — closed-form Maxwellian rate coefficients for electron-impact ionization/(de)excitation from *excited* levels of **atomic hydrogen**, plus the matching three-body recombination; the practical multilevel $S_\mathrm{CR}$ ingredients (§2.1b). Stated accuracy ~10% for $k_BT_e\gtrsim0.05\,E_{pn}$ — i.e. valid across the chromospheric $T_e$ for the excited-level ladder ($E_{pn}=\chi_H/n^2$). [`supporting-papers/VriensSmeets1980.pdf`]
- **Johnson (1972)**, *ApJ* 174, 227 — level-resolved H collisional ionization rate $S_e(n,c)$ (eq. 39) and excitation rates; collisional recombination from detailed balance. The level-resolved input for the multilevel ladder (used by Leenaarts et al. 2007). [`supporting-papers/Johnson1972.pdf`]

**Neutral third body / molecular channel**
- **Bates (1965)**, "Recombination of positive ions and electrons in a dense neutral gas" — neutral atoms as third bodies. [`supporting-papers/Bates1965.pdf`]
- **Drawin (1968)**, *Z. Phys.* 211, 404 — detailed-balance rate for e-ion three-body recombination with a *neutral* third body, into ground and excited states. [`supporting-papers/Drawin1968.pdf`]
- **Gousset et al. (1977)**, *Phys. Rev. A* 16, 2517 — neutral-stabilized *atomic*-ion recombination measured $<10^{-28}\,\mathrm{cm^6\,s^{-1}}$ (small). *(to download)*
- **Pérez-Ríos et al. (2018)**, *Phys. Rev. A* — universal $T$-dependence of ion-neutral-neutral three-body recombination $A^+{+}A{+}A\to A_2^+{+}A$. *(to download)*
- **Cretu et al. (2022)**, *Phys. Rev. A* — ion-atom-atom recombination in cold H/D plasmas; $\mathrm{H^+}{+}\mathrm{H}{+}\mathrm{H}\to\mathrm{H_2^+}{+}\mathrm{H}$. [`supporting-papers/Cretu2022.pdf`]
- **Nóbrega-Siverio et al. (2019)**, *A&A* — non-equilibrium $\mathrm{H_2}$ formation/dissociation in the chromosphere (Bifrost); the molecular channel that matters in cool pockets. [`supporting-papers/NobregaSiverio2019.pdf`] *(most chromosphere-relevant)*

**Ionization coefficients & chromospheric context**
- **Voronov (1997)**, *ADNDT* 65, 1 — collisional ionization $S_i$. [`supporting-papers/Voronov1997.pdf`]
- **Hummer (1994)**, *MNRAS* 268, 109 — radiative recombination $\alpha_r$. [`supporting-papers/Hummer1994.pdf`]
- **Carlsson & Stein (2002)**, *ApJ* 572, 626 — chromospheric ionization via Ly$\alpha$+Balmer; non-equilibrium timescales; photoionizing field varies only ~1–2% (the basis for the frozen-field $P_\mathrm{phot}$). [`supporting-papers/CarlssonStein2002.pdf`]
- **Chae (2021)**, *J. Astron. Space Sci.* 38, 83 — frozen-field photoionization rate $R_{ik}(z)$ inverted from a model atmosphere; the $P_\mathrm{phot}$ closure adopted here (Table 1, FAL-C; §2.2). [`supporting-papers/Chae2021.pdf`]
- **Sollum (1999)**, MSc thesis, Univ. Oslo — "Dynamic hydrogen ionization"; the fixed radiative-rate recipe (frozen radiation field) later adopted by Leenaarts et al. (2007). *(thesis; not in repo)*
- **Leenaarts et al. (2007)**, *A&A* 473, 625 — non-equilibrium H ionization, 2D; Johnson (1972) + Sollum (1999) rates. [`supporting-papers/Leenaarts2007.pdf`]
- **Leenaarts (2020)**, *Living Rev. Solar Phys.* 17, 3 — review of radiation–gas energy-exchange approximations (fixed-rate photoionization, coronal-irradiation heating). [`supporting-papers/Leenaarts2020.pdf`]
- **Carlsson & Leenaarts (2012)**, *A&A* 539, A39 — recipes for chromospheric radiative cooling and for heating by incident coronal radiation (the coronal-EUV channel $P_\mathrm{phot}$ omits). [`supporting-papers/CarlssonLeenaarts2012.pdf`]
- **Golding et al. (2016)**, *ApJ* 817, 125 — non-equilibrium H/He ionization with self-consistent coronal-EUV + Ly$\alpha$ photoionization (binned). [`supporting-papers/Golding2015.pdf`]
- **Leenaarts & Wedemeyer-Böhm (2006)**, *A&A* 460, 301 — time-dependent H ionization, 3D. [`supporting-papers/LeenaartsWedemeyer2006.pdf`]
- **Snow, Druett & Hillier (2023)**, *MNRAS* 525, 4717 — two-fluid chromospheric shocks, multilevel H, collisional + radiative ionization/recombination via the LTE Saha population ratio. [`supporting-papers/Snow2023.pdf`]
- **Le, Cambier et al. (2016)**, *Phys. Plasmas* 23 — multifluid ionization + three-body recombination with detailed balance, recovering Saha in the thermal limit. [`supporting-papers/LeCambier2016.pdf`]

**Framework & reference atmosphere**
- **Meier & Shumlak (2012)**, *Phys. Plasmas* 19, 072508; **Khomenko et al. (2014)**, *Phys. Plasmas* 21, 092901 — two-fluid source-term form and electron-folded reduction.
- **Vernazza, Avrett & Loeser (1981)**, *ApJS* 45, 635 (VAL III) [`supporting-papers/Vernazza1981.pdf`]; **Fontenla, Avrett & Loeser (1993)**, *ApJ* 406, 319 (FAL-C — the model underlying Chae's $R_{ik}$); **Avrett & Loeser (2008)**, *ApJS* 175, 229 (Model C7). [`supporting-papers/AvrettLoeser2008.pdf`]
