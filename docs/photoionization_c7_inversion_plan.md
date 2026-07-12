# Plan: close the C7 over-ionization gap with Chae (2021)'s own method

## What Chae (2021) actually does (verified from the paper, §2, Eqs. 11–29)

**Notation.** The photoionization rate we implement is $P_\mathrm{phot}$ (physics.hpp; formulation §2.2). Chae writes the same quantity as $R_{ik}$, so $P_\mathrm{phot}\equiv R_{ik}$ — we use $P_\mathrm{phot}$ from here on, and keep Chae's $R_{ki}, C_{ik}, C_{ki}$ as his symbols for the recombination/collisional rates (which we map to $\alpha_r, S_i, \alpha_c$ below).

$P_\mathrm{phot}$ is *physically* the radiation-field photoionization rate — Chae's Eq. 14, $P_\mathrm{phot}=\int_{\nu_0}^{\infty}(4\pi\alpha_\nu/h\nu)\,J_\nu\,d\nu$, a pure functional of the radiation field $J_\nu$. But he does **not** compute it by a radiative-transfer solve for $J_\nu$; instead he *infers* it by *inverting the hydrogen ionization-equilibrium balance* against a tabulated 1-D atmosphere. His two-level (ground $i$ + continuum $k$) equilibrium is Eq. 24:

$$n_i\,(P_\mathrm{phot} + C_{ik}) = n_k\,(R_{ki} + C_{ki})$$

with four rates, all *analytic functions of $T,n_e$ except* $P_\mathrm{phot}$:

- $P_\mathrm{phot}$ — photoionization (the unknown), per ground-state atom.
- $R_{ki} = R_{ki,\mathrm{sp}} = 2.5\times10^{-6}\,T^{-3/2}\,e^{h\nu_0/kT} E_1(h\nu_0/kT)\,n_e$ — **radiative recombination to ground** (cgs; his Eq. 20; stimulated term Eq. 23 dropped because $e^{-h\nu/kT}\ll1$).
- $C_{ik} = 5.9\times10^{-11}\,T^{1/2}\,e^{-h\nu_0/kT}\,n_e$ — **collisional (ground-state) ionization** (Cox 2000; his Eq. 21).
- $C_{ki} = 2.45\times10^{-26}\,n_e^2/T$ — **three-body recombination** (detailed balance with $C_{ik}$; his Eq. 22).

He solves for the one unknown (Eq. 25):

$$\boxed{\;P_\mathrm{phot}(z) = \frac{n_k(z)}{n_i(z)}\big[R_{ki}(T,n_e) + C_{ki}(T,n_e)\big] - C_{ik}(T,n_e)\;}$$

evaluated on FAL-C's tabulated $(T(z), n_e(z), n_k(z)/n_i(z))$. Table 1 of the paper is this inversion. For $z>z_t$ (above chromosphere top) he extrapolates by geometric dilution (Eq. 26).

## What $P_\mathrm{phot}$ physically is — and why an equilibrium-derived rate is valid in NEQ gas

**Terminology — which "equilibrium" we mean.** Throughout this plan "equilibrium" means **ionization (statistical) equilibrium**: ionization rate equals recombination rate, so $f$ is steady, $dn_i/dt=0$,

$$\frac{dn_i}{dt}=\underbrace{(P_\mathrm{phot}+S_i n_e+S_\mathrm{CR}n_e)\,n_n}_{\text{ionization}}-\underbrace{(\alpha_r+\alpha_c)\,n_e\,n_i^+}_{\text{recombination}}=0.$$

This is *not* thermal equilibrium (Maxwellian velocities at $T$) and *not* LTE/Saha ($f$ a function of $T$ alone — Chae notes the Saha equation does not apply to the chromosphere). It is *weaker* than Saha: it only requires the rates to balance. In the chromosphere proper that balance is **radiation-dominated** (photoionization $\approx$ radiative recombination), not Saha; only at the deep, dense base does it approach Saha.

**C7 is the same kind of model as FAL-C.** Both are static, 1-D, semi-empirical *non-LTE statistical-equilibrium* models of the average quiet Sun (C7 = Avrett & Loeser 2008 is the modern successor to FAL-C = Fontenla et al. 1993; Chae himself calls them near-twins, "intended to describe the same atmosphere"). So Chae's inversion is atmosphere-agnostic — running Eq. 25 on C7's tabulated $(T,n_e,n_k/n_i)$ is exactly *"Chae's method applied to C7,"* with the same validity as his published FAL-C table, **provided the coefficients we invert with are right** (the consistency requirement below).

**$P_\mathrm{phot}$ is a radiation-field rate, not a gas-state quantity.** Eq. 14 defines it as a pure functional of $J_\nu$ — it depends only on the radiation field, never on the local ionization state. Because the full non-LTE RT solve for $J_\nu$ is expensive, Chae instead *reads off* $P_\mathrm{phot}$ by inverting the equilibrium balance on a static model. **The equilibrium assumption is only a measuring device**; the quantity it extracts is still the radiation-field rate.

**This dissolves the apparent contradiction** of using an "equilibrium-derived" rate in a region we call non-equilibrium. The frozen $P_\mathrm{phot}$ is the radiation-field *drive*; the non-equilibrium lives in the gas's *response* — its slow recombination lag — not in $P_\mathrm{phot}$. The same $P_\mathrm{phot}$ legitimately drives both equilibrium gas (dense base: relaxes to balance) and NEQ gas (upper chromosphere: lags, stays over-ionized), because it is a property of the photons, not of $f$. Chae does exactly this himself: §3.2 applies the FAL-C $P_\mathrm{phot}$ to plasma features with *arbitrary* $p,T$ (not the equilibrium state), calling the result (Eq. 29) "the non-LTE version of the Saha equation." Its justification (his intro; Carlsson & Stein 2002) is that the photoionizing field is *non-local* — set by the slowly-varying deep layers — so it varies only ~1–2% even through shocks and can be frozen as a function of height alone.

**The honest caveat.** The frozen-field *assumption itself* is an approximation, marginal in strongly dynamic gas. Chae's relaxation-time analysis (Eq. 33, Fig. 6) gives $\tau=1/(P_\mathrm{phot}+C_{ik}+R_{ki}+C_{ki})\sim80$–$100$ s in the upper chromosphere — comparable to dynamical times — and he concludes ionization equilibrium "is not fully justified" for dynamic features around 12,000 K (Carlsson & Stein 2002). So the genuine limitation is the *frozen-field/equilibrium-extraction* in shocked plasma, **not** any inconsistency in using $P_\mathrm{phot}$ as a NEQ drive.

## Why our current state has a gap, and what the old "unphysical" profile really was

Our **rejected** C7 calibration was $P=\alpha_r n_e^2/n_n - S_i n_e$. In Chae's notation (with $n_k=n_i^+$, $n_i=n_n$, all neutrals in ground state so $n_i\approx n_n$):

$$P_\text{old} = \frac{n_i^+}{n_n}\,\big(\alpha_r n_e\big) - S_i n_e \;=\; \text{Eq. 25 with } C_{ki}=0,\ R_{ki}\to\alpha_r n_e,\ C_{ik}\to S_i n_e.$$

So **the old calibration was Chae's Eq. 25 with the three-body term dropped and our (Hummer case-B / Voronov) coefficients substituted.** Two reasons it looked unphysical:

1. **Top of the column (ill-conditioned inversion + incomplete network).** At C7's domain top ($h\approx1989$ km, $T\approx6674$ K, $n_e/n_n\approx1.4$) the old inversion returned $P_\mathrm{phot}\approx(n_e/n_n)\alpha_r n_e\approx2.4\times10^{-2}\ \mathrm{s^{-1}}$ — ~10× Chae's FAL-C value (~$2.5\times10^{-3}$), hitting the old $10^{-2}$ cap. The cause is *not* "C7 being out of equilibrium" — C7 is a static SE model, it has no history. Two real causes: (i) **numerical conditioning** — as $f\to1$, $n_n\to0$, so $P_\mathrm{phot}\propto n_e^2/n_n$ is ill-conditioned and amplifies any error in C7's tabulated $n_e/n_n$; (ii) the **incomplete/mismatched network** (case-B $\alpha_r$, no three-body — see point 2). The old code then *froze* this inflated rate while $n_n$ kept dropping, driving $f\to1$ (runaway). Independently, the *real* upper chromosphere is genuinely NEQ (Chae Eq. 33 / Fig. 6: $\tau\sim10^2$ s $\gtrsim$ dynamical time; Carlsson & Stein 2002), so we don't *want* to pin the gas to C7 there anyway — which is why the fix is to use the independent FAL-C rate there, not to invert C7.
2. **Coefficient mismatch.** Chae's $R_{ki}$ is recombination *to the ground level only*; our $\alpha_r$ is Hummer **case-B** (sum over $n\ge2$), ~2× larger — so even in equilibrium our inversion runs ~2× high.

## The honest verdict

- **Lower/mid chromosphere (dense, in equilibrium):** C7 *is* in ionization equilibrium by construction, so inverting the **complete** Eq. 25 (i.e. adding back the three-body $C_{ki}$ + collisional $C_{ik}$ the old version dropped) yields a physical $P_\mathrm{phot}$ and makes C7 an *exact fixed point of whatever rate network we invert with*. **The gap closes here, by construction, physically.**
- **Upper chromosphere (rarefied, gas is NEQ):** the *real* gas here is out of ionization equilibrium (slow recombination lag), so the residual drift is **real over-ionization physics**, not an error — and we *want* the gas free to depart from C7 here. Forcing closure (pinning C7 as a fixed point) is exactly what blew up the old calibration. **The gap should *not* be closed here.** Keep the frozen FAL-C **radiation-field** rate (the current Chae-table approach) — i.e. the radiation-field $P_\mathrm{phot}$ extracted via Chae's equilibrium inversion of a static reference model, used here as the frozen ionizing *drive* while the gas responds out of equilibrium (see "What $P_\mathrm{phot}$ physically is" above).

So: **yes, re-evaluating $P_\mathrm{phot}$ in the C7 context with Chae's method closes the gap — but only in the equilibrium region, which is where it should close.**

## The consistency requirement (the crux)

An inversion makes C7 a fixed point **only of the rate network used in the inversion.** Our Stage-E solver currently uses $\{S_i^\text{Voronov},\ \alpha_r^\text{Hummer},\ P_\mathrm{phot}\}$ with **no three-body**. So to actually close the gap in the *running* code, the Stage-E network and the inversion network must be the same set. Two routes:

### Route A — Chae-faithful (replicate his 2-level model exactly)
Adopt Chae Eqs. 20–22 ($R_{ki}, C_{ik}, C_{ki}$) as the Stage-E coefficients, and set $P_\mathrm{phot}$ by inverting Eq. 25 on C7. C7 becomes an *exact* fixed point of the 2-level equilibrium; inverted $P_\mathrm{phot}$ should match Chae's FAL-C Table 1 up to the FAL-C↔C7 atmosphere difference. Cost: replaces Voronov/Hummer with Cox-2000/ground-recomb coefficients (a step *down* in recombination fidelity — ground-state only, no Rydberg ladder, no case-B).

### Route B — our-network-consistent (recommended)
Keep Voronov $S_i$ + Hummer $\alpha_r$, **add three-body $\alpha_c$** (Hinnov/Stevefelt — already coded in `visualization/plot_c7_ioniz_recomb_rates.py`) and its detailed-balance partner $S_\mathrm{CR}=\alpha_c\Phi/n_e$ (formulation §2.1b / §3.3 option 2), then set $P_\mathrm{phot}$ by inverting our network on C7:

$$P_\mathrm{phot}(z) = \frac{n_i^+}{n_n}\Big[\alpha_r(T)\,n_e + \alpha_c(T,n_e)\,n_e\Big]\;-\;\Big[S_i(T) + S_\mathrm{CR}(T,n_e)\Big]\,n_e$$

This closes the gap at the dense base by construction **and** keeps the better recombination physics (CR-enhanced three-body → correct Saha limit), which the formulation doc independently argues for. It is "Chae's method" (his Eq. 25 inversion + equilibrium assumption) expressed in our coefficient convention. Note: at Saha the $\alpha_c$/$S_\mathrm{CR}$ pair cancels in the inversion, so the base $P_\mathrm{phot}$ reduces to the radiative residual $\alpha_r\Phi - S_i n_e$ — a small, well-defined number, *not* the inflated old value.

**Recommendation: Route B**, because the three-body + $S_\mathrm{CR}$ addition is already on the roadmap (formulation §3.3), it keeps Voronov/Hummer/Hinnov, and it is the physically complete network.

## Proposed method (Route B), step by step

1. **Add three-body recombination to Stage-E.** Port `alpha_c_hinnov(T, n_e)` from the plot script into `physics.hpp` as `recombination_rate_alpha_c(T, n_e)`. Add the detailed-balance collisional ionization $S_\mathrm{CR}=\alpha_c\,\Phi(T)/n_e$ (needs the Saha function $\Phi(T)$ — add `saha_phi(T)` to `physics.hpp`).
2. **Define the C7 inversion in the equilibrium region.** In `model_c7.cpp`, replace `photoionization_rate_chae(h_cell)` with the Eq.-25 inversion above, evaluated on the cell-centered C7 table values `T_E, ne_E, nn_E` (ion = `ne_E`, neutral = `nn_E`). Clip at $P_\mathrm{phot}\ge0$ only.
3. **Blend to the frozen FAL-C rate in the NEQ region.** Define the equilibrium↔NEQ crossover by the relaxation time $\tau=1/(P_\mathrm{phot}+C_{ik}+R_{ki}+C_{ki})$ vs the cell sound-crossing time (or just use Chae's $z\approx1500$ km marker). Below it: use the inversion. Above it: use `photoionization_rate_chae(h)` (the frozen radiation-field rate). Smooth blend (e.g. $\tanh$ in $\log\tau$) across the transition. **Document that the upper-chromosphere residual is retained on purpose** (real NEQ over-ionization).
4. **Update the Stage-E solver** (`src/integrators.cpp::apply_ionization_stage`). Three-body recomb $\propto f^3$ turns the backward-Euler $f$-update from a quadratic into a **cubic** $A f^3 + B f^2 + C f - D = 0$; solve with one or two Newton iterations seeded from the current quadratic root (cheap, point-local). $S_\mathrm{CR}$ ionization enters as a $(1-f)$ term like photoionization. The $\chi_H$ energy drain must include the $S_\mathrm{CR}$ collisional branch and the three-body super-elastic heating (formulation §5 / line 377 — they cancel at LTE).
5. **Validate** (this is what makes it physical, not a fudge):
   - Plot inverted-C7 $P_\mathrm{phot}(z)$ against Chae's FAL-C Table 1 on the existing figure. Expect agreement to ~10–30% through the equilibrium chromosphere; quantify and document the discrepancy (FAL-C↔C7 + case-B-vs-ground-recomb).
   - Confirm the inflated top is gone (no cap needed) and that the blend hands off smoothly to the frozen rate.
6. **Re-run the smoke test:** C7 IC should now be a near-fixed-point in the lower/mid chromosphere (drift → 0 there), with the small, *intended* NEQ relaxation surviving only in the upper chromosphere.
7. **Update the fixed-point test** (`tests/chromo_tests.cpp::test_model_c7_photoionization_uses_chae_frozen_field_rate`): in the equilibrium region assert $|f_a-f_b|\to0$ (now a legitimate fixed point); in the NEQ region assert only the bounded-drift / no-runaway condition.

## What this does and does NOT buy

- **Does:** make C7 an (approximate) fixed point where it physically *is* in equilibrium; replace the borrowed FAL-C-on-C7 rate with a self-consistent C7-context rate; complete the rate network to reach the correct Saha limit at the base; give a validatable $P_\mathrm{phot}$ (compare to Chae Table 1).
- **Does NOT:** make the model "correct" in dynamic/shocked gas — the frozen $P_\mathrm{phot}(z)$ is itself an approximation there ($\tau\sim10^2$ s $\approx$ dynamical time; Chae Eq. 33 / Fig. 6 / §3.2; see "What $P_\mathrm{phot}$ physically is" above). And it deliberately leaves the upper-chromosphere NEQ over-ionization in place, because that is real physics, not a gap to close.

## Files touched (for execution)
- `physics.hpp` — add `recombination_rate_alpha_c`, `saha_phi`, `S_CR`; keep `photoionization_rate_chae` (now used only above the blend).
- `scenarios/model_c7.cpp` — replace the `photoionization_rate_chae` population with the Eq.-25 inversion + blend.
- `src/integrators.cpp` — cubic Stage-E update + $\chi_H$ energy bookkeeping for the new collisional channels.
- `tests/chromo_tests.cpp` — split the assertion by equilibrium vs NEQ region.
- `visualization/plot_c7_ioniz_recomb_rates.py` — overplot inverted-C7 $P_\mathrm{phot}$ vs Chae FAL-C Table 1 (validation).
- `docs/ionization_recombination_formulation.md` §2.2 + writeup `main.tex`/`paper.tex` §3.1 — document the C7-context inversion and the equilibrium/NEQ split.
