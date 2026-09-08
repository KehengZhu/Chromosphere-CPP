# Publication strategy — terminating the AWSoM-R thread

**Written 2026-09-06. Rewritten the same day** around the correct axis: this work exists to serve
Igor Sokolov's AWSoM-R, and the gap it fills is stated in AWSoM-R's own paper. The first draft
organized the field into four literature families and treated the coupling target as one of three
options; that buried the actual story. Everything else is context around this axis.

**Second revision, 2026-09-06 (same day).** Four questions were put to the plan — *does the missing
Alfvén-wave physics break it; why does the host want separate electron and ion temperatures; why is
evaporation the requested run; and if this is "an AWSoM-R-consistent module for later coupling," is
that publishable at all* — and answering them from the literature changed three things that are now
woven through the document rather than appended to it:

1. **A harder motivation than §1.4's quoted sentence, now verified against primary sources.** AWSoM
   pins its chromospheric inner boundary a factor of ten above the same group's own chromospheric-top
   value, and the group's stated purpose for that number ranges across its papers from *preventing*
   chromospheric evaporation to *enabling* it to *replenishing what it removes* (§1.4, "The harder
   evidence"). The host does not merely lack a chromosphere model — it carries a numerical
   placeholder in place of one, and the placeholder has a reported observable cost. This is the
   strongest material available for the Introduction, and it is stronger *because* of the
   inconsistency, not in spite of it. **And the sole authority for the placeholder being harmless —
   Lionello et al. (2009) — turns out to rest on a factor-2 test used to license a factor-10-to-250
   extrapolation, and to limit its own validity explicitly to T > 500,000 K, an order of magnitude
   above AWSoM's inner boundary and above the temperatures where the error is reported** (§1.4
   F1–F5). That reduces the programme's largest risk: the case no longer requires Tier 0 to return a
   large number.
2. **An eighth assumption, and a required companion to §5.1.** `T_e = T_i` is the host's own
   documented efficiency shortcut in exactly our region (§1.3 row 8, §2.7); and a Carlsson &
   Leenaarts radiative sink introduced *without* a background heating source is not a
   simplification but a wrong model (§2.6c, §5.1).
3. **An explicit statement of what the contribution is** and the test that separates it from a
   software deliverable (§2.5), with venues attached to deliverables rather than to papers.

This is a *forward-looking planning document*, not a study recap and not a statement of code
behaviour. Current code is the source of truth for implementation; `paper.tex` is the source of
truth for what is currently written down.

**Plain-language companion: `project_context_zh.md`** — a Chinese explanatory document covering what
AWSoM-R is, where the opening is, and what this work delivers, written for a reader who has not read
Sokolov et al. (2021). It is derived from this plan and is *not* authoritative: all planning,
sequencing and risk decisions live here.

**Animated companion:
`/Users/zkeheng/Developer/CHOMBO/HYPERS-Shock/learning/06-chromosphere-awsom-gap/`** — an 8.5-minute
narrated 1080p lecture animation of the same argument, in English, built from this plan and
`project_context_zh.md`. Its panels are plotted live from `scenarios/model_c7.cpp` and from
`outputs/model_column/twall_sweep/twall_22kK_cfl0.50.txt.gamma_diag`, so it cannot drift from the
code; `chromosphere_gap_NOTES.md` in that folder records the real-vs-schematic split, every literature
number with its source, and **ten places where the primary PDFs differ from how this document
paraphrases them** — including that the "≈1 s AWSoM-R timestep" §3.3 asks to verify against secondary
sources is in fact stated in the 2021 paper itself (p. 8). It is derived from this plan and is *not*
authoritative.

Companion documents: `writeup_standards.md`; `../studies/numerics/state_precision_release_cutover.md`;
`../studies/eos-ionization/ionization_recombination_formulation.md`; `../studies/boundaries/`.
Primary reference: `../supporting-papers/Sokolov_2021_ApJ_908_172.pdf`. Also in
`../supporting-papers/`: `Al Shidi et al. - 2019 - Time-Dependent Two-Fluid ...pdf` (§2.6d),
`FisherCanfieldMcClymont1985.pdf` (§2.5 Bin 1), `CarlssonLeenaarts2012.pdf` + `CL2012 Figures`
(§5.1), `Johnston2019.pdf` / `Johnston2020.pdf` and `NobregaSiverio2019.pdf`.

---

## 1. The gap, stated by the host model

**Sokolov et al. 2021, ApJ 908, 172** (`10.3847/1538-4357/abc000`; submitted 2017, revised 2020,
published 2021). AWSoM-**R** — the R is *realtime*, a design goal rather than a version number.

### 1.1 What AWSoM-R does

The low solar corona (R☉ < R < R_b ≈ 1.1 R☉, containing the transition region at
R☉ < R < 1.03 R☉) is removed from the 3D MHD grid and replaced by **threads**: magnetic field lines
traced through *every grid point of the global coronal model's lower boundary*, along each of which
1D transport equations are solved, interfaced to the full 3D MHD model at R_b.

The motivation is explicitly computational:

> "the greatest number of computational resources is spent in maintaining the numerical solution in
> the **low SC and in the transition region**, where the temperature gradients are sharp and the
> magnetic field topology is complicated. The degraded computational efficiency is caused by the need
> for the highest resolution as well as the use of a fully three-dimensional implicit solver for
> electron heat conduction."

> "the low SC appears to be **a bottleneck limiting the computational efficiency and performance**."

Result: "faster-than-real-time performance … on ∼200 cores."

### 1.2 The 1D justification is the host's, not ours

We do not have to construct an argument that a 1D field-aligned column is adequate here. Sokolov
already made it, and built AWSoM-R on it:

> "although the simulations of the low SC are computationally intense, **the physical nature of the
> processes involved is rather simple as long as the heat fluxes and slow plasma motional velocities
> are mostly aligned with the magnetic field** … the plasma state at any point within the low SC is
> controlled by the plasma, particle, and Alfvén wave transport **along the magnetic field line which
> passes through this point**."

This converts "isn't 1D too simple?" from a question we must defend into one we **inherit**. We make
the same approximation as the host, one step further down.

### 1.3 How the thread currently ends — eight assumptions

The TFLM lower boundary is set at the *top of the transition region* and matched to the chromosphere
through an **analytical TR model** (Lionello et al. 2001, 2009; Downs et al. 2010). The equations, as
printed on p. 6 of the paper:

```
(19)   −∂/∂s( κ₀ T^{5/2} ∂T/∂s ) = −N_e N_i Λ_R(T)

(20)   ½ κ₀² T⁵ (∂T/∂s)² |_{T_ch}^{T}  =  {N_i T}² ∫_{T_ch}^{T} κ₀ T₁^{1/2} Λ_R(T₁) dT₁

(21)   {N_i k_B T} = (1/L_TR) ∫_{T_ch}^{T_TR} κ₀T₁^{5/2} dT₁ / U_heat(T₁),
       U_heat(T)  = sqrt( (2/k_B²) ∫_{T_ch}^{T} κ₀ (T')^{1/2} Λ_R(T') dT' )

(22)   κ₀ T_TR^{5/2} (∂T/∂s)|_{T=T_TR} = {N_i k_B T} · U_heat(T_TR)
```

**Note the square in Equation (20).** Pulling {N_i T} out of the integrand turns N_e N_i into
{N_i T}²/T², which requires the ion density to track the mass density — **full ionization is
hard-coded into the closure.** At the lower limit T_ch ∼ (1÷2)×10⁴ K, Model C7 (Avrett & Loeser
2008, Table 26) gives a hydrogen ionization fraction of **x = 0.75 at 10⁴ K rising to 0.97 at
2×10⁴ K**, so the identity is in error by **3–25 % across the T_ch range** — an approximation that
degrades sharply toward the lower limit rather than a wholesale failure.

*An earlier draft of this document said "x ∼ 10⁻²–10⁻¹, so that identity is simply false there."
That was wrong: 10⁻²–10⁻¹ is C7's value at h ≈ 1000–1400 km, the mid chromosphere, not at T_ch.
Corrected 2026-09-06; see §2.7 for the table and the checking.* **The interesting consequence is not
the error's size but its variation:** a factor of ~8 in the assumption-6 error is swept out by
T_ch's own unconstrained factor of two, which makes it a second Tier 0 quadrature alongside the Λ_R
one (§5.4). There is also **no ionization energy anywhere in the system**: the energy balance is
conduction versus optically thin radiation only, and the reservoir χ_H n_H x does not appear. The
whole construction is **static** — no time derivative, no enthalpy flux, no mass flux.

Unpacked as assumptions:

| # | Assumption | Source in Sokolov (2021) |
| --- | --- | --- |
| 1 | **L_TR ≈ 1 Mm is prescribed** as the TR width along the thread | "we choose at each magnetic thread a short section of length ∼1 Mm to be a width of the TR … an important distinction from previous works, in which the temperature on top of the TR had been set, rather than the TR width" |
| 2 | **{N_i T} is constant** through the TR | "assumed to be constant; therefore, it is separated from the integrand, as the temperature gradient scale within the TR is much shorter than the barometric scale" |
| 3 | **Steady state**: conduction balanced by optically thin radiation, ∂/∂s(κ₀T^{5/2} ∂T/∂s) = −N_e N_i Λ_R | "In the steady-state version of Equation (18), we keep only the terms that dominate within the TR" |
| 4 | **The TR is heated solely from the corona; the base heat flux is zero** | "we postulate that the TR is heated solely by the heat transfer from the corona and the lower boundary of the transition region is where the heat flux from the corona turns to zero" |
| 5 | **T_ch ∼ (1÷2)×10⁴ K is given; Λ_R is a CHIANTI table** (optically thin, ionization equilibrium, coronal approximation) | Eqs. (21)–(22) and Figure 1 |
| 6 | **N_e = N_i — full ionization**, implied by the square in Eq. (20) | Eq. (20) |
| 7 | **No ionization-energy reservoir; no time derivative, enthalpy flux or mass flux** | Eqs. (19)–(22) as printed |
| 8 | **T_e = T_i — one temperature**, by choice and for cost, in exactly this region | "At higher densities as in the low corona, we assume Te=Ti and use the single-temperature energy equation for electron and ions **to improve the computational efficiency**" |

**Assumption 8 was missing from the first draft of this table and is a free addition to the audit's
scope.** Every term in Eqs. (19)–(22) is an *electron* quantity written as if there were a single
temperature: κ₀T^{5/2} is Spitzer **electron** conduction; Λ_R is a function of T_e (the 3D model's
own sink is written Q_rad = N_e N_i Λ_R(T_e)); and N_e N_i is electron × ion density. The 3D model
carries separate electron and ion energy equations and partitions turbulence dissipation between
them with an ion fraction f_p — and then reverts to one temperature precisely where we model. §2.7
gives the reasons the host wants the split back, and the measurement that says it does not yet
matter inside our 553 km domain.

**T_ch is a free parameter with a factor-of-two range**, and it is the lower limit of *both* integrals
in Eq. (21) — including U_heat, which near its lower limit is dominated by Λ_R exactly where the
optically thin CHIANTI table is least reliable. The sensitivity of {N_i k_B T} to T_ch is therefore
computable **by quadrature alone, with no simulation at all** (§5.4, Tier 0).

The output is a Robin-type condition relating T and ∂T/∂s at the top of the TR — "BCs, which set
neither temperature nor its gradient, but which relate the temperature and its gradient at the
boundary."

### 1.4 The sentence — and the harder evidence behind it

> **"There is still a minor uncertainty in the way to distinguish the TR from the top of the
> chromosphere, originating from the fact that _we do not include a consistent chromosphere model_."**
>
> "…the TR solution at some point should be merged with the chromosphere solution with no jump in
> pressure at a location that depends both on {N_I k_B T}_TR and on the pressure barometric
> distribution in the chromosphere. However, the uncertainty in this location, which also results in
> some uncertainty in L_TR, is **negligible, because the barometric scale in the chromosphere is
> small**."

**Read the "negligible" claim precisely.** It is a correct statement about the *location* of the
merge point — static and geometric. It does **not** cover:

- the time-dependent **mass flux** across that interface;
- the share of deposited energy that goes into **ionization** rather than translation;
- whether assumption 2 ({N_i T} = const) survives dynamic conditions;
- whether assumption 4 (zero base heat flux) survives a chromosphere being drained or loaded;
- whether assumption 5's optically thin CHIANTI Λ_R is valid at 1–2×10⁴ K — it is not; that is the
  optically thick chromospheric regime, and the entire reason Carlsson & Leenaarts (2012) exists.

This is a deliberately flagged opening, not an oversight. **It is the opening this work fills.**

#### The harder evidence: the chromospheric boundary is a fixed reservoir, not a model

*Verified against primary sources 2026-09-06. Quotes below are verbatim; sources and remaining gaps
are given at the end of the subsection.*

The quoted sentence says "we did not include a chromosphere model." The *stronger* statement — and
the one to lead the Introduction with — is what AWSoM does instead.

**The canonical source is van der Holst et al. (2014), the AWSoM model paper** (ApJ 782, 81;
arXiv:1311.4093, §2.5) — **not** the application papers:

> "The temperatures are all set to the same value Te = Ti = Ti∥ = T⊙ = 50,000 K uniformly at the
> inner boundary. […] **We overestimate the density for this temperature by an order of magnitude**
> with the value of Ne = Ni = N⊙ = 2×10¹⁷ m⁻³ at the inner boundary."

and immediately after:

> "**This overestimate prevents chromospheric evaporation** and extends the upper chromosphere to
> reach the correct lower density, but does not significantly change the global solution as shown in
> Lionello et al. (2009)."

**The "order of magnitude" is exactly calibrated, and the reference point is Sokolov's own.**
Sokolov et al. (2013, ApJ 764, 23; arXiv:1208.3141 §3.2.1 Eq. 36) sets the analytic chromospheric-top
boundary at `T_ch = 2×10⁴ K, N_ch ≈ 2×10¹⁶ m⁻³` — precisely a factor of ten below AWSoM's
2×10¹⁷ m⁻³. **The gap between the physical chromospheric top and the number AWSoM actually imposes
is a factor of 10, stated by the same group, in two of its own papers.** That is the sentence.

Mechanically the boundary is an effectively infinite mass reservoir standing in for the chromosphere
the model does not have. Note also that the boundary state is written with three temperatures —
T_e, T_i, T_i∥ — which is assumption 8 of §1.3 seen from the interface side.

##### The stated purpose is not consistent across the group's own papers — and that is the better argument

The first draft of this document asserted flatly that the reservoir exists "to prevent chromospheric
evaporation." **That is the majority phrasing but not the only one, and the split tracks lead author
rather than epoch.** Four distinct rationales appear in primary text for the *same* boundary value:

| Framing | Verbatim | Papers |
| --- | --- | --- |
| **Prevents / avoids** | "This overestimate prevents chromospheric evaporation"; "is required to avoid chromospheric evaporation" | van der Holst et al. 2014; Sachdeva et al. 2019, 2021; Lloveras/van der Holst 2020; and two 2025–2026 applications |
| **Allows / enables** | "This procedure **allows** chromospheric evaporation to self-consistently populate the corona with an appropriately high plasma density" | Jin et al. 2017a (arXiv:1611.08897) and 2017b (arXiv:1605.05360); still in use — Panesar et al. 2026 (arXiv:2606.17361) |
| **Replenishes** | "overestimated to provide a ready source to **replenish** the plasma, which may be depleted due to chromospheric evaporation" | Sachdeva et al. 2023 (arXiv:2212.05138) |
| **Hybrid** | "an overestimate of density at the inner boundary to prevent the **excessive** chromospheric evaporation" | Shi et al. 2022 |
| **No evaporation language at all** — purely radiative | "With this density overestimation, radiative cooling will be significant. By maintaining a temperature floor of 50,000 K the density will fall radially outward until radiative cooling is no longer able to drop the temperature below 50,000 K." | van der Holst et al. 2026, multi-ion model (arXiv:2606.19232) |

**Do not write "AWSoM suppresses evaporation" as though it were unanimous.** The accurate — and
stronger — claim is:

> **A fixed inner-boundary density, overestimated by a factor of ten relative to the same group's own
> chromospheric-top value, is carried by every AWSoM paper, and the group's own descriptions of what
> it is *for* range from "prevents evaporation" through "enables evaporation" to "replenishes plasma
> lost to evaporation." That is the signature of a numerical placeholder rather than a physical
> boundary condition — and it is precisely the object a real chromospheric column replaces.**

That formulation is safer in front of an AWSoM author, is harder to rebut, and does not require
picking a side in an inconsistency that is not ours to adjudicate.

##### The cost is reported, and the authors attribute it to two causes at once

- **Sachdeva et al. (2019)**, ApJ 887, 83; arXiv:1910.08110 §4.2 — DEMT comparison, naming *both*
  mechanisms three paragraphs apart: "from r = 1.025 R☉ up to about 1.055 R☉, the model electron
  density is overestimated compared to the DEMT results. **This overestimate is the result of
  artificial broadening of the transition region** to be consistent with our limited numerical
  resolution," and then "if the chromospheric density is too low, the transition region may
  evaporate. […] the density at the inner boundary is taken to be an overestimate, which ensures that
  the base is not affected by chromospheric evaporation. […] Thus, at low radial distances […] the
  AWSoM predicted density is still an overestimate."
- **Shi et al. (2022)** §3.4 names both in a single paragraph: "the modeled transition region still
  extends to higher altitudes and pushes the corona outwards. In addition, we use an overestimate of
  density at the inner boundary…"
- **Shi et al. (2024)**, ApJ 961, 60 states the resulting validity floor outright: the coronal
  solution "converges to physically meaningful values **above 1.01–1.03 R☉**."
- The AWSoM-R stellar paper's density overestimate at T < 10^6.2 K (§1.6).

**Consequence, and it is a real constraint on the story:** the low-height density excess has *two*
authors-acknowledged causes — the reservoir **and** Lionello's artificial TR broadening — and the
literature does not separate them. Only one of them is ours. **This cannot be resolved by reading; it
has to be measured, and Tier 0 is where.** See §2.8.

##### The "it is harmless" authority, read — and it scopes itself out of our region

**Lionello, Linker & Mikić (2009)**, ApJ 690, 902, `10.1088/0004-637X/690/1/902`, obtained
2026-09-06 and filed at `../supporting-papers/LionelloLinkerMikic2009.pdf`. This is the sole
authority every AWSoM paper cites for "the overestimate does not change the global solution," and
none of them re-derives it. It says what they say it says — **but the evidence behind it is far
narrower than the use made of it, and the paper explicitly limits its own validity to temperatures
above 500,000 K.** Four findings, all verbatim from §2.1 and Figures 1 and 3.

**F1 — The claim, and its actual evidential base (p. 904).**

> "The density at the lower boundary is set to a fixed chromospheric value that is large enough (for
> the specified heating) to produce a temperature plateau. […] **If the chromospheric density is too
> low, the transition region may evaporate. It is not crucial to estimate the required chromospheric
> density accurately; it is best to overestimate it. An overestimation will simply produce a slightly
> thicker plateau region, without affecting the overall solution.** Figure 1 demonstrates this
> property by showing one-dimensional loop solutions for different base densities. **The coronal
> parts of the solution are essentially identical.**"

**But Figure 1 tests ρ₀ = 2×10¹² cm⁻³ against 4×10¹² cm⁻³ — a factor of two.** That is the entire
demonstration. AWSoM cites it for a factor of ten (2×10¹⁷ m⁻³ against Sokolov's own 2×10¹⁶ m⁻³
chromospheric top), and Sachdeva et al. (2021) for 5×10¹⁸ m⁻³, a factor of **250** above that same
reference point. Lionello's own base is `n₀ = 2×10¹² cm⁻³, T₀ = 20,000 K` — a *different* value at a
*different* temperature from AWSoM's — and even for that he attributes the pressure-overestimate
claim onward, to Mok et al. (2005), rather than demonstrating it himself.

> **The chain is: a factor-2 test in a 1D loop → "it is best to overestimate it" → a factor-10
> boundary condition in a global model → a factor-250 one. Nobody re-tested at the larger factors.**

**F2 — Both approximations are explicitly validated only above ~500,000 K, and the paper says so
in the same breath as the claim (p. 905).** On the artificial TR broadening:

> "To broaden the transition region we need to increase κ and decrease Q at low temperatures. We have
> modified κ(T) to be constant for temperatures below T_c = 500,000 K. Accordingly, we reduce Q to
> keep κQ unchanged. […] we find that the transition region is broadened significantly (**by a factor
> of ∼400**), and yet the coronal solution remains virtually unchanged. […] We find that the emission
> of the loop in EUV is not significantly modified for coronal temperatures (above ∼500,000 K).
> **More accurate calculations at lower temperatures can be achieved by lowering the temperature at
> which the modification occurs, with the consequence of a requirement for higher resolution.**"

and the Figure 3 caption: "the emissivity is accurately reproduced by the modified model **at coronal
temperatures (above 500,000 K)**."

**AWSoM's inner boundary is at 50,000 K — one order of magnitude below Lionello's stated validity
floor — and the reported density excess is at r = 1.025–1.055 R☉ and T < 10^6.2 K, i.e. inside
exactly the range the founding paper declines to claim accuracy for.** The last sentence is not a
caveat we are imposing on them; it is theirs, and it names the remedy (lower T_c) and its price
(resolution).

**F3 — The mechanism is visible in Lionello's own Figure 3(a).** In the inset over ~5,000–20,000 km
the modified (broadened) solution sits *above* the true one: the approximation makes the transition
region both broader **and denser at a given height**. The low-height density excess AWSoM reports is
therefore not a side effect or a surprise — it is the demonstrated behaviour of the approximation,
plotted in the paper cited to establish that the approximation is harmless. The two effects also act
in the same direction, which is why they are hard to separate.

**F4 — What this does and does not settle.** It does **not** separate the reservoir from the
broadening as the cause of the density excess (§2.8) — both act the same way and Lionello tests them
in the same figures. It does something more useful: **it establishes that neither approximation was
ever validated in the region where the excess is reported, and that the founding paper says as
much.** That converts §2.8's framing from "we suspect the closure is wrong here" to "the closure's
own authority declines to cover here," which is a much easier claim to defend and does not require
attributing blame between the two mechanisms.

**F5 — And it names our deliverable.** "More accurate calculations at lower temperatures can be
achieved by lowering the temperature at which the modification occurs, with the consequence of a
requirement for higher resolution" is a cost/accuracy trade-off stated but never quantified. §4 notes
that no wall-clock cost is published anywhere for the resolved 1D alternative either. **A paper that
answers "what T_c do you actually need, and what does it cost" is directly responsive to a sentence
in the founding paper of this entire closure family** — and §3 already has the cost numbers.

##### Two corrections to the first draft

- Shi et al. (2022) is a restatement, not the origin, and writes only `Te = Tp` — cite van der Holst
  (2014) or Sachdeva (2019) for the three-temperature form.
- Sachdeva et al. (2021)'s 5×10¹⁸ m⁻³ is **25×** the canonical 2×10¹⁷, i.e. ~1.4 orders of magnitude
  — do **not** describe that particular value as "an order of magnitude."
- Sachdeva et al. (2021) cites Sokolov et al. (2013) for the evaporation rationale, but **Sokolov
  (2013) contains no statement about evaporation.** The rationale traces to van der Holst (2014)
  alone. Cite accordingly.

### 1.5 Our domain sits exactly there

The release column spans 1600–2153 km, T ≈ 6 → 22 kK. That is the top of the chromosphere and the
base of the TR — precisely Sokolov's uncertain merge region, neither more nor less.

### 1.6 Is the gap still open? — checked 2026-09-06, **yes**

**This is a time-sensitive judgement. Re-check before every submission**, and at minimum every six
months. The entire programme is void if the closure is revised upstream.

Evidence as of 2026-09-06:

- **No Sokolov first-author paper from 2022–2026 revises the chromospheric boundary condition or the
  TR treatment.** Sokolov et al. (2021) remains the canonical reference. His subsequent direction is
  **SA-MHD** (stream-aligned MHD, Sokolov et al. 2022), a different axis entirely.
- Recent applications **use the closure unchanged**: 2024 total-eclipse forecasting
  (arXiv 2503.10974), far-side magnetic structures (arXiv 2509.02911), flux-rope eruptions
  (arXiv 2606.19546).
- The 2024 eclipse work sets the threaded region to **1.0–1.01 R☉** — roughly 7000 km, *thinner* than
  the 1.03 R☉ of the 2021 paper. The closure is not only unrevised, it is being applied under more
  compression.
- AWSoM_R 1.0 is a **community model hosted at NASA CCMC** (runs-on-request), and Koban et al. (2025)
  describe AWSoM and MAS as "two of the most widely used coronal models."
- **The inner-boundary reservoir of §1.4 is unchanged as of 2026.** van der Holst et al. (2026), the
  multi-ion solar wind model (arXiv:2606.19232), still sets `n_p = 2×10¹⁷ m⁻³` "overestimate[d] […]
  by an order of magnitude." Worth noting for the watch: **that paper drops the evaporation language
  entirely** and justifies the same number on purely radiative grounds — "with this density
  overestimation, radiative cooling will be significant. By maintaining a temperature floor of
  50,000 K the density will fall radially outward until radiative cooling is no longer able to drop
  the temperature below 50,000 K." The *number* is not moving; the *stated reason* is drifting. This
  is the thing to re-check every six months, and a shift toward a radiatively-argued boundary would
  be the first sign the gap is being closed upstream.

**A published symptom already exists.** The AWSoM-R stellar extension (2025, ApJ, "High-Resolution
Modeling of Coronae and Winds in Solar-Type Stars. I.") reports **difficulty capturing the steep
density gradient between transition region and corona, with density overestimated at T < 10^6.2 K**.
That is the observable signature of a mis-set {N_i k_B T} handed to the thread, reported
independently, recently, and from inside the AWSoM-R family. **We do not have to argue that the
problem is real — someone has already seen the symptom without diagnosing the cause.** Check in
Tier 0 whether the T_ch sensitivity is of the right magnitude to explain it.

*Adjacent but distinct — do not confuse them in writing:* Szente et al. (2023, ApJS) implement
non-equilibrium ionization in AWSoM for the **charge states of heavy elements, for spectral
synthesis in the corona**. That is not hydrogen ionization in the chromospheric EOS. A referee may
conflate the two; pre-empt it.

---

## 2. The story

> **AWSoM-R's threads currently terminate in an analytical, steady-state, prescribed-width,
> constant-{N_i T} transition region with zero base heat flux and no chromosphere below it — the
> paper's own words are "we do not include a consistent chromosphere model." This work is the object
> that terminates the thread properly: a cheap, equilibrium-ionization, thermodynamically consistent
> field-aligned chromospheric column that replaces the static bottom closure with a dynamic supply of
> mass and enthalpy, and that audits the eight assumptions the present closure rests on. What it
> replaces is a fixed inner-boundary density, overestimated by a factor of ten relative to the same
> group's own chromospheric-top value, whose stated purpose varies across their papers from
> preventing evaporation to enabling it — a numerical placeholder rather than a physical boundary
> condition.**

Every design choice now has an external justification and stops needing self-defence:

| Design choice | Justification |
| --- | --- |
| 1D field-aligned | Sokolov's own low-β / field-aligned transport argument (§1.2) |
| Three conserved rows (ρ, ρu, ℰ), no carrier storage | The BATS-R-US state-vector shape |
| Truncation at 1600 km | Exactly the merge point Sokolov calls uncertain |
| Equilibrium Saha rather than NEQ | The host is named **AWSoM-R**; cost is the binding constraint |
| Double precision; SWMF exact-Riemann flux; passive E₀ offset | SWMF-native by construction |
| Physical conduction only | The TR closure *is* a conduction-versus-radiation competition; conduction is the input side |
| **No Alfvén waves** | The host applies {S_A/B} as a *boundary condition at the top of the chromosphere*, i.e. it already treats our entire domain as a black box delivering one calibrated constant (§2.6a) |
| **One temperature** | Assumption 8 — the host's own choice in this region — and measured τ_eq = 0.09–3.9 ms across 1600–2153 km, below every competing timescale in every cell (§2.7) |
| **Dynamic mass supply at all** | The host substitutes a fixed reservoir 10× over-dense relative to its own chromospheric-top value, with a reported observable cost (§1.4) |

### 2.1 The audience is larger than one model

The closure is **not unique to AWSoM-R**. Sokolov cites Lionello et al. (2001, 2009) and Downs et al.
(2010) — Predictive Science Inc. papers underlying **MAS**. The same analytical-TR family is the
ancestor of the entire **TRAC** lineage. So the assumption structure spans:

```
Lionello 2001/2009, Downs 2010  ──┬──►  AWSoM-R  (Sokolov 2021, CCMC community model)
   analytical TR closure          ├──►  MAS / PSI
                                  └──►  TRAC family (Johnston, Zhou, Iijima)
```

**Writing consequence: frame the finding as an audit of a class of closure, not a critique of one
code.**

- ❌ "AWSoM-R's boundary condition is wrong."
- ✅ "The constant-{N_i T}, fully ionized, optically thin, steady TR closure — widely adopted —
  departs by X% in the partially ionized merge region."

The second reaches AWSoM-R *and* MAS *and* every TRAC user. The first reaches one group.

### 2.2 Extend, do not correct

This matters because Sokolov is the person this work serves.

**Sokolov flagged the opening himself** — "we do not include a consistent chromosphere model." Doing
the work he invited is not finding his mistake. **Quote that sentence in the Introduction**: it is
simultaneously the motivation, the authorization and the courtesy. Pair it with the stellar-paper
symptom (§1.6) and the paper writes itself:

> A reported symptom of this closure is density overestimation at T < 10^6.2 K [stellar paper]. This
> work identifies one source and supplies a correction.

Diagnosis + correction + an independently observed symptom is a complete paper shape, and it makes
Sokolov an ally or co-author rather than a target.

**The genre is mainstream and well cited.** Bradshaw & Klimchuk (2013) offered no new model — it
showed that a resolution everyone was using gave a 2× density error, and explained the mechanism (88
citations). Johnston et al. (2019, 2020) supplied the fix; Howson et al. (2023) then showed the fix
has frequency-dependent side effects. "We tested assumption X inside widely used model class Y" is a
standard, welcome contribution in this field.

### 2.3 The scope claim, correctly narrowed

The strongest objection is not about dimensionality but about **channel**. Type II spicules are
plausibly the dominant coronal mass supply and are driven by ambipolar diffusion, magnetic tension
and multidimensional reconnection (Martínez-Sykora et al. 2017). A column with no magnetic evolution
and equilibrium ionization misses that channel entirely — not approximately, but structurally.

The only honest response is to scope the claim:

> This model describes the **thermally driven mass-loading channel** (conduction-driven, gentle
> evaporation). It does not claim to cover total coronal mass supply.

That narrowing does **not** damage the coupling story, because AWSoM's coronal heating is Alfvén-wave
turbulence dissipation, which produces exactly the downward heat flux that drives the thermal
channel. We are closing AWSoM's own energy-budget loop, not the Sun's mass budget.

### 2.4 Three claims to retract from earlier framing

1. **Do not use Matsuoka (2024) / chromospheric ambipolar diffusion to motivate the release.** That
   is a non-ideal MHD / two-fluid result and belongs to `src/two_fluid/`, not to the single-fluid
   release. Using it as motivation is a bait-and-switch a referee will catch.
2. **Do not imply this improves the Poynting-flux parameter calibration.** That is a wave-energy
   injection problem. The release has no waves, no induction equation, no magnetic evolution. It does
   not touch that problem at all.
3. **Do not say this replaces TRAC.** Say it provides a *calibratable reference solution* for
   TRAC-like methods. Without radiative losses and impulsive heating we cannot yet even engage the
   problem TRAC solves.

The residue after these retractions is one sharp claim rather than four loose ones, and a single
sharp claim is what survives review.

### 2.5 What the contribution actually is — and the "no new physics" question

The framing in §2 invites an obvious worry: *if the point is a model consistent with AWSoM-R that
can be coupled later, and there is no new physics, is there a paper here at all?*

**Answer the blunt half first. "Consistent with AWSoM-R, couplable later" is a software deliverable,
not a paper.** No referee accepts a manuscript whose claim is that something will plug in one day.
Papers of that shape get rejected or absorbed into someone else's code-description appendix. If the
abstract's second sentence reads *"We present a model consistent with [framework] that can be
coupled to it,"* the work is not ready — and no amount of validation fixes that, because the missing
ingredient is a number, not more evidence.

But that is not what is in hand. Sorted into three bins:

**Bin 1 — genuinely not new. Concede each in one sentence; never defend.**
1D field-aligned hydrodynamics (RADYN, HYDRAD, TFLM itself); Saha/LTE equilibrium EOS (MURaM,
Bifrost, Stagger, CO⁵BOLD all default to it — §6.4); MUSCL-Hancock with an exact-Riemann Godunov
flux (textbook); the physics of conduction-driven chromospheric evaporation (Fisher, Canfield &
McClymont 1985). Claiming novelty in this bin is the fastest available way to lose a referee.

**Bin 2 — new quantitative knowledge about an approximation that is in production use. The paper
lives here.**

- the departure of {N_i k_B T} from the fully ionized, optically thin, steady closure in the
  partially ionized merge region (§5.4 Tier 1/2);
- the fraction of conductive input absorbed by the χ_H n_H x reservoir rather than radiated or
  advected — across 6→22 kK, where x climbs from **0.24 to 1.00** (measured, §2.7), this cannot be
  zero by construction;
- the T_ch sensitivity of Eq. (21) — pure quadrature, days (§5.4 Tier 0);
- the frozen (5/3) versus equilibrium (Γ₁ ≈ 1.09) acoustic contrast: a ~25 % sound-speed difference
  inside the partial-ionization layer, plus a Damköhler map saying where each limit applies (§6.3②).

**None of this is new fundamental physics. All of it is new numbers that change what someone else
does. In this subfield that is the unit of contribution.**

**Bin 3 — new method or criterion; small, but genuinely first.**

- the frozen-composition exact-Riemann flux carrying the Saha ionization energy in SWMF's passive
  E₀ offset — no prior chromospheric art found (§7);
- the storage-precision criterion τ/Δt against 1/u — a general, quotable statement nobody appears to
  have written down (§5.3);
- the boundary-condition failure catalogue, above all that the external hyperbolic datum must live
  in the *pressure*: swapping it for temperature produces catastrophic collapse **while reporting
  clean termination and zero Riemann fallbacks** (§5.4). Anyone building a chromosphere–corona
  interface needs that documented, and nobody has.

**The genre check settles the anxiety.** The load-bearing citations of this exact subfield are
approximation-and-recipe papers, not discovery papers: Lionello et al. (2009) is an artificial TR
broadening; Carlsson & Leenaarts (2012) is a fit to expensive calculations; Bradshaw & Klimchuk
(2013) contains **no new model at all** — it says "the resolution everyone is using gives a 2×
density error, and here is the mechanism" (88 citations); Johnston et al. (2019, 2020) is a
conduction modification. §2.2 already notes this; the thing to internalize is that *"we tested
assumption X inside widely used model class Y and it is off by Z %"* is not a lesser paper in this
field — **it is the modal paper.**

**The one-sentence test.** Before drafting anything, the abstract's second sentence must have the
form *"We find that [quantity] departs by [number] under [conditions]."* If it does not, stop
writing and go get the number. That test is the entire difference between Bin 2 and a software
deliverable, and it is why §5.6 puts Tier 0 first: Tier 0 is the cheapest way to find out whether
Bin 2 has a number in it.

**Venues, attached to deliverables rather than to papers.**

| Deliverable | First choice | Why |
| --- | --- | --- |
| Tier 0+1+2 audit **and the corrected recipe** | **ApJS** | Model grids, tables and detailed numerical descriptions are what ApJS is for; a corrected {N_i k_B T}(T_TR, L_TR, T_ch) fit *is* an ApJS object. Sokolov publishes there |
| The audit as a measurement | **ApJ** | If the headline is a number rather than a table |
| Method and numerics (Angle A) | **A&A**, numerical-methods section | Direct precedent: Johnston et al. (2020) and Zhou et al. (2021) TRAC papers are A&A, and the European TR-closure community reads it |
| Coupling readiness **plus the §3 cost budget as the result** | **JSWSC** (open access) | The natural home for "a component for a real-time framework, with its measured cost." The SWMF community publishes there. **This is the one venue that fits the "AWSoM-R-consistent" framing — and only if §3's numbers are presented as the finding** |
| Angle B precision criterion | **RNAAS** or **A&A Research Note** | §5.3 |
| Model description, low risk | **Solar Physics** | Al Shidi et al. (2019) precedent, same department |

**Do not send the numerics to JCP or CPC.** The method is an assembly of known parts and a
computational-physics referee will say so in one line. The novelty is the *application region*, not
the scheme.

### 2.6 Alfvén waves: inherited on transport, scoped out on momentum, and one real exposure

The release has no waves, no induction equation and no magnetic evolution. That single omission
plays out three different ways and only one of them is a problem.

**(a) As an energy-transport channel through the domain — inherited, exactly as §1.2 inherits the
1D argument.** Sokolov et al. (2021) §2.1 applies the wave boundary condition

```
{S_A/B} = 1.1 × 10⁶ W m⁻² T⁻¹      at the top of the chromosphere, to the outgoing wave only
{L_⊥ √B} = 300 km · √T             (Hollweg 1986 correlation-length scaling)
```

and the thread then carries the w± transport equations and feeds (Γ₋w₋ + Γ₊w₊) into its energy
equation (15). **Everything below that boundary — i.e. this entire domain — is, in the host's own
formulation, a black box whose only job is to deliver that one calibrated constant.** Reflection and
dissipation inside our heights are not computed by AWSoM-R either; they are absorbed into the
number. So the correct reply to "why no Alfvén waves?" is inheritance, not self-defence.

*Corollary, and it confirms §2.4 retraction #2:* the column therefore cannot say anything about the
value of that constant. Its calibration is a separate, active problem — Huang et al. (2022) find the
optimal Poynting-flux parameter correlates with open-field area at ρ_s = 0.96 — and we do not touch
it.

**(b) As a momentum source — also inherited, and scoped out.** Sokolov explicitly drops the Alfvén
wave pressure gradient from Eqs. (14)–(15): it scales as √ρ against the thermal pressure's ρ, and
"keeping this term in Equation (15) would be inconsistent." We inherit that omission. But note what
lives in it: this is the ponderomotive / type-II-spicule / FIP channel. Murabito et al. (2024, PRL)
observed chromospheric Alfvén wave reflection directly and established the ponderomotive force
acting on chromospheric ions but not neutrals as the FIP-effect mechanism; Sakaue & Shibata (2020)
show in 1D MHD that once chromospheric Alfvén nonlinearity exceeds a critical value, chromospheric
dynamics and the wind mass-loss rate **decouple from the photospheric energy input entirely.** The
column cannot represent any of this. Scope per §2.3; do not argue.

**(c) The one real exposure — and it lands on §5.1.**

Wave flux crossing the domain top, from the host's own constant: S_A = 1.1 × 10⁶ · B gives
**1.1 kW m⁻² at B = 10 G and 11 kW m⁻² at B = 100 G**, against chromospheric radiative losses of
≈4 kW m⁻² (quiet Sun) to ≈20 kW m⁻² (active region) (Withbroe & Noyes 1977). So the wave flux is
**25–55 % of the total chromospheric radiative loss** — not a negligible term. The mitigating factor
is that most of that loss occurs *below* our domain, in the H⁻ / Ca II lower chromosphere.

And the literature puts the deposition at our heights specifically:

- **McMurdo et al. (2023)** — phase-mixed Alfvén waves in a single-fluid partially ionized plasma at
  ionization degrees μ = 0.518–0.657, corresponding to **1916–2150 km**, provide enough heating to
  balance quiet-Sun chromospheric radiative losses. That is the upper half of this domain.
- **Soler et al. (2018)** — Ohmic dissipation dominates the low and middle chromosphere,
  **ion–neutral friction dominates the high chromosphere**; the dissipated high-frequency fraction
  balances chromospheric radiative losses.
- **Tu & Song (2013)** — self-consistent 1D plasma + neutral + Maxwell solution; Joule dissipation
  alone balances the quiet chromospheric radiative loss, and the flux transmitted to the corona is
  an order of magnitude below the driver.

> **Consequence: a radiative sink with no source is not a simplification, it is a wrong model.** The
> present release is self-consistent precisely *because* it has neither — the thermal structure is
> carried by the initial condition. Introducing Carlsson & Leenaarts losses without a companion
> heating term makes the chromosphere cool away on the 10²–10³ s timescale, which is the timescale
> of every long run we care about. §5.1 now carries the required companion change.

**(d) What the other 1D codes actually do.** Four families; we are in the first, which is the
majority family of the 1D field-aligned literature.

| Family | Representatives | How waves enter | Relation to us |
| --- | --- | --- | --- |
| Pure hydro, no B, prescribed heating | RADYN (no magnetic field at all), HYDRAD, Johnston TRAC series, Bradshaw & Klimchuk (2013) | not at all; heating is an ad hoc Q_h(s) | **this is us**, and it is the mainstream |
| Prescribed wave-transport equations, no resolved oscillation | Cranmer & van Ballegooijen (2005), ZEPHYR, **Sokolov TFLM / AWSoM** | one extra pair of w± equations; cheap | **this is the host**, and the natural path if waves are ever added |
| Resolved 1D/1.5D single-fluid MHD waves | Suzuki & Inutsuka, Sakaue & Shibata (2020), Soler (2017, 2018, 2026) | solved directly; needs ~250 m and sub-second steps | cost-incompatible with a real-time target |
| Resolved two-fluid waves | JOANNA (Kuźma, Pelekhata, Kraśkiewicz), PIP (Snow, Hillier), Tu & Song (2013), **Al Shidi et al. (2019)** | ion+electron and neutral equations plus collisions | this is `src/two_fluid/`, not the release |

Two notes worth keeping. **Al Shidi et al. (2019)** is a UMich CLaSP + UMass Lowell two-fluid
chromosphere paper (Cohen, Song, Tu) and is already in `../supporting-papers/` — this is Igor's
immediate intellectual neighbourhood and he will know the line of work. And **Soler (2026)** supplies
a ready-made defence against a "you need two fluids" referee: for a broadband photospheric driver
over 0.1–300 mHz, single-fluid MHD and multi-fluid results are almost identical — net coronal energy
flux ~5 % higher single-fluid, and the heating rate underestimated by about 2× only in a narrow band
near 500 km, well below our domain.

### 2.7 Two temperatures: why the host wants them, and the height at which they start to matter

The request to separate T_e and T_i has four distinct justifications. Rank them this way; the first
two are the ones that belong in a paper.

**R1 — the host's state vector is multi-temperature, so the interface demands it.** AWSoM has been
two-temperature (electrons + protons) since van der Holst et al. (2010) and **three-temperature**
since van der Holst et al. (2014): isotropic T_e plus parallel and perpendicular ion temperatures,
with turbulence dissipation apportioned among the three by linear wave damping and nonlinear
stochastic heating. The inner-boundary state of §1.4 is literally written T_e = T_i = T_i∥ = 5×10⁴ K.
A column that hands over a single T forces the host to split it arbitrarily — which hands back the
very uncertainty the work exists to remove.

**R2 — the closure is already an electron closure written as if single-temperature.** This is
assumption 8 (§1.3) and it is the publishable form of the argument: κ₀T^{5/2} is Spitzer *electron*
conduction, Λ_R is Λ_R(T_e), N_e N_i is electron × ion density — three electron quantities collapsed
into one symbol, with the host's stated reason being computational efficiency, in exactly the region
we model. Unpicking it in the thread-termination closure is a natural extension of the §5.4 audit
and costs no new physics argument, only bookkeeping.

**R3 — in a partially ionized layer the source and the sink act on different species.** Conduction,
radiation and electron-impact ionization/recombination are electron channels; ion–neutral friction
and ambipolar heating are ion-and-neutral channels. One temperature therefore averages a heating
term and a cooling term that do not act on the same particles. Martínez-Sykora et al. (2019) find
exactly this asymmetry matters — with NEI and ion–neutral effects together, ambipolar heating in
shock wakes is substantially *less* efficient than under LTE and more efficient in the shock fronts;
Song et al. (2023) go further and attribute the *formation* of the transition region to the
different scale heights of the two fluids. This is the deepest physical reason and the one closest
to the local intellectual context (§2.6d).

**R4 — SWMF history, worth knowing but not a lead argument.** Two temperatures entered the framework
because a single-temperature CME develops an unphysical conductive heat precursor ahead of the shock
(Manchester et al. 2012; Jin et al. 2013 compared 1T and 2T directly for the 2011-03-07 event).
Relevant only if the roadmap includes flux-rope eruptions.

**The counter-measurement, now measured rather than estimated.** Electron–ion energy-equilibration
time from the NRL formulary rate ν_ε = 3.2×10⁻⁹ Z² n_i lnΛ / (μ T_e^{3/2}) (n in cm⁻³, T_e in eV),
evaluated per cell on the canonical release baseline — N=500/R4, 661 cells, CFL 0.50, double
precision, physical conduction only, 22 kK reservoir, t = 4000 s. `T`, `x` and `n_e` are read
directly from the run's own `.gamma_diag` EOS diagnostics, so nothing re-derives Saha. NRL branch 1
(lnΛ = 23 − ln(n_e^{1/2} Z T_e^{−3/2})) applies in every cell; lnΛ = 9.27–12.18, nowhere small or
negative, so the formula is in regime throughout. Script `util/tau_ei_equilibration.py`; scalars in
`outputs/model_column/tau_ei_equilibration.json`; figure
`visualization/model_column/tau_ei_equilibration.png`. Cross-checked against a second independent
double-precision run to four significant figures.

| Location | T [K] | x | n_e [cm⁻³] | lnΛ | τ_eq |
| --- | --- | --- | --- | --- | --- |
| Domain base, 1600.5 km | 6631 | 0.273 | 1.58×10¹¹ | 9.27 | **0.0920 ms** ← minimum |
| Mid-domain, 1850.2 km | 6650 | 0.424 | 8.67×10¹⁰ | 9.57 | **0.163 ms** |
| Domain top, 2152.9 km | 21930 | 1.00 | 1.70×10¹⁰ | 12.18 | **3.92 ms** ← maximum |
| *TR at 10⁶ K — extrapolation, not measured* | *10⁶* | *1* | *~10⁹* | *18.09* | ***13.8 s*** |

**Two corrections to the first-revision hand estimate, and the second one matters.**

1. **The base was wrong by 32×.** The estimate assumed n_e ≈ 5×10⁹ cm⁻³ and x ~ 10⁻²; the model's
   own EOS gives n_e = 1.58×10¹¹ cm⁻³ and **x = 0.273** at 1600 km. τ_eq is therefore 33× *shorter*
   than claimed. (The domain top estimate was right.) **See the flag at the end of this subsection —
   x = 0.27 at 6631 K is worth checking against the rest of the document.**
2. **"Four to five orders of margin against any dynamical timescale" was wrong and must not be
   repeated.** It holds only against the slow, science-relevant timescale. Worst-case ratios over all
   661 cells at t = 4000 s:

   | Competing timescale | worst ratio to τ_eq | where |
   | --- | --- | --- |
   | evaporation / mass loading, 10² s | **2.5×10⁴** (4.4 orders) | 2152.9 km |
   | cell conduction time Δs²/χ | **20.7×** | 2152.6 km |
   | cell sound crossing Δs/c_s | **2.85×** | 2152.6 km |
   | solver Δt = 5.617 ms | **1.43×** | 2152.9 km |

   The margin is strongly height dependent — 1120× against sound crossing at the base, order unity at
   the R4-refined top cell, which is also the CFL-limiting cell.

> **The correct statement: τ_eq is shorter than every competing timescale in every one of the 661
> cells — by ≥2.9× against the fastest resolved physical timescale and by 4.4 orders against the
> evaporation timescale the science claims rest on. So T_e = T_i never fails inside 1600–2153 km, but
> the margin is comfortable rather than overwhelming at the top, and it becomes a poor approximation
> above roughly 10^5.5 K.**

**Electron–neutral collisions are not the controlling channel, contrary to the first revision.** With
σ_en = 10⁻¹⁹ m² (standard e–H elastic momentum transfer), τ_ε^{en} = 43 ms at the base — 468× longer
than τ_eq^{ei}, and that is the *smallest* such ratio in the column. Electron–ion collisions dominate
electron energy relaxation everywhere, including the partially ionized part, because x = 0.27–1.00
rather than 10⁻². Neutrals neither strengthen nor weaken the T_e = T_i argument.

**But they flag the real weak link, which is a different pair.** Ion–neutral energy exchange
(σ ~ 10⁻¹⁸ m²) is 0.14 ms at the base and **exceeds the cell sound-crossing time above 2138 km**,
where the neutral fraction is already down to 4.5 %. So if a multi-temperature treatment is ever
built here, the electron–proton split is the *less* urgent one — the neutral temperature decouples
first. That is a design note for Tier 2, not for the release.

**Time dependence is negligible and the late-time value is the one to quote.** Base τ_eq moves
+0.07 % over 4000 s; top τ_eq moves 4.24 → 3.92 ms (−7.5 %) as the top cell relaxes onto the wall.
The worst-case margin is set at the top and *tightens* with time, so a t = 0 measurement would
overstate it by ~8 %.

**Therefore two temperatures belong in Tier 2, with the upward domain extension, not now** — and the
measurement above is the *evidence* for saying so rather than an excuse for not doing it. "We
measured it; τ_eq is below every timescale in the model with ≥2.9× worst-case margin, and we know
the height at which it stops being true" is a stronger position than implementing it prematurely or
leaving it unaddressed.

##### The x = 0.273 flag, resolved 2026-09-06 — and it turned into a validation result

The measurement raised a flag: the column's base ionization is **x = 0.273 at 6631 K**, more than an
order of magnitude above the "x ~ 10⁻²" this document asserted in three places. Three candidate
explanations were checked. **(a) is correct, and the finding is better than "the number is fine."**

**(b) ruled out — the diagnostic means what it says.** `.gamma_diag` column 3 is `x_eq`, and pure-
hydrogen LTE Saha evaluated analytically at the cell's own (ρ, T) — `(2πm_e kT/h²)^{3/2} =
1.305×10²⁷ m⁻³`, `χ/kT = 23.80`, `x²/(1−x)·n_H = 6.00×10¹⁶` — gives **x = 0.274** against the
diagnostic's 0.273.

**(c) ruled out — the initial condition is right.** `model_column_ic` takes its temperature from
Avrett & Loeser (2008) Table 26, Model C7 (`c7_full_temperature_pchip`), and its density from a
Saha-HSE integration. C7's own row at h = 1617 km is `T = 6633 K` — the column's 6631 K — and
`n_e + n_HI = 5.28×10¹⁷ m⁻³` against the column's `n_H,tot = 5.79×10¹⁷ m⁻³`, agreeing to 10 %.

**(a) confirmed, and it is an asset rather than an explanation.** C7's own ionization at that height
is `x = n_e/(n_e + n_HI) = 1.267/(1.267+4.012) = ` **0.240**.

> **Equilibrium LTE Saha at the column's own (ρ, T) reproduces the semi-empirical NLTE ionization of
> Model C7 to 14 % at the domain base (0.273 vs 0.240), and n_e to 25 %.** That is a validation of the
> release closure against a detailed NLTE model at exactly the height where the closure is most
> suspect, and it belongs in Angle A's asset list (§5.2) and in §6.4's precedent argument.

*Caveat to state when quoting it:* C7's `n_e` includes electrons donated by ionized low-FIP metals.
At 1617 km that contribution is ~10⁻⁴ n_H ≈ 5×10¹³ m⁻³ against n_e = 1.27×10¹⁷, i.e. negligible, so
the comparison is clean **here**. It would not be clean near the temperature minimum, where metals
dominate n_e outright. Do not extend the comparison downward without redoing that estimate.

##### Where "x ~ 10⁻²" actually came from, and the two document corrections it forces

The figure is real but was applied at the wrong heights. From the same C7 table:

| h [km] | T [K] | x = n_e/(n_e+n_HI) |
| --- | --- | --- |
| 1003 | 6225 | **0.0070** |
| 1214 | 6576 | **0.027** |
| 1398 | 6610 | 0.071 |
| 1520 | 6623 | 0.152 |
| **1617 — domain base** | **6633** | **0.240** |
| 1820 | 6652 | 0.447 |
| 2147 | 10 980 | 0.838 |
| 2150 | 15 760 | 0.953 |
| 2152 | 20 510 | 0.973 |

**Correction 1 — about our own domain.** `x ~ 10⁻²–10⁻¹` is C7's value at **h ≈ 1000–1400 km**, the
mid chromosphere, *below* the truncation. At 1600 km it is 0.24. Every statement of the form "x
climbs from ~10⁻² to ~1 across 6–22 kK" (§2.5 Bin 2, §5.4 risk mitigation, §9) should read
**0.24 → 1.00**. The ionization-energy-reservoir argument is unaffected in substance — Δx = 0.73
rather than 0.99, a 27 % reduction in a term that was never marginal.

**Correction 2 — about the audit, and this one cuts against us.** §1.3 argues that Sokolov's implicit
full-ionization identity is "simply false" at T_ch ~ (1÷2)×10⁴ K because x ~ 10⁻²–10⁻¹ there.
**C7 says x = 0.75–0.97 over that range.** Full ionization at T_ch is therefore a **3–25 % error, not
an order-of-magnitude one**, and §1.3's wording is too strong. This is now fixed there.

**But it opens a second, independent Tier 0 measurement.** Across T_ch's unconstrained factor-of-two
freedom, x runs 0.75 → 0.97, so **the assumption-6 error itself varies by a factor of ~8 with a
parameter nobody constrains** — a T_ch sensitivity in the *ionization* channel, entirely separate
from the Λ_R channel Tier 0 already tests, and computable from the same C7 table with no simulation.
Tier 0 should now produce both curves.

### 2.8 Why evaporation is the requested run

Evaporation is not being studied as physics in its own right — the conduction-driven mechanism has
been understood since Fisher, Canfield & McClymont (1985). It is the *instrument*. Five reasons; the
first is decisive and external.

**E1 — the host's boundary is a placeholder whose whole subject is evaporation (§1.4).** AWSoM fixes
the base density a factor of ten above its own chromospheric-top value, and evaporation is what every
one of the group's four different justifications for that number is *about* — preventing it, enabling
it, replenishing what it removes, or limiting it to non-excessive amounts. A request for evaporation
runs is a request for the quantity that replaces the placeholder: a dynamic mass and enthalpy flux
that responds to coronal conditions instead of a constant. Nothing else in the model's output does
that.

**E2 — it is the one test Sokolov's "negligible" argument does not cover.** That argument is about
the *location* of the merge point and is static and geometric (§1.4). The time-dependent mass flux
*across* the merge is untouched by it, and no static computation can reach it (§5.5).

**E3 — it breaks assumptions 3 and 4 simultaneously.** Assumption 4 is that the base heat flux
vanishes; evaporation is what happens when it does not. Assumption 3 is steady state; evaporation is
its dynamic violation. One run audits two rows of the §5.4 table.

**E4 — it is the field's standard benchmark, i.e. the language a referee already speaks, and it
supplies the cross-code comparison §5.2 lists as a blocker.** Bradshaw & Klimchuk (2013), Johnston
et al. (2019, 2020) and Zhou et al. (2021) all report "peak coronal density error" from an
evaporation problem. More usefully, **Longcope (2014)** gives a closed-form scaling for the maximum
evaporation velocity,

```
v_e ≈ 0.38 (F / ρ_co,0)^{1/3}
```

and states that it fits simulations across a range of previously published transition-region
treatments. **A heat-flux sweep testing the 1/3 exponent is a cross-code validation that requires
running nobody else's code** — this is the cheapest available route to §5.2's "add one comparison
against a published result from another code."

**E5 — it is the most sensitive probe of the interface being designed.** The boundary studies
already show the evaporation signal moving by ~19–20× across defensible top-boundary choices while
the wall-temperature exponent barely moves (2.80 vs 2.91) — see §5.4. A diagnostic with 20× dynamic
range against interface design is exactly the diagnostic to design the interface with.

**The attribution gap — reframed by reading Lionello, and no longer the load-bearing question.**
AWSoM's low-height density overestimate has two candidate causes, and the 2026-09-06 source check
established that **the authors themselves name both, in the same paragraphs, without separating
them**: Sachdeva et al. (2019) §4.2 attributes it first to "artificial broadening of the transition
region" and then, three paragraphs later, to the inner-boundary overestimate; Shi et al. (2022) §3.4
lists both in one sentence. Lionello et al. (2009) does not separate them either — its Figures 1 and
3 exhibit both, acting in the same direction (§1.4 F3). **This is not resolvable by reading, and the
attempt should stop.**

**But §1.4's F2 makes the separation much less important than it looked.** The stronger and simpler
claim does not require attributing blame at all: *both* approximations are validated by their
founding paper only above ~500,000 K, and the excess is reported below that. So the sentence to write
is not "the reservoir causes the density error" but **"neither approximation is claimed to be valid
where the error appears, and the founding paper says so."** That survives whichever mechanism
dominates.

Tier 0 still measures the T_ch sensitivity, and a large signal still strengthens the case — but a
small one no longer voids the programme, because the scope-of-validity argument stands on its own.
**This is a material reduction in the plan's single largest risk** (§9, "Tier 0 could return a small
signal"). Ask the SWMF group anyway (action 8); their answer is cheap and would sharpen the story.

---

## 3. Cost — the numbers, measured

### 3.1 This column, measured 2026-09-06

Release configuration (N=500/R4, 661 cells, CFL 0.50, double precision, `build_omp`, this
workstation). Differential timing at T_end = 5 s and 40 s to remove startup cost:

| T_end | wall | steps |
| --- | --- | --- |
| 5 s | 1.067 s | 889 |
| 40 s | 8.017 s | 7110 |
| **difference** | **6.950 s / 35 physical s** | **6221** |

```
marginal cost = 0.199 core-s per physical second   (≈ 5× faster than real time, one core)
              = 1.117 ms per step
              = 1.69 μs per cell-step
startup       = 0.074 s (negligible)
dt = 5.626 ms; 178 steps per physical second
```

**Twelve threads give only 1.25×** (3.990 → 3.201 s over 20 physical s). At 661 cells there is not
enough work to parallelize. This is *good news* for the coupling architecture: threads should not be
spent inside one column but across many columns, one core per thread — which is exactly what TFLM
needs. Note the OpenMP recap's 36.3 s / 2.170× figures are for a 2638-cell mesh and predate the
double-precision cutover; do not quote them for the release configuration.

### 3.2 The AWSoM-R budget

TFLM traces one thread per grid point of the GCM lower boundary. AWSoM's background coronal angular
resolution is ≈1.4° (Shi et al. 2022 describe 0.35° AR refinement as "four times higher than the
background corona"):

```
N_threads ≈ 4π / (1.4° in rad)² ≈ 12.566 / 5.97e-4 ≈ 2.1 × 10⁴
```

**This is an estimate, not a published number — verify the actual TFLM thread count with the SWMF
group before quoting it anywhere.**

```
2.1×10⁴ threads × 0.199 core-s/s  =  4180 core-s per physical second
on 200 cores                      =  ~21 s wall per physical second  →  ~21× slower than real time
real time would need              ≈  4200 cores, or a ~21× cost reduction
```

**The gap is ≈1.3 orders of magnitude, not the six or seven a first guess suggests.** That is a
tractable engineering target, not a fantasy. Also note this is conservative: TFLM's "faster than
real time on ~200 cores" includes the full 3D MHD, not just the threads.

### 3.3 The lever

Cost is set entirely by the **acoustic CFL**: dt = 5.6 ms from c_s ≈ 10 km s⁻¹, while |u| ≈ 5–10 m
s⁻¹ — **Mach ≈ 10⁻³**. Up to three orders of magnitude are being spent on sound waves that carry no
signal of interest.

The tool already surfaced in the well-balanced literature survey:

> **Thomann et al. (2019), "An all speed second order well-balanced IMEX relaxation scheme for the
> Euler equations with gravity"** (`10.1016/j.jcp.2020.109723`) — flux split into a linear implicit
> and a nonlinear explicit part giving a **scale-independent time step**, asymptotic preserving,
> well-balanced, positivity preserving.

Recovering even a fraction of the 21× closes the budget, and doing so is itself a publishable
numerical result.

**Interface subcycling is not a workaround but the natural design.** Secondary sources report that
AWSoM-R advances with a timestep of **≈1 s** (versus milliseconds for the original AWSoM — "a gain of
a factor of several hundred"). Our dt is 5.6 ms, a factor ≈180 smaller. The host neither needs nor
wants chromospheric updates at our cadence, so an interface that exchanges time-averaged {N_i k_B T}
and heat flux once per host step is the correct coupling, and the cost comparison should be made
against *that* exchange rate rather than against lockstep advance. This softens the 21× from a
numerical problem into an interface-design question. **Verify the ≈1 s figure with the SWMF group —
it comes from secondary sources, not the 2021 paper.**

Other levers: only a subset of threads need dynamics at any time; 1.69 μs per cell-step is itself
slow for a three-equation update (Saha Newton inversion plus implicit conduction), so there is
constant-factor headroom.

---

## 4. Context: why the existing 1D codes do not close this gap

Organized by *why each one cannot terminate a TFLM thread*, not by literature family. **§2.6d gives
the complementary cut** — the same codes sorted by *how Alfvén waves enter*, which is the axis a
referee will use when asking why this column has none.

**RADYN — the reference solution, at its cost ceiling.** 1D field-aligned, implicit, adaptive mesh,
full NLTE radiative transfer plus non-equilibrium populations for H (6 levels), Ca II and He. The
standard configuration is **191 grid points**, and Rubio da Costa et al. (2015) state plainly:
"191 grid points (**a denser grid would be computationally too expensive**)." Kowalski et al. (2017):
"The high beam flux simulations are computationally expensive in 1D." **No wall-clock cost is
published anywhere I could find** — which is itself a finding, and an opportunity: a paper that
reports a cost/accuracy trade-off honestly is doing something the field currently does not.

**Bifrost — no separate 1D fork, and it does not matter.** Bifrost (Gudiksen et al. 2011) supports
1D/2D/3D Cartesian configurations, but nobody publishes 1D Bifrost science. The Oslo 1D workhorse is
RADYN: implicit, adaptive, full NLTE/NEQ, **no magnetic field**.

**The decisive precedent is how these two already relate.** Bifrost does not pay RADYN's cost; it
uses the **Carlsson & Leenaarts (2012)** cheap radiative-loss recipes, which were *derived by fitting
detailed RADYN calculations* precisely so that the large model would not have to. Bifrost's own code
paper validates its approximate cooling against a detailed RADYN run.

> **"Distil the expensive 1D physics into a cheap recipe for the big model" is already the community's
> accepted methodology.** This work does the same thing one level up: for ionization energetics
> instead of radiative losses, and for AWSoM-R instead of Bifrost. That is the methodological anchor
> for the whole programme — and it is also the thing we still have to add (§5.1).

**TRAC / LTRAC / jump conditions — solving the opposite problem.** Johnston et al. (2019, 2020)
broaden the TR so coarse grids get the coronal density response right, explicitly sacrificing TR
accuracy; Bradshaw & Klimchuk (2013) established the underlying resolution failure. Howson et al.
(2023) then showed the thermodynamic modification changes the rate of energy injection into the
corona in a **frequency-dependent** way. These methods aim to *avoid* modelling the chromosphere; we
aim to terminate a thread in one. Complementary, not competing — and Howson is the community's own
statement that the avoidance has consequences.

**Equilibrium-ionization EOS in production codes.** Rempel (2016) uses "an equilibrium ionization
equation of state" in the MURaM coronal extension; Bifrost, MURaM, Stagger and CO⁵BOLD all run LTE
EOS tables by default with NEQ hydrogen as an optional module. We are not an outlier (§6.4).

### Citations missing from `reference.bib`

`Sokolov2021` (ApJ 908, 172) — **the anchor**, and currently absent; Lionello et al. (2001, 2009);
Downs et al. (2010); Rempel (2016); Howson et al. (2023); Przybylski et al. (2022); Kerr et al.
(2019); Brown (1973); Ballester et al. (2021, `10.1051/0004-6361/202141851`, Γ₁ in partially ionized
plasma — near-duplicate of our §2.3 framing); Varonov et al. (2023, `10.3847/1538-4357/ad1a12`,
1.1 < γ ≤ 5/3 from Saha in LTE).

**Added by the second revision** (all cited in §§1.4, 2.5–2.8):

- **van der Holst et al. (2014)**, ApJ 782, 81, `10.1088/0004-637X/782/2/81`, arXiv:1311.4093 —
  **the primary source for §1.4**: the inner-boundary values, the "order of magnitude" overestimate
  and the "prevents chromospheric evaporation" rationale, *and* the AWSoM three-temperature
  formulation behind assumption 8 (§2.7 R1). **The single most important addition in this list.**
- **Sokolov et al. (2013)**, ApJ 764, 23, arXiv:1208.3141 §3.2.1 Eq. (36) — `T_ch = 2×10⁴ K,
  N_ch ≈ 2×10¹⁶ m⁻³`, the chromospheric-top reference point that calibrates §1.4's factor of ten.
  Note: it does **not** itself mention evaporation, despite being cited for it.
- **Sachdeva et al. (2019)**, ApJ 887, 83, arXiv:1910.08110 §4.2 — the DEMT density comparison and
  the dual attribution (broadening *and* reservoir) of §2.8's open question.
- **Sachdeva et al. (2021)**, ApJ 923, 176, `10.3847/1538-4357/ac307c` — the 5×10¹⁸ m⁻³ variant
  (**25×**, not one order — see §1.4).
- **Shi et al. (2022)**, ApJ 928, 34, `10.3847/1538-4357/ac52ab` — restatement of the boundary and
  the "excessive evaporation" hybrid framing; **Shi et al. (2024)**, ApJ 961, 60 — the "above
  1.01–1.03 R☉" validity floor.
- **Jin et al. (2017a)**, ApJ 834, 172, arXiv:1611.08897 (and 2017b, ApJ 834, 173, arXiv:1605.05360);
  **Sachdeva et al. (2023)**, arXiv:2212.05138; **Panesar et al. (2026)**, arXiv:2606.17361 — the
  *opposite-valence* and *replenishment* framings of the same reservoir (§1.4 table).
- **van der Holst et al. (2026)**, multi-ion model, arXiv:2606.19232 — the reservoir unchanged in
  2026, with the evaporation rationale replaced by a radiative one (§1.6).
- **Lionello, Linker & Mikić (2009)**, ApJ 690, 902, `10.1088/0004-637X/690/1/902` — **obtained and
  read 2026-09-06**, `../supporting-papers/LionelloLinkerMikic2009.pdf`. The founding paper of the
  whole closure family and the sole authority for "the overestimate does not change the global
  solution." **Cite it for its own scope limit, not only for the claim** (§1.4 F1–F5): the factor-2
  Figure 1 test, the κ(T) = const below T_c = 500,000 K broadening with its ~400× TR widening, and
  the sentence "more accurate calculations at lower temperatures can be achieved by lowering the
  temperature at which the modification occurs, with the consequence of a requirement for higher
  resolution." Also **Mok et al. (2005)**, the onward citation Lionello gives for the
  pressure-overestimate claim — one more hop, not yet checked, probably not worth a trip.
- **Manchester et al. (2012)**; **Jin et al. (2013)**, `10.1088/0004-637X/773/1/50` — the 1T vs 2T
  CME heat-precursor argument (§2.7 R4).
- **Longcope (2014)**, `10.1088/0004-637X/795/1/10` — v_e ≈ 0.38 (F/ρ_co,0)^{1/3}; the cross-code
  comparison of §2.8 E4 and §5.2.
- **Fisher, Canfield & McClymont (1985)** — gentle vs explosive evaporation; the Bin 1 concession
  that the evaporation *mechanism* is not new (§2.5). Already in `../supporting-papers/`.
- **McMurdo et al. (2023)**, `10.3847/1538-4357/ad0364`; **Soler et al. (2018)**,
  `10.3847/1538-4357/aaf64c`; **Tu & Song (2013)**, `10.1088/0004-637X/777/1/53` — chromospheric
  Alfvén-wave heating at *our heights*; the evidence behind the §5.1 sink-needs-a-source change.
- **Soler (2026)**, `10.1051/0004-6361/202659417` — single-fluid ≈ multi-fluid for wave transport;
  the pre-emption of a "you need two fluids" referee (§2.6d).
- **Murabito et al. (2024)**, `10.1103/PhysRevLett.132.215201`; **Sakaue & Shibata (2020)**,
  `10.3847/1538-4357/ababa0` — the ponderomotive/FIP and nonlinear-decoupling channels the column
  structurally cannot represent (§2.6b); cite where §2.3 scopes the claim.
- **Huang et al. (2022)**, `10.3847/2041-8213/acc5ef` — Poynting-flux parameter calibration; cite
  when stating §2.4 retraction #2, to show the problem is active and separate.
- **Song et al. (2023)**, `10.3847/2041-8213/accc27`; **Martínez-Sykora et al. (2019)**,
  `10.3847/1538-4357/ab643f` — species-asymmetric energetics in a partially ionized layer (§2.7 R3).
- **Al Shidi et al. (2019)**, `10.1007/s11207-019-1513-8` — local two-fluid chromosphere precedent
  from the same department (§2.6d). Already in `../supporting-papers/`.
- **Withbroe & Noyes (1977)** — the chromospheric radiative-loss budget the §2.6c comparison rests
  on.

---

## 5. Publication plan

### 5.1 Radiative losses are now the first necessity, not the first risk

Previously listed as the largest referee risk. With AWSoM-R as the axis it is more than that:
**Λ_R is the core of Sokolov's Eqs. (21)–(22).** Without a radiative sink the column cannot engage
the actual TFLM boundary condition at all. The direction is also fixed: at 1–2×10⁴ K an optically
thin CHIANTI table is the wrong instrument, and **Carlsson & Leenaarts (2012)** — already in
`reference.bib` and in `../supporting-papers/CL2012 Figures` — is the right one.

**This single item is simultaneously the missing physics, the interface to the host, and the
methodological precedent of §4.** It should be scheduled first.

#### It cannot be implemented alone — a sink needs a source

**Do not add the radiative sink without a companion heating term in the same change.** The present
release is energetically self-consistent precisely because it has *neither* source nor sink: the
thermal structure is carried by the initial condition, and the only energy input is conduction from
the outer reservoir. Adding Λ_R alone gives the chromosphere a net loss at every step and it cools
away on the 10²–10³ s timescale — which is the timescale of every long run the programme depends on.
The missing source is real chromospheric heating, and §2.6c shows the literature puts a plausible
one (phase-mixed and ion–neutral-damped Alfvén waves) at 1900–2150 km, i.e. inside this domain, at
25–55 % of the total chromospheric radiative loss.

Three options; **take the first.**

1. **Background heating locked to the initial state, Q_h(s) = Λ_R(ρ₀, T₀).** Smallest change, no new
   free parameter, and it is what the 1D field-aligned literature does — HYDRAD, the Johnston TRAC
   series and Bradshaw & Klimchuk all carry a prescribed Q_h. Defensible in one sentence: *we do not
   claim to explain chromospheric heating; we require the initial state to be in energy balance so
   that the measured departure is caused by the dynamics and not by an unbalanced initial condition.*
   This also makes the measurement cleaner, because the audited quantity becomes a *difference* from
   a balanced reference.
2. Deposit a fraction of the host's {S_A/B} with an assumed damping length. More "physical," but it
   introduces two free parameters and walks straight into §2.4 retraction #2.
3. Apply Λ_R only above T_ch. Cheapest, but it abandons the audit of assumption 5 — whether the
   optically thin CHIANTI table is valid at 1–2×10⁴ K — which is one of the four Tier 1 targets.
   **Not acceptable.**

Option 1 changes the §5.6 sequence not at all; it changes the *content* of action 7 in §8.

### 5.2 Angle A — Method paper

**Claim.** A thermodynamically consistent equilibrium-ionization closure for a field-aligned
chromospheric column: frozen-composition exact-Riemann flux carrying the Saha ionization energy in a
passive specific-energy offset, one EOS closure applied consistently across initialization,
reconstruction, flux, conduction and boundaries.

**Assets in hand.** **Equilibrium Saha reproduces Model C7's semi-empirical NLTE ionization at the
domain base to 14 % (x = 0.273 vs 0.240) and n_e to 25 %** (§2.7) — a validation of the closure
against a detailed NLTE model at the height where it is most suspect, and the strongest single answer
to "why is equilibrium ionization good enough here"; Saha to 2×10⁻¹³; caloric round trips; Γ₁ table
audited to 3.1×10⁻⁴; isentropic
perturbation sound speed to 5×10⁻⁴; Toro's five reference Riemann problems; the shock-tube test that
the E₀ offset is inert; 2.4×10⁷ faces with zero Rusanov fallbacks; mesh-generator tests; the
second-order hydrostatic ladder over `ISO_NS` = 48…2000; the measured cost of §3.1.

**Targets.** A&A (numerical methods) first — see the deliverable/venue table in §2.5, which supersedes
this line. Note that Angle A on its own is the weakest of the three against §2.5's one-sentence test:
"we implemented a consistent closure" is a Bin 1/Bin 3 claim, which is why the "must add" item below
is not optional.

**Must add.** "Frozen and equilibrium fluxes agree to a few parts in 10⁴" is a **null result**. Add a
case where the closure choice matters: an acoustic or weak-shock test inside the partial-ionization
layer, where Γ₁ = 1.09 versus 5/3 changes the sound speed by ~25%. Quantify wave travel time, shock
strength, deposition height. Cheap with the existing code; converts "we implemented X" into "the
closure changes Y by Z%." **Highest return on effort in this document.**

**Other blockers.** Finish the matched N=1000 long-duration run; add one comparison against a
published result from another code. **The cheapest route to that comparison is now identified:**
Longcope (2014)'s v_e ≈ 0.38 (F/ρ_co,0)^{1/3} evaporation-velocity scaling, fitted across a range of
published transition-region treatments — a heat-flux sweep testing the 1/3 exponent needs nobody
else's code (§2.8 E4).

### 5.3 Angle B — Research note on storage precision

**Claim.** A resolved chromosphere manufactures a timescale separation that breaks single-precision
conservative updates. Small Δt (5.6 ms, acoustic CFL on a fine mesh) with long mass-loading time
(τ ≈ 1.9×10⁵ s) gives τ/Δt ≈ 3.5×10⁷ against a reciprocal unit round-off of 1.7×10⁷ for `float` and
9.0×10¹⁵ for `double`. In single precision the per-step continuity increment fell below half a ULP
of ρ in **619 of 660 cells** — permanent stagnation, since there is no compensated summation.

**Supporting evidence, already measured.** The float32/float64 ladder contrast (clean order ≈2 across
41× refinement in double; order 0.79 by N=1000 and *reversal* at N=2000 in float); column mass budget
−214.9 → 0.994; the SWMF audit (all 23 `share/build/Makefile.<OS>.<compiler>` templates set
`PRECISION = ${DOUBLEPREC}`; BATS-R-US accumulates U^{n+1} = U^n + dt·R in float64 with no compensated
summation); and the consequence that **any evaporation number from a float32 run is high by order
30%**.

**Why separate.** A general, quotable criterion that nobody appears to have written down, directly
relevant to a double-precision host framework, and it gives Angle A something of our own to cite.
RNAAS or an A&A Research Note; 3–4 pages, a few weeks.

### 5.4 Angle C — Audit the closure. **Do not couple.**

**Full SWMF coupling is not required and should not be attempted.** A component wrapper, MPI
coupling, build integration and a regression suite is software engineering, not science, and it is
not what convinces a referee. What convinces a referee is evidence that the existing closure is
quantifiably wrong. **All of Sokolov's assumptions are testable offline, with the column as a
reference solution and zero SWMF integration.**

Reframe the deliverable: **audit, then supply a recipe** — the Carlsson & Leenaarts (2012) pattern of
§4. Sokolov's side would replace one function, not adopt our code.

#### Tier 0 — zero simulation, days. Do this first.

**The T_ch sensitivity of Eq. (21) is pure quadrature.** T_ch is a free parameter over
(1÷2)×10⁴ K — a factor of two — and it is the lower limit of both integrals. Take a CHIANTI table,
sweep T_ch, and plot the resulting {N_i k_B T} and heat flux.

If the closure's output is strongly sensitive to an unconstrained parameter, **that alone is a
finding**: the number handed to every thread depends on something nobody constrains, and our column
is the instrument that can constrain it. Cross-check the magnitude against the stellar-paper density
overestimate (§1.6).

**Tier 0 also decides whether the rest is worth doing.** Run it before committing.

**Tier 0 now has a second quadrature, added 2026-09-06 (§2.7).** Sweep T_ch over the same
(1÷2)×10⁴ K range against **Model C7's ionization fraction**, which runs x = 0.75 → 0.97 there, and
plot the resulting error in Eq. (20)'s implicit full-ionization identity. This is a T_ch sensitivity
in the **ionization** channel, independent of the Λ_R channel, and it needs only the C7 table already
compiled into `scenarios/model_c7.cpp` — no simulation and no CHIANTI. Two curves from one afternoon:
if *either* is steep, the closure's output depends on an unconstrained parameter, which is the Tier 0
finding either way.

**Tier 0 gains a third item, added 2026-09-06 — the cheapest real contact with the closure anywhere
in this plan.** Evaluate Eq. (21) for a realistic `(T_TR, L_TR = 1 Mm, T_ch)` and compare its
`{N_i k_B T}` against **the pressure the column actually has at its top face**: `p_top ≈ 0.0103 Pa`
at 22 kK, fully ionized, so `{N_i k_B T} = p_top/2 ≈ 0.0051 Pa`.

This is not circular, because the column's top pressure is **not** an output of Eq. (21) — it comes
from the Saha–HSE integration of the Model C7 atmosphere (§2.7). So the comparison reads:

> **Does the TFLM analytical closure reproduce the pressure that a detailed semi-empirical
> atmosphere has at the merge height?**

**Zero simulation, no domain extension, one afternoon** — and if the two differ by a factor of a few
it is a Bin 2 number that can be written down immediately. It also puts a real value on `p_ref` for
the Tier 1 pressure sweep below, replacing "whatever the IC happened to give."

**Tier 0 carries a second decision criterion, added per §2.8.** AWSoM's reported low-height density
overestimate has two entangled causes — the deliberately over-dense inner-boundary reservoir (§1.4)
and Lionello's artificial TR broadening, since the inner boundary sits at the base of the broadened
TR. Compare the magnitude of the T_ch sensitivity against the reported overestimate: **if it is of
the right size, the merge region is the leading cause and Tiers 1–2 are aimed correctly; if it is
much too small, the broadening dominates, this column addresses a second-order effect, and the plan
must be re-ranked before any significant investment.** Answering this costs nothing extra — it is
the same quadrature, read against a different number.

#### Tier 1 — existing 553 km domain + CL2012, weeks. The body of the paper.

Once radiative losses are in, **four of the assumptions are testable without extending the domain**,
because they all fail near T_ch, and T_ch sits at the top of our existing column:

| Assumption | Measurement | Domain extension? |
| --- | --- | --- |
| 6 — implicit full ionization | Column supplies the true x; compare N_e N_i against N_i². **Now known to be a 3–25 % error at T_ch, not an order of magnitude** (§1.3, §2.7) — quote the range, not a single number | no |
| 7 — no ionization-energy reservoir | Fraction of conductive input stored in χ_H n_H x rather than radiated or advected | no |
| 4 — base heat flux → 0 | Measure the base conductive flux directly | no |
| 5 — optically thin CHIANTI Λ_R at 10⁴ K | Compare against Carlsson & Leenaarts (2012) along the run's own (ρ,T) trajectory | no |
| 2 — {N_i T} = const | Measure its variation over 6–22 kK and its dynamic departure | partial |
| 1 — L_TR ≈ 1 Mm | — | **yes** |
| 3 — steady state | Departure from conduction–radiation balance over time | partial |
| 8 — T_e = T_i | **Report τ_eq(ρ, T) from release output** (§2.7): quantify that the single-temperature closure is safe over 6–22 kK and name the height where it fails | no (the *measurement* is Tier 1; the *implementation* is Tier 2) |

Assumption 8 is the cheapest row in the table — τ_eq is a post-processing formula applied to
existing output, no new code, no new run — and it converts an unaddressed weakness into a bounded,
quantified statement. Do it alongside the Damköhler map of §6.3②, which is the same kind of
post-processing on the same output.

**The frozen-x bracketing run (§6.3③) is mandatory at this tier**, not optional: the headline number
depends on the ionization treatment, so it must carry an interval.

##### The outer pressure sweep is also mandatory at this tier — for two independent reasons

Every number in the table above is measured *given* the thread's handoff, i.e. given `kOuterPRef`.
**That value has never been varied.** It is captured from the initial condition on the first
`update_bc` call (`scenarios/model_column.cpp`) and then held for the whole run. So the audit
currently produces "the departure is X%" with no answer to "under what conditions."

**Reason 1 — the audit numbers need a sensitivity.** "Departs by 40%" is not a Bin 2 result until it
carries the conditions it was measured under. A 3–5 point sweep supplies that at ~11 min per run.

**Reason 2 — the recipe cannot be inverted without it** (Tier 2, above). The pressure is the axis
along which the four-quantity relation is entered from our side; a single value gives a
one-dimensional slice that cannot be inverted into Sokolov's parameterization.

**Implementation is small:** add an environment override (e.g. `ISO_OUTER_PREF`) so the datum can be
set explicitly, falling back to the present IC-derived behaviour when unset. Tier 0's third item
(above) supplies a physically motivated value to centre the sweep on.

**Do not respond to this by loosening the boundary instead.** The fixed external pressure datum is
correct and is defended in the boundary-studies block below; a `dp/ds = −ρg` closure leaves the
pressure *level* unanchored, `p_top` drifts upward without bound (`~ln t`, no steady state), and the
column drains from the base to shed weight. The measured identity is clean — `Δ(Mg)` matches
`Δ(p_base − p_top)` to four digits and `−(1/g)·dp_top/dt` matches `dM/dt` to 0.6% — but every
absolute velocity comes out ~20× too small and the reference quantity is itself moving. **For a paper
whose deliverable is a pressure, an unanchored pressure is disqualifying.**

#### Tier 2 — extend the domain upward, months. Converts audit into correction.

Extend to span T_TR, then deliver a corrected fit or table for
**{N_i k_B T}(T_TR, L_TR, T_ch)** including partial ionization and the ionization-energy sink. This
is the CL2012 pattern exactly: distil the expensive computation into a cheap recipe. **Zero coupling
engineering — the deliverable is a function.**

##### How far up — specify a temperature, not a height

**The requirement is a temperature.** Equation (21) integrates to `T_TR`, not to a height:
`{N_i k_B T} = (1/L_TR) ∫_{T_ch}^{T_TR} κ₀T^{5/2}/U_heat(T) dT`. Write the spec as "the domain top
reaches T_TR" and let the height be an *output* of the run. Two reasons this matters: a real thread —
especially in an active region — has a steeper, thinner TR than the quiet-Sun stratification the
numbers below are taken from, and `T_TR` is precisely the control parameter the deliverable is a
function of.

**The target is T_TR ≈ 4×10⁵ K**, fixed by self-consistency with assumption 1 rather than chosen:
`T_TR` is a *given* in Sokolov's formulation ("for the known width L_TR of the TR and for a given
temperature T_TR on top of the TR"), but it is tied to `L_TR ≈ 1 Mm` — the closure covers exactly the
1 Mm above `T_ch`. On Model C7's stratification, `T_ch ≈ 1.6×10⁴ K` sits at 2150 km, and 1 Mm above
that is 3150 km, where T ≈ 3.95×10⁵ K.

Heights from the C7 table already compiled into `scenarios/model_c7.cpp` (present domain top:
2153 km, 22 kK):

| Target T | C7 height | Extension | What it buys |
| --- | --- | --- | --- |
| 5×10⁴ K | 2161 km | **+8 km** | **AWSoM's own inner-boundary temperature.** Essentially free |
| 1×10⁵ K | 2202 km | +49 km | Through the hydrogen Λ_R bump; where assumption 5's CL2012→CHIANTI crossover lives |
| 2×10⁵ K | 2378 km | +225 km | |
| 3×10⁵ K | 2688 km | +535 km | Domain length doubles |
| **4×10⁵ K — the self-consistent target** | **~3150 km** | **+1.0 Mm** | **Completes the Eq. (21) integral → the deliverable** |
| 5×10⁵ K | 4296 km | +2.1 Mm | |
| 1×10⁶ K | ~17 Mm | +17 Mm | **Do not.** That is corona, not TR, and contributes nothing to the closure |

**Note the first row.** Reaching 50,000 K — the temperature at which AWSoM applies its inner boundary
condition (§1.4) — costs **8 km**, i.e. moving the outer conductive reservoir from 22 kK to 50 kK.
That is plausibly the cheapest publishable point of contact with the host anywhere in this plan, and
it should be checked before Step A below.

##### Do it in two steps

- **Step A — +50 km, top at 10⁵ K. A de-risking step, not a science step.** It (i) crosses the
  hydrogen radiative peak and measures the CL2012 → optically-thin CHIANTI crossover temperature,
  **which is itself the assumption-5 audit**; (ii) subjects the implicit conduction solve to a ~44×
  jump in κ before the full 1414× of Step B; and (iii) forces the mesh question, since 22→100 kK is
  only 49 km and is the steepest stretch in the whole domain — it needs ≲1 km cells on the Bradshaw &
  Klimchuk resolution criterion.
- **Step B — +1.0 Mm, top at ~4×10⁵ K.** The full integral and the deliverable.

##### Cost is not the constraint; three pieces of physics are

**Compute barely moves.** The 5.6 ms step is set by the fine cells low down, and the new cells can be
stretched: at 4×10⁵ K, `c_s ≈ 105 km s⁻¹`, so 5 km cells give `dt = 24 ms` — *larger* than the
present step. Cost therefore scales with cell count, roughly 661 → 1200–1500, i.e. **~2×**: a 4000 s
run goes from ~11 min to ~20-odd min. Negligible.

The real cost is that three things become mandatory together:

1. **Two radiative-loss tables and a blend** — CL2012 below ~30 kK, optically thin CHIANTI above. Not
   a burden: **where the crossover falls is the assumption-5 result.**
2. **Separate electron and ion temperatures.** §2.7 measures τ_eq = 0.09–3.9 ms over the present
   domain, below every competing timescale in every cell, so the split cannot change the Tier 1
   answer. At 4×10⁵ K and n_e ~ 10⁹ cm⁻³ the estimate is τ_eq ≈ 3.7 s, leaving only ~2 orders against
   conductive and evaporative timescales rather than the 4.4 orders available now — and the recipe's
   output then becomes a multi-temperature object the three-temperature host can consume directly
   (§2.7 R1). **The extension and the temperature split are one change; do not schedule them apart.**
3. **Conduction stiffness.** `κ ∝ T^{5/2}` rises by 1414× from 22 kK to 4×10⁵ K. Backward Euler
   handles it, but the conditioning degrades. **This is the most likely failure point of Step B,
   which is exactly what Step A exists to hit first.**

**The boundary architecture does not change.** The existing outer conductive reservoir is already the
right object — extension moves it up and changes its temperature to `T_TR`.

##### The column cannot produce the recipe "directly" — it enters the same relation from the other side

An earlier draft of this section said the run series "produces `{N_i k_B T}(T_TR, L_TR, T_ch)`
directly." **That was wrong, and the error is worth stating explicitly because it looks like a
circularity and is not one.**

The apparent problem: the release imposes a pressure datum `kOuterPRef` at the top face, and
`{N_i k_B T}` *is* a pressure at the top of the column. Imposing the answer and reporting it back
computes nothing.

The resolution is that **Eq. (21) is one relation among four quantities** — `(p, T_TR, L_TR, T_ch)` —
so three must be specified and the fourth follows. *Which* three is a modelling choice, and Sokolov
says as much (§1.3 assumption 1: prescribing the width rather than the top temperature is "an
important distinction from previous works"):

```
Sokolov:   specify (T_TR, L_TR, T_ch)  ->  solve for  p
this work: specify (T_TR, p_top)       ->  read off   L_TR    [for each chosen T_ch]
```

The column is a PDE solver, not a function. `T_TR` and `p_top` are **boundary data**; `T_ch` and
`L_TR` are **read out of the converged profile** (`T_ch` is a post-processing choice of which
isotherm to call the merge point; `L_TR` is then the distance from it to `T_TR`). Same four-tuple,
same manifold, different parameterization.

**Two consequences, both operational:**

1. **The pressure dimension must be swept.** Sweep `(T_TR, p_top)` on a grid, read `(L_TR, T_ch)`
   out of each converged profile, then **invert the resulting point cloud** into Sokolov's
   parameterization. Sweeping `T_TR` alone gives a one-dimensional slice that cannot be inverted.
   This is ordinary numerical calibration of an implicit relation, not a workaround.
2. **`p_ref` must become settable.** Its value is currently a by-product of the initial condition,
   captured on the first `update_bc` call. See the Tier 1 action below — the same change is required
   for an independent reason.

**At the *present* domain the column does not compute `{N_i k_B T}` at all, and must not claim to.**
The top face sits at 22 kK, i.e. essentially at `T_ch` (2146–2152 km against a 2153 km domain top),
while Eq. (21) integrates from there up to `T_TR ≈ 4×10⁵ K` — a full Mm above. **The current column
lies entirely below the integration range: it *consumes* the thread's handoff and audits what happens
underneath.** That is Tier 1 and it is a complete deliverable on its own; producing the recipe is
Tier 2 and requires the extension. Do not let the two blur — §5.6's sequence depends on the
distinction.

##### A fourth argument for the recipe, deferred

Once the extension works, the natural next parameter is the footpoint field strength **B**, because
AWSoM applies the correction *per thread* and each thread has its own B: the area expansion
`A(s) ∝ 1/B(s)` changes the mass flux, and `L_TR` is prescribed per thread. The release solver
already carries a prescribed flux tube (`grid.B_i` / `grid.B_iph`, area factors in the flux
differencing); `model_column` simply never sets it. **The cheap version is an ensemble spanning the B
range a real magnetogram produces** — quiet Sun ~5–10 G, network ~50 G, plage/AR ~500 G — with
PFSS-derived `A(s)`, turning the recipe into `{N_i k_B T}(T_TR, L_TR, T_ch, B)`. A magnetogram is
then useful for telling you *which B range to span*, which is a histogram, not a simulation.

**Do not run magnetogram-driven threads for this paper.** The Bin 2 numbers are local thermodynamics;
flux-tube geometry changes their size but not whether the assumptions hold. It would broaden the
parameter study without sharpening the claim, open a large referee surface (which magnetogram? which
PFSS? open versus closed? source-surface height?), and it is the classic route by which a methods
project becomes a software project — the §2.5 failure mode exactly.

#### Tier 3 — actual SWMF coupling: **not on the path.**

One paragraph of future work, or the next person's project.

#### Why this persuades

The paper shape is: *"The closure underlying [widely used model class] departs by X% under conditions
Y; here is the corrected prescription."* The target is a published, used, cited set of equations
rather than a straw man; the result is a number rather than "our model is more self-consistent"; the
deliverable is directly adoptable; and **only this work can produce it**, because it needs a resolved
column with ionization thermodynamics at exactly that height.

**Risk:** if the audit finds a small discrepancy, the result is null. Mitigation — the
ionization-energy term cannot be zero across 6–22 kK, where the column runs **x = 0.24 → 1.00**
(measured, §2.7), and Tier 0 reveals the signal size in *two* independent channels before any
significant investment. **The dominant mitigation is now §1.4 F2:** the closure's founding paper
limits its own validity to T > 500,000 K, so the scope argument does not depend on the size of the
discrepancy at all.

**First, what the outer datum *is*, because this determines how the boundary studies are written up.**
`kOuterPRef` is not an initial-condition residue that happens to be convenient — **it is precisely the
quantity Sokolov's closure emits.** Our top face is the TFLM thread's bottom face, and the thread
hands down exactly two things, which is exactly what the release consumes:

| AWSoM-R side | What is handed down | Our side |
| --- | --- | --- |
| Assumption 2: `{N_i T}` constant across the TR, pulled out of the integrand "as the temperature gradient scale within the TR is much shorter than the barometric scale" | one pressure value | `kOuterPRef` |
| Eq. (22): `κ₀T_TR^{5/2}(∂T/∂s) = {N_i k_B T}·U_heat` — §1.3 calls it "a Robin-type condition relating T and ∂T/∂s" | one conduction relation, *not* a temperature | `outer_conduction_temperature` (22 kK wall), hydro ghost temperature left free |
| Assumption 3: steady state — no time derivative, no enthalpy flux, no mass flux | both are constants | `p_ref` fixed, `T_wall` fixed |

Two consequences for the write-up:

- **This is not over-specified.** Lie splitting separates the operators: the subsonic-outflow
  hyperbolic operator admits one incoming characteristic and takes one datum (the pressure); physical
  conduction is a separate parabolic operator and takes its own (the wall). The count is exact. That
  the hydro ghost temperature floats is not a compromise — under Eq. (22) the top temperature is not
  independent data.
- **"Fixed" is the internally consistent choice, not the lazy one.** The closure being audited is
  *static* (assumption 3). Holding the datum constant reproduces it faithfully; letting the pressure
  drift replaces the audited object with something that does not exist in the host. The honest
  limitation to state in the paper is narrower and should be phrased on the *consistency* side:
  *the overlying TR does not respond to the column within the run, which matches the static character
  of the closure under audit; relaxing it belongs to real coupling (Tier 3).*

**Boundary studies are the other core asset.** These results are unusual in the literature and are
exactly what a coupling paper needs:

- The external hyperbolic datum belongs in the **pressure**, not the temperature. Swapping them
  withholds the top cell's hydrostatic support: the column collapses at the Mach cap, loses half its
  mass by 1000 s, and settles into a stationary 23.7 kK fully ionized downdraft — **while reporting
  clean termination and zero Riemann fallbacks.** A silent, plausible-looking catastrophic failure
  mode is precisely what a coupling referee wants documented.
- A hydrostatically consistent Neumann condition (dp/ds = −ρg on the live top cell) is stable for
  4000 s but leaves the pressure level unpinned: +5.8% drift still rising, evaporation signal down
  ~19×.
- The anchor is a **measurement** requirement, not a stability one: the two-point wall-temperature
  exponent survives (2.80 vs. 2.91) while every absolute velocity is ~20× too small and the top and
  base mass fluxes take opposite signs.
- Lower boundary: a gravity-aware characteristic condition changes the hydrostatic residual not at
  all and the evaporation velocity by ≤0.4%, **but recovers ~17% of the 1000 s column mass loss** —
  nearly irrelevant to the velocity field, material to the mass budget.

**Prerequisites.** Radiative losses (§5.1); domain extended upward far enough to overlap a TFLM
thread; matched N=1000 long-duration convergence. Target ApJ/ApJS.

### 5.5 Two measurement cautions

The per-assumption programme is the Tier 1 and Tier 2 tables in §5.4. Two traps to avoid there:

- **L_TR.** Do not compare our 15.7 km thermal boundary-layer thickness against the prescribed
  ≈1 Mm. They are different quantities: our domain tops at 22 kK and never reaches coronal
  temperatures, while L_TR spans the full TR. **Do not quote the 64× ratio.** The comparison only
  becomes legitimate after the Tier 2 domain extension.
- **The merge-location mass flux.** Sokolov's "negligible" argument covers the *location* of the
  merge point and nothing else (§1.4). Measuring the mass flux *across* that merge is the one test
  the argument does not address, and it should be reported separately from the assumption audit
  rather than folded into it.

### 5.6 Sequence

```
Tier 0: quadrature only (days)   ← FIRST. Decides whether the rest is worth doing.
        (i)  T_ch sensitivity of Eq.21   — is the signal large?
        (ii) T_ch sensitivity of the full-ionization error (C7)
        (iii) Eq.21 output vs our own top pressure — does the closure reproduce C7 at the merge?
        also reads out: reservoir or TR broadening?
   │
   ├─ signal large ──► Tier 1  (weeks)      p_ref is an INPUT here; the column sits BELOW
   │   AND merge-      │                    Eq.21's integration range and audits underneath
   │   region is       │  CL2012 radiative losses + option-1 background heating  (§5.1 — ONE change)
   │   the cause       │  audit assumptions 4 / 5 / 6 / 7
   │                   │  assumption 8 (tau_eq) + Damkohler map  ← post-processing, same pass
   │                   │  frozen-x bracketing run          (mandatory)
   │                   │  ISO_OUTER_PREF + 3-5 pt pressure sweep  (mandatory — two reasons)
   │                   │  Longcope 1/3-exponent sweep  ← the cross-code comparison
   │                   │
   │                   └──► Tier 2  (months)   only here does p become an OUTPUT
   │                          extend domain through T_TR + split T_e/T_i  (one change)
   │                          sweep (T_TR, p_top) in 2-D, read off (L_TR, T_ch), INVERT
   │                          deliver corrected {N_i k_B T}(T_TR, L_TR, T_ch)  ──► submit
   │
   └─ signal small ──► fall back to Angle A; the audit becomes one section of it
       OR broadening
       dominates

in parallel, independent of the above:
   Action 1a: verify the §1.4 reservoir numbers in the primary PDFs   hours, highest value/cost
   Angle B research note (precision criterion)                       ~2–4 weeks, low risk
   Angle A acoustic contrast case                                    ~1 week, decides Angle A viability
```

**Tier 3 (real coupling) is not on the path.**

**The gate at every stage is §2.5's one-sentence test**, not completeness: the question is never "is
the model more self-consistent yet" but "do I have a number of the form *[quantity] departs by
[amount] under [conditions]*." Tier 0 exists to answer that as cheaply as possible; Tier 1 exists to
make the number defensible; Tier 2 exists to turn it into something adoptable.

---

## 6. The instantaneous-Saha question

See also the standing note in the root `CLAUDE.md`. **The new axis makes this substantially easier to
defend, for three reasons that did not exist under the old framing.**

### 6.0 Why the new axis helps

**① The comparison baseline changed — this is the big one.** Under the old framing the implicit
comparison was RADYN/Bifrost, which have full NEQ; against them equilibrium Saha is a *degradation*
and we were defending "why not use the better thing." The comparison is now **Sokolov's analytical
closure, which contains no ionization physics whatsoever** — N_e = N_i hard-coded, no x, no
ionization energy (§1.3). Against that baseline equilibrium Saha is a large *improvement*. The
referee question shifts from

> "Why Saha instead of NEQ?"  →  "Does adding equilibrium ionization to a closure with zero ionization
> physics improve it?"

which is a far easier yes.

**② Cost becomes a first-class design argument rather than an excuse.** The host is named
AWSoM-**R**ealtime and Sokolov's stated motivation is that the TR is the computational bottleneck
(§1.1). "We chose the cheapest closure that retains the dominant ionization energetics" is then the
*requirement*, not a compromise. NEQ would be disqualifying, not better.

**③ Bracketing becomes decisive rather than defensive.** Previously the frozen-x bracket put an error
bar on an admitted weakness. Now it brackets **the correction to Sokolov's closure**. If the
correction to {N_i k_B T} is, say, +40% under Saha and +25% under frozen-x, then it is 25–40% —
**both bounds far from zero, so the finding survives regardless of the ionization treatment.** That
demotes the Saha uncertainty from a threat to a nuisance parameter.

**What got harder in exchange:** the headline number now depends on the ionization treatment, so
**bracketing is mandatory, not optional** (promoted in §5.4 Tier 1). The exposure is at least
concentrated rather than diffuse: Sokolov's TR runs to 10⁶ K where hydrogen is fully ionized and Saha
is trivially exact; the question bites only in the 6–30 kK merge region, which is exactly our domain
and exactly what we can measure.

**A further argument available only for the steady deliverable:** Eqs. (19)–(22) are *static*. If
what we deliver is a corrected quasi-steady relation {N_i k_B T}(T_TR, L_TR), then equilibrium
ionization is the **internally consistent** choice — using NEQ would contradict the steady-state
assumption of the very thing being corrected. **This holds only for the steady recipe, not for
time-dependent evaporation claims.** Keep the honest §6.1 table for those.

### 6.1 Where it actually fails

The criterion is **not** τ_ion versus Δt but τ_ion versus the timescale on which a fluid element's own
(ρ, T) change.

| Region | τ_ion (Carlsson & Stein 2002) | Local evolution timescale here | Verdict |
| --- | --- | --- | --- |
| Top of domain, fully ionized | — | — | Safe: x → 1 leaves no freedom |
| Transition-region base | ∼10² s | conductive heating ∼10²–10³ s | Marginal |
| Mid-chromosphere 1600–1850 km | 10³–10⁵ s | ∼10³–10⁴ s | **Fails** |

"Our problem is smooth and subsonic" evades the shock form of the critique but not the long
mid-chromospheric τ_ion. Note that this failure region coincides with the ~13% monotone flux spread
over 1600–1850 km still relaxing at 4000 s; the two may not be independent, and that is worth
checking.

### 6.2 Reframe — free, and the highest-value writing change

`paper.tex` currently presents the frozen Riemann solve and the per-step Saha decode as **two physical
choices in mutual tension**. Honest, but it hands a referee two targets. The accurate framing:

> **The frozen Riemann solve is a numerical approximation, not a physical assertion.** It approximates
> the *wave structure* at a face, not composition transport. Solving the Riemann problem at a local
> effective γ and decoding through the full EOS is standard in general-EOS Godunov schemes — Rempel
> (2016) uses an equilibrium-ionization EOS in exactly this spirit.
>
> **The only physical closure assertion is x = x_Saha(ρ, T).**

The e_int = p/(5/3−1) + e_ion split supports this directly: an exact algebraic identity, not an
assumption.

### 6.3 Four routes to justification, ranked

**① The error is one-sided and conservative.** Leenaarts (2007) and Przybylski (2022) both find that
under NEQ conditions ionization is a *less effective* energy buffer. In our configuration — monotone
conductive heating, no shocks — LTE lets x follow T immediately, energy goes into ionization rather
than translation, T rises less, pressure is lower, **evaporation is weaker**. Under NEQ, x lags and
evaporation would be stronger.

> **Instantaneous Saha biases our evaporation rate low — the conservative direction.**

Opposite in sign to the float32 bias (~30% high); state both together. **Caveat:** Nóbrega-Siverio
et al. (2019) find LTE *underestimates* ionization in flux emergence, so the sign is configuration
dependent. Always qualify with "in the monotone-heating, shock-free configuration reported here."

**② A posteriori Damköhler diagnostic — best return on effort, pure post-processing.**
`../studies/eos-ionization/ionization_recombination_formulation.md` already contains the full rate
formulation *including three-body recombination* — the inverse process required to reach the Saha
limit at the dense base. Voronov (1997) and Hummer (1994) are in `../writeup-overleaf/papers/` and in
`reference.bib`. Rate functions are implemented in `src/two_fluid/integrators.cpp` (the release
single-fluid path carries no rate coefficients). Compute τ_ion(ρ,T) per cell from release output,
pair with the fluid-element |D ln x / Dt|, and plot **Da = τ_ion · |D ln x / Dt|** as a height–time
map. Da ≪ 1 is a quantitative exemption; Da > 1 is an honest error map. Kerr et al. (2019) did exactly
this for flare Mg II — **the methodological template.**

**③ Bracket the two limits.** Run a frozen-x variant (x locked at IC values) as the τ_ion → ∞ bound;
the existing Saha run is τ_ion → 0. The true NEQ answer lies between; the spread is the error bar.
Small implementation: bypass Saha in the EOS decode, use stored x. Wedemeyer-Böhm et al. (2011)
endorse exactly this framing for calcium — timescales "too long for instantaneous equilibrium and, on
the other hand, not long enough to warrant an assumption of a constant ionization fraction."

**④ Finite-rate ionization — not in this paper.** `apply_ionization_stage` exists in the two-fluid
path and the physics is written down, but it turns three conserved rows into four and breaks the
release architecture. Per the root `CLAUDE.md`, "a modelling project, not a correction." Defer.

### 6.4 Precedent

*Codes using instantaneous equilibrium ionization by default:* Rempel (2016) MURaM coronal extension
("an equilibrium ionization equation of state", 173 citations); Bifrost, MURaM, Stagger, CO⁵BOLD all
default to LTE EOS tables with NEQ hydrogen optional; the MHD control runs in the Snow & Hillier PIP
studies.

*Papers arguing affirmatively that equilibrium is adequate:* **Brown (1973)** — in optical flares,
spontaneous recombination to the second level plus photoionization dominate, so ionization adjusts on
the recombination timescale, *shorter* than the flare heating timescale, hence "the ionisation is
given by a simple LTE-modification of Saha's equation at the instantaneous electron temperature."
**Kerr et al. (2019)** — τ_relax < 0.1 s for most of a flare, so "the equilibrium solution is an
adequate approximation." **Varonov et al. (2023)** — "the quiet solar atmosphere is in LTE, justifying
this theoretical study" (but active regions are not; do not quote half of it).

*Opposing citations that must be engaged in §2 where the closure is introduced, not deferred to a
closing paragraph:* Carlsson & Stein (2002); Leenaarts (2006, 2007); Martínez-Sykora (2019); Anusha
et al. (2021); Osborne et al. (2026).

### 6.5 Writing consequences

- Reframe per §6.2.
- Replace the soft Introduction phrase "a controlled intermediate model" with a quantitative
  applicability domain, the one-sided error direction from ①, and the bracketing interval from ③.
- Move the NEQ counter-citations into §2.

---

## 7. Honest position assessment

| Dimension | This release | Verdict against the AWSoM-R gap |
| --- | --- | --- |
| Position in the atmosphere | 1600–2153 km, 6→22 kK | ✅ Exactly Sokolov's uncertain merge region |
| Ionization thermodynamics | Saha equilibrium, ionization energy in ℰ, CRASH Γ₁ table | ✅ The physics the analytical closure lacks |
| Riemann solver | Frozen-composition exact Riemann, ionization energy in a passive E₀ | ✅ No prior chromospheric art found |
| Cost | 0.199 core-s per physical second | ⚠️ ~21× over the AWSoM-R real-time budget; lever identified (§3.3) |
| Radiative losses | **Absent** | ❌ Blocks engagement with Eqs. (21)–(22) — schedule first |
| Domain top | 22 kK, never reaches coronal T | ❌ Must extend upward to overlap a TFLM thread |
| Long-duration convergence | N=1000 stopped at 2080 s; <10% not established | ❌ Blocking for any quantitative rate |
| Temporal order | First order (Lie split) | Defensible given \|u\| ≪ c; keep it stated |
| Mass-supply channel | Thermal only | ⚠️ Scope-limited by construction (§2.1) |
| Alfvén waves — transport | Absent | ✅ Inherited: the host applies {S_A/B} as a BC at our domain top and computes nothing below it (§2.6a) |
| Alfvén waves — momentum (ponderomotive, wave pressure) | Absent | ⚠️ Sokolov also drops it, but this is the type-II-spicule/FIP channel; scope explicitly (§2.6b) |
| Alfvén waves — as the chromospheric heat source | Absent, and currently harmless | ❌ **Becomes blocking the moment Λ_R is added** — a sink with no source (§2.6c, §5.1) |
| Temperatures | Single T | ✅ for this domain: τ_eq ≈ 1–10 ms vs 10²–10³ s dynamics (§2.7). ❌ at the Tier 2 interface, where the host wants three |
| Nature of the contribution | Bin 2 (quantified assessment) + Bin 3 (method/criterion) | ✅ the modal paper shape in this subfield — **provided §2.5's one-sentence test is passed.** A "coupling-ready module" claim alone is ❌ |

---

## 8. Immediate actions

| # | Action | Cost | Unblocks |
| --- | --- | --- | --- |
| **1** | **Tier 0, now three items over the same T_ch sweep: (i) Eq. (21) with a CHIANTI table → {N_i k_B T} and heat flux; (ii) the full-ionization error in Eq. (20) against Model C7's x = 0.75→0.97 (§2.7) — C7 table already in `scenarios/model_c7.cpp`, no CHIANTI needed; (iii) Eq. (21)'s output against the column's own top pressure ({N_i k_B T} ≈ 0.0051 Pa) — does the closure reproduce a semi-empirical atmosphere at the merge height? Compare (i) against the stellar-paper overestimate; read all three against the §2.8 attribution question.** | **~1 day, no simulation** | **sizes the signal in three independent channels; (iii) also sets the centre of action 4b's sweep** |
| ~~1a~~ | ~~Verify the §1.4 reservoir numbers in the primary PDFs~~ — **DONE 2026-09-06.** All four claims resolved against primary text; van der Holst (2014) established as the source, four competing rationales catalogued, two corrections applied (§1.4) | — | done |
| ~~1b~~ | ~~Obtain Lionello, Linker & Mikić (2009)~~ — **DONE 2026-09-06**, filed and read. It *was* a finding in itself: factor-2 evidence base for a factor-10-to-250 extrapolation, and an explicit self-imposed validity floor at 500,000 K that excludes the region where the error is reported (§1.4 F1–F5). **This is now the second-strongest item in the Introduction after the reservoir itself** | — | done |
| 2 | Add `vanderHolst2014` (**the §1.4 anchor**), `Sokolov2013`, `Sokolov2021` + Lionello/Downs + the §4 second-revision list to `reference.bib`; rewrite the Introduction around §1, leading with the §1.4 factor-of-ten reservoir and the inconsistent-rationale table, then quoting the §1.4 sentence | hours | everything |
| ~~3~~ | ~~τ_eq(ρ,T) map of §2.7~~ — **DONE 2026-09-06**, `util/tau_ei_equilibration.py`; §2.7 now carries measured values and one claim was retracted. **Damköhler map (§6.3②) still outstanding** — same output, same kind of pass | days | §6 |
| ~~3a~~ | ~~Resolve the x = 0.273 flag~~ — **DONE 2026-09-06.** Code and IC both correct; LTE Saha reproduces Model C7's NLTE ionization at the domain base to 14 % (**new Angle A asset**). The "x ~ 10⁻²" was a mid-chromospheric figure misapplied in three places, now fixed. One headline weakened (assumption 6 is a 3–25 % error at T_ch, not orders of magnitude), one new Tier 0 quadrature gained (§2.7) | — | done |
| 4 | Frozen-x bracketing run (τ_ion → ∞ bound) — **now mandatory** | days | §5.4 Tier 1, §6 |
| **4b** | **Add an `ISO_OUTER_PREF` override so the outer pressure datum can be set explicitly (falls back to the present IC-derived capture when unset), then run a 3–5 point pressure sweep. Mandatory at Tier 1 for two independent reasons: the audit numbers need a sensitivity, and the Tier 2 recipe cannot be inverted from a single-value slice.** Do **not** loosen the boundary instead (§5.4) | **hours to add, ~1 h to sweep** | **§5.4 Tier 1 *and* Tier 2; the only unswept axis of the deliverable** |
| 5 | Acoustic / weak-shock contrast case, frozen vs. equilibrium flux | ~1 week | **Angle A viability** |
| 6 | Reframe `paper.tex` §2.4 and Discussion per §6.2, §6.5; apply the §2.4 retractions; add the §2.6b scope sentence on the ponderomotive/spicule channel | hours | Angle A |
| 7 | Implement Carlsson & Leenaarts (2012) radiative losses **together with the option-1 background heating of §5.1 — one change, not two** | weeks | **Tier 1** |
| 7a | Longcope (2014) 1/3-exponent heat-flux sweep | days, after 7 | §5.2's cross-code comparison blocker |
| 8 | Verify with the SWMF group: real TFLM thread count, per-thread cost, and the ≈1 s host timestep. **Add: does the low-height density overestimate trace to the reservoir or to the TR broadening?** | conversation | §3.2, §3.3, §2.8 |
| 8a | **Move the outer conductive reservoir from 22 kK to 50 kK (+8 km on C7)** — reaches AWSoM's own inner-boundary temperature, the cheapest point of contact with the host in this plan (§5.4 Tier 2) | days | a direct comparison against the host's boundary state |
| 9 | **Step A** — extend +50 km to 10⁵ K: CL2012→CHIANTI crossover (= the assumption-5 audit), 44× κ jump, ≲1 km TR mesh. De-risks Step B (§5.4 Tier 2) | weeks | Tier 2 |
| 9b | **Step B** — extend +1.0 Mm to T_TR ≈ 4×10⁵ K **and split T_e/T_i in the same change**; sweep T_TR as the control parameter | months | the deliverable |
| 10 | Draft the precision research note | ~2 weeks | Angle B |
| 11 | Finish the matched N=1000 4000 s run | compute | Angle A |

Items 1, 1b, 2–6, 3a and 8 need no long compute and can proceed immediately. **Item 1 comes first** —
it is a day of quadrature that determines whether Tiers 1–2 are worth funding.

**Item 3a is the cheapest item with real downside risk.** The measured base ionization
(x = 0.273 at 6631 K) is more than an order of magnitude above the "x ~ 10⁻²" this document asserts
in §1.3, §2.5 and §5.4, and the ionization-energy-reservoir argument — the mitigation for the whole
programme's null-result risk (§9) — is stated in terms of that figure. Either the figure is being
misapplied to our layer, or something is off in the IC or the diagnostic. **An hour now, versus a
retraction after Tier 1.**

Note that items 3 and 8 both got *cheaper* in this revision rather than more expensive: assumption 8
is answered by post-processing already-scheduled output, and the cross-code comparison (7a) needs no
external code.

---

## 9. Standing risks

- **The gap could close upstream.** The entire programme rests on §1.6, a judgement with a date on
  it. **Re-check before every submission and at least every six months.** If Sokolov or the PSI group
  publishes a revised chromospheric closure, re-plan rather than press on.
- **Tier 0 could return a small signal**, making the audit a null result. **Substantially reduced on
  2026-09-06 by reading Lionello et al. (2009).** The programme no longer depends on Tier 0 returning
  a large number, because the scope-of-validity argument stands independently: both approximations
  are validated by their own founding paper only above ~500,000 K, and the reported error is below
  that (§1.4 F2, §2.8). A large Tier 0 signal strengthens the paper; a small one now costs a section
  rather than the thesis. Secondary mitigation now on measured footing: the ionization-energy term
  cannot vanish across 6–22 kK, where the column runs **x = 0.24 → 1.00** (§2.7; the earlier
  "~10⁻² → ~1" was a mid-chromospheric figure misapplied, corrected 2026-09-06). Δx = 0.73 rather
  than 0.99 — a 27 % reduction in a term that was never marginal. **Third, independent mitigation
  added the same day:** Tier 0 now also sweeps the ionization channel, where the assumption-6 error
  varies ~8× across T_ch's own factor-of-two freedom.
- **Framing drift.** If the paper reads as "AWSoM-R is wrong" rather than "this closure class needs a
  partial-ionization correction," it reaches one group instead of three and antagonizes the group
  this work serves (§2.1, §2.2).
- **No radiative sink — and, once there is one, no source.** The first scheduled physics item rather
  than a background worry (§5.1). The second revision adds the harder half: **Λ_R must ship with a
  background heating term in the same change**, or the chromosphere cools away on exactly the
  10²–10³ s timescale the long runs depend on (§2.6c). Implementing the sink alone would look like
  progress and be a regression.
- **Cost is ~21× over the host's real-time budget**, with the lever identified but not implemented.
  Do not claim coupling readiness until this is addressed or explicitly scoped.
- **The thermal channel may not be the dominant mass supply** (spicules, §2.1). Scope the claim
  explicitly; do not let it drift.
- **N=500 long-duration mass flux is not established as <10% converged.** A *resolution* limitation,
  entirely separate from storage precision. Do not let recaps merge the two.
- **The localized first-interior-face artefact is precision independent and open**
  (|f_diff|/f_total ≈ 16% at face 0, cell-0 ρV overshoot ≈ 1.28). Never describe it as solved.
- **No positive noise floor established for m/s-scale chromospheric velocities.**
- **The CFL acceptance matrix has not been re-run in double.** CFL 0.50 stands on its original
  float32-era evidence; the residual-velocity convergence mechanism was retired, not the acceptance.
- **The 2.1×10⁴ thread count is my estimate from 1.4°, not a published number** (action 8).

Added by the second revision:

- **The contribution could collapse into a software deliverable.** The failure mode is writing the
  paper around "consistent with AWSoM-R and couplable" instead of around a measured departure. It is
  detectable in one line — §2.5's one-sentence test on the abstract's second sentence — and it is
  detectable *early*, which is why Tier 0 comes first. **Re-apply the test at every draft.**
- ~~The §1.4 reservoir evidence is second-hand.~~ **Closed 2026-09-06** — verified verbatim against
  van der Holst et al. (2014), Lionello et al. (2009) and five other primary sources; see §1.4. The
  framing must now carry the group's four inconsistent rationales rather than asserting one.
  **Residual:** Mok et al. (2005), the onward citation Lionello gives for the pressure-overestimate
  claim, is unchecked — one more hop down a chain that is already three deep, and probably not worth
  the trip unless a referee pulls on it.
- **Reservoir versus TR broadening cannot be attributed from the literature.** Checked and settled as
  a *reading* question on 2026-09-06: the authors name **both** causes themselves, in adjacent
  paragraphs, without separating them (Sachdeva et al. 2019 §4.2; Shi et al. 2022 §3.4). Only one is
  ours. If broadening dominates, the headline shrinks. **Tier 0 is now the only route** — plus a
  direct question to the SWMF group (action 8).
- ~~The τ_eq table of §2.7 is an order-of-magnitude estimate.~~ **Closed 2026-09-06** — measured per
  cell from release `.gamma_diag` output. **It changed the claim:** the "four to five orders of
  margin against any dynamical timescale" phrasing was wrong and is retracted; the worst-case margin
  is 1.43× against Δt and 2.85× against cell sound crossing, both at the top cell, with the 4.4
  orders holding only against the evaporation timescale. The conclusion (defer to Tier 2) survives;
  the sentence supporting it did not.
- ~~The column's base ionization does not match what this document says elsewhere.~~ **Closed
  2026-09-06** (§2.7). The code and the IC are both right — LTE Saha reproduces Model C7's NLTE
  ionization at the domain base to 14 %, which is now an Angle A asset. The "x ~ 10⁻²" figure was
  C7's mid-chromospheric value (h ≈ 1000–1400 km) misapplied to two other places. **One correction
  cuts against us:** Sokolov's full-ionization identity is a 3–25 % error at T_ch, not the
  order-of-magnitude failure §1.3 claimed — a weaker headline, offset by a new Tier 0 quadrature
  (the error varies ~8× across T_ch's own freedom). **Residual:** the Saha-vs-C7 comparison is clean
  at 1617 km but not near the temperature minimum, where low-FIP metals dominate n_e; do not extend
  it downward without redoing that estimate.
- **The Alfvén-wave omission is safe only while the model has no chromospheric energy budget.** The
  §2.6a inheritance argument covers wave *transport*; it does not license running a chromosphere with
  a radiative sink and no heating. Do not let the "we inherit this from the host" sentence migrate
  from §2.6a to §2.6c — they are different claims and only the first one is free.
