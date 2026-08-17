# The coronal conductive-flux boundary condition $q(T)$

This note explains the upper-boundary closure $q(T)$ used by `model_c7` / `model_gentle`, and in particular what "imposed Neumann flux" means. All symbols are those of the writeup (`paper.tex`, §Boundary and initial conditions).

## The term it sets

In the charged-fluid energy row of the six-equation system the field-aligned heat conduction appears as $B\,\partial_s\!\big[(\kappa_e+\kappa_i)/B\,\partial_s T_i\big]$. The quantity being differentiated is the conductive heat flux carried along the field, $F_c = -(\kappa_e+\kappa_i)\,\partial_s T_i$. At the transition-region top the electron conductivity dominates and is Spitzer, $\kappa_e\propto T^{5/2}$, so the flux there is $F_c = -\kappa_0 T^{5/2}\,dT/ds$.

The outer boundary sits at $h\approx2153$ km, $T\approx2.3\times10^4$ K (just above the $\kappa_e/\kappa_n$ crossover) — the natural top of a chromosphere-only column. Across this constant-pressure transition-region base, pressure is continuous and $T$ jumps to the coronal value, so $\rho$ and $T$ there are taken Neumann, and the velocity is a mass-flux outflow capped at $M\le0.05$ of the transition-region sound speed $c_s=\sqrt{2\gamma k_B T/m_p}$. The corona itself is not in the box; its only influence is the heat it conducts down — and that enters as $q$.

## What "imposed Neumann flux" means

A heat-conduction equation needs exactly one boundary condition at each end of the domain — the solver has to be told *something* at the top face, because there is no cell above it to conduct with. There are two standard kinds of condition, and they differ only in *which* quantity you pin down at the edge:

- **Dirichlet** — pin the **temperature** $T_i$ at the face (e.g. "the top is at $10^6$ K"). The heat flowing through the face is then whatever the resulting gradient produces; you do not control it directly.
- **Neumann** — pin the **gradient** $\partial_s T_i$ at the face instead. For heat conduction the flux is $F_c=-(\kappa_e+\kappa_i)\,\partial_s T_i$ — directly proportional to that gradient — so pinning the gradient *is the same as pinning the heat flux through the face*. The face temperature is then left free, to be whatever is consistent with that flux.

A physical picture for the top cell of the column:

- A **Dirichlet** condition is like clamping the end of a metal rod to a thermostat held at a fixed temperature: the temperature is forced, and the heat current into the rod adjusts to whatever maintains it.
- A **Neumann (fixed-flux)** condition is like bolting a fixed-wattage heater to the end of the rod: you control the **heat delivered** (say $90$ W m$^{-2}$), and the end temperature settles to whatever results.

"**Imposed** Neumann flux" therefore means: each timestep we *set the conductive heat flux through the outer face to a chosen value* $q$ — the fixed-wattage heater — and let the top-cell temperature $T_i$ solve to whatever is consistent. We do *not* tell the solver a coronal temperature, and we do *not* read the flux off the interior; we hand it the heat current directly. The corona is not inside the box at all, so its entire influence on the column is summarised by this one number: the heat it conducts down through the top face. Numerically it enters as the outer entry of the implicit conduction tridiagonal — a known flux $q$, with $T_i$ at the top cell solved for.

## Why a fixed flux and not a fixed temperature

Because $\kappa_e\propto T^{5/2}$, a Dirichlet hot temperature on the coarse grid would conduct far too strongly into the top cell (the $T^{5/2}$ makes the conductive flux explode), over-pressurise it, and drive a spurious supersonic downflow. Fixing the flux $q$ instead sets exactly the energy the corona delivers, independent of the under-resolved gradient — which is why the writeup imposes a Neumann *flux*, "a fixed flux, not a fixed hot temperature."

This $q$ is re-radiated at and below the transition-region base, so it is imposed **together with the radiative-loss sink** — the two form the heating/cooling pair that holds the chromosphere top against ionization cooling; without the sink the imposed flux would have nowhere to go and the top would run away.

## The $q(T)$ actually used in `gentle_v3_evolution`

That run extends the column into a **resolved corona** (top at $h\approx11.6$ Mm), so $q$ is *not* the RTV static-loop value $(2/7)\kappa_0 T_{\rm cor}^{7/2}/L$ used by the standard chromosphere-only `model_c7`. Instead $q$ is the flux the Model C7 atmosphere itself carries across the chosen top, taken from its own temperature gradient there:

$$q_0 \;=\; \kappa_0\,T_{\rm top}^{5/2}\,\left.\frac{dT}{ds}\right|_{\rm top},$$

with $\kappa_0=10^{-11}$ W m$^{-1}$ K$^{-7/2}$, and $T_{\rm top}$, $dT/ds$ read from the two highest Avrett & Loeser (2008) Table-26 rows used — at $h\approx8.97$ and $11.6$ Mm, where $T\approx0.78$ and $0.86$ MK:

$$q_0 \;=\; 10^{-11}\,\big(8.65\times10^5\big)^{5/2}\;\frac{(8.65-7.78)\times10^5}{(11.60-8.97)\times10^6}\;\approx\; 2.3\times10^{2}\ \mathrm{W\,m^{-2}}.$$

So the steady (pre-ramp) flux is $q\approx230$ W m$^{-2}$ directed downward. The gentle-evaporation driver then ramps it by a factor 3 with a cosine in time:

$$q(t) \;=\; q_0\,\big[\,1 + 2\,w(t)\,\big],\qquad w(t)=\begin{cases}0, & t<t_{\rm on}\\[2pt]\tfrac12\!\left(1-\cos\dfrac{\pi\,(t-t_{\rm on})}{\Delta}\right), & t_{\rm on}\le t\le t_{\rm on}+\Delta\\[4pt]1, & t>t_{\rm on}+\Delta\end{cases}$$

with onset $t_{\rm on}=3000$ s and ramp width $\Delta=300$ s (the run's `GENTLE_Q_ENHANCE=3`, `GENTLE_Q_T_ON=3000`, `GENTLE_Q_RAMP=300`). Thus $q$ holds at $\approx230$ W m$^{-2}$ through relaxation, then rises smoothly to $3\times230\approx690$ W m$^{-2}$ over $3000$–$3300$ s and holds — the larger downward conductive flux that drives the evaporation.

**Provenance of the ramp.** This time profile is *not* from the literature — it is a code construct ([scenarios/model_c7.cpp](../../../scenarios/model_c7.cpp), `model_c7_update_bc`). The envelope $w(t)=\tfrac12(1-\cos\pi x)$ is a raised-cosine (Hann) smoothstep, chosen only because its value *and* slope are continuous at both ends, so the flux turns on without a kink and launches no spurious acoustic transient (a step or a linear ramp would). The shape is arbitrary; the *physics* requirement is merely that the ramp be **slow** — $\Delta$ much longer than the TR/sound-crossing adjustment time — so the upflow stays subsonic and the response is gentle (A&S 1978; an impulsive onset would push toward the explosive regime of Fisher, Canfield & McClymont 1985). Only $q_0$ (the magnitude) and the slow-ramp requirement come from physics; the cosine envelope is a modeling choice.

### References
- **Spitzer (1962)** — the $\kappa\propto T^{5/2}$ field-aligned electron conductivity giving $F_c=-\kappa_0 T^{5/2}\,dT/ds$.
- **Johnston et al. (2020)** — the SI coefficient $\kappa_0=10^{-11}$ W m$^{-1}$ K$^{-7/2}$.
- **Avrett & Loeser (2008), Table 26** — the Model C7 corona whose $T_{\rm top}$ and $dT/ds$ set $q_0$ (file `docs/supporting-papers/AvrettLoeser2008.pdf`).
- **Antiochos & Sturrock (1978)** — the conduction-driven, subsonic ("gentle") evaporation mechanism the ramped $q(t)$ targets. Verbatim, the abstract: *"conductive losses dominate radiative cooling … the evaporative velocities are small compared with the sound speed"*; and §I (p. 1138): evaporation *"driven by the large conductive heat flux from the high-temperature flare plasma … would involve only small pressure gradients and velocities … this later, 'gentle' evaporation process."* (A&S model the flare *decay* phase; we instead ramp the coronal flux, but the mechanism — conduction-driven, $v\ll c_s$ — is the same.)
- **Fisher, Canfield & McClymont (1985)** — the explicit *gentle-vs-explosive threshold* and its velocity bound ($v_{\max}\lesssim0.3\,v^\*$ for gentle, vs up to $\sim\!2.35\,c_s$ for explosive); this is the source for the "regime" classification, not A&S.
- **Rosner, Tucker & Vaiana (1978)** — the alternative static-loop magnitude $q=(2/7)\kappa_0 T_{\rm cor}^{7/2}/L\approx90$ W m$^{-2}$ used by the standard `model_c7`.
