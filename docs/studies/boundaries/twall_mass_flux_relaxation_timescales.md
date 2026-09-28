# How the evaporation mass flux relaxes and decays in the $T_\mathrm{wall}$ sweep

## Notation

**Abbreviations.**
- **TR:** transition region, the thin layer at the top of the domain where $T$ rises from the $\sim 6.7$ kK chromosphere to the wall temperature.
- **Chromosphere top:** the base of the TR.
- **IC:** the initial condition, the C7 atmosphere at $t = 0$.
- **N=500/R4:** the release mesh, 661 cells with $4\times$ refinement near the top.

Typical values are for the sweep runs at $t = 4000$ s unless stated.

**The flux and the column**

| Symbol | Meaning | Units | Typical value |
|---|---|---|---|
| $F$ | Mass flux $\rho u$ through a cell face, positive upward. This is the "NU" of the question: `f_total` in the `.faceflux` sidecars and the mass-flux panel of the movie. | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | $2.7\times10^{-10}$ (22 kK) to $4.9\times10^{-9}$ (100 kK) |
| $F_\mathrm{top}$, $F_\mathrm{base}$ | $F$ at the top face (2153 km) and at the first interior face (1601 km). $F$ with no subscript means $F_\mathrm{top}$. | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | |
| $F(z)$, $F(z,t)$ | $F$ through the fixed height $z$ | $\mathrm{kg\,m^{-2}\,s^{-1}}$ | |
| $u$ | Flow velocity, positive upward (`V` in the figures) | $\mathrm{m\,s^{-1}}$ | |
| $s$, $z$ | $s$ is distance along the field line and $z$ is height. They are the same here: a vertical, constant-area column. | m or km | 1600–2153 km |
| $z_b$, $z_f$ | Height of the base (1600 km) and of the chromosphere top | km | $z_f \approx 2140$ km at $t = 0$, then it falls |
| $D$ | Column length, $2153 - 1600$ km | km | 553 |
| $m$ | Column mass above a point, $m = \int_z^{\mathrm{top}} \rho\,dz'$. It is used as a label that moves with the gas. | $\mathrm{kg\,m^{-2}}$ | 0 at the top, $M$ at the base |
| $M$, $\Delta M$ | Total column mass, and its deficit $\Delta M = M(0) - M$ | $\mathrm{kg\,m^{-2}}$ | $M \approx 2.1\times10^{-4}$ |
| $\Phi(t)$ | Evaporated mass: everything that has left through the top face since $t = 0$, $\Phi = \int_0^t F_\mathrm{top}\,dt'$ | $\mathrm{kg\,m^{-2}}$ | $1.6\times10^{-6}$ (22 kK) to $3.1\times10^{-5}$ (100 kK) |
| $L \equiv L_\mathrm{TR}$ | TR length: distance from the 7.5 kK isotherm up to the top-cell centre | km | 13 at $t = 0$; 15 (22 kK) to 49 (100 kK) at 4000 s |
| $L_{T_1}$ | The same length measured from a general marker isotherm $T_1$ | km | |
| $L_1$, $t_1$ | TR length at a reference time $t_1$ after the TR has relaxed; the starting point of the decay law (5) | km, s | |
| $M_\mathrm{TR}$ | Mass in the TR, above the 7.5 kK isotherm | $\mathrm{kg\,m^{-2}}$ | |
| $\mu$ | $M_\mathrm{TR}/(\rho_f L)$: the TR mass as a fraction of a layer of chromospheric density and the same length | – | 0.42 (22 kK) to 0.08 (100 kK) |

**Thermodynamics**

| Symbol | Meaning | Units | Typical value |
|---|---|---|---|
| $\rho$, $p$, $T$ | Density, pressure, temperature | $\mathrm{kg\,m^{-3}}$, Pa, K | |
| $x$ | Hydrogen ionization fraction (Saha equilibrium) | – | 0.27 at 1600 km, 0.63 at the chromosphere top, 1 in the TR |
| $e$ | Specific internal energy, thermal plus ionization | $\mathrm{J\,kg^{-1}}$ | |
| $h$ | Specific enthalpy, $h = \tfrac52\,p/\rho + x\chi_H/m_H$: the energy carried by each kg of flowing gas | $\mathrm{J\,kg^{-1}}$ | |
| $h_c$, $h_\mathrm{top}$, $\Delta h$ | $h$ of the chromospheric gas entering the TR; $h$ in the top cell; $\Delta h = h_\mathrm{top} - h_c$ | $\mathrm{J\,kg^{-1}}$ | $h_c \approx 1\times10^{9}$; $h_\mathrm{top}$ $2.2\times10^{9}$ (22 kK) to $5.4\times10^{9}$ (100 kK) |
| $T_\mathrm{wall}$ | Imposed conduction temperature at the top face (`ISO_T_TOP`); the swept parameter | K | 8–100 kK |
| $T_\mathrm{top}$ | Top-cell temperature, held within 0.1–0.3 % of $T_\mathrm{wall}$ by conduction | K | |
| $T_c$ | Chromospheric temperature just below the TR | K | $\approx 6.7$ kK |
| $T_1$ | Marker isotherm that defines the TR base | K | 7.5 kK |
| $T_0(m)$ | Initial ($t = 0$) temperature at column mass $m$ | K | |
| $p_\mathrm{top}$, $p_f$ | Pressure at the top, fixed by the back-pressure `kOuterPRef`, and at the chromosphere top | Pa | $p_\mathrm{top} = 0.0102856$ |
| $\rho_f$, $\rho_b$ | Density at the chromosphere top and at the base | $\mathrm{kg\,m^{-3}}$ | $1.1\times10^{-10}$, $9.7\times10^{-10}$ |
| $V_f$ | $F/\rho_f$: upflow speed at the chromosphere top | $\mathrm{m\,s^{-1}}$ | 2.5 (22 kK) to $\sim 40$ (100 kK) |
| $S$ | Specific entropy | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | |
| $c_p$ | Specific heat at constant pressure, ionization included | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | |
| $R$ | Specific gas constant $(1+x)k_B/m_H$, so that $p = \rho R T$ | $\mathrm{J\,kg^{-1}\,K^{-1}}$ | |
| $\Gamma_1$ | Equilibrium adiabatic index, $(\partial\ln p/\partial\ln\rho)_S$ | – | $\approx 1.09$ in the chromosphere, $5/3$ when fully ionized |
| $\nabla$ | $d\ln T/d\ln p$ along the actual column: how $T$ changes with pressure going up | – | $\approx -0.007$ on the IC ($T$ rises slightly with height) |
| $\nabla_\mathrm{ad}$ | $(\partial\ln T/\partial\ln p)_S$: the same slope for gas compressed or expanded adiabatically | – | $\approx 0.05$; ideal gas 0.4. It is small because ionization absorbs the energy. |
| $\delta$ | $-(\partial\ln\rho/\partial\ln T)_p$: fractional density rise per fractional cooling at fixed pressure | – | $\approx 4$; ideal gas 1. It is large because cooling also recombines hydrogen. |
| $\delta(\nabla_\mathrm{ad} - \nabla)$ | The combination that sets how much a lifted parcel shrinks | – | 0.20–0.23 |
| $\kappa$, $\kappa_0$, $\kappa_c$ | Physical conductivity (Spitzer electrons plus neutrals); the Spitzer coefficient in $\kappa = \kappa_0 T^{5/2}$; $\kappa$ at $T_c$ | $\mathrm{W\,m^{-1}\,K^{-1}}$ ($\kappa_0$: $\mathrm{W\,m^{-1}\,K^{-7/2}}$) | $\kappa_0 = 9.2\times10^{-12}$ |
| $q$ | Conductive heat flux $q = -\kappa\,\partial T/\partial s$, positive upward | $\mathrm{W\,m^{-2}}$ | |
| $q_\mathrm{wall}$ | Conductive heat flux into the domain through the top face (`q_face` in `.outercond`) | $\mathrm{W\,m^{-2}}$ | 0.33 (22 kK) to 23 (100 kK) |
| $\chi_H$, $m_H$, $k_B$ | Hydrogen ionization energy (13.6 eV), hydrogen mass, Boltzmann constant | | |
| $g$, $\phi$ | Gravity ($273.95\ \mathrm{m\,s^{-2}}$) and gravitational potential | | |

**Transition-region balance (stage 2)**

| Symbol | Meaning | Units | Typical value |
|---|---|---|---|
| $J$ | $F L_\mathrm{TR}$: the steady-TR invariant. It is also computed from the boundary state alone, by the integral (1). | $\mathrm{kg\,m^{-1}\,s^{-1}}$ | $4.4\times10^{-6}$ (22 kK) to $2.4\times10^{-4}$ (100 kK) |
| $F L/J$ | Measured $F L_\mathrm{TR}$ divided by the theoretical $J$; 1 means the TR is in steady balance | – | |
| $\mathrm{ER}$ | Energy ratio $F(h_\mathrm{top} - h_c)/q_\mathrm{wall}$: enthalpy carried out over heat conducted in. 1 is balance; $> 1$ means the TR is releasing stored energy. | – | |
| $\mathrm{ER}_\infty$, $A$ | Fit parameters of $\mathrm{ER}(t) = \mathrm{ER}_\infty + A\,e^{-t/\tau_\mathrm{TR}}$: the final value and the initial excess | – | |
| $\Delta x$ | Change of $x$ when chromospheric gas is heated into the TR | – | |

**Contraction (stage 3)**

| Symbol | Meaning | Units | Typical value |
|---|---|---|---|
| $W(z)$ | How far down the level at height $z$ moves per $\mathrm{kg\,m^{-2}}$ evaporated: $W(z) = \int_{z_b}^{z} \delta(\nabla_\mathrm{ad} - \nabla)\,(g/p)\,dz'$ | $\mathrm{m^3\,kg^{-1}}$ | |
| $K$ | $W$ at the chromosphere top: how far the TR base descends per unit evaporated mass | $\mathrm{m^3\,kg^{-1}}$ | $1.61\times10^{9}$ on the IC |
| $K_\mathrm{cur}$ | $K$ recomputed on the current state, which includes the large-displacement effect | $\mathrm{m^3\,kg^{-1}}$ | $2.3\times10^{9}$ (15 kK) to $0.9\times10^{9}$ (100 kK) |
| $K_\mathrm{eff}$ | $K/(1-\mu)$: the actual TR growth rate $dL/d\Phi$, including the mass the TR absorbs | $\mathrm{m^3\,kg^{-1}}$ | |
| $r$ | $-d\ln F/dt$: fractional decay rate of the flux | $\mathrm{s^{-1}}$ | $6\times10^{-5}$ to $1.2\times10^{-4}$, i.e. 3–6 % per 500 s |

**Waves and timescales**

| Symbol | Meaning | Units | Typical value |
|---|---|---|---|
| $c$, $c_\mathrm{eq}$, $c_0$ | Sound speed; equilibrium sound speed $\sqrt{\Gamma_1 p/\rho}$; $c$ at the base | $\mathrm{m\,s^{-1}}$ | $c_\mathrm{eq}$: 8.7 km/s (base) to 10.8 km/s (chromosphere top). The frozen $\sqrt{\tfrac53 p/\rho}$ is $\sim$20–25 % larger. |
| $t_\mathrm{ac}$ | One-way sound-crossing time of the column, $\int ds/c$ | s | 59 ($c_\mathrm{eq}$), 48 (frozen) |
| $H$, $H_f$ | Pressure scale height $p/(\rho g)$, and its value at the chromosphere top | km | 254 (base), 383 (top) |
| $\omega_\mathrm{ac}$ | Acoustic cutoff frequency $c/2H$: waves slower than this cannot propagate freely | $\mathrm{s^{-1}}$ | period $2\pi/\omega_\mathrm{ac} \approx 367$ s |
| $v$, $\psi$ | Velocity perturbation of a wave, and $v$ rescaled as $\psi = e^{-z/2H}\,v$ | $\mathrm{m\,s^{-1}}$ | |
| $k$ | Vertical wavenumber of the fundamental (quarter-wave) mode, $\pi/2D$ | $\mathrm{m^{-1}}$ | |
| $P_1$ | Period of the fundamental ringing mode | s | $\approx$ 190–200 |
| $\tau_d$ | E-folding time of the ringing amplitude | s | 130–170 |
| $A_0$, $\varepsilon$ | Initial relative size of the top/base flux mismatch ($\sim 1$), and the tolerance for calling the flux coherent (5 %) | – | |
| $t_w$ | Time after which $\lvert F_\mathrm{base}/F_\mathrm{top} - 1\rvert$ stays below $\varepsilon$: the end of the wave phase | s | 520–840 |
| $u^*$, $\Delta p$ | Face velocity from the Riemann solver at the base, and the pressure deficit at the base | $\mathrm{m\,s^{-1}}$, Pa | |
| $\tau_M$ | Column-mass response time, $\Delta M = F\tau_M$ | s | 45–50 |
| $\tau_\mathrm{TR}$ | E-folding time of the TR relaxation, from the fit to $\mathrm{ER}$ | s | 530–810 at 15–22 kK |
| $t_\mathrm{TR}$ | Time after which $\mathrm{ER}$ stays within 0.05 of its 4000 s value: the TR has relaxed | s | 2260 (15 kK) to $\lesssim 220$ ($\geq 45$ kK) |
| $t_\mathrm{ready}$ | $\max(t_w, t_\mathrm{TR})$: the earliest time at which the flux can be read | s | |

**Question.** In `visualization/model_column/twall_sweep_N500_4000s_overlay_p_abs.mp4` the mass flux $F$ evolves quickly at first, looks quasi-steady, and then decays slowly. Two guesses were on the table:
1. the early evolution is set by wave propagation through the domain;
2. the late decay is a numerical effect of the fixed upper-boundary pressure.

How long must a run last before $F$ is steady, and how long before it decays? The goal is to read a steady $F$ for each boundary $(T, p)$.

**Answer.**

1. **Guess 1 is right.** The first stage is an acoustic pulse bouncing between the two reservoir faces:
   - it takes $t_\mathrm{ac} \approx$ 50–60 s ($\int ds/c$) to reach the base;
   - it rings with a period of about 190–200 s and damps with an e-folding time of 130–170 s;
   - $F$ is height-coherent (base/top within 5 %) after **520–840 s**, about 10–14 $t_\mathrm{ac}$, at every wall temperature.
2. **Guess 2 is not right.** The boundary data are constant: the back-pressure is fixed, the top-cell pressure moves by at most 0.013 % in 4000 s, and $T_\mathrm{top}$ is thermostatted. The decay is physical within the release equations, and it is quantitatively explained:
   - the upflow lifts the chromosphere;
   - the release has no chromospheric energy source, so every lifted parcel cools **exactly adiabatically** (verified to $\leq 3$ % at fixed column mass);
   - the cooled chromosphere contracts, so the TR lengthens beneath the fixed-height wall;
   - $F$ falls as $1/L_\mathrm{TR}$.

   The fixed pressure only sets the rate, through the density at the TR base. Because $p$ cannot rise, the corona never fills up. That is what lets the evaporation go on indefinitely.
3. **There is a second relaxation, and it depends strongly on $T_\mathrm{wall}$.** The C7 initial TR is not a conduction–enthalpy TR. It takes about 2200 s to reach that balance at 15–18 kK, about 1600 s at 22 kK and about 520 s at 30 kK. At $\geq 45$ kK it gets there in $\lesssim 200$ s, inside the wave phase.
4. **There is no plateau.** The flux $F$, which is what the mass-flux panel draws, starts decaying as soon as it is ready:

   $$F(t)=\frac{J}{\sqrt{L_1^2+2K_\mathrm{eff}J\,(t-t_1)}},$$

   which drifts by 2–6 % per 500 s at 4000 s. What *is* steady is the invariant $J = F L_\mathrm{TR}$. The first-principles $J$ reproduces it to 0.95–0.98 at $t_\mathrm{ready}$ and 0.89–1.01 at 4000 s. From $t_\mathrm{ready}$ to 4000 s, $F$ falls by 10–57 % while $F L_\mathrm{TR}$ changes by only 4–12 %. **Report $J(T, p)$ (and $F$ with its $L_\mathrm{TR}$ and reading time), not $F$ from a "plateau".**

Two by-products correct existing documentation:
- The documented "~13 % monotone spread" of `f_total` over 1600–1850 km is the Eulerian signature of this contraction. It is predicted to 0.2 percentage points on the current release run (12.1 % measured, 11.9 % predicted), and it is not a duration or resolution question.
- The "softening" of the local exponent of $V_\mathrm{top}(T_\mathrm{wall})$ in the sweep recap is mostly an artefact of reading every wall at the same 3500–4000 s.

Evidence:
- script: `util/twall_relaxation_diag.py`
- metrics: `outputs/model_column/twall_sweep/twall_relaxation_metrics.json`
- figures: `visualization/model_column/twall_relaxation_{phases,spacetime,timescales}.png` (commands in `visualize_commands.md` §16)

Post-processing only; no new run.

---

## 0. Setup

Symbols are defined in the Notation tables at the top.

- **Runs.** The 13 accepted runs of the wall-temperature sweep (`upper_wall_temperature_sensitivity_recap.md`): N=500/R4, double state, CFL 0.50 (100 kK at 0.25), 4000 s, `termination=end_time`, zero Godunov fallbacks. They use the previous decoupled outer closure (`ISO_HYDRO_T_DECOUPLE=1`). The one test that was repeated on the current single-T release run (`outputs/model_column/release_singleT_N500_4000s`) reproduces (§5a). The mechanism does not involve the outer temperature closure.
- **Boundaries.** Both are fixed ghost reservoirs:
  - base (1600 km): the IC state with $u = 0$;
  - top (2153 km): the back-pressure `kOuterPRef` $= 0.0102856$ Pa and the extrapolated, Mach-capped velocity. The conduction wall $T_\mathrm{wall}$ sits at the physical face.
- **$F$** is the face flux `f_total`. The top-face value equals the top-cell $\rho u$ to 0.1 %.

The thermodynamics is the pure-hydrogen Saha closure of `eos.hpp`, re-implemented in the script. It reproduces the runs' $x$ and $p$ exactly and the tabulated $\Gamma_1$ to $2\times10^{-4}$.

**Height coordinate.** The `.gamma_diag` height column is each cell's *upper-face* height (it equals `face_km`), not the cell centre. The new script rebuilds the centres. Scripts that treat it as the centre place every quantity half a cell (0.14–0.55 km) too high.

## 1. Stage 1 — the wave phase ($0 \to$ ~600–800 s)

**Mechanism.** At $t = 0$ the wall and the C7 top cell disagree (22.55 kK against a 22 kK wall, and far more for other walls), so the TR starts evaporating immediately. The mass drawn out of the top launches a rarefaction down the column. Both faces hold fixed reservoir states, so the pulse reflects and the column rings until the base reservoir resupplies the mass and the flux is the same at every height.

`twall_relaxation_spacetime.png` shows $F(z,t)/F_\mathrm{top}(t)$: the front and its reflections follow the $\pm c$ characteristics drawn from the top at $t = 0$.

**One-way acoustic time.** From the initial column,

$$t_\mathrm{ac}=\int_{1600}^{2153}\frac{ds}{c}=58.8\ \mathrm{s}\ \ (c=c_\mathrm{eq}),\qquad 47.7\ \mathrm{s}\ \ \left(c=\sqrt{\tfrac53\,p/\rho}\right),$$

with $c_\mathrm{eq}$ rising from 8.7 km/s at the base to 10.8 km/s at the top of the chromosphere. The measured front reaches the base between 53 and 65 s, depending on wall and detection threshold, between the two bounds.

**Ringing period.** Linear vertical acoustic–gravity waves in a stratified column obey a Klein–Gordon equation. With $v = e^{z/2H}\psi$,

$$\psi_{tt}=c^2\psi_{zz}-\omega_\mathrm{ac}^2\psi,\qquad \omega_\mathrm{ac}=\frac{c}{2H}.$$

The initial column gives $H = 254 \to 383$ km (base $\to$ chromosphere top) and an acoustic-cutoff period $2\pi/\omega_\mathrm{ac} \approx 367$ s. The fundamental quarter-wave mode ($k = \pi/2D$, $D = 553$ km) has

$$P_1=\frac{2\pi}{\sqrt{\omega_\mathrm{ac}^2+(\pi c/2D)^2}}\approx 186\ \mathrm{s}\ \text{(WKB phase integral on the IC)},\quad 193\ \mathrm{s}\ \text{(median } c,\ \omega_\mathrm{ac}\text{)}.$$

Measured: **191–199 s** at 15–22 kK and 181 s at 45 kK, from the zero crossings of $F_\mathrm{base}/F_\mathrm{top} - 1$. The 30–100 kK values (134–262 s) are contaminated by the hot-wall numerical ripple. The quarter-wave assignment is inferred from the period. The reflection coefficients of the two ghost closures were not derived.

**Damping and settling.** The ringing envelope decays with an e-folding time $\tau_d \approx$ 130–170 s, less than one period. The residual non-uniformity therefore falls below $\varepsilon$ at

$$t_w\simeq t_\mathrm{ac}+\tau_d\ln\frac{A_0}{\varepsilon}\approx 59+150\ln 20\approx 510\ \mathrm{s}\qquad(A_0\sim1,\ \varepsilon=5\,\%).$$

Measured $\lvert F_\mathrm{base}/F_\mathrm{top} - 1\rvert < 5$ % after **520–840 s**, the same at every wall temperature, because the chromospheric sound speed does not depend on the wall.

**The column mass is a fast variable.** Once the flux is coherent, the base reservoir supplies $F$ through a small pressure deficit:
- the column carries a mass deficit $\Delta M = M(0) - M = F\tau_M$ with $\tau_M =$ 45–50 s measured ($\geq 22$ kK);
- the linear Riemann estimate at the base face ($u^* \approx \Delta p/2\rho c$ with $\Delta p = g\,\Delta M$) gives $\tau_M = 2c_0/g = 64$ s.

So the column mass tracks $F$ within about a minute, and the fixed-pressure reservoirs keep $M$ within $\approx 0.1$ % of $M(0)$ in every run.

## 2. Stage 2 — the TR relaxes to the conduction–enthalpy balance

**Steady TR theory.** The release energy equation has no radiation and no heating:

$$\partial_t\!\left(\rho e+\tfrac12\rho u^2+\rho\phi\right)+\partial_s\!\left[F\left(h+\tfrac12u^2+\phi\right)+q\right]=0,\qquad q=-\kappa\,\partial_sT.$$

For a steady, subsonic, thin TR (so $g\,\Delta z \lesssim 0.5$ % of $\Delta h$), integrate from the chromosphere top, where $q \approx 0$:

$$\kappa(T)\,\frac{dT}{ds}=F\,[h(T)-h_c]\quad\Longrightarrow\quad \boxed{F\,L_{T_1}=J\equiv\int_{T_1}^{T_\mathrm{top}}\frac{\kappa(T,p)}{h(T,p)-h_c}\,dT}\qquad(1)$$

Here $L_{T_1} = s_\mathrm{top} - s(T_1)$. The integral diverges logarithmically as $T_1 \to T_c$: the foot is an exponential tail of length $\kappa_c/(F c_p)$. A marker is therefore needed; $T_1 = 7.5$ kK here.

Three consequences:

- **The steady TR is a one-parameter family.** $F$ and $L$ are not separately fixed by TR physics; their product is. Without a radiative sink, nothing in the TR sets $L$ (§3).
- **The energy balance at the top is $F(h_\mathrm{top} - h_c) = q_\mathrm{wall}$.** The sweep recap's one-line estimate used $h_c = 0$ and was off by 1.2–2$\times$. With $h_c$, the ratio $F(h_\mathrm{top} - h_c)/q_\mathrm{wall}$ is 0.89–1.0 once relaxed.
- **$J$ is a function of the boundary state only.** For $T \gg \chi_H/5k_B$, $h \approx 5k_BT/m_H$ and $J \to (2\kappa_0 m_H/25k_B)\,T_\mathrm{top}^{5/2}$. The ionization term steepens it at lower $T$.

**Verification.** $J$ is computed from first principles: Spitzer plus neutral $\kappa$ (`util/conductive_flux.py`) and the Saha $h$ at the run's $p_\mathrm{top}$, $T_\mathrm{top}$ and measured $h_c$. It is compared with the measured $F L$:

| $T_\mathrm{wall}$ [kK] | $F L/J$ at 300 s | 1000 s | 2000 s | 4000 s |
|---:|---:|---:|---:|---:|
| 15 | 1.64 | 1.23 | 0.98 | 0.89 |
| 18 | 2.09 | 1.23 | 1.00 | 0.92 |
| 22 | 1.82 | 1.08 | 0.95 | 0.92 |
| 30 | 1.25 | 0.95 | 0.94 | 0.94 |
| 45 | 0.98 | 0.98 | 0.94 | 0.96 |
| 70 | 0.97 | 0.97 | 0.97 | 0.99 |
| 100 | 0.95 | 1.01 | 0.99 | 1.01 |

At $\geq 45$ kK the TR is in balance from $\sim 100$ s on. At 15–22 kK it starts far out of balance: the ratio of 1.6–2.1 means the TR is releasing stored energy, mostly hydrogen ionization energy of the C7 TR gas, which leaves as extra enthalpy flux. The late values of 0.89–0.92 are the energy cost of ionizing the chromospheric gas that the lengthening TR absorbs (§3). That cost is a fraction $\approx \delta(\nabla_\mathrm{ad} - \nabla)\,\Delta x\,(\chi_H/m_H)/\Delta h$, about 0.08 at 22 kK and 0.02 at 100 kK.

**Relaxation time.** Two measures:
- the e-folding $\tau_\mathrm{TR}$ of a fit $\mathrm{ER}(t) = \mathrm{ER}_\infty + A\,e^{-t/\tau_\mathrm{TR}}$ to the energy ratio;
- the time $t_\mathrm{TR}$ after which $\mathrm{ER}$ stays within 0.05 of its 4000 s value.

| $T_\mathrm{wall}$ [kK] | 15 | 18 | 22 | 30 | $\geq 45$ |
|---|---:|---:|---:|---:|---:|
| $\tau_\mathrm{TR}$ [s] | 807 | 709 | 534 | 167 | already relaxed at 300 s |
| $t_\mathrm{TR}$ [s] | 2260 | 2160 | 1560 | 520 | $\lesssim 220$ |

This stage is **not** derived a priori. The obvious estimates overestimate it by 2–12$\times$ and do not reproduce its weak dependence on $T_\mathrm{wall}$ between 15 and 22 kK:
- the TR-base flush time $L/V_f$;
- the TR mass over $F$;
- the local conduction time at the top.

Treat these values as measured properties of the C7 initial condition. A pre-relaxed IC would change this stage only.

## 3. Stage 3 — the secular decay: adiabatic contraction of the lifted chromosphere

**(a) Parcels are lifted and cool adiabatically.** Below the TR, $\kappa$ is negligible, and the release has no radiative loss or heating, so $DS/Dt = 0$. Label each parcel by the column mass above it, $m$. Mass leaves through the top at $F_\mathrm{top}$ and nothing enters from above, so every label falls by $\Phi(t)$: the parcel now at $m$ started at $m + \Phi$. Hydrostatic balance, $p = p_\mathrm{top} + g m$ with $p_\mathrm{top}$ fixed, means the parcel has decompressed from $p + g\Phi$ to $p$. Along its Saha adiabat,

$$T(m,t)=T_0(m+\Phi)\left[\frac{p}{p+g\Phi}\right]^{\nabla_\mathrm{ad}}.\qquad(2)$$

**Tested at fixed column mass at 4000 s, from 1620 km to 10 km below the TR base (`..._phases.png` (e)):**
- 22 kK: every chromospheric cell matches to $\pm$0.1–0.3 K out of up to $-17$ K; the $\leq 6$ % worst case is where $\Delta T$ is only a few K.
- 30–100 kK: the maximum relative error over the chromosphere is $\leq 3$ % (0.9, 0.8, 2.2 and 3.2 % at 30, 45, 70 and 100 kK).

The displacements are large: at 100 kK the gas at 2100 km at 4000 s started 228 km lower and has cooled by 600 K.

**(b) Cooling at fixed pressure contracts the column.** From $(\partial\ln\rho/\partial S)_p = -\delta/c_p$ and $dS/dm = -(g c_p/p)(\nabla_\mathrm{ad} - \nabla)$ along the IC,

$$\Delta\ln(1/\rho)\big|_m=-\delta(\nabla_\mathrm{ad}-\nabla)\,\frac{g\Phi}{p},\qquad
\Delta z(m)=-\Phi\,W(z),\quad W(z)\equiv\int_{z_b}^{z}\delta(\nabla_\mathrm{ad}-\nabla)\,\frac{g}{p}\,dz'.\qquad(3)$$

At the chromosphere top this gives the **contraction coefficient**

$$K=W(z_f)\approx\frac{\delta(\nabla_\mathrm{ad}-\nabla)}{\rho_f}\left(1-\frac{\rho_f}{\rho_b}\right).$$

Partially ionized hydrogen has $\nabla_\mathrm{ad} \approx 0.05$ but $\delta \approx 4$, so $\delta(\nabla_\mathrm{ad} - \nabla) \approx$ 0.20–0.23 throughout the chromosphere. The ionization amplifies the contraction: recombination adds to the thermal shrinkage. On the initial column $K = 1.61\times10^{9}\ \mathrm{m^3\,kg^{-1}}$. Evaluated on the current state, $K_\mathrm{cur}$ carries the finite-displacement correction, falling to $0.90\times10^{9}$ at 100 kK by 4000 s, where $g\Phi/p \approx 0.8$.

**(c) The TR lengthens under the fixed-height wall.** The TR top is pinned at 2153 km, while its base follows the chromosphere top down. The TR must also absorb chromospheric mass as it grows. With $M_\mathrm{TR} = \mu\rho_f L$,

$$dL=K\,d\Phi+\mu\,dL\quad\Longrightarrow\quad \frac{dL}{d\Phi}=K_\mathrm{eff}=\frac{K}{1-\mu}.\qquad(4)$$

| $T_\mathrm{wall}$ [kK] | 22 | 30 | 45 | 70 | 100 |
|---|---:|---:|---:|---:|---:|
| $\mu$ (4000 s) | 0.42 | 0.30 | 0.19 | 0.11 | 0.08 |
| $K_\mathrm{cur}/(1-\mu)$ [$10^{9}\ \mathrm{m^3\,kg^{-1}}$] | 2.95 | 2.11 | 1.63 | 1.26 | 0.98 |
| measured $dL/d\Phi$ ($t \geq t_\mathrm{TR}$) [$10^{9}\ \mathrm{m^3\,kg^{-1}}$] | 2.63 | 2.13 | 1.71 | 1.37 | 1.11 |

The model agrees to $\leq 12$ % at $\geq 22$ kK. At 15–18 kK, $\mu$ exceeds 0.5 and the self-similar TR assumption behind (4) overpredicts by 1.4–2$\times$. Use the measured $dL/d\Phi \approx 2.9\times10^{9}\ \mathrm{m^3\,kg^{-1}}$ there.

**(d) $F$ decays.** Combine (1) and (4): $dL/dt = K_\mathrm{eff}F = K_\mathrm{eff}J/L$, so

$$L^2=L_1^2+2K_\mathrm{eff}J\,(t-t_1),\qquad F=\frac{J}{\sqrt{L_1^2+2K_\mathrm{eff}J\,(t-t_1)}}\ \propto t^{-1/2}\ \text{late},\qquad(5)$$

$$r\equiv-\frac{d\ln F}{dt}=\frac{K_\mathrm{eff}F}{L}=\frac{K_\mathrm{eff}F^2}{J}\approx\frac{\delta(\nabla_\mathrm{ad}-\nabla)}{1-\mu}\,\frac{V_f}{L}.\qquad(6)$$

The decay e-folding time is $\approx$ 3–5 **TR-base flow crossings $L/V_f$**. Measured at 4000 s, the $L$-growth part of $r$ against (6):

| $T_\mathrm{wall}$ [kK] | 22 | 30 | 45 | 70 | 100 |
|---|---:|---:|---:|---:|---:|
| $d\ln L/dt$ measured [$10^{-5}\ \mathrm{s^{-1}}$] | 4.65 | 6.96 | 9.22 | 10.5 | 10.2 |
| (6) predicted [$10^{-5}\ \mathrm{s^{-1}}$] | 5.36 | 6.95 | 9.00 | 9.95 | 9.92 |
| total $F$ drift $r$ [$10^{-5}\ \mathrm{s^{-1}}$] | 6.15 | 7.22 | 11.6 | 10.4 | 11.4 |
| $F$ drift per 500 s | 3.1 % | 3.6 % | 5.8 % | 5.2 % | 5.7 % |

The difference between $r$ and the $L$-growth rate is the residual drift of $F L$: 0.8–1.2 % per 500 s at 15–22 kK, where the TR is still settling, and 0.1–1.2 % per 500 s, dominated by noise, at $\geq 30$ kK.

At 100 kK, (5) with $L_1 = 14$ km at 100 s and $K_\mathrm{eff} \approx 1.2\times10^{9}\ \mathrm{m^3\,kg^{-1}}$ gives a crossover time $L_1^2/(2K_\mathrm{eff}J) \approx 330$ s. The flux is on its $t^{-1/2}$ branch before the wave phase is even over. The prediction is $L = 27$ km at 1000 s and 50 km at 4000 s, against 28.3 and 48.7 km measured, and $F(1000)/F(4000) = 1.85$ against 1.79 measured. **On the movie's symlog axis a $t^{-1/2}$ decline looks like a plateau.**

## 4. Why the fixed upper pressure is not the cause, and what it does do

- **The boundary does not change.** `kOuterPRef` is a constant, and the top-cell pressure drifts by $\leq 0.013$ % over 4000 s (100 kK: $1.02877\times10^{-2} \to 1.02864\times10^{-2}$ Pa). $T_\mathrm{top}$ stays within 0.1–0.3 % of the wall. Nothing at the boundary evolves on the decay timescale.
- **The decay is fully accounted for without it.** At $\geq 22$ kK, §3 closes the decay rate to $\leq 15$ %, the TR lengthening to $\leq 12$ %, the temperature change to $\leq 6$ % ($\leq 3$ % at $\geq 30$ kK) and the flux profile to $\leq 1.5$ % (§5a). It is also not numerical: it is a resolved process, with displacements of 3–230 km against 0.28–1.1 km cells, obeying the continuum adiabat (2). It is a physical consequence of a release energy equation with no chromospheric source. In the Sun, radiative losses and heating hold the chromospheric temperature. `main.tex` already lists the missing radiative sink as the reason quasi-steady evaporation is not claimed.
- **What the fixed pressure does:**
  - It sets $\rho_f \approx p_f/(R\,T_c)$, the denominator of $K$, so the decay time scales with $p_\mathrm{top}$.
  - Together with the fixed-height wall, it pins the TR top while the chromosphere contracts beneath it.
  - Because $p_\mathrm{top}$ cannot rise, there is **no coronal-filling feedback**, which is what normally chokes off evaporation in a closed loop. The evaporation therefore never self-limits and keeps lifting the chromosphere.
- **Loosening the pressure is not the fix.** Without it the model has no steady state at all: the hydrostatic Neumann closure drifts 5.8 % in $p_\mathrm{top}$ and cuts the flux by about $19\times$ (`outer_bc_pressure_temperature_swap_recap.md`).

## 5. Two corrections to existing documentation

**(a) The ~12–13 % `f_total` spread over 1600–1850 km is the contraction's Eulerian signature.** Mass conservation with (3) gives

$$\frac{F(z)}{F_\mathrm{top}}=1-\rho(z)\,W(z),\qquad\text{equivalently}\qquad F(z)-F_\mathrm{top}=\frac1g\,\frac{\partial p}{\partial t}\Big|_z.\qquad(7)$$

As the column contracts, the pressure at fixed height falls, and the flux through a fixed height is reduced by $\rho$ times the downward speed of the contracting levels. The TR returns that deficit, since its base sweeps through dense gas at $dL/dt$, so the top and base fluxes still agree (to 0.3 % at 22 kK).

Measured over 3800–4000 s against (7) evaluated on the 4000 s state:
- 1600–1850 km spread: **0.119–0.121 measured, 0.114–0.120 predicted** at 15–45 kK; 0.093/0.092 at 100 kK.
- **Current single-T release run: 0.121/0.119**, with a maximum error over the chromosphere of 0.005.

The spread does not depend on resolution, because (7) is a continuum statement. It is essentially the same at every $T_\mathrm{wall}$ (size $\approx \delta(\nabla_\mathrm{ad} - \nabla)(1 - \rho/\rho_b)$). It persists as long as the column evaporates; its slow 14.3 → 13.0 % change reflects the slowly changing $\delta(\nabla_\mathrm{ad} - \nabla)$. It is therefore **not** "a duration/resolution question", as `validation.md` and `state_precision_release_cutover.md` had it. The localized first-interior-face artifact is separate and still open.

**(b) The softening of the $V_\mathrm{top}(T_\mathrm{wall})$ exponent is mostly a reading-time effect.** $F = J\,(F L/J)/L$, and every wall was read at 3500–4000 s, by which time hotter walls had contracted more. Local exponents between neighbouring walls, at 4000 s:

| interval [kK] | 15→18 | 18→22 | 22→30 | 30→45 | 45→70 | 70→100 |
|---|---:|---:|---:|---:|---:|---:|
| $J = F L$ (steady invariant) | 3.17 | 2.98 | 2.74 | 2.65 | 2.60 | 2.57 |
| $F$ at 4000 s | 3.23 | 2.72 | 2.37 | 2.00 | 1.74 | 1.66 |
| $L_\mathrm{TR}$ at 4000 s | 0.13 | 0.25 | 0.44 | 0.73 | 0.93 | 0.96 |

The exponent of $J$ tends to the Spitzer $5/2$. With $V_\mathrm{top} = F/\rho_\mathrm{top}$ and $\rho_\mathrm{top} \propto 1/T_\mathrm{wall}$, the fixed-$L$ exponent of $V_\mathrm{top}$ is **$\approx$ 3.6–3.7 across 22–100 kK**. The recap's fall from 3.73 to 2.63 is almost entirely the growth of $L$ with $T_\mathrm{wall}$ at a fixed reading time.

## 6. How long to run, and what to report

| $T_\mathrm{wall}$ [kK] | wave phase over (5 %) | TR relaxed ($t_\mathrm{TR}$) | **$t_\mathrm{ready}$** | $F$ change, $t_\mathrm{ready} \to$ 4000 s | $F L$ change | $F L/J$ at $t_\mathrm{ready}$ |
|---:|---:|---:|---:|---:|---:|---:|
| 15 | 835 s | 2260 s | **~2300 s** | −10 % | −7 % | 0.95 |
| 18 | 619 s | 2160 s | **~2200 s** | −13 % | −7 % | 0.98 |
| 22 | 580 s | 1560 s | **~1600 s** | −21 % | −9 % | 0.98 |
| 30 | 710 s | 520 s | **~700 s** | −37 % | −12 % | 0.96 |
| 45 | 535 s | $\lesssim 220$ s | **~540 s** | −47 % | −6 % | 0.96 |
| 70 | 522 s | $\lesssim 220$ s | **~520 s** | −55 % | −6 % | 0.96 |
| 100 | 520 s | $\lesssim 220$ s | **~520 s** | −57 % | −5 % | 0.96 |

Recipe for the current C7 initial condition:

1. **Run to $t_\mathrm{ready}$ plus an averaging window**, e.g. 500 s: about 1100 s for $\geq 45$ kK, 1200 s at 30 kK, 2100 s at 22 kK and 2800 s at 15–18 kK. Near the flow reversal ($\leq 14$ kK) $F \to 0$ and relative criteria mean nothing; use absolute ones.
2. **Report $J = F L_\mathrm{TR}$** as the steady quantity. Its first-principles value (1) needs no run at all, only $T_\mathrm{top}$, $p_\mathrm{top}$ and $h_c$. The simulated $F L_\mathrm{TR}$ is 0.95–0.98 of it at $t_\mathrm{ready}$ and 0.89–1.01 at 4000 s.
3. If an absolute $F$ is needed, report it **with its $L_\mathrm{TR}$ and reading time**, and compare runs at equal $L_\mathrm{TR}$ rather than equal time. In this release $L_\mathrm{TR}$ is set by the domain geometry and the drifting chromosphere, not by TR physics. Making $L$ physically determined needs a chromospheric/TR energy sink, which is outside release physics and a design decision.

**Planning a $p_\mathrm{top}$ sweep.** These are derived, not measured.

- There is **no knob** for it today. `kOuterPRef` is captured from the IC top cell; one would have to be added.
- **Unchanged by $p_\mathrm{top}$:** $t_\mathrm{ac}$ and the wave phase, which are set by the chromospheric sound speed.
- **Weakly changed:** $J$, through the ionization terms of $h(T, p)$ and $\kappa(n_e, n_\mathrm{HI})$.
- **Strongly changed:** $L_\mathrm{TR}$, which is geometric. The chromosphere ends where its hydrostatic pressure reaches $p_\mathrm{top}$, so $\Delta L \approx H_f\,\Delta\ln p_\mathrm{top}$ with $H_f \approx 380$ km. A 10 % change in $p_\mathrm{top}$ moves the TR by about 36 km, against $L_\mathrm{TR} \approx$ 13–15 km at 22 kK. That changes $F$ by a factor of order 3.
  - Lowering $p_\mathrm{top}$ below the IC value puts the would-be chromosphere top above 2153 km and squeezes the TR into the top cells.
  - **At fixed $T_\mathrm{wall}$, an $F(p)$ sweep on this IC would measure mostly this geometric effect.** Rebuild the IC consistently with each $p_\mathrm{top}$, or report $J$.
- **Both slow timescales** (TR relaxation and decay) scale like the TR-base flow crossing $L/V_f = \rho_f L^2/J$, i.e. $\propto p_\mathrm{top} L^2/J$.

## 7. Limits

- Only N=500/R4 was analysed. The standing caveat applies: long-duration evaporation mass flux at N500 is not established as $< 10$ % grid-converged. Relative timescales are better determined than absolute $F$.
- The sweep used the decoupled outer closure. Only the flux-profile test (§5a) was repeated on the single-T release run. The outer-$T$ closure moves $V_\mathrm{top}$ by 0.8 % and does not enter the mechanism.
- $\tau_\mathrm{TR}$ (§2) is measured, not derived.
- The 70 and 100 kK face fluxes carry a numerical ripple, $\pm 30$ % spatially in the lower chromosphere at 100 kK. Profiles are therefore time-averaged over 200 s, and the hot-wall ringing periods are not quoted.
- Walls $\leq 14$ kK, near the flow reversal, are excluded from every ratio-based metric. The mechanism still applies with the signs reversed: downflow compresses and warms the chromosphere, and the TR shortens.

## Reproduce

```bash
.venv/bin/python util/twall_relaxation_diag.py \
    --run 100000=outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt
```

This needs the sweep dumps with their `.faceflux` sidecars (`visualize_commands.md` §16 prerequisites) and the release run `outputs/model_column/release_singleT_N500_4000s`. Runtime is about 10 s on 6 workers.
