# Steady transition-region length $L_\mathrm{TR}(T)$ with the release closures

**Status: CURRENT**

**Question.** In a steady state ($\partial_t = 0$), what is the length of the transition region as a function of temperature? The starting point is the field-aligned energy equation with radiative loss $Q(T)$. The model is the release closure set: pure-hydrogen Saha EOS with the ionization energy included, plus the Spitzer-plus-neutral conductivity. Gravity is dropped. The kinetic term $u^2/2$ is kept, even though it is small.

**Notation.** Only the fluid variables are used: $\rho$, $u$, $p$, $T$, the conductive flux $q=-\kappa\,\partial_sT$ and the loss $Q$. The three shorthands are $F=\rho u$, $h=(e_\mathrm{int}+p)/\rho$ and $H=h+\tfrac12u^2$. In the TR $q<0$ (heat flows down), so $-q=|q|$ is the positive downward flux.

**Answer.**

1. The energy equation turns into a first-order ODE in $T$ for the conductive flux $q$ ((8) below):
   $$\frac{dq}{dT}=-F\,\frac{dH}{dT}+\frac{\kappa\,Q}{q},\qquad q(T_c)=0,\qquad H=h+\tfrac12u^2 .$$
   The length then follows from $ds/dT=-\kappa/q$ by quadrature:
   $$L_\mathrm{TR}(T)=-\int_{T_1}^{T}\frac{\kappa(T')}{q(T')}\,dT'.$$
2. **Release ($Q = 0$): exact closed form.**
   $$L_\mathrm{TR}(T)=\frac{1}{F}\int_{T_1}^{T}\frac{\kappa\,dT'}{h(T',p)-h_c+\tfrac12\left(u^2-u_c^2\right)},\qquad u=\frac{F}{\rho(T',p)} .$$
   This is the existing invariant $F L_\mathrm{TR} = J$ (`twall_mass_flux_relaxation_timescales.md`, eq. 1), with one change: $J$ now depends weakly on $F$ through $u^2$ and through the dynamic-pressure drop.
3. **Size of the $u^2$ term.** In the fully ionized limit
   $$F\,L_\mathrm{TR}=J_0\left(1-\frac17\,\frac{\rho_\mathrm{top}u_\mathrm{top}^2}{p_\mathrm{top}}+\dots\right)=J_0\left(1-\tfrac{5}{21}M_\mathrm{top}^2+O(M^4)\right).$$
   On the release column ($T_\mathrm{top} = 22$ kK, $p_\mathrm{top} = 0.0103$ Pa, $F = 2.68\times10^{-10}$ kg m⁻² s⁻¹), the correction to $L$ is $-7.9\times10^{-8}$.
4. **With radiation ($Q > 0$): bracketed in closed form.** Let $q_e = -F\,\Delta H$ (pure enthalpy) and $q_s = -\left(2\int\kappa Q\,dT\right)^{1/2}$ (pure static radiative). Then
   $$\sqrt{q_e^2+q_s^2}\;\le\;|q|\;\le\;|q_e|+|q_s| .$$
   So $L_\mathrm{TR}$ is pinned to within a factor $\le\sqrt2$ without solving the ODE. A loss $Q$ that switches on at $T_c$ also removes the logarithmic divergence at the TR foot, so $L_\mathrm{TR}$ measured from $T_c$ itself becomes finite.

5. **Radiation-dominated TR ($Q>0$ in the TR, $Q=0$ in the chromosphere): $F$-independent to first order.** At $F=0$ the length $L_s(T;p,T_c)=\int_{T_c}^{T}\kappa\,dT'/|q_s|$ depends only on $T$, $p$, $T_c$ and $\Lambda$. For $F>0$ the correction is linear in $F/F_*$ with $F_*=|q_s|/\Delta H$; for a fully ionized $\Lambda\propto T^{-1/2}$ it is $L_\mathrm{TR}\simeq L_s\left(1-\tfrac47F/F_*\right)$. $F$ is then fixed by the heat flux at the top, not by the TR length (§6e).

Radiation is **not** release physics: the release energy equation has no $Q$. Item 2 is the release result. Items 1, 4 and 5 are the general result, which reduces to item 2 at $Q = 0$.

Numerical check: `util/release/tr_length_steady.py` (§8). Post-processing only; no simulation.

---

## 1. Model and assumptions

Release settings used (`src/single_fluid/mixture.hpp`, `eos.hpp`, `src/single_fluid/integrator.cpp`):

| Item | Release choice |
|---|---|
| Fluid | one fluid, one velocity $u$, one temperature $T$ |
| Geometry | constant cross-section, so the $p\,\partial_s\ln A$ source vanishes |
| EOS | pure hydrogen, Saha equilibrium, $e_\mathrm{int}=\tfrac32p+x\,n_H\chi_H$, $p=(1+x)n_Hk_BT$, $n_H=\rho/m_H$ |
| Conduction | $q=-\kappa\,\partial_sT$, with $\kappa=\kappa_e+\kappa_n$ (Spitzer with e–H collisions, plus neutral H); no TRAC, no artificial diffusivity |
| Release top state | Dirichlet $T_\mathrm{top}=22$ kK, back-pressure $p_\mathrm{top}=k_\mathrm{OuterPRef}\approx0.0103$ Pa |
| Chromosphere top | $T_c\approx6.7$ kK, where $q\approx0$ |

Assumptions used here:
- steady state, $\partial_t=0$;
- gravity dropped ($\phi=0$);
- $u^2$ kept;
- $Q(T)\ge0$, with $Q=0$ for $T\le T_c$ (chromosphere) and $Q>0$ allowed for $T>T_c$ (TR);
- $s$ increases upward, the TR has $dT/ds>0$, so $q<0$ (downward), and evaporation has $F>0$.

**Saha in steady state.** In the release, frozen versus equilibrium composition matters only for the acoustic response inside the Riemann solver. A steady state has no waves, and the EOS decode always imposes Saha. Assuming Saha equilibrium everywhere is therefore exactly the release's steady state.

## 2. Steady conservation laws

With $\partial_t=0$ and constant area, the three release equations become exact ODEs in $s$.

**Mass:**
$$\frac{d}{ds}(\rho u)=0\quad\Longrightarrow\quad \rho u=F=\text{const}.\tag{1}$$

**Momentum** (no gravity):
$$\frac{d}{ds}\left(p+\rho u^2\right)=0\quad\Longrightarrow\quad p+\rho u^2=p_\mathrm{top}+\rho_\mathrm{top}u_\mathrm{top}^2 .\tag{2}$$

**Energy.** With $E=e_\mathrm{int}+\tfrac12\rho u^2$, the energy flux is $(E+p)u=\rho u\left(h+\tfrac12u^2\right)$, where $h=(e_\mathrm{int}+p)/\rho$. Adding the radiative sink:
$$\frac{d}{ds}\Big[F\,H+q\Big]=-Q,\qquad H\equiv h+\tfrac12u^2 .\tag{3}$$

Equation (2) is what "keep $u^2$" requires for consistency. The kinetic term in (3) and the dynamic pressure $\rho u^2$ in (2) are the same order, $O(M^2)$. The pressure drop matters because, in the partially ionized TR, $h$ and $\kappa$ depend on $p$ through $x(T,p)$.

## 3. First integral

Integrate (3) from the chromosphere top $s_c$ (where $T=T_c$ and $q=q_c$) to $s$:
$$F\,H(s)+q(s)=F\,H_c+q_c-\int_{s_c}^{s}Q\,ds' .\tag{4}$$

**Boundary condition at the chromosphere top.** $Q=0$ in the chromosphere and the chromosphere is nearly isothermal at $T_c$, so $q_c=-\kappa\,dT/ds|_c=0$. This is the same condition as the existing $J$ table: the release diagnostic reads $h_c$ where $|q|<1\%\,|q_\mathrm{top}|$. Then
$$q(s)=-F\,\big[H(s)-H_c\big]-\int_{s_c}^{s}Q\,ds'.\tag{5}$$

The downward conductive flux $|q|$ at any height supplies two sinks: the enthalpy (plus kinetic energy) carried up by evaporation, and everything radiated below that height.

## 4. Change of variable to $T$

From $q=-\kappa\,dT/ds$:
$$\frac{ds}{dT}=-\frac{\kappa}{q}.\tag{6}$$

**The loss term in $T$.** In (5) the enthalpy term $F(H-H_c)$ is already a function of $T$, but the column loss $\int Q\,ds$ is an integral over height. To get an ODE in $T$ alone, that term must be rewritten in $T$ too. This is allowed because $T(s)$ is strictly monotonic in the TR: $q<0$ there, so $dT/ds=-q/\kappa>0$, and $s(T)$ is single-valued. Substituting $ds=(ds/dT)\,dT$ with (6):
$$\int_{s_c}^{s}Q\,ds'=\int_{T_c}^{T}Q\,\frac{ds}{dT'}\,dT'=-\int_{T_c}^{T}\frac{\kappa\,Q}{q}\,dT',
\qquad
\frac{d}{dT}\int_{s_c}^{s}Q\,ds'=-\frac{\kappa\,Q}{q}.\tag{7}$$

Physically, the radiation lost over a temperature step $dT$ is $Q$ times the thickness of the layer that spans $dT$, and that thickness is $\kappa/|q|$ per kelvin.

Two remarks:
- At the base $q(T_c)=0$, so $ds/dT$ diverges there. The integrand $\kappa Q/q$ stays finite: $Q=0$ for $T\le T_c$, and just above $T_c$ with $Q>0$, $|q|$ grows like $(T-T_c)^{1/2}$ (§6c), so the integrand $\propto(T-T_c)^{-1/2}$ is integrable.
- If $T(s)$ were not monotonic (a temperature inversion), $s(T)$ would be multivalued and $s$ would have to stay the coordinate. The steady TR here has $q<0$ throughout, so it is monotonic.

Differentiate (5) with respect to $T$ and use (7):
$$\boxed{\;\frac{dq}{dT}=-F\,\frac{dH}{dT}+\frac{\kappa(T)\,Q(T)}{q},\qquad q(T_c)=0\;}\tag{8}$$

Equivalently, multiplying by $q$:
$$\frac{d}{dT}\left(\tfrac12q^2\right)=-F\,q\,\frac{dH}{dT}+\kappa\,Q .\tag{9}$$

Once $q(T)$ is known, the TR length to temperature $T$, measured from a base isotherm $T_1\ge T_c$, follows from (6):
$$\boxed{\;L_\mathrm{TR}(T)\equiv s(T)-s(T_1)=-\int_{T_1}^{T}\frac{\kappa(T')}{q(T')}\,dT'\;}\tag{10}$$

$L_\mathrm{TR}$ quoted without an argument means $L_\mathrm{TR}(T_\mathrm{top})$.

**Interpretation of $dH/dT$.** It is the total derivative along the solution path, so $p=p(T)$ from (2). The momentum equation gives $u\,du=-dp/\rho$, and $dh=d(e_\mathrm{int}/\rho)+p\,d(1/\rho)+dp/\rho$. Hence
$$dH=dh+u\,du=d\!\left(\frac{e_\mathrm{int}}{\rho}\right)+p\,d\!\left(\frac1\rho\right).\tag{11}$$

With $F\,d(1/\rho)/ds=du/ds$ and $F\,d(e_\mathrm{int}/\rho)/ds=d(e_\mathrm{int}u)/ds$, (3) becomes the steady internal-energy equation
$$\frac{d}{ds}\left(e_\mathrm{int}\,u\right)+p\,\frac{du}{ds}=-\frac{dq}{ds}-Q ,$$
and (8) is this equation written with $T$ as the coordinate. This is a check on (2)–(3), and it shows why the kinetic term cannot be added to (3) without (2): only the combination removes $u$ from the internal-energy balance.

## 5. Release closures along the path

### Ionization fraction

With
$$S(T)=\left(\frac{m_ek_BT}{2\pi\hbar^2}\right)^{3/2}e^{-\chi_H/k_BT},$$
pure-hydrogen Saha ($x^2/(1-x)=S/n_H$) combined with $p=(1+x)n_Hk_BT$ gives
$$\frac{x^2}{1-x^2}=\frac{S(T)\,k_BT}{p}\equiv A(T,p),\qquad x=\left(1+A^{-1}\right)^{-1/2}.\tag{12}$$

Its partial derivatives:
$$\left(\frac{\partial x}{\partial T}\right)_p=\frac{x(1-x^2)}{2T}\left(\frac52+\frac{\chi_H}{k_BT}\right),\qquad
\left(\frac{\partial x}{\partial p}\right)_T=-\frac{x(1-x^2)}{2p}.\tag{13}$$

### Density, enthalpy, heat capacity

Density:
$$\rho=\frac{m_H\,p}{(1+x)\,k_BT}.\tag{14}$$

Enthalpy:
$$h=\frac{\tfrac52(1+x)\,k_BT+x\,\chi_H}{m_H}.\tag{15}$$

Heat capacity at constant pressure:
$$c_p\equiv\left(\frac{\partial h}{\partial T}\right)_p=\frac{k_B}{m_H}\left[\frac52(1+x)+\frac{x(1-x^2)}{2}\left(\frac52+\frac{\chi_H}{k_BT}\right)^{2}\right].\tag{16}$$

Pressure derivative of the enthalpy:
$$h_p\equiv\left(\frac{\partial h}{\partial p}\right)_T=-\frac{x(1-x^2)}{2\,m_H\,p}\left(\frac52k_BT+\chi_H\right)\le0.\tag{17}$$

### Conductivity

The conductivity is evaluated with $n_e=x\,n_H$, $n_{HI}=(1-x)\,n_H$ and $n_H=p/[(1+x)k_BT]$, in SI units [W m⁻¹ K⁻¹]:
$$\kappa_e=\frac{9.2048\times10^{-12}\,n_e\,T^{5/2}}{n_e+2.836\times10^{-11}\,n_{HI}\,T^2},\qquad
\kappa_n=\frac{0.0342006\,n_{HI}\,T}{1.20613\,n_e\sqrt{2T}+1.70573\,n_{HI}\sqrt{T}} .\tag{18}$$

Both terms depend on the densities only through the ratio $n_{HI}/n_e=(1-x)/x$. So $\kappa=\kappa(T,x(T,p))$, and it depends on $p$ only through Saha. For $x\to1$, $\kappa\to\kappa_0T^{5/2}$ with $\kappa_0=9.2048\times10^{-12}$ W m⁻¹ K⁻⁷ᐟ².

### Pressure and $dH/dT$ along the path

From (2) with $u=F/\rho$:
$$p(T)=p_\mathrm{top}+F^2\left[\frac{1}{\rho_\mathrm{top}}-\frac{1}{\rho\big(T,p(T)\big)}\right].\tag{19}$$

Differentiating (2) gives $dp/dT=u^2\,d\rho/dT$, so
$$\frac{d\rho}{dT}=\frac{(\partial\rho/\partial T)_p}{1-u^2/c_T^2},\qquad c_T^2\equiv\left(\frac{\partial p}{\partial\rho}\right)_T=\frac{p/\rho}{1+\tfrac12x(1-x)} .\tag{20}$$

The logarithmic derivatives of $\rho$ follow from (13) and (14):
$$\left(\frac{\partial\ln\rho}{\partial\ln T}\right)_p=-1-\frac{x(1-x)}{2}\left(\frac52+\frac{\chi_H}{k_BT}\right),\qquad
\left(\frac{\partial\ln\rho}{\partial\ln p}\right)_T=1+\frac{x(1-x)}{2}.$$

The subsonic branch requires $u<c_T$, the Saha-isothermal sound speed. This is the usual sonic singularity of a steady conductive flow.

Then, with $dh=c_p\,dT+h_p\,dp$ and $u\,du=-dp/\rho$:
$$\frac{dH}{dT}=c_p-u^2\left(\frac1\rho-h_p\right)\frac{d\rho}{dT}.\tag{21}$$

The term $-(u^2/\rho)\,d\rho/dT$ is the kinetic energy gain. The term $u^2h_p\,d\rho/dT$ is the change of $h$ caused by the dynamic pressure drop shifting $x$; it is non-zero only where $0<x<1$. Equation (21) is (11) written explicitly. On the subsonic branch $d\rho/dT<0$ and $1/\rho-h_p>0$, so $dH/dT\ge c_p>0$.

## 6. Solutions

### 6a. Release: $Q=0$ (exact)

The structure is: (10) needs $q(T)$; with $Q=0$ the energy equation gives $q(T)$ algebraically; the closures of §5 turn every quantity in it into a function of $T$ alone. Four steps.

**Step 1: the flux $q(T)$.** With $Q=0$ the first integral (5) has no loss term:
$$q=-F\,\big[H-H_c\big].$$
The same result follows from the ODE (8): $dq/dT=-F\,dH/dT$, integrated from $T_c$ with $q(T_c)=0$. No ODE has to be solved, because without $Q$ the right side of (8) does not depend on $q$. Writing out $H=h+\tfrac12u^2$:
$$q(T)=-F\,\Delta H(T),\qquad \Delta H\equiv H-H_c=h(T,p)-h_c+\tfrac12\left(u^2-u_c^2\right).\tag{22}$$

**Step 2: insert $q$ into (10).** With (22), $-\kappa/q=\kappa/(F\,\Delta H)$. $F$ is constant by (1), so it comes out of the integral:
$$L_\mathrm{TR}(T)=-\int_{T_1}^{T}\frac{\kappa}{q}\,dT'=\int_{T_1}^{T}\frac{\kappa}{F\,\Delta H}\,dT'=\frac1F\int_{T_1}^{T}\frac{\kappa}{\Delta H}\,dT'.$$

**Step 3: make the integrand a function of $T'$ only.** $\kappa$, $h$ and $u$ depend on $(T,p)$, so the integrand is a function of $T'$ only once the pressure along the path is known. That comes from momentum (19), and the closures give the rest. At each $T'$:

| Quantity | From | Uses |
|---|---|---|
| $p(T')$ | momentum (19), $p=p_\mathrm{top}+F^2\left[\rho_\mathrm{top}^{-1}-\rho(T',p)^{-1}\right]$ | $\rho$ below |
| $x(T',p)$ | Saha (12) | $T'$, $p$ |
| $\rho(T',p)$ | EOS (14) | $x$ |
| $u=F/\rho$ | mass (1) | $\rho$ |
| $h(T',p)$ | enthalpy (15) | $x$ |
| $\kappa(T',p)$ | conductivity (18) | $x$ |

$p$ appears on both sides of (19), through $\rho$ and $x$. Iterating $p\leftarrow p_\mathrm{top}+F^2[\rho_\mathrm{top}^{-1}-\rho(T',p)^{-1}]$ from $p=p_\mathrm{top}$ converges: the derivative of the right side with respect to $p$ is $u^2/c_T^2<1$ on the subsonic branch (20). The same chain at $T_c$ gives the base values $p_c=p(T_c)$, $h_c=h(T_c,p_c)$ and $u_c=F/\rho(T_c,p_c)$. The top values $\rho_\mathrm{top}=\rho(T_\mathrm{top},p_\mathrm{top})$ come from (12) and (14).

**Step 4: the result.** Substituting Step 3 into Step 2:
$$\boxed{\;L_\mathrm{TR}(T)=\frac{1}{F}\int_{T_1}^{T}\frac{\kappa\big(T',p(T')\big)}{h\big(T',p(T')\big)-h_c+\tfrac12\left[u(T')^2-u_c^2\right]}\,dT'\equiv\frac{J(T;F)}{F}\;}\tag{23}$$

$J$ is defined as the integral, so (23) is $F\,L_\mathrm{TR}=J$. $F$ also appears inside $J$, through $u=F/\rho$ and through $p(T')$ in (19). That is the only way the kinetic term enters.

**Check: drop $u^2$.** Set $u\to0$ in (22) and $F^2\to0$ in (19). Then $p=p_\mathrm{top}$ everywhere and the denominator is $\Delta h_0=h(T',p_\mathrm{top})-h(T_c,p_\mathrm{top})$. $J$ becomes
$$J_0(T)=\int_{T_1}^{T}\frac{\kappa(T',p_\mathrm{top})}{h(T',p_\mathrm{top})-h(T_c,p_\mathrm{top})}\,dT',$$
independent of $F$. This is exactly the existing release invariant $F\,L_\mathrm{TR}=J(T_\mathrm{top},p_\mathrm{top},T_c)$ of `twall_mass_flux_relaxation_timescales.md` (1) and `analytic_evaporation_J_table.md`. Keeping $u^2$ changes one thing: $J$ is no longer independent of $F$ (§6b).

**The foot.** Near $T_c$, $\Delta H\simeq c_{p,c}(T-T_c)$, so the integrand behaves like $\kappa_c/[F\,c_{p,c}(T-T_c)]$ and $L_\mathrm{TR}$ diverges logarithmically as $T_1\to T_c$. The radiation-free TR foot is an exponential tail of e-folding length $\kappa_c/(F\,c_{p,c})$. For this reason a marker $T_1>T_c$ is required in the release; the existing tables use $T_1=7.5$ kK.

### 6b. Small-Mach expansion of the $u^2$ correction

Expand (23) in $F^2$ at fixed $p_\mathrm{top}$. To leading order (19) gives $\delta p\equiv p-p_\mathrm{top}=F^2\left(\rho_\mathrm{top}^{-1}-\rho^{-1}\right)$, evaluated at $p_\mathrm{top}$. Let $\Delta h_0=h(T,p_\mathrm{top})-h(T_c,p_\mathrm{top})$ be the $u\to0$ denominator. Then
$$J(T;F)=J_0(T)+\delta J(T)+O(F^4),\tag{24}$$
with
$$J_0=\int_{T_1}^{T}\frac{\kappa}{\Delta h_0}\,dT',\qquad
\delta J=\int_{T_1}^{T}\left[\left(\frac{\partial\kappa}{\partial p}\right)_T\frac{\delta p}{\Delta h_0}-\frac{\kappa\,\big(\Delta H-\Delta h_0\big)}{\Delta h_0^2}\right]dT',$$
$$\Delta H-\Delta h_0=\tfrac12\left(u^2-u_c^2\right)+h_p\,\delta p-h_{p,c}\,\delta p_c+O(F^4).$$

Both $\delta p$ and $\Delta H-\Delta h_0$ are $O(F^2)$, so $\delta J$ is too. Therefore
$$L_\mathrm{TR}=\frac{J_0+\delta J}{F}+O(F^3).\tag{25}$$

**Fully ionized closed form.** For $x=1$: $h=5k_BT/m_H$, $h_p=\partial\kappa/\partial p=0$, $\kappa=\kappa_0T^{5/2}$ and $\rho=m_Hp/(2k_BT)$, so $\delta J$ comes only from $\tfrac12u^2$. For $T\gg T_c$, $h_c$ and $u_c$ are negligible. Then
$$J_0=\frac{2\kappa_0m_H}{25k_B}\,T^{5/2},\qquad
\delta J=-\frac{4\kappa_0F^2}{175\,p^2}\,T^{7/2},$$
so
$$\frac{\delta J}{J_0}=-\frac17\,\frac{\rho u^2}{p}=-\frac{5}{21}\,M^2,\qquad M^2=\frac{u^2}{\tfrac53\,p/\rho}.\tag{26}$$

This gives
$$\boxed{\;F\,L_\mathrm{TR}\simeq J_0\left(1-\frac17\,\frac{\rho_\mathrm{top}u_\mathrm{top}^2}{p_\mathrm{top}}\right)=J_0\left(1-\tfrac{5}{21}M_\mathrm{top}^2\right)\;}$$

The $u^2$ term **shortens** the TR. Part of the conducted flux goes into accelerating the expanding gas, so less of it is available to heat the gas, and the temperature must rise over a shorter distance. In the partially ionized TR the measured coefficient is about 1.5–2.2 times $5/21$ (§8): $\Delta h$ is smaller than $h$ there, and the $h_p\,\delta p$ term in $\Delta H-\Delta h_0$ adds.

### 6c. Static radiative TR: $F=0$

From (9), $q^2=2\int_{T_c}^{T}\kappa\,Q\,dT'\equiv q_s^2$, with $q_s<0$. Hence
$$L_\mathrm{TR}(T)=\int_{T_1}^{T}\frac{\kappa\,dT'}{\left(2\int_{T_c}^{T'}\kappa Q\,dT''\right)^{1/2}} .\tag{27}$$

If $Q(T_c^+)>0$, that is, if the loss switches on at the chromosphere top as specified, then $|q_s|\propto(T-T_c)^{1/2}$ near the foot. The integrand $\propto(T-T_c)^{-1/2}$ is then integrable, so **$L_\mathrm{TR}$ from $T_1=T_c$ is finite**.

**Closed form** for an optically thin loss $Q=n_en_H\Lambda_0T^{\alpha}$ in the fully ionized limit, where $n_e=n_H=p/(2k_BT)$:
$$q_s^2=\frac{\kappa_0\Lambda_0p^2}{2k_B^2(\alpha+\tfrac32)}\left(T^{\alpha+3/2}-T_c^{\alpha+3/2}\right).$$

For $T\gg T_c$ and $\alpha<\tfrac{11}{2}$:
$$L_\mathrm{TR}\simeq\frac{4}{11-2\alpha}\left[\frac{2(\alpha+\tfrac32)\,\kappa_0}{\Lambda_0}\right]^{1/2}\frac{k_B}{p}\,T^{(11-2\alpha)/4}.\tag{28}$$

For $\alpha=-\tfrac12$ this is $L_\mathrm{TR}=(k_B/3p)\,(2\kappa_0/\Lambda_0)^{1/2}\,T^3$, the RTV-type scaling $T\propto(pL)^{1/3}$. In the partially ionized general case use $Q=x\,n_H^2\Lambda(T)$ with $n_H$ from (12)–(14) and integrate (27) numerically.

### 6d. General case: $F>0$ and $Q>0$

Equation (8) has no closed form in general. Integrate it from $T_c$, starting from the local solution $-q\simeq F\,c_{p,c}\,\epsilon+\left(2\kappa_cQ_c\,\epsilon\right)^{1/2}$ at $T=T_c+\epsilon$, and then evaluate (10).

**Bracket.** Let $q_e\equiv-F\,\Delta H$ (the $Q=0$ solution (22)) and $q_s$ (the $F=0$ solution of §6c); all three of $q$, $q_e$, $q_s$ are $\le0$. Assume $dH/dT>0$ and $Q\ge0$; both always hold here (§5, after (21)).

- *Lower bound.* (8) gives $d|q|/dT\ge F\,dH/dT$, so $|q|\ge|q_e|$. Then (9) gives $d(q^2)/dT\ge d(q_e^2)/dT+d(q_s^2)/dT$.
- *Upper bound.* Since $|q|\ge|q_s|$, $d|q|/dT\le d\big(|q_e|+|q_s|\big)/dT$.

Together:
$$\sqrt{q_e^2+q_s^2}\;\le\;|q|\;\le\;|q_e|+|q_s| ,\tag{29}$$
and therefore
$$\int_{T_1}^{T}\frac{\kappa\,dT'}{|q_e|+|q_s|}\;\le\;L_\mathrm{TR}(T)\;\le\;\int_{T_1}^{T}\frac{\kappa\,dT'}{\sqrt{q_e^2+q_s^2}} .\tag{30}$$

The two bounds differ by at most $\sqrt2$, and they coincide when either mechanism dominates. In particular, $L_\mathrm{TR}\le\min\left(J/F,\,L_\mathrm{static}\right)$: evaporation and radiation each shorten the TR.

### 6e. Radiation-dominated TR: an $F$-independent length

**Setup.** $Q>0$ in the TR ($T>T_c$) and $Q=0$ in the chromosphere. Mass conservation (1) makes $F$ the same in both layers: a steady state cannot carry $F\ne0$ in the chromosphere and $F=0$ in the TR. In the isothermal chromosphere $H$ is constant, so $F$ does not affect the temperature there whatever its value. The question is when $L_\mathrm{TR}$ stops depending on $F$.

**$F=0$: exactly $F$-independent.** With $F=0$, (2) gives $p=$ const exactly (gravity dropped), so constant pressure is a result here, not an assumption. Then (27), taken from $T_1=T_c$ (finite by §6c), is
$$L_s(T;p,T_c)=\int_{T_c}^{T}\frac{\kappa(T',p)\,dT'}{|q_s(T')|},\qquad |q_s(T)|=\left(2\int_{T_c}^{T}\kappa\,Q\,dT''\right)^{1/2}.\tag{31}$$
$L_s$ depends only on $T$, $p$, $T_c$ and the loss function $\Lambda$. With Saha, $\kappa$ and $Q=x\,n_H^2\Lambda(T)$ depend on $p$ through $x$ and $n_H$, so $L_s$ is a general function of $(T,p)$. In the fully ionized limit $L_s\propto1/p$, as in (28). This is the AWSoM-R lower boundary (Sokolov et al. 2021, ApJ 908, 172, Eqs. 20–21) with the release closures in place of pure Spitzer conduction.

**Small $F$: a correction linear in $F$.** Integrate (9) from $T_c$ with $q(T_c)=0$. Since $q<0$, $-F\,q\,dH/dT=F|q|\,dH/dT$, so exactly
$$q^2=q_s^2+2F\int_{T_c}^{T}|q|\,\frac{dH}{dT'}\,dT'.\tag{32}$$

The second term is already $O(F)$, so to first order $|q|$ inside it can be replaced by $|q_s|$. Also to first order, $p=$ const and $dH/dT=c_p$: the dynamic pressure and $u^2$ are $O(F^2)$ (§6b). With
$$I(T)\equiv\int_{T_c}^{T}|q_s|\,c_p\,dT' ,$$
(32) gives $|q|=|q_s|\left(1+F\,I/q_s^2\right)+O(F^2)$. Inserting this into (10):
$$\boxed{\;L_\mathrm{TR}(T)=L_s(T)-F\int_{T_c}^{T}\frac{\kappa\,I}{|q_s|^3}\,dT'+O(F^2)\;}\tag{33}$$

Near the foot $|q_s|\propto\epsilon^{1/2}$ and $I\propto\epsilon^{3/2}$, with $\epsilon=T-T_c$, so the correction integrand is finite there.

The correction is controlled by
$$F_*(T)\equiv\frac{|q_s(T)|}{\Delta H(T)},$$
the mass flux at which the enthalpy flux $F\,\Delta H$ equals the static radiative flux. The local relative change of $|q|$ is $F\,I/q_s^2\sim F/F_*$. So $L_\mathrm{TR}$ is independent of $F$ to the extent that $F\ll F_*$. In the opposite limit, $F\gg F_*$, radiation is negligible and $L_\mathrm{TR}\to J/F$ (§6a). In between, use the bracket (30) or integrate (8).

**Power-law closed form.** Take a fully ionized gas with $Q=n_en_H\Lambda_0T^{\alpha}$, $\kappa=\kappa_0T^{5/2}$ and $\Delta H=c_pT$ with $c_p=5k_B/m_H$, at $T\gg T_c$. Then $|q_s|\propto T^{(2\alpha+3)/4}$ and every integral in (33) is a power of $T$. For $-\tfrac32<\alpha<3$:
$$L_\mathrm{TR}\simeq L_s\left(1-C_\alpha\,\frac{F}{F_*(T)}\right),\qquad C_\alpha=\frac{11-2\alpha}{(2\alpha+7)(3-\alpha)} ,\tag{34}$$
with $C_{-1/2}=4/7$. Evaporation shortens the radiative TR, consistent with the upper bound in (30).

**What then fixes $F$.** For $F\ll F_*$, $L_\mathrm{TR}$ is set by $(T,p)$ and not by $F$, so the TR length can no longer select $F$ as it does in the release (§7). $F$ instead follows from the heat flux the corona delivers to the TR top. Evaluating $|q|=|q_s|+F\,I/|q_s|$ at the top:
$$F\simeq\frac{|q_s|\,\big(|q_\mathrm{top}|-|q_s|\big)}{I}\bigg|_{T_\mathrm{top}} .$$

$F$ is set by the heat flux in excess of what the static TR radiates. The denominator is $I/|q_s|$ and not $\Delta H$. In the power-law case $I/|q_s|=\Delta H/(1+\tfrac{2\alpha+3}{4})$, which is $\tfrac23\Delta H$ for $\alpha=-\tfrac12$. The reason: evaporation also thins the TR, which reduces the radiated power $\int Q\,ds$, so the conducted flux grows by less than $F\,\Delta H$. If $|q_\mathrm{top}|<|q_s|$, this gives $F<0$, a steady downflow, which is outside the $F>0$ assumption of this note.

## 7. Energy balance at the top and what fixes $F$

Evaluating (5) at the top:
$$-q_\mathrm{top}=F\,\big[H_\mathrm{top}-H_c\big]+\int_{s_c}^{s_\mathrm{top}}Q\,ds .\tag{35}$$

The steady family is one-parameter in $F$, both with and without $Q$. Equations (23) and (30) give $L_\mathrm{TR}(F)$; they do not select $F$. In the release, $Q=0$, $L_\mathrm{TR}$ is set by geometry (the chromosphere ends where its hydrostatic pressure reaches $p_\mathrm{top}$), and $F=J/L_\mathrm{TR}$ follows. With $Q>0$, $L_\mathrm{TR}(F)$ stays finite as $F\to0$ and tends to the static value (31), instead of diverging like $J/F$. For $F\ll F_*$ the length depends on $F$ only through the first-order correction (33), and $F$ is fixed instead by (35) through the heat flux at the top (§6e). In this sense a radiative sink is what makes $L_\mathrm{TR}$ a property of the TR rather than of the domain.

## 8. Numerical check

`util/release/tr_length_steady.py` evaluates (8)–(23) with exactly the release closures: the Saha constants of `eos.hpp` and $\kappa$ from `util/lib/conductive_flux.py`. Inputs: $T_\mathrm{top}=22$ kK, $p_\mathrm{top}=0.0103$ Pa, $T_c=6.7$ kK, $T_1=7.5$ kK.

**1. ODE (8) against the closed form (23)** at $F=2.683\times10^{-10}$ kg m⁻² s⁻¹ (the measured 4000 s top mass flux of the single-T release): both give $L_\mathrm{TR}=16.389$ km (relative difference $1.4\times10^{-7}$). The measured release $L_\mathrm{TR}$ is 13–15 km. That is consistent with the measured $F L/J\approx0.92$ at 4000 s.

**2. Effect of keeping $u^2$** (with the dynamic-pressure drop (19)), for $Q=0$:

| $T_\mathrm{top}$ | $F$ [kg m⁻² s⁻¹] | $u_\mathrm{top}$ [m/s] | $M_\mathrm{top}^2$ | $L_\mathrm{TR}$ [km] | $\Delta L/L$ from $u^2$ | $-\tfrac{5}{21}M_\mathrm{top}^2$ |
|---:|---:|---:|---:|---:|---:|---:|
| 22 kK | $2.68\times10^{-10}$ | 9.5 | $1.5\times10^{-7}$ | 16.39 | $-7.9\times10^{-8}$ | $-3.5\times10^{-8}$ |
| 22 kK | $10^{-8}$ | 353 | $2.1\times10^{-4}$ | 0.440 | $-1.1\times10^{-4}$ | $-4.9\times10^{-5}$ |
| 22 kK | $10^{-7}$ | 3530 | $2.1\times10^{-2}$ | 0.0435 | $-1.1\times10^{-2}$ | $-4.9\times10^{-3}$ |
| 100 kK | $2.68\times10^{-10}$ | 43 | $6.7\times10^{-7}$ | 960 | $-2.5\times10^{-7}$ | $-1.6\times10^{-7}$ |
| 100 kK | $10^{-8}$ | 1600 | $9.3\times10^{-4}$ | 25.8 | $-3.4\times10^{-4}$ | $-2.2\times10^{-4}$ |
| 1 MK | $10^{-8}$ | 16000 | $9.3\times10^{-3}$ | 8810 | $-2.4\times10^{-3}$ | $-2.2\times10^{-3}$ |

At 1 MK the fully ionized prediction (26) holds to 6 %. In the partially ionized TR the coefficient is 1.5–2.2 times larger, for the reasons given in §6b. At the release state the term is $\sim10^{-7}$, far below the discretization error. Keeping it is formally consistent but numerically inert. (Rows with $M^2\gtrsim0.1$ leave the range of the first-order expansion; the script also prints a near-sonic 1 MK row, which is not a valid test of (26).)

**3. Static radiative limit.** For a fully ionized loss with $\Lambda=10^{-32}\,T^{-1/2}$ W m³, a solver check only, the ODE matches the closed form (27)–(28) to $5\times10^{-8}$ ($L_\mathrm{TR}=1.947\times10^{4}$ km from 20 kK to 1 MK).

**4. Bracket (30)** with the same loss, $F>0$:

| $F$ [kg m⁻² s⁻¹] | lower bound [km] | ODE [km] | upper bound [km] |
|---:|---:|---:|---:|
| $10^{-8}$ | 6208 | 7097 | 8256 |
| $3\times10^{-8}$ | 2644 | 2872 | 3034 |
| $10^{-7}$ | 880 | 890 | 923 |

**5. First-order $F$ correction (33)–(34)** with the same loss, from 20 kK to 1 MK at $p=0.0103$ Pa: $L_s=1.947\times10^4$ km, $|q_s(T_\mathrm{top})|=158$ W m⁻², $F_*=3.92\times10^{-9}$ kg m⁻² s⁻¹. The last column recovers $F$ from the ODE's top flux with the formula at the end of §6e.

| $F$ [kg m⁻² s⁻¹] | $F/F_*$ | $\Delta L/L_s$ (ODE) | $-\tfrac47F/F_*$ (34) | $F$ from $q_\mathrm{top}$, rel. error |
|---:|---:|---:|---:|---:|
| $2.68\times10^{-10}$ | 0.068 | $-3.78\times10^{-2}$ | $-3.91\times10^{-2}$ | $+1.1\times10^{-2}$ |
| $10^{-10}$ | 0.026 | $-1.44\times10^{-2}$ | $-1.46\times10^{-2}$ | $+4.3\times10^{-3}$ |
| $3\times10^{-11}$ | 0.0077 | $-4.34\times10^{-3}$ | $-4.38\times10^{-3}$ | $+1.3\times10^{-3}$ |
| $10^{-11}$ | 0.0026 | $-1.45\times10^{-3}$ | $-1.46\times10^{-3}$ | $+4.3\times10^{-4}$ |

Both errors fall linearly with $F$, as expected for an $O(F^2)$ remainder. The 0.6 % residual in $\Delta L/L_s$ at the smallest $F$ is the finite $T_c/T_\mathrm{top}=0.02$, which (34) neglects. The first row uses the release mass flux but not a release state: this loss is a solver check, and $F_*$ for a real $\Lambda(T)$ at the release $(T,p)$ has not been evaluated.

```bash
cd util && ../.venv/bin/python tr_length_steady.py     # ~5 s
```

## 9. Limits

- **Gravity dropped.** This is valid for a thin TR, where $g\,\Delta z\lesssim0.5\,\%$ of $\Delta h$ (`twall_mass_flux_relaxation_timescales.md`). Restoring it adds $F\,\phi$ to $H$ in (3) and $-\rho g$ to the right side of (2).
- **$q_c=0$ at the chromosphere top** is the same idealization as in the existing $J$ table. $T_c$ remains a hidden variable to which $h_c$ is exponentially sensitive at low $p$ (`analytic_evaporation_J_table.md`, point 3).
- **Steady only.** The release TR relaxes over $t_\mathrm{ready}=520$–2300 s, and $F$ decays afterwards by chromospheric contraction.
- **$Q$ is not release physics.** Any quantitative use of §6c–6e needs a chosen $\Lambda(T)$ and is a model extension, not a property of the release.
