# An Approximate Linearised Riemann Solver for the Euler Equations for Real Gases*

**P. GLAISTER**

*Department of Mathematics, University of Reading, Reading, Berkshire RG6 2AX, United Kingdom*

Received September 26, 1986; revised February 20, 1987

Journal of Computational Physics 74, 382–408 (1988)

An approximate (linearised) Riemann solver is presented for the solution of the Euler equations of gas dynamics in one dimension with a general convex equation of state. The scheme is applied to a standard shock reflection test problem for some specimen equations of state. © 1988 Academic Press, Inc.

---

## 1. Introduction

The linearised approximate Riemann solver of Roe [1] was proposed in 1981 for the solution of the Euler equations of gas dynamics where the properties of the fluid are represented by the ideal equation of state. We seek here to extend this scheme to the solution of the Euler equations in one dimension for real gases. At each stage we shall as far as possible draw a parallel with Roe's scheme for the ideal equation of state. Results for the extended scheme are presented for a particular problem of shock reflection for three different equations of state.

In Section 2 we look at the Jacobian matrix of the flux function for the Euler equations with a general convex equation of state, and in Section 3 derive an approximate Riemann solver for the solution of these equations. In Section 4 we give some particular examples of nonideal equations of state, and in Section 5 we describe a standard test problem involving shock reflection. Finally, in Section 6 we display the numerical results achieved for this test problem with three different equations of state.

A construction with similar objectives has been proposed by Roe [2] which, however, differs in both procedure and final form.

---

## 2. Equations of Flow and State

In this section we state the equations of motion for an inviscid compressible fluid in one dimension for any equation of state, and derive the eigenvalues and eigenvectors of the Jacobian of the corresponding flux function.

### 2.1. Equations

The Euler equations governing the flow of an inviscid, compressible fluid in one dimension may be written in conservation form as

$$\mathbf{w}_{,t} + \mathbf{F}_{,x} = \mathbf{0}, \tag{2.1}$$

where

$$\mathbf{w} = (\rho, \rho u, e)^{\mathrm{T}} \tag{2.2}$$

$$\mathbf{F}(\mathbf{w}) = (\rho u,\ p + \rho u^2,\ u(e+p))^{\mathrm{T}}, \tag{2.3}$$

together with

$$e = \rho i + \tfrac{1}{2}\rho u^2, \tag{2.4}$$

where $\rho = \rho(x,t)$, $u = u(x,t)$, $p = p(x,t)$, $i = i(x,t)$, and $e = e(x,t)$ represent the density, velocity, pressure, specific internal energy, and the total energy, respectively, at a position $x$ and time $t$. Equations (2.1) represent conservation of mass, momentum, and energy. In addition, there is an equation of state which is a macroscopic, thermodynamic relationship specific to each particular fluid, and we assume here that this can be written in the form

$$p = p(\rho, i). \tag{2.5}$$

The function $p(\cdot,\cdot)$ will be assumed to satisfy conditions which ensure that the system (2.1) is hyperbolic and the corresponding Riemann problem always possesses a unique solution (see [3]). Furthermore, we shall assume that the first derivatives $\partial p/\partial\rho|_i$ and $\partial p/\partial i|_\rho$ are available. In the case of an ideal gas, Eq. (2.5) becomes

$$p = (\gamma - 1)\rho i, \tag{2.6}$$

where $\gamma$ is the ratio of specific heat capacities of the fluid: this is sometimes called a $\gamma$-gas law. The relationship given in Eq. (2.5) will usually be determined by experimental considerations.

### 2.2. Jacobian

We now construct the Jacobian, $A$, of the flux function, $\mathbf{F}(\mathbf{w})$ given by

$$A = \partial\mathbf{F}/\partial\mathbf{w}, \tag{2.7}$$

and find its eigenvalues and (right) eigenvectors since this will form the basis for our approximate Riemann solver.

Defining the momentum $m$ as $m = \rho u$ we may rewrite Eqs. (2.2), (2.3), and (2.5) in the form

$$\mathbf{w} = (\rho, m, e)^{\mathrm{T}} \tag{2.8a}$$

$$\mathbf{F}(\mathbf{w}) = \left(m,\ p + \frac{m^2}{\rho},\ \frac{me}{\rho} + \frac{mp}{\rho}\right)^{\mathrm{T}} \tag{2.8b}$$

and

$$p = p(\rho, i), \tag{2.8c}$$

where

$$i = \frac{e}{\rho} - \frac{1}{2}\frac{m^2}{\rho^2}. \tag{2.8d}$$

Now,

$$\frac{\partial\mathbf{F}}{\partial\mathbf{w}} = \left(\left.\frac{\partial\mathbf{F}}{\partial\rho}\right|_{m,e},\ \left.\frac{\partial\mathbf{F}}{\partial m}\right|_{\rho,e},\ \left.\frac{\partial\mathbf{F}}{\partial e}\right|_{\rho,m}\right) \tag{2.9}$$

and, in particular, we will need to find, $(\partial p/\partial\rho)(\rho, i(\rho,m,e))|_{m,e}$, $(\partial p/\partial m)(\rho, i(\rho,m,e))|_{\rho,e}$, and $(\partial p/\partial e)(\rho, i(\rho,m,e))|_{\rho,m}$. By the chain rule for partial derivatives, however, we have

$$\frac{\partial p}{\partial\rho}\big(\rho, i(\rho,m,e)\big)\Big|_{m,e} = \frac{\partial p}{\partial\rho}(\rho,i)\Big|_i + \frac{\partial i}{\partial\rho}(\rho,m,e)\Big|_{m,e}\ \frac{\partial p}{\partial i}(\rho,i)\Big|_\rho \tag{2.10a}$$

$$\frac{\partial p}{\partial m}\big(\rho, i(\rho,m,e)\big)\Big|_{\rho,e} = \frac{\partial i}{\partial m}(\rho,m,e)\Big|_{\rho,e}\ \frac{\partial p}{\partial i}(\rho,i)\Big|_\rho \tag{2.10b}$$

$$\frac{\partial p}{\partial e}\big(\rho, i(\rho,m,e)\big)\Big|_{\rho,m} = \frac{\partial i}{\partial e}(\rho,m,e)\Big|_{\rho,m}\ \frac{\partial p}{\partial i}(\rho,i)\Big|_\rho, \tag{2.10c}$$

where

$$i = i(\rho,m,e) = \frac{e}{\rho} - \frac{1}{2}\frac{m^2}{\rho^2}. \tag{2.11}$$

This leads to the following expression for the Jacobian

$$A = \begin{pmatrix}
0 & 1 & 0 \\[4pt]
a^2 - u^2 - \dfrac{p_i}{\rho}(H-u^2) & 2u - \dfrac{u p_i}{\rho} & \dfrac{p_i}{\rho} \\[10pt]
u(a^2 - H) - \dfrac{u p_i}{\rho}(H-u^2) & H - \dfrac{u^2 p_i}{\rho} & u + \dfrac{u p_i}{\rho}
\end{pmatrix}, \tag{2.12}$$

where the enthalpy, $H$, is defined by

$$H = \frac{e+p}{\rho} = \frac{p}{\rho} + i + \frac{1}{2}u^2, \tag{2.13}$$

the "sound speed," $a$, is given by

$$a^2 = \frac{p\,p_i}{\rho^2} + p_\rho, \tag{2.14}$$

and we use the shorthand notation $p_\rho \equiv (\partial p/\partial\rho)(\rho,i)|_i$, $p_i \equiv (\partial p/\partial i)(\rho,i)|_\rho$.

The eigenvalues, $\lambda_i$, and corresponding right eigenvectors, $\mathbf{e}_i$, of $A$ are then found to be

$$\lambda_1 = u+a,\qquad
\mathbf{e}_1 = \begin{pmatrix}1\\ u+a\\ H+ua\end{pmatrix}
= \begin{pmatrix}1\\ u+a\\ \dfrac{p}{\rho}+i+\dfrac{1}{2}u^2+ua\end{pmatrix}, \tag{2.15a}$$

$$\lambda_2 = u-a,\qquad
\mathbf{e}_2 = \begin{pmatrix}1\\ u-a\\ H-ua\end{pmatrix}
= \begin{pmatrix}1\\ u-a\\ \dfrac{p}{\rho}+i+\dfrac{1}{2}u^2-ua\end{pmatrix}, \tag{2.15b}$$

and

$$\lambda_3 = u,\qquad
\mathbf{e}_3 = \begin{pmatrix}1\\ u\\ H - \dfrac{\rho a^2}{p_i}\end{pmatrix}
= \begin{pmatrix}1\\ u\\ i+\dfrac{1}{2}u^2 - \dfrac{\rho p_\rho}{p_i}\end{pmatrix}. \tag{2.15c}$$

We note that in the case of an ideal gas the equation of state (2.8c) becomes

$$p = (\gamma-1)\rho i, \tag{2.16}$$

giving

$$p_i = (\gamma-1)\rho, \qquad p_\rho = (\gamma-1)i \tag{2.17}$$

and thus

$$\frac{a^2}{\gamma-1} = \frac{p}{\rho} + i = H - \frac{1}{2}u^2 = \frac{\gamma p}{\rho(\gamma-1)}. \tag{2.18}$$

In particular, the eigenvectors $\mathbf{e}_1$, $\mathbf{e}_2$, $\mathbf{e}_3$ become

$$\mathbf{e}_1 = \begin{pmatrix}1\\ u+a\\ \dfrac{a^2}{\gamma-1}+\dfrac{1}{2}u^2+ua\end{pmatrix},\qquad
\mathbf{e}_2 = \begin{pmatrix}1\\ u-a\\ \dfrac{a^2}{\gamma-1}+\dfrac{1}{2}u^2-ua\end{pmatrix},\qquad
\mathbf{e}_3 = \begin{pmatrix}1\\ u\\ \dfrac{1}{2}u^2\end{pmatrix}. \tag{2.19a-2.19c}$$

In the next section we develop an approximate Riemann solver using the results in this section.

---

## 3. An Approximate Linearised Riemann Solver

In this section we develop an approximate Riemann solver for the Euler equations in one dimension with a general convex equation of state. We follow a similar course of reasoning as that used by Roe and Pike [4] in the ideal gas case and begin by giving a brief description of their algorithm.

### 3.1. The Approximate (Linearised) Riemann Solver of Roe and Pike for an Ideal Gas

Given two states $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$ (left and right) of a gas close to an average state $\mathbf{w}$, seek coefficients $\alpha_1$, $\alpha_2$, $\alpha_3$ such that, if $\mathbf{e}_1$, $\mathbf{e}_2$, $\mathbf{e}_3$ are the eigenvectors of the Jacobian matrix for the ideal gas flux function (2.19a)–(2.19c)

$$\Delta\mathbf{w} = \sum_{j=1}^{3} \alpha_j \mathbf{e}_j \tag{3.1}$$

to within $O(\Delta^2)$, where $\Delta(\cdot) = (\cdot)_{\mathrm R} - (\cdot)_{\mathrm L}$. This gives the expressions

$$\alpha_1 = \frac{1}{2a^2}(\Delta p + \rho a\,\Delta u) \tag{3.2a}$$

$$\alpha_2 = \frac{1}{2a^2}(\Delta p - \rho a\,\Delta u) \tag{3.2b}$$

$$\alpha_3 = \Delta\rho - \frac{\Delta p}{a^2}, \tag{3.2c}$$

and it can also be shown that with the same values of $\alpha_1,\alpha_2,\alpha_3$

$$\Delta\mathbf{F} = \sum_{j=1}^{3}\lambda_j \alpha_j \mathbf{e}_j, \tag{3.3}$$

where $\lambda_1,\lambda_2,\lambda_3$ are given by Eqs. (2.15a)–(2.15c). The decomposition (3.1) yields exact characteristic fields to $O(\Delta^2)$. The approximate Riemann solver is then constructed by seeking averages $\tilde\rho$, $\tilde u$, $\tilde a$ such that, for states $\mathbf{w}_{\mathrm L}$ and $\mathbf{w}_{\mathrm R}$, not necessarily close,

$$\Delta\mathbf{w} = \sum_{j=1}^{3}\tilde\alpha_j\tilde{\mathbf{e}}_j \tag{3.4}$$

and

$$\Delta\mathbf{F} = \sum_{j=1}^{3}\tilde\lambda_j\tilde\alpha_j\tilde{\mathbf{e}}_j \tag{3.5}$$

(i.e. (3.1) and (3.3) with averaged values) hold, where now

$$\tilde\lambda_{1,2,3} = \tilde u+\tilde a,\ \tilde u-\tilde a,\ \tilde u \tag{3.6a-3.6c}$$

$$\mathbf{e}_1 = \begin{pmatrix}1\\ \tilde u+\tilde a\\ \dfrac{a^2}{\gamma-1}+\dfrac{1}{2}\tilde u^2+\tilde u\tilde a\end{pmatrix},\qquad
\mathbf{e}_2 = \begin{pmatrix}1\\ \tilde u-\tilde a\\ \dfrac{a^2}{\gamma-1}+\dfrac{1}{2}\tilde u^2-\tilde u\tilde a\end{pmatrix},\qquad
\mathbf{e}_3 = \begin{pmatrix}1\\ \tilde u\\ \dfrac{1}{2}\tilde u^2\end{pmatrix}, \tag{3.7a-3.7c}$$

and

$$\tilde\alpha_1 = \frac{1}{2\tilde a^2}(\Delta p + \tilde\rho\tilde a\,\Delta u) \tag{3.8a}$$

$$\tilde\alpha_2 = \frac{1}{2\tilde a^2}(\Delta p - \tilde\rho\tilde a\,\Delta u) \tag{3.8b}$$

$$\tilde\alpha_3 = \Delta\rho - \frac{\Delta p}{\tilde a^2}. \tag{3.8c}$$

The required averages are found to be

$$\tilde\rho = \sqrt{\rho_{\mathrm L}\rho_{\mathrm R}} \tag{3.9}$$

$$\tilde u = \frac{\sqrt{\rho_{\mathrm L}}u_{\mathrm L} + \sqrt{\rho_{\mathrm R}}u_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}} \tag{3.10}$$

and

$$\tilde a^2 = (\gamma-1)\left(\tilde H - \tfrac{1}{2}\tilde u^2\right), \tag{3.11}$$

where

$$\tilde H = \frac{\sqrt{\rho_{\mathrm L}}H_{\mathrm L} + \sqrt{\rho_{\mathrm R}}H_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}} \tag{3.12}$$

and $H = (e+p)/\rho$ is the enthalpy. The approximate Riemann solver can then be implemented in a finite difference scheme as follows (see [4]).

Suppose at time level $n$ the approximate solution consists of a set of piecewise constants

$$\mathbf{w} = \begin{cases}\mathbf{w}_{\mathrm L} & x \in (x_{\mathrm L} - \tfrac{1}{2}\Delta x,\ x_{\mathrm L}+\tfrac{1}{2}\Delta x)\\ \mathbf{w}_{\mathrm R} & x \in (x_{\mathrm R} - \tfrac{1}{2}\Delta x,\ x_{\mathrm R}+\tfrac{1}{2}\Delta x),\end{cases} \tag{3.13}$$

where $\Delta x = x_{\mathrm R} - x_{\mathrm L}$ represents a constant mesh spacing. Thus at either end of the cell $(x_{\mathrm L}, x_{\mathrm R})$ the data is $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$. The solution $\mathbf{w}$ may be updated to time level $n+1$ in an upwind manner as shown schematically in Fig. 1, where $\Delta t$ is the time interval from level $n$ to level $n+1$. This approximate Riemann solver has the important shock-capturing property guaranteed by Eqs. (3.4)–(3.5) (see [1]).

We now use a similar course of reasoning to construct the linearised approximate Riemann solver for a general convex equation of state.

### 3.2. Wavespeeds for Nearby States

Consider two states $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$ (left and right) close to an average state $\mathbf{w}$, and seek $\alpha_1$, $\alpha_2$, $\alpha_3$ such that

$$\Delta\mathbf{w} = \sum_{j=1}^{3}\alpha_j\mathbf{e}_j \tag{3.14}$$

to within $O(\Delta^2)$, where $\Delta(\cdot) = (\cdot)_{\mathrm R} - (\cdot)_{\mathrm L}$ (cf. (3.1)). Writing Eq. (3.14) out in full we have

$$\Delta\rho = \alpha_1+\alpha_2+\alpha_3 \tag{3.15a}$$

$$\Delta(\rho u) = \alpha_1(u+a) + \alpha_2(u-a) + \alpha_3 u \tag{3.15b}$$

$$\Delta e = \alpha_1\left(\frac{p}{\rho}+i+\frac{1}{2}u^2+ua\right) + \alpha_2\left(\frac{p}{\rho}+i+\frac{1}{2}u^2-ua\right) + \alpha_3\left(i+\frac{1}{2}u^2 - \frac{\rho p_\rho}{p_i}\right). \tag{3.15c}$$

From Eqs. (3.15a)–(3.15b) we have that

$$\Delta(\rho u) - u\,\Delta\rho = a(\alpha_1-\alpha_2) \tag{3.16}$$

and from Eqs. (3.15a) and (3.15c),

$$\Delta(\rho i) - i\,\Delta\rho + \Delta\!\left(\frac{\rho u^2}{2}\right) - \frac{1}{2}u^2\Delta\rho = \frac{p}{\rho}(\alpha_1+\alpha_2) + ua(\alpha_1-\alpha_2) - \alpha_3\frac{\rho p_\rho}{p_i}. \tag{3.17}$$

Using Eq. (3.16) together with $\alpha_1+\alpha_2 = \Delta\rho - \alpha_3$, Eq. (3.17) yields the following equation for $\alpha_3$:

$$\frac{\alpha_3}{p_i}\left(\frac{p p_i}{\rho} + \rho p_\rho\right) = i\,\Delta\rho - \Delta(\rho i) + \frac{p}{\rho}\Delta\rho - \frac{u^2}{2}\Delta\rho - \Delta\!\left(\frac{\rho u^2}{2}\right) + u\,\Delta(\rho u). \tag{3.18}$$

Then, since

$$\rho a^2 = \frac{p p_i}{\rho} + \rho p_\rho, \tag{3.19}$$

$\alpha_3$ is given by

$$\rho a^2 \frac{\alpha_3}{p_i} = i\,\Delta\rho - \Delta(\rho i) + \frac{p}{\rho}\Delta\rho - \frac{u^2}{2}\Delta\rho - \Delta\!\left(\frac{\rho u^2}{2}\right) + u\,\Delta(\rho u). \tag{3.20a}$$

The coefficients $\alpha_1$ and $\alpha_2$ can now be calculated from Eqs. (3.15a) and (3.16), i.e.,

$$\alpha_1+\alpha_2 = \Delta\rho - \alpha_3 \tag{3.20b}$$

$$\alpha_1-\alpha_2 = \frac{\Delta(\rho u) - u\Delta\rho}{a}. \tag{3.20c}$$

We have made the assumption that the left and right states $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$ are close to some average state $\mathbf{w}$ to within $O(\Delta^2)$, so that, to this approximation

$$\Delta(\rho u) = u\,\Delta\rho + \rho\,\Delta u \tag{3.21a}$$

$$\Delta(\rho i) = i\,\Delta\rho + \rho\,\Delta i \tag{3.21b}$$

$$\Delta(\rho u^2) = u^2\Delta\rho + 2\rho u\,\Delta u. \tag{3.21c}$$

In that case Eq. (3.20a) gives

$$\rho a^2 \frac{\alpha_3}{p_i} = \frac{p}{\rho}\Delta\rho - \rho\,\Delta i, \tag{3.22}$$

and using Eq. (3.19) we obtain

$$\alpha_3 = \Delta\rho - \frac{(p_\rho\,\Delta\rho + p_i\,\Delta i)}{a^2}. \tag{3.23}$$

But

$$\Delta p = p_\rho\,\Delta\rho + p_i\,\Delta i \tag{3.24}$$

to within $O(\Delta^2)$, and therefore

$$\alpha_3 = \Delta\rho - \frac{\Delta p}{a^2}. \tag{3.25}$$

Finally, Eqs. (3.20b)–(3.20c) become

$$\alpha_1+\alpha_2 = \Delta p/a^2 \tag{3.26}$$

and

$$\alpha_1-\alpha_2 = \rho\,\Delta u/a, \tag{3.27}$$

to give the following expressions for $\alpha_1$, $\alpha_2$, and $\alpha_3$,

$$\alpha_1 = \frac{1}{2a^2}(\Delta p + \rho a\,\Delta u) \tag{3.28a}$$

$$\alpha_2 = \frac{1}{2a^2}(\Delta p - \rho a\,\Delta u) \tag{3.28b}$$

and

$$\alpha_3 = \Delta\rho - \Delta p/a^2. \tag{3.28c}$$

We have found $\alpha_1$, $\alpha_2$, $\alpha_3$ such that

$$\Delta\mathbf{w} = \sum_{j=1}^{3}\alpha_j\mathbf{e}_j \tag{3.29}$$

to within $O(\Delta^2)$, and a routine calculation verifies that

$$\Delta\mathbf{F} = \sum_{j=1}^{3}\lambda_j\alpha_j\mathbf{e}_j \tag{3.30}$$

to within $O(\Delta^2)$. We are now in a position to construct the new approximate Riemann solver.

### 3.3. Decomposition for General $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$

As in Roe and Pike [4], we consider the algebraic problem of finding average eigenvalues $\tilde\lambda_1$, $\tilde\lambda_2$, $\tilde\lambda_3$ and corresponding average eigenvectors $\tilde{\mathbf{e}}_1$, $\tilde{\mathbf{e}}_2$, $\tilde{\mathbf{e}}_3$ such that the relations (3.29) and (3.30) hold exactly for arbitrary states $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$ not necessarily close. Specifically, we seek averages $\tilde\rho$, $\tilde u$, $\tilde p_i$, $\tilde p_\rho$, $\tilde p$, and $\tilde\imath$ in terms of two adjacent states $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$ such that

$$\Delta\mathbf{w} = \sum_{j=1}^{3}\tilde\alpha_j\tilde{\mathbf{e}}_j \tag{3.31}$$

and

$$\Delta\mathbf{F} = \sum_{j=1}^{3}\tilde\lambda_j\tilde\alpha_j\tilde{\mathbf{e}}_j, \tag{3.32}$$

where

$$\Delta(\cdot) = (\cdot)_{\mathrm R} - (\cdot)_{\mathrm L} \tag{3.33a}$$

$$\mathbf{w} = (\rho,\rho u, e)^{\mathrm T} \tag{3.33b}$$

$$\mathbf{F}(\mathbf{w}) = (\rho u,\ p+\rho u^2,\ u(e+p))^{\mathrm T} \tag{3.33c}$$

$$e = \rho i + \tfrac{1}{2}\rho u^2 \tag{3.33d}$$

$$p = p(\rho, i) \tag{3.33e}$$

$$\tilde\lambda_{1,2,3} = \tilde u+\tilde a,\ \tilde u-\tilde a,\ \tilde u \tag{3.34a}$$

$$\tilde{\mathbf{e}}_{1,2,3} = \begin{pmatrix}1\\ \tilde u+\tilde a\\ \dfrac{\tilde p}{\tilde\rho}+\tilde\imath+\dfrac{1}{2}\tilde u^2+\tilde u\tilde a\end{pmatrix},\quad
\begin{pmatrix}1\\ \tilde u-\tilde a\\ \dfrac{\tilde p}{\tilde\rho}+\tilde\imath+\dfrac{1}{2}\tilde u^2-\tilde u\tilde a\end{pmatrix},\quad
\begin{pmatrix}1\\ \tilde u\\ \tilde\imath+\dfrac{1}{2}\tilde u^2-\dfrac{\tilde\rho\tilde p_\rho}{\tilde p_i}\end{pmatrix} \tag{3.34b}$$

$$\tilde\alpha_1 = \frac{1}{2\tilde a^2}(\Delta p + \tilde\rho\tilde a\,\Delta u) \tag{3.35a}$$

$$\tilde\alpha_2 = \frac{1}{2\tilde a^2}(\Delta p - \tilde\rho\tilde a\,\Delta u) \tag{3.35b}$$

$$\tilde\alpha_3 = \Delta\rho - \frac{\Delta p}{\tilde a^2}, \tag{3.35c}$$

and $\tilde a$ is given by

$$\tilde\rho\tilde a^2 = \frac{\tilde\rho\tilde p_i}{\tilde\rho} + \tilde\rho\tilde p_\rho. \tag{3.36}$$

*(N.B. As printed, Eq. (3.36) reads $\tilde\rho\tilde a^2 = \dfrac{\tilde p\tilde p_i}{\tilde\rho} + \tilde\rho\tilde p_\rho$ — i.e. the first term's numerator is $\tilde p\,\tilde p_i$, mirroring Eq. (2.14)/(3.19); the scanned glyph for $\tilde p$ vs $\tilde\rho$ is faint in the numerator, flagged as [?] and resolved by analogy with (3.19).)*

The problem of finding averages $\tilde\rho$, $\tilde u$, $\tilde p_i$, $\tilde p_\rho$, $\tilde p$, and $\tilde\imath$ subject to Eqs. (3.31)–(3.36) will subsequently be denoted by (\*). (N.B. The quantities $\tilde p_i$ and $\tilde p_\rho$ denote approximations to the partial derivatives $p_i$ and $p_\rho$, respectively.)

The solution of problem (\*) will be sought in a way similar to that adopted by Roe and Pike [4] in the specialised, ideal gas case (see Section 3.1). We note, however, that problem (\*) is equivalent to seeking an approximation $\tilde A$ to the Jacobian $A$ with eigenvalues $\tilde\lambda_i$ and eigenvectors $\tilde{\mathbf{e}}_i$, which is an alternative approach also used in the ideal gas case by Roe [1].

The first step in the analysis of problem (\*) is to write out Eqs. (3.31) and (3.32) explicitly, namely,

$$\Delta\rho = \tilde\alpha_1+\tilde\alpha_2+\tilde\alpha_3 \tag{3.37a}$$

$$\Delta(\rho u) = \tilde\alpha_1(\tilde u+\tilde a) + \tilde\alpha_2(\tilde u-\tilde a) + \tilde\alpha_3\tilde u \tag{3.37b}$$

$$\begin{aligned}
\Delta e = \Delta(\rho i) + \Delta\!\left(\frac{\rho u^2}{2}\right) &= \tilde\alpha_1\left(\frac{\tilde p}{\tilde\rho}+\tilde\imath+\frac{1}{2}\tilde u^2+\tilde u\tilde a\right)\\
&\quad + \tilde\alpha_2\left(\frac{\tilde p}{\tilde\rho}+\tilde\imath+\frac{1}{2}\tilde u^2-\tilde u\tilde a\right) + \tilde\alpha_3\left(\tilde\imath+\frac{1}{2}\tilde u^2-\frac{\tilde\rho\tilde p_\rho}{\tilde p_i}\right)
\end{aligned} \tag{3.37c}$$

$$\Delta(\rho u) = \tilde\alpha_1(\tilde u+\tilde a) + \tilde\alpha_2(\tilde u-\tilde a) + \tilde\alpha_3\tilde u \tag{3.37d}$$

$$\Delta(p+\rho u^2) = \Delta p + \Delta(\rho u^2) = \tilde\alpha_1(\tilde u+\tilde a)^2 + \tilde\alpha_2(\tilde u-\tilde a)^2 + \tilde\alpha_3\tilde u^2 \tag{3.37e}$$

and

$$\begin{aligned}
\Delta(u(e+p)) = \Delta(\rho u i) + \Delta\!\left(\frac{\rho u^3}{2}\right) + \Delta(up) &= \tilde\alpha_1(\tilde u+\tilde a)\left(\frac{\tilde p}{\tilde\rho}+\tilde\imath+\frac{1}{2}\tilde u^2+\tilde u\hat a\right)\\
&\quad + \tilde\alpha_2(\tilde u-\tilde a)\left(\frac{\tilde p}{\tilde\rho}+\tilde\imath+\frac{1}{2}\tilde u^2-\tilde u\tilde a\right) + \tilde\alpha_3\tilde u\left(\tilde\imath+\frac{1}{2}\tilde u^2-\frac{\tilde\rho\tilde p_\rho}{\tilde p_i}\right)
\end{aligned} \tag{3.37f}$$

Equation (3.37a) is satisfied by any average we care to define, while Eq. (3.37b) is the same as Eq. (3.37d). Thus it remains to satisfy Eqs. (3.37c)–(3.37f). From Eq. (3.37d) we have

$$\Delta(\rho u) = \tilde u(\tilde\alpha_1+\tilde\alpha_2+\tilde\alpha_3) + \tilde a(\tilde\alpha_1-\tilde\alpha_2) = \tilde u\,\Delta\rho + \tilde\rho\,\Delta u, \tag{3.38}$$

and from Eq. (3.37e) we obtain

$$\Delta(\rho u^2) = \tilde u^2(\tilde\alpha_1+\tilde\alpha_2+\tilde\alpha_3) + 2\tilde u\tilde a(\tilde\alpha_1-\tilde\alpha_2) = \tilde u^2\Delta\rho + 2\tilde u\tilde\rho\,\Delta u. \tag{3.39}$$

Substituting for $\tilde\rho$ from Eq. (3.38) into Eq. (3.39) yields the quadratic equation for $\tilde u$,

$$\tilde u^2\Delta\rho - 2\tilde u\,\Delta(\rho u) + \Delta(\rho u^2) = 0. \tag{3.40}$$

Only one solution of Eq. (3.40) is productive, namely,

$$\tilde u = \frac{\Delta(\rho u) - \sqrt{(\Delta(\rho u))^2 - \Delta\rho\,\Delta(\rho u^2)}}{\Delta\rho}$$

and a routine calculation yields

$$\tilde u = \frac{\sqrt{\rho_{\mathrm L}}u_{\mathrm L} + \sqrt{\rho_{\mathrm R}}u_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}}, \tag{3.41}$$

which, on substituting $\tilde u$ into Eq. (3.38), gives

$$\tilde\rho = \frac{\Delta(\rho u) - \tilde u\,\Delta\rho}{\Delta u} = \sqrt{\rho_{\mathrm L}\rho_{\mathrm R}}. \tag{3.42}$$

We have now determined $\tilde\rho$ and $\tilde u$, and with these we can show that

$$\Delta\!\left(\frac{\rho u^3}{2}\right) - \frac{\tilde u^3}{2}\Delta\rho - 3\frac{\tilde\rho\tilde u^2}{2}\Delta u = \frac{(\Delta u)^3\tilde\rho^2}{2(\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}})^2} \tag{3.43}$$

$$\Delta(up) - \tilde u\,\Delta p = \tilde\rho\,\Delta u\,\frac{(\sqrt{\rho_{\mathrm L}}(p_{\mathrm L}/\rho_{\mathrm L}) + \sqrt{\rho_{\mathrm R}}(p_{\mathrm R}/\rho_{\mathrm R}))}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}}, \tag{3.44}$$

and

$$\frac{\sqrt{\rho_{\mathrm L}}u_{\mathrm L}^2+\sqrt{\rho_{\mathrm R}}u_{\mathrm R}^2}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}} - \tilde u^2 = \frac{\tilde\rho(\Delta u)^2}{(\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}})^2}, \tag{3.45}$$

all of which will be used later.

We are now left with Eqs. (3.37c) and (3.37f) and begin by rewriting them, using Eqs. (3.35)–(3.36), to give

$$\Delta(\rho i) - \tilde\imath\,\Delta\rho - \frac{\tilde\rho\,\Delta p}{\tilde\rho\tilde a^2} + \tilde\alpha_3\tilde\rho\,\frac{\tilde p_\rho}{\tilde p_i} = 0 \tag{3.46}$$

and

$$\begin{aligned}
&\Delta(\rho u i) - \tilde u\tilde\imath\,\Delta\rho - \tilde\rho\tilde\imath\,\Delta u + \Delta(up) - \tilde u\,\Delta p - \tilde\rho\,\Delta u + \Delta\!\left(\frac{\rho u^3}{2}\right) - \frac{\tilde u^3}{2}\Delta\rho - \frac{3}{2}\tilde\rho\tilde u^2\Delta u\\
&\qquad - \frac{\tilde u\tilde\rho\,\Delta p}{\tilde\rho\tilde a^2} + \tilde\alpha_3\tilde u\tilde\rho\,\frac{\tilde p_\rho}{\tilde p_i} = 0.
\end{aligned} \tag{3.47}$$

Now, subtracting Eq. (3.46) multiplied by $\tilde u$ from Eq. (3.47) and using Eqs. (3.43)–(3.45) together with the identity

$$\Delta(\rho u i) - \tilde u\Delta(\rho i) = \tilde\rho\,\Delta u\,\frac{\sqrt{\rho_{\mathrm L}}i_{\mathrm L}+\sqrt{\rho_{\mathrm R}}i_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}},$$

we obtain, after division by $\tilde\rho\,\Delta u$,

$$\frac{\tilde p}{\tilde\rho}+\tilde\imath+\frac{1}{2}\tilde u^2 = \left(\sqrt{\rho_{\mathrm L}}\left(\frac{p_{\mathrm L}}{\rho_{\mathrm L}}+i_{\mathrm L}+\frac{1}{2}u_{\mathrm L}^2\right) + \sqrt{\rho_{\mathrm R}}\left(\frac{p_{\mathrm R}}{\rho_{\mathrm R}}+i_{\mathrm R}+\frac{1}{2}u_{\mathrm R}^2\right)\right)\Big/(\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}). \tag{3.48}$$

Therefore, if we define a mean enthalpy, $\tilde H$, by

$$\tilde H = \frac{\tilde p}{\tilde\rho}+\tilde\imath+\frac{1}{2}\tilde u^2, \tag{3.49}$$

we find, from Eq. (3.48), that

$$\tilde H = \frac{\sqrt{\rho_{\mathrm L}}H_{\mathrm L}+\sqrt{\rho_{\mathrm R}}H_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}}, \tag{3.50}$$

as in the ideal case. We have now specified $\tilde\rho$, $\tilde u$, $\tilde p/\tilde\rho+\tilde\imath$: thus, in order to specify $\tilde p_i$, $\tilde p_\rho$, $\tilde\imath$ (and hence $\tilde p$), we focus our attention on Eq. (3.46) which can be written as

$$\Delta(\rho i) - \tilde\imath\,\Delta\rho - \tilde\rho\,\Delta i + \frac{\tilde\rho}{\tilde p_i}\left(\tilde p_i\,\Delta i + \tilde p_\rho\,\Delta\rho - \Delta p\right) = 0. \tag{3.51}$$

A number of choices can now be made, but it is clear that the most natural choice is to take

$$\Delta(\rho i) - \tilde\imath\,\Delta\rho - \tilde\rho\,\Delta i = 0, \tag{3.52}$$

i.e.,

$$\tilde\imath = \frac{\Delta(\rho i) - \tilde\rho\,\Delta i}{\Delta\rho} = \frac{\sqrt{\rho_{\mathrm L}}i_{\mathrm L}+\sqrt{\rho_{\mathrm R}}i_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}}, \tag{3.53}$$

in which case (3.51) gives

$$\Delta p = \tilde p_i\,\Delta i + \tilde p_\rho\,\Delta\rho \tag{3.54}$$

as a necessary condition. Therefore, all we need to do to complete the approximate Riemann solver is to choose approximations $\tilde p_i$, $\tilde p_\rho$ to $p_i$, $p_\rho$ such that (3.54) holds. This is a straightforward matter.

We propose approximations $\tilde p_i$, $\tilde p_\rho$ to $p_i$ and $p_\rho$ as follows:

$$\tilde p_i = \begin{cases}
\dfrac{1}{\Delta i}\left(\dfrac{1}{2}\big[p(\rho_{\mathrm R},i_{\mathrm R})+p(\rho_{\mathrm L},i_{\mathrm R})\big] - \dfrac{1}{2}\big[p(\rho_{\mathrm R},i_{\mathrm L})+p(\rho_{\mathrm L},i_{\mathrm L})\big]\right) & \text{if } \Delta i \neq 0\\[8pt]
\dfrac{1}{2}\big[p_i(\rho_{\mathrm L},i) + p_i(\rho_{\mathrm R},i)\big] & \text{if } \Delta i = 0,\ i_{\mathrm L}=i_{\mathrm R}=i
\end{cases} \tag{3.55a-3.55b}$$

and

$$\tilde p_\rho = \begin{cases}
\dfrac{1}{\Delta\rho}\left(\dfrac{1}{2}\big[p(\rho_{\mathrm R},i_{\mathrm R})+p(\rho_{\mathrm R},i_{\mathrm L})\big] - \dfrac{1}{2}\big[p(\rho_{\mathrm L},i_{\mathrm R})+p(\rho_{\mathrm L},i_{\mathrm L})\big]\right) & \text{if } \Delta\rho \neq 0\\[8pt]
\dfrac{1}{2}\big[p_\rho(\rho,i_{\mathrm L}) + p_\rho(\rho,i_{\mathrm R})\big] & \text{if } \Delta\rho = 0,\ \rho_{\mathrm L}=\rho_{\mathrm R}=\rho
\end{cases} \tag{3.56a-3.56b}$$

It is a simple matter to check that, for each of the combinations arising from the approximations given by Eqs. (3.55a)–(3.56b), Eq. (3.54) is satisfied. In particular, if the equation of state is separable, i.e., consists of a series of terms of the form $p = \mathcal{R}(\rho)I(i)$ where $\mathcal{R}$, $I$ depend on $\rho$, $i$, respectively, then Eqs. (3.55a)–(3.56b) become

$$\tilde p_i = \begin{cases}\bar{\mathcal{R}}\dfrac{\Delta I}{\Delta i} & \text{if } \Delta i \neq 0\\[6pt] \bar{\mathcal{R}}I'(i) & \text{if } \Delta i=0,\ i_{\mathrm L}=i_{\mathrm R}=i\end{cases} \tag{3.57a-3.57b}$$

and

$$\tilde p_\rho = \begin{cases}\bar I\dfrac{\Delta\mathcal{R}}{\Delta\rho} & \text{if } \Delta\rho\neq 0\\[6pt] \bar I\mathcal{R}'(\rho) & \text{if } \Delta\rho=0,\ \rho_{\mathrm L}=\rho_{\mathrm R}=\rho,\end{cases} \tag{3.58a-3.58b}$$

where $\Delta(\cdot) = (\cdot)_{\mathrm R}-(\cdot)_{\mathrm L}$ as before, and $\bar\cdot = \tfrac{1}{2}[(\cdot)_{\mathrm L}+(\cdot)_{\mathrm R}]$, the arithmetic mean. Although Eqs. (3.55a)–(3.56b) are not the only choices for $\tilde p_i$, $\tilde p_\rho$, these expressions represent a natural extension of the approximations given by Eqs. (3.57a)–(3.58b). In particular, we note that for any particular equation of state equations (3.55a)–(3.56b) can be simplified and the resulting expressions can be incorporated into a finite difference code in such a way as to avoid function evaluations.

Summarising, we can implement the above one-dimensional Riemann solver for the Euler equations with a general convex equation of state in a finite difference scheme in a similar way to that of Roe and Pike [4] as follows. Suppose at time level $n$ we have data $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$ given at either end of the cell $(x_{\mathrm L}, x_{\mathrm R})$. Then we update $\mathbf{w}$ to time level $n+1$ in an upwind manner (cf. Section 3.1). Schematically, we increment $\mathbf{w}$ as in Fig. 1, where $\Delta x = x_{\mathrm R}-x_{\mathrm L}$, $\Delta t$ is the time interval from level $n$ to $n+1$, and $\tilde\lambda_j$, $\tilde\alpha_j$, $\tilde{\mathbf{e}}_j$ are given by

$$\tilde\lambda_{1,2,3} = \tilde u+\tilde a,\ \tilde u-\tilde a,\ \tilde u$$

$$\tilde{\mathbf{e}}_{1,2,3} = \begin{pmatrix}1\\ \tilde u+\tilde a\\ \dfrac{\tilde p}{\tilde\rho}+\tilde\imath+\dfrac{1}{2}\tilde u^2+\tilde u\tilde a\end{pmatrix},\quad
\begin{pmatrix}1\\ \tilde u-\tilde a\\ \dfrac{\tilde p}{\tilde\rho}+\tilde\imath+\dfrac{1}{2}\tilde u^2-\tilde u\tilde a\end{pmatrix},\quad
\begin{pmatrix}1\\ \tilde u\\ \tilde\imath+\dfrac{1}{2}\tilde u^2-\tilde\rho\dfrac{\tilde p_\rho}{\tilde p_i}\end{pmatrix}$$

$$\tilde\alpha_{1,2,3} = \frac{1}{2\tilde a^2}(\Delta p+\tilde\rho\tilde a\,\Delta u),\ \frac{1}{2\tilde a^2}(\Delta p-\tilde\rho\tilde a\,\Delta u),\ \Delta\rho-\frac{\Delta p}{\tilde a^2}$$

$$\tilde\rho = \sqrt{\rho_{\mathrm L}\rho_{\mathrm R}},\qquad \tilde u = \frac{\sqrt{\rho_{\mathrm L}}u_{\mathrm L}+\sqrt{\rho_{\mathrm R}}u_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}}$$

$$\tilde\imath = \frac{\sqrt{\rho_{\mathrm L}}i_{\mathrm L}+\sqrt{\rho_{\mathrm R}}i_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}},\qquad H = \frac{\sqrt{\rho_{\mathrm L}}H_{\mathrm L}+\sqrt{\rho_{\mathrm R}}H_{\mathrm R}}{\sqrt{\rho_{\mathrm L}}+\sqrt{\rho_{\mathrm R}}}$$

$$\tilde p = \tilde\rho\left(\tilde H-\tilde\imath-\tfrac{1}{2}\tilde u^2\right),\qquad \tilde a^2 = \frac{\tilde p\tilde p_i}{\tilde\rho^2}+\tilde p_\rho,$$

$\tilde p_i$, $\tilde p_\rho$ are given by Eqs. (3.55a)–(3.56b), and $\Delta(\cdot)=(\cdot)_{\mathrm R}-(\cdot)_{\mathrm L}$. In addition, we can use the idea of flux limiters [5] to create a second-order algorithm which is

**FIG. 1.** Schematic representation of the first order algorithm. (Diagram: at time level $n$, cell interface data $\mathbf{w}_{\mathrm L}$, $\mathbf{w}_{\mathrm R}$ update to level $n+1$; contributions $-\dfrac{\Delta t}{\Delta x}\tilde\lambda_j\tilde\alpha_j\tilde{\mathbf{e}}_j$ propagate from the right when $\tilde\lambda_j>0$ and $-\dfrac{\Delta t}{\Delta x}\tilde\lambda_j\tilde\alpha_j\tilde{\mathbf{e}}_j$ propagate from the left when $\tilde\lambda_j<0$, for $j=1,2,3$.)


# Glaister (1988), "An approximate linearised Riemann solver for the Euler equations for real gases", J. Comput. Phys. 74, 382–408 — OCR transcription of PDF pages 15–27 (printed pages 396–408)

---

## [p. 396]

oscillation-free, and we can modify the scheme to disperse entropy violating solutions (see [6]).

The Riemann solver we have constructed in this section is a conservative algorithm and has the important shock-capturing property guaranteed by Eqs. (3.31)–(3.32) (see [1]). In the next section we give examples of different equations of state.

### 4. Equations of State

In this section we give three different forms of the equation of state for a fluid.

**(a) *Ideal gas equation of state.*** This can be written in the general form

$$p = (\gamma - 1)\,\rho i, \tag{4.1}$$

where $\gamma$ is a constant and represents the ratio of specific heat capacities of the fluid. Typical values for $\gamma$ are $\gamma = \tfrac{5}{3}$ for a monatomic gas, e.g., helium, and $\gamma = 1.4$ for a diatomic gas, e.g., air.

**(b) *Stiffened equation of state.*** This is usually written in the form

$$p = B\left(\frac{\rho}{\rho_0} - 1\right) + (\gamma - 1)\,\rho i, \tag{4.2}$$

where $B$ is a constant, and $\rho_0$ represents a reference density. This form of the equation of state is a simple extension of the ideal gas equation, and as such can be used in test problems originally designed for ideal gases.

**(c) *General equation of state.*** A more general equation of state has been developed by R. K. Osborne at the Los Alamos Scientific Laboratory [7], and can be written in the form

$$
p = [1/(E+\phi_0)]\{\zeta(a_1 + a_2|\zeta|)
$$
$$
{}+ E[b_0 + \zeta(b_1 + b_2\zeta) + E(c_0 + c_1\zeta)]\}, \tag{4.3}
$$

where $E = \rho_0 i$, $\zeta = \rho/\rho_0 - 1$ and the constants $\rho_0, a_1, a_2, b_0, b_1, b_2, c_0, c_1, \phi_0$ depend on the material in question. Typical values for the material constants for copper are given in Section 6.

Our algorithm requires knowledge of the derivatives $p_i, p_\rho$ which can be explicitly calculated in each of the three cases (a), (b), (c). The most general equations of state may be presented in tabular form, but provided that data is available for $p$, $p_i$ and $p_\rho$, we can always apply our algorithm as in cases (a)–(c).

In the next section we describe a standard test problem for the Euler equations with a general convex equation of state.

---

## [p. 397]

### 5. A Test Problem

In this section we describe a standard test problem in gas dynamics.

The test problem we consider is concerned with shock reflection in one dimension of a gas governed by the Euler equations with a general equation of state. We consider a region $0 \le x \le 1$ with initial conditions (at $t = 0$),

$$
\begin{aligned}
\rho &= \rho_0\\
u &= -u_0,\\
i &= i_0,
\end{aligned} \tag{5.1}
$$

where $p_0 = p(\rho_0, i_0)$ is given. This represents a gas of constant density and pressure moving towards $x=0$ (see [8]). The boundary $x=0$ is a rigid wall and the exact solution describes shock reflection from the wall. The gas is brought to rest at $x=0$ and, denoting initial values by (0), pre-shocked values by $(-)$, and post-shocked values by $(+)$, we can postulate an exact solution of the form

$$
\rho = \rho^{+}, \quad u = u^{+} = 0, \quad i = i^{+}, \quad (p = p^{+} = p(\rho^{+}, i^{+})) \quad \text{for } x/t < S \tag{5.2a}
$$

$$
\rho = \rho^{-}, \quad u = u^{-} = -u_0, \quad i = i^{-} = i_0, \quad (p = p^{-} = p_0 = p(\rho_0, i_0)) \quad \text{for } x/t > S, \tag{5.2b}
$$

where the shock moves out from the origin with speed $S$, and $S$, $\rho^{+}$, $i^{+}$, $p^{+} = p(\rho^{+}, i^{+})$ are given by the Rankine–Hugoniot shock relations. Thus

$$
S = \frac{[\rho u]}{[\rho]} = \frac{[p + \rho u^2]}{[\rho u]} = \frac{[u(e+p)]}{[e]}, \tag{5.3}
$$

where $[v] = v^{+} - v^{-}$ denotes the jump in $v$ across the shock. The solution of Eqs. (5.3) for $S$, $\rho^{+}$, $i^{+}$, $p^{+}$ subject to the initial conditions given by Eq. (5.1), and a precise form for the equation of state $p = p(\rho, i)$, is given by Glaister [9].

In the next section we give the numerical results obtained for the test problem considered here.

### 6. Numerical Results

In this section we show the numerical results obtained for the test problem given in Section 5 using the Riemann solver described in Section 3. Each of Figs. 2–10 refers to one of the equations of state given in Section 4 with different values of the parameters and initial conditions.

---

## [p. 398] — Figure 2 (plate; figure axes and caption transcribed, not the plotted curves themselves)

**Fig. 2.** Solution of the Euler equations with slab symmetry (shock reflection). Results for the ideal equation of state with ratio $p^{+}/p^{-} = \infty$, at time $t = 0.896$.

Key (applies to all four panels): $\rho$ – Density, $u$ – Velocity, $p$ – Pressure, $i$ – Internal energy. Solid line — Exact solution; circles (○○○○○) — Approximate solution.

Parameters:
- Ideal equation of state
- $\gamma = 5/3$
- 100 Mesh points
- 224 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0040$
- Pressure ratio $= \infty$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 1.000$
- $u = -1.000$
- $p = 0.000$
- $i = 0.000$

(Panels show $\rho$, $u$, $p$, $i$ each plotted against $x \in [0,1]$; axis tick values as printed: $\rho$ panel ticks $-3.99, -2.64, -1.35, 0.1, 1.35, 2.66, 3.99$; $u$ panel ticks $-0.99,-0.64,-0.35,0.1,0.35,0.64,0.99$; $p$ panel ticks $-1.32,-0.88,-0.44,0.1,0.44,0.88,1.32$; $i$ panel ticks $-0.46,-0.32,-0.16,0.1,0.16,0.32,0.48$.)

---

## [p. 399] — Figure 3

**Fig. 3.** Same as Fig. 2 with pressure ratio $p^{+}/p^{-} = 10$, at time $t = 0.578$.

Parameters:
- Ideal equation of state
- $\gamma = 5/3$
- 100 Mesh points
- 221 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0026$
- Pressure ratio $= 10$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 1.000$
- $u = -1.000$
- $p = 0.169$
- $i = 0.2531$

---

## [p. 400] — Figure 4

**Fig. 4.** Same as Fig. 2 with pressure ratio $p^{+}/p^{-} = 2$, at time $t = 0.150$.

Parameters:
- Ideal equation of state
- $\gamma = 5/3$
- 100 Mesh points
- 121 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0012$
- Pressure ratio $= 2$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 1.000$
- $u = -1.000$
- $p = 3.000$
- $i = 4.500$

---

## [p. 401] — Figure 5

**Fig. 5.** Solution of the Euler equations with slab symmetry (shock reflection). Results for the stiffened equation of state with $\gamma = 5/3$, $B = 1$, and the pressure ratio $p^{+}/p^{-} = \infty$, at time $t = 0.344$.

Parameters:
- Stiffened equation of state
- $\gamma = 5/3$, $B = 1.00$
- 100 Mesh points
- 172 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0020$
- Pressure ratio $= \infty$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 1.000$
- $u = -1.000$
- $p = 0.000$
- $i = 0.000$

---

## [p. 402] — Figure 6

**Fig. 6.** Same as Fig. 5 with pressure ratio $p^{+}/p^{-} = 10$, at time $t = 0.295$.

Parameters:
- Stiffened equation of state
- $\gamma = 5/3$, $B = 1.00$
- 100 Mesh points
- 160 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0018$
- Pressure ratio $= 10$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 1.000$
- $u = -1.000$
- $p = 0.224$
- $i = 0.3361$ [digit uncertain; printed as "0.336" with a trailing character partly cut — read as $0.3361$; possibly $0.336$]

---

## [p. 403] — Figure 7

**Fig. 7.** Same as Fig. 5 with pressure ratio $p^{+}/p^{-} = 2$, at time $t = 0.130$.

Parameters:
- Stiffened equation of state
- $\gamma = 5/3$, $B = 1.00$
- 100 Mesh points
- 115 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0011$
- Pressure ratio $= 2$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 1.000$
- $u = -1.000$
- $p = 3.303$
- $i = 4.954$)  [printed as "4.954)" — trailing parenthesis appears to be an OCR/plate artifact, not part of the number]

---

## [p. 404] — Figure 8

**Fig. 8.** Solution of the Euler equations with slab symmetry (shock reflection). Results for the general equation of state for copper with the pressure ratio $p^{+}/p^{-} = \infty$, at time $t = 0.311$.

Parameters:
- General equation of state for Copper due to R. K. Osborne
- 100 Mesh points
- 128 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0024$
- Pressure ratio $= \infty$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 8.900$
- $u = -1.000$
- $p = 0.000$
- $i = 0.000$

---

## [p. 405] — Figure 9

**Fig. 9.** Same as Fig. 8 with pressure ratio $p^{+}/p^{-} = 10$, at time $t = 0.288$.

Parameters:
- General equation of state for Copper due to R. K. Osborne
- 100 Mesh points
- 135 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0021$
- Pressure ratio $= 10$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 8.900$
- $u = -1.000$
- $p = 2.013$
- $i = 0.1371$ [printed value ambiguous between "0.1371" and "0.137)"; the trailing character is likely a stray plate mark, number read as $0.137$]

---

## [p. 406] — Figure 10

**Fig. 10.** Same as Fig. 8 with pressure ratio $p^{+}/p^{-} = 2$, at time $t = 0.139$.

Parameters:
- General equation of state for Copper due to R. K. Osborne
- 100 Mesh points
- 117 Time steps
- $\Delta x = 0.01$
- $\Delta t = 0.0012$
- Pressure ratio $= 2$
- "Superbee" limiter used

Initial Conditions (reflected boundary conditions at $x=0$):
- $\rho = 8.900$
- $u = -1.000$
- $p = 27.996$
- $i = 0.5471$ [printed "0.547)"; trailing character likely a stray plate mark, number read as $0.547$]

---

## [p. 407]

**(a) *Ideal equation of state.*** We take $\gamma = \tfrac{5}{3}$ with the initial data

$$\rho(x,0) = \rho_0 = 1$$
$$u(x,0) = -u_0 = -1$$

and choose $i(x,0) = i_0$ such that the pressure jump across the shock, i.e., $p^{+}/p^{-}$, takes the values $\infty$, 10, or 2.

**(b) *Stiffened equation of state.*** The parameters and initial data are taken to have the same values as for (a) and we choose $B = 1.0$. Three pressure ratios are obtained as for (a).

**(c) *General equation of state for copper.*** We consider the general equation of state given by Eq. (4.3) with values for the parameters corresponding to copper, i.e.,

$$
\rho_0 = 8.90, \qquad a_1 = 4.9578, \qquad a_2 = 3.6884,
$$
$$
b_0 = 7.4727, \qquad b_1 = 11.519, \qquad b_2 = 5.5251,
$$
$$
c_0 = 0.39493, \qquad c_1 = 0.52883, \qquad \phi_0 = 3.6000,
$$

together with the initial data

$$\rho(x,0) = \rho_0 = 8.9$$
$$u(x,0) = -u_0 = -1.$$

Again we choose $i(x,0) = i_0$ such that the pressure ratio $p^{+}/p^{-}$ takes the three values $\infty$, 10, or 2.

In each case we take 100 mesh points in $0 \le x \le 1$, and choose the output time so that the shock has moved a distance of 0.3. All computations have been done using a second order scheme with the "superbee limiter" (see [5]). We can see that the approximate solution gives a good representation of the exact solution, in particular, the correct shock speed has been achieved. The results obtained using the first-order algorithm only are not distinguishable from those given here.

Finally, we compare the c.p.u. time to compute the results obtained for the ideal gas case (a) using (i) Roe's original Riemann solver, and (ii) our general Riemann solver applied to the ideal gas case. The comparison, using an Amdahl V7, is as follows:

(i) Using "superbee" and 100 mesh points takes 0.0142 c.p.u. s to compute one time step, and a total of 1.6 c.p.u. s to reach a real time of 0.9 s using 112 time steps.

(ii) Using "superbee" and 100 mesh points takes 0.0178 c.p.u. s to compute one time step, and a total of 2.0 c.p.u. s to reach a real time of 0.9 s using 112 time steps.

This shows that our general Riemann solver is only slightly more expensive than Roe's original, as was to be expected. If we substitute the form of the ideal equation

*(page footer:)* 581/74/2-10

---

## [p. 408]

of state into Eqs. (3.55a)–(3.56b), however, and incorporate the resulting expressions into the finite difference code, we find that the two Riemann solvers are comparable in execution time.

### 7. Conclusions

We have extended the one-dimensional version of Roe's scheme to incorporate a general convex equation of state and have achieved satisfactory results for the shock reflection problem. In addition, we have seen that the algorithm is computationally efficient. This scheme can be extended to three dimensions incorporating operator splitting. Details of this extension together with a two-dimensional calculation of the flow in a tunnel containing a step involving interacting waves are given by Glaister [10].

There may be scope for improving the efficiency of our scheme using the ideas of Colella and Glaz [11] on efficient solution algorithms for the Riemann problem for real gases and the work of Harten [12] on the symmetrisation of systems of conservation laws which possess entropy functions.

### Acknowledgments

I would like to express my thanks to Dr. M. J. Baines for useful discussions and to D. L. Youngs for suggesting the test problem given in Section 5. I acknowledge the financial support of A.W.R.E., Aldermaston.

### References

1. P. L. Roe, *J. Comput. Phys.* **43**, 357 (1981).
2. P. L. Roe, Cranfield Institute of Technology, Cranfield, U. K., private communication (1985).
3. R. G. Smith, *Trans. Amer. Math. Soc.* **249**, 1 (1979).
4. P. L. Roe and J. Pike, "Efficient Construction and Utilisation of Approximate Riemann Solutions," in *Computing Methods in Applied Science and Engineering VI*, edited by R. Glowinski and J.-L. Lions (North-Holland, Amsterdam, 1984), p. 499.
5. P. K. Sweby, *SIAM J. Numer. Anal.* **21**, 995 (1984).
6. P. K. Sweby, University of Reading Numerical Analysis Report 6–82, 1982 (unpublished).
7. T. D. Riney, "Numerical Evaluation of Hypervelocity Impact Phenomena," in *High-Velocity Impact Phenomena*, edited by R. Kinslow (Academic Press, New York/London, 1970), p. 164.
8. K. P. Stanyukovich, *Unsteady Motion of Continuous Media* (Pergamon Press, London/Oxford/Paris/New York, 1960), p. 221.
9. P. Glaister, University of Reading Numerical Analysis Report 7–86, 1986 (unpublished).
10. P. Glaister, University of Reading Numerical Analysis Report 11–86, 1986 (unpublished).
11. P. Colella and H. M. Glaz, *J. Comput. Phys.* **59**, 264 (1985).
12. A. Harten, *J. Comput. Phys.* **49**, 151 (1983).
