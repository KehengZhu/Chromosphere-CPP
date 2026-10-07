"""Steady transition-region length L_TR(T) from the release closures.

Evaluates the derivation in docs/studies/boundaries/tr_length_steady_derivation.md:
steady, constant-area, gravity-free flow with rho u = F, p + rho u^2 = const,
energy d/ds[F(h + u^2/2) + q] = -Q, q = -kappa dT/ds, with the release
pure-hydrogen Saha EOS and the release Spitzer-plus-neutral kappa.

In temperature space the conductive flux q (< 0 in the TR) obeys
    dq/dT = -F dH/dT + kappa Q / q,   q(T_c) = 0,
and L(T) = -int kappa / q dT. With Q = 0 this integrates to q = -F (H - H_c).
The code integrates Phi = -q > 0 to keep the ODE state positive.

Checks printed:
  1. Q = 0: the ODE route reproduces the closed quadrature.
  2. Q = 0: the kinetic + momentum correction to L and its fully ionized
     prediction -5 M_top^2 / 21.
  3. F = 0, fully ionized power-law loss: the ODE reproduces the closed-form
     static length (solver check only; the release has no radiation).
  4. F > 0 with the same loss: L lies inside the closed-form bracket
     int kappa/(|q_e| + |q_s|) <= L <= int kappa/sqrt(q_e^2 + q_s^2).
  5. F << F* with the same loss: L = L_s (1 - C_alpha F/F*), and F recovered
     from the top flux as F = |q_s| (|q_top| - |q_s|) / I.

Usage: .venv/bin/python util/tr_length_steady.py
"""

import numpy as np
from scipy.integrate import solve_ivp

from conductive_flux import physical_conductivity

KB = 1.380649e-23
ME = 9.1093837015e-31
MH = 1.6726219e-27
HP = 6.62607015e-34
CHI = 2.179872361e-18
KAPPA0 = 9.2048e-12            # release Spitzer coefficient [W m^-1 K^-7/2]


def saha_x(T, p):
    """Pure-H Saha x(T, p): x^2/(1 - x^2) = S(T) k_B T / p."""
    log_s = 1.5 * np.log(2.0 * np.pi * ME * KB * T / HP**2) - CHI / (KB * T)
    a = np.exp(log_s) * KB * T / p
    return 1.0 / np.sqrt(1.0 + 1.0 / a)


def rho_of(T, p):
    return MH * p / ((1.0 + saha_x(T, p)) * KB * T)


def h_of(T, p):
    x = saha_x(T, p)
    return (2.5 * (1.0 + x) * KB * T + x * CHI) / MH


def kappa_of(T, p):
    x = saha_x(T, p)
    n_h = p / ((1.0 + x) * KB * T)
    return physical_conductivity(x * n_h, (1.0 - x) * n_h, T)


def pressure_along(T, F, Pi, iters=60):
    """Subsonic root of p + F^2/rho(T, p) = Pi."""
    p = np.full_like(np.asarray(T, dtype=float), Pi)
    for _ in range(iters):
        p = Pi - F * F / rho_of(T, p)
    return p


def profile(T, F, T_top, p_top, T_c, kinetic=True):
    """(p, H, H_c) along T. kinetic=False is the isobaric, u^2-free limit."""
    if not kinetic:
        return np.full_like(T, p_top), h_of(T, p_top), h_of(T_c, p_top)
    Pi = p_top + F * F / rho_of(T_top, p_top)
    p = pressure_along(T, F, Pi)
    pc = pressure_along(np.array([T_c]), F, Pi)[0]
    H = h_of(T, p) + 0.5 * (F / rho_of(T, p))**2
    H_c = h_of(T_c, pc) + 0.5 * (F / rho_of(T_c, pc))**2
    return p, H, H_c


def length_quadrature(F, T_top, p_top, T_c, T_1, kinetic=True, n=4000):
    """Q = 0: L = (1/F) int_{T_1}^{T_top} kappa / (H - H_c) dT [m]."""
    T = np.geomspace(T_1, T_top, n)
    p, H, H_c = profile(T, F, T_top, p_top, T_c, kinetic)
    return np.trapezoid(kappa_of(T, p) / (F * (H - H_c)), T)


def length_ode(F, T_top, p_top, T_c, T_1, Q=None, kappa=None, H_of=None,
               with_flux=False):
    """General route: integrate Phi(T) and L(T) together from T_c to T_top.

    Q(T, p) is the volumetric loss [W m^-3], zero below T_c. Returns L from T_1,
    or (L, Phi(T_top)) with with_flux=True.
    kappa / H_of override the release closures (used for the analytic check).
    """
    kappa = kappa or kappa_of
    Pi = p_top + (F * F / rho_of(T_top, p_top) if F > 0 else 0.0)

    def H(T):
        p = pressure_along(np.array([T]), F, Pi)[0] if F > 0 else p_top
        if H_of is not None:
            return H_of(T), p
        return h_of(T, p) + 0.5 * (F / rho_of(T, p))**2, p

    def rhs(T, y):
        Phi = y[0]
        dT = 1.0e-4 * T
        dH = (H(T + dT)[0] - H(T - dT)[0]) / (2.0 * dT)
        p = H(T)[1]
        k = kappa(T, p)
        q_loss = Q(T, p) if Q is not None else 0.0
        return [F * dH + k * q_loss / Phi, k / Phi if T >= T_1 else 0.0]

    # Start just above T_c on the local analytic solution of the Phi equation.
    eps = 1.0e-6 * T_c
    k_c = kappa(T_c, p_top)
    Q_c = Q(T_c + eps, p_top) if Q is not None else 0.0
    dH_c = (H(T_c + eps)[0] - H(T_c)[0]) / eps
    Phi0 = F * dH_c * eps + np.sqrt(2.0 * k_c * Q_c * eps)
    sol = solve_ivp(rhs, (T_c + eps, T_top), [Phi0, 0.0], method="LSODA",
                    rtol=1e-10, atol=[1e-30, 1e-12])
    if with_flux:
        return sol.y[1, -1], sol.y[0, -1]
    return sol.y[1, -1]


def main():
    T_top, p_top, T_c, T_1 = 22.0e3, 0.0103, 6700.0, 7500.0
    print(f"release top state: T_top={T_top:.0f} K  p_top={p_top} Pa  "
          f"T_c={T_c:.0f} K  marker T_1={T_1:.0f} K")

    print("\n1. Q = 0: ODE route vs quadrature")
    F = 2.683e-10
    Lq = length_quadrature(F, T_top, p_top, T_c, T_1)
    Lo = length_ode(F, T_top, p_top, T_c, T_1)
    print(f"   F={F:.3e}  L_quad={Lq/1e3:.4f} km  L_ode={Lo/1e3:.4f} km  "
          f"rel={Lo/Lq-1:+.2e}")

    print("\n2. Q = 0: effect of keeping u^2 (and the momentum pressure drop)")
    print("   T_top[kK]       F[kg/m2/s]   u_top[m/s]    M_top^2     "
          "L_full[km]   dL/L        -5M^2/21")
    for T_top_i in (22.0e3, 100.0e3, 1.0e6):
        for F in (2.683e-10, 1.0e-8, 1.0e-7):
            rho_t = rho_of(T_top_i, p_top)
            u_t = F / rho_t
            m2 = u_t**2 / ((5.0 / 3.0) * p_top / rho_t)
            Lf = length_quadrature(F, T_top_i, p_top, T_c, T_1, kinetic=True)
            L0 = length_quadrature(F, T_top_i, p_top, T_c, T_1, kinetic=False)
            print(f"   {T_top_i/1e3:9.0f}  {F:12.3e}  {u_t:11.3e}  {m2:10.3e}  "
                  f"{Lf/1e3:11.4e}  {Lf/L0-1:+.3e}  {-5*m2/21:+.3e}")

    print("\n3. Static (F = 0) fully ionized check, Q = n_e n_H Lambda0 T^alpha")
    lam0, alpha, Tc_s, Tt_s = 1.0e-32, -0.5, 2.0e4, 1.0e6
    kap = lambda T, p: KAPPA0 * T**2.5
    Q = lambda T, p: (p / (2 * KB * T))**2 * lam0 * T**alpha
    Lo = length_ode(0.0, Tt_s, p_top, Tc_s, Tc_s, Q=Q, kappa=kap,
                    H_of=lambda T: 5 * KB * T / MH)
    a = alpha + 1.5
    c = KAPPA0 * lam0 * p_top**2 / (2 * KB**2 * a)
    # w^2 = T^a - T_c^a removes the integrable 1/sqrt singularity at T_c.
    w = np.linspace(0.0, np.sqrt(Tt_s**a - Tc_s**a), 200001)
    Tg = (Tc_s**a + w * w)**(1.0 / a)
    Lc = np.trapezoid(2 * KAPPA0 * Tg**(3.5 - a) / (a * np.sqrt(c)), w)
    print(f"   L_ode={Lo/1e3:.5e} km  L_closed={Lc/1e3:.5e} km  "
          f"rel={Lo/Lc-1:+.2e}")

    print("\n4. Same loss with evaporation: L inside the bracket")
    for F in (1.0e-8, 3.0e-8, 1.0e-7):
        Lo = length_ode(F, Tt_s, p_top, Tc_s, Tc_s, Q=Q, kappa=kap,
                        H_of=lambda T: 5 * KB * T / MH)
        Tg = np.geomspace(Tc_s * (1 + 1e-9), Tt_s, 400001)
        phi_e = F * 5 * KB * (Tg - Tc_s) / MH
        phi_s = np.sqrt(c * (Tg**a - Tc_s**a))
        lo = np.trapezoid(kap(Tg, p_top) / (phi_e + phi_s), Tg)
        hi = np.trapezoid(kap(Tg, p_top) / np.hypot(phi_e, phi_s), Tg)
        print(f"   F={F:.0e}  lower={lo/1e3:.4e}  L_ode={Lo/1e3:.4e}  "
              f"upper={hi/1e3:.4e} km  static={Lc/1e3:.4e} km")

    print("\n5. Same loss, F << F*: first-order F correction to the static TR")
    L_s, q_s = length_ode(0.0, Tt_s, p_top, Tc_s, Tc_s, Q=Q, kappa=kap,
                          H_of=lambda T: 5 * KB * T / MH, with_flux=True)
    dH_top = 5 * KB * (Tt_s - Tc_s) / MH
    F_star = q_s / dH_top
    # I = int |q_s| dH from T_c to T_top, with |q_s| = sqrt(c (T^a - T_c^a)).
    Tg = np.geomspace(Tc_s, Tt_s, 400001)
    I_top = np.trapezoid(np.sqrt(c * (Tg**a - Tc_s**a)) * 5 * KB / MH, Tg)
    c_alpha = (11 - 2 * alpha) / ((2 * alpha + 7) * (3 - alpha))
    print(f"   L_s={L_s/1e3:.4e} km  |q_s(T_top)|={q_s:.4e} W/m2  "
          f"F*={F_star:.4e}  C_alpha={c_alpha:.4f}")
    print("        F       F/F*     dL/L(ODE)   -C F/F*     F from q_top  rel")
    for F in (2.683e-10, 1.0e-10, 3.0e-11, 1.0e-11):
        L, q_top = length_ode(F, Tt_s, p_top, Tc_s, Tc_s, Q=Q, kappa=kap,
                              H_of=lambda T: 5 * KB * T / MH, with_flux=True)
        F_rec = q_s * (q_top - q_s) / I_top
        print(f"   {F:9.3e}  {F/F_star:.3e}  {L/L_s-1:+.4e}  "
              f"{-c_alpha*F/F_star:+.4e}  {F_rec:.4e}  {F_rec/F-1:+.2e}")


if __name__ == "__main__":
    main()
