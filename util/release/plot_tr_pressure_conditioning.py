"""Which way the analytic transition region (TR) can be read without amplifying errors.

Rebuilds the static (no-flow) analytic TR of SWMF's ModTransitionRegion.f90
(check_tr_table, uCh = 0) with the same CHIANTI 8.0 loss table, the same
conduction coefficient kappa_0_e(20) and the same trapezoid-in-ln(T) scheme,
but lets the base temperature T_b and the base heat flux q_b be free inputs
(SWMF fixes T_b = 10 kK and q_b = 0 when uCh = 0). Fully ionized hydrogen,
Z = 1, isobaric, steady, gravity-free:

    d(q^2)/dT = 2 kappa_0 Pe^2 Lambda(T) T^{1/2} / k_B^2,    q(T_b) = q_b,
    L         = int_{T_b}^{T_c} kappa_0 T^{5/2} / q dT,

with Pe = n k_B T the electron pressure (SWMF's Pavr for Z = 1); the total
pressure reported here is p = 2 Pe. q is the conductive heat flux magnitude,
directed down the TR; q_b is the part that leaves the TR into the
chromosphere.

For fixed (T_c, L) the TR fixes the base pressure p as a function of the
coronal input T_c and of the chromospheric inputs q_b and T_b. The script
reports the logarithmic sensitivities S_x = d ln p / d ln x. Reading the
relation forwards (inputs -> p) multiplies an input error by |S_x|; reading
it backwards (p -> x) multiplies a pressure error by 1/|S_x|.

The reference point is the release top boundary: T_b = 22 kK, p = 0.0103 Pa,
q_b = 0.527 W m^-2 (the final q_total of the release run's .outercond
diagnostic, outputs/model_column/godunov_release/prod_godunov.txt.outercond),
L = 1.5 Mm; T_c is solved so that these are consistent.

Usage (from the project root):
    MPLCONFIGDIR=/tmp/chromosphere2026-mpl .venv/bin/python \
        util/release/plot_tr_pressure_conditioning.py \
        --out visualization/model_column/tr_pressure_conditioning.pdf
"""

import argparse
import os

import numpy as np
from scipy.optimize import brentq

KB = 1.380649e-23
ME = 9.1093837015e-31
QE = 1.602176634e-19
EPS0 = 8.8541878128e-12

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_TABLE = os.path.join(os.path.dirname(REPO), "SWMF", "GM", "BATSRUS",
                             "Param", "CORONA", "RadCoolCorona_8.0.dat")

# SWMF check_tr_table grid: 310 log-uniform points over [1e4, 2.9988e6] K
DLNT = np.log(2.9988442312309672e6 / 1.0e4) / 309

# Release reference point
T_B_REF = 22.0e3        # [K]
P_REF = 0.0103          # total pressure [Pa]
Q_B_REF = 0.527         # heat flux into the column top [W m^-2]
L_REF = 1.5e6           # TR length [m]


def kappa_0_e(coulomb_log=20.0):
    """SWMF share/Library/src/ModConst.f90 kappa_0_e [W m^-1 K^-7/2]."""
    return (3.2 * 3.0 * 2 * np.pi / coulomb_log
            * np.sqrt(2 * np.pi * KB / ME) * KB
            * ((EPS0 / QE) * (KB / QE)) ** 2)


KAPPA0 = kappa_0_e()


def load_lambda(path):
    """CHIANTI loss table as SWMF reads it: log10 T, Lambda * 1e22 [erg cm^3 s^-1]."""
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[2].split()[0])
    data = np.array([[float(v) for v in ln.split()[:2]] for ln in lines[5:5 + n]])
    log_t, lam = data[:, 0], data[:, 1] * 1.0e-35  # -> [W m^3]
    return lambda t: np.interp(np.log10(t), log_t, lam)


class StaticTR:
    def __init__(self, lam):
        self.lam = lam

    def integrate(self, pe, t_b, t_c, q_b):
        """Return (L [m], q_c [W m^-2]) for electron pressure pe [Pa]."""
        n = max(int(np.ceil(np.log(t_c / t_b) / DLNT)), 8)
        t = t_b * np.exp(np.linspace(0.0, np.log(t_c / t_b), n + 1))
        d = np.log(t_c / t_b) / n
        src = KAPPA0 * pe**2 / KB**2 * self.lam(t) * t**1.5 * d
        q2 = np.empty_like(t)
        q2[0] = q_b**2
        length = 0.0
        for i in range(1, len(t)):
            semi = np.sqrt(q2[i - 1] + src[i - 1])        # SWMF SemiIntUheat
            q2[i] = q2[i - 1] + src[i - 1] + src[i]
            length += 0.5 * d * KAPPA0 * (t[i - 1]**3.5 + t[i]**3.5) / semi
        return length, np.sqrt(q2[-1])

    def pressure(self, t_b, t_c, q_b, length):
        """Total base pressure p = 2 Pe [Pa] of the TR with the given length."""
        f = lambda lpe: np.log(self.integrate(np.exp(lpe), t_b, t_c, q_b)[0] / length)
        return 2.0 * np.exp(brentq(f, np.log(1e-8), np.log(1e2), xtol=1e-12))

    def q_base(self, p, t_b, t_c, q_c):
        """Base heat flux q_b for given total pressure and top flux (no flow)."""
        _, q_c0 = self.integrate(p / 2.0, t_b, t_c, 0.0)
        return np.sqrt(max(q_c**2 - q_c0**2, 0.0)), q_c0**2 > q_c**2


def log_slope(fun, x, h=1e-4):
    return (np.log(fun(x * (1 + h))) - np.log(fun(x * (1 - h)))) / (np.log(1 + h) - np.log(1 - h))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--table", default=DEFAULT_TABLE, help="SWMF RadCoolCorona_8.0.dat")
    ap.add_argument("--out", default="visualization/model_column/tr_pressure_conditioning.pdf")
    args = ap.parse_args()

    tr = StaticTR(load_lambda(args.table))
    t_c = brentq(lambda tc: tr.pressure(T_B_REF, tc, Q_B_REF, L_REF) / P_REF - 1.0, 1.0e5, 2.9e6)
    _, q_c = tr.integrate(P_REF / 2, T_B_REF, t_c, Q_B_REF)
    print(f"kappa_0_e(20) = {KAPPA0:.4e} W m^-1 K^-7/2")
    print(f"reference: T_b = {T_B_REF:.0f} K, p = {P_REF} Pa, q_b = {Q_B_REF} W/m2, "
          f"L = {L_REF/1e6:.2f} Mm  ->  T_c = {t_c/1e3:.1f} kK, q_c = {q_c:.2f} W/m2, "
          f"(q_c/q_b)^2 = {(q_c/Q_B_REF)**2:.3g}")

    p_of_tc = lambda x: tr.pressure(T_B_REF, x, Q_B_REF, L_REF)
    p_of_qb = lambda x: tr.pressure(T_B_REF, t_c, x, L_REF)
    p_of_tb = lambda x: tr.pressure(x, t_c, Q_B_REF, L_REF)
    s_tc, s_qb, s_tb = (log_slope(p_of_tc, t_c), log_slope(p_of_qb, Q_B_REF),
                        log_slope(p_of_tb, T_B_REF))
    print("\nForward sensitivities S = dln p / dln x at the reference point (fixed L):")
    for name, s in (("T_c (corona)", s_tc), ("q_b (chromosphere heat)", s_qb),
                    ("T_b (base temperature)", s_tb)):
        print(f"  {name:26s} S = {s:+.3e}   backwards: 1% in p -> {100/abs(s):.3g}% in x")

    # Same question with the coronal flux q_c held instead of L
    # q_b^2 = q_c^2 - q_c0(p)^2 with q_c0^2 = q_c^2 - q_b^2 at p_ref and q_c0^2 ∝ p^2,
    # so dln q_b / dln p = -(q_c^2 - q_b^2) / q_b^2 exactly
    s_qb_qc = -(q_c**2 - Q_B_REF**2) / Q_B_REF**2
    print(f"\nHolding q_c instead of L: dln q_b / dln p = {s_qb_qc:+.3e} "
          f"(0.01% in p -> {0.01*abs(s_qb_qc):.3g}% in q_b)")
    p_kill = brentq(lambda pp: tr.integrate(pp / 2, T_B_REF, t_c, 0.0)[1] - q_c, P_REF, 2 * P_REF)
    print(f"  q_b reaches 0 (no TR solution beyond) at p = {p_kill:.6g} Pa, "
          f"i.e. {100*(p_kill/P_REF - 1):+.4f}% above the reference")

    # Range checks quoted in the discussion
    for qb in (0.0, 2 * Q_B_REF, 10 * Q_B_REF):
        print(f"  p(q_b = {qb:5.2f}) = {p_of_qb(qb):.6g} Pa ({100*(p_of_qb(qb)/P_REF-1):+.4f}%)")
    for tb in (10e3, 15e3, 30e3, 50e3):
        print(f"  p(T_b = {tb/1e3:4.0f} kK) = {p_of_tb(tb):.6g} Pa ({100*(p_of_tb(tb)/P_REF-1):+.3f}%)")
    for tc in (3e5, 1e6):
        print(f"  p(T_c = {tc/1e3:5.0f} kK) = {p_of_tc(tc):.4g} Pa")

    # ------------------------------------------------------------------ figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    qb_grid = np.linspace(0.0, 10 * Q_B_REF, 30)
    tb_grid = np.linspace(10e3, 50e3, 30)
    p_qb = np.array([p_of_qb(x) for x in qb_grid])
    p_tb = np.array([p_of_tb(x) for x in tb_grid])

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 3.8),
                                        gridspec_kw={"width_ratios": [1.2, 1, 1], "wspace": 0.4})

    # 1. Energy budget of the TR
    rad = q_c - Q_B_REF
    labels = ["heat coming down\nfrom the corona", "radiated away\ninside the TR",
              "left over, flows into\nthe chromosphere"]
    vals = [q_c, rad, Q_B_REF]
    ax1.barh(labels[::-1], vals[::-1], color=["C3", "C1", "C0"][::-1])
    for i, v in enumerate(vals[::-1]):
        ax1.text(v + 3, i, f"{v:.1f} W/m$^2$", va="center")
    ax1.set_xlim(0, 210)
    ax1.set_xlabel("heat flux [W/m$^2$]")
    ax1.set_title("Where the coronal heat goes", loc="left")

    # 2 and 3. Pressure against the chromosphere-side numbers
    for ax, x, y, xref, xlabel, title, note in (
            (ax2, qb_grid, p_qb, Q_B_REF, "heat flowing into the chromosphere [W/m$^2$]",
             "Pressure vs. heat into chromosphere",
             f"10x more heat:\npressure changes {100*(p_qb[-1]/p_qb[0]-1):+.1f} %"),
            (ax3, tb_grid / 1e3, p_tb, T_B_REF / 1e3, "temperature at the TR base [kK]",
             "Pressure vs. TR base temperature",
             f"10 kK to 50 kK:\npressure changes {100*(p_tb[-1]/p_tb[0]-1):+.1f} %")):
        ax.plot(x, y, color="C0", lw=2.5)
        ax.plot([xref], [P_REF], "o", color="k", zorder=5)
        ax.set_ylim(0, 0.02)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("pressure at the TR base [Pa]")
        ax.set_title(title, loc="left")
        ax.text(0.5, 0.25, note, transform=ax.transAxes, ha="center", fontsize=11)
    ax2.annotate("release", (Q_B_REF, P_REF), xytext=(1.2, 0.0135),
                 arrowprops=dict(arrowstyle="->", color="k"))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
