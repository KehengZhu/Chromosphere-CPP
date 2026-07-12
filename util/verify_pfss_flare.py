#!/usr/bin/env python3
"""Phase-3 de-risk check: did beam heating on a PFSS field line drive explosive
chromospheric evaporation?

Compares a control run (beam off) against an explosive run (beam on) produced by

    PFSS_FLARE=1 FLARE_BEAM_FLUX=0   chromo_main ... pfss_field_line line.dat ... -> control
    PFSS_FLARE=1 FLARE_BEAM_FLUX=5e7 chromo_main ... pfss_field_line line.dat ... -> explosive

Prints peak diagnostics (T, upflow, condensation downflow) and saves a
height-time heatmap of T_i and V for the explosive run with the control beside.

Temperature uses the FULL cons2prim (state.cpp): p_i = 2/3 E_i - 1/3 rho_i V^2
- 2/3 rho_i phi_g, with phi_g = g (h - h_base) reconstructed from the cell-center
height. (plot_output.py omits the phi_g term, which inflates T near the top.)
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m_i = 1.6726219e-27
m_n = m_i
k_b = 1.380649e-23
g_si = 274.0
H_BASE_KM = 1003.0

CNI, CNN, CNV, CNU, CEI, CEN = 0, 1, 2, 3, 4, 5


def read_frames(path):
    with open(path) as f:
        raw = [ln for ln in f.read().splitlines() if ln.strip()]
    ns, neq = (int(x) for x in raw[0].split())
    heights = np.fromstring(raw[1], sep=" ")  # cumulative top-face height [km]
    frames, i = [], 2
    while i < len(raw):
        toks = raw[i].split()
        t, stp = float(toks[3]), int(toks[6])
        i += 1
        data = np.zeros((ns, neq))
        for j in range(ns):
            data[j, :] = np.fromstring(raw[i + j], sep=" ")
        i += ns
        frames.append((t, stp, data))
    return heights, frames


def cell_center_km(heights):
    prev = np.concatenate(([H_BASE_KM], heights[:-1]))
    return 0.5 * (heights + prev)


def primitives(xn, phi_g):
    ni = xn[:, CNI] / m_i
    nn = xn[:, CNN] / m_n
    v = xn[:, CNV] / xn[:, CNI]
    u = xn[:, CNU] / xn[:, CNN]
    p_i = 2.0 / 3.0 * xn[:, CEI] - 1.0 / 3.0 * xn[:, CNI] * v * v - 2.0 / 3.0 * xn[:, CNI] * phi_g
    p_n = 2.0 / 3.0 * xn[:, CEN] - 1.0 / 3.0 * xn[:, CNN] * u * u - 2.0 / 3.0 * xn[:, CNN] * phi_g
    Ti = p_i / (2.0 * ni * k_b)
    Tn = p_n / (nn * k_b)
    return ni, nn, v, u, Ti, Tn


def stack(frames, phi_g):
    """Return (times, T_i[nt,ns], V[nt,ns]) arrays over all frames."""
    ts, Ti_all, V_all = [], [], []
    for t, _stp, xn in frames:
        _ni, _nn, v, _u, Ti, _Tn = primitives(xn, phi_g)
        ts.append(t)
        Ti_all.append(Ti)
        V_all.append(v)
    return np.array(ts), np.array(Ti_all), np.array(V_all)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ctrl = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "outputs", "pfss_lineA_control.txt")
    expl = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "outputs", "pfss_lineA_flare.txt")
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(here, "..", "visualization", "pfss_lineA_flare_verify.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    h, f_c = read_frames(ctrl)
    _h, f_e = read_frames(expl)
    phi_g = g_si * (cell_center_km(h) - H_BASE_KM) * 1.0e3  # [J/kg]

    tc, Tic, Vc = stack(f_c, phi_g)
    te, Tie, Ve = stack(f_e, phi_g)

    def report(tag, t, Ti, V):
        it, ix = np.unravel_index(np.nanargmax(Ti), Ti.shape)
        iu = np.nanargmax(V)
        iud = np.unravel_index(iu, V.shape)
        idn = np.unravel_index(np.nanargmin(V), V.shape)
        print(f"[{tag}]")
        print(f"  peak T_i      = {Ti[it, ix]/1e6:8.3f} MK  at t={t[it]:6.2f}s, h={h[ix]:7.1f} km")
        print(f"  peak upflow   = {V[iud]/1e3:8.2f} km/s  at t={t[iud[0]]:6.2f}s, h={h[iud[1]]:7.1f} km")
        print(f"  peak downflow = {V[idn]/1e3:8.2f} km/s  at t={t[idn[0]]:6.2f}s, h={h[idn[1]]:7.1f} km")

    print("=" * 64)
    report("CONTROL  (beam off)", tc, Tic, Vc)
    report("EXPLOSIVE(beam 5e7)", te, Tie, Ve)
    print("=" * 64)

    # Height-time heatmaps: T_i and V for both runs.
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex="col", sharey=True)
    fig.suptitle("PFSS field line (pfss_lineA): control vs explosive beam — Phase-3 de-risk", fontsize=13)

    def hm(ax, t, Z, title, cmap, vmin=None, vmax=None, log=False):
        from matplotlib.colors import LogNorm
        norm = LogNorm(vmin=vmin, vmax=vmax) if log else None
        im = ax.pcolormesh(t, h, Z.T, shading="auto", cmap=cmap,
                           norm=norm, vmin=None if log else vmin, vmax=None if log else vmax)
        ax.set_title(title)
        ax.set_ylabel("height [km]")
        ax.set_xlabel("t [s]")
        fig.colorbar(im, ax=ax)

    tmax_e = te.max()
    hm(axes[0, 0], tc, np.clip(Tic, 3e3, None), "control  T_i [K]", "inferno", vmin=3e3, vmax=2e7, log=True)
    hm(axes[0, 1], te, np.clip(Tie, 3e3, None), "explosive T_i [K]", "inferno", vmin=3e3, vmax=2e7, log=True)
    vlim = max(abs(np.nanpercentile(Ve, 1)), abs(np.nanpercentile(Ve, 99))) / 1e3
    axes[1, 0].pcolormesh(tc, h, (Vc / 1e3).T, shading="auto", cmap="RdBu_r", vmin=-vlim, vmax=vlim)
    axes[1, 0].set_title("control  V [km/s]"); axes[1, 0].set_xlabel("t [s]"); axes[1, 0].set_ylabel("height [km]")
    im = axes[1, 1].pcolormesh(te, h, (Ve / 1e3).T, shading="auto", cmap="RdBu_r", vmin=-vlim, vmax=vlim)
    axes[1, 1].set_title("explosive V [km/s] (red=upflow)"); axes[1, 1].set_xlabel("t [s]")
    fig.colorbar(im, ax=axes[1, :].tolist())
    for ax in axes[:, 1]:
        ax.set_xlim(0, tmax_e)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
