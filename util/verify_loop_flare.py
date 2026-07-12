#!/usr/bin/env python3
"""Verify the full-loop flare run vs the truncated-line artifact.

Unlike verify_pfss_flare.py, height ≠ arc length on a curved loop, so we read the
true per-cell phi_g (→ height h = phi_g/g) from the .dat rather than the output's
cumulative-arc-length column. Plots T and V vs TRUE HEIGHT for control vs flare,
and the coronal mass vs time (does the evaporated plasma fill the corona?).

Usage: python util/verify_loop_flare.py <loop.dat> <control.txt> <flare.txt>
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m_i = 1.6726219e-27
k_b = 1.380649e-23
g_si = 274.0


def parse_dat(path):
    rows, in_cells = [], False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("[CELLS]"):
            in_cells = True; continue
        if s.startswith("[GHOSTS]"):
            break
        if not in_cells or not s or s.startswith("#"):
            continue
        rows.append([float(x) for x in s.split()])
    a = np.array(rows)
    ds = a[:, 1]
    phi_c = 0.5 * (a[:, 4] + a[:, 5])     # cell-center potential [J/kg]
    return ds, phi_c, phi_c / g_si        # ds, phi_g, height[m]


def read_frames(path):
    raw = [l for l in open(path).read().splitlines() if l.strip()]
    ns, neq = map(int, raw[0].split())
    i, frames = 2, []
    while i < len(raw):
        t = float(raw[i].split()[3]); i += 1
        d = np.zeros((ns, neq)); ok = True
        for j in range(ns):
            tk = raw[i + j].split()
            if len(tk) != neq:
                ok = False; break
            d[j] = [float(x) for x in tk]
        i += ns
        if ok:
            frames.append((t, d))
    return frames


def stack(frames, phi_g):
    ts, Ti, V, rho = [], [], [], []
    for t, xn in frames:
        ni = xn[:, 0] / m_i
        v = xn[:, 2] / xn[:, 0]
        p_i = 2/3*xn[:, 4] - 1/3*xn[:, 0]*v*v - 2/3*xn[:, 0]*phi_g
        ts.append(t); Ti.append(p_i/(2*ni*k_b)); V.append(v)
        rho.append(xn[:, 0] + xn[:, 1])
    return np.array(ts), np.array(Ti), np.array(V), np.array(rho)


def main():
    dat = sys.argv[1] if len(sys.argv) > 1 else "scenarios/data/loop_closed_test.dat"
    ctrl = sys.argv[2] if len(sys.argv) > 2 else "outputs/loop_closed_quiet.txt"
    expl = sys.argv[3] if len(sys.argv) > 3 else "outputs/loop_closed_flare.txt"
    out = "visualization/loop_flare_verify.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)

    ds, phi_g, h = parse_dat(dat)
    hMm = h / 1e6
    tc, Tic, Vc, Rc = stack(read_frames(ctrl), phi_g)
    te, Tie, Ve, Re = stack(read_frames(expl), phi_g)

    def rep(tag, t, Ti, V):
        it = np.unravel_index(np.nanargmax(Ti), Ti.shape)
        iu = np.unravel_index(np.nanargmax(V), V.shape)
        idn = np.unravel_index(np.nanargmin(V), V.shape)
        print(f"[{tag}] peak T={Ti[it]/1e6:7.2f} MK @ t={t[it[0]]:.1f}s h={hMm[it[1]]:.2f}Mm | "
              f"up={V[iu]/1e3:+7.1f} | down={V[idn]/1e3:+7.1f} km/s | t_end={t[-1]:.1f}s")

    print("=" * 70)
    rep("CONTROL ", tc, Tic, Vc)
    rep("FLARE   ", te, Tie, Ve)
    # coronal mass (cells above 3 Mm) vs time, normalized to t=0
    cor = hMm > 3.0
    mc = Rc[:, cor] @ ds[cor]; me = Re[:, cor] @ ds[cor]
    print(f"coronal mass (h>3Mm) flare: t=0 {me[0]:.2e} -> peak {me.max():.2e} "
          f"({me.max()/me[0]:.1f}x) @ t={te[np.argmax(me)]:.1f}s")
    print("=" * 70)

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("Closed loop (25 Mm half-length): control vs flare — does the corona fill?", fontsize=13)
    from matplotlib.colors import LogNorm

    def hm(a, t, Z, title, cmap, **kw):
        im = a.pcolormesh(t, hMm, Z.T, shading="auto", cmap=cmap, **kw)
        a.set_title(title); a.set_ylabel("true height [Mm]"); a.set_xlabel("t [s]")
        fig.colorbar(im, ax=a)

    hm(ax[0, 0], te, np.clip(Tie, 3e3, 6e7), "flare  T_i [K]", "inferno",
       norm=LogNorm(vmin=3e3, vmax=6e7))
    vlim = np.nanpercentile(np.abs(Ve), 99) / 1e3
    hm(ax[0, 1], te, Ve/1e3, "flare  V [km/s] (red=up the loop)", "RdBu_r", vmin=-vlim, vmax=vlim)
    # coronal mass vs time
    ax[1, 0].plot(te, me/me[0], "C3", label="flare")
    ax[1, 0].plot(tc, mc/mc[0], "0.5", label="control")
    ax[1, 0].set_xlabel("t [s]"); ax[1, 0].set_ylabel("coronal mass (h>3 Mm) / initial")
    ax[1, 0].grid(alpha=0.3); ax[1, 0].legend(); ax[1, 0].set_title("coronal mass loading")
    # final velocity profile (is there a uniform 200 km/s drift?)
    ax[1, 1].plot(hMm, Ve[-1]/1e3, "C0", label=f"flare t={te[-1]:.0f}s")
    ax[1, 1].axhline(0, color="k", lw=0.6)
    ax[1, 1].set_xlabel("true height [Mm]"); ax[1, 1].set_ylabel("V [km/s]")
    ax[1, 1].grid(alpha=0.3); ax[1, 1].legend(); ax[1, 1].set_title("final velocity profile")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
