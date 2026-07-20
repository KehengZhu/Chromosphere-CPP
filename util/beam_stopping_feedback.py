#!/usr/bin/env python3
"""Show the thick-target beam's stopping-depth FEEDBACK (the loop becomes opaque).

The self-consistent beam is injected at the apex and attenuates down the loop by
collisional stopping (Emslie 1978): the volumetric heating at column depth N (from
the apex) is Q ∝ n_tot·(1+N/N_c)^{-δ/2}, N_c ∝ E_c². Because N = ∫n ds depends on
the CURRENT density, as evaporation fills the loop the column to the chromosphere
grows, so the beam stops higher up — a feedback a fixed window cannot show.

This recomputes Q(s,t) from each output frame's density (identical formula to
physics.hpp::beam_heating_rate) and tracks where the energy lands over time:
  - deposition-weighted mean HEIGHT (the stopping depth), and
  - the fraction of beam energy reaching the chromosphere (h < chromo top).

Usage:
  python util/beam_stopping_feedback.py [loop.dat] [thick.txt] [out.png] [E_cut_keV] [delta]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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
    phi = 0.5 * (a[:, 4] + a[:, 5])
    arc = np.cumsum(ds) - 0.5 * ds
    return ds, arc, phi, phi / g_si        # ds[m], arc[m], phi_g, height[m]


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


def column_from_apex(ntot, ds, a):
    """Trapezoidal column depth N(s) [m^-2] from the apex cell a, both legs."""
    N = np.zeros_like(ntot)
    for i in range(a - 1, -1, -1):
        N[i] = N[i + 1] + 0.5 * (ntot[i + 1] * ds[i + 1] + ntot[i] * ds[i])
    for i in range(a + 1, len(ntot)):
        N[i] = N[i - 1] + 0.5 * (ntot[i] * ds[i] + ntot[i - 1] * ds[i - 1])
    return N


def thick_target_Q(ntot, ds, phi, E_cut_keV, delta):
    """Deposition profile (unnormalized weight w ∝ Q); matches physics.hpp."""
    a = int(np.argmax(phi))
    N = column_from_apex(ntot, ds, a)
    N_c = 2.0e21 * E_cut_keV ** 2
    w = ntot * (1.0 + N / N_c) ** (-0.5 * delta)
    return w / np.sum(w * ds), a, N, N_c       # normalized so ∫Q ds = 1


def main():
    dat = sys.argv[1] if len(sys.argv) > 1 else "scenarios/data/loop_full_fine.dat"
    flare = sys.argv[2] if len(sys.argv) > 2 else "outputs/loops/loop_full_thick.txt"
    out = sys.argv[3] if len(sys.argv) > 3 else "visualization/loops/beam_stopping_feedback.png"
    E_cut = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0
    delta = float(sys.argv[5]) if len(sys.argv) > 5 else 5.0
    T_ON = 2.0
    T_OFF = float(sys.argv[6]) if len(sys.argv) > 6 else 62.0   # beam-off time
    os.makedirs(os.path.dirname(out), exist_ok=True)

    ds, arc, phi, h = parse_dat(dat)
    hMm, arcMm = h / 1e6, arc / 1e6
    a = int(np.argmax(phi))
    legA = np.arange(len(ds)) <= a
    h_chr = 1.2                                  # chromosphere top [Mm] (C7 ~1 Mm + a bit)
    iTR = int(np.argmin(np.abs(hMm[:a] - 1.5)))  # leg-A cell at h~1.5 Mm (just above TR)
    frames = read_frames(flare)

    ts, hcen, chrfrac, n_cor, colfp = [], [], [], [], []
    Qstack = []
    for t, xn in frames:
        ni = xn[:, 0] / m_i; nn = xn[:, 1] / m_i; ntot = ni + nn
        Q, _, N, N_c = thick_target_Q(ntot, ds, phi, E_cut, delta)
        ts.append(t)
        # deposition-weighted mean height + chromospheric fraction (both legs)
        Wds = Q * ds
        hcen.append(np.sum(Wds * hMm) / np.sum(Wds))
        chrfrac.append(np.sum(Wds[hMm < h_chr]) / np.sum(Wds))
        n_cor.append(ni[a])                       # apex density (loop filling)
        colfp.append(N[iTR] / N_c)                # coronal column apex->TR / N_c (opacity)
        Qstack.append(Q)
    ts = np.array(ts); hcen = np.array(hcen); chrfrac = np.array(chrfrac)
    n_cor = np.array(n_cor); colfp = np.array(colfp); Qstack = np.array(Qstack)

    on = (ts >= T_ON) & (ts <= T_OFF)
    print(f"beam window [{T_ON},{T_OFF}]s: deposition mean height "
          f"{hcen[on][0]:.2f} -> {hcen[on][-1]:.2f} Mm; "
          f"chromospheric fraction {chrfrac[on][0]*100:.0f}% -> {chrfrac[on][-1]*100:.0f}%; "
          f"apex n {n_cor[on][0]:.1e} -> {n_cor[on][-1]:.1e} m^-3; "
          f"footpoint column/N_c {colfp[on][0]:.2f} -> {colfp[on][-1]:.2f}")

    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(f"Thick-target beam stopping-depth feedback (E_c={E_cut:.0f} keV, "
                 f"δ={delta:.0f}): the loop becomes opaque as it fills", fontsize=12)

    # (A) deposition profile vs height on leg A, at start / mid / end of beam
    aA = ax[0, 0]
    for tt, c in [(T_ON + 1, "navy"), (0.5 * (T_ON + T_OFF), "darkorange"), (T_OFF - 1, "crimson")]:
        k = int(np.argmin(np.abs(ts - tt)))
        aA.semilogx(hMm[legA], Qstack[k][legA] / Qstack[k][legA].max(),
                    color=c, lw=1.8, label=f"t={ts[k]:.0f}s")
    aA.axvline(h_chr, color="gray", ls="--", lw=1, alpha=0.6)
    aA.text(h_chr * 1.05, 0.05, "chromo top", color="gray", fontsize=8, rotation=90, va="bottom")
    aA.set_xlabel("height [Mm]  (footpoint → apex)"); aA.set_ylabel("normalized deposition Q")
    aA.set_title("deposition profile migrates UP the loop"); aA.legend(fontsize=9)
    aA.grid(alpha=0.3, which="both")

    # (B) space-time deposition on leg A during the beam
    aB = ax[0, 1]
    onk = np.where(on)[0]
    Z = Qstack[onk][:, legA]
    Z = Z / Z.max(axis=1, keepdims=True)
    im = aB.pcolormesh(ts[onk], hMm[legA], Z.T, shading="auto", cmap="inferno",
                       norm=LogNorm(vmin=1e-3, vmax=1.0))
    aB.axhline(h_chr, color="cyan", ls="--", lw=1, alpha=0.7)
    aB.set_xlabel("t [s]"); aB.set_ylabel("height [Mm]"); aB.set_ylim(0, 8)
    aB.set_title("Q(height, t) on one leg (row-normalized)")
    fig.colorbar(im, ax=aB, label="Q / max")

    # (C) stopping depth + chromospheric fraction vs time (DURING the beam only —
    # the deposition profile is meaningful only while the beam is injecting).
    bm = (ts >= T_ON) & (ts <= T_OFF)
    aC = ax[1, 0]
    aC.plot(ts[bm], hcen[bm], color="C0", lw=2.0, label="deposition mean height")
    aC.axvspan(T_ON, T_OFF, color="orange", alpha=0.15, label="beam on")
    aC.set_xlabel("t [s]"); aC.set_ylabel("deposition mean height [Mm]", color="C0")
    aC.tick_params(axis="y", labelcolor="C0"); aC.set_xlim(0, T_OFF * 1.05)
    aC.grid(alpha=0.3)
    aC2 = aC.twinx()
    aC2.plot(ts[bm], chrfrac[bm] * 100, color="crimson", lw=2.0, ls="--")
    aC2.set_ylabel("% energy reaching chromosphere", color="crimson")
    aC2.tick_params(axis="y", labelcolor="crimson"); aC2.set_ylim(0, 105)
    aC.set_title("while beaming: stopping depth rises, chromospheric heating chokes off")

    # (D) the driver: coronal filling + opacity
    aD = ax[1, 1]
    aD.semilogy(ts, n_cor, color="purple", lw=1.8, label="apex $n_e$")
    aD.axvspan(T_ON, T_OFF, color="orange", alpha=0.15)
    aD.set_xlabel("t [s]"); aD.set_ylabel(r"apex $n_e$ [m$^{-3}$]", color="purple")
    aD.tick_params(axis="y", labelcolor="purple"); aD.set_xlim(0, min(ts[-1], T_OFF + 60))
    aD.grid(alpha=0.3, which="both")
    aD2 = aD.twinx()
    aD2.plot(ts, colfp, color="green", lw=1.8, ls="--")
    aD2.axhline(1.0, color="green", lw=0.8, alpha=0.4)
    aD2.set_ylabel(r"apex→TR coronal column / $N_c$", color="green")
    aD2.tick_params(axis="y", labelcolor="green")
    aD.set_xlim(0, T_OFF * 1.05)
    aD.set_title("loop fills (n rises) → coronal column → $N_c$ → opaque")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
