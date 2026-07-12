#!/usr/bin/env python3
"""Evolution frames for the explosive chromospheric evaporation run.

Reads the refined (500-cell) model_flare output and produces, into
visualization/:

  1. flare_evolution_frames.png — a frame strip: rows = (T_i, v, densities),
     columns = selected times anchored on the physical events of the run
     (initial state, explosive upflow launch, condensation peak, temperature
     peak, late relaxation). This is the "pick some frames to show evolution"
     figure for the writeup.
  2. flare_height_time_500.png — continuous height-time diagrams of T_i and v
     for the same run (the classic evaporation view), regenerated from the
     refined data so it is consistent with the frame strip.

It also prints the diagnostics (peak T, max up/downflow, deposition-layer
temperature ~1 s into heating, gentle-control peaks) used to keep the writeup
numbers in sync with the data.

Usage:
  python3 util/visualize_flare_frames.py [explosive_500.txt] [gentle.txt]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "visualization")
os.makedirs(OUT, exist_ok=True)

m_i = 1.6726219e-27
k_b = 1.380649e-23

# Beam deposition window and on-window (model_flare defaults).
H_LO, H_HI = 1600.0, 2000.0
T_ON, T_OFF = 2.0, 12.0


def read_frames(path):
    raw = [l for l in open(path).read().splitlines() if l.strip()]
    ns, neq = map(int, raw[0].split())
    H = np.fromstring(raw[1], sep=" ")
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
    return H, frames


def prim(xn):
    ni = xn[:, 0] / m_i
    nn = xn[:, 1] / m_i
    v = xn[:, 2] / xn[:, 0]
    pi = 2.0 / 3.0 * xn[:, 4] - 1.0 / 3.0 * m_i * ni * v * v
    Ti = pi / (2.0 * ni * k_b)
    return ni, nn, v, Ti


def main():
    exp_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/flare_explosive_500.txt"
    gen_path = sys.argv[2] if len(sys.argv) > 2 else "outputs/flare_gentle.txt"

    H, fe = read_frames(exp_path)
    ts = np.array([t for t, _ in fe])

    # ---- diagnostics over the whole run ------------------------------------
    pkT = np.array([np.nanmax(prim(d)[3]) for _, d in fe])
    up = np.array([np.nanmax(prim(d)[2]) for _, d in fe])
    dn = np.array([np.nanmin(prim(d)[2]) for _, d in fe])
    i_pkT, i_up, i_dn = np.nanargmax(pkT), np.nanargmax(up), np.nanargmin(dn)

    # deposition-layer temperature ~1 s into the heating (t ~ T_ON + 1)
    beam = (H >= H_LO) & (H <= H_HI)
    i_1s = int(np.argmin(np.abs(ts - (T_ON + 1.0))))
    T_dep_1s = np.nanmax(prim(fe[i_1s][1])[3][beam])

    print("=== refined explosive run diagnostics (%s) ===" % exp_path)
    print(f"  grid: {len(H)} cells, h=[{H[0]:.0f},{H[-1]:.0f}] km, t<={ts[-1]:.1f} s")
    print(f"  peak T_i        = {pkT[i_pkT]:.3e} K   at t={ts[i_pkT]:.2f} s")
    print(f"  deposition T_i  = {T_dep_1s:.3e} K   at t={ts[i_1s]:.2f} s (~1 s into beam)")
    print(f"  max upflow      = {up[i_up] / 1e3:.1f} km/s at t={ts[i_up]:.2f} s")
    print(f"  max downflow    = {dn[i_dn] / 1e3:.1f} km/s at t={ts[i_dn]:.2f} s")

    if os.path.exists(gen_path):
        _, fg = read_frames(gen_path)
        gpkT = np.nanmax([np.nanmax(prim(d)[3]) for _, d in fg])
        gup = np.nanmax([np.nanmax(prim(d)[2]) for _, d in fg]) / 1e3
        print(f"  gentle control  : peak T_i={gpkT:.3e} K, max upflow={gup:.1f} km/s ({gen_path})")

    # ---- 1. frame strip -----------------------------------------------------
    # Anchor frames on the physical events, ordered in time, de-duplicated.
    targets = [0.0, ts[i_up], ts[i_dn], ts[i_pkT], ts[-1]]
    idx = sorted({int(np.argmin(np.abs(ts - tt))) for tt in targets})
    cols = len(idx)

    fig, ax = plt.subplots(3, cols, figsize=(3.1 * cols, 8.4),
                           sharex=True, sharey="row")
    for c, k in enumerate(idx):
        t, d = fe[k]
        ni, nn, v, Ti = prim(d)
        Tplot = np.clip(Ti, 3e3, 5e7)
        a0, a1, a2 = ax[0, c], ax[1, c], ax[2, c]
        a0.semilogy(H, Tplot, "r", lw=1.4)
        a1.plot(H, v / 1e3, "b", lw=1.4)
        a1.axhline(0, color="k", lw=0.6)
        a2.semilogy(H, np.clip(ni, 1, None), "C3", lw=1.3, label="$n_i$")
        a2.semilogy(H, np.clip(nn, 1, None), "C0", lw=1.3, label="$n_n$")
        for a in (a0, a1, a2):
            a.axvspan(H_LO, H_HI, color="orange", alpha=0.12)
            a.grid(True, alpha=0.25)
        # mark which phase this frame is
        on = "beam on" if T_ON <= t <= T_OFF else ("pre-beam" if t < T_ON else "post-beam")
        a0.set_title(f"$t={t:.1f}$ s\n({on})", fontsize=10)
        a2.set_xlabel("height [km]")
    ax[0, 0].set_ylabel("$T_i$ [K]")
    ax[1, 0].set_ylabel("$v$ [km s$^{-1}$]\n(+ up / $-$ down)")
    ax[2, 0].set_ylabel("number density [m$^{-3}$]")
    ax[0, 0].set_ylim(3e3, 6e7)
    ax[2, 0].legend(loc="lower left", fontsize=8)
    fig.suptitle(
        "Explosive chromospheric evaporation (Fisher et al. 1985), refined 500-cell run\n"
        "beam $F_e=5\\times10^{10}$ erg cm$^{-2}$ s$^{-1}$ (above threshold), on 2--12 s, "
        "deposited over 1600--2000 km (shaded)",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(OUT, "flare_evolution_frames.png"), dpi=140)
    plt.close(fig)
    print("wrote flare_evolution_frames.png  (frames at t =",
          ", ".join(f"{ts[k]:.1f}" for k in idx), "s)")

    # ---- 2. height-time diagrams (refined) ---------------------------------
    T = np.array([prim(d)[3] for _, d in fe])
    V = np.array([prim(d)[2] for _, d in fe])
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    TT, HH = np.meshgrid(ts, H, indexing="ij")
    pcm0 = ax[0].pcolormesh(TT, HH, np.log10(np.clip(T, 3e3, 5e7)),
                            shading="auto", cmap="inferno", vmin=3.6, vmax=7.3)
    ax[0].set_title("$\\log_{10} T_i$ [K]")
    fig.colorbar(pcm0, ax=ax[0]).set_label("$\\log_{10} T$ [K]")
    vm = np.nanpercentile(np.abs(V), 99) / 1e3
    pcm1 = ax[1].pcolormesh(TT, HH, V / 1e3, shading="auto", cmap="RdBu_r",
                            vmin=-vm, vmax=vm)
    ax[1].set_title("$v$ [km s$^{-1}$]  (red = evaporation, blue = condensation)")
    fig.colorbar(pcm1, ax=ax[1]).set_label("$v$ [km s$^{-1}$]")
    for a in ax:
        a.set_xlabel("time [s]")
        a.axhspan(H_LO, H_HI, color="cyan", alpha=0.12)
    ax[0].set_ylabel("height [km]")
    fig.suptitle("Explosive evaporation height--time (refined 500-cell run)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUT, "flare_height_time_500.png"), dpi=140)
    plt.close(fig)
    print("wrote flare_height_time_500.png")


if __name__ == "__main__":
    main()
