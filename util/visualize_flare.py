#!/usr/bin/env python3
"""Visualize explosive chromospheric evaporation (model_flare scenario).

Reads chromo_main output for the explosive and gentle-control beam runs and
produces, into util/visualization/:

  1. flare_height_time.png  — height–time diagrams of T_i and v for the
     explosive run (the classic evaporation plot): the conductive/heating front,
     the upward evaporation jet, and the downward chromospheric condensation.
  2. flare_compare_snapshot.png — explosive vs gentle profiles at the peak of
     the impulsive phase (T, v, n_i, beam-region marker).
  3. flare_timeseries.png — peak T, max upflow, max (down)flow vs time for both.

Usage:
  python3 util/visualize_flare.py [explosive.txt] [gentle.txt] [t_max_s]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "visualization")
os.makedirs(OUT, exist_ok=True)

m_i = 1.6726219e-27
k_b = 1.380649e-23
gamma = 5.0 / 3.0


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


def stack(frames, tmax):
    ts = [t for t, _ in frames if t <= tmax]
    T = np.array([prim(d)[3] for t, d in frames if t <= tmax])
    V = np.array([prim(d)[2] for t, d in frames if t <= tmax])
    N = np.array([prim(d)[0] for t, d in frames if t <= tmax])
    return np.array(ts), T, V, N


def main():
    exp_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/flare_explosive.txt"
    gen_path = sys.argv[2] if len(sys.argv) > 2 else "outputs/flare_gentle.txt"
    tmax = float(sys.argv[3]) if len(sys.argv) > 3 else 25.0

    H, fe = read_frames(exp_path)
    _, fg = read_frames(gen_path)
    te, Te, Ve, Ne = stack(fe, tmax)
    tg, Tg, Vg, Ng = stack(fg, tmax)

    # ---- 1. height-time diagrams (explosive) --------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    TT, HH = np.meshgrid(te, H, indexing="ij")
    pcm0 = ax[0].pcolormesh(TT, HH, np.log10(np.clip(Te, 3e3, None)),
                            shading="auto", cmap="inferno", vmin=3.6, vmax=7.4)
    ax[0].set_title("Explosive evaporation — $\\log_{10} T_i$ [K]")
    cb0 = fig.colorbar(pcm0, ax=ax[0]); cb0.set_label("$\\log_{10} T$ [K]")
    vmax = np.nanpercentile(np.abs(Ve), 99) / 1e3
    pcm1 = ax[1].pcolormesh(TT, HH, Ve / 1e3, shading="auto", cmap="RdBu_r",
                            vmin=-vmax, vmax=vmax)
    ax[1].set_title("Field-aligned velocity $v$ [km/s]\n(red = upflow / evaporation, blue = condensation)")
    cb1 = fig.colorbar(pcm1, ax=ax[1]); cb1.set_label("$v$ [km/s]")
    for a in ax:
        a.set_xlabel("time [s]")
        a.axhspan(1600, 2000, color="cyan", alpha=0.12)   # beam deposition window
    ax[0].set_ylabel("height [km]")
    fig.suptitle("Explosive chromospheric evaporation (Fisher et al. 1985) — beam $F=5\\times10^{10}$ erg cm$^{-2}$ s$^{-1}$, on 2–12 s",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(OUT, "flare_height_time.png"), dpi=140)
    plt.close(fig)
    print("wrote flare_height_time.png")

    # ---- 2. snapshot comparison at impulsive peak ---------------------------
    kt = int(np.argmax([np.nanmax(prim(d)[3]) for t, d in fe if t <= tmax]))
    t_pk = te[kt]
    kg = int(np.argmin(np.abs(tg - t_pk)))
    nie, nne, ve, Tie = prim(fe[[i for i, (t, _) in enumerate(fe) if t <= tmax][kt]][1])
    nig, nng, vg, Tig = prim(fg[[i for i, (t, _) in enumerate(fg) if t <= tmax][kg]][1])

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    ax[0].plot(H, Tie, "r", label="explosive")
    ax[0].plot(H, Tig, "b", label="gentle (control)")
    ax[0].set_yscale("log"); ax[0].set_title(f"$T_i$ [K] at t={t_pk:.1f} s"); ax[0].legend()
    ax[1].plot(H, ve / 1e3, "r", label="explosive")
    ax[1].plot(H, vg / 1e3, "b", label="gentle")
    ax[1].axhline(0, color="k", lw=0.6); ax[1].set_title("$v$ [km/s] (+ up / − down)"); ax[1].legend()
    ax[2].plot(H, nie, "r", label="explosive")
    ax[2].plot(H, nig, "b", label="gentle")
    ax[2].set_yscale("log"); ax[2].set_title("$n_i$ [m$^{-3}$]"); ax[2].legend()
    for a in ax:
        a.set_xlabel("height [km]"); a.axvspan(1600, 2000, color="cyan", alpha=0.12)
        a.grid(True, alpha=0.3)
    fig.suptitle("Explosive vs gentle chromospheric response at the impulsive peak", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(OUT, "flare_compare_snapshot.png"), dpi=140)
    plt.close(fig)
    print("wrote flare_compare_snapshot.png")

    # ---- 3. time series -----------------------------------------------------
    def series(ts, T, V):
        return (np.nanmax(T, axis=1),
                np.nanmax(V, axis=1) / 1e3,
                np.nanmin(V, axis=1) / 1e3)
    pTe, upe, dne = series(te, Te, Ve)
    pTg, upg, dng = series(tg, Tg, Vg)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].plot(te, pTe, "r", label="explosive"); ax[0].plot(tg, pTg, "b", label="gentle")
    ax[0].set_yscale("log"); ax[0].set_title("peak $T_i$ [K]"); ax[0].legend()
    ax[1].plot(te, upe, "r", label="explosive up"); ax[1].plot(te, dne, "r--", label="explosive down")
    ax[1].plot(tg, upg, "b", label="gentle up"); ax[1].plot(tg, dng, "b--", label="gentle down")
    ax[1].axhline(0, color="k", lw=0.6); ax[1].set_title("max upflow / downflow [km/s]"); ax[1].legend()
    for a in ax:
        a.set_xlabel("time [s]"); a.axvspan(2, 12, color="orange", alpha=0.12)  # beam on
        a.grid(True, alpha=0.3)
    fig.suptitle("Evaporation time series (shaded = beam on, 2–12 s)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(OUT, "flare_timeseries.png"), dpi=140)
    plt.close(fig)
    print("wrote flare_timeseries.png")


if __name__ == "__main__":
    main()
