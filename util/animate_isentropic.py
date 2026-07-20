#!/usr/bin/env python3
"""Render an MP4 of the time evolution of a model_isentropic run.

The isentropic domain is a thin sub-Mm slab (chromosphere + TR base, no resolved
corona), so a single linear height axis is used — no split x-axis needed.

Usage:
    python util/animate_isentropic.py [input.txt] [output.mp4] [fps]

Defaults to outputs/_archive/iso_stage2b.txt (conduction + radiative cooling: the localized
TR + condensation downflow front) -> visualization/_archive/iso_stage2b_evolution.mp4.
Set ANIM_MAX_FRAMES to cap the number of rendered frames (default 400; the run is
subsampled uniformly in time).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _anim_parallel import save_frames_parallel

MP = 1.6726219e-27
KB = 1.380649e-23
G = 273.95
F_ION = 1e-3  # frozen ionization fraction used in the iso staged runs
# Adiabatic index used to decode pressure/T from the conserved energy E = p/(γ−1)
# + KE + PE. The .txt stores no γ, so it MUST match the run's ISO_GAMMA (default
# 5/3). Set ISO_GAMMA=1.05 in the environment to plot a near-isothermal run.
GAMMA = float(os.environ.get("ISO_GAMMA", 5.0 / 3.0))
GM1 = GAMMA - 1.0


def load(fn):
    """Read frames, stopping at the first non-finite frame (so a run that later
    NaNs — e.g. the resolved-corona drainage — still animates its valid window)."""
    f = open(fn)
    ns, neq = map(int, f.readline().split())
    h = np.array(list(map(float, f.readline().split())))
    frames = []
    cur = None
    t = None
    stop = False

    def flush():
        nonlocal stop
        if cur is None or stop:
            return
        rows = [r for r in cur if len(r) == neq]
        if len(rows) != ns:
            stop = True
            return
        arr = np.array(rows)
        if not np.isfinite(arr).all():
            stop = True       # first NaN/inf frame → end the valid window
            return
        frames.append((t, arr))

    for line in f:
        if stop:
            break
        if line.startswith('#'):
            flush()
            t = float(line.split('t =')[1].split('step')[0])
            cur = []
        else:
            try:
                cur.append(list(map(float, line.split())))
            except ValueError:
                cur.append([])
    flush()
    return h, frames


def primitives(U, phi):
    rho_i, rho_n = U[:, 0], U[:, 1]
    V = U[:, 2] / rho_i
    Un = U[:, 3] / rho_n
    p_i = GM1 * U[:, 4] - 0.5 * GM1 * rho_i * V * V - GM1 * rho_i * phi
    p_n = GM1 * U[:, 5] - 0.5 * GM1 * rho_n * Un * Un - GM1 * rho_n * phi
    rho = rho_i + rho_n
    # Single-fluid T with the LOCAL ionization (p = (2 n_i + n_n) k T); n_i=rho_i/MP,
    # n_n=rho_n/MP. Valid for both the uniform-f isentrope runs and the realistic
    # frozen-ionization C7/corona runs (neutral chromosphere → ionized corona).
    T = (p_i + p_n) * MP / ((2.0 * rho_i + rho_n) * KB)
    n_e = rho_i / MP   # n_e = n_i (single ionization)
    n_n = rho_n / MP
    return V, T, rho, (p_i + p_n), n_e, n_n


def kappa_e(n_e, n_n, T):
    """Electron field-aligned heat conductivity [SI], matching physics.hpp eq 53
    (Spitzer e–i numerator + Vranjes & Krstic e–n denominator term)."""
    return 9.2048e-12 * n_e * T ** 2.5 / (n_e + 2.836e-11 * n_n * T ** 2)


def kappa_n(n_i, n_n, T):
    """Neutral field-aligned heat conductivity [SI], matching physics.hpp eq 59
    (single fluid here, so T_i = T_n = T)."""
    return (0.0342006 * n_n * T) / (
        1.20613 * n_i * np.sqrt(2.0 * T) + 1.70573 * n_n * np.sqrt(T))


def heat_flux(h_km, T, n_e, n_n):
    """Field-aligned conductive heat flux q = -(κ_e + κ_n) dT/ds [W/m²].
    Positive = up the field (toward the corona); the strong negative spike at the
    TR is the downward conductive flux that drives evaporation. Spitzer flux (no
    TRAC broadening — that factor only redistributes the TR while preserving the
    integral)."""
    s = h_km * 1.0e3                       # km → m
    dTds = np.gradient(T, s)
    return -(kappa_e(n_e, n_n, T) + kappa_n(n_e, n_n, T)) * dTds


def _render_frame(k, tmpdir, h, V, T, rho, p, q, mflux, t, nF, ylims, title, xsplit):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 3, figsize=(17, 8))
    fig.suptitle(f"{title}  |  t = {t:7.1f} s   ({k + 1}/{nF})", fontsize=13)

    ax[0, 0].plot(h, T, lw=1.4, color="C3")
    ax[0, 0].set_ylabel("T [K]")
    ax[0, 0].set_title("Temperature")
    ax[0, 0].set_ylim(*ylims["T"])

    ax[0, 1].plot(h, V / 1e3, lw=1.4, color="C0")
    ax[0, 1].axhline(0, c="k", lw=0.5)
    ax[0, 1].set_ylabel("V [km/s]")
    ax[0, 1].set_title("Velocity  (>0 up = evaporation, <0 down = condensation)")
    ax[0, 1].set_ylim(*ylims["V"])

    ax[0, 2].plot(h, q / 1e3, lw=1.4, color="C1")
    ax[0, 2].axhline(0, c="k", lw=0.5)
    ax[0, 2].set_ylabel(r"$q$ [kW/m²]")
    ax[0, 2].set_title("Conductive heat flux  (>0 up; <0 = downward into TR)")
    ax[0, 2].set_ylim(*ylims["q"])

    ax[1, 0].semilogy(h, rho, lw=1.4, color="C2")
    ax[1, 0].set_ylabel(r"$\rho$ [kg/m³]")
    ax[1, 0].set_title("Mass density")
    ax[1, 0].set_ylim(*ylims["rho"])

    ax[1, 1].semilogy(h, p, lw=1.4, color="C4")
    ax[1, 1].set_ylabel("p [Pa]")
    ax[1, 1].set_title("Pressure")
    ax[1, 1].set_ylim(*ylims["p"])

    ax[1, 2].plot(h, mflux, lw=1.4, color="C5")
    ax[1, 2].axhline(0, c="k", lw=0.5)
    ax[1, 2].set_ylabel(r"$\rho V$ [kg/m²/s]")
    ax[1, 2].set_title("Mass flux  (>0 up = evaporation, <0 down = condensation)")
    # ax[1, 2].set_ylim(*ylims["mflux"])

    panels = [ax[0, 0], ax[0, 1], ax[0, 2], ax[1, 0], ax[1, 1], ax[1, 2]]

    # Split x-axis for the resolved-corona domain: zoom the thin sub-Mm chromosphere
    # +TR, compress the extended corona (continuous piecewise-linear data→plot map),
    # else a single linear axis. Dashed line marks the chromosphere/TR→corona split.
    if xsplit is not None:
        kd, kp, z_split, ticks = xsplit
        fwd = lambda x: np.interp(x, kd, kp)
        inv = lambda qq: np.interp(qq, kp, kd)
        for a in panels:
            a.set_xscale("function", functions=(fwd, inv))
            a.set_xlim(h.min(), h.max())
            a.set_xticks(ticks)
            a.axvline(z_split, color="gray", ls="--", lw=0.9, alpha=0.6)
            a.set_xlabel("height [km]  (split: chromosphere+TR | corona)")
            a.grid(True, alpha=0.3)
    else:
        for a in panels:
            a.set_xlim(h.min(), h.max())
            a.set_xlabel("height [km]")
            a.grid(True, alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def main():
    in_path = sys.argv[1] if len(sys.argv) > 1 else "outputs/_archive/iso_stage2b.txt"
    if len(sys.argv) > 2:
        out_path = sys.argv[2]
    else:
        stem = os.path.splitext(os.path.basename(in_path))[0]
        out_path = os.path.join("visualization", f"{stem}_evolution.mp4")
    fps = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    max_frames = int(os.environ.get("ANIM_MAX_FRAMES", "400"))
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    h, frames = load(in_path)
    phi = G * (h - h[0]) * 1000.0

    # Subsample uniformly to at most max_frames.
    nAll = len(frames)
    if nAll > max_frames:
        idx = np.linspace(0, nAll - 1, max_frames).round().astype(int)
        idx = np.unique(idx)
        frames = [frames[i] for i in idx]
    nF = len(frames)
    print(f"read {nAll} frames, using {nF}, ns = {h.size}, "
          f"t in [{frames[0][0]:.1f}, {frames[-1][0]:.1f}] s")

    prims = [primitives(fr[1], phi) for fr in frames]  # (V, T, rho, p, n_e, n_n)
    # Field-aligned conductive heat flux per frame (kW/m²), >0 up the field.
    qflux = [heat_flux(h, pr[1], pr[4], pr[5]) for pr in prims]
    # Field-aligned mass flux ρV [kg/m²/s] per frame, >0 up the field (evaporation).
    mflux = [pr[2] * pr[0] for pr in prims]

    def lim(arrs, log=False, pad=0.05):
        a = np.concatenate(arrs)
        a = a[np.isfinite(a)]
        if log:
            a = a[a > 0]
            lo, hi = float(a.min()), float(a.max())
            return lo / 1.3, hi * 1.3
        lo, hi = float(a.min()), float(a.max())
        d = pad * (hi - lo if hi > lo else max(abs(hi), 1.0))
        return lo - d, hi + d

    ylims = {
        "T": lim([p[1] for p in prims]),
        "V": lim([p[0] / 1e3 for p in prims]),
        "rho": lim([p[2] for p in prims], log=True),
        "p": lim([p[3] for p in prims], log=True),
        "q": lim([qi / 1e3 for qi in qflux]),
        "mflux": lim(mflux),
    }

    stem = os.path.splitext(os.path.basename(in_path))[0]
    title = f"model_isentropic — {stem}"

    # Split x-axis for a tall (resolved-corona) domain (> 3 Mm): zoom the thin
    # chromosphere+TR into CHROMO_FRAC of the width, compress the corona into the
    # rest. The split is at the TR — the first height where the initial temperature
    # exceeds 3e4 K (fallback: 20% up the domain). Thin slabs keep a linear axis.
    CHROMO_FRAC = 0.55
    xsplit = None
    if h[-1] - h[0] > 3000.0:
        T0 = prims[0][1]
        hot = np.where(T0 > 3.0e4)[0]
        z_split = float(h[hot[0]]) if hot.size else float(h[0] + 0.2 * (h[-1] - h[0]))
        kd = np.array([h[0], z_split, h[-1]])
        kp = np.array([0.0, CHROMO_FRAC, 1.0])
        ch = [t for t in (0, 500, 1000, 1500, 2000) if h[0] <= t <= z_split]
        co = list(np.linspace(z_split, h[-1], 5)[1:])
        ticks = np.array(ch + co)
        xsplit = (kd, kp, z_split, ticks)
        print(f"split x-axis: chromosphere+TR [{h[0]:.0f},{z_split:.0f}] km | "
              f"corona [{z_split:.0f},{h[-1]:.0f}] km")

    save_frames_parallel(
        _render_frame,
        [(k, h, prims[k][0], prims[k][1], prims[k][2], prims[k][3], qflux[k],
          mflux[k], frames[k][0], nF, ylims, title, xsplit) for k in range(nF)],
        out_path, fps,
    )


if __name__ == "__main__":
    main()
