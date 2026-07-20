#!/usr/bin/env python3
"""Overlay the time evolution of the γ=1.05 iso_t22k runs at three resolutions
(ns=1000/2000/4000) to show grid convergence of the mass flux ρV.

Mass flux ρV = (ρ_i+ρ_n)·V and velocity V are decoded straight from the
conserved RHO/MOM slots — they carry NO γ dependence — so no ISO_GAMMA is needed
here (unlike animate_isentropic.py, which decodes pressure/T from energy).

Panels (three resolutions overlaid in each):
    (a) mass flux ρV over the FULL domain (h>1800 km region shaded)
    (b) mass flux ρV, ZOOMED to h>1800 km only  (the TR evaporation region)
    (c) velocity V, zoomed to h>1800 km only      (the evaporation upflow)

The three runs have different frame cadence; frames are synced by picking, for
each target time on a common axis [0, min(t_max)], the nearest frame from each run.

Usage:
    python util/animate_iso_gamma_compare.py [out.mp4]
Defaults to the three outputs/model_column/iso_t22k*_gamma105.txt files ->
visualization/model_column/iso_t22k_gamma105_massflux_compare.mp4. Set ANIM_MAX_FRAMES to cap
the rendered-frame count (default 300; sampled uniformly in time).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _anim_parallel import save_frames_parallel

H_ZOOM = 1800.0  # km — TR/upper-chromosphere zoom threshold

RUNS = [
    ("outputs/model_column/iso_t22k_ns1000_gamma105.txt", "ns=1000", "tab:red",    1.6),
    ("outputs/model_column/iso_t22k_ns2000_gamma105.txt", "ns=2000", "tab:orange", 1.3),
    ("outputs/model_column/iso_t22k_gamma105.txt",        "ns=4000", "tab:blue",   1.0),
]


def load(fn):
    """Return (h[km], times[np], rhoV[list of np], V[list of np]). Stops at the
    first non-finite frame so a run that later NaNs still animates its valid window."""
    f = open(fn)
    ns, neq = map(int, f.readline().split())
    h = np.array(list(map(float, f.readline().split())))
    times, rhoV, V = [], [], []
    cur = None
    t = None

    def flush():
        if cur is None or len(cur) != ns:
            return True
        U = np.array(cur)
        rho_i, rho_n = U[:, 0], U[:, 1]
        v = U[:, 2] / rho_i
        mf = (rho_i + rho_n) * v
        if not (np.all(np.isfinite(mf)) and np.all(np.isfinite(v))):
            return False
        times.append(t); rhoV.append(mf); V.append(v)
        return True

    for line in f:
        if line.startswith('# t ='):
            if cur is not None and not flush():
                break
            t = float(line.split('t =')[1].split('step')[0])
            cur = []
        else:
            try:
                cur.append(list(map(float, line.split())))
            except ValueError:
                cur.append([])
    else:
        flush()
    return h, np.array(times), rhoV, V


def _render(k, tmpdir, tt, per_run, ylims):
    """per_run: list of (label,color,lw,h,rhoV,V) at this frame. ylims: dict of axis limits."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    axA, axB, axC = axes
    for label, color, lw, h, mf, v in per_run:
        axA.plot(h, mf, color=color, lw=lw, label=label)
        m = h > H_ZOOM
        axB.plot(h[m], mf[m], color=color, lw=lw, label=label)
        axC.plot(h[m], v[m] / 1e3, color=color, lw=lw, label=label)

    axA.axvspan(H_ZOOM, ylims['hmax'], color='0.85', alpha=0.5, zorder=0)
    axA.axhline(0, color='0.6', lw=0.6)
    axA.set_title("mass flux ρV — full domain")
    axA.set_xlabel("height [km]"); axA.set_ylabel("ρV [kg m⁻² s⁻¹]")
    axA.set_xlim(0, ylims['hmax']); axA.set_ylim(*ylims['A'])
    axA.legend(loc='upper right', fontsize=9)

    axB.axhline(0, color='0.6', lw=0.6)
    axB.set_title(f"mass flux ρV — h > {H_ZOOM:.0f} km (TR)")
    axB.set_xlabel("height [km]"); axB.set_ylabel("ρV [kg m⁻² s⁻¹]")
    axB.set_xlim(H_ZOOM, ylims['hmax']); axB.set_ylim(*ylims['B'])
    axB.legend(loc='upper left', fontsize=9)

    axC.axhline(0, color='0.6', lw=0.6)
    axC.set_title(f"velocity V — h > {H_ZOOM:.0f} km (TR)")
    axC.set_xlabel("height [km]"); axC.set_ylabel("V [km s⁻¹]")
    axC.set_xlim(H_ZOOM, ylims['hmax']); axC.set_ylim(*ylims['C'])
    axC.legend(loc='upper left', fontsize=9)

    fig.suptitle(f"iso_t22k  γ=1.05  grid convergence   t = {tt:7.1f} s", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=110)
    plt.close(fig)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else \
        "visualization/model_column/iso_t22k_gamma105_massflux_compare.mp4"
    runs = [(load(fn), lab, col, lw) for fn, lab, col, lw in RUNS]

    t_common = min(r[0][1][-1] for r in runs)   # min of per-run max time
    n = int(os.environ.get("ANIM_MAX_FRAMES", "300"))
    target_t = np.linspace(0.0, t_common, n)

    # Global y-limits (robust) over the rendered window so the axes don't jump.
    A_vals, B_vals, C_vals, hmax = [], [], [], 0.0
    for (h, times, rhoV, V), *_ in runs:
        hmax = max(hmax, h.max())
        m = h > H_ZOOM
        sel = times <= t_common
        for mf in np.array(rhoV, dtype=object)[sel]:
            A_vals.append(mf)
            B_vals.append(mf[m])
        for v in np.array(V, dtype=object)[sel]:
            C_vals.append(v[m] / 1e3)
    A_all = np.concatenate(A_vals); B_all = np.concatenate(B_vals); C_all = np.concatenate(C_vals)

    def pad(lo, hi, f=0.08):
        d = (hi - lo) * f or abs(hi) * f or 1.0
        return lo - d, hi + d
    ylims = {
        'hmax': hmax,
        'A': pad(np.percentile(A_all, 0.2), np.percentile(A_all, 99.8)),
        'B': pad(np.percentile(B_all, 0.5), np.percentile(B_all, 99.5)),
        'C': pad(min(0.0, C_all.min()), C_all.max()),
    }

    frame_args = []
    for k, tt in enumerate(target_t):
        per_run = []
        for (h, times, rhoV, V), lab, col, lw in runs:
            j = int(np.argmin(np.abs(times - tt)))
            per_run.append((lab, col, lw, h, rhoV[j], V[j]))
        frame_args.append((k, tt, per_run, ylims))

    fps = int(os.environ.get("ANIM_FPS", "20"))
    print("runs: " + ", ".join(f"{lab}({len(r[1])}f, t<={r[1][-1]:.0f}s)"
                               for r, lab, _, _ in runs))
    print(f"common t_max = {t_common:.1f} s, rendering {n} frames")
    save_frames_parallel(_render, frame_args, out, fps)


if __name__ == "__main__":
    main()
