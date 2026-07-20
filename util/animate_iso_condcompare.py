#!/usr/bin/env python3
"""Overlay the mass flux ρV of two BEST-well-balanced iso_t22k runs that differ
ONLY in whether heat conduction is on — a demonstration that, once the spurious
base-drainage truncation is removed by the well-balancing stack
(ISO_INNER_WB + ISO_INNER_WB_RHO + ISO_EQ_WB, plus MC3 / log-recon), the remaining
mass flux is entirely physical (conduction-driven evaporation):

    conduction ON  -> a TR evaporation upflow develops; base ρV stays ~0
    conduction OFF -> ρV stays flat at ~round-off everywhere (equilibrium held)

Mass flux ρV = (ρ_i+ρ_n)·V and V are decoded straight from the conserved RHO/MOM
slots (NO γ dependence), so no ISO_GAMMA is needed here.

Panels (both runs overlaid):
    (a) ρV over the FULL domain (h>1800 km region shaded)
    (b) ρV, ZOOMED to h>1800 km (the TR evaporation region)
    (c) V,  ZOOMED to h>1800 km

Frames are synced by nearest-time matching on a common axis [0, min(t_max)].

Usage:
    python util/animate_iso_condcompare.py [condON.txt] [condOFF.txt] [out.mp4]
Defaults to the best-WB ns=2000 γ=1.05 pair ->
visualization/model_column/iso_t22k_ns2000_gamma105_bestwb_condcompare.mp4.
ANIM_MAX_FRAMES caps the rendered-frame count (default 300); ANIM_FPS sets fps.
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

CONDON  = sys.argv[1] if len(sys.argv) > 1 else "outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt"
CONDOFF = sys.argv[2] if len(sys.argv) > 2 else "outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condoff.txt"
OUT     = sys.argv[3] if len(sys.argv) > 3 else "visualization/model_column/iso_t22k_ns2000_gamma105_bestwb_condcompare.mp4"
LABEL   = os.environ.get("ANIM_LABEL", "best well-balanced (INNER_WB+RHO+EQ_WB, MC3)")

RUNS = [
    (CONDON,  "conduction ON",  "tab:red",  1.5),
    (CONDOFF, "conduction OFF", "tab:blue", 1.5),
]


def load(fn):
    """Return (h[km], times[np], rhoV[list], V[list]); stops at first non-finite frame."""
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

    fig.suptitle(f"iso_t22k  γ=1.05  {LABEL}  "
                 f"— mass flux, conduction ON vs OFF   t = {tt:7.1f} s", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=110)
    plt.close(fig)


def main():
    runs = [(load(fn), lab, col, lw) for fn, lab, col, lw in RUNS]
    t_common = min(r[0][1][-1] for r in runs)
    n = int(os.environ.get("ANIM_MAX_FRAMES", "300"))
    target_t = np.linspace(0.0, t_common, n)

    A_vals, B_vals, C_vals, hmax = [], [], [], 0.0
    for (h, times, rhoV, V), *_ in runs:
        hmax = max(hmax, h.max())
        m = h > H_ZOOM
        sel = times <= t_common
        for mf in np.array(rhoV, dtype=object)[sel]:
            A_vals.append(mf); B_vals.append(mf[m])
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
    print("runs: " + ", ".join(f"{lab}({len(r[1])}f, t<={r[1][-1]:.0f}s)" for r, lab, _, _ in runs))
    print(f"common t_max = {t_common:.1f} s, rendering {n} frames")
    save_frames_parallel(_render, frame_args, OUT, fps)


if __name__ == "__main__":
    main()
