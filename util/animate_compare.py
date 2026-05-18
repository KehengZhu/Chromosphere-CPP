#!/usr/bin/env python3
"""Overlay two runs of chromo_main.cpp on the same 8-panel animation.

Usage:
    python3 animate_compare.py [full.txt] [explicit.txt] [out.mp4] [fps]

Defaults look for build/output_full.txt and build/output_explicit.txt, and
write util/model_c7_compare.mp4 at 15 fps.

The two runs may have different frame counts (different CFL-limited step
sizes). The animation is driven by the first run's frame times; for each
of those times the explicit run is sampled at its nearest frame.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from plot_output import read_frames, primitives


def panel_limits(series_list, log=False):
    arrs = []
    for s in series_list:
        arrs.append(np.concatenate(s))
    a = np.concatenate(arrs)
    if log:
        a = a[a > 0]
    lo, hi = float(np.min(a)), float(np.max(a))
    if log:
        return lo / 1.2, hi * 1.2
    pad = 0.05 * (hi - lo if hi > lo else max(abs(hi), 1.0))
    return lo - pad, hi + pad


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    full_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "build", "output_full.txt")
    expl_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "build", "output_explicit.txt")
    out_path  = sys.argv[3] if len(sys.argv) > 3 else os.path.join(here, "model_c7_compare.mp4")
    fps       = int(sys.argv[4]) if len(sys.argv) > 4 else 15

    xx_f, frames_f = read_frames(full_path)
    xx_e, frames_e = read_frames(expl_path)
    assert xx_f.size == xx_e.size and np.allclose(xx_f, xx_e), "grid mismatch between runs"
    xx = xx_f

    prims_f = [primitives(fr[2]) for fr in frames_f]
    prims_e = [primitives(fr[2]) for fr in frames_e]
    t_f = np.array([fr[0] for fr in frames_f])
    t_e = np.array([fr[0] for fr in frames_e])

    # Animate on the full run's timeline; pick nearest explicit frame per tick.
    e_idx_for_f = np.array([int(np.argmin(np.abs(t_e - tf))) for tf in t_f])

    # Per-field global limits taken over BOTH runs so the eye can compare.
    fields = ["nn", "u", "pn", "Tn", "ni", "v", "pi", "Ti"]
    log_yscale = {"nn": True, "pn": True, "ni": True, "pi": True}
    labels = {
        "nn": r"$n_n$ (m$^{-3}$)",
        "u":  r"$u$ (m/s)",
        "pn": r"$P_n$ (Pa)",
        "Tn": r"$T_n$ (K)",
        "ni": r"$n_i$ (m$^{-3}$)",
        "v":  r"$v$ (m/s)",
        "pi": r"$P_i$ (Pa)",
        "Ti": r"$T_i$ (K)",
    }
    # index in primitives() tuple: ni, nn, v, u, p_i, p_n, Ti, Tn
    idx_in_prim = {"ni": 0, "nn": 1, "v": 2, "u": 3, "pi": 4, "pn": 5, "Ti": 6, "Tn": 7}
    field_to_panel = {
        "nn": (0, 0), "u":  (0, 1), "pn": (0, 2), "Tn": (0, 3),
        "ni": (1, 0), "v":  (1, 1), "pi": (1, 2), "Ti": (1, 3),
    }

    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    title = fig.suptitle("")

    lines_full = {}
    lines_expl = {}
    for fld in fields:
        r, c = field_to_panel[fld]
        ax = axes[r, c]
        ipr = idx_in_prim[fld]
        series_full = [p[ipr] for p in prims_f]
        series_expl = [p[ipr] for p in prims_e]
        lo, hi = panel_limits([series_full, series_expl], log=log_yscale.get(fld, False))
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(lo, hi)
        if log_yscale.get(fld, False):
            ax.set_yscale("log")
        ax.set_title(labels[fld])
        ax.grid(True, alpha=0.3)
        (l_full,) = ax.plot(xx, series_full[0], lw=1.4, color="C0", label="semi-implicit (full)")
        (l_expl,) = ax.plot(xx, series_expl[0], lw=1.4, color="C3", linestyle="--", label=r"explicit ($R_I=0$)")
        lines_full[fld] = l_full
        lines_expl[fld] = l_expl

    for ax in axes[1, :]:
        ax.set_xlabel("height (km)")

    # Single legend at the top.
    handles = [lines_full["nn"], lines_expl["nn"]]
    fig.legend(handles, [h.get_label() for h in handles],
               loc="upper left", bbox_to_anchor=(0.01, 0.985), ncol=2, frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    def update(k):
        kf = k
        ke = int(e_idx_for_f[k])
        pf = prims_f[kf]
        pe = prims_e[ke]
        for fld in fields:
            ipr = idx_in_prim[fld]
            lines_full[fld].set_ydata(pf[ipr])
            lines_expl[fld].set_ydata(pe[ipr])
        title.set_text(
            f"Chromosphere — Model C7  |  semi-implicit t = {t_f[kf]:6.2f} s   vs   explicit t = {t_e[ke]:6.2f} s   ({k+1}/{len(t_f)})"
        )
        return list(lines_full.values()) + list(lines_expl.values()) + [title]

    anim = animation.FuncAnimation(fig, update, frames=len(t_f), interval=1000.0/fps, blit=False)
    writer = animation.FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    anim.save(out_path, writer=writer, dpi=120)
    print(f"wrote {out_path}  ({len(t_f)} frames @ {fps} fps = {len(t_f)/fps:.1f} s)")
    print(f"full     run: {len(t_f)} frames, t in [{t_f.min():.1f}, {t_f.max():.1f}] s")
    print(f"explicit run: {len(t_e)} frames, t in [{t_e.min():.1f}, {t_e.max():.1f}] s")


if __name__ == "__main__":
    main()
