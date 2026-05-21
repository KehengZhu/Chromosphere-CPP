#!/usr/bin/env python3
"""Overlay two runs of chromo_main.cpp on the same 8-panel animation,
comparing ionization-on (Stage E enabled) against ionization-off.

Usage:
    python3 animate_ionization_compare.py [ionz_on.txt] [ionz_off.txt] [out.mp4] [fps]

Defaults look for ../out_ioniz_on.txt and ../out_ioniz_off.txt at the
project root, and write util/model_c7_ionization_compare.mp4 at 15 fps.

The two runs may have different frame counts (different CFL-limited dt).
The animation is driven by the ionization-on run's frame times; for each
of those times the ionization-off run is sampled at its nearest frame.
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
    arrs = [np.concatenate(s) for s in series_list]
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
    root = os.path.abspath(os.path.join(here, ".."))
    on_path  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, "out_ioniz_on.txt")
    off_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "out_ioniz_off.txt")
    out_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(here, "visualization", "model_c7_ionization_compare.mp4")
    os.makedirs(os.path.join(here, "visualization"), exist_ok=True)
    fps      = int(sys.argv[4]) if len(sys.argv) > 4 else 15

    xx_on, frames_on   = read_frames(on_path)
    xx_off, frames_off = read_frames(off_path)
    assert xx_on.size == xx_off.size and np.allclose(xx_on, xx_off), "grid mismatch between runs"
    xx = xx_on

    prims_on  = [primitives(fr[2]) for fr in frames_on]
    prims_off = [primitives(fr[2]) for fr in frames_off]
    t_on  = np.array([fr[0] for fr in frames_on])
    t_off = np.array([fr[0] for fr in frames_off])

    # Drive animation on ionization-on timeline; pick nearest off-frame per tick.
    off_idx_for_on = np.array([int(np.argmin(np.abs(t_off - t))) for t in t_on])

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

    lines_on  = {}
    lines_off = {}
    for fld in fields:
        r, c = field_to_panel[fld]
        ax = axes[r, c]
        ipr = idx_in_prim[fld]
        series_on  = [p[ipr] for p in prims_on]
        series_off = [p[ipr] for p in prims_off]
        lo, hi = panel_limits([series_on, series_off], log=log_yscale.get(fld, False))
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(lo, hi)
        if log_yscale.get(fld, False):
            ax.set_yscale("log")
        ax.set_title(labels[fld])
        ax.grid(True, alpha=0.3)
        (l_on,)  = ax.plot(xx, series_on[0],  lw=1.4, color="C0", label="ionization ON (Stage E)")
        (l_off,) = ax.plot(xx, series_off[0], lw=1.4, color="C3", linestyle="--", label="ionization OFF")
        lines_on[fld]  = l_on
        lines_off[fld] = l_off

    for ax in axes[1, :]:
        ax.set_xlabel("height (km)")

    handles = [lines_on["nn"], lines_off["nn"]]
    fig.legend(handles, [h.get_label() for h in handles],
               loc="upper left", bbox_to_anchor=(0.01, 0.985), ncol=2, frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    def update(k):
        kon  = k
        koff = int(off_idx_for_on[k])
        pon  = prims_on[kon]
        poff = prims_off[koff]
        for fld in fields:
            ipr = idx_in_prim[fld]
            lines_on[fld].set_ydata(pon[ipr])
            lines_off[fld].set_ydata(poff[ipr])
        title.set_text(
            f"Model C7  |  ionization ON  t = {t_on[kon]:6.2f} s   vs   ionization OFF  t = {t_off[koff]:6.2f} s   ({k+1}/{len(t_on)})"
        )
        return list(lines_on.values()) + list(lines_off.values()) + [title]

    anim = animation.FuncAnimation(fig, update, frames=len(t_on), interval=1000.0/fps, blit=False)
    writer = animation.FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    anim.save(out_path, writer=writer, dpi=120)
    print(f"wrote {out_path}  ({len(t_on)} frames @ {fps} fps = {len(t_on)/fps:.1f} s)")
    print(f"ionization-on  : {len(t_on)} frames, t in [{t_on.min():.1f}, {t_on.max():.1f}] s")
    print(f"ionization-off : {len(t_off)} frames, t in [{t_off.min():.1f}, {t_off.max():.1f}] s")


if __name__ == "__main__":
    main()
