#!/usr/bin/env python3
"""Overlay two chromo_main runs (ionization on vs off) on the same 8-panel animation.

Usage:
    python3 animate_ionization_compare.py [ionz_on.txt] [ionz_off.txt] [out.mp4] [fps]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_output import read_frames, primitives
from _anim_parallel import save_frames_parallel


# (row, col, prims_index, panel_title, log_yscale)
_PANELS = [
    (0, 0, 1, r"$n_n$ (m$^{-3}$)", True),
    (0, 1, 3, r"$u$ (m/s)",        False),
    (0, 2, 5, r"$P_n$ (Pa)",       True),
    (0, 3, 7, r"$T_n$ (K)",        False),
    (1, 0, 0, r"$n_i$ (m$^{-3}$)", True),
    (1, 1, 2, r"$v$ (m/s)",        False),
    (1, 2, 4, r"$P_i$ (Pa)",       True),
    (1, 3, 6, r"$T_i$ (K)",        False),
]


def _render_frame(k, tmpdir, xx, pon_k, t_on, poff_k, t_off, nF, ylims):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    fig.suptitle(
        f"Model C7  |  ionization ON  t = {t_on:6.2f} s"
        f"   vs   ionization OFF  t = {t_off:6.2f} s   ({k+1}/{nF})"
    )
    for row, col, idx, title, log in _PANELS:
        ax = axes[row, col]
        lo, hi = ylims[title]
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(lo, hi)
        if log:
            ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.plot(xx, pon_k[idx],  lw=1.4, color="C0", label="ionization ON (Stage E)")
        ax.plot(xx, poff_k[idx], lw=1.4, color="C3", linestyle="--", label="ionization OFF")
    for ax in axes[1, :]:
        ax.set_xlabel("height (km)")

    handles = [axes[0, 0].lines[0], axes[0, 0].lines[1]]
    fig.legend(handles, [h.get_label() for h in handles],
               loc="upper left", bbox_to_anchor=(0.01, 0.985), ncol=2, frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def main():
    root = os.path.abspath(os.path.join(HERE, ".."))
    on_path  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, "out_ioniz_on.txt")
    off_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "out_ioniz_off.txt")
    out_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, "visualization", "model_c7_ionization_compare.mp4")
    fps      = int(sys.argv[4]) if len(sys.argv) > 4 else 15
    os.makedirs(os.path.join(HERE, "visualization"), exist_ok=True)

    xx_on,  frames_on  = read_frames(on_path)
    xx_off, frames_off = read_frames(off_path)
    assert xx_on.size == xx_off.size and np.allclose(xx_on, xx_off), "grid mismatch"
    xx = xx_on

    prims_on  = [primitives(fr[2]) for fr in frames_on]
    prims_off = [primitives(fr[2]) for fr in frames_off]
    t_on  = np.array([fr[0] for fr in frames_on])
    t_off = np.array([fr[0] for fr in frames_off])
    nF = len(t_on)

    off_idx = np.array([int(np.argmin(np.abs(t_off - t))) for t in t_on])

    def lim(series_list, log=False):
        a = np.concatenate([np.concatenate(s) for s in series_list])
        if log:
            a = a[a > 0]
        lo, hi = float(np.min(a)), float(np.max(a))
        if log:
            return lo / 1.2, hi * 1.2
        pad = 0.05 * (hi - lo if hi > lo else max(abs(hi), 1.0))
        return lo - pad, hi + pad

    ylims = {
        title: lim([[p[idx] for p in prims_on], [p[idx] for p in prims_off]], log=log)
        for _, _, idx, title, log in _PANELS
    }

    print(f"ionization-on : {len(t_on)} frames, t in [{t_on.min():.1f}, {t_on.max():.1f}] s")
    print(f"ionization-off: {len(t_off)} frames, t in [{t_off.min():.1f}, {t_off.max():.1f}] s")

    save_frames_parallel(
        _render_frame,
        [(k, xx, prims_on[k], t_on[k], prims_off[off_idx[k]], t_off[off_idx[k]], nF, ylims)
         for k in range(nF)],
        out_path, fps,
    )


if __name__ == "__main__":
    main()
