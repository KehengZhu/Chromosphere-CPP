#!/usr/bin/env python3
"""Plot the last (or a chosen) frame of build/output.txt as an 8-panel snapshot.

Output.txt format (written by chromo_main.cpp):
    line 1: "ns num_of_eq"
    line 2: ns cell heights (km)
    repeated: "# t = T step = S" marker, then ns lines of num_of_eq values
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m_i = 1.6726219e-27
m_n = m_i
k_b = 1.380649e-23

CNI, CNN, CNV, CNU, CEI, CEN = 0, 1, 2, 3, 4, 5


def read_frames(path):
    """Return (heights, frames) where frames is a list of (t, step, xn[ns,6])."""
    with open(path) as f:
        raw = [ln for ln in f.read().splitlines() if ln.strip()]

    hdr = raw[0].split()
    ns, neq = int(hdr[0]), int(hdr[1])
    heights = np.fromstring(raw[1], sep=" ")
    assert heights.size == ns, f"expected {ns} heights, got {heights.size}"

    frames = []
    i = 2
    while i < len(raw):
        line = raw[i]
        if not line.startswith("#"):
            raise ValueError(f"expected '# t = ...' at line {i+1}, got: {line!r}")
        toks = line.split()
        # "# t = T step = S"
        t   = float(toks[3])
        stp = int(toks[6])
        i += 1
        data = np.zeros((ns, neq))
        for j in range(ns):
            data[j, :] = np.fromstring(raw[i + j], sep=" ")
        i += ns
        frames.append((t, stp, data))
    return heights, frames


def primitives(xn):
    ni = xn[:, CNI] / m_i
    nn = xn[:, CNN] / m_n
    v  = xn[:, CNV] / (m_i * ni)
    u  = xn[:, CNU] / (m_n * nn)
    p_i = 2.0/3.0 * xn[:, CEI] - 1.0/3.0 * m_i * ni * v * v
    p_n = 2.0/3.0 * xn[:, CEN] - 1.0/3.0 * m_n * nn * u * u
    Ti = p_i / (2.0 * ni * k_b)
    Tn = p_n / (nn * k_b)
    return ni, nn, v, u, p_i, p_n, Ti, Tn


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    in_path  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "build", "output.txt")
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "model_c7_snapshot.png")

    xx, frames = read_frames(in_path)
    t_final, step_final, xn = frames[-1]
    ni, nn, v, u, p_i, p_n, Ti, Tn = primitives(xn)

    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    fig.suptitle(f"Chromosphere — Model C7 final state (t = {t_final:.1f} s, step {step_final})")

    panels = [
        (axes[0, 0], nn,  r"$n_n$ (m$^{-3}$)", "log"),
        (axes[0, 1], u,   r"$u$ (m/s)",        None),
        (axes[0, 2], p_n, r"$P_n$ (Pa)",       "log"),
        (axes[0, 3], Tn,  r"$T_n$ (K)",        None),
        (axes[1, 0], ni,  r"$n_i$ (m$^{-3}$)", "log"),
        (axes[1, 1], v,   r"$v$ (m/s)",        None),
        (axes[1, 2], p_i, r"$P_i$ (Pa)",       "log"),
        (axes[1, 3], Ti,  r"$T_i$ (K)",        None),
    ]
    for ax, y, title, yscale in panels:
        ax.plot(xx, y, lw=1.2)
        ax.set_title(title)
        if yscale:
            ax.set_yscale(yscale)
        ax.grid(True, alpha=0.3)
    for ax in axes[1, :]:
        ax.set_xlabel("height (km)")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=130)
    print(f"wrote {out_path}")
    print(f"frames in file: {len(frames)}  (final t = {t_final:.2f} s, step {step_final})")
    print(f"height range: [{xx.min():.1f}, {xx.max():.1f}] km")
    print(f"final  T_i: [{Ti.min():.1f}, {Ti.max():.1f}] K")
    print(f"final  T_n: [{Tn.min():.1f}, {Tn.max():.1f}] K")
    print(f"final |v|max = {np.max(np.abs(v)):.3e}, |u|max = {np.max(np.abs(u)):.3e} m/s")


if __name__ == "__main__":
    main()
