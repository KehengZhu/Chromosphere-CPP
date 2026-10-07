#!/usr/bin/env python3
"""Is the upper-chromosphere unsteadiness of the hot, high-pressure corner numerical?

Compares the (70 kK, 0.08 Pa) case of the (T_wall, p_top) sweep over its first
T_END seconds with three 1500 s variants of it (see
docs/studies/boundaries/tp_sweep_evaporation_closure_validation.md):

  base       the sweep run itself, N=500/R4, CFL 0.50, grade at 2079-2101 km
  n1000      N=1000/R4 (every cell halved)
  cfl025     CHROMO_CFL=0.25
  grade1950  ISO_REFINE_S_LO_KM=350: the fine band starts at 1950 km, not 2100 km

A numerical origin shows up as unsteadiness that changes with the resolution or
the CFL, or that moves with the mesh grade; a physical one is insensitive to all
three. Per case, over [T0, T_END]:
  - max_z std_t F(z) / S and its height, with S = median_z |mean_t F(z)| one flux
    scale for the column (a local mean can pass through zero and blow a local ratio up);
  - the fraction of top-face captures with F <= 0;
  - the mesh grade band and the TR base (7.5 kK) height.

Writes <corner-dir>/corner_test.json and the figure
<out-dir>/tp_corner_test.png (unsteadiness profile, top-face flux, and the
space-time map of F(z, t)/mean_t F(z) for every case).

Usage:
    python util/tp_corner_test_compare.py [--sweep-dir outputs/model_column/tp_sweep]
        [--corner-dir outputs/model_column/tp_sweep_corner]
        [--out-dir visualization/model_column] [--t0 500] [--t-end 1500]
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from animate_isentropic import load_face_flux
from twall_sweep_diag import COLUMNS, load_gamma_diag

ROOT = os.path.dirname(HERE)
I_T = COLUMNS.index("T")
T_MARK = 7500.0                # TR-base marker isotherm [K]


def cases(sweep_dir, corner_dir):
    return [
        ("base (N500, CFL 0.50)", os.path.join(sweep_dir, "tp_T70kK_p0.08_N500.txt")),
        ("N1000", os.path.join(corner_dir, "n1000", "tp_T70kK_p0.08_N1000.txt")),
        ("CFL 0.25", os.path.join(corner_dir, "cfl025", "tp_T70kK_p0.08_N500.txt")),
        ("grade at 1950 km", os.path.join(corner_dir, "grade1950", "tp_T70kK_p0.08_N500.txt")),
    ]


def grade_band(face_km):
    """Heights [km] spanned by the cells graded between the coarse and fine spacing."""
    w = np.diff(face_km)
    fine, coarse = w.min(), w.max()
    g = np.where((w > 1.05 * fine) & (w < 0.95 * coarse))[0]
    return (float(face_km[g[0]]), float(face_km[g[-1] + 1])) if len(g) else (np.nan, np.nan)


def analyse(label, path, t0, t_end):
    fh, ft, ff = load_face_flux(path + ".faceflux")
    fh, ft, ff = np.asarray(fh), np.asarray(ft), np.asarray(ff)
    keep = ft <= t_end + 1e-6
    ft, ff = ft[keep], ff[keep]
    w = ft >= t0
    mean, std = ff[w].mean(0), ff[w].std(0)
    scale = float(np.median(np.abs(mean)))       # one flux scale for the whole column
    rel = std / scale
    h, tg, d = load_gamma_diag(path)
    h = np.asarray(h)
    lw = (tg >= t0) & (tg <= t_end)
    tr = [h[np.argmax(T > T_MARK)] for T in d[lw, :, I_T]]
    return dict(
        label=label, path=path, n_faces=len(fh), n_captures=int(w.sum()),
        max_relstd=float(rel.max()), max_relstd_height_km=float(fh[rel.argmax()]),
        base_relstd=float(rel[0]),
        top_reversal_fraction=float(np.mean(ff[w, -1] <= 0.0)),
        flux_scale=scale, min_face_flux_over_scale=float(ff[w].min() / scale),
        grade_band_km=grade_band(fh), tr_base_km=[float(min(tr)), float(max(tr))],
        _fh=fh, _ft=ft, _ff=ff, _rel=rel)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sweep-dir", default=os.path.join(ROOT, "outputs/model_column/tp_sweep"))
    ap.add_argument("--corner-dir", default=os.path.join(ROOT, "outputs/model_column/tp_sweep_corner"))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "visualization/model_column"))
    ap.add_argument("--t0", type=float, default=500.0)
    ap.add_argument("--t-end", type=float, default=1500.0)
    a = ap.parse_args()
    res = [analyse(lab, p, a.t0, a.t_end) for lab, p in cases(a.sweep_dir, a.corner_dir)]

    print(f"{'case':24s} {'faces':>5} {'max relstd':>10} {'at km':>7} {'base':>5} "
          f"{'top F<=0':>8} {'min F/scl':>10}  grade [km]      TR base [km]")
    for r in res:
        print(f"{r['label']:24s} {r['n_faces']:5d} {r['max_relstd']:10.2f} "
              f"{r['max_relstd_height_km']:7.1f} {r['base_relstd']:5.2f} "
              f"{r['top_reversal_fraction']:8.3f} {r['min_face_flux_over_scale']:10.1f}  "
              f"{r['grade_band_km'][0]:.0f}-{r['grade_band_km'][1]:.0f}      "
              f"{r['tr_base_km'][0]:.1f}-{r['tr_base_km'][1]:.1f}")
    json.dump(dict(t0_s=a.t0, t_end_s=a.t_end,
                   cases=[{k: v for k, v in r.items() if not k.startswith("_")} for r in res]),
              open(os.path.join(a.corner_dir, "corner_test.json"), "w"), indent=1)

    colors = ["k", "C0", "C1", "C3"]
    fig = plt.figure(figsize=(17, 9.5))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.15])
    ax0 = fig.add_subplot(gs[0, :2])
    ax1 = fig.add_subplot(gs[0, 2:])
    for r, c in zip(res, colors):
        ax0.plot(r["_fh"], r["_rel"], color=c, lw=1.1, label=r["label"])
        g0, g1 = r["grade_band_km"]
        ax0.axvspan(g0, g1, color=c, alpha=0.08)
        ax0.axvline(r["tr_base_km"][1], color=c, ls=":", lw=0.9)
        ft, Ft = r["_ft"], r["_ff"][:, -1]
        ax1.plot(ft, Ft / np.median(Ft[ft >= a.t0]), color=c, lw=0.8, label=r["label"])
    ax0.set(yscale="log", xlabel="face height [km]",
            ylabel=rf"std$_t F(z)$ / median$_z$|mean$_t F|$, {a.t0:.0f}-{a.t_end:.0f} s",
            title="(a) where the flux is unsteady; shading = mesh grade, dotted = TR base (7.5 kK)")
    ax0.legend(fontsize=8)
    ax1.axhline(0.0, color="0.5", lw=0.6)
    ax1.set(xlim=(0, a.t_end), ylim=(-4, 6), xlabel="t [s]",
            ylabel=r"$F_{top}(t)$ / median$_{t \geq t_0}$",
            title="(b) top-face mass flux")
    ax1.legend(fontsize=8, loc="upper right")
    norm = TwoSlopeNorm(1.0, -3.0, 5.0)
    for j, r in enumerate(res):
        ax = fig.add_subplot(gs[1, j])
        sel = r["_fh"] >= 1900.0
        pc = ax.pcolormesh(r["_ft"], r["_fh"][sel], (r["_ff"][:, sel] / r["flux_scale"]).T,
                           norm=norm, cmap="RdBu_r", shading="auto", rasterized=True)
        g0, g1 = r["grade_band_km"]
        for g in (g0, g1):
            ax.axhline(g, color="k", lw=0.6, ls="--")
        ax.set(title=f"({'cdef'[j]}) {r['label']}", xlabel="t [s]",
               ylabel="face height [km]" if j == 0 else None)
    fig.colorbar(pc, ax=fig.axes[-4:], fraction=0.02, label=r"$F(z,t)$ / median$_z$|mean$_t F$|")
    fig.suptitle("model_column, T_wall = 70 kK, p_top = 0.08 Pa: sensitivity of the "
                 "upper-chromosphere unsteadiness to resolution, CFL and the mesh grade "
                 "(dashed: grade band)", fontsize=10)
    out = os.path.join(a.out_dir, "tp_corner_test.png")
    os.makedirs(a.out_dir, exist_ok=True)
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
