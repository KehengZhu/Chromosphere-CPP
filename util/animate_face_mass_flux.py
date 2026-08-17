#!/usr/bin/env python3
r"""Cell-centred rho*V vs the finite-volume face mass flux for a model_column run.

Read-only diagnostic movie. Nothing here re-derives the numerics: every series is
copied out of the `<run>.faceflux` sidecar, which `chromo_main` writes with
CHROMO_FACE_FLUX_DIAG=1 straight out of the production `mixture_rhs_explicit`
call (release single-fluid solver only).

Quantities plotted, with face index i = the UPPER face i+1/2 of cell i:

    rhoV_cell[i] = rho_cell[i] * v_cell[i]     cell-centred mass flux    [cell_km]
    f_total[i]                                 production numerical face mass flux
                                               (the release Godunov flux, or the
                                               Roe / Rusanov reference solver) -- the
                                               thing the continuity row differences
                                               [face_km]
    f_ref[i]     = f_total[i] at t = 0         t = 0 baseline flux
    f_eff[i]     = f_total[i] - f_ref[i]       CHANGE in the transport since t = 0

Why f_eff is worth plotting: the solver differences `f_total` and nothing else,
but f_ref is exactly the O(ds) truncation-error flux the discretisation carries
at rest -- it is nonzero even though the IC has V == 0 everywhere, which is
visible in frame 1: rhoV = 0 and f_ref = f_total != 0 while f_eff == 0.
Subtracting that fixed offset is therefore the cleanest way to see how the
transport has evolved away from the initial condition. It is a DIAGNOSTIC
baseline, not a term in the update.

Layout (2x2, shared y-scale per column so the two rows are directly comparable):

    A  cell rhoV, full column        B  face f_total / f_ref / f_eff, full column
    C  cell rhoV, top zoom           D  face f_total / f_ref / f_eff, top zoom

Panels B/D also carry rhoV as a faint grey line so the cell-vs-face comparison
needs no eye-tracking between rows. All four series come from the SAME sidecar
record, i.e. the same production RHS evaluation at the same step, so the frames
are exactly synchronous by construction (see --check-output below for the
cadence caveat on the main snapshot file).

Usage:
  python util/animate_face_mass_flux.py <run.txt> <out.mp4> [--fps 20] \
      [--zoom LO_KM HI_KM] [--window LO_KM HI_KM] [--stride N] [--max-frames N] \
      [--faceflux PATH] [--final-figure PATH | --no-final-figure] [--no-check-output]

Env (matching the sibling animation scripts):
  ANIM_STRIDE      take every Nth sidecar record (default 1)
  ANIM_MAX_FRAMES  frame budget, evenly subsampled (default: all)
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _anim_parallel import save_frames_parallel                    # noqa: E402
from face_flux_diag import read_faceflux, roughness                # noqa: E402

C_RHOV = "#2f3640"      # neutral dark  -- cell-centred rho*V
C_TOTAL = "#1f6fb2"     # signal blue   -- f_total
C_REF = "#7f8c8d"       # grey dashed   -- f_ref (t=0 baseline, time-independent)
C_EFF = "#c0392b"       # accent red    -- f_eff = f_total - f_ref
C_GHOST = "#b2bec3"     # faint         -- rhoV echoed in the face panels


def _apply_style():
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "figure.dpi": 130,
    })


def check_reference(ref):
    """Print the sanity check that makes the captured face flux trustworthy."""
    split = np.max(np.abs(ref["f_central"] + ref["f_diff"] - ref["f_total"])
                   / np.maximum(np.abs(ref["f_total"]), 1e-300))
    print(f"  split f_central+f_diff == f_total : max rel err {split:.2e}")
    print("  f_ref is the t=0 face flux: a diagnostic baseline only -- the "
          "solver differences f_total itself.")


def check_main_output(run_path, records, ns):
    """Cross-check rhoV against the main snapshot file (cadences differ)."""
    if not os.path.exists(run_path):
        print(f"  main output {run_path} not present -- skipping rhoV cross-check")
        return
    with open(run_path) as fh:
        raw = [ln for ln in fh.read().splitlines() if ln.strip()]
    n_out, rows = (int(x) for x in raw[0].split())
    if n_out != ns or rows != 3:
        print(f"  main output is {n_out} cells x {rows} rows (sidecar has {ns} "
              f"cells, release state has 3) -- skipping rhoV cross-check")
        return
    frames, i = [], 2
    while i < len(raw):
        if not raw[i].startswith("# t ="):
            i += 1
            continue
        t = float(raw[i].split()[3])
        i += 1
        if i + n_out > len(raw):
            break
        blk = np.array([np.fromstring(raw[i + j], sep=" ") for j in range(n_out)])
        i += n_out
        frames.append((t, blk[:, 1]))          # column 1 = mix::MOM = rho*V
    if not frames:
        print(f"  main output {run_path} holds no frames -- skipping cross-check")
        return
    times = np.array([t for t, _ in frames])
    for label, (t_s, _step, c) in (("first", records[0]), ("last", records[-1])):
        j = int(np.argmin(np.abs(times - t_s)))
        mom = frames[j][1]
        rhoV = c["rho_cell"] * c["v_cell"]
        den = max(np.max(np.abs(mom)), 1e-300)
        print(f"  rhoV cross-check ({label:5s}): sidecar t={t_s:9.3f} s vs "
              f"snapshot t={times[j]:9.3f} s (dt={times[j] - t_s:+.3f} s), "
              f"max|rho_cell*v_cell - mom|/max|mom| = "
              f"{np.max(np.abs(rhoV - mom)) / den:.2e}")


def pick_scale(maxabs):
    """Common y multiplier -> (scale, exponent) with the label reading x10^-exp."""
    if not np.isfinite(maxabs) or maxabs <= 0.0:
        return 1.0, 0
    exp = int(np.floor(np.log10(maxabs)))
    return 10.0 ** (-exp), exp


def build_figure(cell_km, face_km, rhoV, f_total, f_ref, f_eff, t, step,
                 full_x, zoom, window, ylims, exp, label):
    import matplotlib.pyplot as plt

    _apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 6.8), sharex="col",
                             gridspec_kw={"width_ratios": [1.45, 1.0]})

    for col, (xlo, xhi) in enumerate((full_x, zoom)):
        ax_c, ax_f = axes[0, col], axes[1, col]
        ylo, yhi = ylims[col]
        for ax in (ax_c, ax_f):
            ax.axvspan(*window, color="#f0c419", alpha=0.13, lw=0, zorder=0)
            ax.axhline(0.0, color="#95a5a6", lw=0.6, zorder=1)
            ax.set_xlim(xlo, xhi)
            ax.set_ylim(ylo, yhi)

        # cell-centred rho*V, plotted against the CELL CENTRES
        ax_c.plot(cell_km, rhoV, color=C_RHOV, lw=1.1)

        # face fluxes, plotted against the FACE positions. rho*V is REPEATED here
        # as a thick halo behind f_eff -- same array as the panel above, same
        # y-scale: where the two agree it shows as a grey fringe around the red
        # line, where they part the separation is the cell-vs-face difference.
        ax_f.plot(cell_km, rhoV, color=C_GHOST, lw=3.4, solid_capstyle="round",
                  zorder=2, label=r"cell $\rho V$ (repeat of top panel)")
        ax_f.plot(face_km, f_ref, color=C_REF, lw=1.0, ls="--", zorder=3,
                  label=r"$f_{\rm ref}$  ($t=0$ baseline)")
        ax_f.plot(face_km, f_total, color=C_TOTAL, lw=1.1, zorder=4,
                  label=r"$f_{\rm total}$")
        ax_f.plot(face_km, f_eff, color=C_EFF, lw=1.3, zorder=5,
                  label=r"$f_{\rm eff}=f_{\rm total}-f_{\rm ref}$")
        ax_f.legend(fontsize=7.5, loc="upper left", ncol=2, framealpha=0.85)
        ax_f.set_xlabel("height  [km]")

    unit = (rf"[$10^{{{exp}}}$ kg m$^{{-2}}$ s$^{{-1}}$]" if exp
            else r"[kg m$^{-2}$ s$^{-1}$]")
    axes[0, 0].set_ylabel(rf"cell-centred $\rho V$  {unit}")
    axes[1, 0].set_ylabel(rf"face mass flux  {unit}")

    axes[0, 0].set_title(f"{label}   full column", fontsize=9, loc="left",
                         color="#57606f")
    axes[0, 1].set_title(f"zoom {zoom[0]:.0f}-{zoom[1]:.0f} km   "
                         f"(shaded = {window[0]:.0f}-{window[1]:.0f} km window)",
                         fontsize=8, loc="left", color="#57606f")
    axes[0, 0].text(0.985, 0.92, f"t = {t:8.1f} s", transform=axes[0, 0].transAxes,
                    ha="right", va="top", fontsize=12, color="#2f3640",
                    family="monospace")

    m = (cell_km >= window[0]) & (cell_km <= window[1])
    mf = (face_km >= window[0]) & (face_km <= window[1])
    r_cell, r_eff = roughness(rhoV[m]), roughness(f_eff[mf])
    ratio = r_cell / r_eff if r_eff > 0 else float("nan")
    axes[1, 1].text(
        0.985, 0.04,
        f"step {step}\n"
        f"{window[0]:.0f}-{window[1]:.0f} km:\n"
        rf"  mean $\rho V$   = {np.mean(rhoV[m]):+.3e}" "\n"
        rf"  mean $f_{{\rm eff}}$  = {np.mean(f_eff[mf]):+.3e}" "\n"
        rf"  rough $\rho V$/$f_{{\rm eff}}$ = {ratio:.0f}$\times$",
        transform=axes[1, 1].transAxes, ha="right", va="bottom", fontsize=7,
        family="monospace", color="#57606f",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#dfe4ea", alpha=0.9))

    fig.tight_layout(h_pad=0.7, w_pad=1.4)
    return fig


def _render(k, tmpdir, *frame_args):
    import matplotlib.pyplot as plt

    fig = build_figure(*frame_args)
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"))
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="model_column output file (.txt); the sidecar is "
                               "<run>.faceflux unless --faceflux is given")
    ap.add_argument("out", help="output .mp4")
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--faceflux", default=None)
    ap.add_argument("--zoom", nargs=2, type=float, default=(2100.0, None),
                    metavar=("LO_KM", "HI_KM"),
                    help="top-region zoom; HI defaults to the domain top")
    ap.add_argument("--window", nargs=2, type=float, default=(2130.0, 2150.0),
                    metavar=("LO_KM", "HI_KM"),
                    help="shaded release analysis window")
    ap.add_argument("--stride", type=int,
                    default=int(os.environ.get("ANIM_STRIDE", 1)))
    ap.add_argument("--max-frames", type=int,
                    default=int(os.environ.get("ANIM_MAX_FRAMES", 0)) or None)
    ap.add_argument("--final-figure", default=None,
                    help="static last-frame figure (default: <out>_final.png)")
    ap.add_argument("--no-final-figure", action="store_true")
    ap.add_argument("--no-check-output", action="store_true",
                    help="skip the rhoV cross-check against the main snapshot file")
    args = ap.parse_args(argv)

    sidecar = args.faceflux or (args.run + ".faceflux")
    if not os.path.exists(sidecar):
        sys.exit(f"missing face-flux sidecar: {sidecar}\n"
                 f"Rerun the case with CHROMO_FACE_FLUX_DIAG=1 (release solver "
                 f"only) to produce it.")

    meta, records = read_faceflux(sidecar)
    if len(records) < 2:
        sys.exit(f"{sidecar}: only {len(records)} record(s); need the t=0 reference "
                 f"record plus at least one later record")
    if records[0][0] != 0.0:
        sys.exit(f"{sidecar}: first record is at t = {records[0][0]} s, not 0. "
                 f"f_ref is defined as the t=0 face flux, so this sidecar cannot "
                 f"be used (capture the run from its start).")

    ref = records[0][2]
    cell_km, face_km = ref["cell_km"], ref["face_km"]
    print(f"{sidecar}\n  "
          + " ".join(f"{k}={v}" for k, v in meta.items() if k and v))
    print(f"  {len(records)} records, t = {records[0][0]:.2f} -> "
          f"{records[-1][0]:.2f} s, {len(cell_km)} captured cells "
          f"({cell_km[0]:.1f} - {face_km[-1]:.1f} km)")
    check_reference(ref)
    if not args.no_check_output:
        check_main_output(args.run, records, len(cell_km))

    f_ref = ref["f_total"]
    frames = records[::max(1, args.stride)]
    if args.max_frames and len(frames) > args.max_frames:
        sel = np.linspace(0, len(frames) - 1, args.max_frames).astype(int)
        frames = [frames[i] for i in sel]
    print(f"  animating {len(frames)} frames @ {args.fps} fps")

    rhoV = [c["rho_cell"] * c["v_cell"] for _, _, c in frames]
    f_tot = [c["f_total"] for _, _, c in frames]
    f_eff = [ft - f_ref for ft in f_tot]

    # Clamp the zoom to what was captured (CHROMO_FACE_FLUX_TOP can truncate the
    # sidecar to the top N cells, in which case the default zoom is the whole file).
    full_x = (float(cell_km[0]), float(face_km[-1]))
    zoom = (max(args.zoom[0], full_x[0]),
            args.zoom[1] if args.zoom[1] is not None else full_x[1])

    # Fixed limits over the whole movie, shared by both rows of a column so that
    # cell rhoV and the face fluxes are read on the SAME scale. The zoom column
    # gets its own limits so late evaporation structure is not squashed by the
    # early transient elsewhere in the column.
    def lim(mask_c, mask_f, pad=0.06):
        vals = ([a[mask_c] for a in rhoV]
                + [a[mask_f] for a in f_tot] + [a[mask_f] for a in f_eff]
                + [f_ref[mask_f]])
        lo = min(float(a.min()) for a in vals)
        hi = max(float(a.max()) for a in vals)
        span = hi - lo or max(abs(hi), 1e-300)
        return lo - pad * span, hi + pad * span

    all_c = np.ones_like(cell_km, dtype=bool)
    zm_c = (cell_km >= zoom[0]) & (cell_km <= zoom[1])
    zm_f = (face_km >= zoom[0]) & (face_km <= zoom[1])
    raw_lims = [lim(all_c, all_c), lim(zm_c, zm_f)]
    scale, exp = pick_scale(max(abs(v) for pair in raw_lims for v in pair))
    ylims = [(lo * scale, hi * scale) for lo, hi in raw_lims]
    label = os.path.basename(args.run)

    frame_args = [
        (k, cell_km, face_km, rhoV[k] * scale, f_tot[k] * scale, f_ref * scale,
         f_eff[k] * scale, frames[k][0], frames[k][1], full_x, zoom,
         tuple(args.window), ylims, exp, label)
        for k in range(len(frames))
    ]
    save_frames_parallel(_render, frame_args, args.out, args.fps)

    if not args.no_final_figure:
        import matplotlib.pyplot as plt

        path = args.final_figure or (os.path.splitext(args.out)[0] + "_final.png")
        fig = build_figure(*frame_args[-1][1:])
        fig.savefig(path)
        plt.close(fig)
        print(f"wrote {path}  (static last frame, t = {frames[-1][0]:.2f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
