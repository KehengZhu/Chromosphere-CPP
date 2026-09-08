#!/usr/bin/env python3
"""Animate the upper-boundary wall-temperature sweep as one overlaid 6-panel movie.

Every accepted run of ``outputs/model_column/twall_sweep/`` is drawn on the SAME
axes, one coloured line per imposed conduction-wall temperature ``T_wall``, in the
six panels of the release evolution movie
``visualization/model_column/lnp_godunov_N500_4000s_f64_evolution.mp4``:

    T(h)      V(h)      q_parallel(h)
    rho(h)    p(h)      mass flux(h)

so the sweep can be read against the single-run release movie panel for panel.
Run discovery, run acceptance and the per-run colours are the same as in
``util/twall_sweep_diag.py`` (imported, not re-implemented), so the movie and the
four static sweep figures always describe the same set of runs.

MASS FLUX PANEL. As in ``util/animate_isentropic.py``, the panel shows the
CONSERVATIVE NUMERICAL FACE MASS FLUX ``f_total`` from the ``<run>.faceflux``
sidecar whenever every run has one, plotted at the face heights. That is the
quantity the continuity row is actually differenced from; the cell-centred product
rho*V is only a reconstruction of it, carrying the full cell-centred ripple and
systematically misleading in the boundary cells, where the ghost closure fixes the
face flux and not the centred product. Produce the sidecars with
``SWEEP_FACE_FLUX=1 util/run_twall_sweep.sh`` (a read-only diagnostic: the
``.gamma_diag`` output is byte-identical with and without it). Set
ANIM_MASS_FLUX=cell to force rho*V, or =face to require the sidecars. The panel
title always says which one is being drawn.

PRESSURE PANEL. The panel draws p RELATIVE TO THE 22 kK RELEASE BASELINE,
p(h,t)/p_22kK(h,t), not p in Pa. Absolute p cannot separate the sweep at any y
scale: p(h) falls 6.6x over the column while the wall-to-wall spread is at most
12 %, so the 13 lines sit within one line width of each other everywhere, and
zooming y also crops x (p drops 12 % in ~30 km), which leaves the fan packed
against the top of the domain instead of spreading it. Dividing by the release run
removes the shared stratification and the shared early acoustic transient, leaving
a 0.90-1.10 fan that is separated over the whole 1600-2153 km. The release run is
then the flat unit reference line rather than a curve. Set ANIM_P_PANEL=abs to get
the old absolute log-Pa panel back; the panel title always says which one is drawn.

ONE DELIBERATE DEVIATION from the single-run movie, forced by overlaying 13 runs
whose responses span four decades: **V, q and the mass flux use symlog axes, and V
is in m/s, not km/s.** The 8 kK wall drains at -0.4 m/s while the 100 kK wall
evaporates at +792 m/s; on the single-run movie's linear km/s axis every wall below
45 kK is a flat line on zero. T and rho are logarithmic for the same reason
(T_top spans 6.7 -> 100 kK).

Runs are keyed by a DISCRETE colorbar with one equal-width band per run, built from
the actual line colours — not a continuous ``LogNorm`` bar, which would show colours
no line has (``wall_colors`` spaces colours by rank and caps their luminance). The
22 kK release baseline is drawn thicker and its tick label is bold.

Usage:
    python util/animate_twall_sweep.py [--dir DIR] [--cfl 0.50]
        [--run T_WALL_K=path] [--out FILE.mp4] [--fps 20]

Defaults reproduce the catalogued figure. Set ANIM_MAX_FRAMES to cap the number of
rendered frames (the run is subsampled uniformly in time).
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _anim_parallel import save_frames_parallel
from animate_isentropic import load_face_flux
from conductive_flux import physical_heat_flux
from twall_sweep_diag import (COLUMNS, IRHO, IT, IP, IV,
                              console_ok, load_gamma_diag, wall_colors)

I_NE = COLUMNS.index("n_e")
I_NHI = COLUMNS.index("n_HI")
T_WALL_RELEASE = 22.0e3      # the release baseline, drawn thicker


def collect_entries(sweep_dir, cfl, explicit):
    """-> sorted [(T_wall_K, path)], the same discovery rule as twall_sweep_diag."""
    entries = []
    for path in sorted(glob.glob(os.path.join(sweep_dir, "twall_*kK_cfl*.txt"))):
        mm = re.search(r"twall_(\d+)kK_cfl([0-9.]+)\.txt$", os.path.basename(path))
        if not mm:
            continue
        if cfl != "all" and mm.group(2) != cfl:
            continue
        entries.append((float(mm.group(1)) * 1.0e3, path))
    for spec in explicit:
        t_wall, path = spec.split("=", 1)
        entries.append((float(t_wall), path))
    entries.sort()
    return entries


def load_run(t_wall, path):
    """-> dict of the six panel fields for one run, or None if the run is rejected.

    A run that did not reach CHROMO_T_END is not a 4000 s result and must not sit
    on the same axes as one that did (the 6 kK wall aborts at t~380 s below the
    EOS table floor). Reject it here rather than truncating the movie to it.
    """
    term, fallbacks = console_ok(path)
    if term != "end_time":
        print(f"{t_wall/1e3:>5g} kK  EXCLUDED: termination={term}")
        return None
    h, times, data = load_gamma_diag(path)
    q = np.stack([physical_heat_flux(h, data[k, :, IT], data[k, :, I_NE],
                                    data[k, :, I_NHI])
                  for k in range(data.shape[0])])
    # The conservative face flux lives on its own height grid (face heights) and
    # its own time grid (a step-based capture stride, not the snapshot cadence),
    # so it is kept separate rather than merged into the cell-centred arrays.
    face = None
    if os.path.exists(path + ".faceflux"):
        face_h, face_t, face_f = load_face_flux(path + ".faceflux")
        face = dict(h=np.asarray(face_h, dtype=float),
                    times=np.asarray(face_t, dtype=float),
                    f=[np.asarray(x, dtype=float) for x in face_f])
    print(f"{t_wall/1e3:>5g} kK  term={term} fallbacks={fallbacks} "
          f"frames={len(times)} t_end={times[-1]:.1f} s ns={h.size} "
          f"faceflux={'-' if face is None else str(len(face['times'])) + ' captures'}")
    return dict(t_wall=t_wall, h=h, times=times, face=face,
                T=data[:, :, IT], V=data[:, :, IV], rho=data[:, :, IRHO],
                p=data[:, :, IP], q=q,
                f=data[:, :, IRHO] * data[:, :, IV])


def symlim(arrays, linthresh_frac=1.0e-4, pad=1.4):
    """Symmetric symlog limits + linear threshold from the pooled data range."""
    a = np.concatenate([np.asarray(x).ravel() for x in arrays])
    a = a[np.isfinite(a)]
    lo, hi = float(a.min()), float(a.max())
    span = max(abs(lo), abs(hi))
    return (lo * pad if lo < 0 else lo / pad,
            hi * pad if hi > 0 else hi / pad,
            max(span * linthresh_frac, np.finfo(float).tiny))


def loglim(arrays, pad=1.3):
    a = np.concatenate([np.asarray(x).ravel() for x in arrays])
    a = a[np.isfinite(a) & (np.asarray(a) > 0.0)]
    return float(a.min()) / pad, float(a.max()) * pad


def linlim(arrays, pad=0.04):
    """Linear limits padded by a fraction of the pooled span (for the p ratio)."""
    a = np.concatenate([np.asarray(x).ravel() for x in arrays])
    a = a[np.isfinite(a)]
    lo, hi = float(a.min()), float(a.max())
    m = pad * max(hi - lo, np.finfo(float).tiny)
    return lo - m, hi + m


def _render_frame(k, tmpdir, h, mflux_h, series, colors, widths, t_walls, t, nF,
                  ylims, title, mflux_label, p_label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    fig, ax = plt.subplots(2, 3, figsize=(17, 8))
    fig.suptitle(f"{title}  |  t = {t:7.1f} s   ({k + 1}/{nF})", fontsize=13)

    specs = [
        (ax[0, 0], "T",    h,        "T [K]",              "Temperature"),
        (ax[0, 1], "V",    h,        "V [m/s]",
         "Velocity  (>0 up = evaporation, <0 down = condensation)"),
        (ax[0, 2], "q",    h,        r"$q$ [W/m$^2$]",
         r"$q_\parallel=-(\kappa_e+\kappa_n)dT/ds$"),
        (ax[1, 0], "rho",  h,        r"$\rho$ [kg/m$^3$]", "Mass density"),
        (ax[1, 1], "p",    h,        p_label[0],           p_label[1]),
        # The mass flux is on the FACE grid when it is f_total, so it carries its
        # own x array rather than the cell centres of the other five panels.
        (ax[1, 2], "f",    mflux_h,  mflux_label[0],       mflux_label[1]),
    ]
    for a, key, x, ylabel, panel_title in specs:
        for y, c, lw in zip(series[key], colors, widths):
            a.plot(x, y, color=c, lw=lw)
        lim = ylims[key]
        if lim[0] == "log":
            a.set_yscale("log")
            a.set_ylim(lim[1], lim[2])
        elif lim[0] == "lin":
            # The p ratio: unity is the release baseline, so mark it.
            a.axhline(1.0, c="k", lw=0.5)
            a.set_ylim(lim[1], lim[2])
        else:
            a.axhline(0.0, c="k", lw=0.5)
            a.set_yscale("symlog", linthresh=lim[3])
            a.set_ylim(lim[1], lim[2])
        a.set_ylabel(ylabel)
        a.set_title(panel_title, fontsize=10)
        a.set_xlim(x.min(), x.max())
        a.set_xlabel("height [km]")
        a.grid(True, alpha=0.3, which="both")

    fig.tight_layout(rect=[0, 0.10, 1, 0.955])

    # DISCRETE colorbar, one equal-width band per run, built from the ACTUAL line
    # colours. It cannot be a continuous LogNorm bar: wall_colors() spaces the
    # colours by RANK and caps their luminance (see util/twall_sweep_diag.py), so
    # a bar drawn from turbo + LogNorm would show colours that no line has and
    # place them at the wrong T_wall. Equal-width bands are the honest rendering
    # of a rank-spaced assignment — the bar is a legend with an axis, and its
    # width is NOT proportional to T_wall.
    cmap_d = ListedColormap(colors)
    norm_d = BoundaryNorm(np.arange(len(colors) + 1) - 0.5, cmap_d.N)
    cax = fig.add_axes([0.20, 0.055, 0.60, 0.020])
    cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap_d, norm=norm_d),
                        cax=cax, orientation="horizontal",
                        ticks=np.arange(len(colors)))
    cbar.ax.set_xticklabels([f"{tw / 1e3:g}" for tw in t_walls], fontsize=9)
    cbar.outline.set_linewidth(0.6)
    cbar.set_label("imposed conduction-wall temperature $T_{\\rm wall}$ [kK]  "
                   "(equal-width bands, one per run; not linear in "
                   "$T_{\\rm wall}$)", fontsize=9.5)
    for lbl, tw in zip(cbar.ax.get_xticklabels(), t_walls):
        if tw == T_WALL_RELEASE:                 # the release baseline
            lbl.set_fontweight("bold")
            lbl.set_fontsize(10.5)
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="outputs/model_column/twall_sweep")
    ap.add_argument("--cfl", default="0.50",
                    help="only animate runs with this CFL tag; 'all' takes every run")
    ap.add_argument("--run", action="append", default=[],
                    help="explicit T_wall_K=path entry, added to the glob")
    ap.add_argument("--out",
                    default="visualization/model_column/"
                            "twall_sweep_N500_4000s_overlay.mp4")
    ap.add_argument("--fps", type=int, default=20)
    args = ap.parse_args()

    entries = collect_entries(args.dir, args.cfl, args.run)
    if not entries:
        raise SystemExit(f"no runs matched in {args.dir} (cfl={args.cfl})")
    runs = [r for r in (load_run(tw, p) for tw, p in entries) if r is not None]
    if not runs:
        raise SystemExit("every matched run was rejected")

    # All runs must share the mesh; they are the same one-parameter configuration.
    h = runs[0]["h"]
    for r in runs[1:]:
        if r["h"].shape != h.shape or not np.allclose(r["h"], h):
            raise SystemExit(f"mesh mismatch at T_wall={r['t_wall']:g} K")

    # Common time grid = the shortest run's snapshot times; every other run is
    # sampled at its nearest snapshot. The worst mismatch is reported rather than
    # hidden, exactly as the face-flux picker in util/animate_isentropic.py does.
    base = min(runs, key=lambda r: len(r["times"]))
    times = base["times"]
    max_frames = int(os.environ.get("ANIM_MAX_FRAMES", "400"))
    if len(times) > max_frames:
        times = times[np.unique(np.linspace(0, len(times) - 1,
                                            max_frames).round().astype(int))]
    picks, worst = [], 0.0
    for r in runs:
        idx = np.abs(r["times"][None, :] - times[:, None]).argmin(axis=1)
        picks.append(idx)
        worst = max(worst, float(np.abs(r["times"][idx] - times).max()))
    nF = len(times)
    print(f"{len(runs)} runs, {nF} frames, t in [{times[0]:.1f}, {times[-1]:.1f}] s; "
          f"worst run-to-run snapshot time mismatch {worst:.3f} s")

    # ------------------------------------------------------------------
    # Mass-flux panel: the conservative face flux f_total when EVERY run has a
    # .faceflux sidecar, else the cell-centred product rho*V. All-or-nothing on
    # purpose: drawing f_total for some walls and rho*V for others on one axis
    # would be comparing two different quantities across the sweep.
    # ------------------------------------------------------------------
    requested = os.environ.get("ANIM_MASS_FLUX")
    mode = (requested or "face").lower()
    if mode not in ("face", "cell"):
        raise ValueError("ANIM_MASS_FLUX must be 'face' or 'cell'")
    missing = [r["t_wall"] for r in runs if r["face"] is None]
    if mode == "face" and missing:
        walls = ", ".join(f"{t/1e3:g} kK" for t in missing)
        msg = (f"no .faceflux sidecar for {walls}; regenerate the sweep with "
               f"SWEEP_FACE_FLUX=1 util/run_twall_sweep.sh")
        if requested:                       # explicitly asked for the face flux
            raise SystemExit(msg)
        print(f"note: {msg}. Falling back to the cell-centred rho*V panel.")
    use_face = mode == "face" and not missing

    if use_face:
        mflux_h = runs[0]["face"]["h"]
        for r in runs[1:]:
            if not np.allclose(r["face"]["h"], mflux_h):
                raise SystemExit(f"face grid mismatch at T_wall={r['t_wall']:g} K")
        picks_f, worst_f = [], 0.0
        for r in runs:
            ft = r["face"]["times"]
            idx = np.abs(ft[None, :] - times[:, None]).argmin(axis=1)
            picks_f.append(idx)
            worst_f = max(worst_f, float(np.abs(ft[idx] - times).max()))
        # The capture cadence is in STEPS and the snapshot cadence in SECONDS, so
        # the two grids do not coincide and the mass-flux line in a given frame can
        # be up to this far from the other five panels. Reported, not hidden.
        # Capture counts differ between walls: the cadence is a fixed number of
        # STEPS and a hotter wall runs a shorter dt, so it captures more often.
        n_cap = [len(r["face"]["times"]) for r in runs]
        print(f"face flux: {min(n_cap)}-{max(n_cap)} captures/run; worst "
              f"frame-to-capture time mismatch {worst_f:.2f} s")
        pooled_f = [np.stack(r["face"]["f"]) for r in runs]
        mflux_label = (r"$f_{\rm total}$ [kg/m$^2$/s]",
                       "Conservative numerical FACE mass flux  (>0 up, <0 down)")
    else:
        mflux_h = h
        picks_f = picks
        pooled_f = [r["f"] for r in runs]
        mflux_label = (r"$\rho V$ [kg/m$^2$/s]",
                       "Cell-centred rho*V (NOT the face flux)  (>0 up, <0 down)")

    # ------------------------------------------------------------------
    # Pressure panel: p RELATIVE to the 22 kK release run at the same frame.
    # Absolute p separates at no y scale — it falls 6.6x over the column while the
    # wall-to-wall spread is at most 12 %, and it falls that 12 % in ~30 km, so a
    # y zoom crops x with it and only compresses the fan further. Dividing by the
    # release run removes the shared stratification and the shared early acoustic
    # transient and leaves a fan that is separated over the whole domain.
    # ------------------------------------------------------------------
    p_mode = os.environ.get("ANIM_P_PANEL", "rel").lower()
    if p_mode not in ("rel", "abs"):
        raise ValueError("ANIM_P_PANEL must be 'rel' or 'abs'")
    ref = next((j for j, r in enumerate(runs) if r["t_wall"] == T_WALL_RELEASE), None)
    if p_mode == "rel" and ref is None:
        print(f"note: no {T_WALL_RELEASE/1e3:g} kK run in this set; falling back to "
              f"the absolute p panel.")
        p_mode = "abs"

    if p_mode == "rel":
        p_ref = runs[ref]["p"][picks[ref]]                     # (nF, ns)
        pooled_p = [runs[j]["p"][picks[j]] / p_ref for j in range(len(runs))]
        p_index = None                                         # already on frames
        p_label = (r"$p\,/\,p_{22\,\rm kK}$",
                   "Pressure relative to the 22 kK release baseline (dimensionless)")
        p_lim = ("lin",) + linlim(pooled_p)
    else:
        pooled_p = [r["p"] for r in runs]
        p_index = picks
        p_label = ("p [Pa]", "Pressure")
        p_lim = ("log",) + loglim(pooled_p)

    keys = ("T", "V", "q", "rho")
    pooled = {k: [r[k] for r in runs] for k in keys}
    pooled["f"] = pooled_f
    ylims = {
        "T": ("log",) + loglim(pooled["T"]),
        "rho": ("log",) + loglim(pooled["rho"]),
        "p": p_lim,
        "V": ("symlog",) + symlim(pooled["V"], linthresh_frac=1.0e-4),
        # q needs a wider linear band than V and rho*V: its range is set by the
        # TR/boundary spike (~100 W/m^2) while the chromospheric interior sits at
        # ~1e-3, so a 1e-4 threshold puts five decades of physically irrelevant
        # interior detail around zero and the +/- tick labels collide.
        "q": ("symlog",) + symlim(pooled["q"], linthresh_frac=1.0e-3),
        "f": ("symlog",) + symlim(pooled["f"], linthresh_frac=1.0e-4),
    }

    t_walls = [r["t_wall"] for r in runs]
    col, _, _ = wall_colors(t_walls)
    colors = [col[tw] for tw in t_walls]
    widths = [2.4 if tw == T_WALL_RELEASE else 1.2 for tw in t_walls]
    title =("model_column upper-boundary sweep — N=500/R4, 4000 s, double state — "
             "one line per imposed $T_{\\rm wall}$ (22 kK release baseline drawn "
             "thick)")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    frame_args = []
    for k in range(nF):
        series = {key: [runs[j][key][picks[j][k]] for j in range(len(runs))]
                  for key in keys}
        series["f"] = [pooled_f[j][picks_f[j][k]] for j in range(len(runs))]
        series["p"] = [pooled_p[j][k if p_index is None else p_index[j][k]]
                       for j in range(len(runs))]
        frame_args.append((k, h, mflux_h, series, colors, widths, t_walls,
                           float(times[k]), nF, ylims, title, mflux_label,
                           p_label))
    save_frames_parallel(_render_frame, frame_args, args.out, args.fps)


if __name__ == "__main__":
    main()
