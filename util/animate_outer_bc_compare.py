#!/usr/bin/env python3
"""Overlay two (or more) model_column upper-boundary CLOSURES on one 6-panel movie.

The release outer face imposes the ONE condition a subsonic outflow admits as a
fixed coronal back-pressure (``kOuterPRef``) and lets the hydro ghost temperature
free-float with the live top cell (``ISO_HYDRO_T_DECOUPLE=1``), with the 22 kK
conduction wall a separate physical-face datum. The experiment animated here
SWAPS which variable carries that condition: the hydro ghost temperature is pinned
to the wall (``ISO_OUTER_T_HYDRO_WALL=1``, so T_hydro = T_cond = 22 kK) and the
pressure becomes zero-gradient (``ISO_OUTER_P_NEUMANN=1``).

Panels and quantities follow ``util/animate_twall_sweep.py`` --

    T(h)      V(h)      q_parallel(h)
    rho(h)    p(h)      mass flux(h)

-- so this movie reads against both the sweep overlay and the single-run release
movie panel for panel. Differences from the sweep overlay, all forced by comparing
two CLOSURES rather than one parameter:

* Runs are keyed by a LABEL with a legend, not by T_wall with a colorbar: both runs
  here have the same 22 kK wall, so T_wall no longer identifies a line.
* p is ABSOLUTE (log Pa), not relative to the release run. The sweep needed the
  ratio because a 12 % wall-to-wall spread is invisible against a 6.6x
  stratification; a closure swap moves p by factors of several, which absolute log
  Pa shows directly and a ratio panel would only rescale.
* A run that did NOT reach CHROMO_T_END is kept, not excluded: for a boundary-
  closure experiment an early abort is the result. The movie then covers the
  overlap window and the legend reports each run's termination.

Usage:
    python util/animate_outer_bc_compare.py [--run "LABEL=path.txt"]...
        [--out FILE.mp4] [--fps 20]

Defaults reproduce the catalogued figure. Env:
    ANIM_MAX_FRAMES  frame budget (default 400; the run is subsampled in time)
    ANIM_MASS_FLUX   face (default, requires .faceflux sidecars) | cell (rho*V)
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _anim_parallel import save_frames_parallel
from animate_isentropic import load_face_flux
from animate_twall_sweep import loglim, symlim
from conductive_flux import physical_heat_flux
from twall_sweep_diag import COLUMNS, IRHO, IT, IP, IV, console_ok, load_gamma_diag

I_NE = COLUMNS.index("n_e")
I_NHI = COLUMNS.index("n_HI")

# Two catalogued run sets. `closure` varies the outer closure at one wall
# temperature; `twall` varies BOTH, so there colour carries the wall temperature and
# the dash pattern carries the closure — one visual channel per varied quantity,
# which a colour-only scheme cannot do for a 2x2.
#
# In `closure`, order matters: colours are assigned by position, so the
# release/mode-1/mode-2 prefix is kept first and the control last, and a three-run
# subset (the transient movie, which has no control run) draws in the same colours
# as the full set. The release closure is drawn thicker, in the neutral dark grey
# the other release figures use; the experiment is the accent red.
PRESETS = {
    "closure": dict(
        title=("model_column outer-face closure — N=500/R4, double state, CFL 0.50, "
               "22 kK conduction wall: WHICH variable carries the one imposed "
               "condition"),
        runs=[
            ("release: fixed $p_{\\rm ref}$, free $T_{\\rm hydro}$",
             "outputs/model_column/twall_sweep/twall_22kK_cfl0.50.txt"),
            ("swap: $T_{\\rm hydro}$ pinned, $dp/ds=0$",
             "outputs/model_column/outer_bc_swap/tfix_pneumann_4000s.txt"),
            ("swap: $T_{\\rm hydro}$ pinned, $dp/ds=-\\rho g$",
             "outputs/model_column/outer_bc_swap/tfix_phse_4000s.txt"),
            ("control: $T_{\\rm hydro}$ pinned, fixed $p_{\\rm ref}$",
             "outputs/model_column/outer_bc_swap/tfix_pref_4000s.txt"),
        ],
        colors=["#2f3640", "#c0392b", "#1f6fb2", "#7f8c8d"],
        widths=[2.8, 1.8, 1.8, 1.4],
        styles=["-", "-", "-", "-"],
        f_linthresh_frac=1.0e-3,
    ),
    "twall": dict(
        title=("model_column wall temperature x outer pressure closure — N=500/R4, "
               "double state  (colour = $T_{\\rm wall}$, dashed = unanchored $p$)"),
        runs=[
            ("release BC, $T_{\\rm wall}=22$ kK",
             "outputs/model_column/twall_sweep/twall_22kK_cfl0.50.txt"),
            ("release BC, $T_{\\rm wall}=100$ kK",
             "outputs/model_column/twall_sweep/twall_100kK_cfl0.25.txt"),
            ("$dp/ds=-\\rho g$, $T_{\\rm wall}=22$ kK",
             "outputs/model_column/outer_bc_swap/hse_22kK_cfl0.25.txt"),
            ("$dp/ds=-\\rho g$, $T_{\\rm wall}=100$ kK",
             "outputs/model_column/outer_bc_swap/hse_100kK_cfl0.25.txt"),
        ],
        colors=["#1f6fb2", "#c0392b", "#1f6fb2", "#c0392b"],
        widths=[2.4, 2.4, 1.8, 1.8],
        styles=["-", "-", "--", "--"],
        # Wider symlog linear band than the closure preset: the 100 kK release run
        # reaches ~5e-9 while the 22 kK unanchored run works at ~1e-11, and a 1e-3
        # band puts eight labelled decades on one short axis.
        f_linthresh_frac=1.0e-2,
    ),
}


def load_run(label, path):
    """-> dict of the six panel fields for one run (early termination is kept)."""
    term, fallbacks = console_ok(path)
    h, times, data = load_gamma_diag(path)
    q = np.stack([physical_heat_flux(h, data[k, :, IT], data[k, :, I_NE],
                                     data[k, :, I_NHI])
                  for k in range(data.shape[0])])
    face = None
    if os.path.exists(path + ".faceflux"):
        face_h, face_t, face_f = load_face_flux(path + ".faceflux")
        face = dict(h=np.asarray(face_h, dtype=float),
                    times=np.asarray(face_t, dtype=float),
                    f=[np.asarray(x, dtype=float) for x in face_f])
    print(f"{label}\n    term={term} fallbacks={fallbacks} frames={len(times)} "
          f"t_end={times[-1]:.1f} s ns={h.size} "
          f"faceflux={'-' if face is None else str(len(face['times'])) + ' captures'}")
    return dict(label=label, term=term, h=h, times=times, face=face,
                T=data[:, :, IT], V=data[:, :, IV], rho=data[:, :, IRHO],
                p=data[:, :, IP], q=q,
                f=data[:, :, IRHO] * data[:, :, IV])


def _render_frame(k, tmpdir, h, mflux_h, series, labels, style, t, nF, ylims, title,
                  mflux_label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 3, figsize=(17, 8))
    fig.suptitle(f"{title}  |  t = {t:7.1f} s   ({k + 1}/{nF})", fontsize=13)

    specs = [
        (ax[0, 0], "T",   h,       "T [K]",              "Temperature"),
        (ax[0, 1], "V",   h,       "V [m/s]",
         "Velocity  (>0 up = evaporation, <0 down = condensation)"),
        (ax[0, 2], "q",   h,       r"$q$ [W/m$^2$]",
         r"$q_\parallel=-(\kappa_e+\kappa_n)dT/ds$"),
        (ax[1, 0], "rho", h,       r"$\rho$ [kg/m$^3$]", "Mass density"),
        (ax[1, 1], "p",   h,       "p [Pa]",             "Total gas pressure"),
        (ax[1, 2], "f",   mflux_h, mflux_label[0],       mflux_label[1]),
    ]
    handles = None
    for a, key, x, ylabel, panel_title in specs:
        lines = [a.plot(x, y, color=c, lw=lw, ls=ls, label=lab)[0]
                 for y, c, lw, ls, lab in zip(series[key], style["colors"],
                                              style["widths"], style["styles"],
                                              labels)]
        if handles is None:
            handles = lines
        lim = ylims[key]
        if lim[0] == "log":
            a.set_yscale("log")
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

    # Two legend columns from three runs on: four of these labels in one row are
    # wider than the figure and get clipped at both ends.
    ncol = 2 if len(labels) > 2 else len(labels)
    fig.tight_layout(rect=[0, 0.075 if ncol == 2 else 0.05, 1, 0.955])
    fig.legend(handles, labels, loc="lower center", ncol=ncol,
               frameon=False, fontsize=10.5)
    fig.savefig(os.path.join(tmpdir, f"frame_{k:06d}.png"), dpi=120)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", choices=sorted(PRESETS), default="closure",
                    help="catalogued run set: 'closure' varies the outer closure at "
                         "22 kK, 'twall' the 2x2 of closure x wall temperature")
    ap.add_argument("--run", action="append", default=[],
                    help="LABEL=path entry; repeat. Overrides the preset's runs but "
                         "keeps its colours/styles. Split on the LAST '=' so a TeX "
                         "label may contain one.")
    ap.add_argument("--out", default="visualization/model_column/"
                                     "outer_bc_tfix_pneumann_vs_release.mp4")
    ap.add_argument("--fps", type=int, default=20)
    args = ap.parse_args()

    preset = PRESETS[args.preset]
    entries = ([tuple(spec.rsplit("=", 1)) for spec in args.run]
               if args.run else preset["runs"])
    if len(entries) > len(preset["colors"]):
        raise SystemExit(f"at most {len(preset['colors'])} runs in preset "
                         f"{args.preset!r}")
    runs = [load_run(lab, path) for lab, path in entries]

    # Same configuration except the outer closure, so the mesh must match.
    h = runs[0]["h"]
    for r in runs[1:]:
        if r["h"].shape != h.shape or not np.allclose(r["h"], h):
            raise SystemExit(f"mesh mismatch for {r['label']!r}")

    # Common time grid = the OVERLAP window (the shortest run's snapshot times);
    # every other run is sampled at its nearest snapshot, worst mismatch reported.
    base = min(runs, key=lambda r: r["times"][-1])
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
    print(f"{len(runs)} runs, {nF} frames, t in [{times[0]:.1f}, {times[-1]:.1f}] s "
          f"(overlap set by {base['label']!r}); worst snapshot mismatch {worst:.3f} s")

    # Mass-flux panel: the conservative face flux f_total when EVERY run has a
    # sidecar, else the cell-centred product rho*V for all of them. All-or-nothing:
    # the two are different quantities and must not share an axis.
    requested = os.environ.get("ANIM_MASS_FLUX")
    mode = (requested or "face").lower()
    if mode not in ("face", "cell"):
        raise ValueError("ANIM_MASS_FLUX must be 'face' or 'cell'")
    missing = [r["label"] for r in runs if r["face"] is None]
    if mode == "face" and missing:
        msg = (f"no .faceflux sidecar for {', '.join(missing)}; rerun with "
               f"CHROMO_FACE_FLUX_DIAG=1")
        if requested:
            raise SystemExit(msg)
        print(f"note: {msg}. Falling back to the cell-centred rho*V panel.")
    use_face = mode == "face" and not missing

    if use_face:
        mflux_h = runs[0]["face"]["h"]
        for r in runs[1:]:
            if not np.allclose(r["face"]["h"], mflux_h):
                raise SystemExit(f"face grid mismatch for {r['label']!r}")
        picks_f, worst_f = [], 0.0
        for r in runs:
            ft = r["face"]["times"]
            idx = np.abs(ft[None, :] - times[:, None]).argmin(axis=1)
            picks_f.append(idx)
            worst_f = max(worst_f, float(np.abs(ft[idx] - times).max()))
        # The capture cadence is in STEPS and the snapshot cadence in SECONDS, so
        # the mass-flux line can be this far from the other five panels in a frame.
        print(f"face flux: {[len(r['face']['times']) for r in runs]} captures/run; "
              f"worst frame-to-capture mismatch {worst_f:.2f} s")
        pooled_f = [np.stack(r["face"]["f"]) for r in runs]
        mflux_label = (r"$f_{\rm total}$ [kg/m$^2$/s]",
                       "Conservative numerical FACE mass flux  (>0 up, <0 down)")
    else:
        mflux_h = h
        picks_f = picks
        pooled_f = [r["f"] for r in runs]
        mflux_label = (r"$\rho V$ [kg/m$^2$/s]",
                       "Cell-centred rho*V (NOT the face flux)  (>0 up, <0 down)")

    keys = ("T", "V", "q", "rho", "p")
    pooled = {k: [r[k][picks[j]] for j, r in enumerate(runs)] for k in keys}
    pooled["f"] = [pooled_f[j][picks_f[j]] for j in range(len(runs))]
    ylims = {
        "T": ("log",) + loglim(pooled["T"]),
        "rho": ("log",) + loglim(pooled["rho"]),
        "p": ("log",) + loglim(pooled["p"]),
        "V": ("symlog",) + symlim(pooled["V"], linthresh_frac=1.0e-4),
        "q": ("symlog",) + symlim(pooled["q"], linthresh_frac=1.0e-3),
        # Wider linear band than the sweep overlay uses: with a collapse run on the
        # same axis the range is set by ~1e-7 while the release sits at ~1e-10, and a
        # 1e-4 threshold puts the +/- decade labels on top of each other near zero.
        # The exact value is per preset — see `f_linthresh_frac`.
        "f": ("symlog",) + symlim(pooled["f"],
                                  linthresh_frac=preset["f_linthresh_frac"]),
    }

    # The termination tag is shown only when a run did NOT reach its end time: it is
    # then part of the result, and four "[end_time, t_end=4000 s]" suffixes are just
    # legend width.
    labels = [r["label"] if r["term"] == "end_time" else
              f"{r['label']}  [{r['term']}, $t_{{\\rm end}}$={r['times'][-1]:.0f} s]"
              for r in runs]
    style = {k: preset[k] for k in ("colors", "widths", "styles")}

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    frame_args = []
    for k in range(nF):
        series = {key: [pooled[key][j][k] for j in range(len(runs))]
                  for key in keys}
        series["f"] = [pooled["f"][j][k] for j in range(len(runs))]
        frame_args.append((k, h, mflux_h, series, labels, style, float(times[k]),
                           nF, ylims, preset["title"], mflux_label))
    save_frames_parallel(_render_frame, frame_args, args.out, args.fps)


if __name__ == "__main__":
    main()
