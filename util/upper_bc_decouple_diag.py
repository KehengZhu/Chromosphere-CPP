#!/usr/bin/env python3
"""Quantitative upper-boundary diagnostics for the model_column hydro-T /
conduction-T decoupling experiment (ISO_HYDRO_T_DECOUPLE).

Reads the EOS-aware ``<run>.gamma_diag`` sidecar (columns
``rho_total v_cm T x_eq n_e n_HI p_total Gamma1 kappa_physical kappa_solver``)
and reports, per run and per requested time:

  V_last, V_{N-2}                          [m/s]
  (rho V)_last / (rho V)_mid               [-]     mid = cell ns//2
  rel_std(rho V) over 2130-2150 km         [-]     std / |mean|
  roughness over the top N cells           [-]     mean|d2(rho V)| / mean|rho V|
  p_top - p_g0                             [Pa]    p_g0 = --p-back (constant)
  T_hydro_g0 - T_top                       [K]     wall - T_top, or 0 if decoupled
  q_top(t)                                 [W/m^2] kappa_solver,top (T_wall-T_top)/ds
  E_cond(t) = int q_top dt                 [J/m^2]
  M(t) = sum rho ds  and  M(t)/M(0) - 1    [kg/m^2]
  boundary-layer thickness, max ripple     [km], [-]

The boundary layer is the contiguous run of top cells whose rho V departs from
the interior plateau (median of rho V over the top 400..90 cells) by more than
``--bl-tol`` of that plateau; the ripple amplitude is the largest such relative
departure over the top ``--top-cells`` cells.

q_top uses the CELL conductivity of the top cell rather than the outer FACE
average (the wall-side kappa needs the ghost density, which the sidecar does not
carry), so it is a first-order estimate of the discrete Dirichlet flux
kappa_face (T_wall - T_top)/ds. Both runs share the wall and very nearly share
kappa, so the comparison between runs is apples-to-apples.

Usage:
    python util/upper_bc_decouple_diag.py --run label=path.txt [--run ...] \
        [--times 400 500] [--t-wall 22000] [--p-back label=Pa] [--json out.json]
"""
import argparse
import json
import os
import sys

import numpy as np

COLUMNS = ["rho", "v", "T", "x_eq", "n_e", "n_HI", "p", "Gamma1",
           "kappa_physical", "kappa_solver"]


def load_gamma_diag(path):
    """-> (heights_km, [(t, array[ns, 10]), ...]), stopping at the first bad frame."""
    diag = path if path.endswith(".gamma_diag") else path + ".gamma_diag"
    with open(diag) as handle:
        ns, ncol = map(int, handle.readline().split())
        h = np.array([float(x) for x in handle.readline().split()])
        frames, rows, t = [], None, None

        def flush():
            if rows is None or len(rows) != ns:
                return False
            arr = np.array(rows)
            if not np.isfinite(arr).all():
                return False
            frames.append((t, arr))
            return True

        for line in handle:
            if line.startswith("# t ="):
                if rows is not None and not flush():
                    break
                t = float(line.split("t =")[1].split("step")[0])
                rows = []
            elif line.startswith("#"):
                continue
            elif rows is not None:
                parts = line.split()
                if len(parts) == ncol:
                    rows.append([float(x) for x in parts])
        else:
            flush()
    if not frames:
        raise SystemExit(f"no finite frames in {diag}")
    return h, frames


def pick(frames, t_target):
    times = np.array([f[0] for f in frames])
    return frames[int(np.argmin(np.abs(times - t_target)))]


def roughness(mass_flux, n_top):
    """mean |second difference| / mean |value| over the top n_top cells."""
    seg = mass_flux[-n_top:]
    d2 = seg[2:] - 2.0 * seg[1:-1] + seg[:-2]
    return float(np.mean(np.abs(d2)) / np.mean(np.abs(seg)))


def boundary_layer(mass_flux, n_top, tol):
    """(thickness in cells, max relative ripple) against the interior plateau."""
    plateau = np.median(mass_flux[-min(400, len(mass_flux)):-n_top])
    seg = mass_flux[-n_top:]
    rel = (seg - plateau) / abs(plateau)
    thickness = 0
    for k in range(len(rel) - 1, -1, -1):          # walk down from the top
        if abs(rel[k]) > tol:
            thickness = len(rel) - k
        else:
            break
    return thickness, float(np.max(np.abs(rel)))


def analyse(label, path, times, t_wall, p_back, n_top, bl_tol, decoupled):
    h, frames = load_gamma_diag(path)
    ds_km = float(np.mean(np.diff(h)))
    ds = ds_km * 1.0e3
    ns = len(h)
    band = (h >= 2130.0) & (h <= 2150.0)

    # Cumulative conductive boundary input over ALL frames (trapezoid in time).
    q_series, t_series = [], []
    for t, arr in frames:
        kappa = arr[-1, COLUMNS.index("kappa_solver")]
        T_top = arr[-1, COLUMNS.index("T")]
        q_series.append(kappa * (t_wall - T_top) / ds)
        t_series.append(t)
    q_series, t_series = np.array(q_series), np.array(t_series)
    e_cond = np.concatenate([[0.0], np.cumsum(
        0.5 * (q_series[1:] + q_series[:-1]) * np.diff(t_series))])
    mass = np.array([arr[:, COLUMNS.index("rho")].sum() * ds for _, arr in frames])

    out = {"label": label, "path": path, "ns": ns, "ds_km": ds_km,
           "t_end": float(t_series[-1]), "n_frames": len(frames),
           "t_wall_K": t_wall, "p_back_Pa": p_back, "times": {}}
    for t_target in times:
        t, arr = pick(frames, t_target)
        rho = arr[:, COLUMNS.index("rho")]
        v = arr[:, COLUMNS.index("v")]
        T = arr[:, COLUMNS.index("T")]
        p = arr[:, COLUMNS.index("p")]
        mf = rho * v
        idx = int(np.argmin(np.abs(t_series - t)))
        thickness_cells, ripple = boundary_layer(mf, n_top, bl_tol)
        out["times"][f"{t:.1f}"] = {
            "t": float(t),
            "V_last": float(v[-1]),
            "V_second_last": float(v[-2]),
            "massflux_last_over_mid": float(mf[-1] / mf[ns // 2]),
            "massflux_mid": float(mf[ns // 2]),
            "rel_std_massflux_2130_2150": float(np.std(mf[band])
                                                / abs(np.mean(mf[band]))),
            f"roughness_top{n_top}": roughness(mf, n_top),
            "T_top": float(T[-1]),
            # Decoupled mode extrapolates T_hydro_g0 = T_top, so the jump is 0 by
            # construction; baseline pins the ghost at the wall.
            "T_hydro_g0_minus_T_top": (0.0 if decoupled
                                       else float(t_wall - T[-1])),
            "p_top": float(p[-1]),
            "p_top_minus_p_back": (None if p_back is None
                                   else float(p[-1] - p_back)),
            "q_top_Wm2": float(q_series[idx]),
            "E_cond_Jm2": float(e_cond[idx]),
            "mass_kgm2": float(mass[idx]),
            "mass_rel_change": float(mass[idx] / mass[0] - 1.0),
            "boundary_layer_km": float(thickness_cells * ds_km),
            "boundary_layer_cells": int(thickness_cells),
            "max_ripple_rel": ripple,
            "max_abs_V": float(np.max(np.abs(v))),
        }
    out["series"] = {"t": t_series.tolist(), "q_top_Wm2": q_series.tolist(),
                     "E_cond_Jm2": e_cond.tolist(), "mass_kgm2": mass.tolist()}
    return out


def plot_top_region(run_specs, out_path, t_target, t_wall, n_top):
    """Side-by-side mass flux / velocity / temperature over the top n_top cells."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    stamp = None
    for spec in run_specs:
        label, path = spec.split("=", 1)
        h, frames = load_gamma_diag(path)
        t, arr = pick(frames, t_target)
        stamp = t
        hs = h[-n_top:]
        rho, v, T = arr[-n_top:, 0], arr[-n_top:, 1], arr[-n_top:, 2]
        ax[0].plot(hs, rho * v, lw=1.3, marker=".", ms=3, label=label)
        ax[1].plot(hs, v, lw=1.3, marker=".", ms=3, label=label)
        ax[2].plot(hs, T, lw=1.3, marker=".", ms=3, label=label)
    ax[0].axhline(0.0, c="k", lw=0.5)
    ax[0].set_ylabel(r"$\rho V$  [kg m$^{-2}$ s$^{-1}$]")
    ax[1].axhline(0.0, c="k", lw=0.5)
    ax[1].set_ylabel("V  [m/s]")
    ax[2].axhline(t_wall, c="k", ls=":", lw=0.9,
                  label=f"conduction wall {t_wall:g} K")
    ax[2].set_ylabel("T  [K]")
    for a in ax:
        a.set_xlabel("h  [km]")
        a.legend(fontsize=8)
        a.grid(alpha=0.25)
    fig.suptitle(f"Upper-boundary top {n_top} cells   |   t = {stamp:.1f} s")
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=140)
    print(f"wrote {out_path}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True,
                    metavar="LABEL=PATH", help="run to analyse (repeatable)")
    ap.add_argument("--times", nargs="*", type=float, default=[500.0])
    ap.add_argument("--t-wall", type=float, default=22000.0,
                    help="conduction Dirichlet wall temperature [K]")
    ap.add_argument("--p-back", action="append", default=[], metavar="LABEL=PA",
                    help="captured outer back-pressure kOuterPRef for a run")
    ap.add_argument("--decoupled", action="append", default=[], metavar="LABEL",
                    help="label(s) whose hydro ghost T follows the live top cell")
    ap.add_argument("--top-cells", type=int, default=90)
    ap.add_argument("--bl-tol", type=float, default=0.05)
    ap.add_argument("--json", default=None)
    ap.add_argument("--plot", default=None,
                    help="write a top-region rho V / V / T comparison figure here")
    ap.add_argument("--plot-time", type=float, default=None,
                    help="time for --plot (default: last of --times)")
    args = ap.parse_args(argv)

    p_back = {}
    for item in args.p_back:
        label, value = item.split("=", 1)
        p_back[label] = float(value)

    results = []
    for item in args.run:
        label, path = item.split("=", 1)
        results.append(analyse(label, path, args.times, args.t_wall,
                               p_back.get(label), args.top_cells, args.bl_tol,
                               label in args.decoupled))

    keys = ["V_last", "V_second_last", "massflux_last_over_mid",
            "rel_std_massflux_2130_2150", f"roughness_top{args.top_cells}",
            "T_top", "T_hydro_g0_minus_T_top", "p_top", "p_top_minus_p_back", "q_top_Wm2", "E_cond_Jm2",
            "mass_kgm2", "mass_rel_change", "boundary_layer_km",
            "max_ripple_rel", "max_abs_V"]
    for t_target in args.times:
        stamp = f"{pick(load_gamma_diag(results[0]['path'])[1], t_target)[0]:.1f}"
        print(f"\n=== t ~ {t_target:g} s (nearest frame {stamp} s) ===")
        width = max(len(k) for k in keys) + 2
        header = "".ljust(width) + "".join(r["label"].rjust(20) for r in results)
        print(header)
        for k in keys:
            row = k.ljust(width)
            for r in results:
                block = r["times"].get(stamp) or list(r["times"].values())[
                    int(np.argmin([abs(float(s) - t_target)
                                   for s in r["times"]]))]
                value = block.get(k)
                row += ("n/a" if value is None else f"{value:.6g}").rjust(20)
            print(row)

    if args.plot:
        plot_top_region(args.run, args.plot,
                        args.plot_time if args.plot_time is not None
                        else args.times[-1], args.t_wall, args.top_cells)

    if args.json:
        os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
        with open(args.json, "w") as handle:
            json.dump(results, handle, indent=1)
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
