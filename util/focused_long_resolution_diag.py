#!/usr/bin/env python3
"""Focused N500/R4 vs N1000/R4 long-run diagnostic for model_column.

Streams the large release sidecars so the 4000 s comparison can be evaluated
without loading every frame into memory. Metric definitions match the existing
face_flux_diag.py / outer_refine_diag.py diagnostics.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

TARGETS = (500.0, 1000.0, 2000.0, 4000.0)
WINDOW = (2130.0, 2150.0)
FF_COLUMNS = ("cell_km face_km rho_cell v_cell T_cell rho_L rho_R v_L v_R T_L T_R "
              "cs_L cs_R a_face f_central f_diff f_total "
              "r_rho phi_plus_rho r_ip1_rho phi_minus_rho r_v phi_plus_v "
              "r_T phi_plus_T").split()
# RELEASE `.outercond` layout (physical conduction only).
OUTER_COLS = ["t", "step", "T_top", "T_wall", "kappa_face", "ds_face",
              "area_ratio", "q_face"]
# Historical layout, when the conduction operator still carried an artificial
# diffusivity. Rows in this shape are mapped onto the release names so old runs
# stay analysable: kappa_face <- kappa_phys_face, q_face <- q_total.
LEGACY_OUTER_COLS = ["t", "step", "T_top", "T_wall", "kappa_phys_face",
                     "chi_num_face", "kappa_num_face", "ds_face", "area_ratio",
                     "q_phys", "q_num", "q_total"]
LEGACY_TO_RELEASE = [0, 1, 2, 3, 4, 7, 8, 11]


def roughness(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.mean(np.abs(np.diff(x, 2))) / max(np.mean(np.abs(x)), 1e-300))


def relstd(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.std(x) / max(abs(np.mean(x)), 1e-300))


def outer_boundary_layer(h_km: np.ndarray, T: np.ndarray, frac: float = 0.1) -> dict:
    s = np.asarray(h_km) * 1e3
    G = np.abs(np.gradient(np.asarray(T), s))
    gmax = float(G.max())
    i = len(G) - 1
    while i > 0 and G[i - 1] >= frac * gmax:
        i -= 1
    return {
        "bl_thickness_km": float(h_km[-1] - h_km[i]),
        "bl_n_cells": int(len(G) - i),
        "bl_edge_km": float(h_km[i]),
        "G_max_K_per_m": gmax,
        "h_Gmax_km": float(h_km[int(np.argmax(G))]),
    }


def update_best(best: dict, t: float, value) -> None:
    for target in TARGETS:
        d = abs(t - target)
        if target not in best or d < best[target][0]:
            best[target] = (d, t, value)


def stream_gamma(base: Path) -> tuple[int, np.ndarray, dict]:
    path = Path(str(base) + ".gamma_diag")
    best = {}
    with path.open() as fh:
        ns, ncol = map(int, fh.readline().split())
        h = np.fromstring(fh.readline(), sep=" ")
        cur_t = None
        rows = []

        def flush():
            nonlocal rows
            if cur_t is None or len(rows) != ns:
                rows = []
                return
            arr = np.asarray(rows, dtype=float)
            m = (h >= WINDOW[0]) & (h <= WINDOW[1])
            v = arr[m, 1]
            neg_h = h[m][v < 0]
            metrics = {
                "mean_V_ms": float(v.mean()),
                "std_V_ms": float(v.std()),
                "min_V_ms": float(v.min()),
                "max_V_ms": float(v.max()),
                "Q_V": roughness(v),
                "n_negative_V": int(np.count_nonzero(v < 0)),
                "negative_V_heights_km": [float(x) for x in neg_h],
                "T_top_K": float(arr[-1, 2]),
                "p_top_Pa": float(arr[-1, 6]),
                "max_abs_V_global_ms": float(np.max(np.abs(arr[:, 1]))),
            }
            metrics.update(outer_boundary_layer(h, arr[:, 2]))
            update_best(best, cur_t, metrics)
            rows = []

        for line in fh:
            if line.startswith("# t ="):
                flush()
                cur_t = float(line.split("t =")[1].split("step")[0])
            elif line.startswith("#"):
                continue
            elif cur_t is not None:
                vals = np.fromstring(line, sep=" ")
                if vals.size == ncol:
                    rows.append(vals)
        flush()
    return ns, h, {t: {"t_gamma": best[t][1], **best[t][2]} for t in TARGETS}


def parse_ff_meta(line: str, meta: dict) -> None:
    for tok in line[1:].split():
        if "=" in tok and not tok.startswith("columns"):
            k, v = tok.split("=", 1)
            meta.setdefault(k, v)


def stream_faceflux(base: Path) -> tuple[dict, dict, dict]:
    path = Path(str(base) + ".faceflux")
    meta = {}
    best = {}
    ref = None
    mass0 = None
    geometry = None
    cur_t = None
    cur_step = None
    rows = []

    def flush():
        nonlocal rows, ref, mass0, geometry
        if cur_t is None or not rows:
            rows = []
            return
        arr = np.asarray(rows, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != len(FF_COLUMNS):
            rows = []
            return
        c = {name: arr[:, k] for k, name in enumerate(FF_COLUMNS)}
        if ref is None:
            ref = c["f_total"].copy()
            faces = c["face_km"]
            cells = c["cell_km"]
            widths_km = np.concatenate(([2.0 * (faces[0] - cells[0])], np.diff(faces)))
            mass0 = float(np.sum(c["rho_cell"] * widths_km * 1e3))
            geometry = {
                "capture_rows": int(len(cells)),
                "min_width_km": float(widths_km.min()),
                "max_width_km": float(widths_km.max()),
                "top_width_km": float(widths_km[-1]),
            }
        if ref is None or len(ref) != len(c["f_total"]):
            rows = []
            return
        mc = (c["cell_km"] >= WINDOW[0]) & (c["cell_km"] <= WINDOW[1])
        mf = (c["face_km"] >= WINDOW[0]) & (c["face_km"] <= WINDOW[1])
        M = c["rho_cell"] * c["v_cell"]
        f_eff = c["f_total"] - ref
        widths_km = np.concatenate(([2.0 * (c["face_km"][0] - c["cell_km"][0])],
                                    np.diff(c["face_km"])))
        mass = float(np.sum(c["rho_cell"] * widths_km * 1e3))
        metrics = {
            "t_faceflux": float(cur_t),
            "step_faceflux": int(cur_step),
            "M_mean": float(M[mc].mean()),
            "A_M": relstd(M[mc]),
            "Q_M": roughness(M[mc]),
            "F_eff_mean": float(f_eff[mf].mean()),
            "A_eff": relstd(f_eff[mf]),
            "Q_eff": roughness(f_eff[mf]),
            "column_mass_kgm2": mass,
            "column_mass_rel_change": float(mass / mass0 - 1.0),
            "n_cells_window": int(mc.sum()),
            "n_faces_window": int(mf.sum()),
        }
        update_best(best, cur_t, metrics)
        rows = []

    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("# t "):
                    flush()
                    parts = line.split()
                    cur_t, cur_step = float(parts[3]), int(parts[6])
                else:
                    parse_ff_meta(line, meta)
                continue
            vals = np.fromstring(line, sep=" ")
            if vals.size == len(FF_COLUMNS):
                rows.append(vals)
        flush()
    out = {t: best[t][2] for t in TARGETS}
    return meta, geometry, out


def stream_outercond(base: Path) -> dict:
    """Interpolate boundary quantities and integrate q_face to exact target times."""
    path = Path(str(base) + ".outercond")
    out = {}
    prev = None
    integral = 0.0
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            vals = np.fromstring(line, sep=" ")
            if vals.size == len(LEGACY_OUTER_COLS):
                vals = vals[LEGACY_TO_RELEASE]
            elif vals.size != len(OUTER_COLS):
                continue
            row = dict(zip(OUTER_COLS, map(float, vals)))
            if prev is None:
                prev = row
                continue
            t0, t1 = prev["t"], row["t"]
            q0, q1 = prev["q_face"], row["q_face"]
            dt = t1 - t0
            if dt <= 0:
                prev = row
                continue
            for target in TARGETS:
                if target in out or not (t0 <= target <= t1):
                    continue
                a = (target - t0) / dt
                interp = {name: float(prev[name] + a * (row[name] - prev[name]))
                          for name in OUTER_COLS}
                q_target = interp["q_face"]
                interp["t"] = float(target)
                interp["E_phys_Jm2"] = float(
                    integral + 0.5 * (q0 + q_target) * (target - t0))
                out[target] = interp
            integral += 0.5 * (q0 + q1) * dt
            prev = row
    for target in TARGETS:
        if target not in out and prev is not None and abs(prev["t"] - target) <= 1.0:
            tail = dict(prev)
            tail["E_phys_Jm2"] = float(integral)
            out[target] = tail
    missing = [t for t in TARGETS if t not in out]
    if missing:
        raise SystemExit(f"outercond does not bracket targets {missing} in {path}")
    return out


def analyze(label: str, base: str) -> dict:
    p = Path(base)
    ns, _, gamma = stream_gamma(p)
    meta, geometry, ff = stream_faceflux(p)
    oc = stream_outercond(p)
    out = {"label": label, "base": base, "ns_actual": ns, "faceflux_meta": meta,
           "geometry": geometry, "times": {}}
    for t in TARGETS:
        out["times"][str(int(t))] = {**gamma[t], **ff[t], **oc[t]}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True, help="LABEL=base-path")
    ap.add_argument("--json")
    args = ap.parse_args()
    runs = []
    for spec in args.run:
        label, base = spec.split("=", 1)
        runs.append(analyze(label, base))
    result = {"window_km": WINDOW, "targets_s": TARGETS, "runs": runs}
    print(json.dumps(result, indent=2))
    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
