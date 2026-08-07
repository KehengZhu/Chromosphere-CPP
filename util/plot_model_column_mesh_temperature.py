#!/usr/bin/env python3
"""Plot the current refined model_column mesh and its initial temperature.

The default configuration is the production-shaped coarse-equivalent N=2000
column: 1600--2153 km, with a 4x outer refinement beginning 500 km above the
base and a minimum 20 km geometric transition.  Temperature is read from the
EOS-aware ``.gamma_diag`` sidecar at t=0.

Usage:
    .venv/bin/python util/plot_model_column_mesh_temperature.py
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/chromosphere2026-mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_RUN = Path("outputs/model_column/realtime_2000s_cfl050.txt")
DEFAULT_OUTPUT = Path(
    "visualization/model_column/model_column_ns2000_mesh_temperature_initial"
)
MESH_MAX_RATIO = 1.1


def _append_uniform(faces: list[float], start: float, end: float, target_ds: float) -> None:
    """Match scenarios/mesh.cpp::append_uniform (coordinates in metres)."""
    span = end - start
    count = max(1, math.ceil(span / max(target_ds, 1.0e-30) - 1.0e-9))
    width = span / count
    faces.extend(start + width * j for j in range(1, count + 1))
    faces[-1] = end


def _grade(fine_ds: float, coarse_ds: float, minimum_width: float) -> tuple[int, float, float]:
    """Match scenarios/mesh.cpp::make_grade."""
    factor = coarse_ds / fine_ds
    count = max(1, math.ceil(math.log(factor) / math.log(MESH_MAX_RATIO)))

    def width_of(n: int) -> float:
        ratio = factor ** (1.0 / n)
        return fine_ds * (factor - 1.0) / (ratio - 1.0)

    while width_of(count) < minimum_width and count < 1_000_000:
        count += 1
    ratio = factor ** (1.0 / count)
    return count, ratio, width_of(count)


def build_outer_faces(
    length_m: float,
    ns_coarse: int,
    refine_factor: float,
    refine_start_m: float,
    transition_min_m: float,
) -> np.ndarray:
    """Build the exact outer-refinement faces used by the production column."""
    coarse_ds = length_m / ns_coarse
    fine_ds = coarse_ds / refine_factor
    count, ratio, transition_width = _grade(fine_ds, coarse_ds, transition_min_m)
    transition_start = refine_start_m - transition_width
    if not (transition_start > 0.0 and refine_start_m < length_m):
        raise ValueError("outer refinement and transition do not fit inside the domain")

    faces = [0.0]
    _append_uniform(faces, 0.0, transition_start, coarse_ds)

    # Coarse-to-fine geometric grade; equivalent to append_grade(..., false).
    coefficient = fine_ds / (ratio - 1.0)
    for j in range(1, count + 1):
        from_fine_end = coefficient * (ratio ** (count - j) - 1.0)
        faces.append(transition_start + transition_width - from_fine_end)
    faces[-1] = refine_start_m

    _append_uniform(faces, refine_start_m, length_m, fine_ds)
    return np.asarray(faces)


def read_initial_temperature(run_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read cell-centre heights [km] and T [K] from the t=0 gamma sidecar."""
    diag_path = Path(str(run_path) + ".gamma_diag")
    with diag_path.open() as stream:
        ns, ncol = (int(value) for value in stream.readline().split())
        height_km = np.fromstring(stream.readline(), sep=" ")
        rows: list[list[float]] = []
        in_initial_frame = False
        for line in stream:
            if line.startswith("# t ="):
                if in_initial_frame:
                    break
                time_s = float(line.split("t =", 1)[1].split("step", 1)[0])
                if not np.isclose(time_s, 0.0):
                    raise ValueError(f"first diagnostic frame is t={time_s:g} s, not t=0")
                in_initial_frame = True
            elif in_initial_frame and not line.startswith("#"):
                values = [float(value) for value in line.split()]
                if len(values) == ncol:
                    rows.append(values)

    data = np.asarray(rows)
    if len(height_km) != ns or data.shape != (ns, ncol):
        raise ValueError(
            f"incomplete t=0 frame in {diag_path}: heights={len(height_km)}, data={data.shape}"
        )
    return height_km, data[:, 2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="output stem; .png and .pdf are written")
    parser.add_argument("--ns-coarse", type=int, default=2000)
    parser.add_argument("--base-height-km", type=float, default=1600.0)
    parser.add_argument("--domain-length-km", type=float, default=553.0)
    parser.add_argument("--refine-factor", type=float, default=4.0)
    parser.add_argument("--refine-start-km", type=float, default=500.0,
                        help="distance above the domain base")
    parser.add_argument("--transition-km", type=float, default=20.0,
                        help="minimum geometric-transition width")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.with_suffix("")
    output.parent.mkdir(parents=True, exist_ok=True)

    printed_face_km, temperature_k = read_initial_temperature(args.run)
    faces_m = build_outer_faces(
        args.domain_length_km * 1000.0,
        args.ns_coarse,
        args.refine_factor,
        args.refine_start_km * 1000.0,
        args.transition_km * 1000.0,
    )
    ds_m = np.diff(faces_m)
    height_km = args.base_height_km + 0.5 * (faces_m[:-1] + faces_m[1:]) * 1.0e-3
    if len(height_km) != len(temperature_k):
        raise ValueError(
            f"mesh has {len(height_km)} cells but {args.run}.gamma_diag has "
            f"{len(temperature_k)}; check the mesh arguments"
        )
    # The state header is documented as height but stores cumulative UPPER faces,
    # printed at six significant digits. Temperature remains cell-centred, so use
    # reconstructed centres for plotting and faces only for this provenance check.
    # Mirror chromo_main's float32, kilometre-by-kilometre accumulation before
    # comparing with its six-significant-digit text header.
    face_height_km = np.float32(args.base_height_km) + np.cumsum(
        ds_m.astype(np.float32) * np.float32(1.0e-3), dtype=np.float32
    )
    face_error_m = np.max(np.abs(face_height_km - printed_face_km)) * 1000.0
    if face_error_m > 5.1:
        raise ValueError(f"reconstructed mesh disagrees with run header by {face_error_m:.2f} m")

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig, axes = plt.subplots(
        2, 1, figsize=(89 / 25.4, 98 / 25.4), sharex=True,
        gridspec_kw={"height_ratios": [0.9, 1.25], "hspace": 0.12},
    )
    mesh_color = "#3C6E9B"
    temperature_color = "#C44E52"
    # Locate the start of the monotone grade from the first cell narrower than
    # the lower uniform plateau; this avoids exposing an implementation detail in labels.
    first_grade = int(np.flatnonzero(ds_m < ds_m[0] * (1.0 - 1.0e-10))[0])
    transition_lo = args.base_height_km + faces_m[first_grade] * 1.0e-3
    refine_lo = args.base_height_km + args.refine_start_km
    domain_top = args.base_height_km + args.domain_length_km

    for axis in axes:
        axis.axvspan(transition_lo, refine_lo, color="#D9E6F2", alpha=0.65, lw=0)
        axis.axvspan(refine_lo, domain_top, color="#B8D3E8", alpha=0.45, lw=0)
        axis.axvline(refine_lo, color="#6689A6", lw=0.7, ls="--")
        axis.grid(axis="y", color="0.88", lw=0.55)
        axis.tick_params(direction="out", length=2.5, width=0.7)

    axes[0].plot(height_km, ds_m, color=mesh_color, lw=1.35)
    axes[0].set_ylabel(r"Cell width $\Delta s$ (m)")
    axes[0].set_ylim(0, ds_m.max() * 1.12)
    axes[0].text(0.01, 0.90, "a", transform=axes[0].transAxes,
                 fontsize=8, fontweight="bold")
    axes[0].text(0.03, 0.15, f"coarse {ds_m[0]:.1f} m",
                 color=mesh_color, transform=axes[0].transAxes,
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.0})
    axes[0].text(0.97, 0.15, f"fine {ds_m[-1]:.1f} m",
                 color=mesh_color, ha="right", transform=axes[0].transAxes,
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.0})

    axes[1].plot(height_km, temperature_k / 1000.0, color=temperature_color, lw=1.45)
    axes[1].set_ylabel("Initial temperature (kK)")
    axes[1].set_xlabel("Height (km)")
    axes[1].set_xlim(args.base_height_km, domain_top)
    axes[1].text(0.01, 0.90, "b", transform=axes[1].transAxes,
                 fontsize=8, fontweight="bold")
    axes[1].text(0.97, 0.08, "outer refined region", color="#557A96",
                 ha="right", transform=axes[1].transAxes)

    fig.suptitle(
        rf"Model column: $N_{{\mathrm{{coarse}}}}={args.ns_coarse}$ "
        rf"({len(ds_m)} cells after refinement)",
        fontsize=8, y=0.995,
    )
    fig.align_ylabels(axes)
    fig.subplots_adjust(left=0.19, right=0.98, bottom=0.12, top=0.93)
    fig.savefig(output.with_suffix(".png"), dpi=600)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(
        f"wrote {output}.png and {output}.pdf; cells={len(ds_m)}, "
        f"ds=[{ds_m.min():.3f}, {ds_m.max():.3f}] m, "
        f"T=[{temperature_k.min():.2f}, {temperature_k.max():.2f}] K, "
        f"header_face_error_max={face_error_m:.3f} m"
    )


if __name__ == "__main__":
    main()
