#!/usr/bin/env python3
"""Compare split conduction-to-hydro diagnostic CSVs at a common physical time.

The CSV is sampled every N Stage-D calls, so this script linearly interpolates
between the two diagnostic samples which bracket ``--time`` before comparing
the ns=2000/3000/4000 profiles. It plots the time-normalized conduction pressure
change, the total momentum-RHS response, and the total mass flux near 250 km.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

# Keep matplotlib and fontconfig cache files out of the repository in restricted
# environments.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/chromosphere-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/chromosphere-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_FILES = [
    "outputs/_archive/iso_cond_hydro_local4_outer2000.csv",
    "outputs/_archive/iso_cond_hydro_local4_outer3000.csv",
    "outputs/_archive/iso_cond_hydro_local4_outer4000.csv",
]

METRICS = (
    ("dp_cond_rate", r"$\Delta p_{\rm cond}/\Delta t$ [Pa s$^{-1}$]"),
    ("dR_mom_total", r"$\Delta R_{\rm mom,i}+\Delta R_{\rm mom,n}$ [kg m$^{-2}$ s$^{-2}$]"),
    ("mass_flux_total", r"$\rho_i V+\rho_n U$ [kg m$^{-2}$ s$^{-1}$]"),
)


def csv_rows(path: Path):
    """Yield data rows while skipping the documented comment preamble."""
    with path.open(newline="") as stream:
        yield from csv.DictReader(line for line in stream if not line.startswith("#"))


def sample_times(path: Path) -> list[tuple[float, int]]:
    """Return one (time, Stage-D index) pair for each full diagnostic sample."""
    samples: list[tuple[float, int]] = []
    previous_stage: int | None = None
    for row in csv_rows(path):
        stage = int(row["stage"])
        if stage != previous_stage:
            samples.append((float(row["time_s"]), stage))
            previous_stage = stage
    if not samples:
        raise ValueError(f"no diagnostic rows in {path}")
    return samples


def time_bracket(samples: list[tuple[float, int]], target: float) -> tuple[tuple[float, int], tuple[float, int]]:
    """Find the two samples that bracket target, or the exact sample twice."""
    samples = sorted(samples)
    if target < samples[0][0] or target > samples[-1][0]:
        raise ValueError(
            f"target t={target:g} s is outside [{samples[0][0]:.6g}, {samples[-1][0]:.6g}] s"
        )
    for current, following in zip(samples, samples[1:]):
        if current[0] == target:
            return current, current
        if current[0] <= target <= following[0]:
            return current, following
    return samples[-1], samples[-1]


def read_stage_profiles(path: Path, stages: set[int]) -> dict[int, dict[str, np.ndarray]]:
    """Read only the requested full snapshots, not the entire large CSV into RAM."""
    columns = {stage: {"height_km": [], "cell_width_km": [],
                       "dp_cond_rate": [], "dR_mom_total": [], "mass_flux_total": []}
               for stage in stages}
    for row in csv_rows(path):
        stage = int(row["stage"])
        if stage not in columns:
            continue
        dt = float(row["dt_s"])
        if dt <= 0.0:
            raise ValueError(f"non-positive dt at Stage-D sample {stage} in {path}")
        columns[stage]["height_km"].append(float(row["height_km"]))
        # Local cell width [km]. `cell_width_m` is emitted by the refined-mesh
        # diagnostic; fall back to a NaN placeholder for legacy CSVs so a uniform
        # median is used instead (see load_run).
        columns[stage]["cell_width_km"].append(
            float(row["cell_width_m"]) * 1.0e-3 if "cell_width_m" in row else float("nan")
        )
        columns[stage]["dp_cond_rate"].append(float(row["dp_cond_total_Pa"]) / dt)
        columns[stage]["dR_mom_total"].append(
            float(row["dR_mom_i"]) + float(row["dR_mom_n"])
        )
        columns[stage]["mass_flux_total"].append(
            float(row["mass_flux_i"]) + float(row["mass_flux_n"])
        )

    profiles: dict[int, dict[str, np.ndarray]] = {}
    for stage, values in columns.items():
        if not values["height_km"]:
            raise ValueError(f"missing selected Stage-D sample {stage} in {path}")
        profiles[stage] = {name: np.asarray(value, dtype=float) for name, value in values.items()}
    return profiles


def interpolate_time(profile0: dict[str, np.ndarray], profile1: dict[str, np.ndarray], weight: float) -> dict[str, np.ndarray]:
    """Interpolate two snapshots at identical cell centers to a physical time."""
    h0 = profile0["height_km"]
    h1 = profile1["height_km"]
    if not np.allclose(h0, h1, rtol=0.0, atol=1.0e-7):
        raise ValueError("cell-center heights changed within one diagnostic CSV")
    result = {"height_km": h0, "cell_width_km": profile0["cell_width_km"]}
    for metric, _ in METRICS:
        result[metric] = (1.0 - weight) * profile0[metric] + weight * profile1[metric]
    return result


def resolution_from_path(path: Path, count: int) -> int:
    # Coarse-equivalent ISO_NS from the filename: ns#### (legacy) or outer#### (the
    # refined iso_cond_hydro_local4_outer#### naming).
    match = re.search(r"(?:ns|outer)(\d+)", path.stem)
    return int(match.group(1)) if match else count


def load_run(path: Path, target_time: float) -> dict[str, object]:
    samples = sample_times(path)
    low, high = time_bracket(samples, target_time)
    selected = {low[1], high[1]}
    profiles = read_stage_profiles(path, selected)
    weight = 0.0 if low[0] == high[0] else (target_time - low[0]) / (high[0] - low[0])
    profile = interpolate_time(profiles[low[1]], profiles[high[1]], weight)
    h = profile["height_km"]
    # Fit against the LOCAL cell width at 250 km (the diagnostic altitude), which
    # on a refined mesh is the fine spacing there — NOT the median domain spacing
    # (a refined run's median is dominated by the coarse column and would misstate
    # the resolution of the artifact). Legacy CSVs without cell_width_m fall back
    # to the median difference of cell-center heights.
    cw = profile["cell_width_km"]
    if np.all(np.isfinite(cw)):
        ds_km = float(np.interp(250.0, h, cw))
    else:
        ds_km = float(np.median(np.diff(h)))
    return {
        "path": path,
        "ns": resolution_from_path(path, len(h)),
        "ds_km": ds_km,
        "t_low": low[0],
        "t_high": high[0],
        "stage_low": low[1],
        "stage_high": high[1],
        "profile": profile,
    }


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Return intercept, slope, and R² for y = intercept + slope*x."""
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (intercept + slope * x)
    total = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 if total == 0.0 else 1.0 - np.sum(residual ** 2) / total
    return float(intercept), float(slope), float(r_squared)


def plot(runs: list[dict[str, object]], target_time: float, hmin: float, hmax: float, out_path: Path) -> None:
    colors = ["C0", "C1", "C2", "C3", "C4"]
    runs = sorted(runs, key=lambda item: float(item["ds_km"]), reverse=True)
    fig, axes = plt.subplots(len(METRICS), 2, figsize=(13.0, 10.4), constrained_layout=True)
    fig.suptitle(
        f"Split conduction-to-hydro diagnostic at t={target_time:g} s (linear time interpolation)",
        fontsize=14,
    )

    for row, (metric, ylabel) in enumerate(METRICS):
        ax_profile, ax_fit = axes[row]
        ds_values: list[float] = []
        values_250: list[float] = []
        for color, run in zip(colors, runs):
            profile = run["profile"]
            assert isinstance(profile, dict)
            h = profile["height_km"]
            values = profile[metric]
            mask = (h >= hmin) & (h <= hmax)
            ns = int(run["ns"])
            ds_km = float(run["ds_km"])
            ax_profile.plot(h[mask], values[mask], color=color, linewidth=1.6,
                            label=fr"ns={ns}, $\Delta s$={ds_km:.3f} km")
            ds_values.append(ds_km)
            values_250.append(float(np.interp(250.0, h, values)))

        x = np.asarray(ds_values)
        y = np.asarray(values_250)
        intercept, slope, r_squared = linear_fit(x, y)
        xfit = np.linspace(0.0, 1.05 * np.max(x), 100)
        ax_fit.plot(x, y, "o", color="C0", markersize=6, label="interpolated diagnostic")
        ax_fit.plot(xfit, intercept + slope * xfit, "--", color="C0", linewidth=1.3,
                    label=fr"fit: {intercept:.3e} + {slope:.3e}$\Delta s$")
        ax_fit.scatter([0.0], [intercept], marker="*", s=110, color="crimson", zorder=4,
                       label=fr"$\Delta s\to0$: {intercept:.3e}")
        ax_fit.axhline(0.0, color="0.35", linewidth=0.7)
        ax_fit.set_xlim(left=0.0)
        ax_fit.set_xlabel(r"local cell width at 250 km $\Delta s$ [km]")
        ax_fit.set_ylabel(f"value at 250 km\n{ylabel}")
        ax_fit.set_title(fr"250 km linear fit ($R^2$={r_squared:.4f})", fontsize=10)
        ax_fit.legend(fontsize=8, loc="best")
        ax_fit.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

        ax_profile.axvline(250.0, color="0.35", linewidth=0.8, linestyle=":")
        ax_profile.axhline(0.0, color="0.35", linewidth=0.7)
        ax_profile.set_xlim(hmin, hmax)
        ax_profile.set_xlabel("height [km]")
        ax_profile.set_ylabel(ylabel)
        ax_profile.set_title("matched-time lower-chromosphere profiles", fontsize=10)
        ax_profile.legend(fontsize=8, loc="best")
        ax_profile.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

        print(f"{metric} @ 250 km: intercept={intercept:.9e}, slope={slope:.9e} per km, R2={r_squared:.8f}")
        for run, value in zip(runs, y):
            print(f"  ns={int(run['ns'])}: ds={float(run['ds_km']):.9f} km, value={value:.9e}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    print(f"wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="*", type=Path, default=[Path(path) for path in DEFAULT_FILES],
                        help="diagnostic CSV files (default: ns=2000/3000/4000 files under outputs/)")
    parser.add_argument("--time", type=float, default=100.0, help="matched physical time in seconds")
    parser.add_argument("--hmin", type=float, default=150.0, help="profile-panel minimum height in km")
    parser.add_argument("--hmax", type=float, default=350.0, help="profile-panel maximum height in km")
    parser.add_argument("--out", type=Path, default=Path("visualization/_archive/iso_cond_hydro_convergence_t100.png"),
                        help="output PNG path")
    args = parser.parse_args()

    if len(args.files) < 2:
        parser.error("provide at least two diagnostic CSV files")
    missing = [str(path) for path in args.files if not path.is_file() or path.stat().st_size == 0]
    if missing:
        parser.error("missing or empty diagnostic CSV: " + ", ".join(missing))

    runs = [load_run(path, args.time) for path in args.files]
    for run in sorted(runs, key=lambda item: int(item["ns"])):
        print(
            f"ns={int(run['ns'])}: ds={float(run['ds_km']):.9f} km; "
            f"bracket t=[{float(run['t_low']):.6f}, {float(run['t_high']):.6f}] s "
            f"at stages [{int(run['stage_low'])}, {int(run['stage_high'])}]"
        )
    plot(runs, args.time, args.hmin, args.hmax, args.out)


if __name__ == "__main__":
    main()
