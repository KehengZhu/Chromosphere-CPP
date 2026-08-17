#!/usr/bin/env python3
"""Compare a float32 diagnostic model_column run against the double-precision release.

The RELEASE stores the conserved state in double; float32 is the default-off
diagnostic build (-DCHROMO_STATE_FLOAT32=ON), which reverts chromosphere::Real --
and with it the state, the mesh metrics, the Grid physical constants and sim_time --
from double to float. Everything else -- physics, grid, boundary conditions,
reconstruction, Godunov flux, CFL, diagnostics -- is identical. (Before the
cutover the option was the other way round and was named CHROMO_STATE_FLOAT64;
that option is retired, and it promoted only Vec, so several float casts survived
in the update chain.)

The question it answers: is the broad lower-chromosphere gradient in the
conservative face mass flux f_total a physical transient, or is it float32
representation round-off? In a genuine quasi-steady state on a constant-area
tube, continuity forces f_total to be independent of height.

Usage:
    python util/plot_precision_control.py f32_run.txt f64_run.txt [output.png]

The first argument is the float32 diagnostic leg, the second the double release.

Each argument names the base .txt; the .faceflux and .gamma_diag sidecars are
appended. Output defaults to
visualization/model_column/precision_control_N500_4000s.png
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_CELL_KM, C_FACE_KM, C_RHO, C_V, C_FDIFF, C_FTOTAL = 0, 1, 2, 3, 15, 16
DT = float(os.environ.get("DT_SECONDS", "0.0056182"))


def load_faceflux(path):
    times, offsets = [], []
    with open(path) as fh:
        off = 0
        for line in fh:
            if line.startswith("# columns="):
                cols = line.split("=", 1)[1].split()
                if cols[C_FACE_KM] != "face_km" or cols[C_FTOTAL] != "f_total":
                    raise ValueError(f"{path}: unexpected .faceflux column layout")
            elif line.startswith("# t ="):
                times.append(float(line.split("=")[1].split()[0]))
                offsets.append(off)
            off += len(line)

    def frame(i):
        rows = []
        with open(path) as fh:
            fh.seek(offsets[i])
            fh.readline()
            for line in fh:
                if line.startswith("#") or not line.strip():
                    break
                rows.append(line.split())
        return np.asarray(rows, dtype=float)

    return np.asarray(times), frame


def load_density(path):
    """-> (height_km, rho_first, rho_last) from a .gamma_diag sidecar."""
    lines = open(path).read().split("\n")
    ns = int(lines[0].split()[0])
    h = np.asarray([float(x) for x in lines[1].split()])
    marks = [i for i, l in enumerate(lines) if l.startswith("# t =")]

    def rho(m):
        i0 = marks[m] + 1
        return np.asarray([float(lines[i0 + j].split()[0]) for j in range(ns)])

    return h, rho(0), rho(len(marks) - 1)


def integrated_divergence(times, frame):
    """Time-integrate the recorded -d(f_total)/ds over every captured frame.

    Returns the mass density each cell should have gained, in kg m^-3, if the
    stored state had responded to the fluxes the run itself recorded. Trapezoid
    in time; face 0 has no left neighbour and is reported as zero.
    """
    d0 = frame(0)
    h = d0[:, C_FACE_KM] * 1.0e3
    ds = np.empty(h.size)
    ds[0] = h[1] - h[0]
    ds[1:] = np.diff(h)
    acc = np.zeros(h.size)
    prev = None
    for i in range(times.size):
        f = frame(i)[:, C_FTOTAL]
        rhs = np.zeros(f.size)
        rhs[1:] = -(f[1:] - f[:-1]) / ds[1:]
        if prev is not None:
            acc += 0.5 * (rhs + prev) * (times[i] - times[i - 1])
        prev = rhs
    return acc


def flux_ulp_ratio(d, dtype):
    """Per-step mass-density increment measured in ULPs of rho at `dtype`."""
    h = d[:, C_FACE_KM] * 1.0e3
    ds = np.empty(h.size)
    ds[0] = h[1] - h[0]
    ds[1:] = np.diff(h)
    rhs = np.empty(h.size)
    rhs[0] = np.nan
    rhs[1:] = -(d[1:, C_FTOTAL] - d[:-1, C_FTOTAL]) / ds[1:]
    ulp = np.abs(np.spacing(d[:, C_RHO].astype(dtype))).astype(float)
    return np.abs(DT * rhs) / ulp


def main():
    args = [a for a in sys.argv[1:]]
    if len(args) < 2:
        raise SystemExit(__doc__)
    f32_base, f64_base = args[0], args[1]
    out = args[2] if len(args) > 2 else \
        "visualization/model_column/precision_control_N500_4000s.png"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    legs = [("float32 (diagnostic storage)", f32_base, "C3", np.float32),
            ("float64 (release storage)", f64_base, "C0", np.float64)]
    data = {}
    for label, base, color, dtype in legs:
        times, frame = load_faceflux(base + ".faceflux")
        data[label] = dict(times=times, frame=frame, color=color, dtype=dtype,
                           dens=load_density(base + ".gamma_diag"))

    fig, ax = plt.subplots(2, 3, figsize=(17.5, 9))

    # (a),(b) f_total(h) at 2000 s and 4000 s.
    for col, target in enumerate((2000.0, 4000.0)):
        a = ax[0, col]
        for label, dd in data.items():
            i = int(np.abs(dd["times"] - target).argmin())
            d = dd["frame"](i)
            seg = d[:, C_FACE_KM] < 1850.0
            spread = 100.0 * (d[seg, C_FTOTAL].max() - d[seg, C_FTOTAL].min()) \
                / d[seg, C_FTOTAL].mean()
            a.plot(d[:, C_FACE_KM], d[:, C_FTOTAL], lw=1.5, color=dd["color"],
                   label=f"{label}\n  base/top = "
                         f"{d[0, C_FTOTAL] / d[-1, C_FTOTAL]:.3f}, "
                         f"1600-1850 km spread {spread:.0f}%")
        a.set_title(f"conservative face mass flux at t = {target:.0f} s\n"
                    "a quasi-steady column must be FLAT", fontsize=10)
        a.set_xlabel("face height [km]")
        a.set_ylabel(r"$f_{total}$ [kg m$^{-2}$ s$^{-1}$]")
        a.set_ylim(bottom=0.0)
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

    # (c) base zoom at the final time, face flux against cell-centred rho*V.
    a = ax[0, 2]
    for label, dd in data.items():
        d = dd["frame"](len(dd["times"]) - 1)
        a.plot(d[:, C_FACE_KM], d[:, C_FTOTAL], "-o", ms=3, color=dd["color"],
               label=f"{label}: $f_{{total}}$")
        a.plot(d[:, C_CELL_KM], d[:, C_RHO] * d[:, C_V], "--", lw=1.0,
               color=dd["color"], alpha=0.6, label=f"{label}: " + r"$\rho V$")
    a.set_xlim(1600, 1640)
    a.set_title("base zoom, final time\nthe first-face artifact survives in both",
                fontsize=10)
    a.set_xlabel("height [km]")
    a.set_ylabel("mass flux")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)

    # (d) velocity at the final time.
    a = ax[1, 0]
    for label, dd in data.items():
        d = dd["frame"](len(dd["times"]) - 1)
        a.plot(d[:, C_CELL_KM], d[:, C_V], lw=1.4, color=dd["color"], label=label)
    a.set_yscale("log")
    a.set_title("cell-centred velocity, final time", fontsize=10)
    a.set_xlabel("height [km]")
    a.set_ylabel("V [m/s]")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # (e) fractional density change over the whole run, against what the
    #     recorded face-flux divergence says it should have been. This is the
    #     column mass budget: it closes in float64 and fails by ~200x, with the
    #     wrong sign, in float32.
    a = ax[1, 1]
    for label, dd in data.items():
        h, r0, rF = dd["dens"]
        obs = (rF - r0) / r0
        pred = integrated_divergence(dd["times"], dd["frame"])[:h.size] / r0
        band = (h >= 1605.0) & (h < 1850.0)
        ratio = pred[band].mean() / obs[band].mean()
        a.plot(h, obs, lw=1.3, color=dd["color"],
               label=f"{label}: observed\n  1605-1850 km budget pred/obs = {ratio:+.1f}")
        a.plot(h[band], pred[band], ":", lw=1.6, color=dd["color"],
               label=f"{label}: implied by " + r"$\int\!-\partial_s f_{total}\,dt$")
    a.axhline(0, c="k", lw=0.5)
    a.set_yscale("symlog", linthresh=1e-5)
    a.set_title("column mass budget over the run\n"
                "float32 loses mass where its own fluxes say it should gain",
                fontsize=10)
    a.set_xlabel("height [km]")
    a.set_ylabel(r"$\Delta\rho/\rho(0)$")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)

    # (f) the mechanism: per-step mass increment in ULPs of the storage type.
    a = ax[1, 2]
    for label, dd in data.items():
        d = dd["frame"](len(dd["times"]) - 1)
        ratio = flux_ulp_ratio(d, dd["dtype"])
        frozen = int(np.sum(ratio[1:] <= 0.5))
        a.semilogy(d[1:, C_CELL_KM], ratio[1:], lw=1.2, color=dd["color"],
                   label=f"{label}: {frozen}/{d.shape[0] - 1} cells frozen")
    a.axhline(0.5, c="k", ls="--", lw=1.2,
              label="0.5 ULP: below this the\nmass update is quantized away")
    a.set_title("per-step mass increment / ULP" + r"($\rho$)" +
                " of the storage type\nthe mechanism", fontsize=10)
    a.set_xlabel("height [km]")
    a.set_ylabel(r"$|\Delta t\,\partial_t\rho|\,/\,\mathrm{ULP}(\rho)$")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    fig.suptitle("model_column N500, 4000 s, release Godunov numerics — "
                 "float32 diagnostic storage vs the double-precision release",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
