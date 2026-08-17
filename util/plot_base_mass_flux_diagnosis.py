#!/usr/bin/env python3
"""Diagnose the elevated mass flux near the lower boundary of a model_column run.

Compares the conservative numerical face mass flux `f_total` of two `.faceflux`
sidecars -- normally the current release run against an older reference run --
and answers three separate questions that the cell-centred `rho*V` panel of an
evolution movie cannot distinguish:

  1. is the base elevation present at t = 0?  A nonzero t = 0 face flux in a
     hydrostatic column is a scheme defect, not a boundary condition.  The
     signature of the historical missing MUSCL-Hancock predictor source is a
     t = 0 face velocity of exactly 0.5 * dt * g.
  2. is it a single-face artifact or a broad gradient?  Panel `|f_diff|` isolates
     the dissipative part, which spikes at the first interior face when the inner
     ghost closure and cell 0 disagree.
  3. can the gradient relax at all?  The bottom-right panel measures the per-step
     mass-density increment `dt * (-df/ds)` in float32 ULPs of `rho`.  Below half
     a ULP the update `rho += dt*rhs` rounds back to `rho` every step and the
     density is frozen at representation round-off, so a flux gradient there can
     never be worked off no matter how long the run.

Usage:
    python util/plot_base_mass_flux_diagnosis.py \
        LABEL_A=run_a.txt LABEL_B=run_b.txt [output.png]

Each run argument names the base `.txt`; the `.faceflux` sidecar is appended.
Output defaults to visualization/model_column/base_mass_flux_diagnosis.png.
Set DT_SECONDS to the run's timestep (default 0.0056182, the N500/R4 CFL 0.50
value) -- it only scales the ULP panel.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Column indices fixed by the .faceflux `# columns=` header; see the io_formats
# Doxygen page. Checked at read time so a format change fails loudly.
C_CELL_KM, C_FACE_KM, C_RHO, C_V, C_FDIFF, C_FTOTAL = 0, 1, 2, 3, 15, 16


def load_faceflux(path):
    """-> (times, frame_reader). frame_reader(i) returns the i-th capture."""
    times, offsets = [], []
    with open(path) as fh:
        off = 0
        for line in fh:
            if line.startswith("# columns="):
                cols = line.split("=", 1)[1].split()
                expect = {C_FACE_KM: "face_km", C_FTOTAL: "f_total",
                          C_FDIFF: "f_diff", C_RHO: "rho_cell", C_V: "v_cell"}
                for idx, name in expect.items():
                    if cols[idx] != name:
                        raise ValueError(f"{path}: column {idx} is {cols[idx]!r}, "
                                         f"expected {name!r}")
            elif line.startswith("# t ="):
                times.append(float(line.split("=")[1].split()[0]))
                offsets.append(off)
            off += len(line)
    if not times:
        raise ValueError(f"{path}: no capture frames found")

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


def cell_widths(face_km):
    """Cell widths [m] from the upper-face heights, base extrapolated."""
    h = face_km * 1.0e3
    ds = np.empty(h.size)
    ds[0] = h[1] - h[0]          # base cell: mirror the first interior spacing
    ds[1:] = np.diff(h)
    return ds


def main():
    args = [a for a in sys.argv[1:] if "=" in a]
    if len(args) != 2:
        raise SystemExit(__doc__)
    runs = []
    for spec, color in zip(args, ("C0", "C3")):
        label, path = spec.split("=", 1)
        runs.append((label, path + ".faceflux", color))
    rest = [a for a in sys.argv[1:] if "=" not in a]
    out = rest[0] if rest else \
        "visualization/model_column/base_mass_flux_diagnosis.png"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    dt = float(os.environ.get("DT_SECONDS", "0.0056182"))

    fig, axes = plt.subplots(2, 3, figsize=(17, 8.5))
    for col, (label, path, color) in enumerate(runs):
        times, frame = load_faceflux(path)

        ax = axes[0, col]
        for target, style, alpha in ((0, ":", 0.9), (2000, "--", 0.9),
                                     (times[-1], "-", 1.0)):
            i = int(np.abs(times - target).argmin())
            d = frame(i)
            ax.plot(d[:, C_FACE_KM], d[:, C_FTOTAL], style, color=color,
                    alpha=alpha, lw=1.4, label=f"t={times[i]:.0f} s")
        ax.set_title(f"{label}\nconservative face mass flux $f_{{total}}$",
                     fontsize=10)
        ax.set_xlabel("face height [km]")
        ax.set_ylabel(r"$f_{total}$ [kg m$^{-2}$ s$^{-1}$]")
        ax.axhline(0, c="k", lw=0.5)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        d = frame(len(times) - 1)
        ax = axes[1, col]
        ax.plot(d[:, C_FACE_KM], d[:, C_FTOTAL], "-o", ms=3, color=color,
                label=r"$f_{total}$ (face)")
        ax.plot(d[:, C_CELL_KM], d[:, C_RHO] * d[:, C_V], "-s", ms=3, color="C7",
                label=r"$\rho V$ (cell centre)")
        lo = d[0, C_FACE_KM]
        ax.set_xlim(lo - 1.0, lo + 39.0)
        ax.set_xlabel("height [km]")
        ax.set_ylabel("mass flux")
        ax.set_title(f"base zoom, t = {times[-1]:.0f} s", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Right column: dissipation localization and the float32 mass-update floor,
    # both for the FIRST run given (normally the current release run).
    times, frame = load_faceflux(runs[0][1])
    d = frame(len(times) - 1)
    fdiff = np.abs(d[:, C_FDIFF])
    interior = np.median(fdiff[5:-5])
    ax = axes[0, 2]
    ax.semilogy(d[:, C_FACE_KM], fdiff, lw=1.2, color="C2")
    ax.set_title(f"{runs[0][0]}: $|f_{{diff}}|$ at t = {times[-1]:.0f} s\n"
                 f"face 0 is {fdiff[0]/interior:.0f}x the interior median",
                 fontsize=10)
    ax.set_xlabel("face height [km]")
    ax.set_ylabel(r"$|f_{diff}|$")
    ax.grid(alpha=0.3)

    ds = cell_widths(d[:, C_FACE_KM])
    ns = d.shape[0]
    rhs = np.empty(ns)
    rhs[0] = np.nan
    rhs[1:] = -(d[1:, C_FTOTAL] - d[:-1, C_FTOTAL]) / ds[1:]
    ulp = np.abs(np.spacing(d[:, C_RHO].astype(np.float32))).astype(float)
    ratio = np.abs(dt * rhs) / ulp
    frozen = int(np.sum(ratio[1:] <= 0.5))
    ax = axes[1, 2]
    ax.semilogy(d[1:, C_CELL_KM], ratio[1:], lw=1.2, color="C4")
    ax.axhline(0.5, c="r", ls="--", lw=1.2,
               label="0.5 ULP: below this the\nmass update is quantized away")
    ax.set_title(f"{runs[0][0]}: per-step mass increment / float32 ULP"
                 r"($\rho$)" + f"\n{frozen} of {ns - 1} cells cannot accumulate mass",
                 fontsize=10)
    ax.set_xlabel("height [km]")
    ax.set_ylabel(r"$|\Delta t\,\partial_t\rho|\,/\,\mathrm{ULP}(\rho)$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle("model_column base-region mass flux: what changed and what remains",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    print(f"  face-0 |f_diff| / interior median = {fdiff[0]/interior:.0f}x")
    print(f"  cells with per-step mass increment <= 0.5 ULP: {frozen}/{ns - 1}"
          f"  (median ratio {np.median(ratio[1:]):.3f})")


if __name__ == "__main__":
    main()
