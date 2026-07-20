#!/usr/bin/env python3
"""Synthetic IRIS Doppler-shift forward model from a flaring-loop run.

Reads the SAME geometry `.dat` + state `.txt` pair as ``animate_loop_flare.py``
and, for a set of spectral lines spanning the chromosphere -> flare-corona
temperature ladder, computes the intensity-weighted line-of-sight Doppler
velocity and intensity as a function of time. This reproduces the canonical
explosive-evaporation fingerprint seen by IRIS (Fisher 1985; Milligan & Dennis
2009): the HOT lines (Fe XXI) blueshift (chromospheric evaporation, upflow)
while the COOL TR/chromospheric lines (Si IV, C II, Mg II) redshift
(chromospheric condensation, downflow).

Physics (optically-thin approximation)
--------------------------------------
Emissivity per cell:        eps(s) = n_e(s)^2 * G(T_e(s))
Contribution function:      G(T) modeled as a Gaussian in log10 T centered on
                            the line's CHIANTI peak-formation temperature,
                            width 0.3 dex (the typical ion-fraction width).
Line centroid Doppler:      v_Dopp = - sum(eps * v_z * ds) / sum(eps * ds)
                            (redshift > 0; blueshift < 0, IRIS convention)
Vertical (disk-center) LOS: v_z = V * (dh/ds), V the field-aligned flow, dh/ds
                            the vertical component of the loop tangent.

For the cool, optically-thick chromospheric lines (Mg II, C II, H-alpha) this
emission-measure weighting is a *quick-look proxy* only -- a rigorous profile
needs NLTE radiative transfer (export the atmosphere to RH1.5D). The hot,
optically-thin lines (Fe XXI and the coronal Fe channels) are well represented.

LOS note: the `.dat` carries only |B| and phi_g (-> height), so only the
VERTICAL velocity component dh/ds is known. The default is therefore a
disk-center (vertical) line of sight. A general inclined LOS would require the
3D field-line trace (azimuth), which this geometry file does not store.

Outputs
-------
  visualization/events/<stem>_doppler.png      v_Dopp(t) + intensity(t) per line
  visualization/events/<stem>_dopplergram.png  synthetic (v, t) spectrograms
  outputs/events/<stem>_doppler.csv            raw v_Dopp(t), I(t) per line

Usage
-----
  python util/synth_doppler.py <loop.dat> <state.txt> [out_stem]

Example
-------
  python util/synth_doppler.py scenarios/data/event_ensemble/c05_101Mm.dat \
      outputs/events/event_te/c05_101Mm_te.txt c05_101Mm_te
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m_i = 1.6726219e-27       # proton mass [kg]
m_n = m_i
k_b = 1.380649e-23        # Boltzmann [J/K]
g_si = 274.0              # solar surface gravity [m/s^2]
amu = 1.66053907e-27      # atomic mass unit [kg]

# Beam timing window [s], read from the same FLARE_* env vars the run used, so the
# shaded impulsive window on the plots matches the simulated beam.
T_ON  = float(os.environ.get("FLARE_T_ON", 2.0))
T_OFF = T_ON + float(os.environ.get("FLARE_DUR", 10.0))

GLOGT_SIGMA = 0.30        # contribution-function width in log10 T [dex]
VTURB_KMS   = 10.0        # micro-turbulent + instrumental floor for line widths
NV          = 161         # velocity samples for synthetic spectrograms
VWIN_KMS    = 400.0       # spectrogram half-window [km/s]

# Spectral lines: (label, rest wavelength [A], log10 peak-formation T, emitting
# ion mass [amu], optically_thick?). Formation temperatures are CHIANTI ion-
# fraction maxima. `thick=True` flags lines for which the optically-thin EM
# weighting is only a quick-look proxy (NLTE caveat printed in the legend).
LINES = [
    ("Fe XXI 1354", 1354.08, 7.05, 56.0, False),   # flare evaporation (hot, up)
    ("Fe XVI 335",   335.41, 6.43, 56.0, False),    # hot corona (AIA 335)
    ("Fe XII 193",   193.51, 6.20, 56.0, False),    # corona (AIA 193)
    ("Fe IX 171",    171.07, 5.80, 56.0, False),    # corona (AIA 171)
    ("O IV 1401",   1401.16, 5.18, 16.0, False),    # TR
    ("Si IV 1394",  1393.76, 4.90, 28.0, False),    # TR (IRIS)
    ("He II 304",    303.78, 4.90,  4.0, False),    # TR / upper chromo (AIA 304)
    ("C II 1334",   1334.53, 4.40, 12.0, True),     # upper chromosphere (IRIS)
    ("Mg II k 2796", 2796.35, 4.00, 24.0, True),    # chromosphere (IRIS)
    ("H-alpha 6563", 6562.80, 4.00,  1.0, True),    # chromosphere (NLTE)
]

# Lines to render as synthetic (velocity, time) spectrograms: one hot (evaporation
# blueshift) + one cool (condensation redshift) to show the opposed drift.
SPECTROGRAM_LINES = ["Fe XXI 1354", "Si IV 1394"]


# ---------------------------------------------------------------------------
# Parsing (mirrors animate_loop_flare.py)
# ---------------------------------------------------------------------------
def read_meta(path):
    meta, in_meta = {}, False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("[META]"):
            in_meta = True; continue
        if s.startswith("[") and not s.startswith("[META]"):
            in_meta = False
        if in_meta and "=" in s and not s.startswith("#"):
            k, v = s.split("=", 1)
            meta[k.strip()] = v.strip()
    return meta


def parse_dat(path):
    """Return per-cell ds[m], arc-length center[m], phi_g[J/kg], height[m]."""
    rows, in_cells = [], False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("[CELLS]"):
            in_cells = True; continue
        if s.startswith("[GHOSTS]"):
            break
        if not in_cells or not s or s.startswith("#"):
            continue
        rows.append([float(x) for x in s.split()])
    a = np.array(rows)
    ds = a[:, 1]
    phi_c = 0.5 * (a[:, 4] + a[:, 5])
    s_arc = np.cumsum(ds) - 0.5 * ds          # cell-center arc length [m]
    return ds, s_arc, phi_c, phi_c / g_si


def read_frames(path):
    raw = [l for l in open(path).read().splitlines() if l.strip()]
    ns, neq = map(int, raw[0].split())
    i, frames = 2, []
    while i < len(raw):
        t = float(raw[i].split()[3]); i += 1
        d = np.zeros((ns, neq)); ok = True
        for j in range(ns):
            tk = raw[i + j].split()
            if len(tk) != neq:
                ok = False; break
            d[j] = [float(x) for x in tk]
        i += ns
        if ok:
            frames.append((t, d))
    return frames


def primitives(xn, phi_g):
    """Return n_e (=n_i), n_n, V [m/s], T_i, T_e from the conserved state.

    col4 = E_I is the TOTAL charged energy; col6 = E_E the electron internal
    energy (7-var T_e!=T_i runs). 6-var runs: T_e == T_i."""
    ni = xn[:, 0] / m_i
    nn = xn[:, 1] / m_n
    v = xn[:, 2] / xn[:, 0]
    p_i = 2/3*xn[:, 4] - 1/3*xn[:, 0]*v*v - 2/3*xn[:, 0]*phi_g
    if xn.shape[1] >= 7:
        p_e = 2.0/3.0 * xn[:, 6]
        p_p = np.maximum(p_i - p_e, 1e-30)
        Te = p_e / (ni * k_b)
        Ti = p_p / (ni * k_b)
    else:
        Ti = p_i / (2.0 * ni * k_b)
        Te = Ti
    return ni, nn, v, Ti, Te


# ---------------------------------------------------------------------------
# Forward model
# ---------------------------------------------------------------------------
def line_emissivity(ne, Te, logT_peak):
    """eps(s) = n_e^2 * G(T_e), G a Gaussian in log10 T_e (peak-normalized)."""
    logTe = np.log10(np.maximum(Te, 1.0))
    G = np.exp(-0.5 * ((logTe - logT_peak) / GLOGT_SIGMA) ** 2)
    return ne * ne * G


def synth_timeseries(frames, phi_g, ds, dhds):
    """For each frame and line: intensity-weighted Doppler velocity [km/s],
    integrated intensity [arb], and emission-weighted formation height [Mm]."""
    times = np.array([t for t, _ in frames])
    nL = len(LINES)
    vDopp = np.full((len(frames), nL), np.nan)
    Iline = np.zeros((len(frames), nL))
    hform = np.full((len(frames), nL), np.nan)
    for fi, (_, xn) in enumerate(frames):
        ne, nn, v, Ti, Te = primitives(xn, phi_g)
        v_z = v * dhds                              # vertical velocity [m/s], +up
        height = phi_g / g_si
        for li, (_, _lam, logTp, _m, _thk) in enumerate(LINES):
            eps = line_emissivity(ne, Te, logTp)
            w = eps * ds
            W = w.sum()
            Iline[fi, li] = W
            if W > 0:
                # redshift (downflow, away from disk-center observer) positive:
                vDopp[fi, li] = -np.sum(w * v_z) / W / 1.0e3   # [km/s]
                hform[fi, li] = np.sum(w * height) / W / 1.0e6  # [Mm]
    return times, vDopp, Iline, hform


def synth_spectrogram(frames, phi_g, ds, dhds, label):
    """Synthetic (velocity, time) spectrogram for one line: each cell emits a
    Gaussian centered on its v_z with thermal+turbulent width, amplitude eps*ds.
    Returns (v_axis[km/s], times[s], image[nv, nt] normalized per column)."""
    idx = [i for i, L in enumerate(LINES) if L[0] == label][0]
    _, _lam, logTp, mline, _thk = LINES[idx]
    Tpk = 10.0 ** logTp
    w_th = np.sqrt(2.0 * k_b * Tpk / (mline * amu)) / 1.0e3      # [km/s]
    w_tot = np.hypot(w_th, VTURB_KMS)
    v_axis = np.linspace(-VWIN_KMS, VWIN_KMS, NV)
    times = np.array([t for t, _ in frames])
    img = np.zeros((NV, len(frames)))
    for fi, (_, xn) in enumerate(frames):
        ne, nn, v, Ti, Te = primitives(xn, phi_g)
        v_z = v * dhds
        eps = line_emissivity(ne, Te, logTp) * ds
        # redshift positive -> spectral velocity is -v_z (downflow = +)
        vc = -v_z / 1.0e3
        prof = np.zeros(NV)
        for k in range(len(eps)):
            if eps[k] <= 0:
                continue
            prof += eps[k] * np.exp(-0.5 * ((v_axis - vc[k]) / w_tot) ** 2)
        if prof.max() > 0:
            prof /= prof.max()
        img[:, fi] = prof
    return v_axis, times, img


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def _temp_color(logTp):
    """Cool -> warm by formation temperature (blue=cool chromo, red=hot flare)."""
    frac = np.clip((logTp - 4.0) / (7.05 - 4.0), 0.0, 1.0)
    return plt.cm.coolwarm(frac)


def plot_timeseries(times, vDopp, Iline, stem, out_png):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for li, (lab, _lam, logTp, _m, thk) in enumerate(LINES):
        c = _temp_color(logTp)
        tag = lab + (" *" if thk else "")
        ax1.plot(times, vDopp[:, li], color=c, lw=1.8, label=tag)
        I = Iline[:, li]
        Imax = I.max()
        if Imax > 0:
            ax2.plot(times, I / Imax, color=c, lw=1.8, label=tag)
    for ax in (ax1, ax2):
        ax.axvspan(T_ON, T_OFF, color="orange", alpha=0.12, lw=0)
        ax.grid(alpha=0.25)
    ax1.axhline(0.0, color="k", lw=0.8)
    ax1.set_ylabel("Doppler velocity  [km/s]\n(blue<0 up  /  red>0 down)")
    ax1.set_title(f"Synthetic IRIS Doppler shifts — {stem}\n"
                  "evaporation: hot lines blueshift, condensation: cool lines redshift "
                  "(* = optically thick, EM-proxy only)")
    ax1.legend(ncol=2, fontsize=8, loc="upper right")
    ax2.set_ylabel("intensity  [peak-normalized]")
    ax2.set_xlabel("time  [s]")
    ax2.set_yscale("log")
    ax2.set_ylim(1e-4, 1.5)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_spectrograms(frames, phi_g, ds, dhds, stem, out_png):
    labels = [l for l in SPECTROGRAM_LINES
              if l in [L[0] for L in LINES]]
    fig, axes = plt.subplots(1, len(labels), figsize=(6 * len(labels), 5),
                             squeeze=False)
    for ax, lab in zip(axes[0], labels):
        v_axis, times, img = synth_spectrogram(frames, phi_g, ds, dhds, lab)
        pcm = ax.pcolormesh(times, v_axis, img, cmap="inferno",
                            shading="auto", vmin=0, vmax=1)
        ax.axhline(0.0, color="cyan", lw=0.8, alpha=0.7)
        ax.axvspan(T_ON, T_OFF, color="white", alpha=0.10, lw=0)
        ax.set_title(f"{lab}  synthetic spectrogram")
        ax.set_xlabel("time  [s]")
        ax.set_ylabel("Doppler velocity  [km/s]  (red>0 down)")
        fig.colorbar(pcm, ax=ax, label="intensity (col-norm)")
    fig.suptitle(f"Synthetic (v, t) spectrograms — {stem}", y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)


def write_csv(times, vDopp, Iline, out_csv):
    cols = ["time_s"]
    for lab, _lam, _lT, _m, _t in LINES:
        key = lab.replace(" ", "_")
        cols += [f"vDopp_kms[{key}]", f"I[{key}]"]
    with open(out_csv, "w") as f:
        f.write(",".join(cols) + "\n")
        for fi, t in enumerate(times):
            row = [f"{t:.6e}"]
            for li in range(len(LINES)):
                row += [f"{vDopp[fi, li]:.4e}", f"{Iline[fi, li]:.4e}"]
            f.write(",".join(row) + "\n")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    dat_path, txt_path = sys.argv[1], sys.argv[2]
    stem = sys.argv[3] if len(sys.argv) > 3 else \
        os.path.splitext(os.path.basename(txt_path))[0]

    meta = read_meta(dat_path)
    ds, s_arc, phi_g, height = parse_dat(dat_path)
    dhds = np.gradient(height, s_arc)           # vertical component of tangent
    frames = read_frames(txt_path)
    print(f"[synth_doppler] {stem}: topology={meta.get('topology','?')}, "
          f"ns={len(ds)}, frames={len(frames)}, "
          f"t=[{frames[0][0]:.2f}, {frames[-1][0]:.2f}] s")

    times, vDopp, Iline, hform = synth_timeseries(frames, phi_g, ds, dhds)

    vis = "visualization/events"
    out = "outputs/events"
    os.makedirs(vis, exist_ok=True)
    os.makedirs(out, exist_ok=True)
    ts_png = os.path.join(vis, f"{stem}_doppler.png")
    sg_png = os.path.join(vis, f"{stem}_dopplergram.png")
    csv    = os.path.join(out, f"{stem}_doppler.csv")

    plot_timeseries(times, vDopp, Iline, stem, ts_png)
    plot_spectrograms(frames, phi_g, ds, dhds, stem, sg_png)
    write_csv(times, vDopp, Iline, csv)

    # Console summary: peak blueshift of the hot line, peak redshift of cool lines.
    print("[synth_doppler]  line            peak v_Dopp [km/s]   at t [s]")
    for li, (lab, _lam, _lT, _m, _thk) in enumerate(LINES):
        col = vDopp[:, li]
        if np.all(np.isnan(col)):
            continue
        ext_i = np.nanargmax(np.abs(col))
        print(f"               {lab:14s}  {col[ext_i]:+10.1f}        {times[ext_i]:6.2f}")
    print(f"[synth_doppler] wrote {ts_png}\n"
          f"                      {sg_png}\n"
          f"                      {csv}")


if __name__ == "__main__":
    main()
