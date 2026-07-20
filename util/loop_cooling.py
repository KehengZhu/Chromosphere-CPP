#!/usr/bin/env python3
"""Measure the post-beam cooldown of a flare loop run (does the loop cool, how fast?).

The beam is on over [t_on, t_on+dur] (default [2,12] s); afterward the evaporated,
flare-heated loop cools — first by thermal conduction (fast, hot/tenuous), then by
radiation (slow, after evaporation has loaded it dense). This script extracts the
apex / coronal temperature and density vs time, finds the post-beam peak, and
reports e-folding cooling times + point estimates of the conductive and radiative
cooling timescales from the sim's own peak (n, T, L) for comparison with Cargill
et al. (1995): tau_cond = 3 n k L^2 / (kappa0 T^{5/2}), tau_rad = 3 k T / (n Lambda).

Usage:
  python util/loop_cooling.py [loop.dat] [flare.txt] [out.png]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m_i = 1.6726219e-27
k_b = 1.380649e-23
g_si = 274.0
kappa0 = 1.0e-11          # Spitzer kappa = kappa0 T^{5/2} [W m^-1 K^-7/2], SI
T_BEAM_OFF = 12.0         # beam off time [s] (t_on 2 + dur 10)


def read_meta(path):
    meta, in_meta = {}, False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("[META]"):
            in_meta = True; continue
        if s.startswith("[") and not s.startswith("[META]"):
            in_meta = False
        if in_meta and "=" in s and not s.startswith("#"):
            k, v = s.split("=", 1); meta[k.strip()] = v.strip()
    return meta


def parse_dat(path):
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
    return ds, phi_c, phi_c / g_si       # ds[m], phi_g[J/kg], height[m]


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


def lambda_si(T):
    """Optically-thin radiative loss Lambda(T) so loss = n_e n_H Lambda.

    CL2012/RTV-like curve in CGS [erg cm^3 s^-1]: plateau ~10^-21.3 near the 10^5 K
    peak, declining ~T^-0.55 through the 10^6-10^7 K coronal range, flattening into
    the bremsstrahlung regime above 10^7 K. Converted to SI [W m^3] (x 1e-13 = the
    erg cm^3 s^-1 -> W m^3 factor, with n in m^-3)."""
    T = np.atleast_1d(T).astype(float)
    lt = np.log10(np.clip(T, 1e4, 1e8))
    logL_cgs = np.where(
        lt < 5.3, -21.3,
        np.where(lt < 7.0, -21.3 - 0.55 * (lt - 5.3),
                 -22.235 + 0.30 * (lt - 7.0)))         # erg cm^3 s^-1
    return 10.0 ** logL_cgs * 1.0e-13                   # -> W m^3


def main():
    dat = sys.argv[1] if len(sys.argv) > 1 else "scenarios/data/loop_full_fine.dat"
    flare = sys.argv[2] if len(sys.argv) > 2 else "outputs/loops/loop_full_flare.txt"
    out = sys.argv[3] if len(sys.argv) > 3 else "visualization/loops/loop_cooling.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)

    topology = read_meta(dat).get("topology", "closed")
    ds, phi_g, height = parse_dat(dat)
    L_half = float(np.sum(ds)) / (2.0 if topology == "full" else 1.0)   # footpoint->apex [m]
    apex = int(np.argmax(phi_g))
    frames = read_frames(flare)

    ts, T_apex, T_max, T_em, n_apex = [], [], [], [], []
    # corona = cells coronal in the IC
    _, xn0 = frames[0]
    ni0 = xn0[:, 0] / m_i
    v0 = xn0[:, 2] / xn0[:, 0]
    p0 = 2/3*xn0[:, 4] - 1/3*xn0[:, 0]*v0*v0 - 2/3*xn0[:, 0]*phi_g
    Ti0 = p0 / (2 * ni0 * k_b)
    cor = Ti0 > 1.0e5

    for t, xn in frames:
        ni = xn[:, 0] / m_i
        v = xn[:, 2] / xn[:, 0]
        p = 2/3*xn[:, 4] - 1/3*xn[:, 0]*v*v - 2/3*xn[:, 0]*phi_g
        Ti = p / (2 * ni * k_b)
        ts.append(t)
        T_apex.append(Ti[apex]); T_max.append(np.nanmax(Ti))
        w = (ni[cor] ** 2) * ds[cor]
        T_em.append(np.sum(w * Ti[cor]) / np.sum(w))     # EM-weighted coronal T
        n_apex.append(ni[apex])
    ts = np.array(ts); T_apex = np.array(T_apex); T_max = np.array(T_max)
    T_em = np.array(T_em); n_apex = np.array(n_apex)

    # Reference = the loop's hottest-cell temperature (smooth cooldown envelope),
    # searched from a short settle buffer after beam-off — this skips both the
    # impulsive single-cell front spike (~beam-off) and the full loop's dynamic,
    # spiky apex (material sloshes through it), and the EM-weighted T artifact
    # (which mis-peaks at t_off for the full loop as cool-dense evaporate arrives).
    SETTLE = 15.0
    post = ts >= (T_BEAM_OFF + SETTLE)
    ip = np.argmax(T_max * post)
    tp, Tp, np_p = ts[ip], T_max[ip], n_apex[ip]
    after = ts >= tp
    def t_cross(frac):
        target = Tp * frac
        idx = np.where(after & (T_max <= target))[0]
        return ts[idx[0]] - tp if len(idx) else np.nan
    t_e = t_cross(1/np.e); t_half = t_cross(0.5); t_q = t_cross(2.0e6/Tp)

    # point-estimate cooling times from the peak (n, T, L)
    tau_cond = 3.0 * np_p * k_b * L_half**2 / (kappa0 * Tp**2.5)
    tau_rad = 3.0 * k_b * Tp / (np_p * lambda_si(Tp)[0])

    print("=" * 72)
    print(f"{flare}  topology={topology}  L(fp->apex)={L_half/1e6:.1f} Mm  beam off @ {T_BEAM_OFF:.0f}s")
    print(f"settled peak loop T (hottest cell) = {Tp/1e6:6.2f} MK @ t={tp:.0f}s")
    print(f"observed cooldown:  t(->T/e)={t_e:7.0f}s   t(->T/2)={t_half:7.0f}s   "
          f"t(->2MK)={t_q if not np.isnan(t_q) else float('nan'):7.0f}s   (run ends t={ts[-1]:.0f}s)")
    print(f"apex collapse: n_apex {n_apex.max():.1e} -> {n_apex[-1]:.1e} m^-3 (loop drains; "
          f"T_apex end {T_apex[-1]/1e6:.3f} MK = catastrophic/over-drain if no background heating)")
    print(f"point estimates @ peak:  tau_cond={tau_cond:7.0f}s   tau_rad={tau_rad:7.0f}s  "
          f"(Cargill 1995: conduction first, then radiation)")
    print(f"apex T: peak {np.nanmax(T_apex)/1e6:.1f} MK -> end {T_apex[-1]/1e6:.2f} MK")
    print("=" * 72)

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"Flare loop cooldown ({topology}, L={L_half/1e6:.0f} Mm): "
                 f"beam off at t={T_BEAM_OFF:.0f}s", fontsize=12)

    a = ax[0]
    a.semilogy(ts, T_max/1e6, color="0.6", lw=1.2, label="max T (any cell)")
    a.semilogy(ts, T_apex/1e6, color="crimson", lw=1.5, label="apex T")
    a.semilogy(ts, T_em/1e6, color="C0", lw=1.6, label="corona T (EM-weighted)")
    a.axvspan(2, T_BEAM_OFF, color="orange", alpha=0.18, label="beam on")
    a.axvline(tp, color="C0", ls=":", lw=1.0)
    if not np.isnan(t_half):
        a.axvline(tp+t_half, color="green", ls="--", lw=1.0, alpha=0.7,
                  label=f"T/2 at +{t_half:.0f}s")
    a.axhline(0.02, color="none")   # keep collapse (~0.04 MK) on-axis
    a.set_xlabel("t [s]"); a.set_ylabel("T [MK]"); a.grid(alpha=0.3, which="both")
    a.legend(fontsize=8, loc="upper right"); a.set_ylim(0.02, max(60, T_max.max()/1e6*1.2))
    # the cooldown + collapse onset all happen early; don't waste axis on the flat tail
    a.set_xlim(0, min(ts[-1], max(2400.0, 4.0 * (tp + (t_q if not np.isnan(t_q) else 600)))))

    # T vs n trajectory of the apex (Cargill radiative phase ~ T proportional to n^2)
    a = ax[1]
    sc = a.scatter(n_apex, T_apex/1e6, c=ts, cmap="viridis", s=10)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"apex $n_e$ [m$^{-3}$]"); a.set_ylabel("apex T [MK]")
    a.set_title("apex (n, T) trajectory — color = time")
    fig.colorbar(sc, ax=a, label="t [s]"); a.grid(alpha=0.3, which="both")

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
