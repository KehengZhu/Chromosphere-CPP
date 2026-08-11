#!/usr/bin/env python3
"""Space-time loader/analyser for model_column `.gamma_diag` sidecars.

The sidecar holds the EOS-decoded cell-centred state written by `chromo_main`
(columns rho_total v_cm T x_eq n_e n_HI p_total Gamma1 kappa_physical
kappa_solver), one block per output frame.  This module turns a run into dense
(n_time, n_cell) arrays, caches them as .npz next to the source, and offers the
small set of space-time primitives used by the low-chromosphere sloshing
diagnosis: height slicing, time series at a height, and cross-correlation lag
between two heights.

Usage:
    python util/column_spacetime.py <run>.gamma_diag [--cache-only]
"""
from __future__ import annotations

import os
import sys

import numpy as np

COLS = ("rho v T x_eq n_e n_HI p Gamma1 kappa_phys kappa_solver").split()


def load(path, use_cache=True):
    """-> (h_km[ns], t[nt], data[nt, ns, 10])."""
    cache = path + ".npz"
    if use_cache and os.path.exists(cache) and \
            os.path.getmtime(cache) >= os.path.getmtime(path):
        z = np.load(cache)
        return z["h"], z["t"], z["d"]

    with open(path) as fh:
        ns, ncol = map(int, fh.readline().split())
        h = np.fromstring(fh.readline(), sep=" ")
        assert h.size == ns, (h.size, ns)
        times, blocks, cur = [], [], []
        for line in fh:
            if line[0] == "#":
                if line.startswith("# t "):
                    if cur:
                        blocks.append(cur)
                        cur = []
                    times.append(float(line.split("t =")[1].split("step")[0]))
                continue
            cur.append(line)
        if cur:
            blocks.append(cur)
    nt = len(blocks)
    d = np.empty((nt, ns, ncol))
    for k, blk in enumerate(blocks):
        d[k] = np.fromstring(" ".join(blk), sep=" ").reshape(ns, ncol)
    t = np.asarray(times[:nt])
    np.savez_compressed(cache, h=h, t=t, d=d)
    return h, t, d


def col(d, name):
    return d[..., COLS.index(name)]


def at_height(h, d, name, hkm):
    """Time series of `name` at the cell nearest `hkm`."""
    i = int(np.argmin(np.abs(h - hkm)))
    return i, col(d, name)[:, i]


def band(h, d, name, lo, hi):
    m = (h >= lo) & (h <= hi)
    return h[m], col(d, name)[:, m]


def xcorr_lag(t, a, b, max_lag=None):
    """Lag (in units of t) maximising corr(a(t), b(t+lag)); + => b lags a."""
    a = a - a.mean()
    b = b - b.mean()
    n = len(a)
    lags = np.arange(-(n - 1), n)
    c = np.correlate(b, a, mode="full")
    norm = np.sqrt((a**2).sum() * (b**2).sum())
    c = c / max(norm, 1e-300)
    if max_lag is not None:
        dt = t[1] - t[0]
        keep = np.abs(lags) <= int(round(max_lag / dt))
        lags, c = lags[keep], c[keep]
    k = int(np.argmax(c))
    dt = t[1] - t[0]
    return lags[k] * dt, float(c[k])


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    for p in sys.argv[1:]:
        if p.startswith("--"):
            continue
        h, t, d = load(p)
        print(f"{p}: ns={h.size} nt={t.size} h=[{h[0]:.2f},{h[-1]:.2f}] km "
              f"t=[{t[0]:.1f},{t[-1]:.1f}] s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
