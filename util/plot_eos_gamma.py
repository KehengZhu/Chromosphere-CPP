#!/usr/bin/env python3
"""Validate and plot the classical pure-H CRASH table (plan Stages 0/1).

The validation is intentionally fail-closed.  It checks the raw table over the
whole rectangular grid, compares CRASH zAv/Gamma/Cv/GammaS with the analytic
Saha closure, verifies strictly increasing internal energy and positive Gamma1,
and confirms that the stripped three-column runtime table is identical to the
validated GammaS column.

Without --validate-only it writes, under visualization/eos_gamma/:
  gamma_hydrogen.png       Gamma_E, Gamma_1, and ionization at fixed densities
  gamma_hydrogen_2d.png    Gamma_1(T,n_H) with ionization contours
  c7_eos_consistency.png   C7 x, gammas, and heat capacity vs height
"""

import argparse
import os
import re

import numpy as np


RAW_TABLE = "outputs/eos_gamma/gamma_hydrogen.dat"
RUNTIME_TABLE = "data/eos/gamma1_hydrogen_v1.dat"
C7_CPP = "scenarios/model_c7.cpp"
OUTDIR = "visualization/eos_gamma"

# Raw CRASH table columns.
COL = dict(T=0, nH=1, rho=2, x=3, p=4, e=5, gamma_e=6,
           gamma1=7, gamma_electron=8, gamma1_electron=9, cv=10)

K_B = 1.380649e-23
M_E = 9.1093837015e-31
M_H = 1.6726219e-27
H_PLANCK = 6.62607015e-34
CHI_H = 2.179872361e-18
SAHA_LOG_C = 1.5 * np.log(2.0 * np.pi * M_E * K_B / H_PLANCK**2)
GAMMA_IDEAL = 5.0 / 3.0
NA_LINES = [1e14, 1e18, 1e22, 1e26]


def saha_ionization_fraction_nh(n_h, temperature):
    """Stable textbook pure-H Saha fraction for arrays/scalars.

    x^2/(1-x) = (2*pi*m_e*k*T/h^2)^(3/2) exp(-chi/kT) / n_H.
    Hydrogen's ground degeneracy cancels the free-electron spin factor.
    """
    n_h, temperature = np.broadcast_arrays(
        np.asarray(n_h, dtype=float), np.asarray(temperature, dtype=float))
    if np.any(~np.isfinite(n_h)) or np.any(n_h <= 0.0):
        raise ValueError("n_H must be positive and finite")
    if np.any(~np.isfinite(temperature)) or np.any(temperature <= 0.0):
        raise ValueError("temperature must be positive and finite")
    log_a = (SAHA_LOG_C + 1.5 * np.log(temperature)
             - CHI_H / (K_B * temperature) - np.log(n_h))
    out = np.empty_like(log_a)
    high = log_a >= 0.0
    out[high] = 2.0 / (1.0 + np.sqrt(1.0 + 4.0 * np.exp(-log_a[high])))
    u = np.exp(0.5 * log_a[~high])
    out[~high] = 2.0 * u / (u + np.sqrt(u * u + 4.0))
    return out.item() if out.ndim == 0 else out


def analytic_thermo(n_h, temperature):
    """Return x, e, gamma_E, Cv/(n_H k_B), and analytic Gamma1."""
    n_h, temperature = np.broadcast_arrays(
        np.asarray(n_h, dtype=float), np.asarray(temperature, dtype=float))
    x = saha_ionization_fraction_nh(n_h, temperature)
    theta = CHI_H / (K_B * temperature)
    dx_dln_t = x * (1.0 - x) * (1.5 + theta) / (2.0 - x)
    dx_dln_n = -x * (1.0 - x) / (2.0 - x)
    pressure = (1.0 + x) * n_h * K_B * temperature
    energy = 1.5 * pressure + x * n_h * CHI_H
    gamma_e = 1.0 + pressure / energy
    cv_norm = 1.5 * (1.0 + x) + (1.5 + theta) * dx_dln_t
    chi_rho = 1.0 + dx_dln_n / (1.0 + x)
    d_p_over_nk_d_t = 1.0 + x + dx_dln_t
    gamma1 = chi_rho + d_p_over_nk_d_t**2 / (cv_norm * (1.0 + x))
    return dict(x=x, p=pressure, e=energy, gamma_e=gamma_e,
                cv=cv_norm, gamma1=gamma1)


def cv_effective_central(n_h, temperature, delta_ln_t=1e-4):
    """C_V/(n_H k_B) from the required central difference in ln(T)."""
    e_plus = analytic_thermo(n_h, temperature * np.exp(delta_ln_t))["e"]
    e_minus = analytic_thermo(n_h, temperature * np.exp(-delta_ln_t))["e"]
    return (e_plus - e_minus) / (2.0 * delta_ln_t * temperature * n_h * K_B)


def load_raw_grid(path):
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] != 11:
        raise RuntimeError(f"{path}: expected 11 columns, got shape {data.shape}")
    if not np.all(np.isfinite(data)):
        raise RuntimeError(f"{path}: contains non-finite values")
    t_vals = np.unique(data[:, COL["T"]])
    n_vals = np.unique(data[:, COL["nH"]])
    if len(t_vals) < 2 or len(n_vals) < 2 or len(data) != len(t_vals) * len(n_vals):
        raise RuntimeError(f"{path}: table is not a complete rectangular grid")
    if not np.all(np.diff(t_vals) > 0.0) or not np.all(np.diff(n_vals) > 0.0):
        raise RuntimeError(f"{path}: axes must be strictly increasing")
    expected_t = np.tile(t_vals, len(n_vals))
    expected_n = np.repeat(n_vals, len(t_vals))
    if not (np.array_equal(data[:, COL["T"]], expected_t)
            and np.array_equal(data[:, COL["nH"]], expected_n)):
        raise RuntimeError(f"{path}: expected density-major, temperature-minor order")
    grids = {name: data[:, index].reshape(len(n_vals), len(t_vals))
             for name, index in COL.items()}
    return t_vals, n_vals, grids


def _max_rel(actual, expected, floor=1e-300):
    return float(np.max(np.abs(actual - expected) / np.maximum(np.abs(expected), floor)))


def validate_tables(raw_path=RAW_TABLE, runtime_path=RUNTIME_TABLE):
    t_vals, n_vals, grid = load_raw_grid(raw_path)
    tt, nn = np.meshgrid(t_vals, n_vals)
    analytic = analytic_thermo(nn, tt)
    cv_central = cv_effective_central(nn, tt)

    failures = []
    if np.any(np.diff(grid["e"], axis=1) <= 0.0):
        failures.append("CRASH internal energy is not strictly increasing in T")
    if np.any(cv_central <= 0.0):
        failures.append("analytic Cv_eff is not positive")
    if np.any(grid["cv"] <= 0.0):
        failures.append("CRASH Cv is not positive")
    if np.any(grid["gamma1"] <= 0.0):
        failures.append("CRASH Gamma1 is not positive")

    transition = (analytic["x"] >= 1e-10) & (analytic["x"] <= 1.0 - 1e-10)
    x_abs = float(np.max(np.abs(grid["x"] - analytic["x"])))
    x_rel = _max_rel(grid["x"][transition], analytic["x"][transition])
    tiny = (analytic["x"] < 1e-8) & (grid["x"] > 0.0)
    x_log = (float(np.max(np.abs(np.log(grid["x"][tiny])
                                 - np.log(analytic["x"][tiny]))))
             if np.any(tiny) else 0.0)
    e_rel = _max_rel(grid["e"], analytic["e"])
    gamma_e_abs = float(np.max(np.abs(grid["gamma_e"] - analytic["gamma_e"])))
    cv_rel = _max_rel(grid["cv"], cv_central)
    gamma1_abs = float(np.max(np.abs(grid["gamma1"] - analytic["gamma1"])))

    # These tolerances are tight enough to catch a factor-two statistical
    # weight, wrong ionization potential, hidden CRASH floor, or column swap,
    # while allowing CRASH's nonlinear solve and table-format rounding.
    checks = [
        (x_abs <= 1e-6, f"Saha zAv max abs error {x_abs:.3e} > 1e-6"),
        (x_rel <= 2e-3, f"Saha zAv transition max rel error {x_rel:.3e} > 2e-3"),
        (x_log <= 2e-3, f"Saha log(x) max error {x_log:.3e} > 2e-3"),
        (e_rel <= 2e-6, f"analytic energy max rel error {e_rel:.3e} > 2e-6"),
        (gamma_e_abs <= 3e-5,
         f"analytic gamma_E max abs error {gamma_e_abs:.3e} > 3e-5"),
        (cv_rel <= 1e-3, f"central-difference Cv max rel error {cv_rel:.3e} > 1e-3"),
        (gamma1_abs <= 7e-4,
         f"analytic Gamma1 max abs error {gamma1_abs:.3e} > 7e-4"),
    ]
    failures.extend(message for passed, message in checks if not passed)

    runtime = np.loadtxt(runtime_path)
    if runtime.ndim != 2 or runtime.shape != (len(t_vals) * len(n_vals), 3):
        failures.append(f"runtime table has wrong shape {runtime.shape}")
    elif not np.all(np.isfinite(runtime)):
        failures.append("runtime table contains non-finite values")
    else:
        expected = np.column_stack((np.log10(grid["T"].ravel()),
                                    np.log10(grid["nH"].ravel()),
                                    grid["gamma1"].ravel()))
        runtime_err = float(np.max(np.abs(runtime - expected)))
        if runtime_err > 5e-14:
            failures.append(f"runtime table differs from validated Gamma1 ({runtime_err:.3e})")

    metrics = dict(nT=len(t_vals), nN=len(n_vals), Tmin=t_vals[0], Tmax=t_vals[-1],
                   nmin=n_vals[0], nmax=n_vals[-1], x_abs=x_abs, x_rel=x_rel,
                   x_log=x_log, e_rel=e_rel, gamma_e_abs=gamma_e_abs,
                   cv_rel=cv_rel, gamma1_abs=gamma1_abs,
                   min_cv=float(min(np.min(grid["cv"]), np.min(cv_central))),
                   min_gamma1=float(np.min(grid["gamma1"])))
    print("Stage 0/1 whole-grid validation:")
    print("  grid: %d x %d, T=[%.6g, %.6g] K, n_H=[%.6g, %.6g] m^-3" %
          (metrics["nN"], metrics["nT"], metrics["Tmin"], metrics["Tmax"],
           metrics["nmin"], metrics["nmax"]))
    print("  Saha: max_abs=%.3e transition_rel=%.3e tiny_log=%.3e" %
          (x_abs, x_rel, x_log))
    print("  analytic: energy_rel=%.3e gammaE_abs=%.3e Cv_rel=%.3e Gamma1_abs=%.3e" %
          (e_rel, gamma_e_abs, cv_rel, gamma1_abs))
    print("  gates: min Cv/(n_H k_B)=%.6g, min Gamma1=%.6g" %
          (metrics["min_cv"], metrics["min_gamma1"]))
    if failures:
        raise RuntimeError("EOS validation failed:\n  - " + "\n  - ".join(failures))
    print("  PASS: raw table, analytic closure, and runtime table are consistent")
    return t_vals, n_vals, grid, metrics


def _nearest_density(n_vals, target):
    return int(np.argmin(np.abs(np.log10(n_vals) - np.log10(target))))


def plot_lines(t_vals, n_vals, grid, path):
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8.4, 8.0), sharex=True,
                                   gridspec_kw=dict(height_ratios=[2, 1], hspace=0.08))
    cmap = matplotlib.colormaps["viridis"]
    norm = colors.Normalize(np.log10(n_vals[0]), np.log10(n_vals[-1]))
    for target in NA_LINES:
        i = _nearest_density(n_vals, target)
        color = cmap(norm(np.log10(n_vals[i])))
        label = r"$n_H=10^{%.1f}\,\mathrm{m}^{-3}$" % np.log10(n_vals[i])
        ax0.plot(t_vals, grid["gamma_e"][i], color=color, lw=2, label=label)
        ax0.plot(t_vals, grid["gamma1"][i], color=color, lw=1.4, ls="--")
        ax1.plot(t_vals, grid["x"][i], color=color, lw=2)
    ax0.axhline(GAMMA_IDEAL, color="0.4", ls=":", lw=1.2)
    ax0.set_ylabel(r"thermodynamic index")
    ax0.set_ylim(1.0, 1.70)
    ax0.legend(fontsize=8.5, ncol=2)
    ax0.set_title(r"Classical pure-H Saha closure (solid $γ_E$, dashed $Γ_1$)")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"temperature $T$ [K]")
    ax1.set_ylabel(r"ionization $x$")
    ax1.set_ylim(-0.03, 1.03)
    for ax in (ax0, ax1):
        ax.grid(True, which="both", ls=":", alpha=0.35)
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(t_vals, n_vals, grid, path):
    tt, nn = np.meshgrid(t_vals, n_vals)
    fig, ax = plt.subplots(figsize=(8.5, 6.3))
    pcm = ax.pcolormesh(tt, nn, grid["gamma1"], shading="auto", cmap="magma_r",
                        vmin=1.0, vmax=GAMMA_IDEAL)
    fig.colorbar(pcm, ax=ax, label=r"sound-speed index $Γ_1$")
    cs = ax.contour(tt, nn, grid["x"], levels=[0.01, 0.1, 0.5, 0.9, 0.99],
                    colors="cyan", linewidths=0.8)
    ax.clabel(cs, fmt=lambda value: f"x={value:g}", fontsize=8)
    ax.set(xscale="log", yscale="log", xlabel=r"temperature $T$ [K]",
           ylabel=r"hydrogen nuclei density $n_H$ [m$^{-3}$]",
           title=r"CRASH $Γ_1(T,n_H)$, classical pure hydrogen")
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def parse_c7_tables(path=C7_CPP):
    with open(path, encoding="utf-8") as stream:
        source = stream.read()
    def grab(name):
        match = re.search(name + r"\[[^\]]*\]\[4\]\s*=\s*\{(.*?)\};", source, re.S)
        if match is None:
            raise RuntimeError(f"cannot find {name} in {path}")
        rows = re.findall(r"\{([^}]*)\}", match.group(1))
        return np.array([[float(item.strip().rstrip("f")) for item in row.split(",")]
                         for row in rows])
    return grab("MODEL_C7_PHOTO"), grab("MODEL_C7")


def sample_c7(height_km, photo, c7):
    table = np.vstack((photo, c7))
    # One duplicate coronal height exists in the source table; retain its first row.
    _, unique = np.unique(table[:, 0], return_index=True)
    table = table[np.sort(unique)]
    height, ne, nhi, temperature = table.T
    return (np.interp(height_km, height, temperature),
            np.exp(np.interp(height_km, height, np.log(ne))),
            np.exp(np.interp(height_km, height, np.log(nhi))))


def interp_log_grid(t_vals, n_vals, values, temperature, n_h):
    lt, ln = np.log10(t_vals), np.log10(n_vals)
    x, y = np.log10(temperature), np.log10(n_h)
    if np.any((x < lt[0]) | (x > lt[-1]) | (y < ln[0]) | (y > ln[-1])):
        raise RuntimeError("C7 profile lies outside the EOS table")
    jt = np.clip(np.searchsorted(lt, x) - 1, 0, len(lt) - 2)
    jn = np.clip(np.searchsorted(ln, y) - 1, 0, len(ln) - 2)
    ft = (x - lt[jt]) / (lt[jt + 1] - lt[jt])
    fn = (y - ln[jn]) / (ln[jn + 1] - ln[jn])
    return ((1-fn) * ((1-ft)*values[jn, jt] + ft*values[jn, jt+1])
            + fn * ((1-ft)*values[jn+1, jt] + ft*values[jn+1, jt+1]))


def _split_height_scale(ax, split=3000.0, compression=0.04):
    def forward(value):
        value = np.asarray(value)
        return np.where(value <= split, value, split + compression * (value - split))
    def inverse(value):
        value = np.asarray(value)
        return np.where(value <= split, value, split + (value - split) / compression)
    ax.set_xscale("function", functions=(forward, inverse))


def plot_c7(t_vals, n_vals, grid, path):
    photo, c7 = parse_c7_tables()
    top = float(c7[-1, 0])
    height = np.unique(np.r_[np.linspace(0.0, 3000.0, 1200),
                             np.linspace(3000.0, top, 800)])
    temperature, ne, nhi = sample_c7(height, photo, c7)
    n_h = ne + nhi
    x_c7 = ne / n_h
    analytic = analytic_thermo(n_h, temperature)
    gamma_e = interp_log_grid(t_vals, n_vals, grid["gamma_e"], temperature, n_h)
    gamma1 = interp_log_grid(t_vals, n_vals, grid["gamma1"], temperature, n_h)
    cv_crash = interp_log_grid(t_vals, n_vals, grid["cv"], temperature, n_h)

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.2), sharex=True,
                             gridspec_kw=dict(hspace=0.08))
    axes[0].plot(height, np.clip(x_c7, 1e-9, 1-1e-9), lw=2, color="#1f77b4",
                 label=r"C7 $x=n_e/(n_e+n_{HI})$")
    axes[0].plot(height, np.clip(analytic["x"], 1e-9, 1-1e-9), lw=1.8,
                 color="#d62728", ls="--", label=r"textbook Saha $x(ρ,T)$")
    axes[0].set_yscale("logit")
    axes[0].set_ylabel(r"ionization $x$")
    axes[0].legend(fontsize=8.5, loc="best")

    axes[1].plot(height, gamma_e, lw=2, label=r"CRASH $γ_E$")
    axes[1].plot(height, gamma1, lw=2, label=r"CRASH $Γ_1$")
    axes[1].plot(height, analytic["gamma_e"], lw=1.2, ls=":",
                 color="black", label=r"analytic $γ_E$")
    axes[1].set_ylabel("index")
    axes[1].set_ylim(1.0, 1.70)
    axes[1].legend(fontsize=8.5, ncol=3)

    cv_effective = cv_effective_central(n_h, temperature)
    axes[2].plot(height, cv_effective, lw=2, color="#2ca02c",
                 label=r"$C_V^{eff}=∂e_{int}/∂T$")
    axes[2].plot(height, cv_crash, lw=1.5, ls="--", color="#9467bd",
                 label=r"CRASH $C_V$")
    axes[2].set_ylabel(r"$C_V/(n_H k_B)$")
    axes[2].set_xlabel("height above photosphere [km]")
    axes[2].legend(fontsize=8.5)

    for ax in axes:
        _split_height_scale(ax)
        ax.axvline(2153.0, color="0.35", ls="--", lw=0.9, alpha=0.6)
        ax.axvline(2628.0, color="0.35", ls="--", lw=0.9, alpha=0.6)
        ax.grid(True, which="both", ls=":", alpha=0.3)
    axes[0].set_title("Model C7 track through the validated classical-Saha EOS\n"
                      "(chromosphere/TR expanded; coronal height compressed)")
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate-only", action="store_true",
                        help="run strict checks without generating figures")
    args = parser.parse_args()
    t_vals, n_vals, grid, _ = validate_tables()
    if args.validate_only:
        return
    global matplotlib, plt, colors
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    os.makedirs(OUTDIR, exist_ok=True)
    outputs = [
        (plot_lines, os.path.join(OUTDIR, "gamma_hydrogen.png")),
        (plot_heatmap, os.path.join(OUTDIR, "gamma_hydrogen_2d.png")),
        (plot_c7, os.path.join(OUTDIR, "c7_eos_consistency.png")),
    ]
    for function, path in outputs:
        function(t_vals, n_vals, grid, path)
        print("wrote", path)


if __name__ == "__main__":
    main()
