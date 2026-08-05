"""Shared physical field-aligned conductivity diagnostics.

The formulas match ``physics.hpp`` exactly for common-temperature pure hydrogen.
They exclude TRAC broadening, numerical diffusivity, imposed boundary fluxes, and
all other solver-only terms.
"""

import numpy as np


def conductivity_components(n_e, n_hi, temperature):
    """Return unbroadened ``(kappa_e, kappa_n)`` in W m^-1 K^-1."""
    n_e = np.asarray(n_e, dtype=float)
    n_hi = np.asarray(n_hi, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    kappa_e = (9.2048e-12 * n_e * temperature**2.5
               / (n_e + 2.836e-11 * n_hi * temperature**2))
    kappa_n = (0.0342006 * n_hi * temperature
               / (1.20613 * n_e * np.sqrt(2.0 * temperature)
                  + 1.70573 * n_hi * np.sqrt(temperature)))
    return kappa_e, kappa_n


def physical_conductivity(n_e, n_hi, temperature):
    """Unbroadened electron-plus-neutral conductivity in W m^-1 K^-1."""
    kappa_e, kappa_n = conductivity_components(n_e, n_hi, temperature)
    return kappa_e + kappa_n


def physical_heat_flux(h_km, temperature, n_e, n_hi):
    """Cell-centred q_phys = -(kappa_e+kappa_n) dT/ds in W m^-2.

    This common diagnostic uses ``numpy.gradient`` on cell centres. It is not the
    nonlinear solver's face flux and deliberately excludes solver-only terms.
    """
    s_m = np.asarray(h_km, dtype=float) * 1.0e3
    dtds = np.gradient(np.asarray(temperature, dtype=float), s_m)
    return -physical_conductivity(n_e, n_hi, temperature) * dtds
