#!/usr/bin/env python3
"""Regression tests for the shared physical conductive-flux diagnostic."""

import pathlib
import sys
import unittest

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "util"))
from conductive_flux import (  # noqa: E402
    conductivity_components,
    physical_conductivity,
    physical_heat_flux,
)


class ConductiveFluxTest(unittest.TestCase):
    def setUp(self):
        self.h = np.array([0.0, 25.0, 70.0, 140.0])
        self.temperature = np.array([6000.0, 6500.0, 9000.0, 20000.0])
        self.n_e = np.array([1e16, 2e16, 4e16, 8e16])
        self.n_hi = np.array([1e20, 8e19, 4e19, 1e19])

    def test_components_match_physics_hpp_formulas(self):
        ke, kn = conductivity_components(self.n_e, self.n_hi, self.temperature)
        expected_e = (9.2048e-12 * self.n_e * self.temperature**2.5 /
                      (self.n_e + 2.836e-11*self.n_hi*self.temperature**2))
        expected_n = (0.0342006*self.n_hi*self.temperature /
                      (1.20613*self.n_e*np.sqrt(2*self.temperature)
                       + 1.70573*self.n_hi*np.sqrt(self.temperature)))
        np.testing.assert_array_equal(ke, expected_e)
        np.testing.assert_array_equal(kn, expected_n)
        np.testing.assert_array_equal(physical_conductivity(
            self.n_e, self.n_hi, self.temperature), expected_e + expected_n)

    def test_identical_fixed_and_gamma_fields_give_identical_flux(self):
        fixed_path = physical_heat_flux(
            self.h, self.temperature, self.n_e, self.n_hi)
        gamma_path = physical_heat_flux(
            self.h.copy(), self.temperature.copy(), self.n_e.copy(), self.n_hi.copy())
        np.testing.assert_array_equal(fixed_path, gamma_path)


if __name__ == "__main__":
    unittest.main()
