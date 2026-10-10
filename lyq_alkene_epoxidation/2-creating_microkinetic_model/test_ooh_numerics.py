# -*- coding: utf-8 -*-
"""No-CatMAP smoke tests for the independent OOH steady-state quality gate.

Run: python -m unittest -v test_ooh_numerics
"""
import unittest

from run_ooh_analysis import check_steady_state


def gas_fluxes(rate):
    return {
        "C6H12O_g": rate,
        "C6H12_g": -rate,
        "H2O2_g": -rate,
        "H2O_g": rate,
    }


class NumericalGateTests(unittest.TestCase):
    def test_mathematically_consistent_cycle(self):
        rate = 1e-20
        r = check_steady_state(
            [rate] * 8, gas_fluxes(rate), [0.2, 0.8],
            solver_residual="1e-130", numbers_solver=True
        )
        self.assertEqual(r["quality_status"], "validated_numerically")
        self.assertTrue(r["quality_pass"])

    def test_spurious_legacy_desorption_tof_is_rejected(self):
        gas = gas_fluxes(1e-55)
        gas["C6H12O_g"] = 1e-26
        r = check_steady_state(
            [1e-55] * 4 + [1e-26] + [1e-55] * 3,
            gas, [0.2, 0.8], solver_residual="1e-51",
            numbers_solver=True
        )
        self.assertEqual(r["quality_status"], "failed_steady_state_checks")
        self.assertFalse(r["quality_pass"])
        self.assertGreater(r["cycle_max_relative_error"], 0.99)

    def test_tiny_consistent_rate_is_not_used_for_tof_fit(self):
        rate = 1e-60
        r = check_steady_state(
            [rate] * 8, gas_fluxes(rate), [0.2, 0.8],
            solver_residual="1e-140", numbers_solver=True
        )
        self.assertEqual(r["quality_status"], "below_diagnostic_rate_floor")
        self.assertFalse(r["quality_pass"])

    def test_wrong_net_gas_signs_are_rejected(self):
        r = check_steady_state(
            [1e-20] * 8, {
                "C6H12O_g": 1e-20,
                "C6H12_g": 1e-20,
                "H2O2_g": -1e-20,
                "H2O_g": 1e-20
            },
            [0.2, 0.8], solver_residual="1e-130",
            numbers_solver=True
        )
        self.assertFalse(r["quality_pass"])


if __name__ == "__main__":
    unittest.main()
