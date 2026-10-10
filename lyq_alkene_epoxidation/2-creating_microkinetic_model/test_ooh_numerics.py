# -*- coding: utf-8 -*-
"""No-CatMAP smoke tests for the independent OOH steady-state quality gate.

Run: python -m unittest -v test_ooh_numerics
"""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from run_ooh_analysis import (
    build_stoichiometric_matrices, check_steady_state,
    compare_refinement_runs, create_single_surface_setup
)


def gas_fluxes(rate):
    return {
        "C6H12O_g": rate,
        "C6H12_g": -rate,
        "H2O2_g": -rate,
        "H2O_g": rate,
    }


class NumericalGateTests(unittest.TestCase):
    def test_parsed_stoichiometry_has_expected_gas_balance(self):
        reactions = [
            [["s", "H2O2_g"], ["H2O2_s"]],
            [["H2O2_s", "s"], ["Ha_s", "OOH_s"]],
            [["s", "C6H12_g"], ["C6H12_s"]],
            [["OOH_s", "C6H12_s", "s"], ["C6H12O_s", "O_s", "Hb_s"]],
            [["C6H12O_s"], ["C6H12O_g", "s"]],
            [["Hb_s", "O_s"], ["OH_s", "s"]],
            [["Ha_s", "OH_s"], ["s", "H2O_s"]],
            [["H2O_s"], ["H2O_g", "s"]],
        ]
        matrix = build_stoichiometric_matrices(
            reactions,
            gas_species=("C6H12O_g", "C6H12_g", "H2O2_g", "H2O_g"),
        )
        self.assertEqual(matrix["gas_matrix"][0], [0, 0, 0, 0, 1, 0, 0, 0])
        self.assertEqual(matrix["gas_matrix"][1], [0, 0, -1, 0, 0, 0, 0, 0])
        self.assertEqual(matrix["gas_matrix"][2], [-1, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(matrix["gas_matrix"][3], [0, 0, 0, 0, 0, 0, 0, 1])

    def test_mathematically_consistent_cycle(self):
        rate = 1e-20
        r = check_steady_state(
            [rate] * 8, gas_fluxes(rate), [0.2, 0.8],
            solver_residual="1e-130", numbers_solver=True
        )
        self.assertEqual(r["quality_status"], "validated_single_precision")
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

    def test_ultralow_but_consistent_rate_is_not_arbitrarily_rejected(self):
        rate = 1e-60
        r = check_steady_state(
            [rate] * 8, gas_fluxes(rate), [0.2, 0.8],
            solver_residual="1e-140", numbers_solver=True
        )
        self.assertEqual(r["quality_status"], "validated_single_precision")
        self.assertTrue(r["quality_pass"])

    def test_absolute_residual_passes_but_relative_residual_fails(self):
        rate = 1e-100
        r = check_steady_state(
            [rate] * 8, gas_fluxes(rate), [0.2, 0.8],
            solver_residual="1e-120", numbers_solver=True
        )
        self.assertFalse(r["quality_pass"])
        self.assertGreater(r["steady_state_residual_relative_to_max_rate"], 1e20)

    def test_refinement_requires_two_independent_valid_solutions(self):
        a = {"quality_pass": True, "log10_net_C6H12O": -55.0}
        b = {"quality_pass": True, "log10_net_C6H12O": -55.005}
        c = {"quality_pass": True, "log10_net_C6H12O": -54.0}
        self.assertTrue(compare_refinement_runs(a, b)["refinement_stable"])
        self.assertFalse(compare_refinement_runs(a, c)["refinement_stable"])
        self.assertFalse(compare_refinement_runs(a, {"quality_pass": False,
            "log10_net_C6H12O": -55.0})["refinement_stable"])

    def test_mpmath_tolerance_never_underflows_to_float_zero(self):
        with TemporaryDirectory() as folder:
            path = create_single_surface_setup(
                "tiw", Path(folder), precision=460, tolerance="1e-340"
            )
            data = path.read_text(encoding="utf-8")
        self.assertIn("from mpmath import mp", data)
        self.assertIn("tolerance = mp.mpf('1e-340')", data)

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
