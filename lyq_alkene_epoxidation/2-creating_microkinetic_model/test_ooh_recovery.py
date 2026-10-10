# -*- coding: utf-8 -*-
"""333.15 K multi-start recovery bookkeeping tests.

These tests do not require CatMAP.  Real CatMAP roots still require the local
CatMAP environment and are accepted only after the independent numerical
quality gates in ``run_ooh_analysis.py``.
"""
import pickle
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from recover_ooh_solver import (
    INITIAL_DOMINANTS,
    TARGET_PRESSURE,
    TARGET_TEMPERATURE,
    _boolean,
    _compare_independent_seeds,
    _coverage_to_squared_numbers,
    _make_initial_coverages,
    _store_single_verified_seed,
    _target_present,
)


class RecoveryMetadataTests(unittest.TestCase):
    def test_target_detection_requires_both_temperature_and_pressure(self):
        self.assertTrue(_target_present([([333.15, 1.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([([450.0, 1.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([([333.15, 2.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([]))

    def test_boolean_columns_require_explicit_true(self):
        values = pd.Series([True, False, "True", "False", "NaN", None])
        self.assertEqual(
            _boolean(values).tolist(), [True, False, True, False, False, False]
        )

    def test_initial_states_are_target_condition_only(self):
        self.assertEqual(TARGET_TEMPERATURE, 333.15)
        self.assertEqual(TARGET_PRESSURE, 1.0)
        self.assertEqual(len(INITIAL_DOMINANTS), 8)
        self.assertEqual(INITIAL_DOMINANTS[0], ("vacant", "vacant"))

    def test_coverage_seeds_sum_to_one_and_convert_to_numbers(self):
        adsorbates = ["OOH_s", "OH_s", "O_s"]
        for _, dominant in INITIAL_DOMINANTS[:4]:
            if dominant != "vacant" and dominant not in adsorbates:
                continue
            coverage = _make_initial_coverages(adsorbates, dominant)
            numbers = _coverage_to_squared_numbers(coverage)
            self.assertEqual(len(coverage), len(numbers))
            self.assertAlmostEqual(float(sum(coverage)), 1.0, places=12)
            self.assertTrue(all(float(value) >= 0 for value in numbers))

    def test_unknown_initial_species_fails_closed(self):
        with self.assertRaises(ValueError):
            _make_initial_coverages(["OOH_s"], "not_a_species")

    def test_cross_seed_requires_two_valid_matching_roots(self):
        valid = {
            "quality_pass": True,
            "refinement_stable": True,
            "log10_net_C6H12O": -80.0,
        }
        agree = {
            "quality_pass": True,
            "refinement_stable": True,
            "log10_net_C6H12O": -80.005,
        }
        differ = {
            "quality_pass": True,
            "refinement_stable": True,
            "log10_net_C6H12O": -82.0,
        }
        self.assertFalse(
            _compare_independent_seeds([valid])["independent_seed_agreement"]
        )
        self.assertTrue(
            _compare_independent_seeds([valid, agree])["independent_seed_agreement"]
        )
        self.assertFalse(
            _compare_independent_seeds([valid, differ])["independent_seed_agreement"]
        )

    def test_seed_pickle_contains_only_native_solution_maps(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "seed.pkl"
            _store_single_verified_seed(
                [333.15, 1.0], [0.3, 0.7], [0.5, 0.5], path
            )
            with path.open("rb") as fh:
                data = pickle.load(fh)
            self.assertEqual(sorted(data), ["coverage_map", "numbers_map"])
            self.assertEqual(data["coverage_map"][0][0], [333.15, 1.0])
            self.assertEqual(data["numbers_map"][0][1], [0.5, 0.5])

    def test_missing_numbers_vector_must_fail_closed(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "seed.pkl"
            with self.assertRaises(ValueError):
                _store_single_verified_seed([333.15, 1.0], [0.3], None, path)
            self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main()
