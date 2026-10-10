# -*- coding: utf-8 -*-
"""CatMAP recovery bookkeeping and verified warm-start tests (no CatMAP required).

Run:
  python -m unittest -v test_ooh_numerics test_ooh_recovery

Actual CatMAP rootfinding, temperature-map residuals and 333.15 K TOF
must be checked separately in the local CatMAP environment.
"""
import pickle
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from recover_ooh_solver import (
    _boolean,
    _target_present,
    _unique_temperature_candidates,
    _compare_independent_seeds,
    _store_single_verified_seed,
    make_temperature_bridge,
)


class RecoveryMetadataTests(unittest.TestCase):
    def test_target_detection_requires_both_temperature_and_pressure(self):
        self.assertTrue(_target_present([([333.15, 1.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([([450.0, 1.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([([333.15, 2.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([]))

    def test_boolean_columns_require_explicit_true(self):
        values = pd.Series([True, False, "True", "False", "NaN", None])
        self.assertEqual(_boolean(values).tolist(),
                         [True, False, True, False, False, False])

    def test_invalid_temperature_rejected_without_catmap(self):
        with self.assertRaises(ValueError):
            make_temperature_bridge("titi", 333.15)
        with self.assertRaises(ValueError):
            make_temperature_bridge("titi", 500.0, steps=2)

    def test_only_independently_validated_sources_can_be_selected(self):
        choices = [
            {"source_temperature_K": 333.15, "source_quality_pass": True,
             "seed_file": "target.pkl"},
            {"source_temperature_K": 360, "source_quality_pass": False,
             "seed_file": "failed.pkl"},
            {"source_temperature_K": 370, "source_quality_pass": True,
             "seed_file": "a.pkl"},
            {"source_temperature_K": 450, "source_quality_pass": True,
             "seed_file": "b.pkl"},
            {"source_temperature_K": 550, "source_quality_pass": True,
             "seed_file": "c.pkl"},
            {"source_temperature_K": 650, "source_quality_pass": True,
             "seed_file": "d.pkl"},
        ]
        picked = _unique_temperature_candidates(choices, max_seeds=3)
        self.assertEqual(len(picked), 3)
        self.assertEqual([v["source_temperature_K"] for v in picked],
                         [370, 450, 650])
        self.assertNotIn("failed.pkl", [x["seed_file"] for x in picked])

    def test_cross_seed_requires_two_valid_matching_roots(self):
        valid = {"quality_pass": True, "refinement_stable": True,
                 "log10_net_C6H12O": -80.0}
        agree = {"quality_pass": True, "refinement_stable": True,
                 "log10_net_C6H12O": -80.005}
        differ = {"quality_pass": True, "refinement_stable": True,
                  "log10_net_C6H12O": -82.0}
        self.assertFalse(_compare_independent_seeds([valid])["independent_seed_agreement"])
        self.assertTrue(_compare_independent_seeds([valid, agree])["independent_seed_agreement"])
        self.assertFalse(_compare_independent_seeds([valid, differ])["independent_seed_agreement"])

    def test_seed_pickle_contains_only_native_solution_maps(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "seed.pkl"
            _store_single_verified_seed([450., 1.], [0.3, 0.7],
                                        [-0.5, 0.0], path)
            with path.open("rb") as f:
                data = pickle.load(f)
            self.assertEqual(sorted(data), ["coverage_map", "numbers_map"])
            self.assertEqual(data["coverage_map"][0][0], [450., 1.])
            self.assertEqual(data["numbers_map"][0][1], [-0.5, 0.0])

    def test_missing_numbers_vector_must_fail_closed(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "seed.pkl"
            with self.assertRaises(ValueError):
                _store_single_verified_seed([450., 1.], [0.3], None, path)
            self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main()
