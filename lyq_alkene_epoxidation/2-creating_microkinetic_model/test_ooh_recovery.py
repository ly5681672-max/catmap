# -*- coding: utf-8 -*-
"""Pure-Python tests for temperature-continuation bookkeeping.

These tests do not run CatMAP; the physical solver must still be tested locally.
Run: python -m unittest -v test_ooh_numerics test_ooh_recovery
"""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from recover_ooh_solver import _boolean, _target_present, make_temperature_bridge


class RecoveryMetadataTests(unittest.TestCase):
    def test_target_detection_requires_both_temperature_and_pressure(self):
        self.assertTrue(_target_present([([333.15, 1.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([([450.0, 1.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([([333.15, 2.0], [0.2, 0.8])]))
        self.assertFalse(_target_present([]))

    def test_boolean_columns_require_explicit_true(self):
        raw = pd.Series([True, False, "True", "False", "NaN", None])
        self.assertEqual(_boolean(raw).tolist(),
                         [True, False, True, False, False, False])

    def test_invalid_temperature_is_rejected_without_catmap(self):
        with self.assertRaises(ValueError):
            make_temperature_bridge("titi", 333.15)
        with self.assertRaises(ValueError):
            make_temperature_bridge("titi", 500.0, steps=2)


if __name__ == "__main__":
    unittest.main()
