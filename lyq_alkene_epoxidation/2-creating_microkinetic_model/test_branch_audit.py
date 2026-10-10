# -*- coding: utf-8 -*-
"""Tests for the independent 333.15 K root replay layer."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

import branch_audit
from audit_conserved_inventory import verify_network
from run_ooh_analysis import serialize_mp_vector


class BranchAuditUnitTests(unittest.TestCase):
    def test_candidate_filter_keeps_only_cross_precision_rows(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "seed_trials.csv"
            pd.DataFrame([
                {
                    "surface": "tire", "seed_id": "O_dominant",
                    "refinement_stage": 1, "quality_pass": True,
                    "refinement_stable": False, "run_precision": 260,
                    "run_tolerance": "1e-180", "initial_seed_source": "a.pkl",
                },
                {
                    "surface": "tire", "seed_id": "O_dominant",
                    "refinement_stage": 2, "quality_pass": True,
                    "refinement_stable": True, "run_precision": 360,
                    "run_tolerance": "1e-260", "initial_seed_source": "a.pkl",
                },
                {
                    "surface": "tire", "seed_id": "Hb_dominant",
                    "refinement_stage": 2, "quality_pass": False,
                    "refinement_stable": False, "run_precision": 360,
                    "run_tolerance": "1e-260", "initial_seed_source": "b.pkl",
                },
            ]).to_csv(path, index=False)
            rows = branch_audit.read_cross_precision_candidates(path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["seed_id"], "O_dominant")
        self.assertEqual(rows[0]["run_precision"], 360)

    def test_high_precision_vector_serialization_does_not_use_float(self):
        from mpmath import mp

        mp.dps = 100
        text = serialize_mp_vector([mp.mpf("1e-340"), mp.mpf("0.999")])
        self.assertIn("1.0e-340", text)
        self.assertNotIn("0.0", text.split(";")[0])

    def test_root_id_is_deterministic(self):
        self.assertEqual(
            branch_audit._root_id({
                "surface": "titi", "seed_id": "Ha_dominant",
                "refinement_stage": "3",
            }),
            "titi__Ha_dominant__stage_03",
        )

    def test_same_tof_but_different_full_coverage_is_not_merged(self):
        base = {
            "surface": "tire",
            "coverage_order": "O_s;Hb_s;vacant",
            "net_rate_vector_high_precision": "1e-20;1e-20;1e-20",
            "log10_net_C6H12O": "-32.7298003469",
        }
        left = {
            **base,
            "full_coverage_vector_high_precision": "0.999999;1e-12;1e-6",
        }
        right = {
            **base,
            "full_coverage_vector_high_precision": "1e-12;0.999999;1e-6",
        }
        same, metrics = branch_audit._same_branch(left, right)
        self.assertFalse(same)
        self.assertGreater(float(metrics["coverage_abs_max"]), 0.9)

    def test_same_full_state_and_rates_are_merged(self):
        base = {
            "surface": "timo",
            "coverage_order": "O_s;vacant",
            "full_coverage_vector_high_precision": "0.999999999;1e-9",
            "net_rate_vector_high_precision": "1e-40;1e-40",
            "log10_net_C6H12O": "-40.0",
        }
        left = {**base}
        right = {**base, "log10_net_C6H12O": "-40.005"}
        same, _ = branch_audit._same_branch(left, right)
        self.assertTrue(same)

    def test_exact_stoichiometric_rank_removes_structural_modes(self):
        rows = [
            [1, 0, 1],
            [0, 1, 1],
            [-1, -1, -2],
        ]
        columns = branch_audit._independent_column_indices(rows)
        self.assertEqual(columns, [0, 1])
        self.assertEqual(branch_audit._exact_rank(rows), 2)

    def test_compatibility_jacobian_has_stoichiometric_dimension(self):
        from mpmath import mp

        rows = [[1, 0, 1], [0, 1, 1], [-1, -1, -2]]
        jacobian = mp.matrix(3, 3)
        for i in range(3):
            jacobian[i, i] = -1
        reduced, basis = branch_audit._compatibility_class_jacobian(
            jacobian, rows
        )
        self.assertEqual(basis, [0, 1])
        self.assertEqual((reduced.rows, reduced.cols), (2, 2))

    def test_stability_classifier_distinguishes_sign_and_precision_zero(self):
        from mpmath import mp

        self.assertEqual(
            branch_audit._classify_linear_stability(
                [mp.mpf("-1"), mp.mpf("-1e-20")], 100, mp.mpf("1e-80")
            ),
            "linearly_stable",
        )
        self.assertEqual(
            branch_audit._classify_linear_stability(
                [mp.mpf("1e-10"), mp.mpf("-1")], 100, mp.mpf("1e-80")
            ),
            "linearly_unstable",
        )
        self.assertEqual(
            branch_audit._classify_linear_stability(
                [mp.mpf("1e-60"), mp.mpf("-1")], 100, mp.mpf("1e-80")
            ),
            "near_neutral_or_undetermined",
        )

    def test_conserved_modes_are_not_promoted_to_full_stability(self):
        self.assertEqual(
            branch_audit._classify_full_simplex_stability(
                "near_neutral_or_undetermined", "linearly_stable", 2
            ),
            "near_neutral_or_undetermined",
        )
        self.assertEqual(
            branch_audit._classify_full_simplex_stability(
                "near_neutral_or_undetermined", "linearly_stable", 0
            ),
            "linearly_stable",
        )

    def test_actual_catmap_network_has_two_exact_invariants(self):
        network = verify_network(Path("alkene_epoxidation.mkm"))
        self.assertEqual(network["reaction_count"], 8)
        self.assertEqual(network["stoichiometric_rank"], 7)
        self.assertEqual(network["independent_conserved_inventories"], 2)


if __name__ == "__main__":
    unittest.main()
