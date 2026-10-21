"""Regression tests for signed NEB endpoint classification and network gates."""
from __future__ import annotations

import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class AnomalyAndNetworkGateTests(unittest.TestCase):
    def read_rows(self, name: str) -> list[dict[str, str]]:
        with (ROOT / "data" / name).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))

    def test_seven_signed_neb_anomalies_preserve_values_and_endpoint_barriers(self) -> None:
        rows = self.read_rows("neb_endpoint_anomaly_audit.csv")
        expected = {
            ("ti", "NEB2"): ("TS_BELOW_FS", -0.1283235),
            ("tizr", "NEB2"): ("TS_BELOW_FS", -2.8199355),
            ("tiw", "NEB1"): ("TS_BELOW_IS", -0.566651),
            ("tita", "NEB1"): ("TS_BELOW_IS", -1.174636),
            ("timo", "NEB2"): ("TS_BELOW_IS", -0.012820),
            ("tiv", "NEB2"): ("TS_BELOW_IS", -0.024677),
            ("ticr", "NEB2"): ("TS_BELOW_IS", -0.1931555),
        }
        selected = {(r["surface_name"], r["path"]): r for r in rows
                    if (r["surface_name"], r["path"]) in expected}
        self.assertEqual(set(selected), set(expected))
        for key, (classification, signed) in expected.items():
            row = selected[key]
            self.assertEqual(row["TS_endpoint_classification"], classification)
            self.assertAlmostEqual(float(row["signed_TS_minus_IS_eV"] if classification == "TS_BELOW_IS"
                                         else row["signed_TS_minus_FS_eV"]), signed, places=8)
            self.assertAlmostEqual(float(row["recomputed_forward_endpoint_barrier_eV"]),
                                   max(float(row["G_IS_eV"]), float(row["G_TS_eV"]), float(row["G_FS_eV"]))
                                   - float(row["G_IS_eV"]), places=10)
            self.assertAlmostEqual(float(row["recomputed_reverse_endpoint_barrier_eV"]),
                                   max(float(row["G_IS_eV"]), float(row["G_TS_eV"]), float(row["G_FS_eV"]))
                                   - float(row["G_FS_eV"]), places=10)

    def test_original_and_candidate_inventory_dimensions_remain_distinct(self) -> None:
        original = self.read_rows("original_8step_inventory_invariants.csv")
        candidate = self.read_rows("candidate_9reaction_inventory_invariants.csv")
        self.assertEqual(len(original), 2)
        self.assertEqual(len(candidate), 1)
        checks = self.read_rows("network_balance_checks.csv")
        self.assertTrue(all(r["atom_balance"] == "PASS" and r["site_balance"] == "PASS"
                            for r in checks))


if __name__ == "__main__":
    unittest.main()
