from __future__ import annotations

import sys
import unittest
from pathlib import Path


V2_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(V2_ROOT / "scripts"))
import construct_global_gibbs as subject  # noqa: E402


class SourceRowMappingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.crosswalk, cls.meta = subject.build_rows()
        cls.by_anchor = {
            (row["catalyst"], row["NEB_role"]): row for row in cls.crosswalk
        }

    def test_all_84_anchors_have_unique_workbook_source_cells(self):
        self.assertEqual(84, len(self.crosswalk))
        keys = [
            (r["catalyst"], r["source_workbook"], r["source_sheet"], r["source_row"], r["raw_G_cell"])
            for r in self.crosswalk
        ]
        self.assertEqual(84, len(set(keys)))
        self.assertTrue(all(r["formula_audit_match_count"] == "1" for r in self.crosswalk))

    def test_tire_six_anchors_match_source_cells_and_expected_signed_barriers(self):
        expected = {
            "NEB1_IS": ("AC12", "-986.726684", ""),
            "NEB1_TS": ("AC16", "-986.281838", "0.444846"),
            "NEB1_FS": ("AC20", "-988.654116", "2.372278"),
            "NEB2_IS": ("AC24", "-886.36569", ""),
            "NEB2_TS": ("AC28", "-886.147817", "0.217873"),
            "NEB2_FS": ("AC32", "-887.178525", "1.030708"),
        }
        for role, (cell, raw, _) in expected.items():
            row = self.by_anchor[("tire", role)]
            self.assertEqual(cell, row["raw_G_cell"])
            self.assertAlmostEqual(float(raw), float(row["raw_G_cached"]), places=9)
        self.assertEqual("0.444846", self.by_anchor[("tire", "NEB1_TS")]["mapped_TS_minus_IS_eV"])
        self.assertEqual("2.372278", self.by_anchor[("tire", "NEB1_TS")]["mapped_TS_minus_FS_eV"])
        self.assertEqual("0.217873", self.by_anchor[("tire", "NEB2_TS")]["mapped_TS_minus_IS_eV"])
        self.assertEqual("1.030708", self.by_anchor[("tire", "NEB2_TS")]["mapped_TS_minus_FS_eV"])

    def test_source_row_key_keeps_timo_conflicting_generic_role_separate(self):
        row = self.by_anchor[("timo", "NEB1_FS")]
        self.assertEqual("20", row["source_row"])
        self.assertEqual("NEB2_IS", row["generic_state_role"])
        self.assertEqual("ROLE_LABEL_CONFLICT_SOURCE_ROW_RETAINED", row["role_crosswalk_status"])
        self.assertEqual("NEB1_FS", row["NEB_role"])

    def test_timo_raw_and_adjusted_fs_energies_remain_separate(self):
        neb1 = self.by_anchor[("timo", "NEB1_FS")]
        neb2 = self.by_anchor[("timo", "NEB2_FS")]
        self.assertEqual("AR20", neb1["raw_G_cell"])
        self.assertAlmostEqual(-887.832689, float(neb1["raw_G_cached"]), places=9)
        self.assertEqual("AR21", neb1["mapped_neb_cell"])
        self.assertIn("AR20+F60", neb1["mapped_neb_formula_or_adjustment"])
        self.assertAlmostEqual(-988.870460, float(neb1["mapped_neb_local_energy_eV"]), places=9)
        self.assertEqual("AR32", neb2["raw_G_cell"])
        self.assertAlmostEqual(-989.164964, float(neb2["raw_G_cached"]), places=9)
        self.assertIn("spectator removed", neb2["mapped_neb_formula_or_adjustment"])
        self.assertAlmostEqual(-888.127193, float(neb2["mapped_neb_local_energy_eV"]), places=9)

    def test_tita_dynamic_columns_and_timn_index_anomaly_do_not_change_source_key(self):
        tita = self.by_anchor[("tita", "NEB2_TS")]
        self.assertEqual("CI28", tita["raw_G_cell"])
        self.assertEqual("CD28", tita["E0_cell"])
        self.assertEqual("CH28", tita["Gibbs_correction_cell"])
        timn = self.by_anchor[("timn", "NEB2_FS")]
        self.assertEqual("32", timn["source_row"])
        self.assertEqual("N32", timn["raw_G_cell"])
        self.assertEqual("总", timn["raw_state_index_in_workbook"])
        self.assertEqual("NEB2_FS", timn["NEB_role"])

    def test_missing_or_duplicate_source_matches_raise(self):
        with self.assertRaisesRegex(ValueError, "expected exactly one"):
            subject._require_exactly_one([], "missing workbook cell")
        with self.assertRaisesRegex(ValueError, "expected exactly one"):
            subject._require_exactly_one([1, 2], "duplicate workbook cell")
        with self.assertRaisesRegex(ValueError, "duplicate source key"):
            subject._unique_index(
                [{"k": "same"}, {"k": "same"}], lambda row: row["k"], "test source"
            )

    def test_global_formation_energy_stays_blank_for_all_records(self):
        self.assertTrue(all(not row["calculated_formation_energy"] for row in self.rows))
        self.assertTrue(all(not row["slab_reference"] for row in self.rows))

    def test_source_crosswalk_is_deterministic(self):
        rows_again, crosswalk_again, meta_again = subject.build_rows()
        self.assertEqual(self.rows, rows_again)
        self.assertEqual(self.crosswalk, crosswalk_again)
        self.assertEqual(self.meta, meta_again)


if __name__ == "__main__":
    unittest.main()
