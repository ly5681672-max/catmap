import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = load("build_phase1_tables")
extract = load("extract_neb_gibbs")


class MappingNetworkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workbook = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\己烯环氧化.xlsx")
        cls.ts_csv = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\DFT_TS_relative.csv")
        cls.states, cls.gas = builder.read_cached_state_data(cls.workbook)
        cls.relative = builder.load_relative(cls.ts_csv)
        cls.rows = builder.mapping_rows(cls.states, cls.relative, cls.workbook)

    def test_tire_neb1_and_neb2_reproduce_signed_forward_reverse_differences(self):
        rows = {r["path"]: r for r in self.rows if r["catalyst"] == "tire"}
        self.assertAlmostEqual(rows["NEB1"]["raw_TS_minus_IS"], 0.444846, places=6)
        self.assertAlmostEqual(rows["NEB1"]["raw_TS_minus_FS"], 2.372278, places=6)
        self.assertAlmostEqual(rows["NEB1"]["CatMAP_effective_Ea_forward"], 0.444846, places=6)
        self.assertAlmostEqual(rows["NEB1"]["CatMAP_effective_Ea_reverse"], 2.372278, places=6)
        self.assertAlmostEqual(rows["NEB2"]["raw_TS_minus_IS"], 0.217873, places=6)
        self.assertAlmostEqual(rows["NEB2"]["raw_TS_minus_FS"], 1.030708, places=6)
        self.assertAlmostEqual(rows["NEB2"]["CatMAP_effective_Ea_forward"], 0.217873, places=6)
        self.assertAlmostEqual(rows["NEB2"]["CatMAP_effective_Ea_reverse"], 1.030708, places=6)

    def test_mapping_preserves_raw_negative_values_without_clipping(self):
        tire_rows = [r for r in self.rows if r["catalyst"] == "tiw" and r["path"] == "NEB1"]
        self.assertEqual(len(tire_rows), 1)
        self.assertLess(tire_rows[0]["raw_TS_minus_IS"], 0)
        self.assertGreaterEqual(tire_rows[0]["CatMAP_effective_Ea_forward"], 0)

    def test_all_fourteen_workbook_mappings_match_declared_csv_rows(self):
        self.assertEqual(len(self.rows), 28)
        self.assertTrue(all(r["source_energy_crosswalk_status"] == "MATCH" for r in self.rows))
        tire_mo = next(r for r in self.rows if r["catalyst"] == "timo" and r["path"] == "NEB1")
        self.assertIn("F60", tire_mo["FS_source_expression_or_adjustment"])
        tire_mo2 = next(r for r in self.rows if r["catalyst"] == "timo" and r["path"] == "NEB2")
        self.assertIn("spectator removed", tire_mo2["FS_source_expression_or_adjustment"])

    def test_candidate_serial_connection_delta_g_is_not_reported_as_a_barrier(self):
        rows = [r for r in self.rows if r["catalyst"] == "tire"]
        network = builder.network_rows(rows, self.states["tire"], self.gas)
        connection = next(r for r in network if r["step"] == "NEB1_FS_to_NEB2_IS")
        self.assertIn("candidate ΔG=1.250655 eV", connection["energy_or_barrier"])
        self.assertIn("connection TS missing", connection["energy_or_barrier"])
        self.assertIn("priority hypothesis only", connection["status"])

    def test_candidate_compositions_conserve_atoms_across_serial_connection(self):
        def total(*comps):
            result = {}
            for comp in comps:
                for atom, count in comp.items():
                    result[atom] = result.get(atom, 0) + count
            return result

        neb1_fs = {"C": 6, "H": 14, "O": 2}
        epoxide_g = {"C": 6, "H": 12, "O": 1}
        neb2_is = {"H": 2, "O": 1}
        self.assertEqual(total(neb1_fs), total(epoxide_g, neb2_is))


if __name__ == "__main__":
    unittest.main()
