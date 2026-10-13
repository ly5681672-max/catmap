import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
MODEL_DIR = ROOT / "models" / "v1_serial"
sys.path.insert(0, str(MODEL_DIR))
spec = importlib.util.spec_from_file_location("map_14_surfaces", MODEL_DIR / "map_14_surfaces.py")
mapping = importlib.util.module_from_spec(spec)
if spec.loader:
    spec.loader.exec_module(mapping)


class SerialFourteenSurfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = mapping.build_surface_rows()
        cls.species_rows = mapping.build_species_energy_rows()

    def test_exact_fourteen_surfaces_and_two_nebs_are_present(self):
        self.assertEqual([row["catalyst"] for row in self.rows], mapping.CATALYSTS)
        self.assertNotIn("tinb", [row["catalyst"] for row in self.rows])
        self.assertTrue(all(row["NEB1_status"] == "MATCH" and row["NEB2_status"] == "MATCH"
                            for row in self.rows))

    def test_signed_negative_raw_barriers_are_preserved(self):
        by_name = {row["catalyst"]: row for row in self.rows}
        self.assertLess(float(by_name["tiw"]["NEB1_raw_TS_minus_IS_eV"]), 0)
        self.assertLess(float(by_name["tita"]["NEB1_raw_TS_minus_IS_eV"]), 0)
        self.assertLess(float(by_name["ticr"]["NEB2_raw_TS_minus_IS_eV"]), 0)
        self.assertGreaterEqual(float(by_name["tiw"]["NEB1_CatMAP_Ea_forward_eV"]), 0)

    def test_timo_spectator_transform_is_retained_and_global_energy_is_not_fabricated(self):
        by_name = {row["catalyst"]: row for row in self.rows}
        self.assertIn("spectator removed", by_name["timo"]["NEB2_FS_reference_note"])
        self.assertEqual(by_name["timo"]["common_slab_reference_status"], "UNRESOLVED")
        self.assertEqual(by_name["timo"]["formation_energy_status"], "UNRESOLVED_REFERENCE")

    def test_missing_post_neb2_states_are_explicit_not_imputed(self):
        by_name = {row["catalyst"]: row for row in self.rows}
        self.assertEqual(by_name["tiv"]["H2O_state_status"], "MISSING_SOURCE_STATE")
        self.assertEqual(by_name["tini"]["water_desorption_delta_G_eV"], "")

    def test_every_composite_state_has_blank_global_formation_energy(self):
        composite = [row for row in self.species_rows if row["species_key"] not in {
            "H2O2_g", "C6H12_g", "C6H12O_g", "H2O_g"
        }]
        self.assertTrue(composite)
        self.assertTrue(all(row["formation_energy_eV"] == "" for row in composite))
        self.assertTrue(all(row["formation_energy_status"] in {
            "UNRESOLVED_REFERENCE", "MISSING_SOURCE_STATE"
        } for row in composite))


if __name__ == "__main__":
    unittest.main()
