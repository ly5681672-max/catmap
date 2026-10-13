import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
MODEL_DIR = ROOT / "models" / "v1_serial"
sys.path.insert(0, str(MODEL_DIR))
spec = importlib.util.spec_from_file_location("validate_serial_energies", MODEL_DIR / "validate_serial_energies.py")
validate = importlib.util.module_from_spec(spec)
if spec.loader:
    spec.loader.exec_module(validate)


class SerialEnergyTests(unittest.TestCase):
    def test_tire_neb1_replays_both_signed_differences_and_catmap_barriers(self):
        row = validate.source_row("tire", "NEB1")
        result = validate.replay_neb_block(row, "NEB1")
        self.assertAlmostEqual(result["raw_TS_minus_IS_eV"], 0.444846, places=6)
        self.assertAlmostEqual(result["raw_TS_minus_FS_eV"], 2.372278, places=6)
        self.assertAlmostEqual(result["CatMAP_effective_Ea_forward_eV"], 0.444846, places=6)
        self.assertAlmostEqual(result["CatMAP_effective_Ea_reverse_eV"], 2.372278, places=6)

    def test_tire_neb2_replays_both_signed_differences_and_catmap_barriers(self):
        row = validate.source_row("tire", "NEB2")
        result = validate.replay_neb_block(row, "NEB2")
        self.assertAlmostEqual(result["raw_TS_minus_IS_eV"], 0.217873, places=6)
        self.assertAlmostEqual(result["raw_TS_minus_FS_eV"], 1.030708, places=6)
        self.assertAlmostEqual(result["CatMAP_effective_Ea_forward_eV"], 0.217873, places=6)
        self.assertAlmostEqual(result["CatMAP_effective_Ea_reverse_eV"], 1.030708, places=6)

    def test_common_energy_translation_preserves_neb_differences(self):
        row = validate.source_row("tire", "NEB1")
        original = validate.replay_neb_block(row, "NEB1")
        shifted = validate.replay_neb_block(row, "NEB1", common_shift_eV=137.25)
        for key in ("raw_TS_minus_IS_eV", "raw_TS_minus_FS_eV", "raw_FS_minus_IS_eV",
                    "CatMAP_effective_Ea_forward_eV", "CatMAP_effective_Ea_reverse_eV"):
            self.assertAlmostEqual(original[key], shifted[key], places=10)

    def test_epoxide_desorption_adds_gas_epoxide_exactly_once_and_stays_conditional(self):
        result = validate.epoxide_connection("tire")
        self.assertEqual(result["epoxide_gas_reference_count"], 1)
        self.assertAlmostEqual(result["delta_G_eV"], 1.250655, places=6)
        self.assertEqual(result["physical_status"], "CONDITIONAL_COMMON_SLAB_UNVERIFIED")
        self.assertEqual(result["activation_barrier_eV"], None)

    def test_local_block_energies_cannot_be_mislabeled_global_formation_energies(self):
        result = validate.replay_neb_block(validate.source_row("tire", "NEB1"), "NEB1")
        self.assertEqual(result["energy_reference_scope"], "BLOCK_LOCAL_NOT_GLOBAL_FORMATION_ENERGY")
        self.assertFalse(result["valid_for_coupled_network_rate_calculation"])

    def test_catmap_table_parser_accepts_diagnostic_energy_table_without_rate_run(self):
        result = validate.parse_diagnostic_table()
        self.assertEqual(result["catmap_version"], "0.3.1")
        self.assertEqual(result["parser_status"], "PASS_DIAGNOSTIC_TABLE_ONLY")
        self.assertFalse(result["coupled_rate_valid"])


if __name__ == "__main__":
    unittest.main()
