import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
MODEL_DIR = ROOT / "models" / "v1_serial"
sys.path.insert(0, str(MODEL_DIR))
spec = importlib.util.spec_from_file_location("validate_serial_network", MODEL_DIR / "validate_serial_network.py")
validate = importlib.util.module_from_spec(spec)
if spec.loader:
    spec.loader.exec_module(validate)


class SerialNetworkTests(unittest.TestCase):
    def test_catmap_031_parses_six_steps_and_identifies_two_ts(self):
        model = validate.build_model()
        self.assertEqual(len(model.elementary_rxns), 6)
        self.assertIn("TS1_s", model.transition_state_names)
        self.assertIn("TS2_s", model.transition_state_names)
        self.assertNotIn("TS1_s", model.adsorbate_names)
        self.assertNotIn("TS2_s", model.adsorbate_names)

    def test_catmap_balance_passes_for_explicit_single_center_assumption(self):
        model = validate.build_model()
        model.verify()
        report = validate.stoichiometry_report()
        self.assertEqual(report["matrix_rank"], 6)
        self.assertEqual(report["right_nullity"], 0)
        self.assertEqual(set(report["independent_conserved_totals"]), {"C", "H", "O", "site"})
        self.assertEqual(report["conserved_total_rank"], report["left_nullity"])
        self.assertEqual(report["left_nullity"], 4)

    def test_final_net_reaction_is_hexene_plus_peroxide_to_epoxide_plus_water(self):
        report = validate.stoichiometry_report()
        self.assertEqual(
            report["overall_net"],
            {"C6H12_g": -1, "H2O2_g": -1, "C6H12O_g": 1, "H2O_g": 1},
        )

    def test_coarse_grained_n_sites_are_never_reported_as_structure_verified(self):
        report = validate.model_metadata()
        self.assertEqual(report["site_occupancy_evidence"], "ASSUMED_SINGLE_EFFECTIVE_CENTER")
        self.assertFalse(report["kinetics_ready"])

    def test_missing_ts_barrier_uses_catmap_endpoint_assumption_and_is_labeled(self):
        model = validate.build_model()
        d_forward, ea_forward = model.get_rxn_energy(
            [["OHHa_s"], ["H2O_s"]], {"OHHa_s": 0.0, "H2O_s": 0.4}
        )
        d_reverse, ea_reverse = model.get_rxn_energy(
            [["H2O_s"], ["OHHa_s"]], {"OHHa_s": 0.0, "H2O_s": 0.4}
        )
        self.assertAlmostEqual(d_forward, 0.4)
        self.assertAlmostEqual(ea_forward, 0.4)
        self.assertAlmostEqual(d_reverse, -0.4)
        self.assertAlmostEqual(ea_reverse, 0.0)
        self.assertEqual(validate.UNACTIVATED_BARRIER_STATUS, "SOFTWARE_ASSUMPTION_NOT_DFT")


if __name__ == "__main__":
    unittest.main()
