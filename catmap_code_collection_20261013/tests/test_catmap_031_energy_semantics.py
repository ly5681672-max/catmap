import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "audit_catmap_031.py"
spec = importlib.util.spec_from_file_location("audit_catmap_031", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


MODEL_PATH = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\alkene_epoxidation.mkm")


class CatMAP031Tests(unittest.TestCase):
    def test_local_model_parser_separates_ts_placeholders_from_stable_adsorbates(self):
        result = module.audit_original_model(MODEL_PATH)
        self.assertEqual(result["catmap_version"], "0.3.1")
        self.assertEqual(result["reaction_count"], 8)
        self.assertFalse(result["OOH-C6H12_is_stable_adsorbate"])
        self.assertTrue(result["OOH-C6H12_is_transition_state"])
        self.assertFalse(result["Hb-O_is_stable_adsorbate"])
        self.assertTrue(result["Hb-O_is_transition_state"])
        self.assertEqual(result["full_element_and_site_balance_verification"], "PASS")

    def test_catmap_effective_barriers_keep_signed_raw_differences_separate(self):
        from catmap.model import ReactionModel

        model = ReactionModel()
        tire_neb1 = module.barrier_values(model, -986.726684, -986.281838, -988.654116)
        self.assertAlmostEqual(tire_neb1["raw_TS_minus_IS"], 0.444846, places=6)
        self.assertAlmostEqual(tire_neb1["raw_TS_minus_FS"], 2.372278, places=6)
        self.assertAlmostEqual(tire_neb1["catmap_Ea_forward"], 0.444846, places=6)
        self.assertAlmostEqual(tire_neb1["catmap_Ea_reverse"], 2.372278, places=6)

        # Negative raw TS-IS remains visible; CatMAP's effective barrier is governed by the endpoint floor.
        tiw = module.barrier_values(model, 0.0, -0.566651, -2.676)
        self.assertAlmostEqual(tiw["raw_TS_minus_IS"], -0.566651, places=6)
        self.assertAlmostEqual(tiw["catmap_Ea_forward"], 0.0, places=6)
        self.assertAlmostEqual(tiw["catmap_Ea_reverse"], 2.676, places=6)

    def test_explicit_n_sites_is_used_by_catmap_balance_verifier(self):
        rxns = ["C_g + 2* -> C_s"]
        model = module.make_model(rxns)
        model.species_definitions["C_s"]["composition"] = {"C": 1}
        with self.assertRaisesRegex(ValueError, "Site balance"):
            model.verify()
        model.species_definitions["C_s"]["n_sites"] = 2
        model.verify()


if __name__ == "__main__":
    unittest.main()
