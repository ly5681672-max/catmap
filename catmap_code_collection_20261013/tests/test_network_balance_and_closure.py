import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("audit_catmap_031", SCRIPTS / "audit_catmap_031.py")
api = importlib.util.module_from_spec(spec)
spec.loader.exec_module(api)


class ProposedNetworkBalanceTests(unittest.TestCase):
    def test_candidate_serial_edges_parse_and_balance_when_assumed_site_counts_are_explicit(self):
        expressions = [
            "NEB1_IS_s <-> NEB1-TS_s -> NEB1_FS_s",
            "NEB1_FS_s + * -> C6H12O_g + Ha_s + Hb_s + O_s",
            "Ha_s + Hb_s + O_s <-> Ha_s + Hb-O_s + * -> OH_s + Ha_s + *",
            "OH_s + Ha_s -> H2O_s + *",
            "H2O_s -> H2O_g + *",
            "H2_g -> H2_g",
        ]
        definitions = {
            "NEB1_IS_s": {"composition": {"C": 6, "H": 14, "O": 2}, "n_sites": 2},
            "NEB1-TS_s": {"composition": {"C": 6, "H": 14, "O": 2}, "n_sites": 2},
            "NEB1_FS_s": {"composition": {"C": 6, "H": 14, "O": 2}, "n_sites": 2},
            "C6H12O_g": {"composition": {"C": 6, "H": 12, "O": 1}},
            "O_s": {"composition": {"O": 1}},
            "OH_s": {"composition": {"O": 1, "H": 1}},
            "H2O_s": {"composition": {"H": 2, "O": 1}},
            "H2O_g": {"composition": {"H": 2, "O": 1}},
            "H2_g": {"composition": {"H": 2}},
        }
        model = api.make_model(expressions, definitions)
        model.verify()
        self.assertEqual(len(model.elementary_rxns), len(expressions))

    def test_connector_without_vacancy_fails_site_balance(self):
        expressions = ["NEB1_FS_s -> C6H12O_g + Ha_s + Hb_s + O_s", "H2_g -> H2_g", "H2O_g -> H2O_g"]
        definitions = {
            "NEB1_FS_s": {"composition": {"C": 6, "H": 14, "O": 2}, "n_sites": 2},
            "C6H12O_g": {"composition": {"C": 6, "H": 12, "O": 1}},
            "O_s": {"composition": {"O": 1}},
            "H2_g": {"composition": {"H": 2}},
            "H2O_g": {"composition": {"H": 2, "O": 1}},
        }
        model = api.make_model(expressions, definitions)
        with self.assertRaisesRegex(ValueError, "Site balance"):
            model.verify()


if __name__ == "__main__":
    unittest.main()
