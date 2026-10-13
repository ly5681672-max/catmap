import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
MODEL_DIR = ROOT / "models" / "v1_serial"
sys.path.insert(0, str(MODEL_DIR))
spec = importlib.util.spec_from_file_location("serial_species", MODEL_DIR / "serial_species.py")
serial_species = importlib.util.module_from_spec(spec)
if spec.loader:
    spec.loader.exec_module(serial_species)


class SerialSpeciesTests(unittest.TestCase):
    def test_complete_coadsorbate_states_have_explicit_confirmed_compositions(self):
        definitions = serial_species.species_definitions()
        expected = {
            "Rpair_s": {"C": 6, "H": 14, "O": 2},
            "TS1_s": {"C": 6, "H": 14, "O": 2},
            "EpoxideHHO_s": {"C": 6, "H": 14, "O": 2},
            "HHO_s": {"H": 2, "O": 1},
            "TS2_s": {"H": 2, "O": 1},
            "OHHa_s": {"H": 2, "O": 1},
            "H2O_s": {"H": 2, "O": 1},
        }
        for species, composition in expected.items():
            self.assertEqual(definitions[species]["composition"], composition)
        self.assertEqual(definitions["EpoxideHHO_s"]["physical_interpretation"],
                         "adsorbed C6H12O plus surface Ha, Hb, O; no molecular H2O")

    def test_transition_states_are_not_stable_coverage_species(self):
        definitions = serial_species.species_definitions()
        self.assertEqual(definitions["TS1_s"]["role"], "transition_state")
        self.assertEqual(definitions["TS2_s"]["role"], "transition_state")
        self.assertNotIn("TS1_s", serial_species.STABLE_SURFACE_SPECIES)
        self.assertNotIn("TS2_s", serial_species.STABLE_SURFACE_SPECIES)

    def test_single_site_count_is_marked_as_modeling_assumption(self):
        for definition in serial_species.species_definitions().values():
            if definition.get("role") in {"stable_adsorbate", "transition_state"}:
                self.assertEqual(definition["n_sites"], 1)
                self.assertEqual(definition["n_sites_evidence"], "ASSUMED_SINGLE_EFFECTIVE_CENTER")

    def test_energy_conversion_refuses_unknown_common_slab_reference(self):
        converted = serial_species.convert_state_energy(
            raw_g_eV=-988.654116,
            source_cell="AC20",
            source_formula="=Y20+AB20",
            common_reference_eV=None,
            spectator_adjustments_eV=[],
        )
        self.assertIsNone(converted["formation_energy_eV"])
        self.assertEqual(converted["conversion_status"], "UNRESOLVED_REFERENCE")
        self.assertIn("common slab", converted["reason"])

    def test_energy_conversion_records_formula_and_signed_spectator_adjustment(self):
        converted = serial_species.convert_state_energy(
            raw_g_eV=-886.365690,
            source_cell="AC24",
            source_formula="=Y24+AB24",
            common_reference_eV=-900.0,
            spectator_adjustments_eV=[-101.037771],
            spectator_notes=["add C6H12O_g once"],
        )
        self.assertAlmostEqual(converted["adjusted_G_eV"], -987.403461)
        self.assertAlmostEqual(converted["formation_energy_eV"], -87.403461)
        self.assertEqual(converted["source_formula"], "=Y24+AB24")
        self.assertEqual(converted["spectator_notes"], ["add C6H12O_g once"])


if __name__ == "__main__":
    unittest.main()
