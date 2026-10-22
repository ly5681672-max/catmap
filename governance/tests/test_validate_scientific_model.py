import hashlib
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from validate_scientific_model import validate


class ModelContractTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.model_rel = Path("model.mkm")
        self.energy_rel = Path("energies.txt")
        self.reactions = ["* + H2O2_g <-> H2O2*", "H2O2* + * -> Ha* + OOH*"]
        self.model = self.root / self.model_rel
        self.energy = self.root / self.energy_rel
        self.model.write_text(
            "rxn_expressions = " + repr(self.reactions) + "\n"
            "temperature = 333.15\npressure_mode = 'static'\n"
            "species_definitions = {}\n"
            "species_definitions['G_g'] = {'pressure': 1.0}\n", encoding="utf-8")
        self.energy.write_text("immutable input\n", encoding="utf-8")
        self.manifest = {
            "source_model": str(self.model_rel),
            "model_sha256": hashlib.sha256(self.model.read_bytes()).hexdigest(),
            "energy_input": str(self.energy_rel),
            "energy_input_sha256": hashlib.sha256(self.energy.read_bytes()).hexdigest(),
            "reaction_expressions": self.reactions,
            "temperature_K": 333.15,
            "pressure_mode": "static",
            "gas_inputs": {"G_g": {"pressure": 1.0}},
        }

    def tearDown(self):
        self.tmp.cleanup()

    def test_unchanged_baseline_passes(self):
        self.assertEqual(validate(self.root, self.manifest), [])

    def test_reaction_edit_is_blocked(self):
        self.model.write_text("rxn_expressions = " + repr(["* + H2O2_g <-> H2O2*"]) + "\n", encoding="utf-8")
        errors = validate(self.root, self.manifest)
        self.assertTrue(any("SHA256 differs" in error for error in errors))
        self.assertTrue(any("reaction expressions differ" in error for error in errors))

    def test_energy_edit_is_blocked(self):
        self.energy.write_text("changed scientific input\n", encoding="utf-8")
        self.assertTrue(any("energy input SHA256 differs" in error for error in validate(self.root, self.manifest)))

    def test_candidate_model_path_is_blocked_even_if_copy_is_identical(self):
        candidate = self.root / "candidate.mkm"
        candidate.write_bytes(self.model.read_bytes())
        errors = validate(self.root, self.manifest, model_arg=candidate)
        self.assertTrue(any("not the approved formal input" in error for error in errors))

    def test_reservoir_mode_change_is_blocked(self):
        source = self.model.read_text(encoding="utf-8")
        self.model.write_text(source.replace("pressure_mode = 'static'", "pressure_mode = 'concentration'"), encoding="utf-8")
        errors = validate(self.root, self.manifest)
        self.assertTrue(any("pressure_mode differs" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
