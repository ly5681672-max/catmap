import sys
import unittest
from decimal import Decimal
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
from screen_tico_zero_inventory_coverages import trace_decision


class TraceCoverageRuleTests(unittest.TestCase):
    def test_no_trace_species(self):
        candidates, omitted, delta, delta_omit = trace_decision({"A_s": "0.2"}, {"A_s": 1}, True)
        self.assertEqual(candidates, set())
        self.assertEqual(omitted, set())
        self.assertEqual(delta, Decimal(0))
        self.assertEqual(delta_omit, Decimal(0))

    def test_single_trace_species_is_display_only(self):
        _, omitted, _, delta_omit = trace_decision({"A_s": "1e-8"}, {"A_s": 1}, True)
        self.assertEqual(omitted, {"A_s"})
        self.assertEqual(delta_omit, Decimal("1e-8"))

    def test_aggregate_limit_blocks_omission(self):
        # 101 individually trace species with two-site occupancy exceed the
        # aggregate cap while each coverage remains below 1e-6.
        coverages = {f"A{i}_s": "9.9e-7" for i in range(101)}
        occupancies = {species: 2 for species in coverages}
        _, omitted, delta, _ = trace_decision(coverages, occupancies, True)
        self.assertEqual(delta, Decimal("1.9998e-4"))
        self.assertEqual(omitted, set())

    def test_ooh_is_never_display_omitted(self):
        candidates, omitted, _, _ = trace_decision({"OOH_s": "1e-20", "A_s": "1e-20"}, {"OOH_s": 1, "A_s": 1}, True)
        self.assertIn("OOH_s", candidates)
        self.assertNotIn("OOH_s", omitted)
        self.assertIn("A_s", omitted)

    def test_bad_numerical_state_cannot_be_omitted(self):
        candidates, omitted, _, _ = trace_decision({"A_s": "1e-30"}, {"A_s": 1}, False)
        self.assertEqual(candidates, set())
        self.assertEqual(omitted, set())


if __name__ == "__main__":
    unittest.main()
