from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]


def csv_rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


class CandidateOutputTests(unittest.TestCase):
    def test_tife_candidate_steady_gate(self) -> None:
        run = ROOT / "runs" / "tife_candidate"
        gate = json.loads((run / "numerical_gate.json").read_text(encoding="utf-8"))
        self.assertEqual(gate["status"], "NUMERICAL_PASS_PROVISIONAL_PHYSICS")
        self.assertEqual(gate["accepted_roots"], 2)
        self.assertLess(mp.mpf(gate["max_relative_surface_residual"]), mp.mpf("1e-25"))
        self.assertLess(mp.mpf(gate["site_balance_error"]), mp.mpf("1e-100"))
        self.assertLess(mp.mpf(gate["inventory_error"]), mp.mpf("1e-100"))
        rows = csv_rows(run / "candidate_steady_state_coverages.csv")
        values = {r["species"]: mp.mpf(r["coverage"]) for r in rows}
        self.assertLess(abs(sum(values.values()) - 1), mp.mpf("1e-100"))
        self.assertEqual(values["Ha_s"], values["OOH_s"])
        self.assertTrue(all(v >= 0 for v in values.values()))

    def test_tife_nine_net_fluxes_close_to_tof(self) -> None:
        run = ROOT / "runs" / "tife_candidate"
        tof = mp.mpf(next(r["turnover_frequency"] for r in csv_rows(run / "candidate_tof.csv")
                          if r["gas_species"] == "C6H12O_g"))
        rates = [mp.mpf(r["net_rate"]) for r in csv_rows(run / "candidate_reaction_rates.csv")]
        self.assertEqual(len(rates), 9)
        self.assertTrue(all(abs(rate - tof) / tof < mp.mpf("1e-45") for rate in rates))

    def test_detailed_balance_and_catmap_arrow_semantics_recorded(self) -> None:
        result = json.loads((ROOT / "data" / "tife_detailed_balance_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(result["reversible_steps_checked"], 9)
        self.assertEqual(result["reversible_steps_mismatched"], 0)
        self.assertEqual(result["zero_reverse_constant_steps"], 0)
        self.assertIn("single-arrow", result["step2_arrow_vs_rate_constant"])

    def test_all_conditional_sensitivity_cases_are_present(self) -> None:
        manifest = json.loads((ROOT / "data" / "tife_sensitivity_run_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(len(manifest), 8)
        self.assertTrue(all(r["status"] == "NUMERICAL_PASS_PROVISIONAL_PHYSICS" for r in manifest))
        table = csv_rows(ROOT / "data" / "tife_conditional_sensitivity.csv")
        self.assertEqual(len(table), 4)
        self.assertTrue(any(r["kind"] == "ENDPOINT_ENERGY_SENSITIVITY_NOT_DRC" for r in table))

    def test_roster_keeps_target_sites_and_tinb_separate(self) -> None:
        roster = csv_rows(ROOT / "data" / "site_readiness.csv")
        self.assertEqual(len(roster), 14)
        self.assertNotIn("tinb", {r["surface_name"] for r in roster})
        self.assertTrue(all(r["independent_inputs_complete"] == "True" for r in roster))
        self.assertTrue(all(r["neb1_neb2_all_six_values_present"] == "True" for r in roster))


if __name__ == "__main__":
    unittest.main()
