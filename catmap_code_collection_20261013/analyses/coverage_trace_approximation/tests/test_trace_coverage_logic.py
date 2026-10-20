import importlib.util
import sys
import unittest
from decimal import Decimal
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "screen_trace_coverage.py"
SPEC = importlib.util.spec_from_file_location("trace_screen", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class TraceCoverageRuleTests(unittest.TestCase):
    def setUp(self):
        self.n = {"A": Decimal(1), "B": Decimal(2), "C": Decimal(1)}

    def test_no_trace_species_preserves_full_balance(self):
        candidates, omitted, bound, decision = MODULE.decide_omissions(
            {"A": Decimal("1e-3")}, {"A": Decimal(1)})
        self.assertEqual(candidates, [])
        self.assertEqual(omitted, [])
        self.assertEqual(bound, 0)
        self.assertEqual(decision, "NO_TRACE_SPECIES")

    def test_single_trace_species_is_omittable(self):
        candidates, omitted, bound, decision = MODULE.decide_omissions(
            {"A": Decimal("2e-7")}, {"A": Decimal(1)})
        self.assertEqual(candidates, ["A"])
        self.assertEqual(omitted, ["A"])
        self.assertEqual(bound, Decimal("2e-7"))
        self.assertEqual(decision, "APPROXIMATION_ELIGIBLE")

    def test_multiple_trace_species_accumulate_site_weighted_bound(self):
        candidates, omitted, bound, decision = MODULE.decide_omissions(
            {"A": Decimal("2e-7"), "B": Decimal("3e-7")}, self.n)
        self.assertEqual(candidates, ["A", "B"])
        self.assertEqual(omitted, ["A", "B"])
        self.assertEqual(bound, Decimal("8e-7"))
        self.assertEqual(decision, "APPROXIMATION_ELIGIBLE")

    def test_individually_trace_but_cumulative_over_limit_are_kept(self):
        candidates, omitted, bound, decision = MODULE.decide_omissions(
            {"A": Decimal("8e-5"), "B": Decimal("1e-5")}, self.n,
            epsilon=Decimal("1e-4"), delta_max=Decimal("1e-4"))
        self.assertEqual(candidates, ["A", "B"])
        self.assertEqual(omitted, [])
        self.assertEqual(bound, Decimal("1e-4"))
        self.assertEqual(decision, "CUMULATIVE_OCCUPANCY_EXCEEDS_LIMIT")

    def test_any_case_above_single_threshold_excludes_species_globally(self):
        # The function receives the maximum over the complete comparison scope.
        candidates, omitted, _, _ = MODULE.decide_omissions(
            {"A": Decimal("1.0001e-6")}, self.n)
        self.assertEqual(candidates, [])
        self.assertEqual(omitted, [])

    def test_failed_data_quality_forbids_all_omission(self):
        candidates, omitted, _, decision = MODULE.decide_omissions(
            {"A": Decimal("1e-20")}, self.n, quality_ok=False)
        self.assertEqual(candidates, [])
        self.assertEqual(omitted, [])
        self.assertEqual(decision, "DATA_QUALITY_FAILED")

    def test_weighted_occupancy_uses_n_sites(self):
        candidates, omitted, bound, _ = MODULE.decide_omissions(
            {"B": Decimal("4e-5")}, self.n,
            epsilon=Decimal("1e-4"), delta_max=Decimal("1e-4"))
        self.assertEqual(candidates, ["B"])
        self.assertEqual(omitted, ["B"])
        self.assertEqual(bound, Decimal("8e-5"))

    def test_near_zero_vacancy_reports_large_relative_error(self):
        ratio = MODULE.delta_to_vacancy(Decimal("1e-5"), Decimal("1e-10"))
        self.assertEqual(ratio, Decimal("1e5"))
        self.assertTrue(ratio > 1)

    def test_zero_vacancy_relative_error_is_infinite_when_delta_nonzero(self):
        ratio = MODULE.delta_to_vacancy(Decimal("1e-8"), Decimal(0))
        self.assertTrue(ratio.is_infinite())

    def test_no_omitted_species_has_zero_omitted_to_vacancy_ratio(self):
        self.assertEqual(MODULE.delta_to_vacancy(Decimal(0), Decimal("1e-30")), Decimal(0))

    def test_nonfinite_and_negative_coverages_are_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.dec("NaN", "test coverage")
        _, errors = MODULE.validate_scope(
            [{"surface":"x", "species":"A", "coverage":"-1e-9", "is_stable_adsorbate":"True"},
             {"surface":"x", "species":"s", "coverage":"1.000000001", "is_stable_adsorbate":"False"}],
            [{"surface":"x", "status":"CONVERGED", "residual":"0", "tolerance":"1e-70"}],
            [{"surface":"x", "status":"CONVERGED", "total_site_coverage":"1", "vacancy_coverage":"1"}],
            {"A": Decimal(1)})
        self.assertTrue(errors)

    def test_archived_thirteen_site_data_passes_quality_and_balance(self):
        project = Path(__file__).resolve().parents[3]
        data = project / "analyses" / "microkinetic_mechanism_13sites" / "data"
        n_sites_all = MODULE.parse_n_sites(project / "catmap_official_run" / "alkene_epoxidation.mkm")
        cov = MODULE.read_csv(data / "coverage_13sites_long.csv")
        stable = {r["species"] for r in cov if r["is_stable_adsorbate"].lower() == "true"}
        n_sites = {s: n_sites_all[s] for s in stable}
        groups, errors = MODULE.validate_scope(
            cov,
            MODULE.read_csv(data / "convergence_13sites.csv"),
            MODULE.read_csv(data / "coverage_site_balance.csv"),
            n_sites,
        )
        self.assertEqual(len(groups), 13)
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
