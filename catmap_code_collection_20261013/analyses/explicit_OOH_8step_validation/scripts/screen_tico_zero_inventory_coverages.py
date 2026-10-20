"""Screen stored Ti-Co explicit-OOH coverages; never runs a CatMAP solve."""

from __future__ import annotations

import argparse
import csv
import json
import os
from decimal import Decimal
from pathlib import Path

from catmap import ReactionModel


TRACE_THRESHOLD = Decimal("1e-6")
AGGREGATE_THRESHOLD = Decimal("1e-4")
PROTECTED_SPECIES = {"OOH_s"}
INVARIANTS = {
    "I1_Ha_minus_OOH_OH_O": {"Ha_s": 1, "OOH_s": -1, "OH_s": -1, "O_s": -1},
    "I2_Hb_minus_O": {"Hb_s": 1, "O_s": -1},
}


def trace_decision(coverages, site_occupancies, numerical_quality_pass, protected=PROTECTED_SPECIES):
    """Return trace candidates and display-only omissions for a single valid state."""
    candidates = {
        species for species, value in coverages.items()
        if Decimal(value) < TRACE_THRESHOLD
    } if numerical_quality_pass else set()
    delta_all = sum(
        (Decimal(site_occupancies[s]) * Decimal(coverages[s]) for s in candidates),
        Decimal(0),
    )
    omit = set()
    if numerical_quality_pass and delta_all < AGGREGATE_THRESHOLD:
        omit = candidates - set(protected)
    delta_omit = sum(
        (Decimal(site_occupancies[s]) * Decimal(coverages[s]) for s in omit),
        Decimal(0),
    )
    return candidates, omit, delta_all, delta_omit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--root-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()

    model_dir = args.model_dir.resolve()
    root_path = args.root_json.resolve()
    output_path = args.output.resolve()
    root = json.loads(root_path.read_text(encoding="utf-8"))

    old_cwd = Path.cwd()
    try:
        os.chdir(model_dir)
        model = ReactionModel(setup_file=str(model_dir / "alkene_epoxidation.mkm"))
    finally:
        os.chdir(old_cwd)

    adsorbates = list(model.adsorbate_names)
    stable_occupancies = {
        species: int(model.species_definitions[species].get("n_sites", 1))
        for species in adsorbates
    }
    if len(model.elementary_rxns) != 8 or len(model.transition_state_names) != 2:
        raise RuntimeError("CatMAP parser did not return the expected original 8-step model.")
    if set(root["coverage"]) != set(adsorbates):
        raise RuntimeError("Stored root species do not match CatMAP-parsed stable adsorbates.")

    # Build the stable-adsorbate stoichiometric matrix from CatMAP's parsed
    # elementary reactions, excluding transition-state intermediates.
    stoich = {sp: [0] * len(model.elementary_rxns) for sp in adsorbates}
    for j, reaction in enumerate(model.elementary_rxns):
        for sp in reaction[0]:
            if sp in stoich:
                stoich[sp][j] -= 1
        for sp in reaction[-1]:
            if sp in stoich:
                stoich[sp][j] += 1
    for name, coeffs in INVARIANTS.items():
        if any(sum(coeffs.get(sp, 0) * stoich[sp][j] for sp in adsorbates) != 0
               for j in range(len(model.elementary_rxns))):
            raise RuntimeError(f"Inventory invariant {name} failed against parsed stoichiometry.")

    coverages = {sp: Decimal(root["coverage"][sp]) for sp in adsorbates}
    vacancy = Decimal(root["vacant"])
    site_sum = vacancy + sum(
        (Decimal(stable_occupancies[sp]) * coverages[sp] for sp in adsorbates), Decimal(0)
    )
    inventory_values = {name: Decimal(root[name]) for name in ("I1", "I2")}
    residuals = [abs(Decimal(v)) for v in root["surface_residuals"].values()]
    net_rates = [Decimal(v) for v in root["reaction_net_rates"]]
    positive_tof = Decimal(root["C6H12O_net_tof"])
    max_rate = max((abs(v) for v in net_rates), default=Decimal(0))
    cycle_relative = max((abs(v - net_rates[0]) for v in net_rates), default=Decimal("Infinity")) / max_rate if max_rate else Decimal("Infinity")
    numeric_checks = {
        "stored_solver_quality_gate": root.get("solver_converged_and_validated") is True,
        "finite_nonnegative_coverages_and_vacancy": all(v.is_finite() and v >= 0 for v in coverages.values()) and vacancy.is_finite() and vacancy >= 0,
        "site_sum_within_1e-50": abs(site_sum - Decimal(1)) <= Decimal("1e-50"),
        "zero_inventory_within_1e-100": max(abs(v) for v in inventory_values.values()) <= Decimal("1e-100"),
        "surface_residual_within_1e-90": max(residuals) <= Decimal("1e-90"),
        "relative_residual_within_1e-60": Decimal(root["relative_residual_to_max_reaction_rate"]) <= Decimal("1e-60"),
        "eight_net_rates_finite_and_cycle_closed": len(net_rates) == 8 and all(v.is_finite() for v in net_rates) and cycle_relative <= Decimal("1e-60"),
        "gas_and_stoichiometry_reconstruction_within_1e-50": Decimal(root["gas_stoichiometry_reconstruction_abs"]) <= Decimal("1e-50") and Decimal(root["stoichiometry_flux_relative_error"]) <= Decimal("1e-50"),
        "cross_precision_replay_passed": root.get("cross_precision_replay", {}).get("passed_existing_steady_quality_gate") is True,
        "positive_product_net_tof": positive_tof.is_finite() and positive_tof > 0,
    }
    numeric_ok = all(numeric_checks.values())

    candidate_names, omit_names, delta_all, delta_omit = trace_decision(
        {sp: str(v) for sp, v in coverages.items()},
        {sp: str(n) for sp, n in stable_occupancies.items()},
        numeric_ok,
    )
    delta_over_vacancy = delta_omit / vacancy if vacancy else Decimal("Infinity")
    cohort_status = "UNVERIFIED_INCOMPLETE_MULTI_CATALYST_COHORT"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "model_scope", "surface_name", "species_internal_id", "species_role",
        "n_sites", "coverage", "site_occupancy_contribution", "local_min_coverage",
        "local_max_coverage", "local_trace_candidate_lt_1e-6", "local_display_omit",
        "local_screen_status", "local_omitted_delta_all_trace", "local_omitted_delta_display",
        "local_delta_over_vacancy", "global_trace_candidate", "global_display_omit",
        "global_screen_status", "numerical_quality_pass", "physical_readiness",
    ]
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for sp in adsorbates:
            value = coverages[sp]
            writer.writerow({
                "model_scope": "explicit_OOH_8step; provisional Ti-Co zero-inventory root",
                "surface_name": "tico",
                "species_internal_id": sp,
                "species_role": "stable_adsorbate",
                "n_sites": stable_occupancies[sp],
                "coverage": str(value),
                "site_occupancy_contribution": str(value * stable_occupancies[sp]),
                "local_min_coverage": str(value),
                "local_max_coverage": str(value),
                "local_trace_candidate_lt_1e-6": str(sp in candidate_names).upper(),
                "local_display_omit": str(sp in omit_names).upper(),
                "local_screen_status": "PROVISIONAL_NUMERICALLY_QUALIFIED" if numeric_ok else "FAIL_NUMERIC_QUALITY",
                "local_omitted_delta_all_trace": str(delta_all),
                "local_omitted_delta_display": str(delta_omit),
                "local_delta_over_vacancy": str(delta_over_vacancy),
                "global_trace_candidate": "UNVERIFIED",
                "global_display_omit": "NO",
                "global_screen_status": cohort_status,
                "numerical_quality_pass": str(numeric_ok).upper(),
                "physical_readiness": "PROVISIONAL; NEB mapping and global reachability unresolved",
            })

    summary = {
        "catmap_version": "0.3.1",
        "parsed_reaction_count": len(model.elementary_rxns),
        "stable_adsorbates": adsorbates,
        "transition_states_excluded_from_coverage_balance": list(model.transition_state_names),
        "n_sites": stable_occupancies,
        "numerical_quality_pass": numeric_ok,
        "numerical_quality_checks": numeric_checks,
        "site_sum_recomputed": str(site_sum),
        "I1": str(inventory_values["I1"]),
        "I2": str(inventory_values["I2"]),
        "local_trace_candidates": sorted(candidate_names),
        "protected_species": sorted(PROTECTED_SPECIES),
        "local_display_omissions": sorted(omit_names),
        "delta_all_trace_candidates": str(delta_all),
        "delta_display_omissions": str(delta_omit),
        "vacancy": str(vacancy),
        "delta_over_vacancy": str(delta_over_vacancy),
        "global_screen_status": cohort_status,
        "tof_s-1_provisional": root["C6H12O_net_tof"],
        "log10_tof_provisional": root["log10_C6H12O_net_tof"],
    }
    summary_path = args.summary.resolve() if args.summary else output_path.with_name("eight_step_parse_screen_summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
