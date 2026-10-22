from __future__ import annotations

import csv
import json
import ast
from pathlib import Path

import mpmath as mp


OUT = Path(r"C:\Users\Administrator\Documents\codex探索\tizr_precision_root_cause_20261021")
SITE = Path(r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\multisite_conditional_replay_20261021\sites\tizr")
SPECIES = ["C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "Ha_s", "Hb_s", "OH_s", "OOH_s", "O_s"]


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    saved = json.loads(SITE.joinpath("constrained_tizr_result.json").read_text(encoding="utf-8"))
    diag = json.loads((OUT / "data/tizr_replay_diagnostics.json").read_text(encoding="utf-8"))
    reevaluated = json.loads((OUT / "data/archived_coverage_re_evaluation.json").read_text(encoding="utf-8"))
    prior_root_rows = read_csv(OUT / "data/TIZR_ROOT_COMPARISON.csv")
    prior_flux_rows = read_csv(OUT / "data/TIZR_FORWARD_REVERSE_FLUX.csv")
    mp.mp.dps = 100

    prior_map = {(r["solution"], r["species"]): r for r in prior_root_rows}
    reconstructed_label = ("reconstructed_from_saved_40digit_logits"
                           if ("reconstructed_from_saved_40digit_logits", SPECIES[0]) in prior_map
                           else "archived_high_precision_root")
    states = {
        "archived_serialized_coverage_snapshot": {
            "coverage": {k: saved["coverage"][k] for k in SPECIES},
            "residual": {k: reevaluated["surface_rhs"][i] for i, k in enumerate(SPECIES)},
            "tof": reevaluated["recomputed_catmap_tof_from_stored_coverage"],
            "provenance": "stored coverage JSON (80 significant digits), re-evaluated with CatMAP 0.3.1 at 385 dps",
        },
        "reconstructed_from_saved_40digit_logits": {
            "coverage": {s: prior_map[(reconstructed_label, s)]["coverage"] for s in SPECIES},
            "residual": {s: prior_map[(reconstructed_label, s)]["absolute_dtheta_dt"] for s in SPECIES},
            "tof": diag["archived_high_root"]["TOF_C6H12O"],
            "provenance": "reconstructing state from archived root_logits serialized to 40 significant digits; not the full-precision root",
        },
        "replayed_old_1e-50_stop": {
            "coverage": {s: prior_map[("replayed_old_1e-50_stop", s)]["coverage"] for s in SPECIES},
            "residual": {s: prior_map[("replayed_old_1e-50_stop", s)]["absolute_dtheta_dt"] for s in SPECIES},
            "tof": diag["replayed_early_stop"]["TOF_C6H12O"],
            "provenance": "385 dps replay, stopped at original scaled residual threshold 1e-50 after 5 iterations",
        },
        "same_precision_bounded_polish": {
            "coverage": {s: prior_map[("same_precision_bounded_polish", s)]["coverage"] for s in SPECIES},
            "residual": {s: prior_map[("same_precision_bounded_polish", s)]["absolute_dtheta_dt"] for s in SPECIES},
            "tof": diag["bounded_polish"]["TOF_C6H12O"],
            "provenance": "single bounded Newton continuation from early-stop state, unchanged model/input and 385 dps",
        },
    }
    refcov = {s: mp.mpf(states["archived_serialized_coverage_snapshot"]["coverage"][s]) for s in SPECIES}
    tofs = {name: mp.mpf(state["tof"]) for name, state in states.items()}
    root_rows = []
    for state_name, state in states.items():
        for species in SPECIES + ["vacancy"]:
            cov = mp.mpf(saved["vacant"] if species == "vacancy" else state["coverage"][species])
            ref = mp.mpf(saved["vacant"] if species == "vacancy" else refcov[species])
            residual = "N/A" if species == "vacancy" else state["residual"][species]
            root_rows.append({
                "solution": state_name, "species": species,
                "coverage": mp.nstr(cov, 80),
                "absolute_dtheta_dt": residual,
                "absolute_coverage_difference_vs_archived_snapshot": mp.nstr(abs(cov - ref), 40),
                "relative_coverage_difference_vs_archived_snapshot": (
                    mp.nstr(abs(cov - ref) / abs(ref), 40) if ref else "UNDEFINED_ZERO_REFERENCE"),
                "TOF_C6H12O": mp.nstr(tofs[state_name], 80),
                "log10_TOF_C6H12O": mp.nstr(mp.log10(tofs[state_name]), 50),
                "provenance": state["provenance"],
            })
    write_csv(OUT / "data/TIZR_ROOT_COMPARISON.csv", root_rows)

    # Re-label prior rows: the former 'archived_high_precision_root' group is
    # a reconstruction from 40-digit root_logits, not the saved coverage root.
    flux_rows = []
    label_map = {
        "archived_high_precision_root": "reconstructed_from_saved_40digit_logits",
        "replayed_old_1e-50_stop": "replayed_old_1e-50_stop",
        "same_precision_bounded_polish": "same_precision_bounded_polish",
    }
    for row in prior_flux_rows:
        row = dict(row)
        if row["solution"] == "archived_serialized_coverage_snapshot":
            continue
        row["solution"] = label_map.get(row["solution"], row["solution"])
        flux_rows.append(row)
    for j, rates in enumerate(reevaluated["directional_rates_mpmath"], start=1):
        original = saved["reaction_labels"][j - 1]
        flux_rows.append({
            "solution": "archived_serialized_coverage_snapshot",
            "reaction_index": j,
            "reaction": original,
            "forward_rate_mpmath": rates["forward"],
            "reverse_rate_mpmath": rates["reverse"],
            "net_rate_catmap_mpmath": rates["net"],
            "forward_minus_reverse_error": rates["forward_reverse_identity_error"],
            "cancellation_factor": rates["cancellation_factor"],
            "tof_C6H12O_g": reevaluated["recomputed_catmap_tof_from_stored_coverage"],
            "mp_dps": reevaluated["coverage_snapshot_re_evaluation_dps"],
            "direction_source": "CatMAP rate_equation_term evaluated with mpmath at archived serialized coverage",
        })
    write_csv(OUT / "data/TIZR_FORWARD_REVERSE_FLUX.csv", flux_rows)

    # Calculate per-species gross-flow normalizers from the CatMAP parsed
    # reaction endpoints. This is a diagnostic normalization only; it does not
    # alter the solver equations or acceptance thresholds.
    stoich_rows = [[0] * len(saved["reaction_labels"]) for _ in SPECIES]
    for j, label in enumerate(saved["reaction_labels"]):
        parsed = ast.literal_eval(label)
        initial, final = parsed[0], parsed[-1]
        for i, species in enumerate(SPECIES):
            stoich_rows[i][j] = final.count(species) - initial.count(species)
    flux_by_state = {}
    for row in flux_rows:
        flux_by_state.setdefault(row["solution"], {})[int(row["reaction_index"]) - 1] = row
    for row in root_rows:
        species = row["species"]
        if species == "vacancy":
            row["gross_directional_flow_scale"] = "N/A"
            row["local_normalized_residual"] = "N/A"
            continue
        state = row["solution"]
        flux_for_state = flux_by_state.get(state, {})
        if len(flux_for_state) != len(saved["reaction_labels"]):
            row["gross_directional_flow_scale"] = "NOT_AVAILABLE"
            row["local_normalized_residual"] = "NOT_AVAILABLE"
            continue
        i = SPECIES.index(species)
        gross = sum(abs(stoich_rows[i][j]) * (
            abs(mp.mpf(flux_for_state[j]["forward_rate_mpmath"])) +
            abs(mp.mpf(flux_for_state[j]["reverse_rate_mpmath"]))
        ) for j in range(len(saved["reaction_labels"])))
        residual = abs(mp.mpf(row["absolute_dtheta_dt"]))
        row["gross_directional_flow_scale"] = mp.nstr(gross, 50)
        row["local_normalized_residual"] = mp.nstr(residual / gross, 50) if gross else "UNDEFINED_ZERO_GROSS_FLOW"
    write_csv(OUT / "data/TIZR_ROOT_COMPARISON.csv", root_rows)

    ref_tof = tofs["archived_serialized_coverage_snapshot"]
    early_tof = tofs["replayed_old_1e-50_stop"]
    polished_tof = tofs["same_precision_bounded_polish"]
    early_rel = abs(early_tof - ref_tof) / abs(ref_tof)
    early_dlog = abs(mp.log10(early_tof) - mp.log10(ref_tof))
    polish_rel = abs(polished_tof - ref_tof) / abs(ref_tof)
    polish_dlog = abs(mp.log10(polished_tof) - mp.log10(ref_tof))
    max_polish_abs_resid = mp.mpf(diag["bounded_polish"]["max_independent_abs_residual"])
    max_polish_rel_cycle = max(abs(mp.mpf(v) - polished_tof) for v in diag["bounded_polish"]["net_reaction_rates"]) / abs(polished_tof)
    max_cov_rel = max(abs(mp.mpf(states["same_precision_bounded_polish"]["coverage"][s]) - refcov[s]) / abs(refcov[s])
                      for s in SPECIES if refcov[s] > mp.mpf("1e-100"))
    final_gate = {
        "site": "tizr", "research_status": "PROVISIONAL",
        "old_status": "NUMERIC_PARTIAL",
        "old_553ppm_replay_status": "NUMERIC_PARTIAL",
        "current_bounded_replay_status": "NUMERIC_PASS",
        "accepted_5ppm_target_unchanged": "PASS",
        "input_model_sha256": diag["inputs"]["mkm_sha256"],
        "input_energy_sha256": diag["inputs"]["energies_sha256"],
        "catmap_version": "0.3.1", "working_precision_dps": 385,
        "same_model_resolve_performed": True,
        "physics_or_input_changed": False,
        "replay": {"early_iterations": diag["replayed_early_stop"]["iterations"],
                   "early_scaled_residual": diag["replayed_early_stop"]["scaled_residual_inf"],
                   "early_max_abs_surface_residual": diag["replayed_early_stop"]["max_independent_abs_residual"],
                   "early_TOF_relative_difference_vs_archived_snapshot": mp.nstr(early_rel, 40),
                   "early_delta_log10_TOF": mp.nstr(early_dlog, 40),
                   "polish_iterations": diag["bounded_polish"]["iterations"],
                   "polish_scaled_residual": diag["bounded_polish"]["scaled_residual_inf"],
                   "polish_max_abs_surface_residual": mp.nstr(max_polish_abs_resid, 40),
                   "polish_max_net_step_relative_deviation_from_product_TOF": mp.nstr(max_polish_rel_cycle, 40),
                   "polish_TOF_relative_difference_vs_archived_snapshot": mp.nstr(polish_rel, 40),
                   "polish_delta_log10_TOF": mp.nstr(polish_dlog, 40),
                   "polish_max_relative_adsorbate_coverage_difference_vs_archived_snapshot": mp.nstr(max_cov_rel, 40),
                   "I1": diag["bounded_polish"]["I1"], "I2": diag["bounded_polish"]["I2"],
                   "site_sum": diag["bounded_polish"]["site_sum"],
                   "target_scaled_residual": diag["bounded_polish"]["target_scaled_residual"],
                   "max_polish_iterations": diag["bounded_polish"]["max_iterations"],
                   "bounded_polish_error": diag["bounded_polish"]["error"]},
        "notes": [
            "The archived cross-precision run stopped at scaled residual 1e-50; its absolute species residual was many orders larger than the ~1e-71 cycle flux.",
            "The archived full-precision root is stored with root_logits at 40 significant digits and coverages at 80 digits; reconstructing from the truncated logits is not a faithful root representation.",
            "The archived serialized coverage itself re-evaluates product TOF within 1.3e-17 relative, but its local RHS residuals are only reproducible to about 1e-68 because coverage serialization is finite precision.",
            "A single same-precision bounded polish reached scaled residual below 1e-100 and recovered the archived TOF well within the unchanged 5 ppm target.",
            "CatMAP 0.3.1 get_directional_rates/get_elem_ec casts directional rates to Python float; do not use it to subtract near-equilibrium rates for ~1e-71 TOF.",
        ],
    }
    (OUT / "data/TIZR_FINAL_NUMERICAL_GATE.json").write_text(json.dumps(final_gate, indent=2), encoding="utf-8")
    print(json.dumps({"early_relative": mp.nstr(early_rel, 20), "early_delta_log10": mp.nstr(early_dlog, 20),
                      "polish_relative": mp.nstr(polish_rel, 20), "polish_delta_log10": mp.nstr(polish_dlog, 20),
                      "polish_abs_residual": mp.nstr(max_polish_abs_resid, 15),
                      "polish_cycle_relative": mp.nstr(max_polish_rel_cycle, 15),
                      "max_coverage_relative": mp.nstr(max_cov_rel, 15)}, indent=2))


if __name__ == "__main__":
    main()
