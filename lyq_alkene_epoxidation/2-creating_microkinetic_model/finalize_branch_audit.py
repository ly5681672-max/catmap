"""Combine 333.15 K branch candidates without altering earlier audit files."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import branch_audit as audit


ROOT = audit.AUDIT_ROOT
EXPANDED = ROOT / "expanded_multistart"
SOURCES = (
    ("unresolved_five", EXPANDED / "unresolved_five"),
    ("historical_five_resume_20261011", EXPANDED / "historical_five_resume_20261011"),
)
REQUIRED_STATE = (
    "coverage_order", "numbers_order", "full_coverage_vector_high_precision",
    "full_numbers_vector_high_precision", "net_rate_vector_high_precision",
    "gas_tof_vector_high_precision",
)


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _verified(row: dict) -> bool:
    if not all(str(row.get(key, "")).strip() for key in REQUIRED_STATE):
        return False
    return (
        audit._as_bool(row.get("quality_pass"))
        and audit._as_bool(row.get("refinement_stable"))
        and audit._as_bool(row.get("recomputed_quality_pass"))
        and audit._as_bool(row.get("parsed_stoichiometry_pass"))
        and audit._as_bool(row.get("recompute_cache_cleared"))
        and float(row.get("rate_recompute_relative_error", "inf")) <= audit.RECOMPUTE_RELATIVE_TOL
        and float(row.get("tof_recompute_relative_error", "inf")) <= audit.RECOMPUTE_RELATIVE_TOL
        and float(row.get("temperature_K", "nan")) == audit.TARGET_TEMPERATURE
    )


def prepare() -> pd.DataFrame:
    """Select independently checked candidates and cluster their full states."""
    old = _read(ROOT / "validated_roots.csv")
    old["source_tag"] = "original_replay"
    old["source_surface"] = old["source_surface"].astype(str)
    original_trials = _read(audit.TRIALS_CSV)
    original_trials["source_tag"] = "original_recovery"
    original_trials["root_id"] = original_trials.apply(
        lambda row: f"original_recovery__{row['surface']}__{row['seed_id']}__stage_{int(row['refinement_stage']):02d}",
        axis=1,
    )
    all_trials = [original_trials]
    candidates = old.to_dict("records")
    for tag, source_root in SOURCES:
        trials = _read(source_root / "expanded_seed_trials.csv")
        trials["source_tag"] = tag
        trials["root_id"] = trials.apply(
            lambda row: f"{tag}__{row['surface']}__{row['seed_id']}__stage_{int(row['refinement_stage']):02d}",
            axis=1,
        )
        all_trials.append(trials)
        for row in trials.to_dict("records"):
            if not _verified(row):
                continue
            stage = source_root / "target" / row["surface"] / row["seed_id"] / f"stage_{int(row['refinement_stage']):02d}"
            config = stage / row["surface"] / f"baseline_{row['surface']}.mkm"
            if not config.is_file():
                raise FileNotFoundError(config)
            row.update({
                "source_surface": row["surface"],
                "source_seed_id": row["seed_id"],
                "replay_root": str(stage.resolve()),
                "replay_quality_pass": "True",
                "replay_rate_self_consistent": "True",
                "replay_tof_self_consistent": "True",
                "source_precision": row["run_precision"],
                "source_tolerance": row["run_tolerance"],
            })
            candidates.append(row)
    all_rows = pd.concat(all_trials, ignore_index=True)
    all_rows.to_csv(ROOT / "combined_all_seed_trials.csv", index=False)
    validated = pd.DataFrame(candidates)
    if validated["root_id"].duplicated().any():
        raise ValueError("Duplicate combined root ID")
    validated.to_csv(ROOT / "combined_validated_roots.csv", index=False)
    clusters = pd.DataFrame(audit.cluster_roots(candidates))
    clusters["source_tag"] = clusters["root_id"].map(
        validated.set_index("root_id")["source_tag"]
    )
    clusters.to_csv(ROOT / "combined_root_clusters.csv", index=False)
    print(f"validated={len(validated)} branches={clusters[['surface', 'branch_id']].drop_duplicates().shape[0]} surfaces={clusters.surface.nunique()}")
    return clusters


def analyze() -> pd.DataFrame:
    """Evaluate one representative per distinct branch; resume completed rows."""
    roots = _read(ROOT / "combined_validated_roots.csv")
    clusters = _read(ROOT / "combined_root_clusters.csv")
    roots = roots.merge(clusters[["root_id", "branch_id"]], on="root_id", how="inner")
    representatives = roots.sort_values(
        ["surface", "branch_id", "source_tag", "root_id"], kind="stable"
    ).drop_duplicates(["surface", "branch_id"])
    path = ROOT / "combined_jacobian_eigenvalues.csv"
    previous = _read(path).to_dict("records") if path.is_file() else []
    for row in previous:
        if row.get("structural_zero_mode_count", ""):
            row.setdefault("compatibility_stability_status", row["stability_status"])
            if int(row["structural_zero_mode_count"]) > 0:
                row["stability_status"] = row["full_reduced_stability_status"]
    completed = {row["branch_id"] for row in previous}
    existing = _read(ROOT / "jacobian_eigenvalues.csv").set_index("root_id")
    for root in representatives.to_dict("records"):
        if root["branch_id"] in completed:
            continue
        print(f"[combined-branch-audit] jacobian {root['branch_id']} {root['root_id']}", flush=True)
        if root["root_id"] in existing.index:
            result = existing.loc[root["root_id"]].to_dict()
            result.update({"root_id": root["root_id"], "branch_id": root["branch_id"]})
            result["compatibility_stability_status"] = result["stability_status"]
            if int(result["structural_zero_mode_count"]) > 0:
                result["stability_status"] = result["full_reduced_stability_status"]
        else:
            try:
                result = audit.analyze_branch_stability(root)
            except Exception as exc:
                result = {
                    "root_id": root["root_id"],
                    "surface": root["surface"],
                    "branch_id": root["branch_id"],
                    "stability_status": "jacobian_evaluation_failed",
                    "jacobian_error": f"{type(exc).__name__}: {exc}",
                }
        previous.append(result)
        pd.DataFrame(previous).to_csv(path, index=False)
        completed.add(root["branch_id"])
    pd.DataFrame(previous).to_csv(path, index=False)
    return pd.DataFrame(previous)


def summarize() -> pd.DataFrame:
    """Report sampled solutions separately from zero-inventory TOF support."""
    roots = _read(ROOT / "combined_validated_roots.csv")
    clusters = _read(ROOT / "combined_root_clusters.csv")
    inventory = _read(ROOT / "conserved_inventory.csv")
    jacobian = _read(ROOT / "combined_jacobian_eigenvalues.csv")
    joined = clusters.merge(
        inventory[[
            "root_id", "surface", "I1_Ha_minus_OOH_OH_O", "I2_Hb_minus_O",
            "zero_inventory_compatible",
        ]],
        on=["root_id", "surface"], how="left",
    )
    if joined["zero_inventory_compatible"].eq("").any():
        raise ValueError("Missing conserved inventory for a validated root")
    branch_rows = []
    for (surface, branch_id), group in joined.groupby(["surface", "branch_id"], sort=True):
        zero_flags = set(group["zero_inventory_compatible"])
        if len(zero_flags) != 1:
            raise ValueError(f"Inconsistent conserved inventory within {branch_id}")
        eigen = jacobian.loc[jacobian["branch_id"] == branch_id]
        if len(eigen) != 1:
            raise ValueError(f"Missing or duplicate Jacobian for {branch_id}")
        representative = group.iloc[0]
        branch_rows.append({
            "surface": surface,
            "branch_id": branch_id,
            "sampled_root_count": len(group),
            "seed_ids": ";".join(group["seed_id"].tolist()),
            "root_ids": ";".join(group["root_id"].tolist()),
            "dominant_species": representative["dominant_species"],
            "log10_TOF": representative["log10_TOF"],
            "I1_Ha_minus_OOH_OH_O": representative["I1_Ha_minus_OOH_OH_O"],
            "I2_Hb_minus_O": representative["I2_Hb_minus_O"],
            "zero_inventory_compatible": zero_flags.pop(),
            "full_simplex_stability": eigen.iloc[0]["stability_status"],
            "compatibility_class_stability": eigen.iloc[0]["compatibility_stability_status"],
            "structural_zero_mode_count": eigen.iloc[0]["structural_zero_mode_count"],
            "jacobian_representative_root": eigen.iloc[0]["root_id"],
        })
    branches = pd.DataFrame(branch_rows)
    branches.to_csv(ROOT / "branch_inventory_summary.csv", index=False)
    summaries = []
    for surface, group in branches.groupby("surface", sort=True):
        zero = group.loc[group["zero_inventory_compatible"] == "True"]
        stable = zero.loc[zero["compatibility_class_stability"] == "linearly_stable"]
        zero_roots = int(zero["sampled_root_count"].sum())
        if len(zero) == 0:
            classification = "unresolved_zero_inventory"
        elif len(zero) > 1:
            classification = "multiple_zero_inventory_candidates"
        elif len(stable) == 1 and zero_roots >= 2:
            classification = "conditional_single_supported_branch"
        else:
            classification = "insufficient_evidence"
        summaries.append({
            "surface": surface,
            "validated_algebraic_roots": int((roots["source_surface"] == surface).sum()),
            "sampled_stationary_solution_clusters": len(group),
            "compatibility_stable_clusters": int((group["compatibility_class_stability"] == "linearly_stable").sum()),
            "full_simplex_neutral_clusters": int((group["full_simplex_stability"] == "near_neutral_or_undetermined").sum()),
            "zero_inventory_candidate_clusters": len(zero),
            "zero_inventory_reproduced_roots": zero_roots,
            "zero_inventory_branch_ids": ";".join(zero["branch_id"]),
            "zero_inventory_log10_TOF": ";".join(zero["log10_TOF"]),
            "zero_inventory_classification": classification,
            "official_TOF_action": "retain_prior_outputs_as_historical; no_automatic_update",
        })
    summary = pd.DataFrame(summaries)
    summary.to_csv(ROOT / "surface_branch_summary.csv", index=False)
    failure_frames = []
    for tag, path in (
        ("original_recovery", audit.TRIALS_CSV),
        *[(tag, source / "expanded_seed_trials.csv") for tag, source in SOURCES],
    ):
        frame = _read(path)
        failed = frame.loc[~frame["quality_pass"].map(audit._as_bool)].copy()
        failed["source_tag"] = tag
        failure_frames.append(failed)
    pd.concat(failure_frames, ignore_index=True).to_csv(
        ROOT / "numerical_failures.csv", index=False
    )
    print(summary[["surface", "sampled_stationary_solution_clusters", "zero_inventory_candidate_clusters", "zero_inventory_classification"]].to_string(index=False))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("prepare", "stability", "summarize"))
    args = parser.parse_args()
    if args.phase == "prepare":
        prepare()
    elif args.phase == "stability":
        rows = analyze()
        print(rows["stability_status"].value_counts().to_string())
    else:
        summarize()


if __name__ == "__main__":
    main()
