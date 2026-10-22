from __future__ import annotations

import csv
import json
from pathlib import Path

BASE = Path(r"C:\Users\Administrator\Documents\codex探索\catmap_original8_zero_inventory_multisite")
TICO = Path(r"C:\Users\Administrator\Documents\codex探索\catmap_tico_zero_inventory_static_replay\constrained_tico_result.json")
SITES = ["ti", "tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "titi", "tiw", "tita", "ticr", "tini"]


def load_site(site: str) -> tuple[dict | None, Path]:
    candidates = [BASE / site, BASE / f"retry_{site}"]
    if site == "tico":
        candidates = [TICO.parent]
    for folder in candidates:
        p = folder / f"constrained_{site}_result.json"
        if p.exists():
            return json.loads(p.read_text(encoding="utf-8")), folder
    return None, candidates[0]


def main() -> None:
    rows = []
    sites = ["tico", *SITES]
    for site in sites:
        result, folder = load_site(site)
        drc = folder / "conditional_drc.csv"
        attempts = folder / "constrained_solver_attempts.json"
        if result is None:
            attempts_data = json.loads(attempts.read_text(encoding="utf-8")) if attempts.exists() else []
            rows.append({
                "surface": site, "run_state": "NUMERICAL_FAILED_OR_UNRESOLVED",
                "TOF_C6H12O_s-1": "", "log10_TOF": "", "theta_OOH": "", "theta_vacancy": "",
                "dominant_adsorbate": "", "dominant_coverage": "", "I1": "", "I2": "",
                "site_sum": "", "relative_residual": "", "cross_precision_dlog10TOF": "",
                "cross_precision_gate": "", "local_stable": "", "DRC_state": "NOT_RUN_NO_ACCEPTED_ROOT",
                "attempt_count": len(attempts_data), "diagnostic": "No accepted result JSON; inspect retained solver attempts/logs",
                "result_directory": str(folder),
            })
            continue
        cov = {k: float(v) for k, v in result["coverage"].items()}
        cross = result.get("cross_precision_replay") or {}
        drc_state = "PROVISIONAL_SAME_INVENTORY_DRC_PRESENT" if drc.exists() else "NOT_PRESENT"
        rows.append({
            "surface": site,
            "run_state": "ROOT_ACCEPTED_CROSS_PRECISION_PASS" if cross.get("passed_existing_steady_quality_gate") else "ROOT_ACCEPTED_CROSS_PRECISION_GATE_NOT_PASSED",
            "TOF_C6H12O_s-1": result.get("C6H12O_net_tof", ""),
            "log10_TOF": result.get("log10_C6H12O_net_tof", ""),
            "theta_OOH": result["coverage"].get("OOH_s", ""),
            "theta_vacancy": result.get("vacant", ""),
            "dominant_adsorbate": max(cov, key=cov.get),
            "dominant_coverage": max(cov.values()),
            "I1": result.get("I1", ""), "I2": result.get("I2", ""),
            "site_sum": result.get("site_sum", ""),
            "relative_residual": result.get("relative_residual_to_max_reaction_rate", ""),
            "cross_precision_dlog10TOF": cross.get("delta_log10_TOF_vs_high_precision", ""),
            "cross_precision_gate": cross.get("passed_existing_steady_quality_gate", ""),
            "local_stable": result.get("locally_stable_in_zero_inventory_subspace", ""),
            "DRC_state": drc_state,
            "attempt_count": "",
            "diagnostic": "CONDITIONAL original-8-step numerical result; NEB physical mapping remains unvalidated",
            "result_directory": str(folder),
        })
    out = BASE / "MULTISITE_CONDITIONAL_SUMMARY.csv"
    with out.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    good_root = sum(r["TOF_C6H12O_s-1"] != "" for r in rows)
    cross_pass = sum(r["cross_precision_gate"] is True for r in rows)
    print(f"Wrote {out}; sites={len(rows)} accepted_roots={good_root} strict_cross_precision_pass={cross_pass}")
    for r in rows:
        print(r["surface"], r["run_state"], r["TOF_C6H12O_s-1"], r["theta_OOH"], r["theta_vacancy"], r["DRC_state"])


if __name__ == "__main__":
    main()
