from __future__ import annotations

import csv
import json
from decimal import Decimal
from pathlib import Path

BASE = Path(r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\multisite_conditional_replay_20261021")
OUT = Path(r"C:\Users\Administrator\Documents\codex探索\numerical_gate_audit_20261021\data")
LOCAL = OUT.parent / "saved_root_local_flux_audit.csv"
SITES = ["tico", "ti", "tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "titi", "tiw", "tita", "ticr", "tini"]
LOG10_TOL = Decimal("0.00000217147")  # log10(1 + 5e-6), six-significant-digit reproducibility target
REL_COVERAGE_TOL = Decimal("1e-4")
EQUATION_TOL = Decimal("1e-12")


def read_csv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return {r["surface"]: r for r in csv.DictReader(f)}


def dec(value: str | None) -> Decimal | None:
    if not value or value in {"", "None", "nan"}:
        return None
    return Decimal(value)


def main() -> None:
    old = read_csv(BASE / "data" / "MULTISITE_CONDITIONAL_SUMMARY.csv")
    local = read_csv(LOCAL)
    reevaluate_rows = []
    status_rows = []
    for site in SITES:
        prior = old.get(site, {})
        audit = local.get(site, {})
        if audit.get("status") != "EVALUATED_FROM_SAVED_ROOT_NO_RESOLVE":
            status = "NUMERIC_FAIL" if prior.get("run_state", "").startswith("NO_ACCEPTED_ROOT") or prior.get("run_state", "").startswith("NUMERICAL_FAILED") else "NOT_EVALUATED"
            status_rows.append({"surface": site, "old_status": prior.get("run_state", "NOT_IN_SUMMARY"), "old_cross_precision_gate": prior.get("cross_precision_gate", ""), "new_status": status, "new_equation_gate": "NOT_EVALUATED", "cross_precision_log10_delta": prior.get("cross_precision_dlog10TOF", ""), "coverage_relative_error_upper_bound_major_and_OOH": "", "reason": "No archived accepted-root JSON for independent CatMAP rate reevaluation."})
            continue

        max_local = dec(audit.get("max_species_normalized_residual"))
        max_abs_cov = None
        cross_path = BASE / "sites" / site / "cross_precision_replay.json"
        cross = json.loads(cross_path.read_text(encoding="utf-8")) if cross_path.exists() else {}
        max_abs_cov = dec(cross.get("max_absolute_coverage_difference"))
        result_path = BASE / "sites" / site / f"constrained_{site}_result.json"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        cov = {k: Decimal(v) for k, v in result["coverage"].items()}
        relevant = [v for v in cov.values() if v >= Decimal("1e-8")]
        relevant.append(cov.get("OOH_s", Decimal(0)))
        min_relevant = min((v for v in relevant if v > 0), default=None)
        cov_bound = (max_abs_cov / min_relevant) if max_abs_cov is not None and min_relevant else None
        dlog = dec(cross.get("delta_log10_TOF_vs_high_precision"))
        cross_scaled = dec(cross.get("scaled_residual_inf"))
        i1, i2, site_sum = dec(cross.get("I1")), dec(cross.get("I2")), dec(cross.get("site_sum"))

        equation_ok = max_local is not None and max_local <= EQUATION_TOL and cross_scaled is not None and cross_scaled <= EQUATION_TOL
        conservation_ok = i1 is not None and i2 is not None and site_sum is not None and abs(i1) <= Decimal("1e-70") and abs(i2) <= Decimal("1e-70") and abs(site_sum - 1) <= Decimal("1e-70")
        coverage_ok = cov_bound is not None and cov_bound <= REL_COVERAGE_TOL
        tof_ok = dlog is not None and abs(dlog) <= LOG10_TOL
        if equation_ok and conservation_ok and coverage_ok and tof_ok:
            status = "NUMERIC_PASS"
        elif equation_ok and conservation_ok:
            status = "NUMERIC_PARTIAL"
        else:
            status = "NUMERIC_PARTIAL"
        reason = "All gates satisfied." if status == "NUMERIC_PASS" else "; ".join(k for k, ok in [("cross-precision log10(TOF) exceeds target", tof_ok), ("coverage upper-bound exceeds target", coverage_ok), ("equation residual criterion", equation_ok), ("inventory/site-balance criterion", conservation_ok)] if not ok)

        status_rows.append({
            "surface": site,
            "old_status": prior.get("run_state", "NOT_IN_SUMMARY"),
            "old_cross_precision_gate": prior.get("cross_precision_gate", ""),
            "new_status": status,
            "new_equation_gate": "PASS" if equation_ok else "FAIL",
            "cross_precision_log10_delta": str(dlog) if dlog is not None else "",
            "coverage_relative_error_upper_bound_major_and_OOH": str(cov_bound) if cov_bound is not None else "",
            "reason": reason,
        })

        if site in {"ti", "tihf", "tizr"}:
            theta_ooh = cov["OOH_s"]
            reevaluate_rows.append({
                "surface": site,
                "archived_root_scaled_residual_inf": result.get("scaled_residual_inf", ""),
                "archived_root_relative_residual_to_max_reaction_rate": result.get("relative_residual_to_max_reaction_rate", ""),
                "old_cross_relative_residual_to_max_net": cross.get("relative_residual_to_max_net_rate", ""),
                "old_cross_scaled_residual": cross.get("scaled_residual_inf", ""),
                "old_cross_abs_residual": cross.get("max_abs_surface_residual_independent", ""),
                "reevaluated_max_species_gross_normalized_residual": audit.get("max_species_normalized_residual", ""),
                "reevaluated_max_abs_species_residual": audit.get("max_abs_residual_independent", ""),
                "site_sum_abs_error": audit.get("site_sum_error", ""),
                "I1_abs_error": audit.get("I1", ""),
                "I2_abs_error": audit.get("I2", ""),
                "TOF_s-1": audit.get("TOF_C6H12O", ""),
                "log10_TOF_high_precision": prior.get("log10_TOF", ""),
                "cross_precision_delta_log10_TOF": str(dlog) if dlog is not None else "",
                "max_abs_coverage_delta": str(max_abs_cov) if max_abs_cov is not None else "",
                "OOH_coverage": str(theta_ooh),
                "conservative_relative_coverage_delta_bound_OOH": str(max_abs_cov / theta_ooh) if max_abs_cov is not None else "",
                "major_and_OOH_relative_coverage_delta_upper_bound": str(cov_bound) if cov_bound is not None else "",
                "new_status": status,
                "interpretation": reason,
            })

    for filename, rows in [("THREE_SITE_REEVALUATION.csv", reevaluate_rows), ("OLD_VS_NEW_QUALITY_STATUS.csv", status_rows)]:
        if not rows:
            raise RuntimeError(f"No output rows for {filename}")
        with (OUT / filename).open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    print(f"wrote {OUT / 'THREE_SITE_REEVALUATION.csv'} rows={len(reevaluate_rows)}")
    print(f"wrote {OUT / 'OLD_VS_NEW_QUALITY_STATUS.csv'} rows={len(status_rows)}")
    for r in status_rows:
        print(r["surface"], r["old_cross_precision_gate"], r["new_status"], r["reason"])


if __name__ == "__main__":
    main()
