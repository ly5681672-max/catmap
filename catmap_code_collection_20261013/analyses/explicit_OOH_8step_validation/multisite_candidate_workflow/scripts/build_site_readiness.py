"""Build a provenance-based readiness roster; this does not run CatMAP."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAP = ROOT / "inputs" / "source_snapshot"
OUT = ROOT / "data"
TARGETS = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "tico", "titi", "tiw", "tita", "ticr", "ti", "tini"]
REQUIRED = {"H2O2", "Ha", "OOH", "C6H12", "H2O", "C6H12O"}


def rows(name: str, delimiter: str = ",") -> list[dict[str, str]]:
    with (SNAP / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    formation = rows("DFT_formation_catmap.txt", "\t")
    neb_rows = {row["surface"]: row for row in rows("DFT_TS_relative.csv")}
    by_surface: dict[str, set[str]] = {}
    for row in formation:
        if row["site_name"] != "gas" and row["surface_name"] not in ("", "None"):
            by_surface.setdefault(row["surface_name"], set()).add(row["species_name"])
    assessment = []
    for surface in TARGETS:
        have = by_surface.get(surface, set())
        neb = neb_rows.get(surface)
        neb_fields = ["G_IS1", "G_TS1", "G_FS1_total", "G_IS2", "G_TS2", "G_FS2_adjusted"]
        missing_single = sorted(REQUIRED - have)
        missing_neb = [] if neb else neb_fields
        if neb:
            missing_neb = [field for field in neb_fields if not neb.get(field, "").strip()]
        signed = []
        if neb and not missing_neb:
            for prefix in ("dG_TS1_minus_IS1", "dG_TS2_minus_IS2"):
                signed.append(float(neb[prefix]))
        status = "SOURCE_ROWS_PRESENT_PHYSICS_MAPPING_UNVERIFIED" if not missing_single and not missing_neb else "INCOMPLETE_SOURCE_ROWS"
        assessment.append({"surface_name": surface, "independent_inputs_complete": not missing_single,
                           "missing_independent_species": ";".join(missing_single),
                           "neb1_neb2_all_six_values_present": bool(neb is not None and not missing_neb),
                           "missing_neb_fields": ";".join(missing_neb),
                           "NEB1_TS_minus_IS_signed_eV": signed[0] if len(signed) == 2 else "",
                           "NEB2_TS_minus_IS_signed_eV": signed[1] if len(signed) == 2 else "",
                           "calculation_status": "NOT_RUN_EXCEPT_SEPARATE_TIFE_CANDIDATE",
                           "physical_gate": "PROVISIONAL_UNVERIFIED_COMPOSITE_SITE_MAPPING",
                           "readiness": status})
    with (OUT / "site_readiness.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(assessment[0]))
        writer.writeheader(); writer.writerows(assessment)
    dataset_fields = ["surface_name", "G_OOH", "G_H2O2", "G_NEB1_IS", "G_NEB1_TS", "G_NEB1_FS",
                      "G_NEB2_IS", "G_NEB2_TS", "G_NEB2_FS", "theta_OOH", "theta_Ha", "theta_vacancy",
                      "TOF_C6H12O_s-1", "log10_TOF", "net_flux_step_1", "net_flux_step_2", "net_flux_step_3",
                      "net_flux_step_4", "net_flux_step_5", "net_flux_step_6", "net_flux_step_7", "net_flux_step_8",
                      "net_flux_step_9", "DRC_NEB1_TS", "DRC_NEB2_TS", "conditional_link_sensitivity",
                      "physical_gate", "numerical_gate", "source_hash_manifest"]
    with (OUT / "multisite_ooh_tof_drc_schema.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        csv.writer(stream).writerow(dataset_fields)
    print(f"target_sites={len(TARGETS)}; complete_source_rows={sum(r['readiness'].startswith('SOURCE_ROWS') for r in assessment)}; TiNb_additional={('tinb' in by_surface)}")


if __name__ == "__main__":
    main()
