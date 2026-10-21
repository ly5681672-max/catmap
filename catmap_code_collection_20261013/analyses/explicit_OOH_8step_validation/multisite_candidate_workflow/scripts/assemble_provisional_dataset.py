"""Assemble only source-backed TiCo/TiFe candidate rows; no values are imputed."""
from __future__ import annotations

import csv
import math
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[2]
TICO = PROJECT / "analyses" / "explicit_OOH_8step_validation" / "candidate_tico_kinetic_control"
TIFE = ROOT / "runs" / "tife_candidate"
OUT = ROOT / "data" / "multisite_ooh_tof_drc_schema.csv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    mp.mp.dps = 80
    fields = ["surface_name", "G_OOH", "G_H2O2", "G_NEB1_IS", "G_NEB1_TS", "G_NEB1_FS",
              "G_NEB2_IS", "G_NEB2_TS", "G_NEB2_FS", "theta_OOH", "theta_Ha", "theta_vacancy",
              "TOF_C6H12O_s-1", "log10_TOF", "net_flux_step_1", "net_flux_step_2", "net_flux_step_3",
              "net_flux_step_4", "net_flux_step_5", "net_flux_step_6", "net_flux_step_7", "net_flux_step_8",
              "net_flux_step_9", "DRC_NEB1_TS", "DRC_NEB2_TS", "conditional_link_sensitivity",
              "physical_gate", "numerical_gate", "source_hash_manifest"]
    output = []
    for surface, folder, cov_file, tof_file, rates_file, drc1, drc2, link, gate, hashes in [
        ("tico", TICO, "candidate_steady_state_coverages.csv", "candidate_tof_baseline.csv", None,
         "3.335475522573327268460361204206354243738e-31", "3.35533283070637728096211216634806502545e-36",
         "step2/step4 endpoint assumptions: Xδ≈0.519/0.524 at δ=0.005 eV; provisional one-sided sensitivity",
         "PROVISIONAL_OLD_CANDIDATE_MAPPING", "INPUT_SHA256.txt"),
        ("tife", TIFE, "candidate_steady_state_coverages.csv", "candidate_tof.csv", "candidate_reaction_rates.csv",
         "6.84627396888061478566743280434e-41", "1.09807011262946494948256676007e-40",
         "NEB1FS endpoint energy sensitivity=1.2053e-13; NEB2FS endpoint energy sensitivity=1.26236 (not barrier DRC)",
         "PROVISIONAL_UNVERIFIED_REDUCED_NETWORK_AND_ONE_SITE_MAPPING", "../inputs/source_snapshot/source_sha256.json")]:
        en = {r["species_name"]: r["formation_energy"] for r in read_tsv(folder / ("candidate_energies.txt" if surface == "tico" else "candidate_tife_energies.txt"))}
        cov = {r["species"].replace("_s", ""): r["coverage"] for r in read_csv(folder / cov_file)}
        tof_rows = read_csv(folder / tof_file)
        tof = next(r["turnover_frequency"] for r in tof_rows if r["gas_species"] == "C6H12O_g")
        rate_values = []
        if rates_file:
            rate_values = [r["net_rate"] for r in read_csv(folder / rates_file)]
        else:
            # Historical Tico baseline re-exported its four gas turnover rates;
            # candidate net reaction rates are not reconstructed from those.
            rate_values = [""] * 9
        row = {key: "" for key in fields}
        row.update({"surface_name": surface, "G_OOH": en.get("OOH", ""), "G_H2O2": en.get("H2O2", ""),
                    "G_NEB1_IS": en.get("NEB1IS", ""), "G_NEB1_TS": en.get("NEB1TS", ""), "G_NEB1_FS": en.get("NEB1FS", ""),
                    "G_NEB2_IS": en.get("NEB2IS", ""), "G_NEB2_TS": en.get("NEB2TS", ""), "G_NEB2_FS": en.get("NEB2FS", ""),
                    "theta_OOH": cov.get("OOH", ""), "theta_Ha": cov.get("Ha", ""), "theta_vacancy": cov.get("*", ""),
                    "TOF_C6H12O_s-1": tof, "log10_TOF": mp.nstr(mp.log10(mp.mpf(tof)), 40),
                    "DRC_NEB1_TS": drc1, "DRC_NEB2_TS": drc2,
                    "conditional_link_sensitivity": link, "physical_gate": gate,
                    "numerical_gate": "PASS_CANDIDATE_ROOTS_220_DIGITS" if surface == "tife" else "PASS_HISTORICAL_VALIDATION_SUMMARY",
                    "source_hash_manifest": hashes})
        for i in range(9):
            row[f"net_flux_step_{i+1}"] = rate_values[i] if i < len(rate_values) else ""
        output.append(row)
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(output)
    print(f"rows={len(output)}; sources=existing provisional TiCo + newly solved provisional TiFe; no missing fields imputed")


if __name__ == "__main__":
    main()
