"""Re-solve constrained CatMAP candidate under small energy perturbations.

TS perturbations are conditional numerical DRC estimates. NEB endpoint-state
perturbations are reported separately as endpoint-model sensitivities, not as
missing-transition-state DRCs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "runs" / "tife_candidate"
SOLVER = ROOT / "scripts" / "run_site_candidate.py"
KB_EV_K = 8.617333262145e-5


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def run_variant(species: str, delta: float, precision: int) -> dict:
    mp.mp.dps = 100
    tag = f"{species.replace('_s','').lower()}_{'plus' if delta > 0 else 'minus'}"
    out = ROOT / "runs" / "tife_sensitivity" / tag
    out.mkdir(parents=True, exist_ok=True)
    source = BASE / "candidate_tife_energies.txt"
    with source.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
        fields = list(rows[0])
    found = 0
    for row in rows:
        if row["surface_name"] == "tife" and row["species_name"] == species:
            row["formation_energy"] = f"{float(row['formation_energy']) + delta:.12f}"
            found += 1
    if found != 1:
        raise RuntimeError(f"expected one energy row for {species}, found {found}")
    energy_path = out / "candidate_energies.txt"
    write_tsv(energy_path, rows, fields)
    model_text = (BASE / "candidate_tife.mkm").read_text(encoding="utf-8")
    model_text = model_text.replace("candidate_tife_energies.txt", energy_path.name)
    model_text = model_text.replace("candidate_tife.pkl", "candidate_sensitivity.pkl")
    model_path = out / "candidate_sensitivity.mkm"
    model_path.write_text(model_text, encoding="utf-8")
    log_path = out / "solver.log"
    env = os.environ.copy(); env["PYTHONIOENCODING"] = "utf-8"
    with log_path.open("w", encoding="utf-8") as log:
        proc = subprocess.run([sys.executable, str(SOLVER), "--model", str(model_path), "--precision", str(precision)],
                              cwd=out, env=env, text=True, stdout=log, stderr=subprocess.STDOUT)
    gate_path = out / "numerical_gate.json"
    if proc.returncode != 0 or not gate_path.is_file():
        return {"species": species, "delta_eV": delta, "status": "FAILED", "exit_code": proc.returncode,
                "log": str(log_path)}
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    with (out / "candidate_tof.csv").open(encoding="utf-8-sig", newline="") as stream:
        tof_rows = list(csv.DictReader(stream))
    product = next(mp.mpf(r["turnover_frequency"]) for r in tof_rows if r["gas_species"] == "C6H12O_g")
    cov = {}
    with (out / "candidate_steady_state_coverages.csv").open(encoding="utf-8-sig", newline="") as stream:
        cov = {r["species"]: float(r["coverage"]) for r in csv.DictReader(stream)}
    return {"species": species, "delta_eV": delta, "status": "NUMERICAL_PASS_PROVISIONAL_PHYSICS",
            "TOF": mp.nstr(product, 80), "log10_TOF": mp.nstr(mp.log10(product), 40) if product > 0 else "",
            "coverages": cov, "gate": gate, "log": str(log_path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=220)
    parser.add_argument("--delta-ts", type=float, default=0.005)
    parser.add_argument("--delta-endpoint", type=float, default=0.005)
    args = parser.parse_args()
    base_gate = json.loads((BASE / "numerical_gate.json").read_text(encoding="utf-8"))
    mp.mp.dps = 100
    baseline_tof = mp.mpf(base_gate["candidate_T_formula_s-1"])
    baseline_cov = {}
    with (BASE / "candidate_steady_state_coverages.csv").open(encoding="utf-8-sig", newline="") as stream:
        baseline_cov = {r["species"]: float(r["coverage"]) for r in csv.DictReader(stream)}
    results = []
    for species in ("NEB1TS", "NEB2TS", "NEB1FS", "NEB2FS"):
        delta = args.delta_ts if species.endswith("TS") else args.delta_endpoint
        for sign in (-1, 1):
            result = run_variant(species, sign * delta, args.precision)
            if result["status"] == "NUMERICAL_PASS_PROVISIONAL_PHYSICS":
                result["max_abs_coverage_change_vs_baseline"] = max(abs(result["coverages"].get(k, 0.0) - v)
                                                                       for k, v in baseline_cov.items())
                result["relative_branch_diagnostic"] = "SAME_NEARBY_BRANCH" if result["max_abs_coverage_change_vs_baseline"] < 0.05 else "BRANCH_CONTINUITY_UNCERTAIN"
            results.append(result)
    summary = []
    for species in ("NEB1TS", "NEB2TS", "NEB1FS", "NEB2FS"):
        pair = [r for r in results if r["species"] == species]
        if len(pair) == 2 and all(r["status"] == "NUMERICAL_PASS_PROVISIONAL_PHYSICS" for r in pair):
            minus = next(r for r in pair if r["delta_eV"] < 0); plus = next(r for r in pair if r["delta_eV"] > 0)
            derivative = (mp.log(mp.mpf(plus["TOF"])) - mp.log(mp.mpf(minus["TOF"]))) / (2 * abs(plus["delta_eV"]))
            kind = "CONDITIONAL_CatMAP_DRC_FINITE_DIFFERENCE" if species.endswith("TS") else "ENDPOINT_ENERGY_SENSITIVITY_NOT_DRC"
            val = -mp.mpf(str(KB_EV_K)) * 333.15 * derivative if species.endswith("TS") else derivative
            summary.append({"perturbed_species": species, "kind": kind, "delta_eV": abs(plus["delta_eV"]),
                            "baseline_TOF_s-1": mp.nstr(baseline_tof, 60), "TOF_minus": minus["TOF"], "TOF_plus": plus["TOF"],
                            "sensitivity": mp.nstr(val, 30),
                            "branch_continuity": ";".join(r["relative_branch_diagnostic"] for r in pair)})
        else:
            summary.append({"perturbed_species": species, "kind": "FAILED_OR_INCOMPLETE", "delta_eV": delta,
                            "baseline_TOF_s-1": baseline_tof, "TOF_minus": "", "TOF_plus": "",
                            "sensitivity": "", "branch_continuity": "see per-run logs"})
    out = ROOT / "data" / "tife_conditional_sensitivity.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0])); writer.writeheader(); writer.writerows(summary)
    (ROOT / "data" / "tife_sensitivity_run_manifest.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(out), "variants": len(results), "passed": sum(r["status"] == "NUMERICAL_PASS_PROVISIONAL_PHYSICS" for r in results)}, indent=2))


if __name__ == "__main__":
    main()
