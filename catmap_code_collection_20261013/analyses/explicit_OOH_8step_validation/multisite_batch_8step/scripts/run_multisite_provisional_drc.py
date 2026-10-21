"""Finite-difference sensitivities for explicit NEB TSs on accepted branches.

These are conditional CatMAP candidate DRC estimates, not formal physical DRCs.
Each perturbation re-solves the same stoichiometry-derived zero-inventory class.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
BATCH = HERE.parent
WORKFLOW = BATCH.parent / "multisite_candidate_workflow"
SOLVER = WORKFLOW / "scripts" / "run_site_candidate.py"
KB = mp.mpf("8.617332478e-5")
TEMP = mp.mpf("333.15")
DELTA = mp.mpf("0.005")


def rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_rows(path: Path, data: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in data for k in r))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(data)


def candidate_folder(surface: str) -> Path:
    return BATCH / "runs" / f"{surface}_candidate"


def execute_variant(surface: str, species: str, delta: mp.mpf, precision: int, timeout: int) -> dict:
    base = candidate_folder(surface)
    tag = f"{species.lower()}_{'plus' if delta > 0 else 'minus'}"
    out = BATCH / "runs" / "sensitivity" / surface / tag
    out.mkdir(parents=True, exist_ok=True)
    source_energy = base / f"candidate_{surface}_energies.txt"
    data = rows(source_energy, "\t")
    selected = [r for r in data if r["surface_name"] == surface and r["species_name"] == species]
    if len(selected) != 1:
        return {"surface_name": surface, "species": species, "delta_eV": str(delta), "status": "INPUT_ENERGY_NOT_UNIQUE"}
    selected[0]["formation_energy"] = mp.nstr(mp.mpf(selected[0]["formation_energy"]) + delta, 30)
    fields = list(data[0])
    with (out / "candidate_energies.txt").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(data)
    model_text = (base / f"candidate_{surface}.mkm").read_text(encoding="utf-8")
    model_text, n_energy = re.subn(r"input_file\s*=\s*'[^']+'", "input_file = 'candidate_energies.txt'", model_text)
    model_text, n_pickle = re.subn(r"data_file\s*=\s*'[^']+'", "data_file = 'candidate_sensitivity.pkl'", model_text)
    if (n_energy, n_pickle) != (1, 1):
        return {"surface_name": surface, "species": species, "delta_eV": str(delta), "status": "MODEL_ANCHOR_NOT_UNIQUE"}
    model_path = out / "candidate_sensitivity.mkm"
    model_path.write_text(model_text, encoding="utf-8")
    env = os.environ.copy(); env["PYTHONIOENCODING"] = "utf-8"
    log = BATCH / "logs" / f"drc_{surface}_{tag}.log"
    try:
        proc = subprocess.run([sys.executable, str(SOLVER), "--model", str(model_path), "--precision", str(precision)],
                              cwd=out, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(f"TIMEOUT {timeout}s\n{exc.stdout!r}\n{exc.stderr!r}\n", encoding="utf-8")
        return {"surface_name": surface, "species": species, "delta_eV": str(delta), "status": "TIMEOUT", "log": str(log)}
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(proc.stdout + proc.stderr, encoding="utf-8")
    gate_path = out / "numerical_gate.json"
    if proc.returncode or not gate_path.is_file():
        return {"surface_name": surface, "species": species, "delta_eV": str(delta), "status": "NO_ACCEPTED_ROOT",
                "exit_code": proc.returncode, "log": str(log)}
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    tof = next(r["turnover_frequency"] for r in rows(out / "candidate_tof.csv") if r["gas_species"] == "C6H12O_g")
    base_cov = {r["species"]: mp.mpf(r["coverage"]) for r in rows(base / "candidate_steady_state_coverages.csv")}
    test_cov = {r["species"]: mp.mpf(r["coverage"]) for r in rows(out / "candidate_steady_state_coverages.csv")}
    branch_dist = max(abs(test_cov.get(k, mp.mpf(0)) - v) for k, v in base_cov.items())
    return {"surface_name": surface, "species": species, "delta_eV": str(delta), "status": "NUMERICAL_PASS",
            "TOF": tof, "relative_residual": gate.get("max_relative_surface_residual", ""),
            "branch_max_abs_coverage_change": mp.nstr(branch_dist, 16),
            "branch_status": "SAME_NEARBY_BRANCH" if branch_dist < mp.mpf("0.05") else "BRANCH_UNCERTAIN",
            "inventory_error": gate.get("inventory_error", ""), "site_balance_error": gate.get("site_balance_error", ""),
            "log": str(log), "run_path": str(out)}


def main() -> None:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--precision", type=int, default=120)
    ap.add_argument("--timeout-seconds", type=int, default=180)
    ap.add_argument("--delta-eV", default="0.005")
    args = ap.parse_args()
    mp.mp.dps = max(80, args.precision)
    status_rows = rows(BATCH / "MULTISITE_CALCULATION_STATUS.csv")
    targets = [r["surface_name"] for r in status_rows if r["calculation_status"] == "COMPUTED_PROVISIONAL"]
    records, summary = [], []
    for surface in targets:
        base = candidate_folder(surface)
        for species in ("NEB1TS", "NEB2TS"):
            variants = [execute_variant(surface, species, -mp.mpf(args.delta_eV), args.precision, args.timeout_seconds),
                        execute_variant(surface, species, mp.mpf(args.delta_eV), args.precision, args.timeout_seconds)]
            records.extend(variants)
            good = [x for x in variants if x.get("status") == "NUMERICAL_PASS"]
            if len(good) == 2:
                minus = next(x for x in good if mp.mpf(x["delta_eV"]) < 0)
                plus = next(x for x in good if mp.mpf(x["delta_eV"]) > 0)
                base_tof = mp.mpf(next(r["turnover_frequency"] for r in rows(base / "candidate_tof.csv") if r["gas_species"] == "C6H12O_g"))
                derivative = (mp.log(mp.mpf(plus["TOF"])) - mp.log(mp.mpf(minus["TOF"]))) / (2 * mp.mpf(args.delta_eV))
                drc = -KB * TEMP * derivative
                kind = "CONDITIONAL_CatMAP_TS_DRC_FINITE_DIFFERENCE"
                quality = ";".join(x["branch_status"] for x in good)
            else:
                base_tof, drc, kind, quality = "", "", "FAILED_OR_INCOMPLETE", "see variant logs"
            energy_rows = rows(base / f"candidate_{surface}_energies.txt", "\t")
            emap = {r["species_name"]: mp.mpf(r["formation_energy"]) for r in energy_rows}
            path_name = species.replace("TS", "")
            prefix = "NEB1" if species == "NEB1TS" else "NEB2"
            below_endpoint = emap[species] < max(emap[prefix + "IS"], emap[prefix + "FS"])
            summary.append({"surface_name": surface, "perturbed_state": species, "kind": kind,
                            "delta_eV": args.delta_eV, "baseline_TOF_s-1": mp.nstr(base_tof, 40) if base_tof else "",
                            "TOF_minus": next((x.get("TOF", "") for x in good if mp.mpf(x["delta_eV"]) < 0), ""),
                            "TOF_plus": next((x.get("TOF", "") for x in good if mp.mpf(x["delta_eV"]) > 0), ""),
                            "conditional_DRC": mp.nstr(drc, 24) if drc != "" else "",
                            "branch_continuity": quality,
                            "CatMAP_TS_controls_endpoint_max": not below_endpoint,
                            "interpretation_note": "finite-difference candidate sensitivity only; not formal DRC; explicit TS below endpoint means CatMAP endpoint-max clipping controls" if below_endpoint else "finite-difference candidate sensitivity only; physical mapping remains provisional"})
            print(f"{surface} {species}: {kind} DRC={summary[-1]['conditional_DRC']}", flush=True)
            write_rows(BATCH / "data" / "multisite_drc_variant_runs.csv", records)
            write_rows(BATCH / "MULTISITE_DRC_PROVISIONAL.csv", summary)
    print(f"variants={len(records)}; passed={sum(r.get('status') == 'NUMERICAL_PASS' for r in records)}", flush=True)


if __name__ == "__main__":
    main()
