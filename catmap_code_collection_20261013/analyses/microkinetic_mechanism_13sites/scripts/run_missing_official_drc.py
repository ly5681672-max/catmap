"""Run CatMAP 0.3.1 official rate_control only for missing strict-baseline sites.

Existing steady-state/coverage/rate CSVs and five existing DRC caches are
read-only. This script creates isolated per-site DRC runs under this analysis.
It does not alter the source project, energies, mechanism, or reservoir values.
"""

from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


OUT = Path(__file__).resolve().parents[1]
WORKSPACE = OUT.parents[1]
STRICT = WORKSPACE / "catmap_official_run" / "results" / "official_strict_tolerance"
EXISTING_DRC = WORKSPACE / "catmap_official_run" / "results" / "activity_analysis_round1" / "official_drc"
TARGET = OUT / "data" / "official_drc_runs_supplemental"
SITES = ["ticr", "tihf", "timn", "timo", "tita", "titi", "tiv", "tiw"]
OUTPUTS = "output_variables = ['coverage', 'rate', 'turnover_frequency', 'directional_rates']"

WORKER = r'''import csv
import json
from catmap import ReactionModel, __version__

if __version__ != "0.3.1":
    raise RuntimeError("Expected CatMAP 0.3.1, got " + str(__version__))
model = ReactionModel(setup_file="alkene_epoxidation.mkm")
model.run(recalculate=True)
labels = model.output_labels["rate_control"]
points = model.rate_control_map
if len(points) != 1:
    raise RuntimeError("Expected exactly one fixed-condition descriptor point")
matrix = points[0][1]
if len(matrix) != len(labels[0]) or any(len(row) != len(labels[1]) for row in matrix):
    raise RuntimeError("Official DRC dimensions do not match output_labels")
with open("rate_control.csv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["gas_species", "energy_parameter", "degree_of_rate_control"])
    for gas, row in zip(labels[0], matrix):
        for parameter, value in zip(labels[1], row):
            writer.writerow([gas, parameter, repr(value)])
metadata = {
    "catmap_version": __version__,
    "surface": model.surface_names[0],
    "temperature_K": float(model.temperature),
    "pressure_mode": model.pressure_mode,
    "solver_tolerance": str(model.solver.tolerance),
    "output_labels": labels,
    "product_gas": "C6H12O_g",
    "product_DRC": matrix[labels[0].index("C6H12O_g")],
    "method": "CatMAP official get_rate_control; numerical Jacobian of turnover_frequency",
}
with open("drc_run_metadata.json", "w", encoding="utf-8") as handle:
    json.dump(metadata, handle, indent=2, ensure_ascii=False, default=str)
print(json.dumps(metadata, ensure_ascii=False, default=str))
'''


def replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"Expected exactly one {label} assignment; found {count}")
    return updated


def prepare_site(site: str) -> dict:
    source = STRICT / site
    run_dir = TARGET / site
    run_dir.mkdir(parents=True, exist_ok=True)
    setup = (source / "alkene_epoxidation.mkm").read_text(encoding="utf-8")
    if setup.count(OUTPUTS) != 1:
        raise RuntimeError(f"{site}: source output_variables anchor is not unique")
    setup = setup.replace(OUTPUTS, "output_variables = ['rate_control']", 1)
    setup = replace_once(setup, r"^data_file\s*=\s*['\"][^'\"]+['\"]$",
                         f"data_file = 'catmap_{site}_mechanism_drc.pkl'", "data_file")
    setup = replace_once(setup, r"^log_file\s*=\s*['\"][^'\"]+['\"]$",
                         f"log_file = 'catmap_{site}_mechanism_drc.log'", "log_file")
    (run_dir / "alkene_epoxidation.mkm").write_text(setup, encoding="utf-8")
    shutil.copy2(source / "energies.txt", run_dir / "energies.txt")
    (run_dir / "run_drc.py").write_text(WORKER, encoding="utf-8")
    proc = subprocess.run([sys.executable, "-X", "utf8", "run_drc.py"], cwd=run_dir,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")
    log = proc.stdout + ("\n" if proc.stdout and proc.stderr else "") + proc.stderr
    (run_dir / "catmap_console.log").write_text(log, encoding="utf-8")
    status = "SUCCESS" if proc.returncode == 0 and (run_dir / "rate_control.csv").is_file() else "FAILED"
    return {"surface": site, "status": status, "returncode": proc.returncode,
            "run_dir": str(run_dir), "log_file": str(run_dir / "catmap_console.log"),
            "source_setup": str(source / "alkene_epoxidation.mkm"),
            "source_energies": str(source / "energies.txt"),
            "rate_control_csv": str(run_dir / "rate_control.csv") if status == "SUCCESS" else ""}


def main() -> None:
    if not STRICT.is_dir() or not EXISTING_DRC.is_dir():
        raise FileNotFoundError("Strict baseline or existing DRC cache directory is missing")
    missing = [site for site in SITES if not (EXISTING_DRC / site / "rate_control.csv").is_file()]
    if len(missing) != len(SITES):
        raise RuntimeError(f"Expected only eight missing sites; found {missing}")
    TARGET.mkdir(parents=True, exist_ok=True)
    results = [prepare_site(site) for site in missing]
    with (TARGET / "supplemental_drc_run_index.json").open("w", encoding="utf-8") as handle:
        json.dump({"missing_sites": missing, "runs": results}, handle, indent=2, ensure_ascii=False)
    print(json.dumps(results, indent=2, ensure_ascii=False))
    failed = [item for item in results if item["status"] != "SUCCESS"]
    if failed:
        raise SystemExit(f"Official CatMAP DRC failed for {len(failed)} sites; see per-site logs")


if __name__ == "__main__":
    main()
