"""Bounded, explicitly hypothetical TS-height sensitivity for endpoint links.

The tested extra barrier span is anchored to the two known Ti-Co NEB forward
barriers (0.225994 and 0.556344 eV) and reverse barriers (0.706574 and
0.796030 eV); it is a scenario span, not evidence or a confidence interval
for any missing transition state.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


LINKS = {
    1: ("* + H2O2_g <-> H2O2*", {"H": 2, "O": 2}, 1),
    2: ("H2O2* + * -> Ha* + OOH*", {"H": 2, "O": 2}, 2),
    3: ("* + C6H12_g <-> C6H12*", {"C": 6, "H": 12}, 1),
    4: ("Ha* + OOH* + C6H12* <-> NEB1IS* + 2*", {"C": 6, "H": 14, "O": 2}, 3),
    6: ("NEB1FS* <-> NEB2IS* + C6H12O_g", {"C": 6, "H": 14, "O": 2}, 1),
    8: ("NEB2FS* <-> H2O*", {"H": 2, "O": 1}, 1),
    9: ("H2O* <-> H2O_g + *", {"H": 2, "O": 1}, 1),
}
EXTRA_BARRIERS_EV = (0.225994, 0.796030)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--base-run", type=Path, required=True)
    p.add_argument("--solver-script", type=Path, required=True)
    args = p.parse_args()
    base = args.base_run.resolve()
    root = base.parent / "endpoint_barrier_sensitivity"
    root.mkdir(parents=True, exist_ok=True)
    energy_source = base / "candidate_energies.txt"
    model_source = base / "candidate_tico.mkm"

    from catmap import ReactionModel
    import os
    cwd = Path.cwd()
    os.chdir(base)
    try:
        model = ReactionModel(setup_file=model_source.name)
        model.run()
        energies = model.scaler.get_free_energies([model.temperature, 1.0])
        reactions = model.elementary_rxns
        if len(reactions) != 9:
            raise RuntimeError("base model must have nine CatMAP reactions")
        endpoint_g = {}
        for step in LINKS:
            rxn = reactions[step - 1]
            gi = model.get_state_energy(rxn[0], energies)
            gf = model.get_state_energy(rxn[-1], energies)
            endpoint_g[step] = max(gi, gf)
    finally:
        os.chdir(cwd)

    with energy_source.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)
    header, energy_rows = rows[0], rows[1:]
    results = []
    for step, (old_expression, composition, n_sites) in LINKS.items():
        for delta in EXTRA_BARRIERS_EV:
            case = root / f"step{step}_extra_{delta:.6f}eV"
            case.mkdir(parents=True, exist_ok=True)
            shutil.copy2(energy_source, case / "candidate_energies.txt")
            ts_internal = f"Link{step}TS"
            energy = endpoint_g[step] + delta
            with (case / "candidate_energies.txt").open("a", encoding="utf-8", newline="") as f:
                f.write("\t".join(["tico", "111", ts_internal,
                                   f"{float(energy):.12f}", "[]",
                                   "HYPOTHETICAL_SENSITIVITY_ONLY_NOT_DFT"]) + "\n")
            config = model_source.read_text(encoding="utf-8")
            ts_species = ts_internal + "_s"
            new_expression = old_expression.replace("<->", f"<-> {ts_internal}* ->", 1) if "<->" in old_expression else old_expression.replace("->", f"<-> {ts_internal}* ->", 1)
            # Preserve the original reversible/irreversible semantics on the
            # reactant side; only insert a scenario TS before the same product.
            if old_expression == "H2O2* + * -> Ha* + OOH*":
                new_expression = f"H2O2* + * <-> {ts_internal}* -> Ha* + OOH*"
            config += (f"\n# Hypothetical TS-height scenario for step {step}; not DFT.\n"
                       f"rxn_expressions[{step - 1}] = {new_expression!r}\n"
                       f"species_definitions[{ts_species!r}] = {{'composition': {composition!r}, 'n_sites': {n_sites}}}\n"
                       "input_file = 'candidate_energies.txt'\n"
                       f"data_file = 'sensitivity_step{step}_{delta:.6f}.pkl'\n")
            model_path = case / "candidate_tico.mkm"
            model_path.write_text(config, encoding="utf-8")
            proc = subprocess.run([sys.executable, str(args.solver_script.resolve()),
                                   "--model", str(model_path), "--precision", "160",
                                   "--seed-coverages", str(base / "candidate_steady_state_coverages.csv")],
                                  cwd=case, capture_output=True, text=True, encoding="utf-8", errors="replace")
            (case / "solver_stdout.log").write_text(proc.stdout, encoding="utf-8")
            (case / "solver_stderr.log").write_text(proc.stderr, encoding="utf-8")
            result = {"step": step, "extra_barrier_eV": f"{delta:.6f}",
                      "hypothetical_TS_eV": f"{float(energy):.12f}",
                      "scenario_only": "TRUE_NOT_DFT",
                      "process_exit_code": proc.returncode,
                      "solve_status": "FAILED"}
            attempts_path = case / "candidate_solver_attempts.json"
            tof_path = case / "candidate_tof.csv"
            if proc.returncode == 0 and attempts_path.exists():
                attempts = json.loads(attempts_path.read_text(encoding="utf-8"))
                result["solve_status"] = attempts["status"]
                accepted = attempts.get("accepted_roots", [])
                if accepted:
                    result["product_TOF"] = accepted[0]["product_TOF"]
                result["accepted_root_count"] = len(accepted)
            result["case_directory"] = str(case)
            result["log_path"] = str(case / "solver_stdout.log")
            results.append(result)
            print(json.dumps(result, ensure_ascii=False))

    output = root / "endpoint_barrier_sensitivity.csv"
    fields = list(dict.fromkeys(k for row in results for k in row))
    with output.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    print(f"summary={output}")


if __name__ == "__main__":
    main()
