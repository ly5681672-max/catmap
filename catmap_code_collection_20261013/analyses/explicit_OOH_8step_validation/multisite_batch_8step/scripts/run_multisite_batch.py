"""Run the remaining explicit-OOH candidate sites with the existing CatMAP workflow.

All outputs are provisional endpoint-model results. The script reuses the source
snapshot and solver from multisite_candidate_workflow and never edits scientific
source tables or the CatMAP installation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BATCH = HERE.parent
WORKFLOW = BATCH.parent / "multisite_candidate_workflow"
SNAP = WORKFLOW / "inputs" / "source_snapshot"
PREPARE = WORKFLOW / "scripts" / "prepare_site_candidate.py"
SOLVER = WORKFLOW / "scripts" / "run_site_candidate.py"
TICO = BATCH.parent / "candidate_tico_kinetic_control"
TARGETS = ["ti", "timn", "tihf", "tire", "timo", "tiv", "tizr", "tico",
           "titi", "tiw", "tita", "ticr", "tife", "tini"]
ALREADY = {"tico", "tife"}
NEGATIVE_TS = {"tiw": "NEB1", "tita": "NEB1", "timo": "NEB2", "tiv": "NEB2", "ticr": "NEB2"}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def directional_from_saved_state(model_path: Path, coverage_map: dict[str, str]) -> dict[str, str]:
    """Call CatMAP's directional-rate API at a stored accepted coverage vector.

    CatMAP 0.3.1 returns these directional fluxes as Python floats; preserve
    that API precision limitation in the output labels and keep the solver's
    high-precision net flux separately.
    """
    old_cwd = Path.cwd()
    try:
        os.chdir(model_path.parent)
        from catmap import ReactionModel
        model = ReactionModel(setup_file=model_path.name)
        model.run()
        solver = model.solver
        solver.compile()
        cv = [float(coverage_map[name]) for name in solver.adsorbate_names]
        params = list(model.scaler.get_rxn_parameters([333.15, 1.0]))
        solver._coverage = cv
        solver.get_rate_constants(params, cv)
        values = list(solver.get_directional_rates(params))
        nrxn = len(solver.elementary_rxns)
        out = {f"forward_flux_step_{i+1}": repr(values[i]) for i in range(nrxn)}
        out.update({f"reverse_flux_step_{i+1}": repr(values[nrxn+i]) for i in range(nrxn)})
        return out
    finally:
        os.chdir(old_cwd)


def tico_reused() -> dict:
    gate = json.loads((TICO / "validation_summary.json").read_text(encoding="utf-8"))
    cov = {r["species"]: r["coverage"] for r in read_csv(TICO / "candidate_steady_state_coverages.csv")}
    tof_rows = read_csv(TICO / "candidate_tof_baseline.csv")
    tof = next(r["turnover_frequency"] for r in tof_rows if r.get("gas_species") == "C6H12O_g")
    flux_rows = read_csv(TICO / "FORWARD_REVERSE_RATE_AUDIT.csv")
    flux = {f"net_flux_step_{i+1}": r["net_flux_high_precision"] for i, r in enumerate(flux_rows)}
    directional = {f"forward_flux_step_{i+1}": r["official_directional_api_forward_float"] for i, r in enumerate(flux_rows)}
    directional.update({f"reverse_flux_step_{i+1}": r["official_directional_api_reverse_float"] for i, r in enumerate(flux_rows)})
    en_rows = read_csv(TICO / "candidate_energies.txt", "\t")
    energies = {r["species_name"]: r["formation_energy"] for r in en_rows}
    signed = {"NEB1_TS_minus_IS_eV": float(energies["NEB1TS"]) - float(energies["NEB1IS"]),
              "NEB1_TS_minus_FS_eV": float(energies["NEB1TS"]) - float(energies["NEB1FS"]),
              "NEB2_TS_minus_IS_eV": float(energies["NEB2TS"]) - float(energies["NEB2IS"]),
              "NEB2_TS_minus_FS_eV": float(energies["NEB2TS"]) - float(energies["NEB2FS"])}
    source_hash = sha(TICO / "INPUT_SHA256.txt") + ":" + sha(TICO / "candidate_tico.mkm") + ":" + sha(TICO / "candidate_energies.txt")
    return {"surface_name": "tico", "status": "ALREADY_COMPUTED_PROVISIONAL",
            "numerical_status": "PASS_HISTORICAL_23_CASE_VALIDATION", "catmap_version": "0.3.1",
            "direction_semantics": "REVERSIBLE_NATIVE_CATMAP", "TOF": tof,
            "coverages": cov, "fluxes": flux, "directional_fluxes": directional, "energies": energies,
            "NEB_signed": signed,
            "max_relative_residual": gate.get("maximum_relative_surface_residual", ""),
            "source_hash": source_hash, "physical_gate": "PROVISIONAL_UNVERIFIED_COMPOSITE_SITE_MAPPING",
            "run_path": str(TICO), "reused": True}


def tife_reused() -> dict:
    folder = WORKFLOW / "runs" / "tife_candidate"
    gate = json.loads((folder / "numerical_gate.json").read_text(encoding="utf-8"))
    cov_rows = read_csv(folder / "candidate_steady_state_coverages.csv")
    cov = {r["species"]: r["coverage"] for r in cov_rows}
    tof = next(r["turnover_frequency"] for r in read_csv(folder / "candidate_tof.csv")
               if r.get("gas_species") == "C6H12O_g")
    flux_rows = read_csv(folder / "candidate_reaction_rates.csv")
    flux = {f"net_flux_step_{r['reaction_index']}": r["net_rate"] for r in flux_rows}
    directional = {f"forward_flux_step_{r['reaction_index']}": r["forward_rate"] for r in flux_rows if r.get("forward_rate")}
    directional.update({f"reverse_flux_step_{r['reaction_index']}": r["reverse_rate"] for r in flux_rows if r.get("reverse_rate")})
    if len(directional) != 18:
        directional = directional_from_saved_state(folder / "candidate_tife.mkm", cov)
    en_rows = read_csv(folder / "candidate_tife_energies.txt", "\t")
    energies = {r["species_name"]: r["formation_energy"] for r in en_rows}
    return {"surface_name": "tife", "status": "ALREADY_COMPUTED_PROVISIONAL",
            "numerical_status": gate["status"], "catmap_version": gate["catmap_version"],
            "direction_semantics": "REVERSIBLE_NATIVE_CATMAP", "TOF": tof,
            "coverages": cov, "fluxes": flux, "directional_fluxes": directional, "energies": energies,
            "max_relative_residual": gate.get("max_relative_surface_residual", ""),
            "source_hash": sha(folder / "candidate_tife.mkm") + ":" + sha(folder / "candidate_tife_energies.txt") + ":" + sha(SNAP / "source_sha256.json"),
            "physical_gate": "PROVISIONAL_UNVERIFIED_COMPOSITE_SITE_MAPPING",
            "run_path": str(folder), "reused": True}


def load_batch_candidate(surface: str, folder: Path) -> dict | None:
    gate_path = folder / "numerical_gate.json"
    if not folder.exists():
        return None
    if not gate_path.is_file():
        attempt_path = folder / "candidate_solver_attempts.json"
        metadata_path = folder / "preparation_metadata.json"
        energy_path = folder / f"candidate_{surface}_energies.txt"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.is_file() else {}
        signed = metadata.get("NEB_signed_differences_eV", {})
        negative = [name for name, val in signed.items() if "_TS_minus_IS_" in name and float(val) < 0]
        attempts = json.loads(attempt_path.read_text(encoding="utf-8")) if attempt_path.is_file() else {}
        if attempts.get("status") == "NO_ACCEPTED_ROOT":
            numerical = "NO_ACCEPTED_ROOT_AT_P120"
            if (folder / "attempts" / "precision_120" / "candidate_solver_attempts.json").is_file():
                if (BATCH / "logs" / f"{surface}_p240.log").is_file():
                    numerical += ";NO_ACCEPTED_ROOT_AT_P240"
                elif (BATCH / "logs" / f"{surface}_timeout.log").is_file():
                    numerical += ";P240_TIMEOUT"
                else:
                    numerical += ";P240_STATUS_UNRESOLVED"
            err = "; ".join(f"{a.get('seed')}:relative_residual={a.get('relative_surface_residual')}" for a in attempts.get("attempts", []))
        else:
            numerical = "INTERRUPTED_NO_GATE"
            err = "run was interrupted before solver attempt diagnostics were written"
        energies = ({r["species_name"]: r["formation_energy"] for r in read_csv(energy_path, "\t")}
                    if energy_path.is_file() else {})
        source_hash = (sha(folder / f"candidate_{surface}.mkm") + ":" + sha(energy_path) + ":" + sha(SNAP / "source_sha256.json")) if energy_path.is_file() else ""
        return {"surface_name": surface, "status": "NUMERICAL_FAILED", "numerical_status": numerical,
                "error": err, "run_path": str(folder), "energies": energies,
                "NEB_signed": signed, "negative_TS_flags": ";".join(negative),
                "physical_gate": "PROVISIONAL_UNVERIFIED_COMPOSITE_SITE_MAPPING" + (":" + ";".join(negative) if negative else ""),
                "source_hash": source_hash}
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    coverages = {r["species"]: r["coverage"] for r in read_csv(folder / "candidate_steady_state_coverages.csv")}
    tof = next(r["turnover_frequency"] for r in read_csv(folder / "candidate_tof.csv")
               if r.get("gas_species") == "C6H12O_g")
    flux_rows = read_csv(folder / "candidate_reaction_rates.csv")
    flux = {f"net_flux_step_{r['reaction_index']}": r["net_rate"] for r in flux_rows}
    directional = {f"forward_flux_step_{r['reaction_index']}": r.get("forward_rate", "") for r in flux_rows}
    directional.update({f"reverse_flux_step_{r['reaction_index']}": r.get("reverse_rate", "") for r in flux_rows})
    en = {r["species_name"]: r["formation_energy"] for r in read_csv(folder / f"candidate_{surface}_energies.txt", "\t")}
    meta_path = folder / "preparation_metadata.json"
    metadata = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.is_file() else {}
    signed = metadata.get("NEB_signed_differences_eV", {})
    negative = [name for name, val in signed.items() if "_TS_minus_IS_" in name and float(val) < 0]
    directional = {f"forward_flux_step_{r['reaction_index']}": r.get("forward_rate", "") for r in flux_rows}
    directional.update({f"reverse_flux_step_{r['reaction_index']}": r.get("reverse_rate", "") for r in flux_rows})
    if len(directional) != 18 or not all(directional.values()):
        directional = directional_from_saved_state(folder / f"candidate_{surface}.mkm", coverages)
    return {"surface_name": surface, "status": "COMPUTED_PROVISIONAL",
            "numerical_status": gate["status"], "catmap_version": gate.get("catmap_version", "unknown"),
            "direction_semantics": "REVERSIBLE_NATIVE_CATMAP", "TOF": tof,
            "coverages": coverages, "fluxes": flux, "directional_fluxes": directional,
            "energies": en, "NEB_signed": signed, "negative_TS_flags": ";".join(negative),
            "max_relative_residual": gate.get("max_relative_surface_residual", ""),
            "site_balance_error": gate.get("site_balance_error", ""), "inventory_error": gate.get("inventory_error", ""),
            "source_hash": sha(folder / f"candidate_{surface}.mkm") + ":" + sha(folder / f"candidate_{surface}_energies.txt") + ":" + sha(SNAP / "source_sha256.json"),
            "physical_gate": "PROVISIONAL_UNVERIFIED_COMPOSITE_SITE_MAPPING" + (":" + ";".join(negative) if negative else ""),
            "run_path": str(folder), "reused": True}


def assemble_outputs(results: list[dict]) -> None:
    status_rows, value_rows, barrier_rows = [], [], []
    for r in results:
        surface = r["surface_name"]
        status_rows.append({"surface_name": surface, "calculation_status": r.get("status", "NOT_RUN"),
                            "numerical_status": r.get("numerical_status", ""),
                            "catmap_version": r.get("catmap_version", "0.3.1"),
                            "direction_semantics": r.get("direction_semantics", "REVERSIBLE_NATIVE_CATMAP"),
                            "physical_gate": r.get("physical_gate", "PROVISIONAL_UNVERIFIED_COMPOSITE_SITE_MAPPING"),
                            "negative_signed_TS_flags": r.get("negative_TS_flags", ""),
                            "max_relative_residual": r.get("max_relative_residual", ""),
                            "site_balance_error": r.get("site_balance_error", ""),
                            "inventory_error": r.get("inventory_error", ""),
                            "source_hash": r.get("source_hash", ""), "reused_existing_result": r.get("reused", False),
                            "run_path": r.get("run_path", ""), "error": r.get("error", "")})
        if r.get("TOF"):
            try:
                logtof = math.log10(float(r["TOF"]))
            except (ValueError, OverflowError):
                logtof = ""
            row = {"surface_name": surface, "TOF_C6H12O_s-1": r["TOF"], "log10_TOF": logtof,
                   "theta_OOH": r.get("coverages", {}).get("OOH_s", r.get("coverages", {}).get("OOH", "")),
                   "theta_Ha": r.get("coverages", {}).get("Ha_s", r.get("coverages", {}).get("Ha", "")),
                   "theta_vacancy": r.get("coverages", {}).get("*", ""),
                   "calculation_status": r.get("status", ""), "numerical_status": r.get("numerical_status", ""),
                   "physical_gate": r.get("physical_gate", ""), "source_hash": r.get("source_hash", ""),
                   "run_path": r.get("run_path", "")}
            for state in ("H2O2_s", "C6H12_s", "H2O_s", "C6H12O_s", "NEB1IS_s", "NEB1FS_s", "NEB2IS_s", "NEB2FS_s"):
                row[f"theta_{state}"] = r.get("coverages", {}).get(state, "")
            for i in range(1, 10):
                row[f"forward_flux_step_{i}"] = r.get("directional_fluxes", {}).get(f"forward_flux_step_{i}", "")
                row[f"reverse_flux_step_{i}"] = r.get("directional_fluxes", {}).get(f"reverse_flux_step_{i}", "")
                row[f"net_flux_step_{i}"] = r.get("fluxes", {}).get(f"net_flux_step_{i}", "")
            value_rows.append(row)
        en = r.get("energies", {})
        signed = r.get("NEB_signed", {})
        for neb, names in (("NEB1", ("NEB1IS", "NEB1TS", "NEB1FS")),
                           ("NEB2", ("NEB2IS", "NEB2TS", "NEB2FS"))):
            if all(name in en for name in names):
                ei, et, ef = (float(en[name]) for name in names)
                barrier_rows.append({"surface_name": surface, "path": neb, "G_IS_eV": ei, "G_TS_eV": et, "G_FS_eV": ef,
                                     "signed_TS_minus_IS_eV": signed.get(f"{neb}_TS_minus_IS_eV", et-ei),
                                     "signed_TS_minus_FS_eV": signed.get(f"{neb}_TS_minus_FS_eV", et-ef),
                                     "CatMAP_endpoint_scheme_forward_barrier_eV": max(ei, et, ef)-ei,
                                     "CatMAP_endpoint_scheme_reverse_barrier_eV": max(ei, et, ef)-ef,
                                     "barrier_status": "CATMAP_MAX_OF_IS_TS_FS; raw signed NEB differences retained",
                                     "physical_gate": r.get("physical_gate", "PROVISIONAL_UNVERIFIED_COMPOSITE_SITE_MAPPING")})
    write_csv(BATCH / "MULTISITE_CALCULATION_STATUS.csv", status_rows)
    write_csv(BATCH / "MULTISITE_TOF_COVERAGE.csv", value_rows)
    write_csv(BATCH / "MULTISITE_REACTION_BARRIERS.csv", barrier_rows)


def run_one(surface: str, precision: int, python: str, timeout_seconds: int) -> dict:
    if surface == "tico":
        return tico_reused()
    if surface == "tife":
        return tife_reused()
    run_root = BATCH / "runs"
    existing = load_batch_candidate(surface, run_root / f"{surface}_candidate")
    if existing is not None:
        return existing
    log_dir = BATCH / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    try:
        child_env = os.environ.copy()
        child_env["PYTHONIOENCODING"] = "utf-8"
        prep = subprocess.run([python, str(PREPARE), "--surface", surface], cwd=WORKFLOW,
                              env=child_env, text=True, capture_output=True, encoding="utf-8", errors="replace",
                              timeout=timeout_seconds)
        (log_dir / f"{surface}_prepare.log").write_text(prep.stdout + prep.stderr, encoding="utf-8")
        if prep.returncode:
            return {"surface_name": surface, "status": "ENERGY_MAPPING_FLAGGED",
                    "numerical_status": "NOT_RUN_PREPARATION_FAILED", "error": prep.stderr[-2000:],
                    "run_path": str(run_root / f"{surface}_candidate")}
        prepared_dir = WORKFLOW / "runs" / f"{surface}_candidate"
        destination_dir = run_root / f"{surface}_candidate"
        if destination_dir.exists():
            return {"surface_name": surface, "status": "ENERGY_MAPPING_FLAGGED",
                    "numerical_status": "ISOLATED_DESTINATION_ALREADY_EXISTS",
                    "error": f"refusing to overwrite {destination_dir}", "run_path": str(destination_dir)}
        destination_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(prepared_dir), str(destination_dir))
        model = destination_dir / f"candidate_{surface}.mkm"
        run_dir = model.parent
        run_precision = precision
        run = subprocess.run([python, str(SOLVER), "--model", str(model), "--precision", str(run_precision)],
                             cwd=run_dir, text=True, capture_output=True,
                             env=child_env, encoding="utf-8", errors="replace", timeout=timeout_seconds)
        log = log_dir / f"{surface}_p{run_precision}.log"
        log.write_text(run.stdout + run.stderr, encoding="utf-8")
        gate_path = model.parent / "numerical_gate.json"
        if run.returncode or not gate_path.is_file():
            # Preserve first-pass diagnostics before the higher-precision retry.
            attempt_dir = run_dir / "attempts" / f"precision_{run_precision}"
            attempt_dir.mkdir(parents=True, exist_ok=True)
            for path in run_dir.iterdir():
                if path.is_file() and path.name not in {model.name, f"candidate_{surface}_energies.txt", "run_manifest.json"}:
                    shutil.copy2(path, attempt_dir / path.name)
            retry_precision = max(220, 2 * run_precision)
            retry = subprocess.run([python, str(SOLVER), "--model", str(model), "--precision", str(retry_precision)],
                                   cwd=run_dir, text=True, capture_output=True,
                                   env=child_env, encoding="utf-8", errors="replace", timeout=timeout_seconds)
            run_precision = retry_precision
            run = retry
            log = log_dir / f"{surface}_p{run_precision}.log"
            log.write_text(run.stdout + run.stderr, encoding="utf-8")
        manifest = {"surface_name": surface, "python": python, "precision_digits": precision,
                    "actual_precision_digits": run_precision, "catmap_version_expected": "0.3.1", "model_sha256": sha(model),
                    "energies_sha256": sha(model.parent / f"candidate_{surface}_energies.txt"),
                    "source_manifest_sha256": sha(SNAP / "source_sha256.json"),
                    "run_exit_code": run.returncode, "log": str(log), "model": str(model)}
        (model.parent / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        if run.returncode or not gate_path.is_file():
            return {"surface_name": surface, "status": "NUMERICAL_FAILED",
                    "numerical_status": "SOLVER_FAILED", "error": run.stdout[-2000:] + run.stderr[-2000:],
                    "run_path": str(model.parent), "source_hash": manifest["energies_sha256"]}
        gate = json.loads(gate_path.read_text(encoding="utf-8"))
        coverages = {r["species"]: r["coverage"] for r in read_csv(model.parent / "candidate_steady_state_coverages.csv")}
        tof = next(r["turnover_frequency"] for r in read_csv(model.parent / "candidate_tof.csv")
                   if r.get("gas_species") == "C6H12O_g")
        flux_rows = read_csv(model.parent / "candidate_reaction_rates.csv")
        flux = {f"net_flux_step_{r['reaction_index']}": r["net_rate"] for r in flux_rows}
        energies = {r["species_name"]: r["formation_energy"]
                    for r in read_csv(model.parent / f"candidate_{surface}_energies.txt", "\t")}
        signed = {k: float(v) for k, v in json.loads((model.parent / "preparation_metadata.json").read_text(encoding="utf-8"))["NEB_signed_differences_eV"].items()}
        flags = []
        if signed["NEB1_TS_minus_IS_eV"] < 0:
            flags.append("NEB1_TS_BELOW_IS")
        if signed["NEB2_TS_minus_IS_eV"] < 0:
            flags.append("NEB2_TS_BELOW_IS")
        directional = {f"forward_flux_step_{r['reaction_index']}": r.get("forward_rate", "") for r in flux_rows}
        directional.update({f"reverse_flux_step_{r['reaction_index']}": r.get("reverse_rate", "") for r in flux_rows})
        return {"surface_name": surface, "status": "COMPUTED_PROVISIONAL",
                "numerical_status": gate["status"], "catmap_version": gate["catmap_version"],
                "direction_semantics": "REVERSIBLE_NATIVE_CATMAP", "TOF": tof,
                "coverages": coverages, "fluxes": flux, "directional_fluxes": directional, "energies": energies,
                "NEB_signed": signed, "negative_TS_flags": ";".join(flags),
                "max_relative_residual": gate.get("max_relative_surface_residual", ""),
                "site_balance_error": gate.get("site_balance_error", ""),
                "inventory_error": gate.get("inventory_error", ""),
                "source_hash": manifest["model_sha256"] + ":" + manifest["energies_sha256"] + ":" + manifest["source_manifest_sha256"],
                "physical_gate": "PROVISIONAL_UNVERIFIED_COMPOSITE_SITE_MAPPING" + (":" + ";".join(flags) if flags else ""),
                "run_path": str(model.parent), "reused": False}
    except subprocess.TimeoutExpired as exc:
        (BATCH / "logs").mkdir(parents=True, exist_ok=True)
        timeout_log = BATCH / "logs" / f"{surface}_timeout.log"
        timeout_log.write_text(f"TIMEOUT after {timeout_seconds}s\nstdout={exc.stdout!r}\nstderr={exc.stderr!r}\n", encoding="utf-8")
        return {"surface_name": surface, "status": "NUMERICAL_FAILED",
                "numerical_status": f"TIME_LIMIT_EXCEEDED_{timeout_seconds}s", "error": str(exc),
                "run_path": str(BATCH / "runs" / f"{surface}_candidate"), "log": str(timeout_log)}
    except Exception as exc:
        return {"surface_name": surface, "status": "NUMERICAL_FAILED",
                "numerical_status": "BATCH_EXCEPTION", "error": f"{type(exc).__name__}: {exc}",
                "run_path": str(BATCH / "runs" / f"{surface}_candidate")}


def main() -> None:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--surfaces", nargs="+", default=TARGETS)
    ap.add_argument("--precision", type=int, default=120)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--timeout-seconds", type=int, default=300)
    args = ap.parse_args()
    BATCH.mkdir(parents=True, exist_ok=True)
    results = []
    for surface in args.surfaces:
        print(f"START {surface}", flush=True)
        r = run_one(surface, args.precision, args.python, args.timeout_seconds)
        results.append(r)
        print(json.dumps({k: r.get(k) for k in ("surface_name", "status", "numerical_status", "TOF", "negative_TS_flags", "run_path", "error")}, ensure_ascii=False), flush=True)
        (BATCH / "batch_run_progress.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (BATCH / "batch_run_progress.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    assemble_outputs(results)
    readiness = WORKFLOW / "data" / "site_readiness.csv"
    if readiness.is_file():
        rows = read_csv(readiness)
        by_result = {r["surface_name"]: r for r in results}
        for row in rows:
            res = by_result.get(row["surface_name"])
            if res:
                row["calculation_status"] = res.get("status", "NOT_RUN")
                row["numerical_status"] = res.get("numerical_status", "")
                row["physical_gate"] = res.get("physical_gate", row.get("physical_gate", ""))
        write_csv(readiness, rows)
    print(f"complete={len(results)}; numerical_pass={sum(r.get('numerical_status') == 'NUMERICAL_PASS_PROVISIONAL_PHYSICS' for r in results)}", flush=True)


if __name__ == "__main__":
    main()
