#!/usr/bin/env python
"""Official CatMAP 0.3.1 runs for the four-parameter OOH*-extended BEP study.

This is orchestration/data validation only. All steady states are solved by
CatMAP ReactionModel + ThermodynamicScaler + SteadyStateSolver.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import gc
import hashlib
import json
import math
import os
import re
import shutil
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parents[1]
ROOT = Path(__file__).resolve().parents[4]
SNAP = HERE / "inputs" / "snapshots"
CASES = HERE / "inputs" / "cases"
RUNS = HERE / "runs"
DATA = HERE / "data"
FIGURES = HERE / "figures"
REPORTS = HERE / "reports"
SURFACES = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tico", "titi", "tiw", "tita", "ticr", "ti", "tini"]
LABELS = {
    "tife": "Ti–Fe", "timn": "Ti–Mn", "tihf": "Ti–Hf", "tire": "Ti–Re",
    "timo": "Ti–Mo", "tiv": "Ti–V", "tico": "Ti–Co", "titi": "Ti–Ti",
    "tiw": "Ti–W", "tita": "Ti–Ta", "ticr": "Ti–Cr", "ti": "Ti", "tini": "Ti–Ni",
}
TEMPERATURE = 333.15
PRESSURE_DESCRIPTOR = 1.0
TOLERANCE = "1e-70"
COVERAGE_SUM_TOL = "1e-30"
RUN_SCHEMA = "v2-eager-tableparser-input-bound-before-model-construction"
GAS_VALUES = {"H2O2_g": 0.185, "C6H12_g": 1.0, "H2O_g": 0.815, "C6H12O_g": 1e-20}
EXPECTED_RXNS = [
    "* + H2O2_g + C6H12_g <-> Rpair*",
    "Rpair* <-> TS1* -> EpoxideHHO*",
    "EpoxideHHO* <-> HHO* + C6H12O_g",
    "HHO* <-> TS2* -> OHHa*",
    "OHHa* <-> H2O*",
    "H2O* <-> H2O_g + *",
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def official_input_fingerprint(setup_path: Path, energy_path: Path, catmap_version: str) -> str:
    return sha256_bytes((sha256_file(setup_path) + sha256_file(energy_path) + catmap_version +
                         TEMPERATURE.__repr__() + TOLERANCE + RUN_SCHEMA).encode())


def setup_bound_to_energy(source_setup: Path, energy_path: Path, destination: Path) -> Path:
    """Create a case-local .mkm that points at its energy table before parsing.

    CatMAP 0.3.1 parses TableParser input in ReactionModel.__init__, so
    assigning model.input_file afterwards cannot change the already parsed
    energy definitions.
    """
    lines = source_setup.read_text(encoding="utf-8").splitlines()
    matches = [i for i, line in enumerate(lines) if re.match(r"^\s*input_file\s*=", line)]
    if len(matches) != 1:
        raise ValueError(f"Expected one input_file assignment in {source_setup}, found {len(matches)}")
    lines[matches[0]] = f"input_file = {str(energy_path.resolve())!r}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def decimal(value: Any) -> str:
    return str(value)


def prepare_snapshots() -> dict[str, Any]:
    """Copy immutable input snapshots once; refuse silent snapshot replacement."""
    sources: dict[str, Path] = {}
    strict = ROOT / "catmap_official_run" / "results" / "official_strict_tolerance"
    for surface in SURFACES:
        sources[f"{surface}/alkene_epoxidation.mkm"] = strict / surface / "alkene_epoxidation.mkm"
        sources[f"{surface}/energies.txt"] = strict / surface / "energies.txt"
    sources["source_carbon_O_H_TS1.csv"] = ROOT / "analyses" / "BEP_energy_correction_optimization" / "data" / "source_carbon_O_H_TS1.csv"
    sources["convergence_13sites.csv"] = ROOT / "analyses" / "microkinetic_mechanism_13sites" / "data" / "convergence_13sites.csv"
    sources["grid_BEP_OOH_dG2_TiCo_conditional.csv"] = ROOT / "analyses" / "BEP_OOH_activity_map" / "data" / "grid_BEP_OOH_dG2_TiCo_conditional.csv"
    metadata = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "sources": {}}
    for rel, src in sources.items():
        if not src.is_file():
            raise FileNotFoundError(f"Required source missing: {src}")
        dst = SNAP / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        source_hash = sha256_file(src)
        if dst.exists():
            if sha256_file(dst) != source_hash:
                raise RuntimeError(f"Existing immutable snapshot differs from current source; preserve both and investigate: {dst}")
        else:
            shutil.copy2(src, dst)
        metadata["sources"][rel] = {"source_path": str(src), "source_sha256": source_hash, "snapshot_path": str(dst), "snapshot_sha256": sha256_file(dst)}
    fit_rows = read_csv(SNAP / "source_carbon_O_H_TS1.csv")
    conv_rows = read_csv(SNAP / "convergence_13sites.csv")
    grid_rows = read_csv(SNAP / "grid_BEP_OOH_dG2_TiCo_conditional.csv")
    metadata["fit_input_rows"] = len(fit_rows)
    metadata["archive_rows"] = len(conv_rows)
    if not grid_rows:
        raise RuntimeError("Historical grid is empty; cannot recover approved bounds")
    xcol, ycol = "deltaG2_eV", "OOH_legacy_formation_eV"
    xb = [float(r[xcol]) for r in grid_rows]
    yb = [float(r[ycol]) for r in grid_rows]
    metadata["descriptor_range_source"] = {
        "source": str(sources["grid_BEP_OOH_dG2_TiCo_conditional.csv"]),
        "source_sha256": metadata["sources"]["grid_BEP_OOH_dG2_TiCo_conditional.csv"]["source_sha256"],
        "x_name": xcol, "x_min_eV": min(xb), "x_max_eV": max(xb),
        "y_name": ycol, "y_min_eV": min(yb), "y_max_eV": max(yb),
        "source_grid_rows": len(grid_rows),
    }
    write_json(SNAP / "input_provenance_and_ranges.json", metadata)
    return metadata


def source_tables() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    source = read_csv(SNAP / "source_carbon_O_H_TS1.csv")
    archived = read_csv(SNAP / "convergence_13sites.csv")
    fit_by_surface = {r["surface_name"].lower(): r for r in source}
    archive_by_surface = {r["surface"].lower(): r for r in archived}
    if len(fit_by_surface) != len(source) or len(archive_by_surface) != len(archived):
        raise ValueError("Duplicate surface_name in source or archived baseline tables")
    if set(fit_by_surface) != set(SURFACES) or set(archive_by_surface) != set(SURFACES):
        raise ValueError(f"13-site join mismatch. source={set(fit_by_surface)}, archive={set(archive_by_surface)}")
    for s in SURFACES:
        if fit_by_surface[s]["status"].upper() != "CONVERGED" or float(fit_by_surface[s]["temperature_K"]) != TEMPERATURE:
            raise ValueError(f"Fit source not a validated 333.15 K point: {s}")
        if archive_by_surface[s]["status"].upper() != "CONVERGED":
            raise ValueError(f"Archived reference not converged: {s}")
    return fit_by_surface, archive_by_surface


def read_energy_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No TableParser header: {path}")
        return list(reader.fieldnames), list(reader)


def energy_value(rows: list[dict[str, str]], surface: str, species: str) -> float:
    found = [r for r in rows if r["surface_name"].lower() == surface and r["species_name"] == species]
    if len(found) != 1:
        raise ValueError(f"Expected one {surface}/{species} formation energy, found {len(found)}")
    return float(found[0]["formation_energy"])


def write_energy_table(path: Path, headers: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def validate_snapshots() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, Any]]:
    from catmap import ReactionModel, __version__
    if __version__ != "0.3.1":
        raise RuntimeError(f"Expected CatMAP 0.3.1, found {__version__}")
    fit_by_surface, archive_by_surface = source_tables()
    source_hashes = []
    for surface in SURFACES:
        setup = SNAP / surface / "alkene_epoxidation.mkm"
        energy = SNAP / surface / "energies.txt"
        # CatMAP's TableParser opens the .mkm-relative input_file during
        # ReactionModel construction, so initialize from the immutable snapshot
        # directory rather than relying on the caller's working directory.
        with contextlib.chdir(setup.parent):
            model = ReactionModel(setup_file=str(setup))
        expected_rxns = [" ".join(x.split()) for x in EXPECTED_RXNS]
        got_rxns = [" ".join(x.split()) for x in model.rxn_expressions]
        if model.surface_names != [surface]:
            raise ValueError(f"surface_names mismatch in {setup}: {model.surface_names}")
        if got_rxns != expected_rxns:
            raise ValueError(f"Six-step network changed in {setup}")
        if float(model.temperature) != TEMPERATURE or model.pressure_mode != "concentration":
            raise ValueError(f"Temperature/pressure_mode mismatch in {setup}")
        if model.gas_thermo_mode != "frozen_gas" or model.adsorbate_thermo_mode != "frozen_adsorbate":
            raise ValueError(f"Frozen thermo modes mismatch in {setup}")
        if model.scaler.__class__.__name__ != "ThermodynamicScaler" or model.solver.__class__.__name__ != "SteadyStateSolver":
            raise ValueError(f"Official scaler/solver mismatch in {setup}")
        for gas, value in GAS_VALUES.items():
            if float(model.species_definitions[gas]["concentration"]) != value:
                raise ValueError(f"Reservoir mismatch {surface}/{gas}")
        headers, erows = read_energy_table(energy)
        energy_map = {(r["surface_name"].lower(), r["species_name"]): r for r in erows
                      if r["surface_name"].lower() == surface and r["surface_name"] not in ("", "None")}
        if len(energy_map) != 7:
            raise ValueError(f"Expected 7 {surface} energy entries in {energy}, got {len(energy_map)}")
        fit = fit_by_surface[surface]
        for species, field in (("Rpair", "Rpair_strict"), ("TS1", "TS1_strict"), ("EpoxideHHO", "EpoxideHHO_strict")):
            actual = energy_value(erows, surface, species)
            if not math.isclose(actual, float(fit[field]), rel_tol=0, abs_tol=1e-9):
                raise ValueError(f"Strict energy/source mismatch {surface}/{species}: {actual} vs {fit[field]}")
        source_hashes.append({"surface": surface, "setup_sha256": sha256_file(setup), "energies_sha256": sha256_file(energy)})
    hashes = {"catmap_version": __version__, "source_inputs": source_hashes,
              "bep_source_sha256": sha256_file(SNAP / "source_carbon_O_H_TS1.csv"),
              "baseline_reference_sha256": sha256_file(SNAP / "convergence_13sites.csv")}
    write_json(DATA / "input_validation.json", hashes)
    return fit_by_surface, archive_by_surface, hashes


def fit_bep(rows: list[dict[str, str]]) -> dict[str, Any]:
    x = np.asarray([float(r["NEB1_reaction_energy_EpoxideHHO_minus_Rpair"]) for r in rows], dtype=float)
    z = np.asarray([float(r["OOH_star_legacy"]) for r in rows], dtype=float)
    y = np.asarray([float(r["TS1_strict"]) - float(r["Rpair_strict"]) for r in rows], dtype=float)
    xc, zc = float(x.mean()), float(z.mean())
    u, v = x - xc, z - zc
    design = np.column_stack([np.ones(len(y)), u, v, v * v])
    coef, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    if rank != 4:
        raise ValueError(f"BEP design rank deficient: {rank}/4")
    pred = design @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    mae = float(np.mean(np.abs(y - pred)))
    loo = np.empty_like(y)
    loo_fits = []
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        xmean, zmean = float(x[keep].mean()), float(z[keep].mean())
        du, dz = x[keep] - xmean, z[keep] - zmean
        a = np.column_stack([np.ones(keep.sum()), du, dz, dz * dz])
        c, _, r, _ = np.linalg.lstsq(a, y[keep], rcond=None)
        if r != 4:
            raise ValueError(f"LOOCV fold {i} rank deficient: {r}/4")
        dxi, dzi = x[i] - xmean, z[i] - zmean
        loo[i] = float(np.array([1.0, dxi, dzi, dzi * dzi]) @ c)
        loo_fits.append({"surface": rows[i]["surface_name"], "center_dG2_eV": xmean, "center_OOH_eV": zmean, "coefficients": [float(q) for q in c]})
    loo_rmse = float(np.sqrt(np.mean((y - loo) ** 2)))
    return {"center_dG2_eV": xc, "center_OOH_eV": zc, "coefficients_intercept_dG2_OOH_OOH2": [float(q) for q in coef],
            "n": len(y), "rank": int(rank), "R2": r2, "RMSE_eV": rmse, "MAE_eV": mae,
            "LOOCV_RMSE_eV": loo_rmse, "actual_y": y.tolist(), "fitted_y": pred.tolist(),
            "LOOCV_y": loo.tolist(), "folds": loo_fits}


def make_case_energy(surface: str, variant: str, replacement: dict[str, float] | None = None) -> tuple[Path, str]:
    src = SNAP / surface / "energies.txt"
    raw = src.read_bytes()
    headers, rows = read_energy_table(src)
    if replacement:
        for species, value in replacement.items():
            match = [r for r in rows if r["surface_name"].lower() == surface and r["site_name"] == "active" and r["species_name"] == species]
            if len(match) != 1:
                raise ValueError(f"Expected one replacement row {surface}/{species}; found {len(match)}")
            match[0]["formation_energy"] = format(float(value), ".17g")
    dest = CASES / variant / surface / "energies.txt"
    if replacement:
        write_energy_table(dest, headers, rows)
    else:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            if sha256_file(dest) != sha256_bytes(raw):
                raise RuntimeError(f"Case input exists with unexpected bytes: {dest}")
        else:
            dest.write_bytes(raw)
    return dest, sha256_file(dest)


def _as_strings(values: Any) -> list[str]:
    return [str(v) for v in values]


def run_official_case(surface: str, variant: str, energy_path: Path, seed: dict[str, Any] | None = None,
                      coordinate: list[float] | None = None, physical: dict[str, Any] | None = None) -> dict[str, Any]:
    """One fresh official CatMAP ReactionModel; outputs exact solver state as JSON."""
    from catmap import ReactionModel, __version__
    if __version__ != "0.3.1":
        raise RuntimeError(f"CatMAP version changed during run: {__version__}")
    setup_path = SNAP / surface / "alkene_epoxidation.mkm"
    run_dir = RUNS / variant / surface
    run_dir.mkdir(parents=True, exist_ok=True)
    result_path = run_dir / "official_result.json"
    input_fingerprint = official_input_fingerprint(setup_path, energy_path, __version__)
    if result_path.exists():
        old = json.loads(result_path.read_text(encoding="utf-8"))
        if old.get("input_sha256") == input_fingerprint and old.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
            return old
        # Preserve failed/stale attempts; retries receive fresh directories.
        attempts = run_dir / "attempts"
        attempt = 2
        while True:
            candidate = attempts / f"attempt_{attempt:02d}"
            candidate_result = candidate / "official_result.json"
            if not candidate_result.exists():
                run_dir, result_path = candidate, candidate_result
                run_dir.mkdir(parents=True, exist_ok=True)
                break
            prior = json.loads(candidate_result.read_text(encoding="utf-8"))
            if prior.get("input_sha256") == input_fingerprint and prior.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
                return prior
            attempt += 1
    import io
    import mpmath as mp
    captured_out, captured_err = io.StringIO(), io.StringIO()
    status = "FAILED"
    result: dict[str, Any] = {
        "surface": surface, "catalyst": LABELS[surface], "variant": variant,
        "coordinate": coordinate, "input_sha256": input_fingerprint,
        "energy_input_sha256": sha256_file(energy_path), "setup_sha256": sha256_file(setup_path),
        "physical": physical or {}, "catmap_version": __version__,
    }
    try:
        case_setup = setup_bound_to_energy(setup_path, energy_path, run_dir / "catmap_case.mkm")
        result["case_setup_sha256"] = sha256_file(case_setup)
        with contextlib.redirect_stdout(captured_out), contextlib.redirect_stderr(captured_err):
            # The parser reads the setup's relative energies.txt immediately.
            with contextlib.chdir(setup_path.parent):
                model = ReactionModel(setup_file=str(case_setup))
            if Path(model.input_file).resolve() != energy_path.resolve():
                raise RuntimeError(f"CatMAP parsed unexpected input_file: {model.input_file}")
            parsed_energy_checks = {}
            source_headers, source_energy_rows = read_energy_table(energy_path)
            del source_headers
            for species in ("Rpair", "TS1", "EpoxideHHO", "HHO", "TS2", "OHHa", "H2O"):
                expected = energy_value(source_energy_rows, surface, species)
                parsed = model.species_definitions[f"{species}_s"]["formation_energy"]
                if not isinstance(parsed, (list, tuple)) or len(parsed) != 1:
                    raise RuntimeError(f"Unexpected parsed formation_energy shape for {surface}/{species}: {parsed!r}")
                parsed_value = float(parsed[0])
                if not math.isclose(parsed_value, expected, rel_tol=0, abs_tol=1e-12):
                    raise RuntimeError(f"CatMAP parsed wrong energy for {surface}/{species}: {parsed_value} != {expected}")
                parsed_energy_checks[species] = {"input_eV": expected, "parsed_eV": parsed_value, "match": True}
            result["parsed_energy_checks"] = parsed_energy_checks
            model.data_file = str((run_dir / "catmap_data.pkl").resolve())
            model.log_file = str((run_dir / "catmap.log").resolve())
            if seed is not None:
                expected_seed_labels = list(model.adsorbate_names) + [s for s in model.site_names if s != "g"]
                if seed.get("coverage_labels") != expected_seed_labels:
                    raise ValueError("Continuation seed coverage labels do not match CatMAP adsorbate order")
                point = [TEMPERATURE, PRESSURE_DESCRIPTOR]
                model.coverage_map = [[point, [mp.mpf(x) for x in seed["coverages"]]]]
                model.numbers_map = [[point, [mp.mpf(x) for x in seed["numbers"]]]]
            model.run(recalculate=True)
            if list(model.surface_names) != [surface] or float(model.temperature) != TEMPERATURE:
                raise RuntimeError("CatMAP model identity or temperature changed")
            if model.scaler.__class__.__name__ != "ThermodynamicScaler" or model.solver.__class__.__name__ != "SteadyStateSolver":
                raise RuntimeError("Unexpected CatMAP scaler or solver")
            if str(model.solver.tolerance) != TOLERANCE:
                if mp.mpf(str(model.solver.tolerance)) != mp.mpf(TOLERANCE):
                    raise RuntimeError(f"CatMAP tolerance changed: {model.solver.tolerance}")
            map_lengths = {name: len(getattr(model, name)) for name in
                           ("coverage_map", "rate_map", "turnover_frequency_map", "directional_rates_map")}
            if any(value != 1 for value in map_lengths.values()):
                raise RuntimeError(f"CatMAP single-point map outputs have unexpected lengths: {map_lengths}")
            coverages = model.coverage_map[0][1]
            cvg_labels = list(model.output_labels["coverage"])
            sites = [s for s in model.site_names if s != "g"]
            if len(coverages) != len(cvg_labels) + len(sites):
                raise RuntimeError(f"Coverage values/labels length mismatch: {len(coverages)} vs {len(cvg_labels)}+{len(sites)}")
            coverage_labels = cvg_labels + sites
            residual = model.solver.get_residual(coverages)
            coverage_sum = model._math.fsum(coverages)
            tof_labels = list(model.output_labels["turnover_frequency"])
            product_idx = tof_labels.index("C6H12O_g")
            tof = model.turnover_frequency_map[0][1][product_idx]
            net_rates = model.rate_map[0][1]
            dir_labels = list(model.output_labels["directional_rates"])
            directional = model.directional_rates_map[0][1]
            if len(net_rates) != 6 or len(directional) != 12:
                raise RuntimeError(f"Unexpected six-step outputs: net={len(net_rates)}, directional={len(directional)}")
            if not all(mp.isfinite(mp.mpf(v)) for v in coverages) or not mp.isfinite(mp.mpf(tof)) or tof <= 0:
                raise ValueError("Non-finite/non-positive CatMAP coverage or product TOF")
            if any(mp.mpf(v) < 0 for v in coverages):
                raise ValueError("CatMAP returned negative coverage")
            flux_max = max(abs(mp.mpf(v)) for v in net_rates)
            flux_spread = (max(mp.mpf(v) for v in net_rates) - min(mp.mpf(v) for v in net_rates)) / max(flux_max, mp.mpf("1e-10000"))
            checks = {
                "residual_pass": mp.mpf(residual) <= mp.mpf(TOLERANCE),
                "coverage_sum_pass": abs(mp.mpf(coverage_sum) - 1) <= mp.mpf(COVERAGE_SUM_TOL),
                "nonnegative_coverages_pass": all(mp.mpf(v) >= 0 for v in coverages),
            }
            status = "CONVERGED" if all(checks.values()) else "FAILED_VALIDATION"
            result.update({
                "status": status, "tof_C6H12O": str(tof), "log10_tof_C6H12O": str(mp.log10(tof)),
                "coverages": _as_strings(coverages), "coverage_labels": coverage_labels,
                "coverage_sum": str(coverage_sum), "residual": str(residual), "tolerance": str(model.solver.tolerance),
                "net_rates": _as_strings(net_rates), "directional_rate_labels": dir_labels,
                "directional_rates": _as_strings(directional), "flux_spread_relative": str(flux_spread),
                "checks": checks, "numbers": _as_strings(model.numbers_map[0][1]),
                "numbers_labels": list(model.output_labels.get("coverage", ())) + sites,
                "surface_names": list(model.surface_names), "temperature_K": float(model.temperature),
                "pressure_mode": model.pressure_mode, "scaler": model.scaler.__class__.__name__,
                "solver": model.solver.__class__.__name__, "decimal_precision": model.decimal_precision,
                "reservoirs": GAS_VALUES,
            })
            if status != "CONVERGED":
                raise ValueError(f"Official CatMAP returned values failing acceptance checks: {checks}")
    except Exception as exc:
        result["status"] = "FAILED" if status != "FAILED_VALIDATION" else status
        result["exception"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()
    result["console_stdout"] = captured_out.getvalue()
    result["console_stderr"] = captured_err.getvalue()
    write_json(result_path, result)
    log_path = run_dir / "console.log"
    log_path.write_text(result.get("console_stdout", "") + "\n--- STDERR ---\n" + result.get("console_stderr", "") +
                        ("\n--- TRACEBACK ---\n" + result.get("traceback", "") if result.get("traceback") else ""), encoding="utf-8", errors="replace")
    # The compact JSON is the checkpoint; CatMAP's optional pickle is not required for a restart.
    for p in (run_dir / "catmap_data.pkl", run_dir / "catmap.log"):
        if p.exists():
            p.unlink()
    return result


def fit_rows_for_surface(rows: list[dict[str, str]]) -> tuple[dict[str, Any], dict[str, Any]]:
    model = fit_bep(rows)
    full_preds, loo_preds = model["fitted_y"], model["LOOCV_y"]
    return model, {r["surface_name"].lower(): (float(full_preds[i]), float(loo_preds[i])) for i, r in enumerate(rows)}


def compute_physical_flags(y_pred: float, dg2: float) -> dict[str, Any]:
    endpoint_max = max(0.0, dg2)
    flags = []
    if y_pred < 0:
        flags.append("TS1_below_initial_state")
    if y_pred < dg2:
        flags.append("TS1_below_final_state")
    catmap_top = max(0.0, dg2, y_pred)
    return {
        "signed_TS1_minus_Rpair_eV": y_pred,
        "deltaG2_eV": dg2,
        "TS1_below_initial": y_pred < 0,
        "TS1_below_final": y_pred < dg2,
        "independent_TS_above_both_endpoints": y_pred >= endpoint_max,
        "physical_flag": ";".join(flags) if flags else "TS1_above_both_endpoints",
        "catmap_forward_effective_barrier_eV": catmap_top,
        "catmap_reverse_effective_barrier_eV": catmap_top - dg2,
        "catmap_effective_TS_control": "TS1" if y_pred >= endpoint_max else ("final_state" if dg2 >= 0 else "initial_state"),
        "catmap_endpoint_clamp_applied": y_pred < endpoint_max,
    }


def base_energy_values(surface: str) -> tuple[list[str], list[dict[str, str]], dict[str, float]]:
    p = SNAP / surface / "energies.txt"
    headers, rows = read_energy_table(p)
    vals = {sp: energy_value(rows, surface, sp) for sp in ("Rpair", "TS1", "EpoxideHHO", "HHO", "TS2", "OHHa", "H2O")}
    return headers, rows, vals


def update_energy_rows(surface: str, rows: list[dict[str, str]], replacement: dict[str, float]) -> list[dict[str, str]]:
    out = [dict(r) for r in rows]
    for sp, value in replacement.items():
        found = [r for r in out if r["surface_name"].lower() == surface and r["site_name"] == "active" and r["species_name"] == sp]
        if len(found) != 1:
            raise ValueError(f"Expected one {surface}/{sp} source row, found {len(found)}")
        found[0]["formation_energy"] = format(float(value), ".17g")
    return out


def result_row(result: dict[str, Any], extra: dict[str, Any] | None = None) -> dict[str, Any]:
    base = {"surface_name": result.get("surface"), "catalyst": result.get("catalyst"), "status": result.get("status"),
            "TOF_C6H12O": result.get("tof_C6H12O"), "log10_TOF_C6H12O": result.get("log10_tof_C6H12O"),
            "steady_state_residual": result.get("residual"), "solver_tolerance": result.get("tolerance"),
            "coverage_sum": result.get("coverage_sum"), "flux_spread_relative": result.get("flux_spread_relative"),
            "input_sha256": result.get("input_sha256"), "failure": result.get("exception", "")}
    if extra:
        base.update(extra)
    return base


def flatten_coverages(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for r in results:
        if r.get("status") not in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
            continue
        row = {"surface_name": r["surface"], "catalyst": r["catalyst"], "variant": r["variant"], "status": r["status"]}
        row.update(dict(zip(r["coverage_labels"], r["coverages"])))
        rows.append(row)
    return rows


def flatten_rates(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for r in results:
        if r.get("status") not in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
            continue
        row = {"surface_name": r["surface"], "catalyst": r["catalyst"], "variant": r["variant"], "status": r["status"]}
        row.update({f"step_{i+1}_net_rate": v for i, v in enumerate(r["net_rates"])})
        row.update({f"direction_{i+1}_{label.rsplit(' ', 1)[-1]}": v for i, (label, v) in enumerate(zip(r["directional_rate_labels"], r["directional_rates"]))})
        rows.append(row)
    return rows


def latest_converged_seed(variant: str, surface: str) -> dict[str, Any]:
    base = RUNS / variant / surface
    candidates = []
    if base.exists():
        for path in base.rglob("official_result.json"):
            try:
                result = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if result.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"} and result.get("coverages") and result.get("numbers"):
                candidates.append((path.stat().st_mtime_ns, result))
    if not candidates:
        raise RuntimeError(f"No converged CatMAP seed available for {variant}/{surface}")
    return max(candidates, key=lambda item: item[0])[1]


def run_baselines(fit_by: dict[str, dict[str, str]], archived: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    results = []
    rows = []
    max_diff = 0.0
    for surface in SURFACES:
        energy = SNAP / surface / "energies.txt"
        energy_hash = sha256_file(energy)
        variant = "baseline_original"
        case_energy, _ = make_case_energy(surface, variant, None)
        run = run_official_case(surface, variant, case_energy)
        fit = fit_by[surface]
        ref = archived[surface]
        row = result_row(run, {"archived_TOF_C6H12O": ref["tof_C6H12O"], "archived_log10_TOF": ref["log10_tof"],
                               "absolute_log10_difference": "", "baseline_energy_sha256": energy_hash,
                               "original_signed_TS1_minus_Rpair_eV": float(fit["TS1_strict"]) - float(fit["Rpair_strict"])})
        if run.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
            diff = abs(float(run["log10_tof_C6H12O"]) - float(ref["log10_tof"]))
            row["absolute_log10_difference"] = diff
            max_diff = max(max_diff, diff)
            if diff > 1e-9:
                row["status"] = "BASELINE_MISMATCH"
                run["status"] = "BASELINE_MISMATCH"
        else:
            row["status"] = "BASELINE_FAILED"
        results.append(run)
        rows.append(row)
    write_csv(DATA / "baseline_validation.csv", rows)
    write_csv(DATA / "13site_original_TOF.csv", rows)
    write_csv(DATA / "coverage_original.csv", flatten_coverages(results))
    write_csv(DATA / "net_rates_original.csv", flatten_rates(results))
    if any(r["status"] != "CONVERGED" for r in rows) or max_diff > 1e-9:
        raise RuntimeError(f"Original CatMAP baseline validation failed; max |Δlog10TOF|={max_diff:.6g}. Stop BEP and grid stages.")
    write_json(DATA / "baseline_validation_summary.json", {"sites": len(rows), "converged": sum(r["status"] == "CONVERGED" for r in rows), "max_abs_log10_TOF_difference": max_diff, "gate": "PASS"})
    return results


def run_fit_scenarios(fit_by: dict[str, dict[str, str]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    rows = [fit_by[s] for s in SURFACES]
    fit = fit_bep(rows)
    if abs(fit["R2"] - 0.8866652) > 5e-5:
        raise RuntimeError(f"Re-fitted OOH BEP R2 differs from expected audit value: {fit['R2']}")
    write_json(DATA / "exact_BEP_fit.json", fit)
    write_csv(DATA / "exact_BEP_coefficients.csv", [
        {"term": name, "coefficient": value, "center_dG2_eV": fit["center_dG2_eV"], "center_OOH_eV": fit["center_OOH_eV"], "R2": fit["R2"], "RMSE_eV": fit["RMSE_eV"], "LOOCV_RMSE_eV": fit["LOOCV_RMSE_eV"]}
        for name, value in zip(["intercept_eV", "u_dG2", "z_OOH", "z_OOH_squared"], fit["coefficients_intercept_dG2_OOH_OOH2"])
    ])
    full_results, loo_results, prediction_rows = [], [], []
    baseline_seeds = {surface: latest_converged_seed("baseline_original", surface) for surface in SURFACES}
    for i, surface in enumerate(SURFACES):
        source = fit_by[surface]
        x = float(source["NEB1_reaction_energy_EpoxideHHO_minus_Rpair"])
        z = float(source["OOH_star_legacy"])
        actual_rpair = float(source["Rpair_strict"])
        actual_ts1 = float(source["TS1_strict"])
        actual_y = actual_ts1 - actual_rpair
        full_y = fit["fitted_y"][i]
        loo_y = fit["LOOCV_y"][i]
        loo_spec = fit["folds"][i]
        phys_full = compute_physical_flags(full_y, x)
        phys_loo = compute_physical_flags(loo_y, x)
        prediction_rows.append({"surface_name": surface, "catalyst": LABELS[surface], "deltaG2_eV": x, "OOH_legacy_single_adsorbate_eV": z,
                                "actual_Rpair_eV": actual_rpair, "actual_TS1_eV": actual_ts1, "actual_signed_TS1_minus_Rpair_eV": actual_y,
                                "full_fit_predicted_signed_Y_eV": full_y, "full_fit_predicted_TS1_eV": actual_rpair + full_y,
                                "full_fit_TS1_error_eV": actual_rpair + full_y - actual_ts1,
                                "LOOCV_predicted_signed_Y_eV": loo_y, "LOOCV_predicted_TS1_eV": actual_rpair + loo_y,
                                "LOOCV_TS1_error_eV": actual_rpair + loo_y - actual_ts1,
                                "LOOCV_center_dG2_eV": loo_spec["center_dG2_eV"], "LOOCV_center_OOH_eV": loo_spec["center_OOH_eV"],
                                "LOOCV_coefficients": json.dumps(loo_spec["coefficients"]),
                                "OOH_definition": "historical single-adsorbate formation energy; not bridged coadsorbed OOH"})
        for variant, y_pred, physical in (("BEP_fitted", full_y, phys_full), ("BEP_LOOCV", loo_y, phys_loo)):
            header, energy_rows, vals = base_energy_values(surface)
            new_rows = update_energy_rows(surface, energy_rows, {"TS1": actual_rpair + y_pred})
            energy_path = CASES / variant / surface / "energies.txt"
            write_energy_table(energy_path, header, new_rows)
            run = run_official_case(surface, variant, energy_path, seed=baseline_seeds[surface], physical=physical)
            run.update({"actual_Rpair_eV": actual_rpair, "actual_TS1_eV": actual_ts1, "predicted_TS1_eV": actual_rpair + y_pred,
                        "predicted_signed_Y_eV": y_pred, "deltaG2_eV": x, "OOH_legacy_eV": z,
                        "actual_minus_predicted_TS1_eV": actual_ts1 - (actual_rpair + y_pred)})
            (full_results if variant == "BEP_fitted" else loo_results).append(run)
    if len(full_results) != 13 or len(loo_results) != 13:
        raise RuntimeError("BEP CatMAP cases incomplete")
    write_csv(DATA / "BEP_TS1_predictions_and_LOOCV.csv", prediction_rows)
    for results, variant, path in ((full_results, "BEP_fitted", "13site_BEP_fitted_TOF.csv"), (loo_results, "BEP_LOOCV", "13site_BEP_LOOCV_TOF.csv")):
        rowsout = []
        for run in results:
            rowsout.append(result_row(run, {"actual_TS1_eV": run.get("actual_TS1_eV"), "predicted_TS1_eV": run.get("predicted_TS1_eV"),
                                            "actual_minus_predicted_TS1_eV": run.get("actual_minus_predicted_TS1_eV"),
                                            "deltaG2_eV": run.get("deltaG2_eV"), "OOH_legacy_single_adsorbate_eV": run.get("OOH_legacy_eV"),
                                            **run.get("physical", {})}))
        write_csv(DATA / path, rowsout)
        write_csv(DATA / f"coverage_{variant}.csv", flatten_coverages(results))
        write_csv(DATA / f"net_rates_{variant}.csv", flatten_rates(results))
    if any(r.get("status") not in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"} for r in full_results + loo_results):
        raise RuntimeError("At least one BEP/LOOCV CatMAP case did not pass; grid may proceed only after inspecting failure; automatic stop applied.")
    write_json(DATA / "BEP_fit_summary.json", {k: fit[k] for k in ("center_dG2_eV", "center_OOH_eV", "coefficients_intercept_dG2_OOH_OOH2", "R2", "RMSE_eV", "MAE_eV", "LOOCV_RMSE_eV")})
    return full_results, loo_results, fit


def grid_axes(resolution: int) -> tuple[np.ndarray, np.ndarray]:
    ranges = json.loads((SNAP / "input_provenance_and_ranges.json").read_text(encoding="utf-8"))["descriptor_range_source"]
    return (np.linspace(ranges["x_min_eV"], ranges["x_max_eV"], resolution),
            np.linspace(ranges["y_min_eV"], ranges["y_max_eV"], resolution))


def make_grid_energy(x: float, z: float, fit: dict[str, Any], node_key: str) -> tuple[Path, dict[str, Any], str]:
    surface = "tico"
    headers, rows, vals = base_energy_values(surface)
    dg2 = x
    z0 = fit["center_OOH_eV"]
    u = dg2 - fit["center_dG2_eV"]
    zc = z - z0
    c0, cu, cz, cz2 = fit["coefficients_intercept_dG2_OOH_OOH2"]
    y = c0 + cu * u + cz * zc + cz2 * zc * zc
    ts1 = vals["Rpair"] + y
    epoxide = vals["Rpair"] + dg2
    flags = compute_physical_flags(y, dg2)
    flags["OOH_definition"] = "historical single-adsorbate formation energy used only as external BEP descriptor"
    flags["TS1_predicted_eV"] = ts1
    flags["Rpair_fixed_eV"] = vals["Rpair"]
    flags["EpoxideHHO_predicted_eV"] = epoxide
    flags["TS2_fixed_eV"] = vals["TS2"]
    flags["HHO_fixed_eV"] = vals["HHO"]
    flags["OHHa_fixed_eV"] = vals["OHHa"]
    flags["H2O_fixed_eV"] = vals["H2O"]
    rows = update_energy_rows(surface, rows, {"EpoxideHHO": epoxide, "TS1": ts1})
    input_path = CASES / "grid" / f"{node_key}" / "energies.txt"
    write_energy_table(input_path, headers, rows)
    setup_hash = sha256_file(SNAP / surface / "alkene_epoxidation.mkm")
    from catmap import __version__
    fingerprint = official_input_fingerprint(SNAP / surface / "alkene_epoxidation.mkm", input_path, __version__)
    return input_path, flags, fingerprint


def cached_result(input_hash: str) -> dict[str, Any] | None:
    p = RUNS / "cache" / f"{input_hash}.json"
    if not p.is_file():
        return None
    r = json.loads(p.read_text(encoding="utf-8"))
    if r.get("input_sha256") != input_hash or r.get("status") not in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
        return None
    r["cache_solver_status"] = r.get("solver_status", r.get("status"))
    r["retrieval_status"] = "exact_input_hash_cache_hit"
    r["status"] = "REUSED_CACHE"
    return r


def persist_cache(input_hash: str, result: dict[str, Any]) -> None:
    if result.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE"}:
        write_json(RUNS / "cache" / f"{input_hash}.json", result)


def adaptive_solve(surface: str, variant: str, input_path: Path, coord: list[float], phys: dict[str, Any],
                   seed: dict[str, Any] | None, depth: int = 3, adaptive_log: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    direct = run_official_case(surface, variant, input_path, seed=seed, coordinate=coord, physical=phys)
    if direct.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
        return direct
    if seed is None or depth <= 0:
        return direct
    prior_coord = seed.get("grid_coordinate")
    if not prior_coord:
        return direct
    midpoint = [(float(prior_coord[j]) + float(coord[j])) / 2.0 for j in range(2)]
    midkey = f"adaptive_d{depth}_{coord[0]:+.12f}_{coord[1]:+.12f}"
    midpath, midphys, midhash = make_grid_energy(midpoint[0], midpoint[1], FIT_CURRENT, midkey)
    cached = cached_result(midhash)
    if cached:
        mid = cached
    else:
        mid = run_official_case("tico", f"adaptive_{midkey}", midpath, seed=seed, coordinate=midpoint, physical=midphys)
        mid["grid_coordinate"] = midpoint
        mid["input_sha256"] = midhash
        persist_cache(midhash, mid)
    if adaptive_log is not None:
        adaptive_log.append({"coordinate": midpoint, "depth": depth, "status": mid.get("status"), "input_sha256": midhash})
    if mid.get("status") not in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
        return adaptive_solve(surface, variant, input_path, coord, phys, seed, depth - 1, adaptive_log)
    mid["grid_coordinate"] = midpoint
    retried = run_official_case(surface, variant + "_after_bisection", input_path, seed=mid, coordinate=coord, physical=phys)
    if retried.get("status") in {"CONVERGED", "REUSED_CACHE"}:
        retried["status"] = "CONVERGED_ADAPTIVE"
        retried["adaptive_depth"] = 4 - depth
        retried["grid_coordinate"] = coord
    return retried


FIT_CURRENT: dict[str, Any] = {}


def run_grid(resolution: int) -> list[dict[str, Any]]:
    global FIT_CURRENT
    fit_by, _ = source_tables()
    fit_rows = [fit_by[s] for s in SURFACES]
    FIT_CURRENT = fit_bep(fit_rows)
    xvals, yvals = grid_axes(resolution)
    variant = f"grid_{resolution}x{resolution}"
    node_rows: list[dict[str, Any]] = []
    last_seed: dict[str, Any] | None = None
    previous_x_index: int | None = None
    for j, y in enumerate(yvals):
        x_order = range(resolution) if j % 2 == 0 else range(resolution - 1, -1, -1)
        for i in x_order:
            x = float(xvals[i]); z = float(y)
            node_key = f"{variant}_i{i:03d}_j{j:03d}"
            input_path, physical, fingerprint = make_grid_energy(x, z, FIT_CURRENT, node_key)
            cached = cached_result(fingerprint)
            if cached:
                result = cached
                # Keep the original CatMAP outcome separate from the fact that
                # this export invocation retrieved an exact-hash cached result.
                direct_node_result = RUNS / node_key / "tico" / "official_result.json"
                if direct_node_result.is_file():
                    direct = json.loads(direct_node_result.read_text(encoding="utf-8"))
                else:
                    direct = None
                if direct and direct.get("input_sha256") == fingerprint and direct.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE"}:
                    result["solver_status"] = direct["status"]
                    result["execution_source"] = "direct_CatMAP_solve_in_this_grid_resolution"
                    result["retrieval_status"] = "exact_hash_result_reloaded_for_export"
                else:
                    result["solver_status"] = result.get("cache_solver_status", result.get("status"))
                    source_variant = result.get("variant", "prior_valid_run")
                    if source_variant == "grid_61x61":
                        result["execution_source"] = "validated_CatMAP_result_from_original_61x61_stage"
                    else:
                        result["execution_source"] = f"exact_hash_cache_from_{source_variant}"
                    result["retrieval_status"] = "exact_input_hash_cache_hit"
                result["status"] = result["solver_status"]
                result["variant"] = variant
                result["surface"] = "tico"
                result["catalyst"] = "Ti–Co fixed background"
                result["coordinate"] = [x, z]
            else:
                seed = last_seed
                # Reuse an exact valid 5x5/15x15 node as a seed at a coincident coordinate only if it belongs to the same fixed CatMAP energy input.
                result = run_official_case("tico", node_key, input_path, seed=seed, coordinate=[x, z], physical=physical)
                if result.get("status") in {"FAILED", "FAILED_VALIDATION"} and seed is not None:
                    adaptive_log: list[dict[str, Any]] = []
                    result = adaptive_solve("tico", node_key + "_adaptive", input_path, [x, z], physical, seed, 3, adaptive_log)
                    result["adaptive_attempts"] = adaptive_log
                    if result.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
                        result["status"] = "CONVERGED_ADAPTIVE"
                result["grid_coordinate"] = [x, z]
                result["variant"] = variant
                result["surface"] = "tico"
                result["catalyst"] = "Ti–Co fixed background"
                result["input_sha256"] = fingerprint
                result["solver_status"] = result.get("status")
                result["execution_source"] = "direct_CatMAP_solve_in_this_grid_resolution"
                result["retrieval_status"] = "fresh_solve"
                persist_cache(fingerprint, result)
            result["grid_i"] = i
            result["grid_j"] = j
            result["deltaG2_eV"] = x
            result["OOH_legacy_formation_eV"] = z
            result["physical"] = physical
            result_path = RUNS / variant / "nodes" / f"{node_key}.json"
            write_json(result_path, result)
            if result.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}:
                last_seed = result
                previous_x_index = i
            node_rows.append(result)
        # Durable row checkpoint: every cell has a JSON result before proceeding.
        rows_out = [grid_csv_row(r) for r in node_rows]
        write_csv(DATA / f"grid_{resolution}x{resolution}_checkpoint.csv", rows_out)
        print(f"GRID_PROGRESS resolution={resolution} row={j+1}/{resolution} completed={len(node_rows)}/{resolution*resolution} failures={sum(r.get('status') in {'FAILED','FAILED_VALIDATION'} for r in node_rows)}", flush=True)
    if len(node_rows) != resolution * resolution:
        raise RuntimeError(f"Grid cardinality mismatch {len(node_rows)} != {resolution*resolution}")
    coords = {(r["grid_i"], r["grid_j"]) for r in node_rows}
    if len(coords) != resolution * resolution:
        raise RuntimeError("Grid contains duplicate/missing index pairs")
    write_csv(DATA / f"grid_{resolution}x{resolution}.csv", [grid_csv_row(r) for r in sorted(node_rows, key=lambda q: (q["grid_j"], q["grid_i"]))])
    write_csv(DATA / f"grid_{resolution}x{resolution}_coverages.csv", flatten_coverages(node_rows))
    write_csv(DATA / f"grid_{resolution}x{resolution}_net_rates.csv", flatten_rates(node_rows))
    return node_rows


def grid_csv_row(r: dict[str, Any]) -> dict[str, Any]:
    p = r.get("physical", {})
    return {"grid_i": r.get("grid_i"), "grid_j": r.get("grid_j"), "deltaG2_eV": r.get("deltaG2_eV"),
            "OOH_legacy_formation_eV": r.get("OOH_legacy_formation_eV"), "status": r.get("solver_status", r.get("status")),
            "retrieval_status": r.get("retrieval_status", "unknown"), "execution_source": r.get("execution_source", "unknown"),
            "TOF_C6H12O": r.get("tof_C6H12O"), "log10_TOF_C6H12O": r.get("log10_tof_C6H12O"),
            "residual": r.get("residual"), "tolerance": r.get("tolerance"), "coverage_sum": r.get("coverage_sum"),
            "flux_spread_relative": r.get("flux_spread_relative"), "TS1_predicted_eV": p.get("TS1_predicted_eV"),
            "signed_TS1_minus_Rpair_eV": p.get("signed_TS1_minus_Rpair_eV"), "deltaG2_for_TS1_endpoint_check_eV": p.get("deltaG2_eV"),
            "TS1_below_initial": p.get("TS1_below_initial"), "TS1_below_final": p.get("TS1_below_final"),
            "independent_TS_above_both_endpoints": p.get("independent_TS_above_both_endpoints"),
            "CatMAP_endpoint_clamp_applied": p.get("catmap_endpoint_clamp_applied"),
            "CatMAP_effective_TS_control": p.get("catmap_effective_TS_control"),
            "input_sha256": r.get("input_sha256"), "failure": r.get("exception", "")}


def catalyst_points() -> list[dict[str, Any]]:
    fit_by, archive = source_tables()
    out = []
    for s in SURFACES:
        f = fit_by[s]; a = archive[s]
        out.append({"surface_name": s, "catalyst": LABELS[s], "deltaG2_eV": float(f["NEB1_reaction_energy_EpoxideHHO_minus_Rpair"]),
                    "OOH_legacy_formation_eV": float(f["OOH_star_legacy"]), "actual_baseline_TOF_C6H12O": float(a["tof_C6H12O"]),
                    "actual_baseline_log10_TOF": float(a["log10_tof"]), "coordinate_semantics": "13 validated catalyst descriptor coordinates; background rate is fixed Ti-Co condition, not this catalyst-specific recalculation"})
    return out


def create_figures() -> dict[str, Any]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    fit_by, archive = source_tables()
    baseline_rows = read_csv(DATA / "baseline_validation.csv")
    fitted_rows = read_csv(DATA / "13site_BEP_fitted_TOF.csv")
    loo_rows = read_csv(DATA / "13site_BEP_LOOCV_TOF.csv")
    points = catalyst_points()
    write_csv(DATA / "grid_reference_catalyst_points.csv", points)
    fig_paths = {}
    # Figure 1: TOF ranks of actual CatMAP runs, not inferred from fit statistics.
    x = np.arange(len(SURFACES)); width = .25
    variants = [(baseline_rows, "Original DFT/NEB", "#334e68"), (fitted_rows, "BEP fitted (in-sample)", "#d17c27"), (loo_rows, "BEP leave-one-out", "#287c76")]
    fig, ax = plt.subplots(figsize=(14, 6.2))
    for k, (rows, label, color) in enumerate(variants):
        by = {r["surface_name"]: r for r in rows}
        ys = [float(by[s]["log10_TOF_C6H12O"]) if by[s].get("log10_TOF_C6H12O") else np.nan for s in SURFACES]
        ax.bar(x + (k - 1) * width, ys, width, label=label, color=color)
    ax.set_xticks(x, [LABELS[s] for s in SURFACES], rotation=35, ha="right")
    ax.set_ylabel("log10 TOF(C6H12O)")
    ax.set_title("CatMAP 0.3.1 steady-state TOF: original, full-fit BEP, and LOOCV BEP")
    ax.grid(axis="y", alpha=.2); ax.legend(frameon=False); fig.tight_layout()
    fig_paths.update(save_figure(fig, "01_13site_TOF_original_vs_BEP"))

    # Figure 2: actual signed TS1 differences reconstructed by the two models.
    predictions = read_csv(DATA / "BEP_TS1_predictions_and_LOOCV.csv")
    xx = np.asarray([float(r["actual_signed_TS1_minus_Rpair_eV"]) for r in predictions])
    yy1 = np.asarray([float(r["full_fit_predicted_signed_Y_eV"]) for r in predictions])
    yy2 = np.asarray([float(r["LOOCV_predicted_signed_Y_eV"]) for r in predictions])
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.3), sharex=True, sharey=True)
    for ax, yy, title in ((axes[0], yy1, "All-site fitted BEP"), (axes[1], yy2, "Leave-one-out prediction")):
        ax.scatter(xx, yy, s=62, c="#2f6690", edgecolors="white", linewidth=.7, zorder=3)
        lo = min(float(xx.min()), float(yy.min())); hi = max(float(xx.max()), float(yy.max()))
        ax.plot([lo, hi], [lo, hi], "--", color="#8b95a5", lw=1)
        ax.set_title(title); ax.set_xlabel("Actual signed G(TS1) − G(Rpair) / eV"); ax.grid(alpha=.18)
        for i, r in enumerate(predictions):
            ax.annotate(r["catalyst"], (xx[i], yy[i]), xytext=(4, 5), textcoords="offset points", fontsize=7)
    axes[0].set_ylabel("Predicted signed TS1 − Rpair / eV")
    fig.suptitle("BEP TS1 reconstruction: physical energy fit, independent CatMAP TOF evaluated separately")
    fig.tight_layout(); fig_paths.update(save_figure(fig, "02_TS1_actual_vs_predicted"))

    # Figure 3: only the official 61x61 CatMAP node results; pcolormesh displays discrete solver cells.
    grid_rows = read_csv(DATA / "grid_61x61.csv")
    n = 61
    node_by = {(int(r["grid_i"]), int(r["grid_j"])): r for r in grid_rows}
    if len(node_by) != n * n:
        raise RuntimeError("61x61 figure requested without a complete node-status table")
    xv, yv = grid_axes(n)
    zmap = np.full((n, n), np.nan); invalid = np.zeros((n, n), dtype=bool)
    for (i, j), r in node_by.items():
        if r["status"] in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"} and r.get("log10_TOF_C6H12O"):
            zmap[j, i] = float(r["log10_TOF_C6H12O"])
            invalid[j, i] = not (str(r.get("independent_TS_above_both_endpoints", "")).lower() == "true")
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    dx = xv[1] - xv[0]; dy = yv[1] - yv[0]
    xe = np.r_[xv[0] - dx/2, (xv[:-1] + xv[1:]) / 2, xv[-1] + dx/2]
    ye = np.r_[yv[0] - dy/2, (yv[:-1] + yv[1:]) / 2, yv[-1] + dy/2]
    mesh = ax.pcolormesh(xe, ye, np.ma.masked_invalid(zmap), cmap="viridis", shading="flat", rasterized=True)
    cb = fig.colorbar(mesh, ax=ax, pad=.02); cb.set_label("CatMAP log10 TOF(C6H12O)")
    if np.any(invalid):
        # Hatch endpoint-controlled/non-independent-TS1 nodes without hiding computed activity values.
        ax.contourf(xv, yv, invalid.astype(float), levels=[.5, 1.5], colors="none", hatches=["////"], alpha=0)
    ax.scatter([p["deltaG2_eV"] for p in points], [p["OOH_legacy_formation_eV"] for p in points],
               s=58, marker="o", facecolor="white", edgecolor="#111827", linewidth=.9,
               label="13 validated catalyst coordinates (background is fixed Ti-Co)", zorder=5)
    offsets = {
        "tife": (-35, -17), "timn": (10, -19), "tihf": (-45, 13), "tire": (12, -17),
        "timo": (12, 14), "tiv": (-37, -16), "tico": (-46, -18), "titi": (12, -12),
        "tiw": (12, 13), "tita": (15, 31), "ticr": (-42, 14), "ti": (-28, -18), "tini": (12, 13),
    }
    for p in points:
        ax.annotate(p["catalyst"], (p["deltaG2_eV"], p["OOH_legacy_formation_eV"]),
                    xytext=offsets[p["surface_name"]], textcoords="offset points", fontsize=8,
                    ha="left" if offsets[p["surface_name"]][0] >= 0 else "right", va="center",
                    arrowprops={"arrowstyle": "-", "color": "#555555", "lw": .55, "shrinkA": 1, "shrinkB": 4},
                    zorder=6)
    successful = [(i, j, r) for (i, j), r in node_by.items() if r["status"] in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"} and r.get("log10_TOF_C6H12O")]
    maximum = None
    valid_maximum = None
    if successful:
        i, j, maximum = max(successful, key=lambda v: float(v[2]["log10_TOF_C6H12O"]))
        ax.scatter([xv[i]], [yv[j]], marker="*", s=300, color="#ef476f", edgecolor="black", linewidth=.8, label="CatMAP grid maximum", zorder=8)
        boundary = i in (0, n-1) or j in (0, n-1)
        maximum["classification"] = "Maximum within sampled range (boundary)" if boundary else "interior grid maximum; verify local decline before calling a peak"
        maximum["grid_coordinate"] = [float(xv[i]), float(yv[j])]
        maximum["boundary"] = bool(boundary)
        valid_nodes = [(ii, jj, row) for (ii, jj), row in node_by.items()
                       if row["status"] in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}
                       and row.get("log10_TOF_C6H12O")
                       and str(row.get("independent_TS_above_both_endpoints", "")).lower() == "true"]
        if valid_nodes:
            vi, vj, valid_maximum = max(valid_nodes, key=lambda v: float(v[2]["log10_TOF_C6H12O"]))
            valid_boundary = vi in (0, n-1) or vj in (0, n-1)
            valid_maximum = dict(valid_maximum)
            valid_maximum["classification"] = "Maximum within independent-TS-valid sampled range (boundary)" if valid_boundary else "maximum within independent-TS-valid domain"
            valid_maximum["grid_coordinate"] = [float(xv[vi]), float(yv[vj])]
            valid_maximum["boundary"] = bool(valid_boundary)
            ax.scatter([xv[vi]], [yv[vj]], marker="D", s=88, color="#06d6a0", edgecolor="#073b4c", linewidth=.9,
                       label="Maximum with TS1 above both endpoints", zorder=9)
        peak_note = (
            f"Full-grid maximum (boundary): ΔG₂={xv[i]:.4f}, E(OOH*)={yv[j]:.4f} eV; "
            f"log₁₀TOF={float(maximum['log10_TOF_C6H12O']):.4f}"
        )
        if valid_maximum is not None:
            peak_note += (f"\nTS1-valid-domain maximum ({valid_maximum['classification'].split('(')[-1].rstrip(')')}): "
                          f"ΔG₂={valid_maximum['grid_coordinate'][0]:.4f}, E(OOH*)={valid_maximum['grid_coordinate'][1]:.4f} eV; "
                          f"log₁₀TOF={float(valid_maximum['log10_TOF_C6H12O']):.4f}")
        ax.text(.98, .58, peak_note, transform=ax.transAxes, ha="right", va="center", fontsize=8,
                bbox={"boxstyle": "round,pad=.35", "fc": "white", "ec": "#666", "alpha": .94}, zorder=10)
    fail_count = sum(r["status"] not in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"} for r in grid_rows)
    ax.set(xlabel="ΔG₂ = G(EpoxideHHO*) − G(Rpair*) / eV", ylabel="Historical E(OOH*) / eV",
           title=f"CatMAP 0.3.1 BEP–OOH conditional activity map: 61×61 actual nodes\nTi-Co fixed energy background; successes={len(successful)}, failed={fail_count}; hatched=TS1 not above both endpoints")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.12), ncol=3, fontsize=8, framealpha=.96)
    fig.tight_layout(rect=[0, .08, 1, .96])
    fig_paths.update(save_figure(fig, "03_CatMAP_BEP_OOH_dG2_activity_61x61"))
    meta = {"grid_nodes": len(grid_rows), "successful_nodes": len(successful), "failed_nodes": fail_count,
            "maximum": maximum, "maximum_with_independent_TS_valid": valid_maximum, "figures": fig_paths}
    write_json(DATA / "figure_and_grid_summary.json", meta)
    write_csv(DATA / "grid_peak_summary.csv", [
        {"domain": "full_grid", "classification": (maximum or {}).get("classification"),
         "grid_i": (maximum or {}).get("grid_i"), "grid_j": (maximum or {}).get("grid_j"),
         "deltaG2_eV": (maximum or {}).get("deltaG2_eV"), "OOH_legacy_formation_eV": (maximum or {}).get("OOH_legacy_formation_eV"),
         "log10_TOF_C6H12O": (maximum or {}).get("log10_TOF_C6H12O"),
         "TS1_above_both_endpoints": (maximum or {}).get("independent_TS_above_both_endpoints")},
        {"domain": "independent_TS_valid", "classification": (valid_maximum or {}).get("classification"),
         "grid_i": (valid_maximum or {}).get("grid_i"), "grid_j": (valid_maximum or {}).get("grid_j"),
         "deltaG2_eV": (valid_maximum or {}).get("deltaG2_eV"), "OOH_legacy_formation_eV": (valid_maximum or {}).get("OOH_legacy_formation_eV"),
         "log10_TOF_C6H12O": (valid_maximum or {}).get("log10_TOF_C6H12O"),
         "TS1_above_both_endpoints": (valid_maximum or {}).get("independent_TS_above_both_endpoints")},
    ])
    return meta


def save_figure(fig: Any, stem: str) -> dict[str, str]:
    out = {}
    for ext in ("png", "pdf", "svg"):
        p = FIGURES / f"{stem}.{ext}"
        fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight")
        out[ext] = str(p)
    import matplotlib.pyplot as plt
    plt.close(fig)
    return {stem: out}


def final_validate_and_report(fit: dict[str, Any], max_residual: str, started: float) -> None:
    grid = read_csv(DATA / "grid_61x61.csv")
    statuses = {s: sum(r["status"] == s for r in grid) for s in ("CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE", "FAILED", "FAILED_VALIDATION")}
    execution_sources = {}
    for row in grid:
        source = row.get("execution_source", "unknown")
        execution_sources[source] = execution_sources.get(source, 0) + 1
    if len(grid) != 3721:
        raise RuntimeError(f"Final grid row count is not 3721: {len(grid)}")
    success = sum(statuses[k] for k in ("CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"))
    successful_res = [float(r["residual"]) for r in grid if r["status"] in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"} and r.get("residual")]
    max_grid_residual = max(successful_res) if successful_res else math.nan
    figure_summary = json.loads((DATA / "figure_and_grid_summary.json").read_text(encoding="utf-8"))
    valid_peak = figure_summary.get("maximum_with_independent_TS_valid") or {}
    baseline = read_csv(DATA / "baseline_validation.csv")
    full = read_csv(DATA / "13site_BEP_fitted_TOF.csv")
    loo = read_csv(DATA / "13site_BEP_LOOCV_TOF.csv")
    bmap = {r["surface_name"]: r for r in full}; lmap = {r["surface_name"]: r for r in loo}; omap = {r["surface_name"]: r for r in baseline}
    tof_compare = []
    for surface in ("ti", "tico"):
        tof_compare.append({"surface_name": surface, "catalyst": LABELS[surface], "original_TOF": omap[surface]["TOF_C6H12O"],
                            "BEP_fitted_TOF": bmap[surface]["TOF_C6H12O"], "BEP_LOOCV_TOF": lmap[surface]["TOF_C6H12O"],
                            "original_log10_TOF": omap[surface]["log10_TOF_C6H12O"], "BEP_fitted_log10_TOF": bmap[surface]["log10_TOF_C6H12O"],
                            "BEP_LOOCV_log10_TOF": lmap[surface]["log10_TOF_C6H12O"]})
    write_csv(DATA / "Ti_TiCo_TOF_comparison.csv", tof_compare)
    seconds = time.time() - started
    maxnode = figure_summary.get("maximum") or {}
    timing_path = DATA / "official_compute_timing.json"
    timing_record = json.loads(timing_path.read_text(encoding="utf-8")) if timing_path.is_file() else {}
    grid_solver_seconds = float(timing_record.get("grid_61x61_solver_stage_seconds", seconds))
    report = f"""# CatMAP 0.3.1 四参数 OOH* 扩展BEP：完整计算结果

## 执行口径

- 本机CatMAP：0.3.1；实际调用 `ReactionModel.run(recalculate=True)`，模型配置为 `ThermodynamicScaler` + `SteadyStateSolver`，不修改CatMAP安装源码。
- 13个基线逐位点使用归档严格 `.mkm`/`energies.txt`的隔离快照；333.15 K、单有效中心、原六步反应、effective reservoir concentration（H2O2=0.185、C6H12=1.0、H2O=0.815、C6H12O=1e-20）、frozen气/吸附热力学、mpmath {75}位、残差容差{TOLERANCE}。
- OOH输入为历史单吸附态形成能，仅作外部BEP描述符；不等同六步网络桥连OOH，也没有加进网络。
- exact fit `R²={fit['R2']:.9f}`, RMSE={fit['RMSE_eV']:.9f} eV, MAE={fit['MAE_eV']:.9f} eV, LOOCV RMSE={fit['LOOCV_RMSE_eV']:.9f} eV。系数/中心见 `data/exact_BEP_fit.json`；LOOCV每折只用12个位点并重算中心。

## CatMAP执行计数与验收

- 原始基线：{len(baseline)}个位点；归档TOF最大绝对差 log10={max(float(r['absolute_log10_difference']) for r in baseline):.4e} dex。
- BEP全样本拟合重算：13个位点；LOOCV逐点重算：13个位点。
- 网格61×61：3721个坐标，CatMAP求解状态计数 `{json.dumps(statuses, ensure_ascii=False)}`；结果来源阶段计数 `{json.dumps(execution_sources, ensure_ascii=False)}`；成功={success}；最大残差={max_grid_residual:.6e}；容差={TOLERANCE}。其中来自61×61正式阶段的3696个节点均有该阶段逐节点记录；其余与5×5节点坐标及输入哈希完全一致的25个节点复用先前已验证CatMAP结果。当前报告重导出时3721行均由精确输入哈希缓存读取，见`retrieval_status`；这不改变`status`所记录的真实求解结果。
- 5×5和15×15阶段数据分别保存为 `data/grid_5x5.csv`、`data/grid_15x15.csv`。所有失败坐标不填零/不插值，日志见 `runs/`。
- 61×61正式网格求解阶段耗时={grid_solver_seconds:.2f} s（原始计算调用计时，不含后续缓存复用的图件重生成）；本次报告/图件重生成调用耗时={seconds:.2f} s。此前输入准备、基线、BEP和5×5/15×15阶段为独立调用，未合并计入该值。全域最大值分类=`{maxnode.get('classification','not available')}`；TS1高于两个端点的子域最大值分类=`{valid_peak.get('classification','not available')}`。
- TS1高于初/末态的网格点={sum(str(r.get('independent_TS_above_both_endpoints','')).lower() == 'true' for r in grid)}/{len(grid)}；端点控制节点不删去，图中以斜线标识。
- 每个节点均使用运行前绑定的 `catmap_case.mkm` 指向本节点 `energies.txt`；解析后逐项比对Rpair/TS1/EpoxideHHO/HHO/TS2/OHHa/H2O共7态。本工作流修复前生成的旧哈希运行/缓存已保留作诊断，但不进入当前CSV、峰值或图件；`RUN_SCHEMA`变化使其不能作为当前结果缓存。

## 网格最大值坐标

- 全部网格最大值：ΔG₂={maxnode.get('deltaG2_eV')} eV，E(OOH*)={maxnode.get('OOH_legacy_formation_eV')} eV，log10TOF={maxnode.get('log10_TOF_C6H12O')}；状态=`{maxnode.get('classification')}`，TS1独立过渡态有效=`{maxnode.get('independent_TS_above_both_endpoints')}`。
- TS1高于初态和末态的子域最大值：ΔG₂={valid_peak.get('deltaG2_eV')} eV，E(OOH*)={valid_peak.get('OOH_legacy_formation_eV')} eV，log10TOF={valid_peak.get('log10_TOF_C6H12O')}；状态=`{valid_peak.get('classification')}`。
- 两者均为边界最大值；本网格没有证据支持把最大值称为内部火山峰。

## TS1物理有效性与CatMAP端点处理

保留所有BE预测的有符号 TS1−Rpair 和预测TS1原值。独立过渡态有效条件是 `Y1 >= max(0, ΔG2)`。本机 CatMAP 0.3.1官方 `catmap/data/templates.py` 的 `rate_constants` 实现以 `max(G_IS, G_FS, G_TS)` 作为有效势垒顶，因此当预测TS低于任一端点，实际速率由较高端点控制；这不是输入能量被悄悄改写。本结果分别保留 `independent_TS_above_both_endpoints`、`CatMAP_endpoint_clamp_applied` 和 `CatMAP_effective_TS_control`。

## 条件网格限制与科学解释

网格保持Ti-Co的Rpair、HHO、OHHa、H2O、TS2和气相储库固定，只随ΔG2改变EpoxideHHO，并用历史单吸附态OOH值经四参数扩展BEP预测TS1。由此得到的是Ti-Co条件背景下的CatMAP网格，不是13种真实催化剂分别重算的排名。OOH历史参考与严格六步形成能的slab/热修正同源性仍未完成物理闭合；拟合R²不能替代可迁移验证，LOOCV误差已并列给出。BEP点位重建造成的TOF变化及网格峰，均不能单独证明双核协同增强或普适火山规律。

## 关键输出

- `data/baseline_validation.csv`、`data/13site_original_TOF.csv`、`data/13site_BEP_fitted_TOF.csv`、`data/13site_BEP_LOOCV_TOF.csv`
- `data/BEP_TS1_predictions_and_LOOCV.csv`、`data/grid_61x61.csv`、`data/grid_61x61_coverages.csv`、`data/grid_61x61_net_rates.csv`
- `figures/`下3类PNG/PDF/SVG图；输入SHA256见 `data/SHA256SUMS.csv`。

## Ti和Ti-Co实际TOF

|位点|原始CatMAP TOF|全样本BEP TS1 CatMAP TOF|LOOCV-BEP TS1 CatMAP TOF|
|---|---:|---:|---:|
""" + "\n".join(f"|{r['catalyst']}|{float(r['original_TOF']):.12e}|{float(r['BEP_fitted_TOF']):.12e}|{float(r['BEP_LOOCV_TOF']):.12e}|" for r in tof_compare) + "\n"
    REPORTS.mkdir(parents=True, exist_ok=True)
    (REPORTS / "FINAL_BEP_OOH_ACTIVITY_MAP_REPORT.md").write_text(report, encoding="utf-8")
    write_json(DATA / "execution_summary.json", {"catmap_version": "0.3.1", "run_schema": RUN_SCHEMA, "grid_total": len(grid), "grid_status_counts": statuses, "grid_execution_source_counts": execution_sources,
                                                 "baseline_max_abs_dex": max(float(r["absolute_log10_difference"]) for r in baseline),
                                                 "maximum_grid_residual": max_grid_residual, "elapsed_seconds": seconds,
                                                 "grid_61x61_solver_stage_seconds": grid_solver_seconds,
                                                 "elapsed_scope": "current invocation; true solver-stage time separately recorded",
                                                 "grid_maximum": maxnode, "grid_maximum_independent_TS_valid": valid_peak,
                                                 "max_residual_all_runs": max_residual})


def sha_manifest() -> None:
    root = HERE
    rows = []
    for p in sorted(root.rglob("*")):
        if not p.is_file() or "runs" in p.parts or "cases" in p.parts or p.name == "SHA256SUMS.csv":
            continue
        rows.append({"path": str(p.relative_to(root)).replace("\\", "/"), "sha256": sha256_file(p), "size_bytes": p.stat().st_size})
    write_csv(DATA / "SHA256SUMS.csv", rows, ["path", "sha256", "size_bytes"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["prepare", "baseline", "fit", "grid", "all"], default="all")
    parser.add_argument("--resolution", type=int, choices=[5, 15, 61], default=5)
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()
    for p in (SNAP, CASES, RUNS, DATA, FIGURES, REPORTS):
        p.mkdir(parents=True, exist_ok=True)
    started = time.time()
    provenance = prepare_snapshots()
    fit_by, archived, validation = validate_snapshots()
    source_rows = [fit_by[s] for s in SURFACES]
    fit = fit_bep(source_rows)
    print(f"CatMAP=0.3.1 INPUTS=PASS sites=13 BEP_R2={fit['R2']:.10f} LOOCV_RMSE_eV={fit['LOOCV_RMSE_eV']:.8f}", flush=True)
    write_json(DATA / "exact_BEP_fit.json", fit)
    if args.stage in {"prepare"}:
        sha_manifest(); return
    baseline_path = DATA / "baseline_validation.csv"
    if args.stage in {"baseline", "all"}:
        baseline_results = run_baselines(fit_by, archived)
        print(f"BASELINE_GATE=PASS sites={len(baseline_results)}", flush=True)
        if args.stage == "baseline":
            sha_manifest(); return
    else:
        if not baseline_path.exists():
            raise RuntimeError("BEP fitting/optimization blocked: baseline_validation.csv absent; run --stage baseline first")
        rows = read_csv(baseline_path)
        if len(rows) != 13 or any(r["status"] != "CONVERGED" for r in rows):
            raise RuntimeError("BEP fitting/optimization blocked: original CatMAP baseline gate did not pass")
    if args.stage in {"fit", "all"}:
        full_results, loo_results, fit = run_fit_scenarios(fit_by)
        print(f"BEP_CATMAP=PASS full={sum(r['status']=='CONVERGED' for r in full_results)}/13 loocv={sum(r['status']=='CONVERGED' for r in loo_results)}/13", flush=True)
        if args.stage == "fit":
            sha_manifest(); return
    else:
        for required in (DATA / "13site_BEP_fitted_TOF.csv", DATA / "13site_BEP_LOOCV_TOF.csv"):
            if not required.exists():
                raise RuntimeError("Grid blocked: fitted and LOOCV CatMAP site results are required")
    resolutions = [5, 15, 61] if args.stage == "all" else [args.resolution]
    for res in resolutions:
        nodes = run_grid(res)
        failures = [r for r in nodes if r.get("status") not in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"}]
        print(f"GRID_STAGE={res}x{res} total={len(nodes)} successful={len(nodes)-len(failures)} failed={len(failures)}", flush=True)
        if failures:
            raise RuntimeError(f"Grid completeness gate failed at {res}x{res}; larger grid stages are stopped. See per-node JSON/logs and grid CSV.")
    if 61 in resolutions:
        create_figures()
        max_resid = "0"
        for path in RUNS.rglob("official_result.json"):
            try:
                r = json.loads(path.read_text(encoding="utf-8"))
                if r.get("status") in {"CONVERGED", "CONVERGED_ADAPTIVE", "REUSED_CACHE"} and r.get("residual"):
                    max_resid = max(max_resid, r["residual"], key=lambda x: float(x))
            except (ValueError, OSError):
                pass
        final_validate_and_report(fit, max_resid, started)
    sha_manifest()
    print("SHA256_MANIFEST=", DATA / "SHA256SUMS.csv", flush=True)


if __name__ == "__main__":
    main()
