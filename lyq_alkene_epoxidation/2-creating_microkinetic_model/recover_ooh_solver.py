# -*- coding: utf-8 -*-
"""CatMAP MinResidMapper recovery with *validated* warm starts.

This is a numerical-method change only. The original 14-surface DFT input,
eight elementary reactions, effective activities, and target temperature
333.15 K are never changed. Official mechanisms used:
- CatMAP MinResidMapper neighbor guesses and descriptor-line bisection.
- numbers_map/coverage_map as native CatMAP recalculation seeds.
- ReactionModel.run(recalculate=True).

IMPORTANT: the source point of any cached seed MUST lie on the new
descriptor grid, otherwise MinResidMapper ignores the off-grid guess.
Results become publishable numerical candidates only if all 8 step fluxes
and coverage residuals pass at 333.15 K at two precision levels AND at
least two independent, vetted source temperatures yield matching TOFs.

Usage:
  python run_ooh_analysis.py --recover --recover-surfaces titi timn ticr
  python run_ooh_analysis.py --recover
"""
from __future__ import annotations

import math
import os
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd

import run_ooh_analysis as core


# v2 output directory preserves all earlier recovery/ CSVs and logs.
RECOVERY_ROOT = core.OUTPUT / "recovery_v2"
TARGET_TEMPERATURE = 333.15
TARGET_PRESSURE = 1.0
SOURCE_MIN_SEPARATION_K = 0.01
MAX_SOURCE_SEEDS = 3
CROSS_SEED_LOG10_TOL = core.TOF_STABILITY_LOG10_TOL
REFINED_STAGES = ((260, "1e-180"), (360, "1e-260"), (460, "1e-340"))


def _boolean(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def _target_present(coverage_map) -> bool:
    return any(
        abs(float(pt[0]) - TARGET_TEMPERATURE) < 1e-5
        and abs(float(pt[1]) - TARGET_PRESSURE) < 1e-7
        for pt, _ in (coverage_map or [])
    )


def _find_point(point, entries, temperature_tol=1e-6):
    for candidate_pt, value in entries or []:
        if (abs(float(candidate_pt[0]) - float(point[0])) < temperature_tol
                and abs(float(candidate_pt[1]) - float(point[1])) < 1e-8):
            return value
    return None


def _unique_temperature_candidates(records: list[dict],
                                   max_seeds: int = MAX_SOURCE_SEEDS) -> list[dict]:
    """Pick nearest/middle/farthest *independently checked* temperature seeds.

    Never select a point based on 'mapper_success' alone.
    """
    eligible = [
        r for r in records
        if r.get("source_quality_pass") is True
        and math.isfinite(float(r.get("source_temperature_K", float("nan"))))
        and float(r["source_temperature_K"]) > TARGET_TEMPERATURE + SOURCE_MIN_SEPARATION_K
        and r.get("seed_file")
    ]
    eligible.sort(key=lambda r: (float(r["source_temperature_K"]),
                                 str(r["seed_file"])))
    distinct = []
    for item in eligible:
        if any(abs(float(item["source_temperature_K"])
                   -float(old["source_temperature_K"])) < SOURCE_MIN_SEPARATION_K
               for old in distinct):
            continue
        distinct.append(item)
    if len(distinct) <= max_seeds:
        return distinct
    indices = {0, len(distinct)-1}
    if max_seeds >= 3:
        indices.add((len(distinct)-1)//2)
    # Deterministic, evenly spaced selection for max_seeds>3.
    for idx in np.linspace(0, len(distinct)-1, max_seeds).astype(int):
        if len(indices) >= max_seeds:
            break
        indices.add(int(idx))
    return [distinct[i] for i in sorted(indices)[:max_seeds]]


def _compare_independent_seeds(valid: list[dict],
                               limit: float = CROSS_SEED_LOG10_TOL) -> dict:
    """Disallow acceptance of different stable roots from different seeds."""
    values = [float(v["log10_net_C6H12O"]) for v in valid
              if (v.get("quality_pass") is True
                  and v.get("refinement_stable") is True
                  and math.isfinite(float(v.get("log10_net_C6H12O", float("nan")))))]
    delta = max(values) - min(values) if values else float("nan")
    return {
        "independent_seed_validated_count": len(values),
        "independent_seed_delta_log10_TOF": float(delta),
        "independent_seed_agreement": (
            len(values) >= 2 and math.isfinite(delta) and delta <= limit
        ),
    }


def _store_single_verified_seed(point, coverages, numbers, path: Path) -> None:
    """Write ONLY native CatMAP solution maps; no stale output/metadata.

    One entry per map avoids any sorting mismatch between coverage_map and
    numbers_map in CatMAP MinResidMapper's initial-map handling.
    """
    if numbers is None:
        raise ValueError("A numbers_solver warm start requires numbers_map")
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "coverage_map": [[list(point), coverages]],
        "numbers_map": [[list(point), numbers]],
    }
    with path.open("wb") as fh:
        pickle.dump(payload, fh, protocol=2)


def _inspect_bridge_point(model, point) -> dict:
    """Independently validate a CatMAP map entry at its OWN temperature.

    The map being populated only proves solver acceptance. This function
    checks the eight real net elementary rates and the independent coverage
    residual at the same descriptor point.
    """
    from mpmath import mp

    model.mapper.get_point_output(point)
    cover_names = list(model.output_labels["coverage"])
    coverage = list(model._coverage)
    if len(coverage) == len(cover_names):
        coverage.append(mp.mpf(1) - mp.fsum(coverage))
    if len(coverage) != len(cover_names) + 1:
        raise ValueError("Inconsistent adsorbate and vacant-site coverage vector")

    gas_labels = list(model.output_labels["turnover_frequency"])
    gas = dict(zip(gas_labels, model._turnover_frequency))
    rates = list(model._rate)
    model.solver._descriptors = list(point)
    residual = model.solver.get_residual(
        coverage, validate_coverages=False, refresh_rate_constants=True
    )
    return core.check_steady_state(
        rates, gas, coverage, residual, numbers_solver=True
    )


def make_temperature_bridge(
    surface: str, high_temperature: float,
    steps: int = 17, precision: int = 260, tolerance: str = "1e-180",
) -> tuple[dict, list[dict]]:
    """Build + independently inspect a native CatMAP temperature coverage map."""
    if not math.isfinite(high_temperature) or high_temperature <= TARGET_TEMPERATURE:
        raise ValueError("Bridge high temperature must be finite and exceed 333.15 K")
    if steps < 3:
        raise ValueError("At least three temperature path points are required")
    if surface not in core.SURFACES:
        raise ValueError(f"Unknown DFT catalyst surface: {surface}")

    from catmap import ReactionModel

    work = RECOVERY_ROOT / "bridges" / surface / f"T{high_temperature:g}"
    work.mkdir(parents=True, exist_ok=True)
    config = core.create_single_surface_setup(
        surface, work, precision=precision, tolerance=tolerance,
        numbers_solver=True
    )
    with config.open("a", encoding="utf-8") as fh:
        fh.write("\n# Numerical-only temperature path; final target is 333.15 K.\n")
        fh.write(f"descriptor_ranges = [[333.15, {high_temperature!r}], [1.0, 1.0]]\n")
        fh.write(f"resolution = [{steps}, 1]\n")
        fh.write("max_bisections = 5\n")
        fh.write("max_rootfinding_iterations = 350\n")
    for suffix in (".log", ".pkl"):
        candidate = config.with_suffix(suffix)
        if candidate.is_file():
            candidate.unlink()

    previous = Path.cwd()
    try:
        os.chdir(work)
        model = ReactionModel(setup_file=config.name)
        model.output_variables = ["coverage", "rate", "turnover_frequency"]
        model.run()
        point_map = list(getattr(model, "coverage_map", []) or [])
        numbers_map = list(getattr(model, "numbers_map", []) or [])
        diagnostics = []
        vetted_count = 0

        for point, coverages in point_map:
            temperature = float(point[0])
            source_id = f"T{int(round(temperature*100000)):09d}"
            record = {
                "surface": surface,
                "bridge_temperature_K": float(high_temperature),
                "source_temperature_K": temperature,
                "source_pressure": float(point[1]),
                "source_quality_pass": False,
                "seed_file": "",
            }
            try:
                quality = _inspect_bridge_point(model, point)
                record.update({f"source_{k}": v for k, v in quality.items()
                               if k in ("quality_status", "cycle_max_relative_error",
                                        "balance_max_relative_error",
                                        "steady_state_residual_relative_to_max_rate",
                                        "coverage_sum_error")})
                is_eligible = (
                    quality["quality_pass"]
                    and temperature > TARGET_TEMPERATURE + SOURCE_MIN_SEPARATION_K
                )
                if is_eligible:
                    nums = _find_point(point, numbers_map)
                    if nums is None:
                        raise ValueError("Bridge has coverages but lacks corresponding numbers")
                    # Use the map's own solution, never an unverified extrapolation.
                    path = (RECOVERY_ROOT / "seeds" / surface
                            / f"T{high_temperature:g}_{source_id}.pkl")
                    _store_single_verified_seed(point, coverages, nums, path)
                    record["seed_file"] = str(path.resolve())
                    record["source_quality_pass"] = True
                    vetted_count += 1
            except Exception as exc:
                record["source_quality_status"] = "inspection_failed"
                record["source_error"] = f"{type(exc).__name__}: {exc}"
            diagnostics.append(record)

        meta = {
            "surface": surface,
            "bridge_temperature": float(high_temperature),
            "bridge_point_count": len(point_map),
            "bridge_numbers_point_count": len(numbers_map),
            "bridge_includes_target": _target_present(point_map),
            "bridge_independently_verified_seed_count": vetted_count,
            "bridge_status": "mapped_and_vetted" if vetted_count
                             else "mapped_but_no_vetted_seed",
        }
        return meta, diagnostics
    finally:
        os.chdir(previous)


def _refine_with_verified_seed(surface: str, seed: dict
                               ) -> tuple[dict, list[dict]]:
    """At every precision stage restart from SAME vetted high-T source.

    NEVER use a target-point result which failed independent flux checks
    as a seed for the next stage.
    """
    temperature = float(seed["source_temperature_K"])
    seed_file = Path(seed["seed_file"]).resolve()
    trial_id = "sourceT" + str(int(round(temperature * 100000)))
    prev = None
    attempts = []
    accepted = None
    for i, (precision, tol) in enumerate(REFINED_STAGES, 1):
        root = RECOVERY_ROOT / "target" / surface / trial_id / f"stage_{i:02d}"
        try:
            row = core.run_single_surface(
                surface, precision=precision, tolerance=tol,
                output_root=root, seed_file=seed_file,
                seed_temperature=temperature,
                max_iterations=400, max_bisections=5,
            )
        except Exception as exc:
            row = {
                "surface": surface,
                "status": f"FAILED: {type(exc).__name__}: {exc}",
                "quality_status": "solver_failed",
                "quality_pass": False,
            }
        row["source_temperature_K"] = temperature
        row["seed_file"] = str(seed_file)
        row["recovery_stage"] = i
        row["refinement_stable"] = False
        row["refinement_delta_log10_TOF"] = float("nan")
        if prev is not None:
            row.update(core.compare_refinement_runs(prev, row))
        attempts.append(dict(row))
        if row["refinement_stable"]:
            accepted = dict(row)
            accepted["quality_status"] = "validated_two_precisions"
            accepted["quality_pass"] = True
            break
        prev = row

    if accepted is not None:
        return accepted, attempts
    last = dict(attempts[-1])
    last["quality_pass"] = False
    last["refinement_stable"] = False
    last["quality_status"] = "not_stable_from_this_seed"
    return last, attempts


def recover_one_surface(
    surface: str, bridge_temperatures: list[float],
    steps: int = 17,
) -> tuple[dict, list[dict], list[dict], list[dict]]:
    """Check mapped points, test diverse verified seeds, require agreement."""
    bridge_rows, map_point_rows, seed_rows, candidates = [], [], [], []
    for high_temperature in bridge_temperatures:
        try:
            meta, inspected = make_temperature_bridge(
                surface, high_temperature, steps=steps
            )
            bridge_rows.append(meta)
            map_point_rows.extend(inspected)
        except Exception as err:
            bridge_rows.append({
                "surface": surface, "bridge_temperature": high_temperature,
                "bridge_status": f"FAILED: {type(err).__name__}: {err}",
            })

    selected = _unique_temperature_candidates(map_point_rows, MAX_SOURCE_SEEDS)
    if not selected:
        return (
            {"surface": surface, "status": "no validated map point",
             "quality_status": "no_vetted_source_seed",
             "quality_pass": False, "refinement_stable": False},
            [], bridge_rows, map_point_rows
        )

    for seed in selected:
        res, attempts = _refine_with_verified_seed(surface, seed)
        candidates.append(res)
        seed_rows.extend(attempts)

    valid = [r for r in candidates
             if r.get("quality_pass") is True and
             r.get("refinement_stable") is True]
    agreement = _compare_independent_seeds(valid)
    result = dict(valid[0]) if valid else {"surface": surface}
    result.update(agreement)
    result["vetted_source_seed_count"] = len(selected)
    result["attempted_source_temperatures_K"] = ";".join(
        f"{float(x['source_temperature_K']):.5f}" for x in selected
    )
    if agreement["independent_seed_agreement"]:
        result["quality_pass"] = True
        result["refinement_stable"] = True
        result["quality_status"] = "validated_cross_seed_and_precision"
        result["status"] = "solved"
    else:
        # A single converged root is not enough to demonstrate
        # independent-start convergence; differing roots may represent
        # branch multiplicity. Both cases need user inspection.
        result["quality_pass"] = False
        result["refinement_stable"] = False
        result["quality_status"] = (
            "different_valid_steady_state_branches"
            if len(valid) >= 2
            else "insufficient_independent_valid_roots"
        )
        result.setdefault("status", "not independently validated")
    return result, seed_rows, bridge_rows, map_point_rows


def recover_unvalidated_surfaces(
    surfaces: list[str] | None = None,
    bridge_temperatures: list[float] | None = None,
    steps: int = 17,
) -> None:
    """Update previous baseline ONLY after cross-seed AND cross-precision QC."""
    base_file = core.OUTPUT / "dft_baseline.csv"
    if not base_file.exists():
        raise FileNotFoundError("Please run baseline CatMAP analysis first")
    baseline = pd.read_csv(base_file)
    existing_good = _boolean(baseline["quality_pass"]) & _boolean(
        baseline["refinement_stable"]
    )
    if not surfaces:
        selected = baseline.loc[~existing_good, "surface"].tolist()
    else:
        selected = list(dict.fromkeys(surfaces))
        bad = set(selected) - set(core.SURFACES)
        if bad:
            raise ValueError(f"Unknown catalyst surface names: {sorted(bad)}")
        previous_good = set(baseline.loc[existing_good, "surface"])
        selected = [x for x in selected if x not in previous_good]

    if not selected:
        print("All selected surfaces already meet the existing quality standard.")
        return

    bridge_temperatures = list(bridge_temperatures or (550.0, 750.0))
    RECOVERY_ROOT.mkdir(parents=True, exist_ok=True)
    all_results, all_trials, all_bridges, all_points, updated = [], [], [], [], []
    for surface in selected:
        print(f"\n====== VETTED CatMAP MinResidMapper recovery: {surface} ======")
        try:
            final, trials, bridges, points = recover_one_surface(
                surface, bridge_temperatures=bridge_temperatures, steps=steps
            )
        except Exception as err:
            final = {"surface": surface, "status": f"{type(err).__name__}: {err}",
                     "quality_status": "recovery_exception", "quality_pass": False,
                     "refinement_stable": False}
            trials, bridges, points = [], [], []
        all_results.append(final)
        all_trials.extend(trials)
        all_bridges.extend(bridges)
        all_points.extend(points)
        print(f"[{surface}] {final['quality_status']}; "
              f"validated source seeds: {final.get('vetted_source_seed_count',0)}")

        if final.get("quality_pass") is True and final.get("refinement_stable") is True:
            updated.append(surface)
            for key, value in final.items():
                if key not in baseline.columns:
                    baseline[key] = np.nan
                baseline.loc[baseline["surface"].eq(surface), key] = value

    # Diagnostics are always written; they never corrupt validated baseline.
    pd.DataFrame(all_bridges).to_csv(
        RECOVERY_ROOT / "bridge_diagnostics.csv", index=False
    )
    pd.DataFrame(all_points).to_csv(
        RECOVERY_ROOT / "bridge_point_quality.csv", index=False
    )
    pd.DataFrame(all_trials).to_csv(
        RECOVERY_ROOT / "seed_trials.csv", index=False
    )
    pd.DataFrame(all_results).to_csv(
        RECOVERY_ROOT / "recovery_results.csv", index=False
    )
    if updated:
        backup = RECOVERY_ROOT / "baseline_before_recovery.csv"
        if not backup.exists():
            import shutil
            shutil.copy2(base_file, backup)
        baseline.to_csv(base_file, index=False)
        df = core.read_data().merge(baseline, on="surface", how="left")
        good = (_boolean(df["quality_pass"])
                & _boolean(df["refinement_stable"]))
        log_tof = pd.to_numeric(df["log10_net_C6H12O"], errors="coerce")
        df["eligible_for_tof_fit"] = good & np.isfinite(log_tof)
        df.to_csv(core.OUTPUT / "OOH_TOF_results.csv", index=False)
        fit_data = df.loc[df["eligible_for_tof_fit"]]
        core.plot_linear(
            fit_data, "G_OOH_eV", "log10_net_C6H12O",
            "OOH_logTOF_linear.png", "OOH* formation energy (eV)",
            "log10(net epoxide TOF)", min_points=5
        )
    print("\nRecovered validated surfaces:", len(updated), updated)
    print("Bridge candidate checks:", RECOVERY_ROOT / "bridge_point_quality.csv")
    print("Results and attempts:", RECOVERY_ROOT / "recovery_results.csv")
    print("Only 333.15 K is used for final TOF; all DFT energies remain unchanged.")


if __name__ == "__main__":
    recover_unvalidated_surfaces()
