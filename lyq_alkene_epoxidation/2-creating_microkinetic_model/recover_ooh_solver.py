# -*- coding: utf-8 -*-
"""CatMAP official MinResidMapper coverage/number-map warm-start recovery.

Runs only the not-yet-validated surfaces using *unchanged* DFT energies and
reaction network. Numerical homotopy spans temperature exclusively as a
path-finding device; the reported final TOF is always at 333.15 K.

Scientific basis:
https://catmap.readthedocs.io/en/latest/tutorials/refining_a_microkinetic_model.html
https://catmap.readthedocs.io/en/latest/_modules/catmap/mappers/min_resid_mapper.html

Usage: python run_ooh_analysis.py --recover
       python run_ooh_analysis.py --recover --recover-surfaces titi timn ticr
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

import run_ooh_analysis as core


RECOVERY_ROOT = core.OUTPUT / "recovery"
REFINED_STAGES = (
    (260, "1e-180"),
    (360, "1e-260"),
    (460, "1e-340"),
)


def _target_present(coverage_map) -> bool:
    return any(
        abs(float(point[0]) - 333.15) < 1e-5
        and abs(float(point[1]) - 1.0) < 1e-7
        for point, _ in (coverage_map or [])
    )


def make_temperature_bridge(
    surface: str, high_temperature: float,
    steps: int = 17, precision: int = 260, tolerance: str = "1e-180",
) -> dict:
    """Create a CatMAP MinResidMapper temperature map with reusable numbers_map.

    Descending temperature is handled by the mapper's built-in reverse traversal
    and nearest-neighbor guesses. No changes to actual reaction conditions
    survive in the final single-point model.
    """
    if high_temperature <= 333.15:
        raise ValueError("Bridge temperature must exceed 333.15 K.")
    if steps < 3:
        raise ValueError("At least three temperature path points are needed.")

    from catmap import ReactionModel

    work = RECOVERY_ROOT / "bridges" / surface / f"T{high_temperature:g}"
    work.mkdir(parents=True, exist_ok=True)
    config = core.create_single_surface_setup(
        surface, work, precision=precision, tolerance=tolerance,
        numbers_solver=True
    )
    with config.open("a", encoding="utf-8") as fh:
        fh.write("\n# Numerical temperature continuation only; not an experiment.\n")
        fh.write(f"descriptor_ranges = [[333.15, {high_temperature!r}], [1.0, 1.0]]\n")
        fh.write(f"resolution = [{steps}, 1]\n")
        fh.write("max_bisections = 5\n")
        fh.write("max_rootfinding_iterations = 350\n")

    cache = config.with_suffix(".pkl")
    if cache.exists():
        cache.unlink()
    logfile = config.with_suffix(".log")
    if logfile.exists():
        logfile.unlink()

    original_cwd = Path.cwd()
    try:
        os.chdir(work)
        model = ReactionModel(setup_file=config.name)
        # The MinResidMapper still computes and saves coverage_map + numbers_map.
        model.output_variables = ["coverage"]
        model.run()
        mapped = getattr(model, "coverage_map", None) or []
        seed_map = getattr(model, "numbers_map", None) or []
        mapped_count = len(mapped)
        seed_count = len(seed_map)
        target_mapped = _target_present(mapped)
        best_temperature = max((float(pt[0]) for pt, _ in mapped), default=float("nan"))
        if not cache.exists() or mapped_count == 0 or seed_count == 0:
            raise RuntimeError("CatMAP found no reusable numbers-map seed on bridge")
        return {
            "surface": surface,
            "bridge_temperature": high_temperature,
            "bridge_point_count": mapped_count,
            "bridge_numbers_point_count": seed_count,
            "bridge_includes_target": target_mapped,
            "bridge_highest_converged_temperature": best_temperature,
            "bridge_status": "available",
            "bridge_seed_file": str(cache.resolve()),
        }
    finally:
        os.chdir(original_cwd)


def recover_one_surface(
    surface: str, bridge_temperatures: list[float],
    steps: int = 17,
) -> tuple[dict, list[dict], list[dict]]:
    """Warm-start from a temperature map, then independently refine at T=333.15K."""
    bridge_logs = []
    selected_bridge = None
    for high_t in bridge_temperatures:
        try:
            entry = make_temperature_bridge(surface, high_t, steps=steps)
            bridge_logs.append({k: v for k, v in entry.items()
                                if k != "bridge_seed_file"})
            # Prefer a seed map which already includes the final target.
            if selected_bridge is None or entry["bridge_includes_target"]:
                selected_bridge = entry
            if entry["bridge_includes_target"]:
                break
        except Exception as err:
            bridge_logs.append({
                "surface": surface, "bridge_temperature": high_t,
                "bridge_status": f"FAILED: {type(err).__name__}: {err}",
            })
        if selected_bridge is not None:
            # A partial map still provides valid higher-T initial guesses;
            # retry with a larger high-temperature path only if needed.
            continue

    if selected_bridge is None:
        return (
            {"surface": surface, "quality_status": "no_reusable_temperature_bridge",
             "quality_pass": False, "refinement_stable": False,
             "status": "no temperature bridge converged"},
            [], bridge_logs
        )

    latest_seed = Path(selected_bridge["bridge_seed_file"])
    previous = None
    attempts = []
    accepted = None

    for stage, (prec, tol) in enumerate(REFINED_STAGES, 1):
        output_root = RECOVERY_ROOT / "target" / f"stage_{stage:02d}"
        try:
            row = core.run_single_surface(
                surface, precision=prec, tolerance=tol,
                output_root=output_root, seed_file=latest_seed,
                max_iterations=400, max_bisections=5,
            )
        except Exception as err:
            row = {
                "surface": surface, "quality_status": "solver_failed",
                "quality_pass": False, "refinement_stable": False,
                "status": f"FAILED: {type(err).__name__}: {err}",
            }

        row["recovery_stage"] = stage
        row["recovery_bridge_temperature"] = selected_bridge["bridge_temperature"]
        row["recovery_seed_origin"] = str(latest_seed)
        row["refinement_stable"] = False
        row["refinement_delta_log10_TOF"] = float("nan")
        if previous is not None:
            row.update(core.compare_refinement_runs(previous, row))
        attempts.append(dict(row))

        if row["refinement_stable"]:
            accepted = dict(row)
            accepted["quality_pass"] = True
            accepted["quality_status"] = "validated_refinement_stable"
            break

        # Recycle CatMAP's official numbers_map even if that stage's
        # full eight-step physical validation is not yet satisfactory.
        generated_cache = (output_root / surface / f"baseline_{surface}.pkl")
        if generated_cache.is_file() and row.get("status") == "solved":
            latest_seed = generated_cache.resolve()
        previous = row

    if accepted is not None:
        return accepted, attempts, bridge_logs

    last = dict(attempts[-1]) if attempts else {"surface": surface}
    last["quality_pass"] = False
    last["refinement_stable"] = False
    last["quality_status"] = "unresolved_after_continuation"
    return last, attempts, bridge_logs


def _boolean(series):
    return series.astype(str).str.lower().eq("true")


def recover_unvalidated_surfaces(
    surfaces: list[str] | None = None,
    bridge_temperatures: list[float] | None = None,
    steps: int = 17,
) -> None:
    """Update baseline only for independently validated recovered rows.

    Existing accepted results remain immutable. Failed attempts are recorded
    separately instead of replacing the existing baseline with failure rows.
    """
    base_file = core.OUTPUT / "dft_baseline.csv"
    if not base_file.is_file():
        raise FileNotFoundError(
            "Run the baseline analysis first, then recover unresolved surfaces."
        )
    baseline = pd.read_csv(base_file)
    baseline["surface"] = baseline["surface"].astype(str)
    existing_good = _boolean(baseline["quality_pass"]) & _boolean(baseline["refinement_stable"])

    if surfaces is None or len(surfaces) == 0:
        selected = baseline.loc[~existing_good, "surface"].tolist()
    else:
        selected = list(dict.fromkeys(surfaces))
        unknown = set(selected) - set(core.SURFACES)
        if unknown:
            raise ValueError(f"Unknown catalyst labels: {sorted(unknown)}")
        good_selected = set(baseline.loc[existing_good, "surface"])
        selected = [x for x in selected if x not in good_selected]

    if not selected:
        print("All selected surfaces already passed cross-precision verification.")
        return

    RECOVERY_ROOT.mkdir(parents=True, exist_ok=True)
    bridge_temperatures = bridge_temperatures or [550.0, 750.0]
    final_rows, stages, bridges = [], [], []
    newly_accepted = []

    for surface in selected:
        print(f"\n===== CatMAP temperature-continuation recovery: {surface} =====")
        try:
            result, attempts, bridge_log = recover_one_surface(
                surface, bridge_temperatures=bridge_temperatures, steps=steps
            )
        except Exception as err:
            result, attempts, bridge_log = (
                {"surface": surface, "status": f"FAILED: {type(err).__name__}: {err}",
                 "quality_status": "recovery_exception", "quality_pass": False,
                 "refinement_stable": False},
                [], []
            )
        final_rows.append(result)
        stages.extend(attempts)
        bridges.extend(bridge_log)
        print(f"[恢复结束] {surface}: {result['quality_status']}")

        if result.get("quality_pass") and result.get("refinement_stable"):
            newly_accepted.append(surface)
            for key, val in result.items():
                if key not in baseline.columns:
                    baseline[key] = np.nan
                baseline.loc[baseline["surface"].eq(surface), key] = val

    pd.DataFrame(bridges).to_csv(RECOVERY_ROOT / "bridge_diagnostics.csv", index=False)
    pd.DataFrame(stages).to_csv(RECOVERY_ROOT / "target_refinement.csv", index=False)
    pd.DataFrame(final_rows).to_csv(RECOVERY_ROOT / "recovery_results.csv", index=False)

    if newly_accepted:
        backup = RECOVERY_ROOT / "baseline_before_recovery.csv"
        if not backup.exists():
            import shutil
            shutil.copy2(base_file, backup)
        baseline.to_csv(base_file, index=False)
        descriptors = core.read_data()
        merged = descriptors.merge(baseline, on="surface", how="left")
        good = _boolean(merged["quality_pass"]) & _boolean(merged["refinement_stable"])
        log_tof = pd.to_numeric(merged["log10_net_C6H12O"], errors="coerce")
        merged["eligible_for_tof_fit"] = good & np.isfinite(log_tof)
        merged.to_csv(core.OUTPUT / "OOH_TOF_results.csv", index=False)
        selected_fit = merged.loc[merged["eligible_for_tof_fit"]]
        core.plot_linear(
            selected_fit, "G_OOH_eV", "log10_net_C6H12O",
            "OOH_logTOF_linear.png", "OOH* formation energy (eV)",
            "log10(net epoxide TOF)", min_points=5
        )

    print("\n新增稳定通过数:", len(newly_accepted),
          "; 催化剂:", ", ".join(newly_accepted) or "无")
    print("温度路径仅用于初始覆盖度；正式TOF仍为333.15 K下的DFT模型结果。")
    print("诊断文件:", RECOVERY_ROOT / "recovery_results.csv")


if __name__ == "__main__":
    recover_unvalidated_surfaces()
