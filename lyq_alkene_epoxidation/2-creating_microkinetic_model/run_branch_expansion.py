# -*- coding: utf-8 -*-
"""Run additional 333.15 K starts in an isolated branch-audit directory."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import recover_ooh_solver as recovery
import run_ooh_analysis as core


AUDIT_ROOT = core.OUTPUT / "steady_state_branch_audit"
DEFAULT_SURFACES = ("timn", "tihf", "tiv", "ticr", "tizr")
HISTORICAL_SURFACES = ("tife", "tico", "tini", "tiw", "ti")


def run_surfaces(surfaces: list[str], tag: str) -> Path:
    output_root = AUDIT_ROOT / "expanded_multistart" / tag
    output_root.mkdir(parents=True, exist_ok=True)
    results = []
    trials = []
    seeds = []
    failures = []
    for surface in surfaces:
        print(f"[expanded-branch-audit] {surface} | tag={tag}")
        try:
            final, surface_trials, surface_seeds = recovery.recover_one_surface(
                surface, output_root=output_root, expanded=True
            )
        except Exception as exc:
            final = {
                "surface": surface,
                "status": f"FAILED: {type(exc).__name__}: {exc}",
                "quality_status": "recovery_exception",
                "quality_pass": False,
                "refinement_stable": False,
            }
            surface_trials, surface_seeds = [], []
            failures.append({
                "surface": surface,
                "stage": "expanded_multistart",
                "error": f"{type(exc).__name__}: {exc}",
            })
        results.append(final)
        trials.extend(surface_trials)
        seeds.extend(surface_seeds)
        print(
            f"[{surface}] {final.get('quality_status')}; "
            f"validated={final.get('initial_seed_validated_count', 0)}; "
            f"matching_pairs={final.get('independent_seed_matching_pair_count', 0)}"
        )

    pd.DataFrame(results).to_csv(output_root / "expanded_recovery_results.csv", index=False)
    pd.DataFrame(trials).to_csv(output_root / "expanded_seed_trials.csv", index=False)
    pd.DataFrame(seeds).to_csv(output_root / "expanded_seed_definitions.csv", index=False)
    pd.DataFrame(failures).to_csv(output_root / "expanded_failures.csv", index=False)
    print(f"[expanded-branch-audit] output={output_root}")
    return output_root


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", required=True)
    parser.add_argument("--surfaces", nargs="*", default=None)
    parser.add_argument("--historical", action="store_true")
    args = parser.parse_args()
    surfaces = args.surfaces
    if surfaces is None:
        surfaces = list(HISTORICAL_SURFACES if args.historical else DEFAULT_SURFACES)
    unknown = set(surfaces) - set(core.SURFACES)
    if unknown:
        parser.error(f"unknown surfaces: {sorted(unknown)}")
    run_surfaces(list(dict.fromkeys(surfaces)), args.tag)


if __name__ == "__main__":
    main()
