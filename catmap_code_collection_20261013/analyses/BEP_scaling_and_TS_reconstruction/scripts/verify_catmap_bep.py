"""Static BEP-mode check against the installed CatMAP 0.3.1 scaler.

This script parses and verifies an isolated copy of the existing descriptor
grid setup, asks the official GeneralizedLinearScaler for its BEP coefficients,
and never calls ReactionModel.run().
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import catmap
from catmap import ReactionModel


PACKAGE = Path(__file__).resolve().parents[1]
DATA = PACKAGE / "data"
REPORTS = PACKAGE / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)
SOURCE_SETUP = DATA / "source_grid_setup.mkm"
SOURCE_ENERGIES = DATA / "source_strict_energies.txt"
PROBE_SETUP = REPORTS / "catmap_bep_mode_probe.mkm"
NO_RUN_CACHE = REPORTS / "catmap_bep_mode_probe_not_run.pkl"


def main() -> None:
    if catmap.__version__ != "0.3.1":
        raise RuntimeError(f"Expected CatMAP 0.3.1, found {catmap.__version__}")

    text = SOURCE_SETUP.read_text(encoding="utf-8")
    old_mode = "transition_state_scaling_mode = 'initial_state'"
    if text.count(old_mode) != 1:
        raise RuntimeError("Could not uniquely locate initial_state mode in source setup")
    text = text.replace(old_mode, "transition_state_scaling_mode = 'BEP'", 1)
    old_input = "input_file = 'energies.txt'"
    if text.count(old_input) != 1:
        raise RuntimeError("Could not uniquely locate input_file in source setup")
    text = text.replace(old_input, f"input_file = {str(SOURCE_ENERGIES)!r}", 1)
    old_cache = "data_file = 'catmap_ts1_ts2_5x5.pkl'"
    if text.count(old_cache) != 1:
        raise RuntimeError("Could not uniquely locate data_file in source setup")
    text = text.replace(old_cache, f"data_file = {str(NO_RUN_CACHE)!r}", 1)
    PROBE_SETUP.write_text(text, encoding="utf-8")

    model = ReactionModel(setup_file=str(PROBE_SETUP))
    if len(model.surface_names) != 13:
        raise AssertionError(f"Expected 13 valid surfaces, got {len(model.surface_names)}")
    if len(model.elementary_rxns) != 6:
        raise AssertionError(f"Expected 6 reactions, got {len(model.elementary_rxns)}")

    scaler = model.scaler
    # Follow the scaler's normal static call path: get_coefficient_matrix()
    # initializes default adsorbate constraints before building TS rows.
    scaler.get_coefficient_matrix()
    coeffs = scaler.transition_state_scaling_coefficients
    output = []
    for name, pair in zip(model.transition_state_names, coeffs):
        state = name.removesuffix("_s")
        if state == "TS1":
            path = "NEB1_path2"
            initial, final = "Rpair_s", "EpoxideHHO_s"
        elif state == "TS2":
            path = "NEB2_path4"
            initial, final = "HHO_s", "OHHa_s"
        else:
            continue
        output.append({
            "catmap_version": catmap.__version__,
            "path": path,
            "transition_state": name,
            "initial_state": initial,
            "final_state": final,
            "transition_state_scaling_mode": scaler.transition_state_scaling_mode,
            "official_alpha": float(pair[0]),
            "official_beta_eV": float(pair[1]),
            "coefficient_source": "CatMAP GeneralizedLinearScaler.get_transition_state_scaling_matrix",
            "reaction_model_run_called": False,
        })
    if {r["path"] for r in output} != {"NEB1_path2", "NEB2_path4"}:
        raise AssertionError(f"Expected TS1 and TS2 BEP entries, got {output}")

    csv_path = DATA / "official_BEP_scaler_coefficients.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    metadata = {
        "catmap_version": catmap.__version__,
        "source_setup": str(SOURCE_SETUP),
        "probe_setup": str(PROBE_SETUP),
        "strict_energy_snapshot": str(SOURCE_ENERGIES),
        "surface_count": len(model.surface_names),
        "reaction_count": len(model.elementary_rxns),
        "transition_state_names": list(model.transition_state_names),
        "scaler": type(scaler).__name__,
        "mode": scaler.transition_state_scaling_mode,
        "reaction_model_run_called": False,
        "probe_data_cache_created": NO_RUN_CACHE.exists(),
        "coefficients": output,
    }
    (DATA / "official_BEP_scaler_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
