"""Independent NEB energy replay with CatMAP's 0.3.1 barrier semantics."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from catmap.model import ReactionModel
from catmap.parsers.table_parser import TableParser

from validate_serial_network import build_model


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
MAPPING_FILE = DATA / "neb_to_catmap_energy_mapping.csv"
ENERGY_MANIFEST = DATA / "energy_mapping_manifest.json"
HERE = Path(__file__).resolve().parent
UNIQUE_FORMATION_REFERENCE_STATUS = "UNRESOLVED_COMMON_SLAB_REFERENCE"
LOCAL_REFERENCE_SCOPE = "BLOCK_LOCAL_NOT_GLOBAL_FORMATION_ENERGY"


def _mapping_rows() -> list[dict[str, str]]:
    with MAPPING_FILE.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def source_row(catalyst: str, path: str) -> dict[str, str]:
    matches = [row for row in _mapping_rows()
               if row["catalyst"] == catalyst and row["path"] == path]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one mapping for {catalyst}/{path}, found {len(matches)}")
    return matches[0]


def _energy_dict(row: dict[str, str], path: str, common_shift_eV: float) -> tuple[dict[str, float], list[str]]:
    g_is = float(row["G_IS_workbook_or_normalized"]) + common_shift_eV
    g_ts = float(row["G_TS_workbook_or_normalized"]) + common_shift_eV
    g_fs = float(row["G_FS_workbook_or_normalized"]) + common_shift_eV
    names = ("Rpair_s", "TS1_s", "EpoxideHHO_s") if path == "NEB1" else ("HHO_s", "TS2_s", "OHHa_s")
    return dict(zip(names, (g_is, g_ts, g_fs))), list(names)


def replay_neb_block(
    row: dict[str, str], path: str, common_shift_eV: float = 0.0
) -> dict[str, Any]:
    energy_dict, names = _energy_dict(row, path, common_shift_eV)
    model = ReactionModel()
    forward = model.get_rxn_energy([[names[0]], [names[1]], [names[2]]], energy_dict)
    reverse = model.get_rxn_energy([[names[2]], [names[1]], [names[0]]], energy_dict)
    g_is, g_ts, g_fs = (energy_dict[name] for name in names)
    raw_fwd = g_ts - g_is
    raw_rev = g_ts - g_fs
    dg = g_fs - g_is
    return {
        "catalyst": row["catalyst"],
        "path": path,
        "source_cells": [row["IS_source_cell"], row["TS_source_cell"], row["FS_source_cell"]],
        "source_formulas": [row["IS_source_expression"], row["TS_source_expression"],
                            row["FS_source_expression_or_adjustment"]],
        "source_G_eV": [float(row["G_IS_workbook_or_normalized"]),
                         float(row["G_TS_workbook_or_normalized"]),
                         float(row["G_FS_workbook_or_normalized"])],
        "block_local_G_eV": [g_is, g_ts, g_fs],
        "raw_TS_minus_IS_eV": raw_fwd,
        "raw_TS_minus_FS_eV": raw_rev,
        "raw_FS_minus_IS_eV": dg,
        "CatMAP_effective_Ea_forward_eV": forward[1],
        "CatMAP_effective_Ea_reverse_eV": reverse[1],
        "energy_reference_scope": LOCAL_REFERENCE_SCOPE,
        "formation_energy_status": UNIQUE_FORMATION_REFERENCE_STATUS,
        "valid_for_coupled_network_rate_calculation": False,
    }


def epoxide_connection(catalyst: str) -> dict[str, Any]:
    mapping = {row["path"]: row for row in _mapping_rows() if row["catalyst"] == catalyst}
    if set(mapping) != {"NEB1", "NEB2"}:
        raise ValueError(f"Incomplete NEB pair for {catalyst}")
    with ENERGY_MANIFEST.open("r", encoding="utf-8") as stream:
        manifest = json.load(stream)
    epoxide_g = float(manifest["gas_reference_cached_values"]["C6H12O_g"])
    g_neb1_fs = float(mapping["NEB1"]["G_FS_workbook_or_normalized"])
    g_neb2_is = float(mapping["NEB2"]["G_IS_workbook_or_normalized"])
    # One epoxide gas reference is added to the post-desorption state exactly once.
    g_after = g_neb2_is + epoxide_g
    delta = g_after - g_neb1_fs
    return {
        "catalyst": catalyst,
        "NEB1_FS_source_cell": mapping["NEB1"]["FS_source_cell"],
        "NEB1_FS_source_formula": mapping["NEB1"]["FS_source_expression_or_adjustment"],
        "NEB1_FS_G_eV": g_neb1_fs,
        "NEB2_IS_source_cell": mapping["NEB2"]["IS_source_cell"],
        "NEB2_IS_source_formula": mapping["NEB2"]["IS_source_expression"],
        "NEB2_IS_G_eV": g_neb2_is,
        "C6H12O_g_reference_eV": epoxide_g,
        "epoxide_gas_reference_count": 1,
        "post_desorption_endpoint_G_eV": g_after,
        "delta_G_eV": delta,
        "activation_barrier_eV": None,
        "physical_status": "CONDITIONAL_COMMON_SLAB_UNVERIFIED",
        "reason": "same slab/supercell, raw atom mapping, and relaxation equivalence are not established",
    }


def unactivated_step_barrier(delta_g_eV: float) -> dict[str, Any]:
    model = ReactionModel()
    model.species_definitions = {
        "IS_s": {"type": "adsorbate"},
        "FS_s": {"type": "adsorbate"},
    }
    forward = model.get_rxn_energy([["IS_s"], ["FS_s"]], {"IS_s": 0.0, "FS_s": delta_g_eV})
    reverse = model.get_rxn_energy([["FS_s"], ["IS_s"]], {"IS_s": 0.0, "FS_s": delta_g_eV})
    return {
        "delta_G_eV": delta_g_eV,
        "forward_effective_Ea_eV": forward[1],
        "reverse_effective_Ea_eV": reverse[1],
        "barrier_status": "SOFTWARE_ASSUMPTION_NOT_DFT",
    }


def mapping_barrier_records() -> list[dict[str, Any]]:
    results = []
    for row in _mapping_rows():
        path = row["path"]
        values = replay_neb_block(row, path)
        values["catalyst"] = row["catalyst"]
        results.append(values)
    return results


def parse_diagnostic_table() -> dict[str, Any]:
    """Test the local CatMAP table parser without solving the diagnostic model."""
    model = build_model()
    parser = TableParser(model)
    parser.input_file = str(HERE / "energies_tire_serial_neb.txt")
    parser.parse()
    return {
        "catmap_version": __import__("catmap").__version__,
        "parser_status": "PASS_DIAGNOSTIC_TABLE_ONLY",
        "coupled_rate_valid": False,
        "reason": "NEB1 and NEB2 adsorbate entries use distinct block-local references",
    }


if __name__ == "__main__":
    for path in ("NEB1", "NEB2"):
        print(json.dumps(replay_neb_block(source_row("tire", path), path), ensure_ascii=False, indent=2))
    print(json.dumps(epoxide_connection("tire"), ensure_ascii=False, indent=2))
