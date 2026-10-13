"""Build traceable 14-surface serial-NEB tables without inventing global energies."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from serial_species import convert_state_energy
from validate_serial_energies import (
    DATA,
    ENERGY_MANIFEST,
    MAPPING_FILE,
    epoxide_connection,
    replay_neb_block,
)


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = DATA
H_MODEL_DIR = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model")
GIBBS_STATE_FILE = DATA / "original_gibbs_state_inventory.csv"
OLD_ENERGY_FILE = H_MODEL_DIR / "energies_dft_neb.txt"
CATALYSTS = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "tico",
             "titi", "tiw", "tita", "ticr", "ti", "tini"]
SPECIES_PER_PATH = {
    "NEB1": ("Rpair_s", "TS1_s", "EpoxideHHO_s"),
    "NEB2": ("HHO_s", "TS2_s", "OHHa_s"),
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _states() -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = {name: [] for name in CATALYSTS}
    for row in _read_csv(GIBBS_STATE_FILE):
        if row["catalyst"] in result:
            result[row["catalyst"]].append(row)
    return result


def _post_state(rows: list[dict[str, str]], state_index: str, minimum_source_row: int) -> dict[str, str] | None:
    matches = [row for row in rows if row["state_index"] == state_index
               and int(row["source_row"]) >= minimum_source_row]
    if not matches:
        return None
    # When worksheet layout reuses an index on a NEB anchor and the following
    # state row, choose the post-NEB row (largest source row) and preserve its label.
    return max(matches, key=lambda row: int(row["source_row"]))


def _mapping() -> dict[tuple[str, str], dict[str, str]]:
    return {(row["catalyst"], row["path"]): row for row in _read_csv(MAPPING_FILE)}


def build_species_energy_rows() -> list[dict[str, Any]]:
    map_rows = _mapping()
    all_states = _states()
    rows: list[dict[str, Any]] = []
    for catalyst in CATALYSTS:
        for path, species_names in SPECIES_PER_PATH.items():
            record = map_rows[(catalyst, path)]
            raw_energies = [float(record[key]) for key in (
                "G_IS_workbook_or_normalized", "G_TS_workbook_or_normalized", "G_FS_workbook_or_normalized"
            )]
            local_values = [value - raw_energies[0] for value in raw_energies]
            source_cells = [record["IS_source_cell"], record["TS_source_cell"], record["FS_source_cell"]]
            source_formulas = [record["IS_source_expression"], record["TS_source_expression"],
                               record["FS_source_expression_or_adjustment"]]
            for role, (species, raw_g, local_g, cell, formula) in enumerate(
                zip(species_names, raw_energies, local_values, source_cells, source_formulas)
            ):
                role_name = ("IS", "TS", "FS")[role]
                rows.append({
                    "catalyst": catalyst,
                    "species_key": species,
                    "physical_state": f"{path}_{role_name}",
                    "source_cell": cell,
                    "source_formula_or_adjustment": formula,
                    "source_G_eV": raw_g,
                    "block_local_relative_G_eV": local_g,
                    "local_reference_state": f"{path}_IS",
                    "formation_energy_eV": "",
                    "formation_energy_status": "UNRESOLVED_REFERENCE",
                    "reason": "clean-slab/common formation-energy reference and cross-path offset are unavailable",
                })
        states = all_states[catalyst]
        water = _post_state(states, "9", 36)
        if water:
            water_g = float(water["G_cached_value"])
            neb2 = map_rows[(catalyst, "NEB2")]
            hho_g = float(neb2["G_IS_workbook_or_normalized"])
            conversion = convert_state_energy(
                raw_g_eV=water_g,
                source_cell=water["G_cell"],
                source_formula=water["G_formula_or_value"],
                common_reference_eV=None,
                spectator_adjustments_eV=[],
                spectator_notes=["post-NEB water state label inspected; exact spectator/slab relationship remains unverified"],
            )
            rows.append({
                "catalyst": catalyst,
                "species_key": "H2O_s",
                "physical_state": "post_NEB2_water_adsorbed_candidate",
                "source_cell": water["G_cell"],
                "source_formula_or_adjustment": water["G_formula_or_value"],
                "source_G_eV": water_g,
                "block_local_relative_G_eV": water_g - hho_g,
                "local_reference_state": "NEB2_IS (conditional; state-to-state slab reference unverified)",
                "formation_energy_eV": "",
                "formation_energy_status": conversion["conversion_status"],
                "reason": conversion["reason"],
                "workbook_state_label": water["state_label"],
            })
        else:
            rows.append({
                "catalyst": catalyst,
                "species_key": "H2O_s",
                "physical_state": "post_NEB2_water_adsorbed_candidate",
                "source_cell": "",
                "source_formula_or_adjustment": "",
                "source_G_eV": "",
                "block_local_relative_G_eV": "",
                "local_reference_state": "",
                "formation_energy_eV": "",
                "formation_energy_status": "MISSING_SOURCE_STATE",
                "reason": "no matching post-NEB2 water state at the documented source row",
            })
    return rows


def build_barrier_rows() -> list[dict[str, Any]]:
    return [replay_neb_block(row, row["path"]) for row in _read_csv(MAPPING_FILE)]


def build_surface_rows() -> list[dict[str, Any]]:
    map_rows = _mapping()
    states = _states()
    with ENERGY_MANIFEST.open("r", encoding="utf-8") as stream:
        energy_manifest = json.load(stream)
    water_gas = float(energy_manifest["gas_reference_cached_values"]["H2O_g"])
    rows = []
    for catalyst in CATALYSTS:
        neb1 = map_rows[(catalyst, "NEB1")]
        neb2 = map_rows[(catalyst, "NEB2")]
        connection = epoxide_connection(catalyst)
        state9 = _post_state(states[catalyst], "9", 36)
        state10 = _post_state(states[catalyst], "10", 41)
        if state9:
            water_g = float(state9["G_cached_value"])
            water_formation_dg = water_g - float(neb2["G_FS_workbook_or_normalized"])
            water_status = "CONDITIONAL_SLAB_AND_SPECTATOR_REFERENCE_UNVERIFIED"
            water_cell, water_formula = state9["G_cell"], state9["G_formula_or_value"]
            water_label = state9["state_label"]
        else:
            water_g = None
            water_formation_dg = None
            water_status = "MISSING_SOURCE_STATE"
            water_cell = water_formula = water_label = ""
        if state9 and state10:
            # Workbook state 10 represents the product reference; this expression
            # follows the documented H2O gas reference once. It remains conditional
            # because the clean-slab/spectator basis cannot be matched from labels alone.
            water_desorption_dg = float(state10["G_cached_value"]) + water_gas - water_g
            desorption_status = "CONDITIONAL_SLAB_AND_SPECTATOR_REFERENCE_UNVERIFIED"
            state10_cell, state10_formula = state10["G_cell"], state10["G_formula_or_value"]
        else:
            water_desorption_dg = None
            desorption_status = "MISSING_SOURCE_STATE"
            state10_cell = state10_formula = ""
        rows.append({
            "catalyst": catalyst,
            "NEB1_status": neb1["source_energy_crosswalk_status"],
            "NEB1_IS_G_eV": neb1["G_IS_workbook_or_normalized"],
            "NEB1_TS_G_eV": neb1["G_TS_workbook_or_normalized"],
            "NEB1_FS_G_eV": neb1["G_FS_workbook_or_normalized"],
            "NEB1_raw_TS_minus_IS_eV": neb1["raw_TS_minus_IS"],
            "NEB1_raw_TS_minus_FS_eV": neb1["raw_TS_minus_FS"],
            "NEB1_CatMAP_Ea_forward_eV": neb1["CatMAP_effective_Ea_forward"],
            "NEB1_CatMAP_Ea_reverse_eV": neb1["CatMAP_effective_Ea_reverse"],
            "NEB1_FS_source_formula": neb1["FS_source_expression_or_adjustment"],
            "NEB2_status": neb2["source_energy_crosswalk_status"],
            "NEB2_IS_G_eV": neb2["G_IS_workbook_or_normalized"],
            "NEB2_TS_G_eV": neb2["G_TS_workbook_or_normalized"],
            "NEB2_FS_G_eV": neb2["G_FS_workbook_or_normalized"],
            "NEB2_raw_TS_minus_IS_eV": neb2["raw_TS_minus_IS"],
            "NEB2_raw_TS_minus_FS_eV": neb2["raw_TS_minus_FS"],
            "NEB2_CatMAP_Ea_forward_eV": neb2["CatMAP_effective_Ea_forward"],
            "NEB2_CatMAP_Ea_reverse_eV": neb2["CatMAP_effective_Ea_reverse"],
            "NEB2_FS_reference_note": neb2["FS_source_expression_or_adjustment"],
            "epoxide_desorption_delta_G_eV": connection["delta_G_eV"],
            "epoxide_desorption_activation_barrier_eV": "",
            "epoxide_desorption_status": connection["physical_status"],
            "H2O_state_source_cell": water_cell,
            "H2O_state_formula": water_formula,
            "H2O_state_label": water_label,
            "H2O_state_status": water_status,
            "water_formation_delta_G_eV": "" if water_formation_dg is None else water_formation_dg,
            "water_formation_status": water_status,
            "water_product_state_source_cell": state10_cell,
            "water_product_state_formula": state10_formula,
            "water_desorption_delta_G_eV": "" if water_desorption_dg is None else water_desorption_dg,
            "water_desorption_status": desorption_status,
            "common_slab_reference_status": "UNRESOLVED",
            "formation_energy_status": "UNRESOLVED_REFERENCE",
            "kinetics_ready": False,
        })
    return rows


def _gas_formation_rows() -> list[dict[str, str]]:
    result = []
    with OLD_ENERGY_FILE.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            if row["surface_name"] == "None" and row["species_name"] in {"H2O2", "C6H12", "C6H12O", "H2O"}:
                result.append(row)
    if len(result) != 4:
        raise ValueError("Expected four unchanged gas-reference rows in the original energy table")
    return result


def write_tire_energy_file(path: Path) -> None:
    rows = []
    for gas in _gas_formation_rows():
        rows.append({
            "surface_name": "None", "site_name": "gas", "species_name": gas["species_name"],
            "formation_energy": gas["formation_energy"], "frequencies": "[]",
            "reference": f"inherited unchanged from original gas reference; {gas['reference']}",
        })
    tire = {row["path"]: row for row in _read_csv(MAPPING_FILE) if row["catalyst"] == "tire"}
    for path_name, names in SPECIES_PER_PATH.items():
        audited = replay_neb_block(tire[path_name], path_name)
        energy_values = audited["block_local_G_eV"]
        local_values = [value - energy_values[0] for value in energy_values]
        for name, local_energy in zip(names, local_values):
            species_name = name.removesuffix("_s")
            rows.append({
                "surface_name": "tire", "site_name": "111", "species_name": species_name,
                "formation_energy": f"{local_energy:.12f}", "frequencies": "[]",
                "reference": f"{path_name} IS local zero; diagnostic only; NOT global formation energy; do not load for coupled rates",
            })
    states = _states()["tire"]
    water = _post_state(states, "9", 36)
    if water:
        neb2_is = float(tire["NEB2"]["G_IS_workbook_or_normalized"])
        rows.append({
            "surface_name": "tire", "site_name": "111", "species_name": "H2O",
            "formation_energy": f"{float(water['G_cached_value']) - neb2_is:.12f}", "frequencies": "[]",
            "reference": "NEB2 IS local zero assumed; state-reference equivalence unverified; diagnostic only",
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["surface_name", "site_name", "species_name",
                                                     "formation_energy", "frequencies", "reference"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs() -> dict[str, Any]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    surface_rows = build_surface_rows()
    species_rows = build_species_energy_rows()
    barrier_rows = build_barrier_rows()
    _write_csv(OUTPUT_DIR / "SERIAL_14_SURFACE_MAPPING.csv", surface_rows)
    _write_csv(OUTPUT_DIR / "SERIAL_SPECIES_ENERGIES.csv", species_rows)
    _write_csv(OUTPUT_DIR / "SERIAL_NEB_BARRIER_AUDIT.csv", barrier_rows)
    write_tire_energy_file(ROOT / "models" / "v1_serial" / "energies_tire_serial_neb.txt")
    return {"surface_count": len(surface_rows), "neb_barrier_rows": len(barrier_rows),
            "species_energy_rows": len(species_rows), "tof_calculated": False}


if __name__ == "__main__":
    print(json.dumps(write_outputs(), ensure_ascii=False, indent=2))
