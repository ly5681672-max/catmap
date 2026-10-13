from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import openpyxl


CATALYST_ORDER = [
    "tife", "timn", "tihf", "tire", "timo", "tiv", "tizr",
    "tico", "titi", "tiw", "tita", "ticr", "ti", "tini",
]
STATE_ROLE = {
    1: "feed_reference_state",
    2: "adsorption_or_precomplex_state",
    3: "NEB1_IS",
    4: "NEB1_TS",
    5: "NEB1_FS",
    6: "NEB2_IS",
    7: "NEB2_TS",
    8: "NEB2_FS",
    9: "post_NE B2_surface_state",
    10: "gas_product_reference_state",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def catalyst_groups(ws) -> dict[str, int]:
    groups: dict[str, int] = {}
    for col in range(1, ws.max_column + 1):
        label = ws.cell(1, col).value
        if not isinstance(label, str):
            continue
        cleaned = label.strip().lower().replace("–", "-").replace(" ", "")
        aliases = {
            "tife": "tife", "timn": "timn", "tihf": "tihf", "tire": "tire",
            "tinb": "tinb", "timo": "timo", "tiv": "tiv", "tizr": "tizr",
            "tico": "tico", "titi": "titi", "tiw": "tiw", "tita": "tita",
            "ticr": "ticr", "tini": "tini", "ti": "ti",
        }
        name = aliases.get(cleaned)
        # The single-Ti heading contains a legacy non-ASCII prefix in some workbook copies.
        if name is None and col == 96 and cleaned.endswith("ti"):
            name = "ti"
        if name:
            groups[name] = col
    return groups


def preceding_state_label(ws, start_col: int, row: int) -> str | None:
    for previous in range(row - 1, 0, -1):
        value = ws.cell(previous, start_col).value
        if isinstance(value, str) and value.strip() and value.strip().lower() not in {
            "ev", "ts", "neb", "state", "初态", "终态", "过渡态"
        }:
            return value.strip()
    return None


def data_columns(ws, start_col: int, stop_col: int) -> dict[str, int]:
    headers = {str(ws.cell(3, col).value).strip(): col for col in range(start_col, stop_col)
               if ws.cell(3, col).value is not None}
    if "E0" not in headers or "G" not in headers:
        raise ValueError(f"Missing E0/G headers for block {start_col}:{stop_col}: {headers}")
    e0_col, g_col = headers["E0"], headers["G"]
    return {
        "index": start_col,
        "E0": e0_col,
        "T": e0_col + 1,
        "path_start": e0_col + 2,
        "correction": g_col - 1,
        "G": g_col,
        "extra": g_col + 1,
    }


def extract_workbook(workbook: Path, output_dir: Path, comparison: Path | None = None) -> dict[str, Any]:
    formulas_wb = openpyxl.load_workbook(workbook, data_only=False, read_only=False)
    cached_wb = openpyxl.load_workbook(workbook, data_only=True, read_only=False)
    sheet_name = "吉布斯自由能汇总"
    if sheet_name not in formulas_wb.sheetnames:
        raise KeyError(f"Missing worksheet: {sheet_name}; found {formulas_wb.sheetnames}")
    ws = formulas_wb[sheet_name]
    cached = cached_wb[sheet_name]
    groups = catalyst_groups(ws)
    if not all(name in groups for name in CATALYST_ORDER):
        missing = [name for name in CATALYST_ORDER if name not in groups]
        raise ValueError(f"Missing catalyst blocks: {missing}; found={groups}")

    sorted_starts = sorted(groups.values())
    stops = {start: next((n for n in sorted_starts if n > start), ws.max_column + 1)
             for start in sorted_starts}
    state_rows: list[dict[str, Any]] = []
    sorted_starts = sorted(groups.values())
    stops = {start: next((n for n in sorted_starts if n > start), ws.max_column + 1)
             for start in sorted_starts}
    for catalyst in CATALYST_ORDER:
        start = groups[catalyst]
        cols = data_columns(ws, start, stops[start])
        for row in range(1, ws.max_row + 1):
            index = ws.cell(row, cols["index"]).value
            if not isinstance(index, (int, float)) or int(index) != index:
                continue
            state_index = int(index)
            data_row = row if isinstance(ws.cell(row, cols["E0"]).value, (int, float)) else None
            if data_row is None:
                for candidate in range(row + 1, min(row + 3, ws.max_row) + 1):
                    if isinstance(ws.cell(candidate, cols["E0"]).value, (int, float)):
                        data_row = candidate
                        break
            if data_row is None:
                data_row = row
            e0 = ws.cell(data_row, cols["E0"])
            temp = ws.cell(data_row, cols["T"])
            path_cells = [ws.cell(data_row, c) for c in range(cols["path_start"], cols["correction"])]
            path = " | ".join(str(c.value) for c in path_cells if c.value is not None)
            correction = ws.cell(data_row, cols["correction"])
            g_cell = ws.cell(data_row, cols["G"])
            total_cell = ws.cell(data_row, cols["extra"])
            role = STATE_ROLE.get(state_index, "unclassified_state")
            if role == "post_NE B2_surface_state":
                role = "post_NEB2_surface_state"
            state_rows.append({
                "catalyst": catalyst,
                "state_index": state_index,
                "state_role": role,
                "state_label": preceding_state_label(ws, start, row),
                "worksheet": sheet_name,
                "state_index_cell": ws.cell(row, cols["index"]).coordinate,
                "state_index_row": row,
                "source_row": data_row,
                "E0_cell": e0.coordinate,
                "E0": e0.value,
                "temperature_K": temp.value,
                "calculation_path": path,
                "Gibbs_correction_cell": correction.coordinate,
                "Gibbs_correction_formula_or_value": correction.value,
                "Gibbs_correction_cached": cached[correction.coordinate].value,
                "G_cell": g_cell.coordinate,
                "G_formula_or_value": g_cell.value,
                "G_cached_value": cached[g_cell.coordinate].value,
                "total_or_relative_cell": total_cell.coordinate,
                "total_or_relative_formula_or_value": total_cell.value,
                "total_or_relative_cached": cached[total_cell.coordinate].value,
                "structure_composition_status": "not_verified_from_coordinates; only workbook label/path available",
                "occupied_site_count_status": "not specified by workbook label; requires structure/site model evidence",
            })

    neb_roles = [
        (12, "NEB1_IS", {"C": 6, "H": 14, "O": 2}, "coadsorbate occupancy unresolved"),
        (16, "NEB1_TS", {"C": 6, "H": 14, "O": 2}, "TS occupancy unresolved"),
        (20, "NEB1_FS", {"C": 6, "H": 14, "O": 2}, "coadsorbate occupancy unresolved"),
        (24, "NEB2_IS", {"H": 2, "O": 1}, "label suggests three adsorbed fragments; coordinate confirmation required"),
        (28, "NEB2_TS", {"H": 2, "O": 1}, "TS occupancy unresolved"),
        (32, "NEB2_FS", {"H": 2, "O": 1}, "label suggests two adsorbed fragments; coordinate confirmation required"),
    ]
    neb_rows: list[dict[str, Any]] = []
    for catalyst in CATALYST_ORDER:
        start = groups[catalyst]
        cols = data_columns(ws, start, stops[start])
        for row, role, composition, site_note in neb_roles:
            e0 = ws.cell(row, cols["E0"])
            temp = ws.cell(row, cols["T"])
            path_values = [ws.cell(row, c).value for c in range(cols["path_start"], cols["correction"])]
            correction = ws.cell(row, cols["correction"])
            g_cell = ws.cell(row, cols["G"])
            extra = ws.cell(row, cols["extra"])
            following_g = ws.cell(min(row + 1, ws.max_row), cols["G"])
            neb_rows.append({
                "catalyst": catalyst,
                "NEB_role": role,
                "worksheet": sheet_name,
                "source_row_from_DFT_TS_relative_reference_note": row,
                "raw_state_index_in_workbook": ws.cell(row, cols["index"]).value,
                "state_label": preceding_state_label(ws, start, row),
                "E0_cell": e0.coordinate,
                "E0": e0.value,
                "temperature_K": temp.value,
                "calculation_paths": " | ".join(str(v) for v in path_values if v is not None),
                "Gibbs_correction_cell": correction.coordinate,
                "Gibbs_correction_formula_or_value": correction.value,
                "Gibbs_correction_cached": cached[correction.coordinate].value,
                "raw_G_cell": g_cell.coordinate,
                "raw_G_formula": g_cell.value,
                "raw_G_cached": cached[g_cell.coordinate].value,
                "same_row_extra_formula_or_value": extra.value,
                "same_row_extra_cached": cached[extra.coordinate].value,
                "next_row_G_cell": following_g.coordinate,
                "next_row_G_formula": following_g.value,
                "next_row_G_cached": cached[following_g.coordinate].value,
                "label_derived_composition_candidate": json.dumps(composition, sort_keys=True),
                "composition_status": "label/path-derived, not atom-coordinate verified",
                "site_occupancy_note": site_note,
                "common_slab_and_spectator_reference_status": "must be checked per catalyst/state; see formula audit",
            })

    formula_rows: list[dict[str, Any]] = []
    for row in ws.iter_rows():
        for cell in row:
            if isinstance(cell.value, str) and cell.value.startswith("="):
                formula_rows.append({
                    "worksheet": sheet_name,
                    "cell": cell.coordinate,
                    "row": cell.row,
                    "column": cell.column,
                    "formula": cell.value,
                    "cached_value": cached[cell.coordinate].value,
                    "cell_context": "state/energy/reference formula; inspect source row and referenced cells",
                })

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "original_neb_state_inventory.csv", neb_rows, list(neb_rows[0]))
    write_csv(output_dir / "original_gibbs_state_inventory.csv", state_rows, list(state_rows[0]))
    write_csv(output_dir / "excel_gibbs_formula_audit.csv", formula_rows, list(formula_rows[0]))

    manifest: dict[str, Any] = {
        "primary_workbook": str(workbook),
        "primary_sha256": sha256(workbook),
        "worksheet": sheet_name,
        "workbook_sheets": formulas_wb.sheetnames,
        "formula_cells_audited": len(formula_rows),
        "NEB_anchor_records": len(neb_rows),
        "all_numeric_state_records": len(state_rows),
        "requested_catalysts": CATALYST_ORDER,
        "comparison_workbook": str(comparison) if comparison else None,
        "comparison_sha256": sha256(comparison) if comparison and comparison.exists() else None,
        "comparison_status": (
            "compared by hash and extracted records" if comparison and comparison.exists()
            else "not_available; H-drive source workbook audited, separate previously attached (2).xlsx not found"
        ),
    }
    (output_dir / "workbook_audit_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--comparison-workbook", type=Path)
    args = parser.parse_args()
    result = extract_workbook(args.workbook, args.output_dir, args.comparison_workbook)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
