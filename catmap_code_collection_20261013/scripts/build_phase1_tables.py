from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import openpyxl
from openpyxl import Workbook
from catmap.model import ReactionModel

from audit_catmap_031 import barrier_values
from extract_neb_gibbs import CATALYST_ORDER, catalyst_groups, data_columns, sha256, write_csv


STATE_IDS = {"NEB1_IS": 3, "NEB1_TS": 4, "NEB1_FS": 5,
             "NEB2_IS": 6, "NEB2_TS": 7, "NEB2_FS": 8}
GAS_REFERENCE_CELLS = {"H2O2_g": "F57", "C6H12_g": "F58", "H2O_g": "F59", "C6H12O_g": "F60"}


def read_cached_state_data(workbook: Path) -> tuple[dict[str, dict[int, dict[str, Any]]], dict[str, float]]:
    formulas = openpyxl.load_workbook(workbook, data_only=False, read_only=False)["吉布斯自由能汇总"]
    cached = openpyxl.load_workbook(workbook, data_only=True, read_only=False)["吉布斯自由能汇总"]
    groups = catalyst_groups(formulas)
    sorted_starts = sorted(groups.values())
    stops = {start: next((n for n in sorted_starts if n > start), formulas.max_column + 1)
             for start in sorted_starts}
    gas_refs = {name: cached[cell].value for name, cell in GAS_REFERENCE_CELLS.items()}
    states: dict[str, dict[int, dict[str, Any]]] = {}
    for catalyst in CATALYST_ORDER:
        start = groups[catalyst]
        cols = data_columns(formulas, start, stops[start])
        states[catalyst] = {}
        for row in range(1, formulas.max_row + 1):
            idx = formulas.cell(row, cols["index"]).value
            if not isinstance(idx, (int, float)) or int(idx) != idx:
                continue
            idx = int(idx)
            data_row = row if isinstance(formulas.cell(row, cols["E0"]).value, (int, float)) else None
            if data_row is None:
                for candidate in range(row + 1, min(row + 3, formulas.max_row) + 1):
                    if isinstance(formulas.cell(candidate, cols["E0"]).value, (int, float)):
                        data_row = candidate
                        break
            if data_row is None:
                data_row = row
            g_cell = formulas.cell(data_row, cols["G"])
            total_cell = formulas.cell(data_row, cols["extra"])
            states[catalyst][idx] = {
                "index_row": row,
                "row": data_row,
                "G_cell": g_cell.coordinate,
                "G_formula": g_cell.value,
                "G": cached[g_cell.coordinate].value,
                "extra_cell": total_cell.coordinate,
                "extra_formula": total_cell.value,
                "extra_cached": cached[total_cell.coordinate].value,
                "path": " | ".join(str(formulas.cell(data_row, c).value) for c in range(cols["path_start"], cols["correction"])
                                     if formulas.cell(data_row, c).value is not None),
                "correction": formulas.cell(data_row, cols["correction"]).value,
            }
    return states, gas_refs


def load_relative(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return {row["surface"].strip().lower(): row for row in csv.DictReader(stream)}


def mapping_rows(states, relative, workbook: Path) -> list[dict[str, Any]]:
    rows = []
    model = ReactionModel()
    formulas_ws = openpyxl.load_workbook(workbook, data_only=False, read_only=False)["吉布斯自由能汇总"]
    cached_ws = openpyxl.load_workbook(workbook, data_only=True, read_only=False)["吉布斯自由能汇总"]
    groups = catalyst_groups(formulas_ws)
    starts = sorted(groups.values())
    stops = {start: next((n for n in starts if n > start), formulas_ws.max_column + 1)
             for start in starts}
    layouts = {name: data_columns(formulas_ws, groups[name], stops[groups[name]]) for name in CATALYST_ORDER}
    gas_epoxide = cached_ws["F60"].value
    # DFT_TS_relative.csv explicitly documents source rows 12/16/20 and 24/28/32.
    # The worksheet may put a total/reference adjustment in the following row or require
    # removal of the epoxide spectator for TiMo; therefore semantic row numbers are not
    # inferred from the unrelated state-index column.
    anchors = {
        "NEB1": (12, 16, 20, "G_IS1", "G_TS1", "G_FS1_total", "dG_TS1_minus_IS1", "dG_TS1_minus_FS1"),
        "NEB2": (24, 28, 32, "G_IS2", "G_TS2", "G_FS2_adjusted", "dG_TS2_minus_IS2", "dG_TS2_minus_FS2"),
    }
    for catalyst in CATALYST_ORDER:
        rel = relative.get(catalyst, {})
        for path_name, (r_is, r_ts, r_fs, k_is, k_ts, k_fs, k_fwd, k_rev) in anchors.items():
            values = [rel.get(k) for k in (k_is, k_ts, k_fs)]
            if any(v in (None, "") for v in values):
                rows.append({"catalyst": catalyst, "path": path_name, "status": "missing relative CSV energy"})
                continue
            e_is, e_ts, e_fs = map(float, values)
            barriers = barrier_values(model, e_is, e_ts, e_fs)
            cols = layouts[catalyst]

            def mapped_source(row_number: int, target: float, role: str) -> tuple[str, str, float | None]:
                gcell = formulas_ws.cell(row_number, cols["G"])
                raw = cached_ws[gcell.coordinate].value
                if isinstance(raw, (int, float)) and abs(float(raw) - target) < 1e-6:
                    return gcell.coordinate, str(gcell.value), float(raw)
                # NEB1 FS total may explicitly append gas epoxide one row below the surface state.
                if path_name == "NEB1" and role == "FS":
                    for candidate in range(row_number, min(row_number + 2, formulas_ws.max_row) + 1):
                        c = formulas_ws.cell(candidate, cols["G"])
                        cv = cached_ws[c.coordinate].value
                        if isinstance(cv, (int, float)) and abs(float(cv) - target) < 1e-6:
                            return c.coordinate, str(c.value), float(cv)
                # TiMo NEB2 FS contains an epoxide spectator in the workbook total;
                # the relative-energy table removes precisely that explicit gas reference.
                if path_name == "NEB2" and role == "FS" and isinstance(raw, (int, float)):
                    if abs((float(raw) - float(gas_epoxide)) - target) < 1e-6:
                        return gcell.coordinate, f"({gcell.value}) - $F$60 [C6H12O_g spectator removed]", float(raw) - float(gas_epoxide)
                return gcell.coordinate, str(gcell.value), float(raw) if isinstance(raw, (int, float)) else None

            is_cell, is_expr, is_mapped = mapped_source(r_is, e_is, "IS")
            ts_cell, ts_expr, ts_mapped = mapped_source(r_ts, e_ts, "TS")
            fs_cell, fs_expr, fs_mapped = mapped_source(r_fs, e_fs, "FS")
            source_status = "MATCH" if all(v is not None and abs(v - tgt) < 1e-6 for v, tgt in
                                              ((is_mapped, e_is), (ts_mapped, e_ts), (fs_mapped, e_fs))) else "MISMATCH_REQUIRES_REVIEW"
            rows.append({
                "catalyst": catalyst, "path": path_name, "status": "CSV relative state energies evaluated against source workbook",
                "IS_row": r_is, "TS_row": r_ts, "FS_row": r_fs,
                "IS_source_cell": is_cell, "IS_source_expression": is_expr,
                "TS_source_cell": ts_cell, "TS_source_expression": ts_expr,
                "FS_source_cell": fs_cell, "FS_source_expression_or_adjustment": fs_expr,
                "source_energy_crosswalk_status": source_status,
                "G_IS_workbook_or_normalized": e_is, "G_TS_workbook_or_normalized": e_ts,
                "G_FS_workbook_or_normalized": e_fs,
                "raw_TS_minus_IS": barriers["raw_TS_minus_IS"],
                "raw_TS_minus_FS": barriers["raw_TS_minus_FS"],
                "raw_FS_minus_IS": barriers["deltaG_forward"],
                "CatMAP_effective_Ea_forward": barriers["catmap_Ea_forward"],
                "CatMAP_effective_Ea_reverse": barriers["catmap_Ea_reverse"],
                "CSV_TS_minus_IS": float(rel[k_fwd]), "CSV_TS_minus_FS": float(rel[k_rev]),
                "workbook_mapping_rule": "cell-level cached value/formula matched; any catalyst-specific spectator adjustment is explicit",
                "state_formula_reference_note": rel.get("reference_note", ""),
                "energy_source_note": "signed relative CSV values retained; effective barriers computed with local CatMAP 0.3.1",
            })
    return rows


def definition_rows() -> list[dict[str, Any]]:
    rows = []
    for name, role, comp, evidence, sites in [
        ("NEB1_IS_s", "NEB1 initial coadsorbed state", {"C": 6, "H": 14, "O": 2},
         "workbook umbrella label H2O2+C6H12; Ti-Re path text includes 1-hex+ooh+h; coordinates unavailable", None),
        ("NEB1_TS_s", "NEB1 transition-state configuration", {"C": 6, "H": 14, "O": 2},
         "same-NEB atom inventory inferred from endpoints; transition-state coordinates unavailable", None),
        ("NEB1_FS_s", "NEB1 final coadsorbed product state", {"C": 6, "H": 14, "O": 2},
         "workbook label C6H12O+H2O coadsorbed; coordinates unavailable", None),
        ("NEB2_IS_s", "NEB2 initial separated surface H/O state", {"H": 2, "O": 1},
         "workbook label Hb*+O*+Ha*; gas epoxide is a spectator reference for cross-path comparison", None),
        ("NEB2_TS_s", "NEB2 transition-state configuration", {"H": 2, "O": 1},
         "same-NEB atom inventory inferred from endpoints; transition-state coordinates unavailable", None),
        ("NEB2_FS_s", "NEB2 final OH+H state", {"H": 2, "O": 1},
         "workbook label OH*+Ha*; gas epoxide is a spectator reference for cross-path comparison", None),
    ]:
        rows.append({
            "species_key_candidate": name,
            "role": role,
            "element_composition_candidate": json.dumps(comp, sort_keys=True),
            "composition_evidence": evidence,
            "n_sites": sites,
            "n_sites_status": "unknown; must be derived from original NEB coordinates/local site definition",
            "CatMAP_syntax_status": "representable with explicit composition and n_sites; physical definition not validated",
            "stable_or_TS_status": "TS species must only occur in the explicit TS position; not a stable adsorbate",
        })
    return rows


def network_rows(tire_mapping: list[dict[str, Any]], tire_states: dict[int, dict[str, Any]], gas_refs: dict[str, float]) -> list[dict[str, Any]]:
    lookup = {r["path"]: r for r in tire_mapping if r["catalyst"] == "tire"}
    neb1, neb2 = lookup["NEB1"], lookup["NEB2"]
    neb1_fs_total = float(neb1["G_FS_workbook_or_normalized"])
    neb2_is = float(neb2["G_IS_workbook_or_normalized"])
    g6_total = neb2_is + float(gas_refs["C6H12O_g"])
    connect_dg = g6_total - neb1_fs_total
    g8_total = float(tire_states[8]["G"]) + float(gas_refs["C6H12O_g"])
    # State 9 is an adsorbed-water endpoint plus gas epoxide; state 10 is gas products on clean-slab reference.
    g9_total = float(tire_states[9]["G"]) + float(gas_refs["C6H12O_g"])
    g10_total = float(tire_states[10]["G"]) + float(gas_refs["H2O_g"]) + float(gas_refs["C6H12O_g"])
    return [
        {"scheme": "A_minimal_NEB1", "step": "feed_to_NEB1_IS", "candidate_stoichiometry": "C6H12_g + H2O2_g + n* -> NEB1_IS_s",
         "atom_balance": "nominal C6H14O2 on both sides", "site_balance": "unknown n; NEB1_IS occupancy unverified",
         "energy_or_barrier": "no adsorption/activation TS supplied", "evidence": "reservoir species + NEB1_IS composition inferred from label/path",
         "status": "incomplete; do not assign a rate or TOF"},
        {"scheme": "A_minimal_NEB1", "step": "NEB1_internal", "candidate_stoichiometry": "NEB1_IS_s <-> NEB1-TS_s -> NEB1_FS_s",
         "atom_balance": "C6H14O2 to C6H14O2, label-derived", "site_balance": "same NEB cell likely conserves sites; CatMAP n_sites still unverified",
         "energy_or_barrier": f"Ti-Re raw TS-IS={neb1['raw_TS_minus_IS']:.6f} eV; TS-FS={neb1['raw_TS_minus_FS']:.6f} eV",
         "evidence": "direct NEB1 endpoint/TS Gibbs energies in workbook", "status": "direct NEB step; per-species absolute formation energy not unique"},
        {"scheme": "A_minimal_NEB1", "step": "NEB1_FS_to_products", "candidate_stoichiometry": "NEB1_FS_s -> C6H12O_g + H2O_g + n*",
         "atom_balance": "nominal balanced", "site_balance": "depends on NEB1_FS site occupancy n",
         "energy_or_barrier": f"Ti-Re endpoint total ΔG to state10+gas refs={g10_total-neb1_fs_total:.6f} eV; desorption TS missing",
         "evidence": "workbook state5, state10 and gas reference rows", "status": "thermodynamic endpoint comparison only; no kinetic barrier"},
        {"scheme": "B_serial_NEB1_to_NEB2", "step": "NEB1_FS_to_NEB2_IS", "candidate_stoichiometry": "NEB1_FS_s + n* -> C6H12O_g + NEB2_IS_s",
         "atom_balance": "nominal C6H14O2 conserved if epoxide is included as spectator gas", "site_balance": "unknown; NEB1 complex occupancy vs three separated H/O adsorbates not coordinate-verified",
         "energy_or_barrier": f"Ti-Re candidate ΔG={connect_dg:.6f} eV from G(state6)+G(C6H12O_g)-G(state5); connection TS missing",
         "evidence": "state5 label C6H12O+H2O*, state6 label Hb*+O*+Ha*, Excel gas reference F60", "status": "priority hypothesis only; no zero-barrier assumption"},
        {"scheme": "B_serial_NEB1_to_NEB2", "step": "NEB2_internal", "candidate_stoichiometry": "Ha_s + Hb_s + O_s <-> Ha_s + Hb-O_s + * -> OH_s + Ha_s + *",
         "atom_balance": "H2O conserved, label-derived", "site_balance": "3 sites each side under monodentate label interpretation; must verify coordinates",
         "energy_or_barrier": f"Ti-Re raw TS-IS={neb2['raw_TS_minus_IS']:.6f} eV; TS-FS={neb2['raw_TS_minus_FS']:.6f} eV",
         "evidence": "direct NEB2 endpoint/TS Gibbs energies; extra gas epoxide cancels within this path", "status": "direct NEB step; occupancy interpretation provisional"},
        {"scheme": "B_serial_NEB1_to_NEB2", "step": "NEB2_FS_to_adsorbed_water_state9", "candidate_stoichiometry": "OH_s + Ha_s -> H2O_s + *",
         "atom_balance": "H2O conserved", "site_balance": "2-to-2 under monodentate interpretation", 
         "energy_or_barrier": f"Ti-Re endpoint ΔG={g9_total-g8_total:.6f} eV; no connecting TS in supplied NEB data",
         "evidence": "workbook state8, state9 and gas epoxide reference", "status": "candidate thermodynamic connection; rate barrier missing"},
        {"scheme": "B_serial_NEB1_to_NEB2", "step": "water_desorption", "candidate_stoichiometry": "H2O_s -> H2O_g + *",
         "atom_balance": "H2O conserved", "site_balance": "balanced for one-site H2O_s", 
         "energy_or_barrier": f"Ti-Re endpoint ΔG={g10_total-g9_total:.6f} eV; desorption TS missing",
         "evidence": "workbook state9/state10 + gas references", "status": "endpoint thermodynamics only; no kinetic barrier"},
    ]


def write_manual_xlsx(path: Path, states, gas_refs, mapping_rows_all, tire_network) -> None:
    wb = Workbook()
    ws = wb.active
    ws.title = "TiRe_NEB核算"
    ws.append(["NEB路径", "状态", "工作簿行", "G单元格", "G公式", "缓存G/eV", "来源计算路径"])
    for path_name, state_ids in (("NEB1", [3, 4, 5]), ("NEB2", [6, 7, 8])):
        for idx in state_ids:
            st = states["tire"][idx]
            ws.append([path_name, idx, st["row"], st["G_cell"], st["G_formula"], st["G"], st["path"]])
    ws.append([])
    ws.append(["路径", "G_FS-G_IS", "TS-G_IS", "TS-G_FS", "CatMAP正向Ea", "CatMAP逆向Ea"])
    for r in mapping_rows_all:
        if r.get("catalyst") == "tire":
            ws.append([r["path"], r["raw_FS_minus_IS"], r["raw_TS_minus_IS"], r["raw_TS_minus_FS"],
                       r["CatMAP_effective_Ea_forward"], r["CatMAP_effective_Ea_reverse"]])
    ws2 = wb.create_sheet("连接态与参考"); ws2.append(["量", "数值/eV", "来源/定义"])
    ws2.append(["NEB1 state5 G", states["tire"][5]["G"], states["tire"][5]["G_cell"]])
    ws2.append(["NEB2 state6 G", states["tire"][6]["G"], states["tire"][6]["G_cell"]])
    ws2.append(["C6H12O gas G reference", gas_refs["C6H12O_g"], "吉布斯自由能汇总!F60"])
    ws2.append(["state6+spectator epoxide total G", states["tire"][6]["G"] + gas_refs["C6H12O_g"], "state6 cached G + F60"])
    ws2.append(["candidate ΔG state5→state6+epoxide_g", states["tire"][6]["G"] + gas_refs["C6H12O_g"] - states["tire"][5]["G"], "not a barrier; assumes common slab/reference and complete spectator accounting"])
    ws2.append(["state8+epoxide_g total G", states["tire"][8]["G"] + gas_refs["C6H12O_g"], "state8 cached G + F60"])
    ws2.append(["state9+epoxide_g total G", states["tire"][9]["G"] + gas_refs["C6H12O_g"], "state9 cached G + F60"])
    ws2.append(["state10+epoxide_g+water_g total G", states["tire"][10]["G"] + gas_refs["C6H12O_g"] + gas_refs["H2O_g"], "state10 cached G + F60 + F59"])
    ws2.append(["连接证据限制", "not validated", "raw NEB coordinates unavailable; cell composition/site count/reference equivalence must be verified"])
    for sheet in wb.worksheets:
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for col in sheet.columns:
            letter = col[0].column_letter
            sheet.column_dimensions[letter].width = min(48, max(14, max(len(str(c.value or "")) for c in col) + 2))
    path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)


def build(workbook: Path, ts_csv: Path, output_dir: Path) -> None:
    data = output_dir / "data"
    reports = output_dir / "reports"
    states, gas_refs = read_cached_state_data(workbook)
    relative = load_relative(ts_csv)
    mrows = mapping_rows(states, relative, workbook)
    mapping_path = data / "neb_to_catmap_energy_mapping.csv"
    write_csv(mapping_path, mrows, list(mrows[0]))
    defs = definition_rows()
    write_csv(data / "proposed_adsorbate_definitions.csv", defs, list(defs[0]))
    tire_rows = [r for r in mrows if r.get("catalyst") == "tire"]
    network = network_rows(tire_rows, states["tire"], gas_refs)
    write_csv(data / "proposed_reaction_network.csv", network, list(network[0]))
    write_manual_xlsx(reports / "ti_re_manual_energy_validation.xlsx", states, gas_refs, mrows, network)
    tire_connection = next(r for r in network if r["step"] == "NEB1_FS_to_NEB2_IS")
    summary = {
        "input_workbook_sha256": sha256(workbook),
        "input_neb_csv_sha256": sha256(ts_csv),
        "catmap_version": __import__("catmap").__version__,
        "gas_reference_cached_values": gas_refs,
        "mapping_rows": len(mrows),
        "mapping_csv": str(mapping_path),
        "tire_neb1_neb2_connection_candidate_deltaG_eV": tire_connection["energy_or_barrier"],
        "tof_calculated": False,
    }
    (data / "energy_mapping_manifest.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--ts-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    build(args.workbook, args.ts_csv, args.output_dir)
