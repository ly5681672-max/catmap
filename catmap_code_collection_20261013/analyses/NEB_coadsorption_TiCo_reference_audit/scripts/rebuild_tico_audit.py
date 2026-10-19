#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read-only source audit for Ti-Co's two full-system NEB pathways.

Rebuild CSVs from existing repository DFT/NEB source tables. When the raw
workbook is accessible, additionally locate its Ti-Co G columns dynamically
and verify the actual cached Excel cell values. Does not modify any source
or call/run CatMAP.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[4]
SOURCE_DIR = REPO / "lyq_alkene_epoxidation"
MKM = SOURCE_DIR / "2-creating_microkinetic_model"
WORKBOOK = SOURCE_DIR / "己烯环氧化.xlsx"
CONVERTER = SOURCE_DIR / "1-generating_input_file" / "convert_dft_to_catmap.py"
DATA = HERE / "data"
TOL = Decimal("0.000001")
EXPECTED_ATOM_MU = {"C": Decimal("-9.28"), "H": Decimal("-1.11"), "O": Decimal("-4.37")}
COMPOSITION = {
    "Ha": {"H": 1}, "Hb": {"H": 1}, "OOH": {"O": 2, "H": 1},
    "C6H12": {"C": 6, "H": 12}, "C6H12O": {"C": 6, "H": 12, "O": 1},
    "O": {"O": 1}, "OH": {"H": 1, "O": 1},
}
SINGLE_EXCEL_ROW = {"Ha": 8, "Hb": 12, "OOH": 16, "C6H12": 20,
                    "C6H12O": 24, "O": 29, "OH": 34}
ROUTES = {
    "NEB1": {
        "source_cols": ("G_IS1", "G_TS1", "G_FS1_total"),
        "source_rows": {"IS": 12, "TS": 16, "FS": 20},
        "composition": {"C": 6, "H": 14, "O": 2}, "sites": 4,
        "independent": {"IS": ("Ha", "OOH", "C6H12"),
                        "FS": ("Ha", "C6H12O", "O", "Hb")},
        "formula": {"IS": "BN8+BN16+BN20-2*BN51",
                    "FS": "BN8+BN24+BN29+BN12-3*BN51"},
        "model_ts": "OOH-C6H12", "model_is": ("OOH", "C6H12"),
    },
    "NEB2": {
        "source_cols": ("G_IS2", "G_TS2", "G_FS2_adjusted"),
        "source_rows": {"IS": 24, "TS": 28, "FS": 32},
        "composition": {"H": 2, "O": 1}, "sites": 3,
        "independent": {"IS": ("Ha", "Hb", "O"), "FS": ("Ha", "OH")},
        "formula": {"IS": "BN8+BN12+BN29-2*BN51",
                    "FS": "BN8+BN34-BN51"},
        "model_ts": "Hb-O", "model_is": ("Hb", "O"),
    },
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def decimal(x: object) -> Decimal:
    return Decimal(str(x))


def dstr(x: Decimal | None) -> str:
    return "" if x is None else format(x, "f")


def assert_equal(actual: Decimal, expected: Decimal, label: str) -> None:
    if abs(actual - expected) > TOL:
        raise AssertionError(f"{label}: actual={actual}; expected={expected}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def full_mu(comp: dict[str, int]) -> Decimal:
    return sum((EXPECTED_ATOM_MU[a] * n for a, n in comp.items()), Decimal(0))


def inspect_original_workbook(
    slab: Decimal, formation: dict[str, Decimal], neb: dict[str, dict[str, Decimal]]
) -> tuple[dict[str, dict[str, str]], str]:
    """Use existing standard-library XLSX reader; check cached G, not formulas."""
    if not WORKBOOK.exists():
        return {}, "RAW_WORKBOOK_MISSING: derived CSV/TXT comparison only"
    spec = importlib.util.spec_from_file_location("existing_dft_converter", CONVERTER)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load the repository's existing converter")
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)
    for atom, val in EXPECTED_ATOM_MU.items():
        assert_equal(decimal(converter.MU[atom]), val, f"converter atomic reference {atom}")

    tables = converter.load_sheets(WORKBOOK)
    single_table = tables["catmap数据集"]
    neb_table = tables["吉布斯自由能汇总"]
    single_col = converter.surfaces(single_table)["tico"][1]
    neb_col = converter.surfaces(neb_table)["tico"][1]
    if converter.cell(single_col, 4) != "BN4":
        raise AssertionError("Source single-adsorbate Ti-Co G column is not BN")

    slab_ref = converter.cell(single_col, 51)
    assert_equal(decimal(single_table[slab_ref]), slab, slab_ref)
    for species, row in SINGLE_EXCEL_ROW.items():
        cell = converter.cell(single_col, row)
        expected_full = slab + full_mu(COMPOSITION[species]) + formation[species]
        assert_equal(decimal(single_table[cell]), expected_full, f"catmap数据集!{cell}: {species}")

    verified_cells: dict[str, dict[str, str]] = {}
    for name, route in ROUTES.items():
        verified_cells[name] = {}
        for status, row in route["source_rows"].items():
            cell = converter.cell(neb_col, row)
            assert_equal(decimal(neb_table[cell]), neb[name][status],
                         f"吉布斯自由能汇总!{cell}: {name}/{status}")
            verified_cells[name][status] = cell

    return verified_cells, f"RAW_WORKBOOK_VERIFIED: TiCo NEB G column={converter.cell(neb_col, 1)[:-1]}"


def write_csv(path: Path, headers: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


def main() -> None:
    energy_file = MKM / "energies_dft_neb.txt"
    relative_file = MKM / "DFT_TS_relative.csv"
    adsorption_file = MKM / "DFT_adsorption.csv"
    with energy_file.open("r", encoding="utf-8", newline="") as handle:
        energies = list(csv.DictReader(handle, delimiter="\t"))
    formation = {r["species_name"]: decimal(r["formation_energy"])
                 for r in energies if r["surface_name"] == "tico"}
    relative = next(r for r in read_csv(relative_file) if r["surface"] == "tico")
    adsorption = next(r for r in read_csv(adsorption_file)
                      if r["surface"] == "tico" and r["species"] == "H2O2")
    slab = decimal(adsorption["G_slab_eV"])
    neb = {key: dict(zip(("IS", "TS", "FS"),
                        (decimal(relative[c]) for c in route["source_cols"])))
           for key, route in ROUTES.items()}
    cells, workbook_status = inspect_original_workbook(slab, formation, neb)

    state_headers = ["path", "state", "original_workbook_sheet",
                     "source_excel_semantic_row", "source_excel_G_cell",
                     "NEB_full_composition", "formal_site_count_assumption",
                     "G_NEB_full_eV", "G_clean_slab_eV",
                     "atomic_mu_reference_total_eV", "G_NEB_formation_eV",
                     "G_independent_formation_sum_eV",
                     "G_independent_full_equivalent_eV",
                     "NEB_minus_independent_eV",
                     "independent_abs_G_excel_expression_literally",
                     "raw_workbook_check_status", "physical_site_structure_status"]
    comp_headers = ["path", "mechanism_step", "NEB_TS_minus_IS_eV",
                    "NEB_FS_minus_IS_eV", "independent_FS_minus_IS_eV",
                    "independent_minus_NEB_reaction_eV",
                    "NEB_TS_minus_FS_eV", "current_CatMAP_TS_formation_eV",
                    "current_CatMAP_TS_minus_IS_eV",
                    "current_CatMAP_TS_minus_FS_eV",
                    "formal_CatMAP_forward_barrier_eV", "physical_mapping_status"]
    state_rows: list[list[str]] = []
    comp_rows: list[list[str]] = []
    for name, route in ROUTES.items():
        full_ref = full_mu(route["composition"])
        independent = {
            stage: sum((formation[s] for s in route["independent"][stage]), Decimal(0))
            for stage in ("IS", "FS")
        }
        for stage in ("IS", "TS", "FS"):
            gf = neb[name][stage]
            gform = gf - slab - full_ref
            gind = independent.get(stage)
            original_cell = cells.get(name, {}).get(stage)
            source_cell = (f"吉布斯自由能汇总!{original_cell}" if original_cell
                           else f"TiCo区块G列·第{route['source_rows'][stage]}行（列字母待原工作簿复核）")
            state_rows.append([
                name, stage, "吉布斯自由能汇总", str(route["source_rows"][stage]),
                source_cell, ("C6H14O2" if name == "NEB1" else "H2O"),
                f"{route['sites']} (形式位点模型，非结构实证)",
                dstr(gf), dstr(slab), dstr(full_ref), dstr(gform), dstr(gind),
                dstr(slab + full_ref + gind if gind is not None else None),
                dstr(gform - gind if gind is not None else None),
                "" if stage == "TS" else "'=" + route["formula"][stage],
                workbook_status, "NEB真实局域位点占用/同一slab参考尚未验证",
            ])

        neb_dg = neb[name]["FS"] - neb[name]["IS"]
        independent_dg = independent["FS"] - independent["IS"]
        ts = formation[route["model_ts"]]
        model_is = sum((formation[k] for k in route["model_is"]), Decimal(0))
        model_fs = independent["FS"] - formation["Ha"]
        model_ts_f = ts - model_is
        model_ts_r = ts - model_fs
        forward_barrier = max(Decimal(0), independent_dg, model_ts_f)
        comp_rows.append([
            name, ("4" if name == "NEB1" else "6"),
            dstr(neb[name]["TS"] - neb[name]["IS"]), dstr(neb_dg),
            dstr(independent_dg), dstr(independent_dg - neb_dg),
            dstr(neb[name]["TS"] - neb[name]["FS"]), dstr(ts),
            dstr(model_ts_f), dstr(model_ts_r), dstr(forward_barrier),
            "原始NEB含Ha；CatMAP独立态均值参考尚未等价；不得直接改输入",
        ])

    DATA.mkdir(parents=True, exist_ok=True)
    write_csv(DATA / "tico_neb_state_alignment.csv", state_headers, state_rows)
    write_csv(DATA / "tico_neb_path_comparison.csv", comp_headers, comp_rows)
    source_paths = [energy_file, relative_file, adsorption_file, CONVERTER]
    if WORKBOOK.exists():
        source_paths.append(WORKBOOK)
    manifest = {
        "status": workbook_status,
        "sources_sha256": {str(p.relative_to(REPO)): sha256(p) for p in source_paths},
        "assumptions": [
            "TiCo clean slab reference is shared with NEB full systems (pending cell/supercell verification)",
            "NEB1 and NEB2 full atom inventories and Ha spectator were confirmed by researcher",
            "Formal site counts are bookkeeping, not validated physical active center counts",
            "No separate DFT TS candidate energy was substituted for true coadsorbed NEB TS",
        ],
        "no_catmap_run": True,
        "no_source_modification": True,
    }
    (DATA / "tico_source_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for row in comp_rows:
        print(f"{row[0]}: NEB ΔG={row[3]}, independent ΔG={row[4]}, gap={row[5]} eV")
    print(workbook_status)
    if not WORKBOOK.exists():
        print("WARNING: raw Excel workbook unavailable; single-cell validation NOT completed.")


if __name__ == "__main__":
    main()
