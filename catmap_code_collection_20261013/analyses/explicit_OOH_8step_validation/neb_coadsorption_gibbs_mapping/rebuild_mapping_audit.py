#!/usr/bin/env python3
"""Rebuild a read-only NEB/coadsorption-to-independent-state Gibbs audit.

This script reads frozen project CSV/TXT inputs only. It never invokes a
CatMAP solver and never modifies scientific inputs. All arithmetic uses
Decimal to preserve source decimal strings.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from decimal import Decimal
from pathlib import Path
from datetime import datetime, timezone


TARGET_SURFACES = [
    "ti", "tico", "ticr", "tife", "tihf", "timn", "timo", "tini",
    "tire", "tita", "titi", "tiv", "tiw", "tizr",
]
CATALYST_NAMES = {
    "ti": "Ti", "tico": "Ti-Co", "ticr": "Ti-Cr", "tife": "Ti-Fe",
    "tihf": "Ti-Hf", "timn": "Ti-Mn", "timo": "Ti-Mo", "tini": "Ti-Ni",
    "tire": "Ti-Re", "tita": "Ti-Ta", "titi": "Ti-Ti", "tiv": "Ti-V",
    "tiw": "Ti-W", "tizr": "Ti-Zr",
}
ZERO = Decimal("0")
TOL = Decimal("0.000001")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: str(row.get(k, "")) for k in fields})


def d(value: str | Decimal) -> Decimal:
    return value if isinstance(value, Decimal) else Decimal(value.strip())


def fmt(x: Decimal) -> str:
    return format(x, "f")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()
    root, out = args.project_root.resolve(), args.output_dir.resolve()
    analysis = root / "catmap_code_collection_20261013" / "analyses" / "explicit_OOH_8step_validation"
    neb_csv = root / "lyq_alkene_epoxidation" / "2-creating_microkinetic_model" / "DFT_TS_relative.csv"
    energy_txt = root / "lyq_alkene_epoxidation" / "2-creating_microkinetic_model" / "energies_dft_neb.txt"
    cell_audit = analysis / "data" / "source_audit" / "tico_excel_cell_audit.csv"
    state_alignment = analysis / "data" / "source_audit" / "tico_neb_state_alignment.csv"
    neb_barrier_audit = analysis / "multisite_batch_8step" / "MULTISITE_REACTION_BARRIERS.csv"
    catmap_barriers = analysis / "conditional_DRC_OOH_TOF_14sites" / "energy_rate_DRC_mechanism_audit" / "data" / "EIGHT_STEP_EFFECTIVE_BARRIERS.csv"
    comparison = analysis / "data" / "source_audit" / "tico_neb_path_comparison.csv"
    source_manifest = analysis / "data" / "source_audit" / "tico_excel_source_manifest.json"
    workbook = root / "lyq_alkene_epoxidation" / "己烯环氧化.xlsx"
    formal_model = root / "lyq_alkene_epoxidation" / "2-creating_microkinetic_model" / "alkene_epoxidation.mkm"
    required = [neb_csv, energy_txt, cell_audit, state_alignment, neb_barrier_audit, catmap_barriers, comparison, source_manifest, workbook, formal_model]
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        raise SystemExit("Missing required frozen source(s):\n" + "\n".join(missing))
    frozen_hashes = {
        str(formal_model): "86115cefb7784fe28ab519271cdfaf0e308b6f6e3a8308b1e90e4318cb497cfe",
        str(energy_txt): "55ad4e94d9c603b6705ec8c92630c9a1ceff5859cf69865145152375e203ca30",
    }
    for source, expected_hash in frozen_hashes.items():
        if sha256(Path(source)) != expected_hash:
            raise SystemExit(f"Frozen formal input hash differs from approved baseline: {source}")
    out.mkdir(parents=True, exist_ok=True)

    neb_rel = {r["surface"].strip().lower(): r for r in read_csv(neb_csv)}
    excel_audit = {
        r["label"].strip(): r for r in read_csv(cell_audit)
        if r.get("record_type") in {"neb_state", "independent_adsorbate", "clean_slab"}
    }
    energies: dict[str, dict[str, Decimal]] = {}
    for r in read_tsv(energy_txt):
        s, species = r["surface_name"].strip().lower(), r["species_name"].strip()
        if s in TARGET_SURFACES:
            energies.setdefault(s, {})[species] = d(r["formation_energy"])
    raw_rows = {r["surface_name"].strip().lower() + ":" + r["path"].strip(): r for r in read_csv(neb_barrier_audit)}
    catmap_rows = {(r["surface_name"].strip().lower(), r["reaction_index_1based"].strip()): r for r in read_csv(catmap_barriers)}
    cross_source_checks = 0

    tico_rows = []
    mapping_rows = []
    site_rows = []
    for s in TARGET_SURFACES:
        nr = neb_rel.get(s)
        en = energies.get(s, {})
        has_ind = all(x in en for x in ["OOH", "C6H12", "C6H12O", "O", "Hb", "OH"])
        for path, reaction_index, is_species, fs_species, ts_species, n_is, n_ts, n_fs in [
            ("NEB1", "4", ("OOH", "C6H12"), ("C6H12O", "O", "Hb"), "OOH-C6H12", "G_IS1", "G_TS1", "G_FS1_total"),
            ("NEB2", "6", ("Hb", "O"), ("OH",), "Hb-O", "G_IS2", "G_TS2", "G_FS2_adjusted"),
        ]:
            raw = raw_rows.get(f"{s}:{path}")
            cm = catmap_rows.get((s, reaction_index))
            neb_is = d(nr[n_is]) if nr else None
            neb_ts = d(nr[n_ts]) if nr else None
            neb_fs = d(nr[n_fs]) if nr else None
            neb_dg = neb_fs - neb_is if nr else None
            neb_af = neb_ts - neb_is if nr else None
            neb_ar = neb_ts - neb_fs if nr else None
            independent_dg = sum((en[x] for x in fs_species), ZERO) - sum((en[x] for x in is_species), ZERO) if has_ind else None
            ref_delta = independent_dg - neb_dg if (has_ind and nr) else None
            # A forward-barrier-preserving shift is uniquely determined by G_IS_ind.
            # Its reverse barrier necessarily differs by -ref_delta.
            mapped_reverse = neb_af - independent_dg if (has_ind and nr) else None
            reverse_error = mapped_reverse - neb_ar if mapped_reverse is not None else None
            if not nr or not raw:
                status = "INSUFFICIENT_DATA"
            elif not has_ind:
                status = "INSUFFICIENT_DATA"
            elif abs(ref_delta) > TOL:
                status = "REFERENCE_MISMATCH"
            else:
                status = "MAPPING_CONSISTENT"
            if nr and raw:
                assert abs(d(raw["signed_TS_minus_IS_eV"]) - neb_af) <= TOL, (s, path, "TS-IS difference")
                assert abs(d(raw["signed_TS_minus_FS_eV"]) - neb_ar) <= TOL, (s, path, "TS-FS difference")
                assert abs((d(raw["G_FS_eV"]) - d(raw["G_IS_eV"])) - neb_dg) <= TOL, (s, path, "FS-IS difference")
                assert abs((neb_af - neb_ar) - neb_dg) <= TOL, (s, path, "NEB thermodynamic identity")
                cross_source_checks += 1
            if cm and independent_dg is not None:
                assert abs(d(cm["deltaG_rxn_eV"]) - independent_dg) <= TOL, (s, path, "official CatMAP-derived independent reaction energy")
            if s == "tico":
                neb_labels = {"NEB1": ("NEB1_IS", "NEB1_TS", "NEB1_FS"), "NEB2": ("NEB2_IS", "NEB2_TS", "NEB2_FS")}[path]
                independent_labels = list(is_species) + list(fs_species)
                neb_sources = [excel_audit.get(label, {}) for label in neb_labels]
                independent_sources = [excel_audit.get(label, {}) for label in independent_labels]
                # For tico retain detailed state energies, source cells and actual
                # current CatMAP endpoint/TS semantics, where audited.
                tico_rows.append({
                    "path": path, "mechanism_step": reaction_index,
                    "NEB_IS_G_eV": fmt(neb_is) if neb_is is not None else "",
                    "NEB_TS_G_eV": fmt(neb_ts) if neb_ts is not None else "",
                    "NEB_FS_G_eV": fmt(neb_fs) if neb_fs is not None else "",
                    "NEB_energy_basis": "Excel absolute Gibbs total; reaction differences are reference-invariant within each path",
                    "NEB_deltaG_FS_minus_IS_eV": fmt(neb_dg) if neb_dg is not None else "",
                    "NEB_forward_barrier_TS_minus_IS_eV": fmt(neb_af) if neb_af is not None else "",
                    "NEB_reverse_barrier_TS_minus_FS_eV": fmt(neb_ar) if neb_ar is not None else "",
                    "independent_G_IS_sum_eV": fmt(sum((en[x] for x in is_species), ZERO)) if has_ind else "",
                    "independent_G_TS_current_eV": fmt(en[ts_species]) if ts_species in en else "",
                    "independent_G_FS_sum_eV": fmt(sum((en[x] for x in fs_species), ZERO)) if has_ind else "",
                    "independent_energy_basis": "sum of CatMAP formation_energy values; vacancy/site reference cancels by reaction stoichiometry",
                    "independent_deltaG_eV": fmt(independent_dg) if independent_dg is not None else "",
                    "NEB_excel_sheet": "吉布斯自由能汇总",
                    "NEB_excel_G_cells_IS_TS_FS": ";".join(x.get("G_cell", "") for x in neb_sources),
                    "NEB_excel_G_formulas_IS_TS_FS": ";".join(x.get("G_formula", "") for x in neb_sources),
                    "NEB_excel_G_cached_values_IS_TS_FS": ";".join(x.get("G_cached_eV", "") for x in neb_sources),
                    "NEB_excel_calculation_paths_IS_TS_FS": ";".join(x.get("calculation_path", "") for x in neb_sources),
                    "independent_excel_sheet": ";".join(x.get("workbook_sheet", "") for x in independent_sources),
                    "independent_excel_G_cells_IS_species_then_FS_species": ";".join(x.get("G_cell", "") for x in independent_sources),
                    "independent_excel_G_formulas_IS_species_then_FS_species": ";".join(x.get("conversion_or_reference_definition", "") for x in independent_sources),
                    "independent_excel_G_cached_values_IS_species_then_FS_species": ";".join(x.get("G_cached_eV", "") for x in independent_sources),
                    "independent_species_order": ";".join(independent_labels),
                    "reference_difference_independent_minus_NEB_eV": fmt(ref_delta) if ref_delta is not None else "",
                    "current_CatMAP_effective_forward_barrier_eV": cm.get("deltaG_forward_eff_eV", "") if cm else "",
                    "current_CatMAP_effective_reverse_barrier_eV": cm.get("deltaG_reverse_eff_eV", "") if cm else "",
                    "current_CatMAP_control": cm.get("effective_control", "") if cm else "",
                    "source_NEBCSV": str(neb_csv), "source_independent_energy": str(energy_txt),
                    "mapping_status": status,
                })
            mapping_rows.append({
                "surface_name": s, "catalyst": CATALYST_NAMES[s], "path": path,
                "NEB_deltaG_eV": fmt(neb_dg) if neb_dg is not None else "",
                "independent_state_deltaG_eV": fmt(independent_dg) if independent_dg is not None else "",
                "independent_minus_NEB_deltaG_eV": fmt(ref_delta) if ref_delta is not None else "",
                "NEB_forward_barrier_eV": fmt(neb_af) if neb_af is not None else "",
                "NEB_reverse_barrier_eV": fmt(neb_ar) if neb_ar is not None else "",
                "mapped_TS_forward_barrier_eV": fmt(neb_af) if neb_af is not None else "",
                "mapped_TS_reverse_barrier_eV": fmt(mapped_reverse) if mapped_reverse is not None else "",
                "mapped_reverse_minus_NEB_reverse_eV": fmt(reverse_error) if reverse_error is not None else "",
                "TS_above_independent_IS_and_FS": (str((neb_af > ZERO) and (mapped_reverse > ZERO)).upper() if mapped_reverse is not None else ""),
                "thermodynamic_closure_NEb_af_minus_ar_equals_dG": (str(abs((neb_af-neb_ar)-neb_dg) <= TOL).upper() if neb_dg is not None else ""),
                "mapping_status": status,
                "interpretation": "Forward-preserving TS shift is a diagnostic only; not an approved energy edit." if nr and has_ind else "Required source comparison incomplete.",
            })
        site_rows.append({
            "surface_name": s, "catalyst": CATALYST_NAMES[s],
            "NEB1_data_present": str((nr is not None) and (f"{s}:NEB1" in raw_rows)).upper(),
            "NEB2_data_present": str((nr is not None) and (f"{s}:NEB2" in raw_rows)).upper(),
            "independent_energy_data_present": str(has_ind).upper(),
            "NEB1_deltaG_eV": next((r["NEB_deltaG_eV"] for r in mapping_rows if r["surface_name"] == s and r["path"] == "NEB1"), ""),
            "NEB1_independent_deltaG_eV": next((r["independent_state_deltaG_eV"] for r in mapping_rows if r["surface_name"] == s and r["path"] == "NEB1"), ""),
            "NEB1_difference_eV": next((r["independent_minus_NEB_deltaG_eV"] for r in mapping_rows if r["surface_name"] == s and r["path"] == "NEB1"), ""),
            "NEB2_deltaG_eV": next((r["NEB_deltaG_eV"] for r in mapping_rows if r["surface_name"] == s and r["path"] == "NEB2"), ""),
            "NEB2_independent_deltaG_eV": next((r["independent_state_deltaG_eV"] for r in mapping_rows if r["surface_name"] == s and r["path"] == "NEB2"), ""),
            "NEB2_difference_eV": next((r["independent_minus_NEB_deltaG_eV"] for r in mapping_rows if r["surface_name"] == s and r["path"] == "NEB2"), ""),
            "NEB1_mapping_status": next((r["mapping_status"] for r in mapping_rows if r["surface_name"] == s and r["path"] == "NEB1"), ""),
            "NEB2_mapping_status": next((r["mapping_status"] for r in mapping_rows if r["surface_name"] == s and r["path"] == "NEB2"), ""),
            "physical_cause_assignment": "UNVERIFIED: numerical difference does not identify coadsorption/slab/spectator cause without state-identity and common-reference evidence.",
        })

    # Deterministic source-data integrity checks for the two highlighted Ti-Co paths.
    tico_map = {(r["surface_name"], r["path"]): r for r in mapping_rows if r["surface_name"] == "tico"}
    expected = {
        "NEB1": (Decimal("-0.239686"), Decimal("1.800647500"), Decimal("2.040333500")),
        "NEB2": (Decimal("-0.480580"), Decimal("-2.062067000"), Decimal("-1.581487000")),
    }
    for path, (expected_neb, expected_ind, expected_diff) in expected.items():
        row = tico_map[("tico", path)]
        assert abs(d(row["NEB_deltaG_eV"]) - expected_neb) <= TOL, (path, "NEB endpoint difference")
        assert abs(d(row["independent_state_deltaG_eV"]) - expected_ind) <= TOL, (path, "independent-state reaction energy")
        assert abs(d(row["independent_minus_NEB_deltaG_eV"]) - expected_diff) <= TOL, (path, "reference difference")
        assert row["thermodynamic_closure_NEb_af_minus_ar_equals_dG"] == "TRUE", (path, "NEB detailed-balance identity")
        assert abs(d(row["mapped_reverse_minus_NEB_reverse_eV"]) + expected_diff) <= TOL, (path, "mapped reverse barrier inconsistency")
    assert len(site_rows) == 14 and len(mapping_rows) == 28

    write_csv(out / "TICO_NEB1_NEB2_GIBBS_REFERENCE.csv", tico_rows, list(tico_rows[0]))
    write_csv(out / "TS_MAPPING_CONSISTENCY.csv", mapping_rows, list(mapping_rows[0]))
    write_csv(out / "14_SITE_NEB_MAPPING_STATUS.csv", site_rows, list(site_rows[0]))

    src_hashes = {str(p): sha256(p) for p in required}
    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(root), "analysis_sources": str(analysis),
        "model": "original_explicit_OOH_8step", "temperature_K": 333.15,
        "calculation_scope": "Gibbs bookkeeping only; no CatMAP steady-state/TOF/DRC and no DFT/NEB rerun",
        "arithmetic": "Decimal on frozen source decimal strings; mapping comparison tolerance 1e-6 eV (rounding-scale check only)",
        "source_sha256": src_hashes,
        "validation_checks": {
            "NEB_DFT_TS_relative_vs_multisite_barrier_rows_matched": cross_source_checks,
            "NEB_detailed_balance_identity": "PASS for every matched path",
            "CatMAP_audited_independent_deltaG_vs_frozen_energy_sums": "PASS where archived CatMAP barrier rows are available",
            "tico_mapping_diagnostic_assertions": "PASS",
            "target_surface_count": len(site_rows),
            "mapping_rows": len(mapping_rows),
        },
        "preflight_sha256": {
            "formal_model": "86115cefb7784fe28ab519271cdfaf0e308b6f6e3a8308b1e90e4318cb497cfe",
            "energies_dft_neb_txt": "55ad4e94d9c603b6705ec8c92630c9a1ceff5859cf69865145152375e203ca30",
        },
    }
    (out / "SOURCE_RECORDS_AND_HASHES.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    n1 = mapping_rows[[r["surface_name"] for r in mapping_rows].index("tico")]
    n2 = mapping_rows[[r["surface_name"] for r in mapping_rows].index("tico") + 1]
    tico1, tico2 = tico_rows
    delta1, delta2 = d(tico1["reference_difference_independent_minus_NEB_eV"]), d(tico2["reference_difference_independent_minus_NEB_eV"])
    report_ref = f"""# Ti–Co 第4、6步参考差异

## 结论

从冻结的 CatMAP 独立态能量输入重算后，第4步独立态反应自由能为 **{tico1['independent_deltaG_eV']} eV**，完整 NEB1 端点给出 **{tico1['NEB_deltaG_FS_minus_IS_eV']} eV**，差为 **{fmt(delta1)} eV**。第6步分别为 **{tico2['independent_deltaG_eV']} eV** 与 **{tico2['NEB_deltaG_FS_minus_IS_eV']} eV**，差为 **{fmt(delta2)} eV**。

这些数值确认的是两套已输入端点反应能不同（`CONFIRMED`）；它们**不能单独证明**差异来自共吸附相互作用、slab参考、旁观 Ha 的相互作用或某一特定结构变化（均 `UNVERIFIED`）。研究者确认 Gibbs 已完成热修正；本审计接受该前提，不复查温度或频率修正。

## NEB 原始路径（eV）

表中 NEB 的 `G_IS/G_TS/G_FS` 是 Excel 汇总中的绝对 Gibbs 总能；独立态总和是 CatMAP 的 formation_energy 组合。两者绝对值不直接横向相减，审计只比较各自路径内部的能量差及两种端点定义的反应自由能。

| 路径 | IS | TS | FS | FS−IS | TS−IS | TS−FS |
|---|---:|---:|---:|---:|---:|---:|
| NEB1 / 第4步 | {tico1['NEB_IS_G_eV']} | {tico1['NEB_TS_G_eV']} | {tico1['NEB_FS_G_eV']} | {tico1['NEB_deltaG_FS_minus_IS_eV']} | {tico1['NEB_forward_barrier_TS_minus_IS_eV']} | {tico1['NEB_reverse_barrier_TS_minus_FS_eV']} |
| NEB2 / 第6步 | {tico2['NEB_IS_G_eV']} | {tico2['NEB_TS_G_eV']} | {tico2['NEB_FS_G_eV']} | {tico2['NEB_deltaG_FS_minus_IS_eV']} | {tico2['NEB_forward_barrier_TS_minus_IS_eV']} | {tico2['NEB_reverse_barrier_TS_minus_FS_eV']} |

两条路径均满足有符号恒等式 `ΔG‡f − ΔG‡r = ΔGrxn`，这是各自同一组 IS/TS/FS 数据内部的算术闭合，不证明 NEB 与独立吸附态网络端点物理同一。

## 差异来源的证据边界

- **状态/元素计量：** 既有 Excel 单元格审计把 NEB1 完整状态记为包含旁观 Ha 的 `C6H14O2`，NEB2 状态记为 `H2O`（也含 Ha）；独立 CatMAP 状态由分别输入的吸附物种形成能相加。工作簿记录了 NEB 及独立态各自的数值与 calculation_path。该记录确认状态账目不同，但没有结构一一映射证据。
- **Slab/参考：** 旧审计用 Ti–Co clean-slab 和原子化学势计算了条件性形成能；它也明确指出不同结构的 clean slab/超胞同一性未被工作簿证明。故不能把差值定量归因到 slab 项。
- **Ha 旁观物：** Ha 在 NEB1、NEB2 状态组成中被保留；没有证明其与独立 Ha* 能量可直接相消，故不能从 NEB 总能中扣除 Ha* 并宣称获得“修正势垒”。
- **共吸附：** `ΔG_ind−ΔG_NEB` 只能称参考/端点反应能差。缺少同一 slab、同一微观构型对照时，不称为 `ΔG_coadsorption`。

## 最小结论

差异的**算术事实已确认**；其唯一物理来源尚不能判定。当前证据不支持对 CatMAP TS 作平移或把 NEB 局部势垒直接套到独立态端点。最小所需证据是：可确认 NEB 与独立端点的状态同一性及共同 slab/位点参考；若不是同一状态，则需真实状态连接的自由能和连接 TS/势垒。详见 `COADSORPTION_REFERENCE_AUDIT.md` 与 `SOURCE_RECORDS_AND_HASHES.json`。
"""
    (out / "TICO_REFERENCE_DIFFERENCE.md").write_text(report_ref, encoding="utf-8")

    report_coads = f"""# 共吸附态—独立态参考审计

## 范围与证据等级

本报告只重用已归档的 Ti–Co Gibbs/Excel 审计、当前冻结 `energies_dft_neb.txt`、`DFT_TS_relative.csv` 和 CatMAP 有效势垒表；不重算 DFT/NEB/Gibbs，不复核研究者已确认的热修正，不运行 CatMAP 稳态、TOF 或 DRC。所有原始输入只读。来源哈希见 `SOURCE_RECORDS_AND_HASHES.json`。

- `CONFIRMED`：当前归档表中可直接重算的能量恒等式、reaction-energy 差及当前 CatMAP 使用的势垒数值。
- `INFERRED`：依据状态标签/化学计量记录的路径解释。
- `UNVERIFIED`：结构身份、真实 footprint、slab/超胞等价以及差异的唯一物理成因。

## 独立态能量组合

按原始八步反应计量（空位形成能取 CatMAP 格点参考零）重构：

- 第4步：`G_IS = G(OOH*) + G(C6H12*)`；`G_FS = G(C6H12O*) + G(O*) + G(Hb*)`；故 `ΔG_ind = G_FS−G_IS = {tico1['independent_deltaG_eV']} eV`。TS 当前输入为 `G(TS1*)={tico1['independent_G_TS_current_eV']} eV`。
- 第6步：`G_IS = G(Hb*) + G(O*)`；`G_FS = G(OH*)`；故 `ΔG_ind = {tico2['independent_deltaG_eV']} eV`。TS 当前输入为 `G(TS2*)={tico2['independent_G_TS_current_eV']} eV`。

对应 NEB 独立汇总的单元格来源是 `吉布斯自由能汇总!BM12/BM16/BM20`（NEB1）及 `BM24/BM28/BM32`（NEB2），逐态公式、缓存值及 calculation_path 见项目既有 `data/source_audit/tico_excel_cell_audit.csv`。NEB1、NEB2完整组成中均保留 Ha 旁观组分；此报告不把 Ha 的存在简化为“无相互作用”。

## 当前 CatMAP 有效势垒

直接重用已归档的 CatMAP 解析能量/势垒表（该表按官方模型解析反应及其 IS/TS/FS 状态和端点规则生成）：

| 步骤 | 当前独立态 ΔG | 当前正向有效势垒 | 当前逆向有效势垒 | CatMAP控制态 | 说明 |
|---|---:|---:|---:|---|---|
| 第4步 | {tico1['independent_deltaG_eV']} | {tico1['current_CatMAP_effective_forward_barrier_eV']} | {tico1['current_CatMAP_effective_reverse_barrier_eV']} | {tico1['current_CatMAP_control']} | 原始 NEB TS 在独立态 FS 以下，当前正向有效势垒由 FS 端点控制；不是 NEB1 双向势垒复现 |
| 第6步 | {tico2['independent_deltaG_eV']} | {tico2['current_CatMAP_effective_forward_barrier_eV']} | {tico2['current_CatMAP_effective_reverse_barrier_eV']} | {tico2['current_CatMAP_control']} | 当前输入的独立态 TS 高于两端点，使用该 TS 对应正逆势垒；其逆向值不等于完整 NEB2 的局部逆势垒 |

这里的“有效势垒”描述当前 CatMAP 输入/端点处理；不能把它与真实完整共吸附 NEB 路径的两端势垒混称。

## 归因判断

| 候选来源 | 状态 | 依据 |
|---|---|---|
| 端点反应自由能不一致 | CONFIRMED | NEB 与独立态的数值差分别为 {fmt(delta1)}、{fmt(delta2)} eV |
| NEB状态含旁观Ha，独立反应式未显式计入该共吸附构型 | CONFIRMED（计量标签层面） | Excel审计/状态来源表；不能据此把差值单独归因于Ha相互作用 |
| 共吸附相互作用能等于上述差值 | UNVERIFIED | 需同一slab、相同组成基准与状态映射；当前差值并非已证明的共吸附能 |
| NEB与独立吸附态的clean slab/超胞基准完全相同 | UNVERIFIED | 既有审计指出工作簿不能证明结构同源；不以形成能代数闭合代替结构证据 |
| 差异可由统一常数平移TS解决 | FAILED（数学上） | 平移TS不能改变 NEB 反应能恒等式；详见 `TS_MAPPING_CONSISTENCY.csv` |

## 尚需的最小证据

1. NEB IS/FS 与独立态组合是否为同一微观状态的可追溯对应（或明确是不同构型）。
2. 若不同：独立态进入 NEB IS、NEB FS 离开到后续稳定态的实际连接状态自由能；只有端点差不能定出动力学速率。
3. 共吸附状态的真实格点占据定义若要写入 CatMAP，须有结构支持；不得按片段数推 `n_sites`。

在这些证据齐备前，Ti–Co到其他催化剂不存在可推广的统一能量转换规则；14位点表只报告逐位点已有输入的数值对照，不施加 Ti–Co 偏移。
"""
    (out / "COADSORPTION_REFERENCE_AUDIT.md").write_text(report_coads, encoding="utf-8")
    report_final = f"""# 最终 Gibbs 映射结论

## 结论先行

1. **NEB 与原始 CatMAP 独立态网络目前不能通过单一 TS 平移完成热力学一致映射。** Ti–Co 第4步和第6步的反应自由能差分别为 {fmt(delta1)} eV、{fmt(delta2)} eV，均远大于输入舍入量级。
2. **不能只保留 NEB 前向势垒并声称正逆路径一致。** 前向保持的 TS 映射会使逆向势垒分别为 {n1['mapped_TS_reverse_barrier_eV']} eV、{n2['mapped_TS_reverse_barrier_eV']} eV，与 NEB 逆向值分别差 {n1['mapped_reverse_minus_NEB_reverse_eV']} eV、{n2['mapped_reverse_minus_NEB_reverse_eV']} eV。
3. 当前 CatMAP 第4步的实际有效势垒为 {tico1['current_CatMAP_effective_forward_barrier_eV']}/{tico1['current_CatMAP_effective_reverse_barrier_eV']} eV（前/逆），因显式TS落在独立态FS之下而由FS端点控制；第6步为 {tico2['current_CatMAP_effective_forward_barrier_eV']}/{tico2['current_CatMAP_effective_reverse_barrier_eV']} eV。前者不复现真实 NEB1 双向势垒，后者的热力学驱动力也不同于 NEB2。
4. **差异的数值存在性已确认，差异的唯一结构/物理起因未确认。** 不将其直接命名为共吸附能、slab校正或Ha相互作用能。
5. Ti–Co 的反应能不一致已足以否定“把一套统一转换直接推广到其他13个位点”；14位点表仅归档各位点原始 NEB 与本机独立态能量的现有对照，并逐位点标出差值状态，不套用 Ti–Co 偏移。

## 核心恒等式

任一给定端点的过渡态能量满足：

`ΔG‡f − ΔG‡r = (G_TS−G_IS)−(G_TS−G_FS) = G_FS−G_IS = ΔG_rxn`。

如果NEB与独立态的 `ΔG_rxn` 不同，则不存在一个共同 TS 能量（不论整体加同一常数或从独立 IS 平移）能使前、逆势垒同时等于各自 NEB 势垒。这是代数约束，不是求解器或精度问题。CatMAP 的端点最高态规则只能按当前端点定义有效势垒，不能消除参考不一致，也不能证明真实 NEB路径已被复现。

## 其他位点可推广范围

所有14个位点都能从归档反应势垒表得到 NEB1/NEB2 的 IS/TS/FS 数值；当前 `energies_dft_neb.txt` 亦含用于独立态网络的各物种输入。逐位点 reaction-energy 对照列于 `14_SITE_NEB_MAPPING_STATUS.csv`。这只支持“逐位点计算并报告差异”，不支持统一修正规则。几何参考或端点同一性未证明的位点状态不升级为 `MAPPING_CONSISTENT`。

## 研究者决策前的最小下一步

- 若端点本应是同一微观状态：提供或确认可追溯的共吸附/独立态一一对应及共同 slab/参考映射，然后重新核对端点能量（不动现有原始值）。
- 若端点确为不同状态：需要连接状态的自由能以及相应真实连接势垒/路径；不能用反应能差冒充势垒，也不能设零势垒。
- 共吸附态占据数仅在有几何依据后定义；禁止以化学片段数作为 `n_sites`。

## 执行范围和状态

- 项目预检：`MODEL_CONTRACT=PASS`，CatMAP 0.3.1，原始八步模型，9个稳定吸附态，温度333.15 K。
- 本轮执行：只读源表重算、Decimal恒等式和映射可行性审计。
- 本轮未执行：VASP/DFT、NEB、CatMAP稳态、TOF、DRC或模型修改。
- 结果状态：能量算术对照 `CONFIRMED`；成因/结构映射 `UNVERIFIED`；正式物理能量映射 `REFERENCE_MISMATCH`。
- 所有结果仍是既有数据的审核整理，不构成新的 CatMAP 动力学预测。
"""
    (out / "FINAL_GIBBS_MAPPING_CONCLUSION.md").write_text(report_final, encoding="utf-8")

    print(f"AUDIT_OUTPUT={out}")
    print(f"TICO_STEP4_REF_DIFF={fmt(delta1)} eV; mapped reverse error={n1['mapped_reverse_minus_NEB_reverse_eV']} eV")
    print(f"TICO_STEP6_REF_DIFF={fmt(delta2)} eV; mapped reverse error={n2['mapped_reverse_minus_NEB_reverse_eV']} eV")
    print(f"SITE_ROWS={len(site_rows)}; MAPPING_ROWS={len(mapping_rows)}")
    print("CATMAP_STEADY_STATE=NOT_RUN; DFT_NEB=NOT_RUN")


if __name__ == "__main__":
    main()
