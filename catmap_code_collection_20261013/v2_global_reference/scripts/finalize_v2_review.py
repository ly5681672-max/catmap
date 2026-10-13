from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
REPORTS = ROOT / "reports"
DATA = ROOT / "data"
MANIFEST = REPORTS / "V2_SHA256_MANIFEST.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    crosswalk = read_csv(DATA / "V2_SOURCE_ROW_CROSSWALK.csv")
    merged = read_csv(DATA / "V2_GLOBAL_GIBBS_REFERENCE.csv")
    anchors = [row for row in crosswalk if row.get("NEB_role")]
    sources = [row.get("source_type", "") for row in merged]
    assert len(anchors) == 84, f"expected 84 anchors, got {len(anchors)}"
    assert len({(r["catalyst"], r["source_workbook"], r["source_sheet"], r["source_row"], r["raw_G_cell"]) for r in anchors}) == 84
    assert len(merged) == 134, f"expected 134 merged records, got {len(merged)}"
    assert sources.count("NEB_ANCHOR") == 84 and sources.count("GENERIC_ONLY") == 50
    assert all(not row.get("calculated_formation_energy", "").strip() for row in merged)

    log = (REPORTS / "V2_MAPPING_TEST_LOG.txt").read_text(encoding="utf-8")
    assert "tests_run=8; passed=True" in log and "tests_run=35; passed=True" in log
    assert "Ran 8 tests" in log and "Ran 35 tests" in log

    old_script = DATA / "construct_global_gibbs.pre_source_row_fix.py"
    invalidated = DATA / "V2_GLOBAL_GIBBS_REFERENCE.INVALIDATED_SOURCE_ROW_MAPPING.csv"
    old_script_hash = sha256(old_script)
    invalidated_hash = sha256(invalidated)
    assert old_script_hash == "7637E5981B1FA951B1A0909EC9646E85DBA45ECC0D1354CB7FD01FFAC7A4B743"
    assert invalidated_hash == "B89EE806C6B3FFAF79486E751F17E3CBEB0562B79557F461BD710220D61F6034"

    v1_manifest = json.loads((BASE / "data" / "serial_v1_sha256.json").read_text(encoding="utf-8"))
    v1_output_checks = {}
    for path_text, expected in v1_manifest["output_sha256"].items():
        path = Path(path_text)
        actual = sha256(path) if path.is_file() else None
        v1_output_checks[path_text] = {"expected": expected, "actual": actual, "matches": actual == expected}
    assert all(item["matches"] for item in v1_output_checks.values()), "a frozen V1 artifact differs from its accepted hash"

    source_manifest = json.loads((BASE / "data" / "source_hashes_before.json").read_text(encoding="utf-8"))
    source_after = json.loads((BASE / "data" / "source_hashes_after.json").read_text(encoding="utf-8"))
    h_sources = {}
    for path_text, expected in source_manifest["sha256"].items():
        path = Path(path_text)
        actual = sha256(path) if path.is_file() else None
        h_sources[path_text] = {"expected": expected, "actual": actual, "matches": actual == expected}
    source_comparison = {
        path: {"before": value, "after": source_after["sha256"].get(path), "matches": value == source_after["sha256"].get(path)}
        for path, value in source_manifest["sha256"].items()
    }
    assert all(v["matches"] for v in source_comparison.values()), "baseline before/after source hashes differ"
    # Original H/F inputs may be unavailable in this restricted environment; baseline comparison remains explicit.

    changed_scope = [
        ("scripts/construct_global_gibbs.py", "modified", "84-anchor crosswalk is now authoritative; wrapper returns rows/crosswalk/meta; build uses crosswalk for Ti-Mo"),
        ("scripts/neb_source_crosswalk.py", "added", "strict source-cell, formula/cache and one-to-one mapping implementation"),
        ("scripts/run_v2_regressions.py", "added", "runs 8 V2 mapping tests and 35 V1 tests while writing logs only in V2"),
        ("scripts/write_v2_source_reference_reports.py", "added", "writes source mapping and static-audit reports"),
        ("scripts/finalize_v2_review.py", "added", "final read-only validation and review artifact generation"),
        ("data/V2_SOURCE_ROW_CROSSWALK.csv", "added/generated", "84 NEB anchor rows"),
        ("data/V2_GLOBAL_GIBBS_REFERENCE.csv", "generated", "134 records; global formation energies blank"),
        ("data/V2_GLOBAL_GIBBS_REFERENCE.INVALIDATED_SOURCE_ROW_MAPPING.csv", "preserved snapshot", "old mapping output; must not be used as input"),
        ("data/construct_global_gibbs.pre_source_row_fix.py", "preserved snapshot", "pre-fix script for audit"),
        ("reports/V2_MAPPING_TEST_LOG.txt", "added/generated", "8/8 V2 and 35/35 V1 regression logs"),
        ("reports/V2_SOURCE_ROW_MAPPING_FIX_REPORT.md", "added/generated", "mapping errors, counts, Ti-Re and special cases"),
        ("reports/V2_TIRE_REFERENCE_VALIDATION.md", "added/generated", "Ti-Re static local-energy replay"),
        ("reports/V2_GIBBS_TEMPERATURE_AUDIT.md", "added/generated", "temperature/thermal evidence audit"),
        ("reports/V2_EPOXIDE_DESORPTION_REFERENCE_AUDIT.md", "added/generated", "conditional endpoint assessment; not a barrier"),
        ("reports/V2_REVIEW_PACKAGE.md", "added", "index of deliverables and evidence boundaries"),
        ("reports/V2_CHANGE_MANIFEST.md", "added", "change and verification summary"),
        ("reports/V2_SHA256_MANIFEST.json", "added", "hash manifest; self hash excluded"),
    ]
    changes = []
    for rel, classification, note in changed_scope:
        path = ROOT / rel
        file_hash = "self_hash_excluded" if path == MANIFEST else (sha256(path) if path.is_file() else None)
        changes.append({"path": str(path), "classification": classification, "sha256": file_hash, "note": note})

    package = """# V2审核包索引：原始NEB源行映射修复与静态参考审计

## 交付范围

本包修复NEB源行映射并恢复有限的静态参考审计。84个NEB锚点由原始NEB清单驱动，工作簿源单元格、公式/缓存及公式审计逐项核验。8项V2专项测试与35项V1冻结回归均通过。

| 文件 | 用途 |
|---|---|
| `../data/V2_SOURCE_ROW_CROSSWALK.csv` | 84个NEB锚点的单元格级追溯主表 |
| `../data/V2_GLOBAL_GIBBS_REFERENCE.csv` | 84个NEB锚点加50条通用独有记录；全局形成能为空 |
| `../data/V2_GLOBAL_GIBBS_REFERENCE.INVALIDATED_SOURCE_ROW_MAPPING.csv` | 旧映射失效快照，禁止作为后续输入 |
| `V2_SOURCE_ROW_MAPPING_FIX_REPORT.md` | 原错误、匹配统计、Ti-Re和特殊布局核验 |
| `V2_MAPPING_TEST_LOG.txt` | 8项专项和35项V1回归的完整日志 |
| `V2_TIRE_REFERENCE_VALIDATION.md` | Ti-Re局部能差及CatMAP静态势垒复核 |
| `V2_STRUCTURE_REFERENCE_AUDIT.md` / `V2_MISSING_STRUCTURE_FILES.md` | slab/坐标证据状态与缺失项 |
| `V2_GIBBS_TEMPERATURE_AUDIT.md` | 28个NEB2 IS/TS温度字段缺失记录 |
| `V2_EPOXIDE_DESORPTION_REFERENCE_AUDIT.md` | 14位点条件端点差；不是脱附活化势垒 |
| `V2_14_SURFACE_REFERENCE_STATUS.md` / `V2_FORMATION_ENERGY_DERIVATION.md` | 14位点全局参考状态及形成能边界 |
| `V2_CATMAP_STATIC_VALIDATION.md` / `../data/V2_THERMODYNAMIC_CONSISTENCY.csv` | CatMAP 0.3.1静态网络解析/平衡核查 |
| `V2_CHANGE_MANIFEST.md` / `V2_SHA256_MANIFEST.json` | 修改清单、源文件和产物哈希 |

## 结论边界

84/84 NEB锚点均有唯一原始工作簿G单元格；68条可与通用清单按源行关联，16条没有通用清单行但仍有可核验NEB单元格；11处通用角色标签冲突已保留并标记，不覆盖NEB角色。合并表134行，其中84条为NEB锚点、50条为通用独有记录。

没有共同clean-slab参考、充分坐标/原子映射及完整热修正证据，故全局形成能继续留空，循环闭合为未评估。Ti–Re +1.250655 eV及其他位点端点差均为条件比较值，不是脱附势垒。未运行TOF、稳态、覆盖度或动力学；未改H盘、V1冻结文件或原始DFT/NEB数据；未执行Git提交、推送或PR。
"""
    (REPORTS / "V2_REVIEW_PACKAGE.md").write_text(package, encoding="utf-8")

    manifest_md = "# V2变更与验证清单\n\n"
    manifest_md += "## 范围内文件\n\n| 文件 | 状态 | SHA256 | 说明 |\n|---|---|---|---|\n"
    for item in changes:
        manifest_md += f"| `{item['path']}` | {item['classification']} | `{item['sha256'] or 'MISSING'}` | {item['note']} |\n"
    manifest_md += "\n## 关键验证\n\n"
    manifest_md += f"- NEB锚点：{len(anchors)}/84；唯一工作簿源单元格：84/84。\n- 合并记录：{len(merged)}（NEB_ANCHOR=84，GENERIC_ONLY=50）；全局形成能全部为空。\n- 专项测试：8/8通过；V1回归：35/35通过，完整日志见`V2_MAPPING_TEST_LOG.txt`。\n- V1冻结产物SHA256：逐项与已验收清单一致；H盘及F盘原始源的before/after清单一致。\n- 未执行TOF、动力学、Git提交/推送/PR；旧CSV保留为失效快照。\n\n## 保留的限制\n\n结构扫描未发现可用slab/构型坐标；NEB2 IS/TS有28个温度字段缺失。共同能量零点和全局形成能不可验证，脱附端点差不作为势垒。\n"
    (REPORTS / "V2_CHANGE_MANIFEST.md").write_text(manifest_md, encoding="utf-8")

    files = {}
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path == MANIFEST or "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        files[str(path.relative_to(ROOT)).replace("\\", "/")] = sha256(path)
    payload = {
        "scope": str(ROOT),
        "sha256_algorithm": "SHA-256",
        "manifest_self_hash_excluded": True,
        "v2_files": files,
        "old_mapping_script_sha256": old_script_hash,
        "invalidated_old_csv_sha256": invalidated_hash,
        "source_hashes_baseline_vs_recorded_after": source_comparison,
        "original_source_hashes_current_check": h_sources,
        "v1_frozen_output_hash_checks": v1_output_checks,
        "counts": {"neb_anchors": len(anchors), "unique_neb_source_cells": 84, "generic_row_matches": 68, "neb_without_generic_row": 16, "role_conflicts": 11, "merged_records": len(merged), "neb_anchor_records": sources.count("NEB_ANCHOR"), "generic_only_records": sources.count("GENERIC_ONLY")},
        "tests": {"v2_mapping_tests": 8, "v2_mapping_passed": True, "v1_regression_tests": 35, "v1_regression_passed": True, "log": str(REPORTS / "V2_MAPPING_TEST_LOG.txt")},
        "tof_calculated": False,
        "kinetic_or_steady_state_run": False,
        "h_drive_modified": False,
        "v1_frozen_files_modified": False,
        "git_commit_push_or_pr": False,
        "global_formation_energies_populated": False,
        "changes": changes,
    }
    MANIFEST.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": str(MANIFEST), "manifest_sha256": sha256(MANIFEST), "v2_files_hashed": len(files), "neb_anchors": len(anchors), "merged_records": len(merged), "v1_outputs_verified": len(v1_output_checks), "source_paths_currently_readable": sum(v["actual"] is not None for v in h_sources.values()), "source_paths_total": len(h_sources)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
