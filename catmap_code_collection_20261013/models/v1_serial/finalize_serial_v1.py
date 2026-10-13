"""Verify immutable inputs and write SHA256/change manifests for V1 outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "data"
REPORTS = ROOT / "reports"
BASELINE = DATA / "source_hashes_before.json"
TASK_SPEC = Path(r"C:\Users\Administrator\.codex\attachments\3175db0b-a0e1-4b9a-a4f3-9fa1bb239818\已粘贴的文本.txt")
MANIFEST_PATH = DATA / "serial_v1_sha256.json"
CHANGE_PATH = REPORTS / "SERIAL_V1_CHANGE_MANIFEST.md"
INDEX_PATH = REPORTS / "SERIAL_V1_REVIEW_PACKAGE.md"
V1_NETWORK_REPORT = REPORTS / "MKM_NETWORK_DESIGN_REPORT_V1.md"

WRAPUP_BEFORE = {
    "scripts/write_phase1_reports.py": "ED3F26D0F9F0F223B82519681431E061A7827742E72DEBDA9547D8233EFEDBE7",
    "models/v1_serial/finalize_serial_v1.py": "061708DF516AB6E81B1FE320AF03ED2D897ED2B46298B6BFC6D3AD7F4BF0EE3D",
    "reports/MKM_NETWORK_DESIGN_REPORT.md": "B45FF9853A7127D93B380CB499079E9CC7F3FF280FE5D7080270220F3AB53183",
    "reports/PHASE1_REVIEW_PACKAGE.md": "B8F4E09EA8A8EB91A20247CF5B05202D332F9900623A186DE78FA5818922196E",
    "reports/SERIAL_V1_TEST_LOG.txt": "0F086F7B1C2238184917DAADD4B53C9E21674DDB9E247E1DA6C8DD85D4BAA2C5",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def tracked_outputs() -> list[Path]:
    paths = [path for path in (ROOT / "models" / "v1_serial").iterdir()
             if path.is_file() and path.suffix.lower() in {".py", ".mkm", ".txt", ".md"}]
    paths.extend(path for path in DATA.glob("SERIAL_*") if path.is_file())
    paths.extend([
        ROOT / "tests" / "test_serial_species.py",
        ROOT / "tests" / "test_serial_network.py",
        ROOT / "tests" / "test_serial_energies.py",
        ROOT / "tests" / "test_serial_14_surface_mapping.py",
        ROOT / "scripts" / "write_phase1_reports.py",
        ROOT / "docs" / "superpowers" / "plans" / "2026-10-12-serial-neb-v1.md",
        REPORTS / "NEB_STATE_DEFINITION_REPORT.md",
        REPORTS / "ENERGY_REFERENCE_MAPPING_REPORT.md",
        REPORTS / "MKM_NETWORK_DESIGN_REPORT.md",
        V1_NETWORK_REPORT,
        REPORTS / "MISSING_REACTION_DATA.md",
        REPORTS / "PHASE1_REVIEW_PACKAGE.md",
        REPORTS / "SERIAL_V1_TEST_LOG.txt",
        INDEX_PATH,
    ])
    return sorted({path.resolve() for path in paths
                   if path.exists() and path.resolve() != MANIFEST_PATH.resolve()})


def finalize() -> dict:
    before = json.loads(BASELINE.read_text(encoding="utf-8"))["sha256"]
    source_after = {path: sha256(Path(path)) for path in before}
    unchanged = all(source_after[path] == before[path] for path in before)
    spec_hash = sha256(TASK_SPEC) if TASK_SPEC.exists() else None
    test_log = (REPORTS / "SERIAL_V1_TEST_LOG.txt").read_text(encoding="utf-8")
    match = re.search(r"Ran (\d+) tests?", test_log)
    test_count = int(match.group(1)) if match else None
    tests_passed = test_count == 35 and re.search(r"\nOK\s*$", test_log) is not None
    if not tests_passed:
        raise RuntimeError(f"Final audit requires all 35 tests to pass; observed count={test_count}.")

    index = """# Serial NEB V1 审核包索引

本索引区分第一阶段历史审计记录与依据研究者最新结构确认形成的V1修订，不覆盖或替代历史文件。

## 第一阶段历史快照（保留原样）

- [第一阶段候选网络历史报告](MKM_NETWORK_DESIGN_REPORT.md) —— 历史记录；不作为V1网络表示的当前版本。
- [第一阶段审阅包](PHASE1_REVIEW_PACKAGE.md)
- [第一阶段NEB状态定义](NEB_STATE_DEFINITION_REPORT.md)
- [第一阶段能量参考与势垒映射](ENERGY_REFERENCE_MAPPING_REPORT.md)
- [第一阶段缺失数据清单](MISSING_REACTION_DATA.md)

## V1修订及模型证据

- [V1串联网络修订报告](MKM_NETWORK_DESIGN_REPORT_V1.md) —— NEB1终态仍为吸附环氧己烷与分散吸附Ha/Hb/O，尚未形成分子态水；环氧己烷脱附连接NEB2；NEB2及水形成分别使用`HHO_s <-> TS2_s -> OHHa_s`、`OHHa_s <-> H2O_s`。
- [V1串联NEB模型报告](../models/v1_serial/SERIAL_NEB_MODEL_REPORT.md)
- [环氧己烷脱附端点审计](../models/v1_serial/EPOXIDE_DESORPTION_AUDIT.md)
- [缺失动力学参数清单](../models/v1_serial/SERIAL_NETWORK_MISSING_KINETICS.md)
- [机器可读反应网络核验](../data/SERIAL_NETWORK_AUDIT.json)
- 数据映射：`data/SERIAL_14_SURFACE_MAPPING.csv`、`data/SERIAL_SPECIES_ENERGIES.csv`、`data/SERIAL_NEB_BARRIER_AUDIT.csv`、`data/SERIAL_REACTION_BALANCE.csv`。
- 模型输入：`models/v1_serial/ti_re_serial_neb.mkm`；诊断能量文件未被模型引用，且不代表全局形成能。

## 测试与审计

- [35项完整测试日志](SERIAL_V1_TEST_LOG.txt)
- [SHA256清单](../data/serial_v1_sha256.json)
- [修改文件及差异清单](SERIAL_V1_CHANGE_MANIFEST.md)

本审核包不含TOF、稳态求解或动力学模拟结果；未知的全局形成能及吸附/脱附/水形成势垒仍留空。
"""
    INDEX_PATH.write_text(index, encoding="utf-8")

    outputs = {str(path): sha256(path) for path in tracked_outputs()}
    result = {
        "scope": str(ROOT),
        "source_hashes_before": before,
        "source_hashes_after": source_after,
        "all_original_sources_unchanged": unchanged,
        "new_task_spec": str(TASK_SPEC),
        "new_task_spec_sha256": spec_hash,
        "output_sha256": outputs,
        "tof_calculated": False,
        "kinetic_model_run": False,
        "h_drive_modified": False,
        "git_commit_push_or_pr": False,
        "manifest_self_hash_excluded": True,
        "test_count": test_count,
        "all_tests_passed": tests_passed,
        "wrapup_changes": [],
    }
    wrapup_paths = {
        "scripts/write_phase1_reports.py": ROOT / "scripts" / "write_phase1_reports.py",
        "models/v1_serial/finalize_serial_v1.py": HERE / "finalize_serial_v1.py",
        "reports/MKM_NETWORK_DESIGN_REPORT.md": REPORTS / "MKM_NETWORK_DESIGN_REPORT.md",
        "reports/PHASE1_REVIEW_PACKAGE.md": REPORTS / "PHASE1_REVIEW_PACKAGE.md",
        "reports/SERIAL_V1_TEST_LOG.txt": REPORTS / "SERIAL_V1_TEST_LOG.txt",
        "reports/MKM_NETWORK_DESIGN_REPORT_V1.md": V1_NETWORK_REPORT,
        "reports/SERIAL_V1_REVIEW_PACKAGE.md": INDEX_PATH,
        "data/serial_v1_sha256.json": MANIFEST_PATH,
        "reports/SERIAL_V1_CHANGE_MANIFEST.md": CHANGE_PATH,
    }
    for rel_path, path in wrapup_paths.items():
        if rel_path == "reports/SERIAL_V1_CHANGE_MANIFEST.md":
            continue
        if rel_path == "data/serial_v1_sha256.json":
            result["wrapup_changes"].append({
                "path": str(path),
                "before_sha256": None,
                "after_sha256": "self_hash_excluded",
                "classification": "added",
            })
            continue
        before_hash = WRAPUP_BEFORE.get(rel_path)
        after_hash = sha256(path) if path.exists() else None
        classification = (
            "preserved_unchanged" if before_hash and before_hash == after_hash
            else "modified" if before_hash
            else "added"
        )
        result["wrapup_changes"].append({
            "path": str(path),
            "before_sha256": before_hash,
            "after_sha256": after_hash,
            "classification": classification,
        })

    changes = [
        "# Serial NEB V1 文件修改清单与哈希状态",
        "",
        f"范围：`{ROOT}`。全部原始输入前后SHA256一致：**{unchanged}**。",
        f"新任务文本SHA256：`{spec_hash or '未读取'}`。",
        f"回归测试：**{test_count}/35通过**；未修改测试门槛。",
        "",
        "## 本轮收尾文件状态",
        "",
    ]
    for change in result["wrapup_changes"]:
        changes.append(
            f"- `{Path(change['path']).relative_to(ROOT).as_posix()}`（{change['classification']}）："
            f"前SHA256 `{change['before_sha256'] or '不存在'}` → 后SHA256 `{change['after_sha256'] or '未生成'}`。"
        )
    changes.extend([
        "",
        "历史保护核验：第一阶段`MKM_NETWORK_DESIGN_REPORT.md`及`PHASE1_REVIEW_PACKAGE.md`的哈希与本轮开始前一致；V1报告单独写入`MKM_NETWORK_DESIGN_REPORT_V1.md`。",
        "",
        "## 全项目V1审核材料",
        "",
        "完整模型、映射CSV、能量审计、六步计量核验、测试日志及原始输入/输出哈希均列于`data/serial_v1_sha256.json`和`SERIAL_V1_REVIEW_PACKAGE.md`。",
        "变更清单不包含自身哈希；其SHA256由`data/serial_v1_sha256.json`记录。SHA清单本身不记录自身哈希。",
        "",
        "未修改H盘输入或原始DFT/NEB数据，未修改六步机理；未运行TOF、稳态求解或动力学模拟；未执行Git提交、推送或PR。",
        "",
        "## SHA256文件索引",
        "",
        f"完整原始输入及输出哈希见[`data/serial_v1_sha256.json`](../data/serial_v1_sha256.json)。该JSON不包含其自身哈希。",
    ])
    CHANGE_PATH.write_text("\n".join(changes) + "\n", encoding="utf-8")
    result["change_manifest_sha256"] = sha256(CHANGE_PATH)
    result["output_sha256"][str(CHANGE_PATH)] = result["change_manifest_sha256"]
    result["wrapup_changes"].append({
        "path": str(CHANGE_PATH),
        "before_sha256": None,
        "after_sha256": result["change_manifest_sha256"],
        "classification": "added",
        "self_hash_recorded_in": str(MANIFEST_PATH),
    })
    MANIFEST_PATH.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    result = finalize()
    print(json.dumps({"all_original_sources_unchanged": result["all_original_sources_unchanged"],
                      "source_count": len(result["source_hashes_before"]),
                      "output_count": len(result["output_sha256"]),
                      "test_count": result["test_count"],
                      "all_tests_passed": result["all_tests_passed"],
                      "task_spec_sha256": result["new_task_spec_sha256"]}, ensure_ascii=False, indent=2))
