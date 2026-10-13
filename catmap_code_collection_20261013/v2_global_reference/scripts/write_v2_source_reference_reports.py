"""Write mapping, Ti-Re reference, temperature, and conditional desorption reports."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


V2_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = V2_ROOT.parent
V1_DATA = PROJECT_ROOT / "data"
DATA = V2_ROOT / "data"
REPORTS = V2_ROOT / "reports"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def write_reports() -> dict:
    REPORTS.mkdir(parents=True, exist_ok=True)
    crosswalk = read_csv(DATA / "V2_SOURCE_ROW_CROSSWALK.csv")
    global_rows = read_csv(DATA / "V2_GLOBAL_GIBBS_REFERENCE.csv")
    surfaces = read_csv(V1_DATA / "SERIAL_14_SURFACE_MAPPING.csv")
    old_snapshot = DATA / "V2_GLOBAL_GIBBS_REFERENCE.INVALIDATED_SOURCE_ROW_MAPPING.csv"
    old_sha = sha256(old_snapshot)
    counts = {
        "anchors": len(crosswalk),
        "unique_cells": len({(r["source_workbook"], r["source_sheet"], r["source_row"], r["raw_G_cell"]) for r in crosswalk}),
        "generic_matches": sum(r["generic_source_row_match_count"] == "1" for r in crosswalk),
        "generic_missing": sum(r["generic_source_row_match_count"] == "0" for r in crosswalk),
        "role_conflicts": sum(r["role_crosswalk_status"] == "ROLE_LABEL_CONFLICT_SOURCE_ROW_RETAINED" for r in crosswalk),
        "label_matches": sum(r["role_crosswalk_status"] == "MATCH_SOURCE_ROW_AND_LABEL" for r in crosswalk),
        "global_rows": len(global_rows),
        "global_energies_blank": all(not r["calculated_formation_energy"] for r in global_rows),
    }
    if counts["anchors"] != 84 or counts["unique_cells"] != 84 or counts["generic_matches"] + counts["generic_missing"] != 84:
        raise ValueError(f"crosswalk counts failed report gate: {counts}")

    tire = [row for row in crosswalk if row["catalyst"] == "tire"]
    map_report = [
        "# V2源行映射修复报告", "",
        "## 旧错误与修复", "",
        "旧版`construct_global_gibbs.py`先按`(catalyst, NEB_role)`建索引，再用通用状态表的`(catalyst, state_role)`查找；通用角色并不总与实际工作簿源行对应，造成漏配或角色错配。修复后以84条`original_neb_state_inventory.csv`记录为NEB主表，使用`catalyst + source_workbook + source_sheet + source_row`作为唯一主键，并进一步验证原始G单元格行号、工作簿公式/缓存及`excel_gibbs_formula_audit.csv`的一对一记录。通用状态表仅作补充元数据，不决定NEB角色。没有近似能量或邻近行回退。", "",
        "旧输出已保留为`../data/V2_GLOBAL_GIBBS_REFERENCE.INVALIDATED_SOURCE_ROW_MAPPING.csv`，SHA256：`" + old_sha + "`；不得作为后续计算输入。", "",
        "## 匹配统计", "",
        "| 项目 | 数量 | 判定 |", "|---|---:|---|",
        f"| NEB锚点 | {counts['anchors']} | 全部进入交叉表 |",
        f"| 唯一原始工作簿单元格 | {counts['unique_cells']} | 一对一 |",
        f"| 与通用清单精确源行关联 | {counts['generic_matches']} | 补充信息关联 |",
        f"| 无通用清单源行 | {counts['generic_missing']} | 原始NEB单元格仍已核验并保留 |",
        f"| 通用/NEB角色冲突 | {counts['role_conflicts']} | `ROLE_LABEL_CONFLICT_SOURCE_ROW_RETAINED`，不覆盖NEB角色 |",
        f"| 角色与标签一致 | {counts['label_matches']} | 按精确源行匹配 |",
        f"| 合并状态表记录 | {counts['global_rows']} | 84个NEB锚点 + 通用独有行；来源类型分列 |",
        "", "## Ti–Re六个原始锚点", "",
        "| NEB角色 | 源行 | 原始G单元格 | 原始G公式 | 原始G缓存/eV | 局部比较G/eV | TS−IS/eV | TS−FS/eV | 温度/K |",
        "|---|---:|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in tire:
        map_energy = row["mapped_neb_local_energy_eV"]
        map_ts_is = row["mapped_TS_minus_IS_eV"]
        map_ts_fs = row["mapped_TS_minus_FS_eV"]
        map_report.append(
            f"| {row['NEB_role']} | {row['source_row']} | {row['raw_G_cell']} | `{row['raw_G_formula']}` | {row['raw_G_cached']} | {map_energy} | {map_ts_is} | {map_ts_fs} | {row['temperature_K'] or '未知'} |"
        )
    map_report.extend([
        "", "Ti–Re NEB1：TS−IS=0.444846 eV、TS−FS=2.372278 eV；NEB2：TS−IS=0.217873 eV、TS−FS=1.030708 eV。映射值与`DFT_TS_relative.csv`一致。", "",
        "## 特殊布局与旁观气相处理", "",
        "- Ti–Mo NEB1 FS：原始AR20=-887.832689 eV；比较值AR21=`AR20+F60`=-988.870460 eV。两列分别记录，不回写原始G。",
        "- Ti–Mo NEB2 FS：原始AR32=-989.164964 eV；按既有V1映射减去F60后的局部比较值=-888.127193 eV。F60公式/缓存及处理说明保留在交叉表。",
        "- Ti–Ta NEB2 TS：使用该催化剂实际列CI28，E0为CD28、G修正为CH28；不按其他位点列号推算。",
        "- Ti–Mn NEB2 FS：工作簿原状态编号文本为`总`，但NEB身份由源行32/单元格N32及NEB清单角色确认，不把编号文字强制转为数值。",
        "", "## 参考能边界", "",
        f"合并表的全局formation_energy空值检查：`{counts['global_energies_blank']}`。没有clean-slab能量、经坐标确认的原子/位点映射及完整共同热修正，因而不生成全局形成能。条件端点差不作为活化势垒。", "",
        "来源：原始工作簿、`DFT_TS_relative.csv`、V1冻结NEB/Gibbs清单、公式审计表及V1路径映射。逐行追溯见`../data/V2_SOURCE_ROW_CROSSWALK.csv`。",
    ])
    (REPORTS / "V2_SOURCE_ROW_MAPPING_FIX_REPORT.md").write_text("\n".join(map_report) + "\n", encoding="utf-8")

    tire_lines = [
        "# Ti–Re NEB全局参考与源行复核", "",
        "此报告区分NEB路径内部可复现量与跨路径全局参考。下表逐项来自工作簿原始G单元格，并由公式审计表和V1映射交叉核验。", "",
        "| 路径 | 状态 | 源行/单元格 | 原始公式 | 原始G/eV | 局部映射G/eV | T/K |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for row in tire:
        tire_lines.append(
            f"| {row['path']} | {row['state']} | {row['source_row']} / {row['raw_G_cell']} | `{row['raw_G_formula']}` | {row['raw_G_cached']} | {row['mapped_neb_local_energy_eV']} | {row['temperature_K'] or '缺失'} |"
        )
    tire_lines.extend([
        "", "NEB1局部能差：TS−IS=0.444846 eV，TS−FS=2.372278 eV，FS−IS=−1.927432 eV。NEB2局部能差：TS−IS=0.217873 eV，TS−FS=1.030708 eV，FS−IS=−0.812835 eV。CatMAP 0.3.1 `get_rxn_energy`静态正/逆势垒复核见`V2_TIRE_ENERGY_CYCLE.md`。", "",
        "NEB1 FS到NEB2 IS的环氧己烷脱附端点比较仍为条件值`+1.250655 eV`（依赖跨NEB共同slab和组成参考），不是脱附势垒。NEB1与NEB2的绝对能量零点不同；目前没有clean-slab及原子映射证据，因此全局formation_energy留空、跨路径闭合为UNRESOLVED/NOT_EVALUABLE。",
        "", "NEB1三态温度记录为333.15 K；NEB2 IS与TS温度单元格缺失，FS记录333.15 K。温度和热修正缺口不由局部能差补齐。未执行TOF、稳态或动力学计算。",
    ])
    (REPORTS / "V2_TIRE_REFERENCE_VALIDATION.md").write_text("\n".join(tire_lines) + "\n", encoding="utf-8")

    neb_rows = read_csv(V1_DATA / "original_neb_state_inventory.csv")
    temps = {"333.15": 0, "missing": 0, "other": 0}
    for row in neb_rows:
        temps["333.15" if row["temperature_K"] == "333.15" else "missing" if not row["temperature_K"] else "other"] += 1
    missing = [row for row in neb_rows if not row["temperature_K"]]
    temp_lines = [
        "# Gibbs温度与热修正字段审计", "",
        f"84个NEB状态中，温度字段为333.15 K的有{temps['333.15']}个，缺失{temps['missing']}个，其他温度{temps['other']}个。缺失项全部来自NEB2 IS/TS（14个位点×2），不得自动填成333.15 K。", "",
        "| 催化剂 | NEB角色 | 源行 | G单元格 | E0单元格 | Gibbs修正单元格 | 温度字段 |",
        "|---|---|---:|---|---|---|---|",
    ]
    for row in missing:
        temp_lines.append(
            f"| {row['catalyst']} | {row['NEB_role']} | {row['source_row_from_DFT_TS_relative_reference_note']} | {row['raw_G_cell']} | {row['E0_cell']} | {row['Gibbs_correction_cell']} | 缺失 |"
        )
    temp_lines.extend([
        "", "Excel中虽有E0与Gibbs修正公式/缓存，不代表缺失的温度元数据已被证实。原频率输出、振动自由能与计算温度设置尚未提供，故不能独立重建热修正。该限制影响NEB2热化学解释，不改变其源表有符号代数差值的复现。",
    ])
    (REPORTS / "V2_GIBBS_TEMPERATURE_AUDIT.md").write_text("\n".join(temp_lines) + "\n", encoding="utf-8")

    surface_lines = [
        "# 环氧己烷脱附条件端点审计", "",
        "以下数值沿用V1 14位点映射表，仅是NEB1终态与NEB2初态加气相环氧己烷后的条件端点差。共同slab、同一表面组成/占位及跨路径热修正尚未验证，因此全部保持CONDITIONAL；这些数值不是脱附活化势垒。", "",
        "| 位点 | 条件脱附ΔG/eV | 活化势垒/eV | 状态 | 全局形成能/eV |",
        "|---|---:|---:|---|---|",
    ]
    for row in surfaces:
        barrier = row["epoxide_desorption_activation_barrier_eV"]
        delta = row["epoxide_desorption_delta_G_eV"]
        surface_lines.append(
            f"| {row['catalyst']} | {delta} | {barrier or '未提供'} | {row['epoxide_desorption_status']} |  |"
        )
    surface_lines.extend([
        "", "Ti–Re条件端点差为+1.250655 eV；不得将其解释为势垒。全局formation_energy列留空。没有额外假设或人工平移。",
    ])
    (REPORTS / "V2_EPOXIDE_DESORPTION_REFERENCE_AUDIT.md").write_text("\n".join(surface_lines) + "\n", encoding="utf-8")
    return {**counts, "generic_matches": counts["generic_matches"], "temperature_counts": temps, "old_csv_invalidated_sha256": old_sha}


if __name__ == "__main__":
    import json
    print(json.dumps(write_reports(), ensure_ascii=False, indent=2))
