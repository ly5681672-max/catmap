"""Read-only structural reference audit; never opens a structure for writing."""

from __future__ import annotations

import csv
import os
from pathlib import Path


V2_ROOT = Path(__file__).resolve().parents[1]
ROOT = V2_ROOT.parent
V1_DATA = ROOT / "data"
REPORTS = V2_ROOT / "reports"
DATA = V2_ROOT / "data"
SOURCE_ROOT = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation")
STRUCTURE_NAMES = {"POSCAR", "CONTCAR", "OUTCAR", "vasprun.xml", "XDATCAR"}
STRUCTURE_SUFFIXES = {".vasp", ".xyz", ".traj", ".cif", ".freq", ".freqs", ".hess"}
TIRE_ROLES = ("NEB1_IS", "NEB1_TS", "NEB1_FS", "NEB2_IS", "NEB2_TS", "NEB2_FS")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def find_structure_files(roots: list[Path]) -> tuple[list[str], list[str]]:
    found: list[str] = []
    unavailable: list[str] = []
    for root in roots:
        if not root.exists():
            unavailable.append(str(root))
            continue
        for current, dirs, files in os.walk(root, onerror=lambda error: unavailable.append(str(error.filename))):
            # Avoid recursing into the output tree while auditing the C workspace.
            dirs[:] = [name for name in dirs if name not in {".git", "__pycache__", ".pytest_cache"}]
            for name in files:
                path = Path(current) / name
                if name.upper() in STRUCTURE_NAMES or path.suffix.lower() in STRUCTURE_SUFFIXES:
                    found.append(str(path))
    return sorted(set(found)), sorted(set(unavailable))


def audit() -> dict:
    rows = read_csv(V1_DATA / "original_neb_state_inventory.csv")
    tire_rows = [row for row in rows if row["catalyst"] == "tire" and row["NEB_role"] in TIRE_ROLES]
    found, unavailable = find_structure_files([ROOT, SOURCE_ROOT])
    findings: list[dict[str, str]] = []
    for row in tire_rows:
        linux_path = row["calculation_paths"]
        path_available = bool(linux_path) and Path(linux_path).exists()
        findings.append({
            "comparison_item": row["NEB_role"],
            "catalyst": "tire",
            "source_sheet": row["worksheet"],
            "source_row": row["source_row_from_DFT_TS_relative_reference_note"],
            "source_cell": row["raw_G_cell"],
            "state_label": row["state_label"],
            "recorded_calculation_path": linux_path,
            "recorded_path_accessible_here": str(path_available).upper(),
            "matching_structure_files_found": " | ".join(found),
            "slab_atom_count": "",
            "cell_and_surface_termination": "",
            "constraints_and_calculation_settings": "",
            "atom_mapping_to_other_neb": "UNRESOLVED_NO_COORDINATES",
            "same_clean_slab_reference": "UNRESOLVED_NO_CLEAN_SLAB_OR_COORDINATES",
            "composition_evidence": "WORKBOOK_LABEL_OR_RESEARCHER_CONFIRMED_ONLY; NOT_COORDINATE_VERIFIED",
            "effective_site_occupancy": "UNRESOLVED_NO_COORDINATES",
        })
    findings.extend([
        {
            "comparison_item": "NEB1_vs_NEB2_common_slab",
            "catalyst": "tire", "source_sheet": "", "source_row": "", "source_cell": "",
            "state_label": "NEB1 and NEB2 blocks", "recorded_calculation_path": "",
            "recorded_path_accessible_here": "FALSE", "matching_structure_files_found": " | ".join(found),
            "slab_atom_count": "", "cell_and_surface_termination": "",
            "constraints_and_calculation_settings": "", "atom_mapping_to_other_neb": "UNRESOLVED",
            "same_clean_slab_reference": "UNRESOLVED", "composition_evidence": "",
            "effective_site_occupancy": "UNRESOLVED",
        },
        {
            "comparison_item": "clean_slab", "catalyst": "tire", "source_sheet": "", "source_row": "",
            "source_cell": "", "state_label": "clean Ti-Re slab reference", "recorded_calculation_path": "",
            "recorded_path_accessible_here": "FALSE", "matching_structure_files_found": " | ".join(found),
            "slab_atom_count": "", "cell_and_surface_termination": "",
            "constraints_and_calculation_settings": "", "atom_mapping_to_other_neb": "UNRESOLVED",
            "same_clean_slab_reference": "MISSING_REFERENCE_ENERGY_AND_STRUCTURE", "composition_evidence": "",
            "effective_site_occupancy": "UNRESOLVED",
        },
        {
            "comparison_item": "post_NEB2_water_state", "catalyst": "tire", "source_sheet": "吉布斯自由能汇总",
            "source_row": "36/41", "source_cell": "AC36/AC41", "state_label": "water adsorbed/product candidate",
            "recorded_calculation_path": "see original_gibbs_state_inventory.csv", "recorded_path_accessible_here": "FALSE",
            "matching_structure_files_found": " | ".join(found), "slab_atom_count": "",
            "cell_and_surface_termination": "", "constraints_and_calculation_settings": "",
            "atom_mapping_to_other_neb": "UNRESOLVED_NO_COORDINATES", "same_clean_slab_reference": "UNRESOLVED",
            "composition_evidence": "WORKBOOK_LABEL_ONLY", "effective_site_occupancy": "UNRESOLVED",
        },
    ])
    DATA.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    with (DATA / "V2_SLAB_COMPARISON.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(findings[0]))
        writer.writeheader()
        writer.writerows(findings)

    linux_paths = [row["calculation_paths"] for row in tire_rows]
    report = [
        "# V2结构与共同slab参考审计（Ti–Re）", "",
        "## 结论", "",
        "**结构对应与共同slab参考：UNRESOLVED。**研究者确认的结构解释作为输入事实保留，但当前C盘工作区和H盘原始项目中未发现POSCAR、CONTCAR、OUTCAR、vasprun.xml、XDATCAR、VASP/XYZ/TRJ/CIF或振动文件。扫描发现结构文件数为0。Excel记录的Linux计算路径在本机不可访问，路径名本身不能证明结构或slab一致。", "",
        "NEB1终态仍为吸附C6H12O，Ha/Hb/O分散吸附，尚未形成分子态H2O；C6H12O脱附后进入NEB2。此结构解释来自研究者确认，不等价于坐标核验。", "",
        "## Ti–Re Excel记录的六个计算路径", "",
    ]
    report.extend(f"- `{row['NEB_role']}`（源行{row['source_row_from_DFT_TS_relative_reference_note']}，{row['raw_G_cell']}）：`{row['calculation_paths']}`" for row in tire_rows)
    report.extend([
        "", "## 逐项判定", "",
        "| 审计项 | 状态 | 依据 |", "|---|---|---|",
        "| NEB1/NEB2催化剂与slab是否相同 | UNRESOLVED | 无原子坐标、晶胞、slab原子数及洁净slab能量 |",
        "| 晶胞/表面终止/约束/计算设置 | UNRESOLVED | 仅Excel Linux路径文字可用，计算目录未挂载 |",
        "| Ha/Hb/O原子来源和逐原子映射 | UNRESOLVED | 标签及机理确认不能替代原子索引/坐标 |",
        "| NEB1 FS吸附环氧己烷及NEB2 IS脱附状态 | 研究者结构解释已确认；坐标核验未完成 | 缺NEB1 FS、NEB2 IS构型文件 |",
        "| 局部表面重构和有效活性中心占位 | UNRESOLVED | 缺优化前后坐标及位点几何 |",
        "| 单有效中心表示 | CONDITIONAL_MODEL_ASSUMPTION | V1 `n_sites=1`仅为粗粒化建模假设 |",
        "", "缺失路径扫描错误/不可访问目录：" + ("；".join(unavailable) if unavailable else "无扫描异常记录。"),
        "完整逐状态矩阵见`../data/V2_SLAB_COMPARISON.csv`；需提供的文件及原Excel路径见`V2_MISSING_STRUCTURE_FILES.md`。",
    ])
    (REPORTS / "V2_STRUCTURE_REFERENCE_AUDIT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    missing = [
        "# V2待补结构与计算文件清单（Ti–Re优先）", "",
        "当前本机对Excel中以下Linux计算路径均返回不可访问；H盘原始项目与C盘工作区的结构文件扫描结果为0。请从原计算节点导出文件并保留原目录/作业标识。", "",
        "| 结构/计算态 | Excel工作表与源行 | 已记录的Linux计算路径 | 需要提供 |", "|---|---:|---|---|",
    ]
    for row in tire_rows:
        missing.append(f"| {row['NEB_role']} | 吉布斯自由能汇总!{row['source_row_from_DFT_TS_relative_reference_note']} ({row['raw_G_cell']}) | `{row['calculation_paths']}` | POSCAR/CONTCAR或等价坐标、OUTCAR/vasprun.xml、完整输入参数及频率热化学输出 |")
    missing.extend([
        "| NEB1/NEB2共同clean slab | 无源单元格/路径 | 未记录 | 与两条NEB相同超胞、终止面和约束的clean-slab结构、总能与热修正 |",
        "| NEB1 FS→NEB2 IS连接 | 状态源行20→24 | 两条各自路径如上 | 原子编号映射、含/不含C6H12O的配对结构和表面弛豫说明 |",
        "| NEB2后水吸附/水生成态 | 吉布斯自由能汇总!36/41（Ti–Re列AC） | 见原始状态清单 | H2O*与产物态坐标、计算输出、对应slab/振动数据 |",
        "| 独立附件《己烯环氧化(2).xlsx》 | — | 当前F盘下载及附件位置未找到 | 若需版本对照，重新提供附件；不以主工作簿替代副本 |",
        "", "在这些文件到齐并完成逐原子/计算设置核对前，不生成全局CatMAP形成能。",
    ])
    (REPORTS / "V2_MISSING_STRUCTURE_FILES.md").write_text("\n".join(missing) + "\n", encoding="utf-8")
    return {"structure_files_found": found, "unavailable_roots_or_paths": unavailable, "tire_state_rows": len(tire_rows)}


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), ensure_ascii=False, indent=2))
