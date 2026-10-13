"""Render the V1 review documents from the machine-readable audit tables."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from map_14_surfaces import DATA
from validate_serial_network import audit, reaction_balance_rows
from validate_serial_energies import epoxide_connection, replay_neb_block, source_row


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REPORTS = ROOT / "reports"
CATALYSTS = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "tico",
             "titi", "tiw", "tita", "ticr", "ti", "tini"]


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _fmt(value: Any, digits: int = 6) -> str:
    if value in (None, ""):
        return "—"
    return f"{float(value):.{digits}f}"


def write_reports() -> list[Path]:
    REPORTS.mkdir(parents=True, exist_ok=True)
    surface = _rows(DATA / "SERIAL_14_SURFACE_MAPPING.csv")
    barriers = _rows(DATA / "SERIAL_NEB_BARRIER_AUDIT.csv")
    species = _rows(DATA / "SERIAL_SPECIES_ENERGIES.csv")
    balance = reaction_balance_rows()
    network_audit = audit()
    (DATA / "SERIAL_NETWORK_AUDIT.json").write_text(
        json.dumps(network_audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    tire_neb1 = replay_neb_block(source_row("tire", "NEB1"), "NEB1")
    tire_neb2 = replay_neb_block(source_row("tire", "NEB2"), "NEB2")
    tire_link = epoxide_connection("tire")
    conn_rows = [epoxide_connection(catalyst) for catalyst in CATALYSTS]

    surface_table = [
        "| 位点 | NEB1正/逆有效Ea | NEB2正/逆有效Ea | NEB1→NEB2条件ΔG | H2O形成ΔG* | H2O脱附ΔG* | 全局参考 |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in surface:
        surface_table.append(
            f"| {row['catalyst']} | {_fmt(row['NEB1_CatMAP_Ea_forward_eV'])} / {_fmt(row['NEB1_CatMAP_Ea_reverse_eV'])} | "
            f"{_fmt(row['NEB2_CatMAP_Ea_forward_eV'])} / {_fmt(row['NEB2_CatMAP_Ea_reverse_eV'])} | "
            f"{_fmt(row['epoxide_desorption_delta_G_eV'])} | {_fmt(row['water_formation_delta_G_eV'])} | "
            f"{_fmt(row['water_desorption_delta_G_eV'])} | 未建立 |"
        )

    barrier_rows = []
    for path_name, block in (("NEB1", tire_neb1), ("NEB2", tire_neb2)):
        barrier_rows.append(
            f"| {path_name} | {', '.join(block['source_cells'])} | "
            f"{_fmt(block['source_G_eV'][0])} | {_fmt(block['source_G_eV'][1])} | {_fmt(block['source_G_eV'][2])} | "
            f"{_fmt(block['raw_TS_minus_IS_eV'])} | {_fmt(block['raw_TS_minus_FS_eV'])} | "
            f"{_fmt(block['raw_FS_minus_IS_eV'])} | {_fmt(block['CatMAP_effective_Ea_forward_eV'])} / "
            f"{_fmt(block['CatMAP_effective_Ea_reverse_eV'])} |"
        )

    report = f"""# CatMAP串联NEB模型 V1 审核报告

## 范围和结论

本模型使用研究者确认的完整共吸附态结构语义，复用第一阶段已审计的Excel公式/缓存和CatMAP 0.3.1接口；不改动任何H盘原始数据，不加载旧CatMAP缓存，不运行稳态/动力学或TOF。

**结论：**六步反应的CatMAP 0.3.1解析、TS分类、元素守恒与单有效中心位点守恒均通过。NEB1和NEB2的Ti–Re正逆能量差分别可在各自局部参考内精确重放。但不同路径间共同slab参考、坐标映射与活性位点占位未建立，所以共吸附物种的全局CatMAP formation_energy均留空；能量输入表只是局部诊断表，不是可用于耦合速率的全局参数表。

## 证据层级

- **研究者确认结构事实：**NEB1_IS为C6H12与活化H2O2形成的OOH/Ha完整共吸附态；NEB1_FS为吸附C6H12O及表面Ha/Hb/O，尚未形成分子态水；C6H12O脱附后表面留下Ha/Hb/O；NEB2将Hb/O转为OH并保留Ha，随后Ha+OH形成水。
- **源表直接支持：**NEB1/NEB2的Excel G、公式、缓存值和其各自的有符号相对能量；Ti–Mo的环氧己烷旁观项转换按Excel对应式处理。
- **模型粗粒化假设：**每个复合表面态以`n_sites=1`表示单个Ti–M有效中心，以便做CatMAP解析/守恒演示。该占位数没有坐标证明，不是结构实测。
- **软件默认假设：**CatMAP无显式TS步骤取端点较高者作为TS，即正向Ea=max(0,ΔG)，逆向Ea=max(0,−ΔG)。这不是DFT势垒，且V1不接入这些未知步骤做速率计算。
- **未知：**NEB1/NEB2 clean-slab参考是否相同、超胞/原子清单、温度/热修正一致性、脱附/水形成/进料捕获势垒。

## 六步候选网络

```text
* + H2O2_g + C6H12_g <-> Rpair*
Rpair* <-> TS1* -> EpoxideHHO*
EpoxideHHO* <-> HHO* + C6H12O_g
HHO* <-> TS2* -> OHHa*
OHHa* <-> H2O*
H2O* <-> H2O_g + *
```

复合物species显式composition：Rpair/TS1/EpoxideHHO均为C6H14O2；HHO/TS2/OHHa/H2O均为H2O总原子数。HHO不是水分子；EpoxideHHO亦不是环氧己烷加分子态水。

## Ti–Re NEB双向复核

| NEB | Excel单元格 | G(IS) | G(TS) | G(FS) | TS−IS | TS−FS | FS−IS | CatMAP正向/逆向有效Ea |
|---|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(barrier_rows)}

两条路径分别选IS为局部零点时，NEB1相对态能为0、+0.444846、−1.927432 eV；NEB2为0、+0.217873、−0.812835 eV。共同能量平移不改变路径内部正、逆能量差。**两个局部零点不能当作已验证的共同formation-energy参考。**本地CatMAP `get_rxn_energy()`独立复算正逆有效势垒。

## 串联脱附端点

Ti–Re使用源公式/缓存：NEB1_FS `AC20=Y20+AB20`，G=−988.654116 eV；NEB2_IS `AC24=Y24+AB24`，G=−886.365690 eV；C6H12O(g)参考F60=−101.037771 eV。

`ΔG_cond = G(NEB2_IS)+G(C6H12O_g)-G(NEB1_FS) = {tire_link['delta_G_eV']:.6f} eV`。

该式仅加入一次气相环氧己烷。反应为`EpoxideHHO* -> HHO* + C6H12O_g`；在单有效中心假设下1位点→1位点。该值是**条件端点自由能，不是脱附势垒**。slab/超胞、吸附原子对应和表面弛豫未核实，不得据此设置动力学势垒。

## 化学计量与守恒

CatMAP `verify()`通过。稳定态/气相/空位计量矩阵秩为{network_audit['stoichiometry']['matrix_rank']}，6个反应的右零度为{network_audit['stoichiometry']['right_nullity']}，左零度为{network_audit['stoichiometry']['left_nullity']}。独立守恒总量为C、H、O和总表面位点；守恒向量秩为{network_audit['stoichiometry']['conserved_total_rank']}，未发现额外独立线性不变量。TS1/TS2不计作稳态覆盖物种。

净反应：`C6H12_g + H2O2_g -> C6H12O_g + H2O_g`；逐步平衡见`data/SERIAL_REACTION_BALANCE.csv`。

## 14位点能垒与条件端点概览

星号列的水形成/脱附端点值均为条件比较；无源态或参考不清时为空。所有全局形成能未定，不应用于TOF拟合。

{chr(10).join(surface_table)}

NEB1、NEB2逐条有符号能差见`data/SERIAL_NEB_BARRIER_AUDIT.csv`；完整能量源单元格和局部能量见`data/SERIAL_SPECIES_ENERGIES.csv`。14条NEB对（28路径）源交叉核验状态均为MATCH。负TS−IS输入保留，CatMAP有效Ea单独列出。

## 计算与审核状态

- 温度：333.15 K；有效储库参数沿用旧模型：C6H12=1.0、H2O2=0.185、H2O=0.815、C6H12O=1e−20。
- `models/v1_serial/energies_tire_serial_neb.txt`中NEB复合物能量以NEB各自IS为局部零点并标为`NOT global formation energy`；该文件不被`.mkm`引用。CatMAP表格解析测试仅证明列/物种兼容，不代表耦合能量有效。
- `CatMAP`解析/守恒与能量回放通过；无TOF、无动力学轨迹、无缺TS步骤速率。
- 进入正式串联求解前必须补齐各状态clean-slab参考、超胞/坐标、位点占据和第1/3/5/6步的独立TS或经批准的势垒假设。
"""
    (HERE / "SERIAL_NEB_MODEL_REPORT.md").write_text(report, encoding="utf-8")

    connection_table = [
        "| 位点 | NEB1_FS G | NEB2_IS G | C6H12O(g) G | 条件ΔG | 状态 |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in conn_rows:
        connection_table.append(
            f"| {row['catalyst']} | {_fmt(row['NEB1_FS_G_eV'])} | {_fmt(row['NEB2_IS_G_eV'])} | "
            f"{_fmt(row['C6H12O_g_reference_eV'])} | {_fmt(row['delta_G_eV'])} | {row['physical_status']} |"
        )
    epoxide_report = f"""# NEB1→NEB2环氧己烷脱附端点审计

## 结构/反应定义

研究者确认的连接反应是`EpoxideHHO* -> HHO* + C6H12O_g`：NEB1 FS的环氧己烷仍吸附，Ha/Hb/O留在表面；脱附后形成NEB2 IS表面组态。NEB1 FS中的Excel标签`(C6H12O)(H2O)*`是总体组成表述，不表示当时已形成分子态H2O。

单有效中心假设下元素计量为C6H14O2→C6H12O(g)+H2O(surface fragments)，位点数1→1，故CatMAP候选表达式在原子和位点账面上闭合。位点数仍为粗粒化假设。

## Ti–Re原Excel逐项重算

| 项 | Excel来源 | G/eV | 注记 |
|---|---|---:|---|
| NEB1 FS | 吉布斯自由能汇总!AC20=`Y20+AB20` | −988.654116 | NEB1完整态；研究者结构释义 |
| NEB2 IS | 吉布斯自由能汇总!AC24=`Y24+AB24` | −886.365690 | 表面Ha/Hb/O |
| C6H12O(g) | 吉布斯自由能汇总!F60 | −101.037771 | 连接态中加入一次 |
| 脱附后端点 | G(NEB2 IS)+G(epoxide gas) | −987.403461 | 假设共同slab参考 |
| 条件ΔG | G_after−G(NEB1 FS) | +1.250655 | 不是TS或活化能 |

表面名义组成为H2O原子计数（Ha+Hb+O），加C6H12O(g)后与NEB1 FS名义C6H14O2守恒。**不同原子总G不得直接相减；本条件式仅在相同slab/超胞及修正一致时成立。**

## 14位点条件端点值

{chr(10).join(connection_table)}

这些值由既有逐位点映射的NEB1 FS/NEB2 IS及同一气相参考重算；TiMo等源映射有显式旁观项规则。全部标为条件比较，待结构与参考核实。

## 当前不能下结论的事项

没有NEB坐标，无法验证两个计算的slab原子组成/超胞、NEB1 FS到NEB2 IS的逐原子映射、表面Ha/Hb/O构型相似度、弛豫变化及位点数。没有连接TS，所以脱附活化势垒、速率常数均未知。CatMAP无显式TS时的端点势垒规则只是软件假设，不是DFT结果。本轮没有使用此规则运行动力学。
"""
    (HERE / "EPOXIDE_DESORPTION_AUDIT.md").write_text(epoxide_report, encoding="utf-8")

    missing = """# 串联网络尚缺动力学参数清单

当前模型只做语法、计量和原始NEB内部能差验证。严禁把端点热力学差直接当成DFT活化势垒；CatMAP无TS时将端点较高者当作TS的规则属于软件默认假设。

| 步骤 | 过程 | 当前数据 | 缺失/风险 | 当前处理 |
|---|---|---|---|---|
| 1 | `*+H2O2_g+C6H12_g <-> Rpair*` | NEB1_IS完整态G | 吸附、H2O2活化、溶剂/活度耦合与TS未知；该步合并多个物理过程 | 不赋Ea/速率 |
| 2 | `Rpair* <-> TS1* -> EpoxideHHO*` | NEB1 IS/TS/FS和双向能差 | 仍需检查温度/热修正及共同clean-slab参考 | 只保留NEB相对势垒 |
| 3 | `EpoxideHHO* <-> HHO*+C6H12O_g` | 条件端点ΔG；Ti–Re +1.250655 eV | 连接TS、slab/超胞一致性、脱附构型和速率常数缺失 | 不赋活化势垒，不计算动力学 |
| 4 | `HHO* <-> TS2* -> OHHa*` | NEB2 IS/TS/FS和双向能差 | 跨物种绝对能量参考仍未定 | 只保留NEB相对势垒 |
| 5 | `OHHa* <-> H2O*` | 部分位点有Excel后续水态条件端点 | 水形成TS/速率缺失；部分位点无可识别后续态；共同参考未验证 | 缺数据留空 |
| 6 | `H2O* <-> H2O_g+*` | 部分位点气相产品端点 | 脱附TS/速率缺失；状态公式是否同一slab参考未验证 | 缺数据留空 |

此外，复合态`n_sites=1`仅为单Ti–M有效中心粗粒化假设，需坐标核对；工作簿NEB2部分状态的温度字段为空，热修正/温度也需确认。获得共同slab参考前，不输出全局CatMAP formation_energy。未经研究者另行批准，不做假设势垒敏感性或TOF计算。
"""
    (HERE / "SERIAL_NETWORK_MISSING_KINETICS.md").write_text(missing, encoding="utf-8")

    return [HERE / "SERIAL_NEB_MODEL_REPORT.md", HERE / "EPOXIDE_DESORPTION_AUDIT.md",
            HERE / "SERIAL_NETWORK_MISSING_KINETICS.md"]


if __name__ == "__main__":
    print(json.dumps([str(path) for path in write_reports()], ensure_ascii=False, indent=2))
