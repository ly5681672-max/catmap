from __future__ import annotations

import csv
import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(r"C:\Users\Administrator\Documents\codex探索\catmap_neb_coadsorption_model")
DATA = ROOT / "data"
REPORTS = ROOT / "reports"
CATALYSTS = ["tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "tico", "titi", "tiw", "tita", "ticr", "ti", "tini"]
SOURCE_FILES = [
    Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\己烯环氧化.xlsx"),
    Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\alkene_epoxidation.mkm"),
    Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\DFT_TS_relative.csv"),
    Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\energies_dft_neb.txt"),
    Path(r"F:\download\CatMAP_%E9%A1%B9%E7%9B%AE%E4%BA%A4%E6%8E%A5%E4%B8%8ECodex%E6%89%A7%E8%A1%8C%E8%AF%B4%E6%98%8E_2026-10-10.docx"),
    Path(r"F:\download\CatMAP_%E9%A1%B9%E7%9B%AE%E4%BA%A4%E6%8E%A5%E4%B8%8ECodex%E6%89%A7%E8%A1%8C%E8%AF%B4%E6%98%8E_2026-10-10.md"),
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def fmt(value, places=6):
    if value in (None, ""):
        return "—"
    return f"{float(value):.{places}f}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def energy_rows_text(mapping: list[dict[str, str]]) -> str:
    lines = [
        "| 催化位点 | NEB路径 | G(IS) | G(TS) | G(FS) | TS−IS | TS−FS | FS−IS | CatMAP正向Ea | CatMAP逆向Ea | 源表交叉核验 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in mapping:
        lines.append(
            f"| {r['catalyst']} | {r['path']} | {fmt(r['G_IS_workbook_or_normalized'])} | "
            f"{fmt(r['G_TS_workbook_or_normalized'])} | {fmt(r['G_FS_workbook_or_normalized'])} | "
            f"{fmt(r['raw_TS_minus_IS'])} | {fmt(r['raw_TS_minus_FS'])} | {fmt(r['raw_FS_minus_IS'])} | "
            f"{fmt(r['CatMAP_effective_Ea_forward'])} | {fmt(r['CatMAP_effective_Ea_reverse'])} | {r['source_energy_crosswalk_status']} |"
        )
    return "\n".join(lines)


def make_reports(v1_network_only: bool = False) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    states = read_csv(DATA / "original_neb_state_inventory.csv")
    generic_states = read_csv(DATA / "original_gibbs_state_inventory.csv")
    mapping = read_csv(DATA / "neb_to_catmap_energy_mapping.csv")
    formulas = read_csv(DATA / "excel_gibbs_formula_audit.csv")
    network = read_csv(DATA / "proposed_reaction_network.csv")
    manifest = json.loads((DATA / "workbook_audit_manifest.json").read_text(encoding="utf-8"))
    api = json.loads((DATA / "catmap_031_api_audit.json").read_text(encoding="utf-8"))
    refs = json.loads((DATA / "energy_mapping_manifest.json").read_text(encoding="utf-8"))

    tire = {r["path"]: r for r in mapping if r["catalyst"] == "tire"}
    tire1, tire2 = tire["NEB1"], tire["NEB2"]
    tire_group = {r["NEB_role"]: r for r in states if r["catalyst"] == "tire"}
    state9 = next(r for r in generic_states if r["catalyst"] == "tire" and r["state_index"] == "9" and r["source_row"] == "36")
    state10 = next(r for r in generic_states if r["catalyst"] == "tire" and r["state_index"] == "10")
    gas = refs["gas_reference_cached_values"]
    g8_plus_epox = float(tire_group["NEB2_FS"]["raw_G_cached"]) + gas["C6H12O_g"]
    g9_plus_epox = float(state9["G_cached_value"]) + gas["C6H12O_g"]
    g10_products = float(state10["G_cached_value"]) + gas["C6H12O_g"] + gas["H2O_g"]
    delta_5_to_6 = float(tire2["G_IS_workbook_or_normalized"]) + gas["C6H12O_g"] - float(tire1["G_FS_workbook_or_normalized"])

    mapping_text = energy_rows_text(mapping)
    missing_temperature = [r for r in states if r["temperature_K"] in (None, "")]
    missing_cache = [r for r in formulas if r["cached_value"] in (None, "")]
    crosswalk_bad = [r for r in mapping if r["source_energy_crosswalk_status"] != "MATCH"]
    negative = [r for r in mapping if float(r["raw_TS_minus_IS"]) < 0]
    hashes = None
    if not v1_network_only:
        hashes = {str(p): sha256(p) for p in SOURCE_FILES if p.exists()}
        (DATA / "source_hashes_before.json").write_text(
            json.dumps({"sha256": hashes, "capture_stage": "phase1_before_report_generation"}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    state_report = f"""# 原始NEB状态定义与Excel来源审计

## 结构解释更新（以2026-10-12研究者确认信息为准）

NEB1_FS的Excel标签`(C6H12O)(H2O)*`仅表示总体组成记账；本项目后续V1机理按研究者确认的完整构型解释为吸附C6H12O及表面Ha、Hb、O，**此时尚未形成分子态H2O**。以下旧标签行只保留源表字面信息，不再作为物种结构解释。

## 范围与证据等级

读取H盘只读源表《己烯环氧化.xlsx》工作表“吉布斯自由能汇总”，并对照`DFT_TS_relative.csv`中注明的NEB源行（12/16/20及24/28/32）。提取器分别保留单元格公式与Excel缓存值；共审计{len(formulas)}个公式单元格，缓存值缺失{len(missing_cache)}个；NEB锚点共{len(states)}行（14个位点×6个IS/TS/FS锚点）。另保留{len(generic_states)}条表内带编号状态记录，避免把不同催化剂工作簿编号误当成统一状态序号。

**文件版本核验限制：**当前可用的是H盘源表，SHA256为`{manifest['primary_sha256']}`。附件目录和F盘下载目录未找到单独的《己烯环氧化(2).xlsx》，因此无法独立比较“此前提供的副本”与H盘版本；本报告不把两者宣称为同一版本。若要做副本级版本对照，请重新提供该xlsx。

**结构证据限制：**本地项目目录中没有NEB原子坐标/图像文件（如POSCAR/CONTCAR、xyz、traj）。工作簿只保留标签、Linux计算路径、E0、修正值、G公式和缓存值。故以下化学式为标签/路径支持的原子总量候选，不等同于已从坐标逐原子计数；活性中心构型、同一slab/晶胞和位点占据数尚未坐标验证。

## 两条NEB的状态定义

| 路径状态 | 工作簿标签/路径语义 | 元素组成候选 | 位点/旁观分子说明 |
|---|---|---|---|
| NEB1_IS（源行12） | `(H2O2)(C6H12)*`；Ti–Re计算路径含`1-hex+ooh+h` | C6H14O2 | 反应物共吸附/OOH+H+己烯的总原子数相符；具体吸附物拆分与位点数未由坐标验证 |
| NEB1_TS（源行16） | NEB1过渡态计算路径；原子由NEB路径守恒推定 | C6H14O2 | TS完整坐标不可用；不是独立CatMAP稳定吸附态 |
| NEB1_FS（源行20） | Excel标签`(C6H12O)(H2O)*`；按研究者确认的原始构型解释为`C6H12O* + Ha* + Hb* + O*` | C6H14O2 | 该标签只作总体组成记账，不表示分子态H2O已形成；结构事实来自研究者确认，坐标仍待归档 |
| NEB2_IS（源行24） | `H_b*+O*+H_a*` | H2O（表面片段总组成） | 标签列出三个单原子吸附片段；坐标/局域位点仍需确认。跨路径对照需核对环氧化物旁观气相参考 |
| NEB2_TS（源行28） | NEB2过渡态计算路径 | H2O（NEB内原子守恒候选） | 真实TS坐标不可用 |
| NEB2_FS（源行32） | `OH*+H_a*` | H2O（表面片段总组成） | 标签列出两个吸附片段；位点数、空位和重排需结构核实 |

逐催化剂锚点、G公式/缓存值、能量路径及单元格见`data/original_neb_state_inventory.csv`和`data/excel_gibbs_formula_audit.csv`。NEB1与NEB2不因共享分子式就自动视为首尾相接。

## 14个位点能量差复现摘要

| 位点 | NEB1 TS−IS | NEB1 TS−FS | NEB2 TS−IS | NEB2 TS−FS |
|---|---:|---:|---:|---:|
""" + "\n".join(
        f"| {c} | {fmt(next(r for r in mapping if r['catalyst']==c and r['path']=='NEB1')['raw_TS_minus_IS'])} | "
        f"{fmt(next(r for r in mapping if r['catalyst']==c and r['path']=='NEB1')['raw_TS_minus_FS'])} | "
        f"{fmt(next(r for r in mapping if r['catalyst']==c and r['path']=='NEB2')['raw_TS_minus_IS'])} | "
        f"{fmt(next(r for r in mapping if r['catalyst']==c and r['path']=='NEB2')['raw_TS_minus_FS'])} |"
        for c in CATALYSTS
    ) + f"""

## 来源映射及特殊表格规则

`DFT_TS_relative.csv`的28组G_IS/G_TS/G_FS均与H盘工作簿指定源行/缓存值或明确旁观分子变换相符（容差1e-6 eV）：交叉核验未匹配数={len(crosswalk_bad)}。Ti–Ta块在G列前额外插入一列备用路径；提取器按E0/G表头动态定位，没有固定列偏移。Ti–Mn有状态编号与能量表头共行的情况；状态编号与CSV所称NEB状态角色不总是同一编号，故本报告优先依据CSV写明的源行，而不把工作簿内部编号强行统一。

Ti–Mo是明确的参考态例外：NEB1终态使用工作簿下一行总G公式`=AR20+F60`，即把环氧化物气相参考加入水吸附态；NEB2终态CSV的`G_FS2_adjusted`则为工作簿`AR32`减去`F60`环氧化物旁观项。两者的公式来源和缓存值分别列于`neb_to_catmap_energy_mapping.csv`。不得将该调整套用到其它位点，也不得直接把未调整的不同组分总G相减。

Excel的“总G”可能是状态G、状态G加旁观气体参考、或用于相对能差的公式；E0、Gibbs修正、G与总G不是可互换字段。完整公式列表有407项且均有缓存；逐项记录见公式审计CSV。
"""
    if not v1_network_only:
        (REPORTS / "NEB_STATE_DEFINITION_REPORT.md").write_text(state_report, encoding="utf-8")

    tire_table = f"""## Ti–Re逐步手算

### NEB1（源行12/16/20）

- G(IS) = `{fmt(tire1['G_IS_workbook_or_normalized'])}` eV（{tire1['IS_source_cell']}）
- G(TS) = `{fmt(tire1['G_TS_workbook_or_normalized'])}` eV（{tire1['TS_source_cell']}）
- G(FS) = `{fmt(tire1['G_FS_workbook_or_normalized'])}` eV（{tire1['FS_source_cell']}）
- ΔG(FS−IS) = G(FS)−G(IS) = **{fmt(tire1['raw_FS_minus_IS'])} eV**
- 正向原始 TS−IS = **{fmt(tire1['raw_TS_minus_IS'])} eV**；逆向原始 TS−FS = **{fmt(tire1['raw_TS_minus_FS'])} eV**。
- CatMAP有效正向Ea = **{fmt(tire1['CatMAP_effective_Ea_forward'])} eV**；有效逆向Ea = **{fmt(tire1['CatMAP_effective_Ea_reverse'])} eV**。

### NEB2（源行24/28/32）

- G(IS) = `{fmt(tire2['G_IS_workbook_or_normalized'])}` eV（{tire2['IS_source_cell']}）
- G(TS) = `{fmt(tire2['G_TS_workbook_or_normalized'])}` eV（{tire2['TS_source_cell']}）
- G(FS) = `{fmt(tire2['G_FS_workbook_or_normalized'])}` eV（{tire2['FS_source_cell']}）
- ΔG(FS−IS) = **{fmt(tire2['raw_FS_minus_IS'])} eV**
- 正向原始 TS−IS = **{fmt(tire2['raw_TS_minus_IS'])} eV**；逆向原始 TS−FS = **{fmt(tire2['raw_TS_minus_FS'])} eV**。
- CatMAP有效正向Ea = **{fmt(tire2['CatMAP_effective_Ea_forward'])} eV**；有效逆向Ea = **{fmt(tire2['CatMAP_effective_Ea_reverse'])} eV**。

四个目标相对能量均在H盘Excel缓存值中重现。对应单元格公式/缓存和完整算式另见`reports/ti_re_manual_energy_validation.xlsx`。

## CatMAP 0.3.1实际势垒语义

在本机`G:\\anaconda\\Lib\\site-packages\\catmap\\model.py`的`ReactionModel.get_rxn_energy()`上直接核验。对显式TS，返回有符号`ΔG=G_FS−G_IS`，并用`Ea_fwd=max(0,G_TS−G_IS,G_FS−G_IS)`；逆向交换IS/FS后相同计算。因此，原始有符号TS差与有效势垒必须分列。负的TS−IS输入保留：当前28条路径中{len(negative)}条TS−IS为负，位点/路径为`{', '.join(r['catalyst']+'-'+r['path'] for r in negative)}`；有效势垒由TS与两端点共同决定，部分正向有效Ea为0并不等同于修改原始DFT值。

NEB状态内部的直接相对能量对共同参考平移不敏感；把完整NEB总G转成CatMAP各伪物种formation_energy则需要同一slab基准、元素/气相reference和旁观组分核算。当前工作簿某些状态使用不同总G/旁观气体公式，且结构/洁净表面跨状态等价尚未验证，故不能给每个复合NEB态强制赋唯一绝对CatMAP形成能。

## NEB1终态→NEB2初态候选连接（仅热力学差，不是势垒）

工作簿Ti–Re：NEB1_FS `G5={fmt(tire1['G_FS_workbook_or_normalized'])}` eV；NEB2_IS `G6={fmt(tire2['G_IS_workbook_or_normalized'])}` eV；环氧化物气相参考`F60={fmt(gas['C6H12O_g'])}` eV。按“NEB2初态表面H/O片段 + 环氧化物气相旁观项”构造候选全组分端点：

`ΔG_connect = G6 + G(C6H12O_g) − G5 = {fmt(delta_5_to_6)} eV`。

按研究者确认的构型，NEB1_FS是C6H12O*与Ha/Hb/O表面片段，NEB2_IS是Ha/Hb/O表面片段；NEB1_FS的Excel标签`(C6H12O)(H2O)*`仅用于总体组成记账，不表示已形成分子态H2O。计入一次C6H12O_g后，两端名义原子数相符。但ΔG仅在slab/超胞、吸附原子集合、旁观气体计数及G修正一致时可比；连接TS/势垒**缺失**，不得令势垒为0。

作为后续状态的工作簿端点提示：NEB2_FS+环氧化物气相参考到水吸附态state9的端点差为`{fmt(g9_plus_epox-g8_plus_epox)}` eV；state9到气相水+环氧化物state10的差为`{fmt(g10_products-g9_plus_epox)}` eV。两段均没有相应TS，不构成可运行完整微观动力学速率的证据。

## 14个位点逐路径能量与CatMAP有效势垒

{mapping_text}

原始`.mkm`与旧`energies_dft_neb.txt`保持只读。旧模型的`OOH-C6H12`/`Hb-O`是从NEB相对信息映射的TS占位值，不是独立计算的真实TS绝对形成能；本报告用完整NEB端点的相对能量做接口语义核对，没有把新复合态能量写回旧输入。
"""
    if not v1_network_only:
        (REPORTS / "ENERGY_REFERENCE_MAPPING_REPORT.md").write_text(
            f"# Gibbs参考态、NEB能量映射与CatMAP有效势垒\n\n"
            f"**机理释义更新：** NEB1_FS按研究者确认构型为C6H12O*+Ha*+Hb*+O*；Excel水分子标签只作总组成记账，不表示已形成H2O。下文端点比较依然是条件值，不能代替共同slab参考核实。\n\n"
            f"数据版本：H盘Excel SHA256 `{manifest['primary_sha256']}`；相对能CSV SHA256 `{refs['input_neb_csv_sha256']}`。\n\n"
            f"{tire_table}", encoding="utf-8")

    network_heading = (
        "# CatMAP共吸附态候选反应网络设计（V1修订版）"
        if v1_network_only else "# CatMAP共吸附态候选反应网络设计"
    )
    v1_preamble = (
        "**版本说明：** 本文件依据研究者最新确认的NEB共吸附结构形成；第一阶段历史报告 `MKM_NETWORK_DESIGN_REPORT.md` 保留不变。\n\n"
        if v1_network_only else ""
    )
    network_report = f"""{network_heading}

{v1_preamble}## 结构解释更新

按研究者2026-10-12确认，NEB1_FS中环氧己烷仍吸附，Ha/Hb/O仍分散吸附，尚未形成分子态H2O；NEB1→NEB2串联候选通过环氧己烷脱附连接，留下Ha/Hb/O。此事实不补足脱附TS、slab参考或位点占据证据。

## 建模结论

**优先假设：NEB1→NEB2串联，但当前不成立为已验证完整微观动力学循环。** NEB1和NEB2各自内部TS相对能量均可在源表复现；连接能量可做有条件的端点比较，但连接结构、原子坐标、位点占据和连接势垒未提供。CatMAP解析通过只能证明语法/计量可以表示，不能证明所写结构就是实际DFT态。

CatMAP 0.3.1将源`.mkm`解析为8个反应、9个稳定吸附物种和2个显式TS。`OOH-C6H12_s`与`Hb-O_s`被解析到transition_state_names，而非稳定adsorbate_names。原八步网络在本机`verify()`元素/位点检查通过。对候选新态，测试证明可通过species `composition`与`n_sites`定义复合态/多位点占位；但NEB1复合态的`n_sites`仍未知。示范候选串联式在暂定位点数下可通过CatMAP检验；去掉连接式中的一个空位后，位点守恒测试失败。该语法测试不构成物理验证。

官方语法要求显式TS写作`IS <-> TS -> FS`，并在初始化检查质量/位点守恒；非化学式物种需显式给出composition。官方 formation-energy 方法还要求同一模型使用一致的参考态。参考：[CatMAP Code Overview](https://catmap.readthedocs.io/en/stable/topics/code_overview.html)、[Generating an Input File / Formation Energy Approach](https://catmap.readthedocs.io/en/latest/tutorials/generating_an_input_file.html)。本文中的具体`n_sites`行为以本机0.3.1 `verify()`和`ParserBase`测试为准。

## 方案A：NEB1为核心的最小网络

| 步骤 | 候选表示 | 证据 | 缺口 |
|---|---|---|---|
| 反应物储库→NEB1_IS | `C6H12_g + H2O2_g + n* -> NEB1_IS_s` | NEB1_IS总元素组成名义为C6H14O2 | 共吸附形成、H2O2活化路径、`n`与真实吸附构型/势垒缺失 |
| NEB1核心 | `NEB1_IS_s <-> NEB1-TS_s -> NEB1_FS_s` | 原始NEB1提供IS/TS/FS相对能 | 每态原子坐标和可用位点数仍需回溯；绝对形成能参考未唯一确定 |
| NEB1_FS→NEB2_IS | `EpoxideHHO_s -> HHO_s + C6H12O_g` | 研究者确认的串联构型；连接端点能仅条件可比 | 环氧化物脱附TS及共同slab参考未验证；水尚未形成 |

方案A不是闭合的可定量速率网络；它可用于保留NEB1路径能量，却不能单独给出可信TOF。

## 方案B：NEB1→NEB2串联扩展网络

| 步骤 | 候选表示 | Ti–Re能量/计量 | 状态 |
|---|---|---|---|
| 进料→NEB1_IS | 同方案A | C6H12+H2O2与NEB1_IS名义组成一致 | 吸附/活化TS缺失 |
| NEB1核心 | 同方案A | 正向/逆向NEB势垒0.444846/2.372278 eV | 直接NEB支持相对能 |
| NEB1_FS→NEB2_IS+epoxide_g | `EpoxideHHO_s -> HHO_s + C6H12O_g` | 条件候选ΔG=**{fmt(delta_5_to_6)} eV**；研究者确认NEB1_FS保留Ha/Hb/O且未形成分子水 | 连接TS/势垒、共用slab与位点数均未确认；绝不设0势垒 |
| NEB2内部重排 | `HHO_s <-> TS2_s -> OHHa_s` | NEB2正/逆势垒0.217873/1.030708 eV；名义H2O原子数守恒 | HHO和OHHa是完整共吸附构型物种；真实位点数需坐标确认 |
| NEB2_FS→H2O* | `OHHa_s <-> H2O_s` | 端点ΔG={fmt(g9_plus_epox-g8_plus_epox)} eV（条件值） | 无连接TS；表面Ha+OH形成分子水的路径/势垒未提供 |
| 水脱附 | `H2O_s -> H2O_g+*` | 端点ΔG={fmt(g10_products-g9_plus_epox)} eV | 无脱附TS/速率数据 |

方案B是优先核查路径，不是当前已验证机理。即使将这些平衡方程通过CatMAP `verify()`，CatMAP对没有显式TS的步会采用端点控制的无额外势垒假设；这不是原始NEB证据，未经批准不能用作正式速率计算。

## 复合态可行性

| NEB态 | composition候选 | CatMAP n_sites | 结论 |
|---|---|---|---|
| NEB1_IS/TS/FS | C6H14O2 | 未知 | 需要坐标确定单一局部反应中心是否等价于一个species，以及占位/可用空位数 |
| NEB2_IS/TS/FS | H2O表面片段；NEB2比较时另有C6H12O气相旁观项 | 标签暗示初态3片段、终态2片段，但坐标未证实 | CatMAP可通过`n_sites`与显式空位表示多位点物种，必须根据结构赋值 |

**结论：CatMAP语法上可建，物理参数化尚不充分。**当前不创建或运行新模型，不计算TOF。
"""
    network_report_path = REPORTS / (
        "MKM_NETWORK_DESIGN_REPORT_V1.md" if v1_network_only else "MKM_NETWORK_DESIGN_REPORT.md"
    )
    network_report_path.write_text(network_report, encoding="utf-8")

    missing_report = f"""# 缺失反应数据与证据清单

**机理澄清：** NEB1_FS按研究者确认结构是C6H12O*+Ha*+Hb*+O*，不是已形成分子水的共吸附态；后续Ha+OH成水归入NEB2之后。

第一阶段确认28/28组NEB状态相对能量可从工作簿/CSV交叉核对，但完整机理定量求解仍受以下数据缺口阻止。

1. **原始结构坐标及NEB图像**：当前项目本地没有POSCAR/CONTCAR/xyz/traj或NEB image目录。需要提供14体系中IS/TS/FS原子结构、超胞、slab原子数、吸附物原子标记、局域Ti–M中心和图像对应关系，以核对实际化学计量及位点数。
2. **NEB1/NEB2同循环连接**：Ti–Re候选态的名义组成可对上，但state5到state6的原子映射、环氧化物是否从表面转为旁观气相、H2O/Ha/Hb/O片段重排、同一slab/同一总原子数及连接TS均无证据。端点候选ΔG={fmt(delta_5_to_6)} eV不能代替连接活化能。
3. **NEB1前置步骤**：从C6H12/H2O2储库形成NEB1_IS（含OOH/H等局域物种）的吸附、解离/氢迁移步骤及各TS/自由能，结构与障碍均缺失。
4. **NEB1后续产品释放**：NEB1_FS为吸附环氧己烷+Ha/Hb/O表面构型，水尚未形成；环氧化物脱附TS缺失，Ha/Hb/O后续水形成和脱附势垒亦缺失。
5. **NEB2后续水生成/脱附**：NEB2_FS到水吸附state9的构型转换TS，以及H2O*脱附TS缺失。Excel endpoint energy只能提供热力学候选差。
6. **CatMAP复合species占位定义**：NEB1完整复合态、NEB2片段态的`n_sites`、位点类型、空位计量与同一催化位点状态映射，需从坐标/表面模型确认。CatMAP代码支持显式composition和n_sites并不意味着这些参数可猜。
7. **状态能量共同参考**：逐状态clean-slab基准、旁观分子/气相加减及总G含义需随结构配对确认。TiMo已证实NEB1_FS和NEB2_FS使用不同旁观分子处理；其他体系须逐个沿公式而非统一常数转换。
8. **附件版本对照**：独立的《己烯环氧化(2).xlsx》副本未在当前附件或F盘下载目录找到；H盘源表可审计，但不能替代副本哈希比对。
9. **自由能温度与热修正元数据**：14个位点的28组NEB锚点行均未提供可机器读取的温度字段。虽然工作表名为“吉布斯自由能汇总”，当前证据不足以确认这些缓存G值采用的温度、熵/零点能修正方案及其与333.15 K研究条件的一致性；不得直接把本轮相对能量称为已验证的333.15 K自由能。

因此当前证据等级是“NEB内部能量差已复现；串联连接是假设；完整CatMAP动力学网络未闭合/未参数化”。
"""
    if not v1_network_only:
        (REPORTS / "MISSING_REACTION_DATA.md").write_text(missing_report, encoding="utf-8")

    review = f"""# 第一阶段审阅包：原始NEB共吸附CatMAP重建

## 结论先行

- **NEB1→NEB2串联是假设，不是既定连接。** Ti–Re从NEB1_FS到“NEB2_IS表面H/O片段+气相环氧化物旁观项”的条件端点差为+{fmt(delta_5_to_6)} eV；连接TS/势垒和坐标映射缺失。
- **Ti–Re NEB1、NEB2原始正逆相对能量均已复现**，CatMAP 0.3.1显式TS有效势垒与这四个正值一致；计算明细在手算xlsx及映射CSV。
- **CatMAP物种/计量可表示，当前不能物理定参。** 0.3.1实际解析区分stable adsorbates与TS placeholders；模型语法允许显式composition和多位点n_sites，但NEB复合态站点占据缺坐标依据。
- **当前不创建第二阶段模型、不计算新TOF、不修改原始数据、不做Git操作。**

## 已交付

- [NEB状态定义报告](NEB_STATE_DEFINITION_REPORT.md)
- [能量参考与有效势垒报告](ENERGY_REFERENCE_MAPPING_REPORT.md)
- [候选网络设计报告](MKM_NETWORK_DESIGN_REPORT.md)
- [缺失反应数据清单](MISSING_REACTION_DATA.md)
- [Ti–Re手工核算表](ti_re_manual_energy_validation.xlsx)
- 数据：`data/original_neb_state_inventory.csv`（84 NEB锚点）、`original_gibbs_state_inventory.csv`、`excel_gibbs_formula_audit.csv`（407公式）、`neb_to_catmap_energy_mapping.csv`（28路径）、`proposed_adsorbate_definitions.csv`、`proposed_reaction_network.csv`、CatMAP API审计JSON及输入哈希清单。

## 用户审核问题

1. 是否认可把NEB1_FS→NEB2_IS作为优先待证串联假设，连接ΔG只作条件端点差，不赋零势垒？
2. 是否能重新提供《己烯环氧化(2).xlsx》附件，完成与H盘源表哈希和公式级版本比对？
3. 是否可提供NEB原子坐标/图像及其计算路径访问方式，以确认原子映射、slab一致性、旁观物种和n_sites？

## 文件版本与限制

H盘源Excel SHA256：`{manifest['primary_sha256']}`；独立附件副本未发现。H盘`.mkm`、`energies_dft_neb.txt`、`DFT_TS_relative.csv`和交接DOCX/Markdown的当前哈希见`data/source_hashes_before.json`。第一阶段输出只在C盘授权目录。

另有28个NEB锚点未标明Gibbs能温度/热修正元数据，因此Ti–Re逐步核算是对源表相对数值及CatMAP势垒语义的复现，不构成“333.15 K自由能已确认”的结论。

阶段测试：Excel提取4项、CatMAP 0.3.1接口3项、能量映射/连接5项、候选网络计量2项（合计14项；全部测试命令和实测最终状态将在最终交付记录中确认）。
"""
    if not v1_network_only:
        (REPORTS / "PHASE1_REVIEW_PACKAGE.md").write_text(review, encoding="utf-8")

    source_hashes_equal = None
    if not v1_network_only:
        after = {str(p): sha256(p) for p in SOURCE_FILES if p.exists()}
        source_hashes_equal = hashes == after
        (DATA / "source_hashes_after.json").write_text(
            json.dumps({"sha256": after, "matches_before": source_hashes_equal}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    print(json.dumps({"reports": [p.name for p in REPORTS.iterdir()], "formula_count": len(formulas),
                      "cache_missing": len(missing_cache), "mapping_rows": len(mapping),
                      "crosswalk_mismatches": len(crosswalk_bad), "missing_temperature_anchor_rows": len(missing_temperature),
                      "source_hashes_equal": source_hashes_equal, "generated_network_report": network_report_path.name,
                      "historical_report_written": not v1_network_only, "tof_calculated": False}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate phase-1 reports or only the isolated V1 network revision.")
    parser.add_argument(
        "--v1-network-only", action="store_true",
        help="write only reports/MKM_NETWORK_DESIGN_REPORT_V1.md; preserve all phase-1 report files and source hash snapshots",
    )
    make_reports(v1_network_only=parser.parse_args().v1_network_only)
