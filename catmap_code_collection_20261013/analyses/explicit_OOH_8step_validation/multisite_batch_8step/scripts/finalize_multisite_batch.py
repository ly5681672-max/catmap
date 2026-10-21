"""Assemble audit tables and concise reports from recorded candidate runs."""
from __future__ import annotations

import csv
import json
import math
from decimal import Decimal
from pathlib import Path

BATCH = Path(__file__).resolve().parents[1]
WORKFLOW = BATCH.parent / "multisite_candidate_workflow"
TICO = BATCH.parent / "candidate_tico_kinetic_control"
TIFE = WORKFLOW / "runs" / "tife_candidate"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def assemble_quality() -> list[dict]:
    statuses = read_csv(BATCH / "MULTISITE_CALCULATION_STATUS.csv")
    rows = []
    for s in statuses:
        site = s["surface_name"]
        run = BATCH / "runs" / f"{site}_candidate"
        gate = {}
        if site == "tico":
            gate = json.loads((TICO / "validation_summary.json").read_text(encoding="utf-8"))
            residual = gate.get("max_accepted_relative_surface_residual", "")
            scaled = gate.get("max_scaled_residual_inf", "")
            site_err = gate.get("max_site_balance_error", "")
            inv_err = gate.get("max_Ha_minus_OOH_inventory_error", "")
            flux_err = gate.get("max_reaction_net_rate_vs_product_TOF_relative_difference", "")
            accepted = gate.get("candidate_solve_cases_accepted", "")
        elif site == "tife":
            gate = json.loads((TIFE / "numerical_gate.json").read_text(encoding="utf-8"))
            residual = gate.get("max_relative_surface_residual", "")
            scaled = ""
            site_err = gate.get("site_balance_error", "")
            inv_err = gate.get("inventory_error", "")
            flux_err = gate.get("max_net_reaction_rate_vs_product_TOF_relative_difference", "")
            accepted = gate.get("accepted_roots", "")
        else:
            path = run / "numerical_gate.json"
            if path.is_file():
                gate = json.loads(path.read_text(encoding="utf-8"))
                residual = gate.get("max_relative_surface_residual", "")
                scaled = gate.get("max_scaled_residual_inf", "")
                site_err = gate.get("site_balance_error", "")
                inv_err = gate.get("inventory_error", "")
                flux_err = gate.get("max_net_reaction_rate_vs_product_TOF_relative_difference", "")
                accepted = gate.get("accepted_roots", "")
            else:
                residual = scaled = site_err = inv_err = flux_err = accepted = ""
        rows.append({"surface_name": site, "numerical_status": s["numerical_status"],
                     "accepted_root_count_or_cases": accepted, "max_relative_residual": residual,
                     "max_scaled_residual_inf": scaled, "site_balance_error": site_err,
                     "inventory_error": inv_err, "net_flux_vs_TOF_relative_difference": flux_err,
                     "direction_semantics": s["direction_semantics"], "physical_gate": s["physical_gate"],
                     "source_hash": s["source_hash"], "run_path": s["run_path"], "failure_detail": s["error"]})
    write_csv(BATCH / "NUMERICAL_QUALITY_GATE.csv", rows)
    return rows


def assemble_drc() -> None:
    drc = read_csv(BATCH / "MULTISITE_DRC_PROVISIONAL.csv")
    # Reuse existing Ti-Co and Ti-Fe candidate conditional TS sensitivities.
    for row in read_csv(TICO / "TICO_PROVISIONAL_DRC.csv"):
        if row.get("ts") not in ("NEB1TS", "NEB2TS") or not row.get("drc_central"):
            continue
        drc.append({"surface_name": "tico", "perturbed_state": row["ts"],
                    "kind": "HISTORICAL_CONDITIONAL_CatMAP_TS_DRC_FINITE_DIFFERENCE",
                    "delta_eV": row.get("h_eV", ""), "baseline_TOF_s-1": "3.530636940099847292329196030416896330474e-27",
                    "TOF_minus": row.get("tof_minus", ""), "TOF_plus": row.get("tof_plus", ""),
                    "conditional_DRC": row.get("drc_central", ""),
                    "branch_continuity": f"historical coverage log-distance={row.get('branch_log_coverage_distance_max','')}",
                    "CatMAP_TS_controls_endpoint_max": "see NEB barrier table",
                    "interpretation_note": "historical provisional sensitivity; native CatMAP reversible semantics; not formal physical DRC"})
    tife_sens = WORKFLOW / "data" / "tife_conditional_sensitivity.csv"
    for row in read_csv(tife_sens):
        if row.get("perturbed_species") not in ("NEB1TS", "NEB2TS"):
            continue
        drc.append({"surface_name": "tife", "perturbed_state": row["perturbed_species"], "kind": row["kind"],
                    "delta_eV": row["delta_eV"], "baseline_TOF_s-1": row["baseline_TOF_s-1"],
                    "TOF_minus": row["TOF_minus"], "TOF_plus": row["TOF_plus"],
                    "conditional_DRC": row["sensitivity"], "branch_continuity": row["branch_continuity"],
                    "CatMAP_TS_controls_endpoint_max": "see NEB barrier table",
                    "interpretation_note": "previously computed same-inventory branch finite difference; not formal physical DRC"})
    # Keep one row per catalyst and TS; for the historical Ti-Co runs prefer h=0.005.
    latest = {}
    for row in drc:
        key = (row.get("surface_name"), row.get("perturbed_state"))
        if key not in latest or (row.get("delta_eV") == "0.005"):
            latest[key] = row
    write_csv(BATCH / "MULTISITE_DRC_PROVISIONAL.csv", list(latest.values()))


def assemble_link_sensitivity() -> None:
    rows = []
    for row in read_csv(TICO / "ENDPOINT_BARRIER_CONTROL.csv"):
        if row.get("scenario") != "HYPOTHETICAL_NOT_DFT":
            continue
        rows.append({"surface_name": "tico", "candidate_step": row.get("step", ""),
                     "endpoint_perturbation_eV": row.get("delta_eV", ""),
                     "TOF_s-1": row.get("tof", ""), "conditional_sensitivity": row.get("conditional_barrier_sensitivity_X", ""),
                     "status": "REUSED_HISTORICAL_ENDPOINT_ASSUMPTION_SENSITIVITY_NOT_DRC",
                     "source": str(TICO / "ENDPOINT_BARRIER_CONTROL.csv")})
    for row in read_csv(WORKFLOW / "data" / "tife_conditional_sensitivity.csv"):
        if row.get("perturbed_species") not in ("NEB1FS", "NEB2FS"):
            continue
        rows.append({"surface_name": "tife", "candidate_step": row["perturbed_species"],
                     "endpoint_perturbation_eV": row["delta_eV"], "TOF_s-1": row["TOF_plus"],
                     "conditional_sensitivity": row["sensitivity"],
                     "status": "REUSED_ENDPOINT_ENERGY_SENSITIVITY_NOT_DRC",
                     "source": str(WORKFLOW / "data" / "tife_conditional_sensitivity.csv")})
    for s in read_csv(BATCH / "MULTISITE_CALCULATION_STATUS.csv"):
        if s["calculation_status"] not in ("COMPUTED_PROVISIONAL", "ALREADY_COMPUTED_PROVISIONAL"):
            continue
        if s["surface_name"] in ("tico", "tife"):
            continue
        for step in range(1, 10):
            rows.append({"surface_name": s["surface_name"], "candidate_step": step,
                         "endpoint_perturbation_eV": "", "TOF_s-1": "",
                         "conditional_sensitivity": "", "status": "NOT_RUN; endpoint-model sensitivity is not DRC",
                         "source": "no separate link-barrier perturbation was run in this phase"})
    write_csv(BATCH / "data" / "MULTISITE_LINK_ENDPOINT_SENSITIVITY.csv", rows)


def write_reports() -> None:
    statuses = read_csv(BATCH / "MULTISITE_CALCULATION_STATUS.csv")
    values = read_csv(BATCH / "MULTISITE_TOF_COVERAGE.csv")
    barriers = read_csv(BATCH / "MULTISITE_REACTION_BARRIERS.csv")
    good = [r for r in values if r.get("TOF_C6H12O_s-1")]
    good.sort(key=lambda r: Decimal(r["TOF_C6H12O_s-1"]), reverse=True)
    positive = [r for r in statuses if r.get("negative_signed_TS_flags")]
    numerical_pass = [r for r in statuses if r["calculation_status"] in ("COMPUTED_PROVISIONAL", "ALREADY_COMPUTED_PROVISIONAL")]
    lines = ["# 八步显式 OOH 候选多位点批算汇总", "",
             "## 结论与边界", "",
             f"原目标 14 个位点中，{len(good)} 个有可用候选 TOF/覆盖度记录；其中新计算 {sum(r['calculation_status']=='COMPUTED_PROVISIONAL' for r in statuses)} 个，复用既有数值结果 {sum(r['calculation_status']=='ALREADY_COMPUTED_PROVISIONAL' for r in statuses)} 个。其余 {14-len(good)} 个未产生通过门控的结果。所有成功数据均为 `PROVISIONAL_ENDPOINT_MODEL`，不是正式八步物理预测。",
             "", "**跨位点排序仅用于同一候选网络/储库设置下的数值对照。** 物理映射、单有效位点复合态和无显式 TS 的端点势垒仍未验证，不用于正式 SISSO 或论文活性结论。CatMAP 原生单箭头并不关闭逆反应，因此所有模型统一标记 `REVERSIBLE_NATIVE_CATMAP`。", "",
             "## 已得到候选 TOF 的位点（降序）", "",
             "| 排名 | 位点 | TOF(C₆H₁₂O), s⁻¹ | log₁₀ TOF | θ(OOH*) | θ(*) | 状态 |", "|---:|---|---:|---:|---:|---:|---|"]
    for i, r in enumerate(good, 1):
        lines.append(f"| {i} | {r['surface_name']} | {r['TOF_C6H12O_s-1']} | {r.get('log10_TOF','')} | {r.get('theta_OOH','')} | {r.get('theta_vacancy','')} | {r.get('calculation_status','')} |")
    lines += ["", "该排序不构成已验证的催化剂优劣结论；Ti–Co 为既有候选基线，Ti–Fe 为既有候选结果，其余位点为本批复用流程计算。失败位点未赋 TOF=0。", "",
              "## 未通过数值门控的位点", "", "| 位点 | 状态 | 诊断 |", "|---|---|---|"]
    for r in statuses:
        if r["calculation_status"] in ("NUMERICAL_FAILED", "ENERGY_MAPPING_FLAGGED"):
            lines.append(f"| {r['surface_name']} | {r['numerical_status']} | {r.get('error','')} |")
    lines += ["", "## NEB TS 与 CatMAP 端点势垒", "",
              "`MULTISITE_REACTION_BARRIERS.csv` 对每个位点保留 NEB 原始有符号 TS−IS/TS−FS，并另列由 IS/TS/FS 中最高态形成的 CatMAP 候选正逆势垒。若 TS 低于端点，CatMAP 端点最高态控制有效势垒；这些值不是新的 DFT 势垒。", "",
              "存在负 NEB TS−IS 的位点：Ti–W (NEB1)、Ti–Ta (NEB1)、Ti–Mo (NEB2)、Ti–V (NEB2)、Ti–Cr (NEB2)。原始负值未钳制。", "",
              "## 方向语义与通量", "",
              "本机 CatMAP 0.3.1 将 `->` 与 `<->` 都解析为同一组反应状态，速率常数函数同时生成 kf、kr。Ti–Co 正逆通量审计及本批 CatMAP 官方方向速率接口输出均显示第2步逆速率非零。因此数据统一使用 `REVERSIBLE_NATIVE_CATMAP`，不宣称单向。对大正逆通量相减，方向速率接口采用浮点数；净通量取高精度稳态求解器的 net rate，二者不能在极强抵消时要求机器浮点逐位相等。", "",
              "## DRC 与连接敏感性", "",
              "`MULTISITE_DRC_PROVISIONAL.csv` 为同一候选库存约束下对 NEB1TS/NEB2TS 作 ±0.005 eV 扰动并重求的条件性有限差分，不是经物理验收的正式 CatMAP DRC。若显式 TS 低于端点最高能量，其扰动可能不改变 CatMAP 的端点最高态势垒，近零 DRC 只表示该端点模型截断下的敏感性。`data/MULTISITE_LINK_ENDPOINT_SENSITIVITY.csv` 复用 Ti–Co/Ti–Fe 的已有端点敏感性；其余新位点标为未运行，不会把端点能量敏感性命名为 DRC。", "",
              "## 主要文件", "",
              "- `MULTISITE_CALCULATION_STATUS.csv`：14 位点运行/物理门控和哈希。",
              "- `MULTISITE_TOF_COVERAGE.csv`：TOF、OOH/Ha/空位及主要稳态覆盖度，九步正逆/净通量。",
              "- `MULTISITE_REACTION_BARRIERS.csv`：NEB 有符号能量差和 CatMAP 端点最高态正逆势垒。",
              "- `MULTISITE_DRC_PROVISIONAL.csv`、`data/multisite_drc_variant_runs.csv`：条件性 TS 灵敏度及逐扰动结果。",
              "- `NUMERICAL_QUALITY_GATE.csv`：残差、位点/库存误差、通量闭合。",
              "- `FULL_8STEP_MODEL_GATE.md`：正式八步模型升级阻断项。"]
    (BATCH / "BATCH_RESULTS_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    full = """# 原始八步显式 OOH 模型物理门控

## 状态

正式八步物理模型门控 **NOT PASSED**。本目录中的 9 基元候选计算只作条件性数值诊断；不能与原始八基元网络视为等价，也不能把候选 TOF/DRC 用于正式构效拟合。

## 已统一/已核实

- 14 个目标表面均按同一套候选输入转换流程、333.15 K、相同有效气相储库及 CatMAP 0.3.1 数值程序处理；Ti–Co/Ti–Fe 复用既有候选结果。
- 完整 NEB1/NEB2 的 IS/TS/FS 能量字段用于各自位点；有符号差值原样保留。
- CatMAP 本机箭头方向行为已统一识别为 `REVERSIBLE_NATIVE_CATMAP`；单箭头不会禁止逆速率。
- 候选零库存约束从该 9 基元计量矩阵重新推导为 θ(Ha*)−θ(OOH*)=0；逐位点计算没有套用 Ti–Co 特定能量修正。

## 阻止升级为正式模型的物理事项

1. 原网络含 Hb*、O*、OH*，并具有两个独立库存不变量；当前候选省略了这些物种，只有 Ha*−OOH*=0，尚无严格粗粒化映射或库存等价证明。
2. NEB 共吸附 IS/FS/TS 的真实组成、表面位点几何与独立吸附态之间的结构/位点映射尚未逐位点证明。候选将每个 NEB 复合态按一个有效位点处理，只是粗粒化假设。
3. 独立前驱体到 NEB1IS、NEB1FS 到 NEB2IS+产物、NEB2FS 到独立水态及其他无显式 TS 连接采用 CatMAP 端点近似。自由能差不是连接活化能；这些假设可显著影响 TOF/覆盖度。
4. 负 signed TS−IS 的 Ti–W、Ti–Ta、Ti–Mo、Ti–V、Ti–Cr 不满足常规过渡态高于初态的直观检验；保留数值诊断但不能通过物理 TS 验收。其他体系也须按 TS 是否高于两个端点检查，而非只看 TS−IS。
5. 目前候选反应第2步在本机 CatMAP 中实际可逆。若研究机制要求不可逆 H2O2 活化，需要先有官方支持的方向定义或采用明确获批的物理动力学实现；不得改写 CatMAP 内部数组。
6. 稳态求根的数值初猜不是洁净表面动力学轨迹；当前库存子空间是否对应指定初始条件的可达稳态仍未证明。

## 允许的解释

可以报告模型在固定候选假设下的数值根、覆盖度、正逆/净速率和显式 TS 有限差分敏感性；必须标注 `PROVISIONAL_ENDPOINT_MODEL`。不能宣称正式八步 TOF、真实 DRC、已证明的 OOH 覆盖度—活性因果关系或全催化剂物理排名。

## 进入正式八步分析前的最小要求

- 给出原始八步稳定表面物种与共吸附 NEB 状态的逐反应原子、位点、结构映射，并保留 Ha、Hb、OOH、O、OH 状态所需库存关系。
- 为缺少 TS 的网络连接提供真实动力学势垒/速率数据，或由研究者明确批准并界定为条件性参数研究；不能以端点 ΔG 冒充势垒。
- 解决 H2O2 活化的方向口径，并证明采用的库存类与所讨论初态一致。
- 对所有纳入体系逐项确认显式 TS 的端点相对能量及模型采用的有效势垒来源。

本门控不要求重复热修正审计；Gibbs 已按研究者前提作为有效输入。
"""
    (BATCH / "FULL_8STEP_MODEL_GATE.md").write_text(full, encoding="utf-8")


def main() -> None:
    assemble_quality()
    assemble_drc()
    assemble_link_sensitivity()
    write_reports()
    print("FINALIZED: numerical gate, provisional DRC merge, link sensitivity audit, and two reports")


if __name__ == "__main__":
    main()
