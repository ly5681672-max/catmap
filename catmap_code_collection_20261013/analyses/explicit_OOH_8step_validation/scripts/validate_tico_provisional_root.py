"""Validate archived Ti-Co explicit-OOH zero-inventory root without solving CatMAP.

This checks provenance and numerical invariants in the archived JSON only. It
does not upgrade the result to a physically validated eight-step prediction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from decimal import Decimal
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diagnostic-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    diagnostic = args.diagnostic_dir
    root_path = diagnostic / "data" / "tico_zero_inventory_root.json"
    replay_path = diagnostic / "data" / "cross_precision_replay.json"
    manifest_path = diagnostic / "data" / "INPUT_SHA256.txt"
    root = json.loads(root_path.read_text(encoding="utf-8"))
    replay = json.loads(replay_path.read_text(encoding="utf-8"))

    checks: list[dict[str, str]] = []

    def check(name: str, status: str, observed: str, criterion: str, evidence: str) -> None:
        checks.append({
            "check": name,
            "status": status,
            "observed": observed,
            "criterion_or_scope": criterion,
            "evidence": evidence,
        })

    input_hashes: list[tuple[str, str, str]] = []
    for line in manifest_path.read_text(encoding="utf-8-sig").splitlines():
        if not line or line.startswith("SHA256"):
            continue
        expected, source = line.split(maxsplit=1)
        source_path = Path(source.strip())
        actual = sha256(source_path) if source_path.is_file() else "MISSING"
        input_hashes.append((source, expected.upper(), actual))
    scientific_inputs = [item for item in input_hashes if item[0].lower().endswith((".mkm", ".txt", ".csv")) and "2026" not in item[0].lower()]
    hash_ok = all(expected == actual for _, expected, actual in scientific_inputs)
    check("original_model_and_energy_input_hashes", "PASS" if hash_ok else "FAIL",
          f"{sum(expected == actual for _, expected, actual in scientific_inputs)}/{len(scientific_inputs)} matched",
          "SHA256 equals archived diagnostic manifest", "data/INPUT_SHA256.txt")

    check("catmap_version_and_conditions", "PASS" if root["catmap_version"] == "0.3.1" and root["temperature_K"] == "333.15" else "FAIL",
          f"CatMAP {root['catmap_version']}; T={root['temperature_K']} K; reservoirs={root['effective_gas_reservoirs']}",
          "reuse archived original eight-step conditions; no new solve", str(root_path))

    coverages = {key: Decimal(value) for key, value in root["coverage"].items()}
    vacant = Decimal(root["vacant"])
    site_sum = vacant + sum(coverages.values(), Decimal(0))
    nonnegative = all(value >= 0 for value in coverages.values()) and vacant >= 0
    check("coverage_nonnegative", "PASS" if nonnegative else "FAIL", str(nonnegative),
          "all saved stable-species coverages and vacancy >= 0", str(root_path))
    check("site_balance", "PASS" if abs(site_sum - Decimal(1)) < Decimal("1e-300") else "FAIL",
          str(site_sum), "vacancy + nine stable adsorbate coverages = 1", str(root_path))

    inventory_ok = abs(Decimal(root["I1"])) < Decimal("1e-300") and abs(Decimal(root["I2"])) < Decimal("1e-300")
    replay_inventory_ok = abs(Decimal(replay["I1"])) < Decimal("1e-300") and abs(Decimal(replay["I2"])) < Decimal("1e-300")
    check("zero_inventory_class", "PASS" if inventory_ok and replay_inventory_ok else "FAIL",
          f"root I1={root['I1']}, I2={root['I2']}; replay I1={replay['I1']}, I2={replay['I2']}",
          "both archived root and precision replay satisfy I1=I2=0 to stored precision", str(root_path))

    rates = [Decimal(value) for value in root["reaction_net_rates"]]
    tof = Decimal(root["C6H12O_net_tof"])
    tof_relative_difference = abs(rates[0] - tof) / abs(tof) if tof else Decimal("Infinity")
    flux_ok = (len(rates) == 8 and all(rate == rates[0] for rate in rates)
               and tof_relative_difference < Decimal("1e-30")
               and Decimal(root["cycle_max_relative_error"]) < Decimal("1e-60")
               and Decimal(root["stoichiometry_flux_relative_error"]) < Decimal("1e-60"))
    check("eight_step_net_flux_matches_product_tof", "PASS" if flux_ok else "FAIL",
          f"n_rates={len(rates)}, first_rate={rates[0]}, product_TOF={tof}, relative_difference={tof_relative_difference}, cycle_error={root['cycle_max_relative_error']}, stoichiometry_flux_error={root['stoichiometry_flux_relative_error']}",
          "all eight archived net rates equal each other; rate/TOF relative difference <1e-30; archived cycle and stoichiometric flux errors <1e-60", str(root_path))

    numeric_ok = (root.get("solver_converged_and_validated") is True
                  and replay.get("passed_existing_steady_quality_gate") is True
                  and Decimal(root["max_independent_residual_abs"]) < Decimal("1e-90")
                  and Decimal(replay["max_abs_surface_residual_independent"]) < Decimal("1e-50"))
    check("archived_numerical_steady_quality", "PASS" if numeric_ok else "FAIL",
          f"root residual={root['max_independent_residual_abs']}; replay residual={replay['max_abs_surface_residual_independent']}; replay Δlog10TOF={replay['delta_log10_TOF_vs_high_precision']}",
          "recorded existing gates and independent precision replay; no thresholds changed", str(replay_path))

    check("neb_to_independent_state_physical_mapping", "UNVERIFIED",
          "NEB1 ΔG=-0.239686 vs independent-state ΔG=+1.8006475 eV; NEB2 ΔG=-0.480580 vs independent-state ΔG=-2.0620670 eV",
          "must establish state identity or provide real link-state free energies/barriers; differences are not TS corrections",
          "NEB_MAPPING_VALIDATION.md; TICO_COADSORPTION_MAPPING_UPDATE.md")
    check("clean_surface_global_reachability", "UNVERIFIED",
          "root is locally stable within the constrained zero-inventory subspace; no clean-surface time trajectory / basin proof in archived evidence",
          "local stability is not global reachability", str(root_path))
    check("official_catmap_drc_same_inventory_branch", "BLOCKED",
          "official MeanFieldSolver.get_rate_control finite-differences get_turnover_frequency; the latter calls get_rate(... verify_coverages=True), which obtains coverage through the standard SteadyStateSolver.get_coverage path",
          "no evidence that perturbation solves retain the custom seven-coordinate I1=I2=0 constraint and same branch; do not report DRC values",
          "CatMAP 0.3.1 mean_field_solver.py:41-62,139-158; steady_state_solver.py:81")
    check("formal_eight_step_tof_status", "BLOCKED_PROVISIONAL_ONLY",
          f"archived old-input conditional TOF={tof} s^-1 (log10={root['log10_C6H12O_net_tof']})",
          "numerical quality passes, physical NEB mapping and clean-surface reachability do not", str(root_path))

    output = args.output_dir
    data_dir = output / "data"
    reports_dir = output / "reports"
    data_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    csv_path = data_dir / "tico_zero_inventory_provisional_gate.csv"
    report_path = reports_dir / "TICO_ZERO_INVENTORY_PROVISIONAL_REUSE.md"
    write_csv(csv_path, checks)

    report = f"""# Ti–Co 显式 OOH 八步 provisional 根复用与计算门控

日期：2026-10-20  
范围：只复核已归档的 Ti–Co 八步零库存根、其输入哈希、数值不变量及 CatMAP 官方 DRC 调用链；**本次没有重新运行 CatMAP 求解器，也没有重算热修正**。研究者已确认原始相关 Gibbs 输入可直接使用，本报告接受该前提。

## 结论

已复用的旧八步输入下 Ti–Co 根在数值层面通过既有稳态质量门控：九个稳定表面物种与空位非负、位点和为 1、零库存不变量满足、八步净通量与产物 TOF 一致，并由独立高精度重放复现。结果仍只能称为 **provisional 数值根**，不能作为正式物理 TOF：完整共吸附 NEB 与独立吸附态反应端点不一致，物理映射尚未建立；从洁净表面出发的全局可达性也未证明。

因此本轮没有生成或修改候选 `.mkm`，没有正式 TOF/DRC，也没有尝试任意能量平移或连接势垒。CatMAP 官方 DRC 不自动保持此自定义零库存子空间，直接给出的数值无法通过同库存/同分支门控。

## 复用根的关键数值

- CatMAP 版本：{root['catmap_version']}；温度：{root['temperature_K']} K；浓度模式及有效储库沿用原输入。
- 产物 C₆H₁₂O 净 TOF：`{root['C6H12O_net_tof']} s⁻¹`；`log10(TOF)={root['log10_C6H12O_net_tof']}`。
- 覆盖度主导态：烯烃吸附态 `θ={root['coverage']['C6H12_s']}`；显式 OOH 中间体 `θ={root['coverage']['OOH_s']}`；空位 `θ*={root['vacant']}`。OOH 覆盖很低不意味着可从八步网络或动力学中删除。
- `I1={root['I1']}`，`I2={root['I2']}`；位点总和 `{root['site_sum']}`。
    - 八步净速率：均为 `{root['reaction_net_rates'][0]}`（与产物净 TOF 一致）。最大独立表面残差 `{root['max_independent_residual_abs']}`；相对净速率残差 `{root['relative_residual_to_max_reaction_rate']}`。
- 独立精度重放：工作精度 {replay['working_dps']} 位；覆盖度最大差 `{replay['max_absolute_coverage_difference']}`；`Δlog10(TOF)={replay['delta_log10_TOF_vs_high_precision']}`；已有质量门控通过。
- 零库存子空间局部稳定性在已有诊断中通过；此项不证明洁净表面全局可达。

覆盖度、反应速率、残差及其精度重放均直接来自 `data/tico_zero_inventory_root.json` 与 `data/cross_precision_replay.json`，没有重算。

## NEB 与独立态映射：正式物理门控未通过

| 路径 | 完整共吸附 NEB ΔG (eV) | 独立吸附态端点组合 ΔG (eV) | 差值 (eV) | 完整 NEB 局部正/逆能垒 (eV) |
|---|---:|---:|---:|---:|
| NEB1 / 候选氧转移步 | −0.239686 | +1.8006475 | +2.0403335 | +0.556344 / +0.796030 |
| NEB2 / 候选后续氧氢转移步 | −0.480580 | −2.0620670 | −1.5814870 | +0.225994 / +0.706574 |

这些差值只能说明当前两套端点定义不相同，不能当作共吸附修正、TS 修正或连接势垒。NEB1/NEB2 都含 Ha 旁观物；必须保留 Ha 的显式计量，但不得假定它与其他吸附物能量可加和或无相互作用。当前可访问的导出结构与原始计算目录尚未一一对应；结构状态同一性和 CatMAP 位点占据映射仍待验证。

### 最少物理缺口

1. 将完整 NEB1/NEB2 IS、TS、FS 与各自 clean slab 的真实计算结构映射到八步中间体，确定哪些端点与独立吸附态是同一微观状态，并确认包含 Ha 的结构/组成和有效位点占据。
2. 若状态并非相同，补齐独立吸附态网络与 NEB 共吸附态网络之间真实连接步骤的自由能和势垒；至少要覆盖进入 NEB1-IS、NEB1-FS 离开至后续网络、环氧物脱附连接 NEB2-IS、NEB2-FS 接回水生成支路。不得设零势垒或从上述端点差值臆造。

研究者已确认 Gibbs 热修正完整，因此热修正/温度栏不列为本轮缺失项。

## 官方 DRC 与可达性

本机 CatMAP 0.3.1 `MeanFieldSolver.get_rate_control` 对 `get_turnover_frequency` 作有限差分；其速率路径调用 `get_rate(..., verify_coverages=True)`，再由 `SteadyStateSolver.get_coverage` 使用常规覆盖度稳态路径。该调用链没有本零库存诊断外层七维重构的 I₁=I₂=0 约束接口证据。因此，直接计算的官方 DRC 不能保证扰动前后仍在相同库存类及同一分支；本轮不输出 DRC 数值。

已有子空间局部特征值全负只表明受约束子空间内局部稳定。洁净空表面到该根的真实动力学轨迹/吸引域尚未验证，故“可达性”保持 UNVERIFIED。

## 分阶段门控

- **数值复用基线：PASS（provisional）**：输入模型/能量哈希匹配历史清单，数值根与重放满足归档门控。
- **物理网络闭合：UNVERIFIED / 阻断正式计算**：NEB 共吸附路径与独立态端点/连接尚未建立同一可追溯能量网络。
- **洁净表面可达性：UNVERIFIED**。
- **正式稳态 TOF：未通过物理验收；仅保留旧输入条件数值**。
- **DRC：BLOCKED；未运行**。
- **本轮新 CatMAP 求解：未运行**；没有新的模型/DFT/TOF 数值。

逐项机器可读门控见 `../data/tico_zero_inventory_provisional_gate.csv`。原始 `.mkm`、能量表、Excel 和历史诊断均未修改。
"""
    report_path.write_text(report, encoding="utf-8")

    print(f"CSV={csv_path}")
    print(f"REPORT={report_path}")
    print(f"CHECKS={len(checks)}; FAIL={sum(row['status']=='FAIL' for row in checks)}")
    return 0 if all(row["status"] != "FAIL" for row in checks) else 2


if __name__ == "__main__":
    raise SystemExit(main())
