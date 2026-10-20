"""Screen existing validated six-step CatMAP coverages for display-only trace omission.

This script never runs CatMAP and never edits mechanism, energy, rate, TOF, or DRC data.
Run from any directory with Python 3.10+: python scripts/screen_trace_coverage.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path


DEFAULT_PROJECT = Path(__file__).resolve().parents[3]
DEFAULT_EPSILON = Decimal("1e-6")
DEFAULT_DELTA_MAX = Decimal("1e-4")
SITE_SUM_TOL = Decimal("1e-24")


def dec(value: str, field: str) -> Decimal:
    try:
        result = Decimal(value)
    except (InvalidOperation, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid decimal in {field}: {value!r}") from exc
    if not result.is_finite():
        raise ValueError(f"Non-finite value in {field}: {value!r}")
    return result


def decide_omissions(max_coverages, n_sites, epsilon=DEFAULT_EPSILON,
                     delta_max=DEFAULT_DELTA_MAX, quality_ok=True):
    """Return trace candidates and globally omittable adsorbates.

    The conservative cumulative bound sums n_i * max(theta_i) over the whole
    comparison scope. If any row failed quality validation, no omission is allowed.
    """
    if not quality_ok:
        return [], [], Decimal(0), "DATA_QUALITY_FAILED"
    candidates = [name for name, maximum in max_coverages.items()
                  if maximum < epsilon]
    upper_bound = sum((n_sites[name] * max_coverages[name]
                       for name in candidates), Decimal(0))
    if upper_bound >= delta_max:
        return candidates, [], upper_bound, "CUMULATIVE_OCCUPANCY_EXCEEDS_LIMIT"
    return candidates, candidates, upper_bound, "APPROXIMATION_ELIGIBLE" if candidates else "NO_TRACE_SPECIES"


def delta_to_vacancy(delta: Decimal, vacancy: Decimal) -> Decimal:
    """Relative omitted occupancy; infinity warns when vacancy is zero."""
    if vacancy < 0:
        raise ValueError("Vacancy coverage cannot be negative")
    if vacancy == 0:
        return Decimal("Infinity") if delta > 0 else Decimal(0)
    return delta / vacancy


def parse_n_sites(model_path: Path) -> dict[str, Decimal]:
    text = model_path.read_text(encoding="utf-8")
    result = {}
    pattern = re.compile(
        r"^species_definitions\[['\"](?P<name>[^'\"]+)['\"]\]\s*=\s*\{(?P<body>[^\r\n]*)",
        re.MULTILINE,
    )
    for match in pattern.finditer(text):
        name, body = match.group("name"), match.group("body")
        if name == "s" or "'composition'" not in body and '"composition"' not in body:
            continue
        n_match = re.search(r"['\"]n_sites['\"]\s*:\s*([0-9.eE+-]+)", body)
        if not n_match:
            raise ValueError(f"Missing explicit n_sites for adsorbate {name} in {model_path}")
        result[name] = Decimal(n_match.group(1))
    if not result:
        raise ValueError(f"No surface adsorbates with explicit n_sites found in {model_path}")
    return result


def read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_scope(coverage_rows, convergence_rows, balance_rows, n_sites):
    errors = []
    conv = {r["surface"]: r for r in convergence_rows}
    balance = {r["surface"]: r for r in balance_rows}
    if len(conv) != len(convergence_rows) or len(balance) != len(balance_rows):
        errors.append("duplicate surface in quality tables")
    if set(conv) != set(balance):
        errors.append("convergence and site-balance surface sets differ")
    groups = defaultdict(list)
    for row in coverage_rows:
        groups[row["surface"]].append(row)
    if set(groups) != set(conv):
        errors.append("coverage and quality surface sets differ")
    required_species = set(n_sites)
    for surface, rows in groups.items():
        if conv.get(surface, {}).get("status") != "CONVERGED":
            errors.append(f"{surface}: convergence status not CONVERGED")
        b = balance.get(surface, {})
        if b.get("status") != "CONVERGED":
            errors.append(f"{surface}: site-balance status not CONVERGED")
        try:
            residual = dec(conv[surface]["residual"], f"{surface}.residual")
            tolerance = dec(conv[surface]["tolerance"], f"{surface}.tolerance")
            if residual < 0 or residual > tolerance:
                errors.append(f"{surface}: residual exceeds recorded tolerance")
            total = dec(b["total_site_coverage"], f"{surface}.total_site_coverage")
            if abs(total - Decimal(1)) > SITE_SUM_TOL:
                errors.append(f"{surface}: total site coverage outside {SITE_SUM_TOL}")
        except (KeyError, ValueError) as exc:
            errors.append(f"{surface}: malformed quality record ({exc})")
        stable_rows = [r for r in rows if r["is_stable_adsorbate"].lower() == "true"]
        observed = {r["species"] for r in stable_rows}
        if observed != required_species:
            errors.append(f"{surface}: stable adsorbate set differs from model; {observed} vs {required_species}")
        if len(stable_rows) != len(required_species):
            errors.append(f"{surface}: duplicate or missing stable adsorbate rows")
        try:
            theta_by_species = {r["species"]: dec(r["coverage"], f"{surface}.{r['species']}")
                                for r in stable_rows}
            vacancy_rows = [r for r in rows if r["species"] == "s"]
            if len(vacancy_rows) != 1:
                errors.append(f"{surface}: expected one vacancy row")
                continue
            vacancy = dec(vacancy_rows[0]["coverage"], f"{surface}.vacancy")
            all_values = list(theta_by_species.values()) + [vacancy]
            if any(value < 0 or value > 1 for value in all_values):
                errors.append(f"{surface}: coverage outside [0,1]")
            site_sum = vacancy + sum((n_sites[s] * theta_by_species[s] for s in required_species), Decimal(0))
            if abs(site_sum - Decimal(1)) > SITE_SUM_TOL:
                errors.append(f"{surface}: recomputed n_sites-weighted balance error exceeds {SITE_SUM_TOL}")
        except (KeyError, ValueError) as exc:
            errors.append(f"{surface}: invalid coverage value ({exc})")
    return groups, errors


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(project: Path, epsilon=DEFAULT_EPSILON, delta_max=DEFAULT_DELTA_MAX):
    analysis = project / "analyses" / "microkinetic_mechanism_13sites"
    data = analysis / "data"
    model = project / "catmap_official_run" / "alkene_epoxidation.mkm"
    cov_path = data / "coverage_13sites_long.csv"
    conv_path = data / "convergence_13sites.csv"
    balance_path = data / "coverage_site_balance.csv"
    coverage_rows = read_csv(cov_path)
    model_n_sites = parse_n_sites(model)
    stable_from_results = {row["species"] for row in coverage_rows
                           if row["is_stable_adsorbate"].lower() == "true"}
    missing_model_species = stable_from_results - set(model_n_sites)
    if missing_model_species:
        raise ValueError(f"Coverage table has stable adsorbates missing from model n_sites: {missing_model_species}")
    # TS definitions may carry n_sites in the .mkm but are not steady-state coverages.
    n_sites = {name: model_n_sites[name] for name in model_n_sites if name in stable_from_results}
    convergence_rows = read_csv(conv_path)
    balance_rows = read_csv(balance_path)
    groups, quality_errors = validate_scope(coverage_rows, convergence_rows, balance_rows, n_sites)
    quality_ok = not quality_errors

    max_coverages = {}
    min_coverages = {}
    max_occupancies = {}
    for species in n_sites:
        values = [dec(row["coverage"], f"{row['surface']}.{species}")
                  for row in coverage_rows
                  if row["species"] == species and row["is_stable_adsorbate"].lower() == "true"]
        if values:
            min_coverages[species] = min(values)
            max_coverages[species] = max(values)
            max_occupancies[species] = n_sites[species] * max(values)
        else:
            quality_errors.append(f"no coverage values for {species}")
            quality_ok = False
    candidates, omitted, bound, decision = decide_omissions(
        max_coverages, n_sites, epsilon, delta_max, quality_ok
    )
    actual_delta = Decimal(0)
    max_delta_over_vacancy = Decimal(0)
    if omitted:
        for surface, rows in groups.items():
            lookup = {r["species"]: dec(r["coverage"], f"{surface}.{r['species']}") for r in rows}
            current = sum((n_sites[s] * lookup[s] for s in omitted), Decimal(0))
            actual_delta = max(actual_delta, current)
            vacancy = lookup.get("s", Decimal(0))
            ratio = delta_to_vacancy(current, vacancy)
            max_delta_over_vacancy = max(max_delta_over_vacancy, ratio)
    vacancy_rows = [(dec(row["vacancy_coverage"], f"{row['surface']}.vacancy_coverage"), row["surface"])
                    for row in balance_rows]
    min_vacancy, min_vacancy_surface = min(vacancy_rows)

    output = project / "analyses" / "coverage_trace_approximation"
    output.mkdir(parents=True, exist_ok=True)
    screening_rows = []
    for species in n_sites:
        max_theta = max_coverages.get(species)
        max_occ = max_occupancies.get(species)
        species_rows = [r for r in coverage_rows if r["species"] == species and r["is_stable_adsorbate"].lower() == "true"]
        min_row = min(species_rows, key=lambda r: dec(r["coverage"], f"{r['surface']}.{species}")) if species_rows else {}
        max_row = max(species_rows, key=lambda r: dec(r["coverage"], f"{r['surface']}.{species}")) if species_rows else {}
        single_ok = bool(quality_ok and max_theta is not None and max_theta < epsilon)
        screening_rows.append({
            "species_id": species,
            "min_steady_coverage": str(min_coverages.get(species, "")),
            "min_coverage_surface": min_row.get("surface", ""),
            "max_steady_coverage": str(max_theta if max_theta is not None else ""),
            "max_coverage_surface": max_row.get("surface", ""),
            "n_sites": str(n_sites[species]),
            "max_site_occupancy_contribution": str(max_occ if max_occ is not None else ""),
            "single_species_trace_rule": "PASS" if single_ok else ("NOT_EVALUABLE_QUALITY_FAIL" if not quality_ok else "FAIL_THRESHOLD"),
            "omitted_from_approximate_balance": "YES" if species in omitted else "NO",
            "scope": "13 converged six-step strict-baseline catalyst sites; 333.15 K; fixed baseline reservoirs",
            "data_quality": "PASS" if quality_ok else "FAIL: " + " | ".join(quality_errors),
            "screening_decision": decision,
        })
    write_csv(output / "data" / "trace_coverage_screening.csv", list(screening_rows[0]), screening_rows)

    hash_rows = []
    for path in [model, cov_path, conv_path, balance_path]:
        hash_rows.append({"file": str(path), "sha256": sha256(path), "role": "read-only screening input"})
    write_csv(output / "data" / "trace_coverage_input_hashes.csv", ["file", "sha256", "role"], hash_rows)

    species_notes = "\n".join(
        f"- `{name}`（模型态；n_sites={n_sites[name]}）：θ范围 {min_coverages[name]} 至 {max_coverages[name]}；最大位点占据贡献 {max_occupancies[name]}；单项阈值 {'满足' if name in candidates else '不满足'}；全局近似省略 {'是' if name in omitted else '否'}。"
        for name in n_sites
    )
    if omitted:
        approx_lines = [f"近似守恒式：θ_* + Σ_(保留 i) n_i θ_i ≈ 1。", f"省略态：{', '.join(omitted)}。",
                        f"逐工况最大实际累计省略位点占据 δ_max={actual_delta}；保守上界 Σ n_i max(θ_i)={bound}（门槛 {delta_max}）。",
                        f"最大 δ/θ_*={max_delta_over_vacancy}。"]
    else:
        approx_lines = ["本数据范围没有可全局省略的稳定吸附态；不生成删项后的近似式。完整守恒式保持原样。",
                        f"通过单项阈值的候选态：{', '.join(candidates) if candidates else '无'}。",
                        f"候选态保守累计位点占据上界 δ_bound={bound}；本次实际省略 δ=0；最大 δ/θ_*={max_delta_over_vacancy}（无省略项时该比值仅为空集合记账值）。"]
    quality_block = "PASS：全部纳入位点通过收敛、非负/有限覆盖度及加权位点守恒核验。" if quality_ok else "FAIL：筛选质量门控未通过；禁止作微量近似。\n\n" + "\n".join(f"- {e}" for e in quality_errors)
    report = f"""# 微量覆盖度近似筛选规则与结果

## 范围与边界

- 实际模型：CatMAP 六步复合态网络；不是显式 OOH 八步模型。模型来源：`catmap_official_run/alkene_epoxidation.mkm`。
- 数据：既有13个位点稳态覆盖度（Ti–Zr未纳入），333.15 K固定基准储库条件；本脚本不调用 CatMAP、不重算稳态，也不改能量、TOF、速率或DRC。
- 阈值是可配置的分析规则，不是 CatMAP 官方标准或普适物理常数：单态 θ_i < {epsilon}；全部候选态的保守累计位点占据上界 Σ_i n_i max(θ_i) < {delta_max}。
- 本次只筛选13个位点的固定严格基准条件；TS1–TS2描述符网格是能量描述符扫描而非催化剂/操作条件对照，其他浓度扰动只覆盖少数位点，故均不并入本次“跨位点全局可忽略”结论。要对这些范围作统一省略判断，应将每个有效网格/工况覆盖度追加到同一筛选范围后重算最大值。
- 只有收敛状态、残差、覆盖度有限且非负、以及加权位点总和检查均通过的数据才参与判断。

## 数据质量门控

{quality_block}

质量表中求解残差按其逐位点记录容差核验；加权总覆盖度与1的差使用1e-24绝对输出核验阈值，区别于稳态残差容差。

## 逐态筛选

{species_notes}

## 位点守恒式

所有六步模型稳定吸附态均为 `n_sites=1`，因此完整式为：

`θ_* + θ_R + θ_P + θ_Q + θ_W + θ_U = 1`

其中 R、P、Q、W、U 仅为公式简写，分别指模型中的反应物共吸附态、环氧化物仍吸附且 H_a/H_b/O 分散的状态、H_a/H_b/O分散状态、含羟基/氢的后续状态、分子水吸附态。公式符号不改变结构定义；具体程序态 ID 见 CSV。

## 近似判断

{chr(10).join(approx_lines)}

被忽略位点占据相对空位的比值为 δ/θ_*。若空位趋近于零，该比值会很大或无定义；因此即使绝对 δ 小，也必须单独评估。本范围最小空位覆盖度为 {min_vacancy}（{min_vacancy_surface}）。由于当前没有省略物种，δ/θ_* = 0仅为无省略项的记账结果，不代表对非空省略集合的影响评估。

## 动力学边界

覆盖度近似仅用于展示或位点守恒的近似书写。任何低覆盖物种仍保留在 CatMAP 机理、求解方程、反应速率/净通量、TOF 与 DRC 中。本轮未改变这些输出。

## 可复现文件

- 筛选表：`data/trace_coverage_screening.csv`
- 输入哈希：`data/trace_coverage_input_hashes.csv`
- 脚本：`scripts/screen_trace_coverage.py`
- 回归测试：`tests/test_trace_coverage_logic.py`
"""
    (output / "reports" / "TRACE_COVERAGE_RULE.md").parent.mkdir(parents=True, exist_ok=True)
    (output / "reports" / "TRACE_COVERAGE_RULE.md").write_text(report, encoding="utf-8")

    approx_report = f"""# 完整与近似位点守恒表达

本分析使用六步模型的五个稳定吸附态和空位；过渡态不作为稳态覆盖度项。所有稳定吸附态的 CatMAP `n_sites=1`。

## 完整式

`θ_* + θ_R + θ_P + θ_Q + θ_W + θ_U = 1`

此式为 CatMAP 完整位点守恒表达，筛选输出的最大绝对偏差见输入 `coverage_site_balance.csv` 及其哈希。

## 微量筛选后的式

当前全局比较范围中满足单态阈值且可安全整体省略的态：{', '.join(omitted) if omitted else '无'}。

因此本轮**不产生删项近似式**；若把完整式原样写成近似式既没有简化也没有意义。省略累计误差 δ=0（因为没有省略项），不是“观测误差为零”的声明。最小空位覆盖度为 {min_vacancy}（{min_vacancy_surface}）；由于没有省略项，δ/θ_* = 0仅为空集合记账值。

单态阈值={epsilon}；累计位点占据上限={delta_max}；全局判定={decision}；数据门控={'PASS' if quality_ok else 'FAIL'}。

当未来获得满足筛选的数据并实际省略态集合 S 时，正确形式为 `θ_* + Σ_(i∉S) n_i θ_i ≈ 1`，误差项 `δ=Σ_(i∈S) n_i θ_i`；应同时报告 `δ/θ_*`，不得重新归一化保留态。

## 公式中模型态的物理命名

- R：反应物共吸附模型态（`Rpair_s`；保留程序ID用于追溯）。
- P：环氧化物仍吸附且 Ha/Hb/O 分散的模型态（`EpoxideHHO_s`）。
- Q：Ha/Hb/O 分散模型态（`HHO_s`；不是分子水）。
- W：后续含羟基与氢的模型态（`OHHa_s`）。
- U：分子水吸附态（`H2O_s`）。

这些是当前六步模型的态标记，不声称未由结构文件证明的几何细节。

## 影响范围

本筛选不改变任何 CatMAP 输入、求解结果、真实覆盖度、反应通量、TOF 或 DRC；它只判断是否可在展示文字中使用近似式。所有动力学计算必须继续使用完整模型态。
"""
    (output / "reports" / "SITE_BALANCE_APPROXIMATION.md").write_text(approx_report, encoding="utf-8")
    return {"project": project, "output": output, "quality_ok": quality_ok,
            "quality_errors": quality_errors, "decision": decision, "candidates": candidates,
            "omitted": omitted, "delta_bound": bound, "max_coverages": max_coverages,
            "max_occupancies": max_occupancies, "site_count": len(groups)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--epsilon", default=str(DEFAULT_EPSILON))
    parser.add_argument("--delta-max", default=str(DEFAULT_DELTA_MAX))
    args = parser.parse_args()
    result = run(args.project_root.resolve(), Decimal(args.epsilon), Decimal(args.delta_max))
    print(f"quality_ok={result['quality_ok']} sites={result['site_count']} decision={result['decision']}")
    print(f"trace_candidates={result['candidates']} globally_omitted={result['omitted']} delta_bound={result['delta_bound']}")
    for species, maximum in result["max_coverages"].items():
        print(f"{species}: max_theta={maximum}, max_n_theta={result['max_occupancies'][species]}")
    for name in ["TRACE_COVERAGE_RULE.md", "SITE_BALANCE_APPROXIMATION.md"]:
        print(result["output"] / "reports" / name)


if __name__ == "__main__":
    main()
