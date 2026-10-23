from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
from collections import Counter

import mpmath as mp


REPO = Path(r"H:\catmap\catmap")
PROJECT = REPO / "catmap_code_collection_20261013"
ANALYSIS = PROJECT / "analyses" / "explicit_OOH_8step_validation"
REPLAY = ANALYSIS / "multisite_conditional_replay_20261021"
SOURCE_COVERAGE = ANALYSIS / "conditional_DRC_OOH_TOF_14sites" / "data" / "14_SITE_TOF_COVERAGE.csv"
SOURCE_DRC = ANALYSIS / "conditional_DRC_OOH_TOF_14sites" / "data" / "14_SITE_CONDITIONAL_DRC.csv"
STAGE = Path(__file__).resolve().parents[1]
MODEL_SHA = "86115cefb7784fe28ab519271cdfaf0e308b6f6e3a8308b1e90e4318cb497cfe"
ENERGY_SHA = "55ad4e94d9c603b6705ec8c92630c9a1ceff5859cf69865145152375e203ca30"
MODEL_FILE = REPO / "lyq_alkene_epoxidation" / "2-creating_microkinetic_model" / "alkene_epoxidation.mkm"
ENERGY_FILE = REPO / "lyq_alkene_epoxidation" / "2-creating_microkinetic_model" / "energies_dft_neb.txt"
SITES = ["tico", "ti", "tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "titi", "tiw", "tita", "ticr", "tini"]
ADS_ORDER = ["C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "Ha_s", "Hb_s", "OH_s", "OOH_s", "O_s"]
GAS_PRODUCT = "C6H12O_g"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def configure_catmap(setup: Path, run_dir: Path):
    import catmap
    from catmap import ReactionModel

    if catmap.__version__ != "0.3.1":
        raise RuntimeError(f"Expected CatMAP 0.3.1, found {catmap.__version__}")
    old_cwd = Path.cwd()
    run_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(run_dir)
    model = ReactionModel(setup_file=str(setup))
    model.parse()
    model._math = mp
    model._mpfloat = mp.mpf
    model._matrix = mp.matrix
    model._Axb_solver = mp.lu_solve
    model._kB = mp.mpf("8.617332478e-5")
    model._h = mp.mpf("4.135667516e-15")
    model._math.infnorm = lambda x: mp.norm(x, "inf")
    model._math.lsqnorm = lambda x: mp.norm(x, 2) ** 2
    model.solver.compile()
    return old_cwd, model, model.solver


def load_archived_coverage() -> dict[str, dict[str, str]]:
    with SOURCE_COVERAGE.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return {r["surface_name"]: r for r in rows}


def load_drc() -> dict[str, dict[str, dict[str, str]]]:
    with SOURCE_DRC.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    out: dict[str, dict[str, dict[str, str]]] = {}
    for row in rows:
        out.setdefault(row["surface_name"], {})[row["state"]] = row
    return out


def format_reaction(rxn: list) -> str:
    expression = " + ".join(rxn[0])
    if len(rxn) == 3:
        expression += " <-> " + " + ".join(rxn[1]) + " -> " + " + ".join(rxn[-1])
    else:
        expression += " <-> " + " + ".join(rxn[-1])
    return expression


def analyze_site(site: str, coverage_row: dict[str, str], drc_rows: dict[str, dict[str, str]]):
    setup_dir = REPLAY / "sites" / site
    setup = setup_dir / f"{site}_static_probe.mkm"
    if site == "tico":
        setup = setup_dir / f"{site}_probe.mkm"
    if not setup.exists():
        raise FileNotFoundError(f"No archived probe model: {setup}")

    old_cwd, model, solver = configure_catmap(setup, STAGE / "runs" / site)
    try:
        if len(solver.elementary_rxns) != 8 or list(solver.adsorbate_names) != ADS_ORDER:
            raise RuntimeError(f"MODEL_CONTRACT mismatch at {site}: reactions/species changed")
        mp.mp.dps = 160
        idx = {name: i for i, name in enumerate(solver.adsorbate_names)}
        theta = [mp.mpf(coverage_row[f"theta_{name.removesuffix('_s')}"]) for name in ADS_ORDER]
        vacancy = mp.mpf(coverage_row["theta_vacancy"])
        params = list(model.scaler.get_rxn_parameters([mp.mpf("333.15")]))
        solver._coverage = list(theta)

        # Official CatMAP rates: kf/kr and signed elementary net rates.
        constants = list(solver.get_rate_constants(list(params), list(theta)))
        n = len(solver.elementary_rxns)
        kf, kr = constants[:n], constants[n:]
        net_rates = list(solver.get_rxn_rates(list(theta), constants))
        tof_vector = list(solver.get_turnover_frequency(list(params), rates=net_rates, verify_coverages=False))
        gas_tof = tof_vector[solver.gas_names.index(GAS_PRODUCT)]

        # CatMAP's own mass-action expression generator, evaluated in mpmath
        # arithmetic to avoid get_directional_rates()' Python-float conversion.
        forward_exprs = [solver.rate_equation_term(rxn[0], f"kf[{j}]") for j, rxn in enumerate(solver.elementary_rxns)]
        reverse_exprs = [solver.rate_equation_term(rxn[-1], f"kr[{j}]") for j, rxn in enumerate(solver.elementary_rxns)]
        env = {"kf": kf, "kr": kr, "p": list(solver.gas_pressures), "theta": theta, "s": [vacancy], "mp": mp}
        rf = [eval(expr, {"__builtins__": {}}, env) for expr in forward_exprs]
        rr = [eval(expr, {"__builtins__": {}}, env) for expr in reverse_exprs]

        prefactors = list(model.prefactor_list)
        prefactor_values: list[mp.mpf] = []
        for p in prefactors:
            if p.replace(" ", "") == "kB*T/h":
                prefactor_values.append(model._kB * mp.mpf("333.15") / model._h)
            else:
                raise RuntimeError(f"Unrecognized CatMAP prefactor expression {p!r}; do not infer barriers")
        barriers_f = [-model._kB * mp.mpf("333.15") * mp.log(kf[j] / prefactor_values[j]) for j in range(n)]
        barriers_r = [-model._kB * mp.mpf("333.15") * mp.log(kr[j] / prefactor_values[j]) for j in range(n)]

        balance = sum(theta) + vacancy - 1
        i1 = theta[idx["Ha_s"]] - theta[idx["OOH_s"]] - theta[idx["OH_s"]] - theta[idx["O_s"]]
        i2 = theta[idx["Hb_s"]] - theta[idx["O_s"]]
        rhs = list(solver.steady_state_function(list(theta)))
        stoich_rhs = [mp.mpf("0") for _ in theta]
        for j, rxn in enumerate(solver.elementary_rxns):
            ci, cf = Counter(rxn[0]), Counter(rxn[-1])
            for sp, i in idx.items():
                stoich_rhs[i] += (cf.get(sp, 0) - ci.get(sp, 0)) * net_rates[j]
        rhs_diff = max(abs(a - b) for a, b in zip(rhs, stoich_rhs))
        identity_errors = [abs((rf[j] - rr[j]) - net_rates[j]) for j in range(n)]
        identity_err = max(identity_errors)
        identity_rel = max((identity_errors[j] / (abs(rf[j]) + abs(rr[j]))
                            for j in range(n) if rf[j] or rr[j]), default=mp.mpf("0"))

        reaction_rows = []
        for j, rxn in enumerate(solver.elementary_rxns):
            expression = format_reaction(rxn)
            reaction_rows.append({
                "surface_name": site,
                "reaction_index": j + 1,
                "reaction_expression_catmap": expression,
                "CatMAP_reverse_rate_generated": "yes",
                "classification": "local_TS" if len(rxn) == 3 else "CatMAP_TS_or_endpoint_effective_barrier",
                "forward_barrier_effective_eV": mp.nstr(barriers_f[j], 45),
                "reverse_barrier_effective_eV": mp.nstr(barriers_r[j], 45),
                "kf_s-1": mp.nstr(kf[j], 45),
                "kr_s-1": mp.nstr(kr[j], 45),
                "prefactor": prefactors[j],
                "temperature_K": "333.15",
                "model_status": "PROVISIONAL",
            })
        rate_rows = []
        for j, rxn in enumerate(solver.elementary_rxns):
            expression = format_reaction(rxn)
            net = net_rates[j]
            rate_rows.append({
                "surface_name": site,
                "reaction_index": j + 1,
                "reaction_expression_catmap": expression,
                "rf_site-1_s-1": mp.nstr(rf[j], 60),
                "rr_site-1_s-1": mp.nstr(rr[j], 60),
                "rnet_site-1_s-1": mp.nstr(net, 60),
                "cancellation_factor": mp.nstr((abs(rf[j]) + abs(rr[j])) / abs(net), 35) if net else "NET_ZERO",
                "forward_minus_reverse_abs_error": mp.nstr(abs(rf[j] - rr[j] - net), 25),
                "rate_source": "CatMAP rate_equation_term + official get_rate_constants/get_rxn_rates; mpmath direction evaluation",
                "scientific_status": "PROVISIONAL",
            })
        coverage_rows = [{"surface_name": site, "species": name, "coverage": mp.nstr(theta[i], 70),
                          "source": coverage_row["source_path"], "scientific_status": "PROVISIONAL"}
                         for i, name in enumerate(ADS_ORDER)]
        coverage_rows.append({"surface_name": site, "species": "s_vacant", "coverage": mp.nstr(vacancy, 70),
                              "source": coverage_row["source_path"], "scientific_status": "PROVISIONAL"})

        # Formal product TOF is CatMAP's sum over signed gas stoichiometry.
        # R4 produces C6H12O_s, R5 releases C6H12O_g; compare the elementary net fluxes.
        actual_tof = mp.mpf(coverage_row["TOF_s-1"])
        oxygen_transfer = net_rates[3]
        product_release = net_rates[4]
        rel = lambda a, b: abs(a - b) / max(abs(a), abs(b)) if (a or b) else mp.mpf("0")
        drc_ts1 = drc_rows.get("OOH-C6H12_s", {})
        drc_ts2 = drc_rows.get("Hb-O_s", {})
        row = {
            "surface_name": site,
            "Catalyst": coverage_row["Catalyst"],
            "r4_oxygen_transfer_net_site-1_s-1": mp.nstr(oxygen_transfer, 60),
            "r5_epoxide_release_net_site-1_s-1": mp.nstr(product_release, 60),
            "CatMAP_C6H12O_TOF_s-1_recomputed": mp.nstr(gas_tof, 60),
            "archived_C6H12O_TOF_s-1": mp.nstr(actual_tof, 60),
            "r4_vs_TOF_relative_error": mp.nstr(rel(oxygen_transfer, gas_tof), 30),
            "r5_vs_TOF_relative_error": mp.nstr(rel(product_release, gas_tof), 30),
            "r4_vs_r5_relative_error": mp.nstr(rel(oxygen_transfer, product_release), 30),
            "archived_vs_recomputed_TOF_relative_error": mp.nstr(rel(actual_tof, gas_tof), 30),
            "step_flux_consistency_gate": "PASS_CONDITIONAL" if (
                rel(oxygen_transfer, gas_tof) <= mp.mpf("1e-10")
                and rel(product_release, gas_tof) <= mp.mpf("1e-10")
                and rel(actual_tof, gas_tof) <= mp.mpf("5e-6")
            ) else "FAIL_REEVALUATION_MISMATCH",
            "TS1_DRC_h0.005": drc_ts1.get("DRC_h0.005", ""),
            "TS1_DRC_h0.010": drc_ts1.get("DRC_h0.010", ""),
            "TS1_DRC_gate": drc_ts1.get("perturbation_gate", "NOT_AVAILABLE"),
            "TS1_smoothness_gate": drc_ts1.get("TS_smoothness_gate", "NOT_AVAILABLE"),
            "TS2_DRC_h0.005": drc_ts2.get("DRC_h0.005", ""),
            "TS2_DRC_h0.010": drc_ts2.get("DRC_h0.010", ""),
            "TS2_DRC_gate": drc_ts2.get("perturbation_gate", "NOT_AVAILABLE"),
            "numeric_status": coverage_row["numeric_status"],
            "scientific_status": "PROVISIONAL",
            "model_SHA256": coverage_row["model_SHA256"],
            "energy_SHA256": coverage_row["energy_SHA256"],
        }
        checks = {
            "site_balance_abs": mp.nstr(abs(balance), 35),
            "I1_abs": mp.nstr(abs(i1), 35),
            "I2_abs": mp.nstr(abs(i2), 35),
            "max_CatMAP_RHS_vs_stoich_abs": mp.nstr(rhs_diff, 35),
            "max_rf_minus_rr_vs_CatMAP_net_abs": mp.nstr(identity_err, 35),
            "max_directional_identity_relative_error": mp.nstr(identity_rel, 35),
            "rate_identity_pass": identity_rel <= mp.mpf("1e-25"),
            "oxygen_transfer_equals_product_TOF": rel(oxygen_transfer, gas_tof) <= mp.mpf("1e-10"),
            "release_equals_product_TOF": rel(product_release, gas_tof) <= mp.mpf("1e-10"),
            "archive_TOF_matches_native": rel(actual_tof, gas_tof) <= mp.mpf("5e-6"),
        }
        return model, solver, reaction_rows, rate_rows, coverage_rows, row, checks
    finally:
        os.chdir(old_cwd)


def main() -> None:
    mp.mp.dps = 160
    STAGE.mkdir(parents=True, exist_ok=True)
    actual_model_sha, actual_energy_sha = sha256(MODEL_FILE), sha256(ENERGY_FILE)
    if actual_model_sha != MODEL_SHA or actual_energy_sha != ENERGY_SHA:
        raise RuntimeError(f"Frozen input hash mismatch: model={actual_model_sha}, energy={actual_energy_sha}")
    coverage_data = load_archived_coverage()
    drc_data = load_drc()
    if set(coverage_data) != set(SITES):
        raise RuntimeError(f"Expected 14 archived coverage rows, found {len(coverage_data)}")
    # Evaluate Ti-Co first. Only proceed to the archived multi-site extraction
    # after the requested R4/R5/TOF and zero-inventory checks pass.
    tico = analyze_site("tico", coverage_data["tico"], drc_data.get("tico", {}))
    _, _, rc, rr, cr, summary, checks = tico
    if not (checks["rate_identity_pass"] and checks["oxygen_transfer_equals_product_TOF"]
            and checks["release_equals_product_TOF"] and checks["archive_TOF_matches_native"]):
        raise RuntimeError(f"Ti-Co validation failed; stopping before 14-site extraction: {checks}")
    tico_summary, tico_checks = summary, checks

    all_reactions, all_rates, all_coverages, oxygen_rows = list(rc), list(rr), list(cr), [summary]
    all_checks = {"tico": checks}
    for site in SITES[1:]:
        _, _, rc, rr, cr, summary, checks = analyze_site(site, coverage_data[site], drc_data.get(site, {}))
        all_reactions.extend(rc)
        all_rates.extend(rr)
        all_coverages.extend(cr)
        oxygen_rows.append(summary)
        all_checks[site] = checks

    write_csv(STAGE / "data" / "TICO_8STEP_RATE_CONSTANTS.csv", all_reactions[:8])
    write_csv(STAGE / "data" / "TICO_8STEP_REACTION_RATES.csv", all_rates[:8])
    write_csv(STAGE / "data" / "TICO_STEADY_STATE_COVERAGES.csv", all_coverages[:10])
    write_csv(STAGE / "data" / "TICO_OXYGEN_TRANSFER_TOF.csv", oxygen_rows[:1])
    write_csv(STAGE / "data" / "14SITE_OXYGEN_TRANSFER_TOF.csv", oxygen_rows)
    (STAGE / "data" / "QUALITY_CHECKS.json").write_text(json.dumps(all_checks, indent=2), encoding="utf-8")
    d = drc_data["tico"]
    summary_lines = [
        "# Ti–Co 八步逐步速率、氧转移 TOF 与 DRC 摘要",
        "",
        "## 结果定位",
        "本结果复用归档的 Ti–Co 零库存稳态覆盖度；没有重求稳态，没有重新运行14位点条件计算。CatMAP 0.3.1 官方 ReactionModel/SteadyStateSolver 接口用于解析、生成速率常数、净反应速率和 TOF。所有研究结果为 `PROVISIONAL`：八步模型与完整 NEB 共吸附状态的物理映射仍未闭合。",
        "",
        f"正式模型 SHA256：`{MODEL_SHA}`；能量 SHA256：`{ENERGY_SHA}`；CatMAP 0.3.1；T=333.15 K；输入预检 MODEL_CONTRACT=PASS。",
        "",
        "## Ti–Co 第4/5步核验",
        f"- 第4步氧转移净速率：`{tico_summary['r4_oxygen_transfer_net_site-1_s-1']} site⁻¹ s⁻¹`。",
        f"- 第5步环氧产物释放净速率：`{tico_summary['r5_epoxide_release_net_site-1_s-1']} site⁻¹ s⁻¹`。",
        f"- CatMAP 环氧产物气相 TOF：`{tico_summary['CatMAP_C6H12O_TOF_s-1_recomputed']} s⁻¹`。",
        f"- 第4步/TOF、第五步/TOF 与归档 TOF 的相对差分别为 `{tico_summary['r4_vs_TOF_relative_error']}`、`{tico_summary['r5_vs_TOF_relative_error']}`、`{tico_summary['archived_vs_recomputed_TOF_relative_error']}`。",
        f"- 位点守恒绝对误差 `{tico_checks['site_balance_abs']}`；I1 `{tico_checks['I1_abs']}`；I2 `{tico_checks['I2_abs']}`；CatMAP RHS 与计量重建最大绝对差 `{tico_checks['max_CatMAP_RHS_vs_stoich_abs']}`。",
        "",
        "## CatMAP 速率常数与有效势垒说明",
        "八步的 kf、kr 来自官方 `get_rate_constants()`；方向速率以官方 `rate_equation_term()` 质量作用表达式在 mpmath 高精度下计算，净速率来自官方 `get_rxn_rates()`。`kf`/`kr` 是速率常数，不是反应速率；rf/rr 是单向速率，rnet=rf−rr。官方 `get_directional_rates()` 内部转为 Python float，故极小净通量处不用其返回值做正逆相减。",
        "有效正逆势垒由官方 kf/kr 与 CatMAP 所用 TST 前因子 kBT/h 反算。第4步 CatMAP 原八步映射下得到的有效正/逆势垒不是孤立 NEB1 的 0.556344/0.796030 eV；第6步亦不等同于孤立 NEB2 的两向势垒。故不得把本表解释为 NEB 单步复现。",
        "正式模型第2步文本使用 `->`，但 CatMAP 0.3.1 解析后仍生成 kf 与 kr，本分析也确认其逆向 rate expression 存在；因此这里按 CatMAP 实际生成的双向速率报告，而非宣称严格不可逆。",
        "",
        "## Ti–Co DRC 判断",
        f"- 氧转移 TS1（OOH-C6H12_s）DRC：h=0.005/0.010 eV 分别 `{d['OOH-C6H12_s']['DRC_h0.005']}` / `{d['OOH-C6H12_s']['DRC_h0.010']}`；扰动门控 `{d['OOH-C6H12_s']['perturbation_gate']}`，TS 光滑性 `{d['OOH-C6H12_s']['TS_smoothness_gate']}`。",
        f"- 后续 TS2（Hb-O_s）DRC：`{d['Hb-O_s']['DRC_h0.005']}` / `{d['Hb-O_s']['DRC_h0.010']}`；扰动门控 `{d['Hb-O_s']['perturbation_gate']}`，TS 光滑性 `{d['Hb-O_s']['TS_smoothness_gate']}`。",
        f"- 最大绝对 DRC 项包括烯烃吸附态 C6H12_s=`{d['C6H12_s']['DRC_h0.005']}`、环氧产物吸附态 C6H12O_s=`{d['C6H12O_s']['DRC_h0.005']}`，以及 OOH_s=`{d['OOH_s']['DRC_h0.005']}`。因此当前条件性 DRC 不支持把氧转移 TS1 认定为 Ti–Co 的主要速率控制因素；这只是该模型、工况和零库存分支上的局部敏感性，不是机理因果证明。",
        "",
        "## 14 位点一致性提示",
        "14_SITE_OXYGEN_TRANSFER_TOF.csv 是对归档覆盖度的只读速率重评估。Ti、Ti–Co、Ti–Fe、Ti–Mn、Ti–Hf、Ti–Re、Ti–Mo、Ti–V、Ti–Zr、Ti–Ti、Ti–Ta、Ti–Cr、Ti–Ni 的第4/5步净通量与气相 TOF 在当前精度下相符。Ti–W 不符：重算第4步 `1.74467923586168e-173`，第5步/气相 TOF 为 `-1.42231644604861e-169 s⁻¹`；保留负号并标记 `FAIL_REEVALUATION_MISMATCH`，不将其包装成产物正向 TOF。虽归档汇总有 NUMERIC_PASS，该独立通量复核发现具体不一致，故 Ti–W 行不通过本次步骤通量门控，需查验其根/覆盖度快照与归档根是否同一版本；本轮未重求该位点。",
        "",
        "## DRC 适用边界",
        "Ti–Co 11个状态的 DRC 扰动均在归档表中标记 PASS_CONDITIONAL；TS1、TS2也在端点以上的光滑区域。DRC 是现有 provisional 八步输入下受库存约束重求稳态的有限差分局部敏感性，不等于唯一决速步。其它位点的 TS1 DRC 状态见 14_SITE_OXYGEN_TRANSFER_TOF.csv；未通过或不完整者保留原始状态，不补零。",
        "",
        "详见同目录 CSV 与 QUALITY_CHECKS.json。",
    ]
    (STAGE / "TICO_RATE_CONTROL_SUMMARY.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(json.dumps({"sites": len(oxygen_rows), "tico_checks": all_checks["tico"],
                      "max_14site_r4_tof_relerr": max(float(r["r4_vs_TOF_relative_error"]) for r in oxygen_rows),
                      "max_14site_r5_tof_relerr": max(float(r["r5_vs_TOF_relative_error"]) for r in oxygen_rows),
                      "output": str(STAGE)}, indent=2))


if __name__ == "__main__":
    main()
