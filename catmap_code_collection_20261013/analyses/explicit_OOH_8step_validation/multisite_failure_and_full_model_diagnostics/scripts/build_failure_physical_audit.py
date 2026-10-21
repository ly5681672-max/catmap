"""Rebuild deterministic NEB anomaly, coverage and network-mapping audits."""
from __future__ import annotations

import ast
import csv
import hashlib
import json
from collections import Counter
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from pathlib import Path

from catmap import ReactionModel

HERE = Path(__file__).resolve().parents[1]
VALIDATION = HERE.parent
BATCH = VALIDATION / "multisite_batch_8step"
ORIGINAL_MKM = HERE / "inputs" / "original_8step_snapshot" / "alkene_epoxidation.mkm"
CANDIDATE_MKM = HERE / "inputs" / "candidate_9reaction_snapshot" / "candidate_ti.mkm"
DATA = HERE / "data"
REPORTS = HERE / "reports"
TOL = 1e-8


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def classify_signed_neb(ts_minus_is: float, ts_minus_fs: float, tol: float = TOL) -> str:
    below_is = ts_minus_is < -tol
    below_fs = ts_minus_fs < -tol
    if below_is and below_fs:
        return "TS_BELOW_BOTH_ENDPOINTS"
    if below_is:
        return "TS_BELOW_IS"
    if below_fs:
        return "TS_BELOW_FS"
    return "TS_NOT_BELOW_ENDPOINTS"


def exact_nullspace(matrix: list[list[int]]) -> list[list[int]]:
    """Integer basis for null(matrix), exact rational RREF."""
    if not matrix:
        return []
    a = [[Fraction(x) for x in row] for row in matrix]
    nr, nc = len(a), len(a[0])
    pivots, r = [], 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        scale = a[r][c]
        a[r] = [x / scale for x in a[r]]
        for i in range(nr):
            if i != r and a[i][c]:
                q = a[i][c]
                a[i] = [x - q * y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    free = [c for c in range(nc) if c not in pivots]
    basis = []
    for f in free:
        v = [Fraction(0) for _ in range(nc)]
        v[f] = 1
        for rr, pc in enumerate(pivots):
            v[pc] = -a[rr][f]
        den = reduce(lcm, (x.denominator for x in v), 1)
        ints = [int(x * den) for x in v]
        common = reduce(gcd, (abs(x) for x in ints if x), 0) or 1
        ints = [x // common for x in ints]
        first = next((x for x in ints if x), 1)
        if first < 0:
            ints = [-x for x in ints]
        basis.append(ints)
    return basis


def parse_catmap_reactions(path: Path) -> tuple[ReactionModel, list[list[list[str]]]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    expr = next(ast.literal_eval(node.value) for node in tree.body
                if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "rxn_expressions" for t in node.targets))
    model = ReactionModel()
    model.parse_elementary_rxns(expr)
    return model, list(model.elementary_rxns)


def stoich_and_invariants(model: ReactionModel, rxns: list[list[list[str]]]) -> tuple[list[dict], list[dict], list[list[int]]]:
    species = list(model.adsorbate_names)
    stoich = [[0] * len(rxns) for _ in species]
    for j, rxn in enumerate(rxns):
        initial, final = Counter(rxn[0]), Counter(rxn[-1])
        for i, name in enumerate(species):
            stoich[i][j] = final.get(name, 0) - initial.get(name, 0)
    null = exact_nullspace([list(row) for row in zip(*stoich)])
    matrix_rows = [{"stable_adsorbate_identifier": name, **{f"r{j+1}": row[j] for j in range(len(rxns))}}
                   for name, row in zip(species, stoich)]
    invariant_rows = []
    for idx, vec in enumerate(null, 1):
        terms = [f"{coef:+d}·θ({name})" for coef, name in zip(vec, species) if coef]
        invariant_rows.append({"invariant_id": f"I{idx}", "expression": " ".join(terms) + " = constant",
                               **{name: coef for name, coef in zip(species, vec)}})
    return matrix_rows, invariant_rows, null


def n_sites_for(model: ReactionModel, name: str) -> int:
    return int(getattr(model, "species_definitions", {}).get(name, {}).get("n_sites", 1))


def site_and_atom_checks(model: ReactionModel, rxns: list[list[list[str]]]) -> list[dict]:
    formulas = {
        "H2O2_g": {"H": 2, "O": 2}, "H2O2_s": {"H": 2, "O": 2},
        "H2O_g": {"H": 2, "O": 1}, "H2O_s": {"H": 2, "O": 1},
        "C6H12_g": {"C": 6, "H": 12}, "C6H12_s": {"C": 6, "H": 12},
        "C6H12O_g": {"C": 6, "H": 12, "O": 1}, "C6H12O_s": {"C": 6, "H": 12, "O": 1},
        "Ha_s": {"H": 1}, "Hb_s": {"H": 1}, "O_s": {"O": 1},
        "OH_s": {"O": 1, "H": 1}, "OOH_s": {"O": 2, "H": 1},
        "Hb-O_s": {"H": 1, "O": 1}, "OOH-C6H12_s": {"C": 6, "H": 13, "O": 2},
        "NEB1IS_s": {"C": 6, "H": 14, "O": 2}, "NEB1TS_s": {"C": 6, "H": 14, "O": 2},
        "NEB1FS_s": {"C": 6, "H": 14, "O": 2},
        "NEB2IS_s": {"H": 2, "O": 1}, "NEB2TS_s": {"H": 2, "O": 1}, "NEB2FS_s": {"H": 2, "O": 1},
    }
    rows = []
    for j, rxn in enumerate(rxns, 1):
        states = []
        for state in rxn:
            atoms: Counter = Counter()
            sites = 0
            for sp in state:
                if sp == "s":
                    sites += 1
                    continue
                if sp not in formulas:
                    raise KeyError(f"No audited formula for CatMAP species {sp}")
                atoms.update({el: count for el, count in formulas[sp].items()})
                if not sp.endswith("_g"):
                    sites += n_sites_for(model, sp)
            states.append((atoms, sites))
        atom_ok = all(states[k][0] == states[k+1][0] for k in range(len(states)-1))
        site_ok = all(states[k][1] == states[k+1][1] for k in range(len(states)-1))
        rows.append({"reaction_index": j, "states": " | ".join(" + ".join(s) for s in rxn),
                     "atom_balance": "PASS" if atom_ok else "FAIL", "site_balance": "PASS" if site_ok else "FAIL",
                     "state_site_counts": ";".join(str(sites) for _, sites in states),
                     "basis": "CatMAP parsed states; all adsorbates n_sites default/declared 1; formula assignments documented in script"})
    return rows


def energy_span(site: str) -> tuple[float, float]:
    path = BATCH / "runs" / f"{site}_candidate" / f"candidate_{site}_energies.txt"
    with path.open(encoding="utf-8-sig", newline="") as f:
        values = [float(r["formation_energy"]) for r in csv.DictReader(f, delimiter="\t") if r["surface_name"] == site]
    return min(values), max(values)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    original, orig_rxns = parse_catmap_reactions(ORIGINAL_MKM)
    candidate, cand_rxns = parse_catmap_reactions(CANDIDATE_MKM)
    orig_matrix, orig_inv, _ = stoich_and_invariants(original, orig_rxns)
    cand_matrix, cand_inv, _ = stoich_and_invariants(candidate, cand_rxns)
    write_csv(DATA / "original_8step_surface_stoichiometric_matrix.csv", orig_matrix)
    write_csv(DATA / "original_8step_inventory_invariants.csv", orig_inv)
    write_csv(DATA / "candidate_9reaction_surface_stoichiometric_matrix.csv", cand_matrix)
    write_csv(DATA / "candidate_9reaction_inventory_invariants.csv", cand_inv)
    orig_checks = site_and_atom_checks(original, orig_rxns)
    cand_checks = site_and_atom_checks(candidate, cand_rxns)
    write_csv(DATA / "network_balance_checks.csv", [
        {"network": "original_8step", **r} for r in orig_checks
    ] + [{"network": "candidate_9reaction", **r} for r in cand_checks])

    mapping = [
        (1, "H2O2 adsorption", "candidate r1", "state correspondence retained"),
        (2, "H2O2 activation to Ha*/OOH*", "candidate r2", "explicit OOH* retained; native CatMAP arrow still reversible"),
        (3, "1-hexene adsorption", "candidate r3", "state correspondence retained"),
        (4, "OOH*/alkene oxygen-transfer stage", "candidate r4+r5", "independent activation-to-NEB1IS link plus full NEB1 IS/TS/FS; Ha spectator folded into composite formulas"),
        (5, "epoxide desorption", "candidate r6", "coupled to NEB1FS→NEB2IS and product release; not an independent-state mapping"),
        (6, "Hb*/O* hydrogen transfer to OH*", "candidate r7 and NEB2 composite", "Hb*, O*, OH* not independently represented; no strict state lumping proof"),
        (7, "Ha*/OH* water formation", "candidate r8", "absorbed into NEB2FS→H2O* link; Ha/Hb/OH inventory distinctions lost"),
        (8, "water desorption", "candidate r9", "state correspondence retained"),
    ]
    write_csv(DATA / "original_to_candidate_reaction_mapping.csv", [
        {"original_step": i, "formal_chemical_stage": stage, "candidate_reactions": candidate_rxns,
         "mapping_status": status} for i, stage, candidate_rxns, status in mapping])

    barrier_rows = read_csv(BATCH / "MULTISITE_REACTION_BARRIERS.csv")
    anomaly_rows = []
    for r in barrier_rows:
        a, b = float(r["signed_TS_minus_IS_eV"]), float(r["signed_TS_minus_FS_eV"])
        category = classify_signed_neb(a, b)
        g = {"IS": float(r["G_IS_eV"]), "TS": float(r["G_TS_eV"]), "FS": float(r["G_FS_eV"])}
        max_state = max(g, key=g.get)
        ef = g[max_state] - g["IS"]
        er = g[max_state] - g["FS"]
        archived_f = float(r["CatMAP_endpoint_scheme_forward_barrier_eV"])
        archived_r = float(r["CatMAP_endpoint_scheme_reverse_barrier_eV"])
        anomaly_rows.append({**r, "TS_endpoint_classification": category,
                             "TS_below_any_endpoint": category != "TS_NOT_BELOW_ENDPOINTS",
                             "endpoint_max_state": max_state, "recomputed_forward_endpoint_barrier_eV": ef,
                             "recomputed_reverse_endpoint_barrier_eV": er,
                             "forward_barrier_match": abs(ef-archived_f) < TOL,
                             "reverse_barrier_match": abs(er-archived_r) < TOL,
                             "interpretation": "raw signed NEB differences retained; endpoint maximum is CatMAP candidate barrier, not a repaired DFT TS"})
    write_csv(DATA / "neb_endpoint_anomaly_audit.csv", anomaly_rows)

    coverage_rows = read_csv(BATCH / "MULTISITE_TOF_COVERAGE.csv")
    species_columns = [
        ("theta_H2O2_s", "adsorbed hydrogen peroxide"), ("theta_Ha", "Ha-labelled surface H state"),
        ("theta_OOH", "OOH*"), ("theta_C6H12_s", "adsorbed 1-hexene-labelled state"),
        ("theta_NEB1IS_s", "NEB1 initial composite state"), ("theta_NEB1FS_s", "NEB1 final composite state"),
        ("theta_NEB2IS_s", "NEB2 initial composite state"), ("theta_NEB2FS_s", "NEB2 final composite state"),
        ("theta_H2O_s", "adsorbed water-labelled state"),
    ]
    dominant_rows = []
    for r in coverage_rows:
        parsed = [(col, label, float(r[col])) for col, label in species_columns if r.get(col, "") not in (None, "")]
        col, label, theta = max(parsed, key=lambda x: x[2])
        dominant_rows.append({"surface_name": r["surface_name"], "TOF_C6H12O_s-1": r["TOF_C6H12O_s-1"],
                              "log10_TOF": r["log10_TOF"], "theta_OOH": r["theta_OOH"], "theta_vacancy": r["theta_vacancy"],
                              "dominant_candidate_state_column": col, "dominant_state_cautious_label": label,
                              "dominant_state_coverage": theta,
                              "forward_flux_step_6": r["forward_flux_step_6"], "reverse_flux_step_6": r["reverse_flux_step_6"],
                              "net_flux_step_6": r["net_flux_step_6"], "interpretation_limit": "model-internal only; high coverage alone does not establish rate control"})
    write_csv(DATA / "successful_site_dominant_state_summary.csv", dominant_rows)

    status_rows = read_csv(BATCH / "MULTISITE_CALCULATION_STATUS.csv")
    fail_rows = []
    for r in status_rows:
        if r["surface_name"] not in {"tizr", "titi", "tini", "tiw"}:
            continue
        site = r["surface_name"]
        minimum, maximum = energy_span(site)
        attempts_path = BATCH / "runs" / f"{site}_candidate" / "candidate_solver_attempts.json"
        attempt_json = json.loads(attempts_path.read_text(encoding="utf-8")) if attempts_path.exists() else {}
        attempts = attempt_json.get("attempts", [])
        failure_log = BATCH / "logs" / f"{site}_p240.log"
        failure_text = failure_log.read_text(encoding="utf-8", errors="replace") if failure_log.exists() else ""
        fail_rows.append({"surface_name": site, "previous_status": r["numerical_status"],
                          "energy_min_eV": minimum, "energy_max_eV": maximum, "input_energy_span_eV": maximum-minimum,
                          "attempt_count_recorded": len(attempts),
                          "recorded_seed_relative_residuals": ";".join(f"{x.get('seed')}={x.get('relative_surface_residual')}" for x in attempts),
                          "observed_solver_error": ";".join(sorted({str(x.get('solver_error')) for x in attempts if x.get('solver_error')})),
                          "log_has_singular_jacobian": "numerically singular" in failure_text.lower(),
                          "interpretation": ("interrupted before solver diagnostics; current cause unknown" if site == "tizr" else
                              "multiple starts converge to stationary/non-root residual; logs do not prove a thermodynamic impossibility" if site in {"titi","tini"} else
                              "initial solver matrix numerically singular and p240 timed out; additionally NEB1 TS lies below IS, endpoint clipping gives zero forward TS barrier")})
    write_csv(DATA / "failed_site_root_cause_diagnostics.csv", fail_rows)

    provenance_rows = []
    for r in status_rows:
        site = r["surface_name"]
        if r["calculation_status"] == "ALREADY_COMPUTED_PROVISIONAL":
            provenance = "REUSED_HISTORICAL_STEADY_RESULT"
            evidence = r["run_path"]
        elif r["calculation_status"] == "COMPUTED_PROVISIONAL":
            provenance = "SOLVED_IN_THIS_BATCH"
            evidence = str(BATCH / "runs" / f"{site}_candidate" / "candidate_solver_attempts.json")
        else:
            provenance = "ATTEMPTED_NO_ACCEPTED_STEADY_RESULT"
            evidence = str(BATCH / "logs" / f"{site}_p120.log")
        provenance_rows.append({"surface_name": site, "calculation_status": r["calculation_status"],
                                "legacy_reused_existing_result_field": r["reused_existing_result"],
                                "corrected_provenance_interpretation": provenance,
                                "evidence_path": evidence,
                                "note": "legacy boolean is ambiguous for batch-computed rows; this audit uses explicit status plus per-site solver evidence"})
    write_csv(DATA / "calculation_provenance_audit.csv", provenance_rows)

    # Join only existing batch outputs by the stable surface key; no numerical
    # recomputation and no mutation of the source batch tables.
    tof_by_site = {r["surface_name"]: r for r in coverage_rows}
    prov_by_site = {r["surface_name"]: r for r in provenance_rows}
    barrier_by_site = {}
    for row in anomaly_rows:
        barrier_by_site.setdefault(row["surface_name"], {})[row["path"]] = row
    drc_rows = read_csv(BATCH / "MULTISITE_DRC_PROVISIONAL.csv")
    drc_by_site = {}
    for row in drc_rows:
        drc_by_site.setdefault(row["surface_name"], []).append(row)
    link_rows = read_csv(BATCH / "data" / "MULTISITE_LINK_ENDPOINT_SENSITIVITY.csv")
    link_by_site = {}
    for row in link_rows:
        link_by_site.setdefault(row["surface_name"], []).append(row)

    def barrier_fields(site: str, path: str, prefix: str) -> dict:
        row = barrier_by_site.get(site, {}).get(path, {})
        return {
            f"{prefix}_G_IS_eV": row.get("G_IS_eV", ""),
            f"{prefix}_G_TS_eV": row.get("G_TS_eV", ""),
            f"{prefix}_G_FS_eV": row.get("G_FS_eV", ""),
            f"{prefix}_signed_TS_minus_IS_eV": row.get("signed_TS_minus_IS_eV", ""),
            f"{prefix}_signed_TS_minus_FS_eV": row.get("signed_TS_minus_FS_eV", ""),
            f"{prefix}_effective_forward_barrier_eV": row.get("CatMAP_endpoint_scheme_forward_barrier_eV", ""),
            f"{prefix}_effective_reverse_barrier_eV": row.get("CatMAP_endpoint_scheme_reverse_barrier_eV", ""),
            f"{prefix}_TS_endpoint_classification": row.get("TS_endpoint_classification", ""),
        }

    master_rows = []
    for status in status_rows:
        site = status["surface_name"]
        tof = tof_by_site.get(site, {})
        drcs = drc_by_site.get(site, [])
        neb1_drc = next((r.get("conditional_DRC", "") for r in drcs
                         if r.get("perturbed_state") == "NEB1TS"), "")
        neb2_drc = next((r.get("conditional_DRC", "") for r in drcs
                         if r.get("perturbed_state") == "NEB2TS"), "")
        valid_drcs = [(abs(float(r["conditional_DRC"])), r["perturbed_state"], r["conditional_DRC"])
                      for r in drcs if r.get("conditional_DRC", "") not in ("", None)]
        dominant_drc = max(valid_drcs, default=(None, "", ""))
        link = link_by_site.get(site, [])
        link_statuses = sorted({r.get("status", "") for r in link if r.get("status", "")})
        master_rows.append({
            "surface_name": site,
            "calculation_status": status.get("calculation_status", ""),
            "numerical_status": status.get("numerical_status", ""),
            "physical_gate": status.get("physical_gate", ""),
            "direction_semantics": status.get("direction_semantics", ""),
            "provenance": prov_by_site.get(site, {}).get("corrected_provenance_interpretation", ""),
            "legacy_reused_existing_result_field": status.get("reused_existing_result", ""),
            "source_hash": status.get("source_hash", ""),
            "TOF_C6H12O_s-1": tof.get("TOF_C6H12O_s-1", ""),
            "log10_TOF": tof.get("log10_TOF", ""),
            "theta_OOH": tof.get("theta_OOH", ""),
            "theta_vacancy": tof.get("theta_vacancy", ""),
            "dominant_candidate_state": (max(
                ((k, float(tof[k])) for k in ("theta_H2O2_s", "theta_Ha", "theta_C6H12_s", "theta_H2O_s",
                   "theta_C6H12O_s", "theta_NEB1IS_s", "theta_NEB1FS_s", "theta_NEB2IS_s", "theta_NEB2FS_s")
                 if tof.get(k, "") not in ("", None)), key=lambda x: x[1])[0]
                if tof else ""),
            "dominant_candidate_coverage": (max(
                (float(tof[k]) for k in ("theta_H2O2_s", "theta_Ha", "theta_C6H12_s", "theta_H2O_s",
                   "theta_C6H12O_s", "theta_NEB1IS_s", "theta_NEB1FS_s", "theta_NEB2IS_s", "theta_NEB2FS_s")
                 if tof.get(k, "") not in ("", None)), default="") if tof else ""),
            **barrier_fields(site, "NEB1", "NEB1"),
            **barrier_fields(site, "NEB2", "NEB2"),
            "conditional_DRC_NEB1TS": neb1_drc,
            "conditional_DRC_NEB2TS": neb2_drc,
            "dominant_conditional_sensitivity_state": dominant_drc[1],
            "dominant_conditional_sensitivity_value": dominant_drc[2],
            "link_endpoint_sensitivity_status": ";".join(link_statuses),
            "interpretation": "candidate-model diagnostics only; not formal validated 8-step TOF/DRC",
        })
    write_csv(DATA / "multisite_candidate_master.csv", master_rows)

    source_hashes = []
    for path in (ORIGINAL_MKM, HERE / "inputs" / "original_8step_snapshot" / "energies_dft_neb.txt",
                 HERE / "inputs" / "original_8step_snapshot" / "DFT_TS_relative.csv", CANDIDATE_MKM,
                 HERE / "inputs" / "candidate_9reaction_snapshot" / "candidate_ti_energies.txt"):
        source_hashes.append({"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "size_bytes": path.stat().st_size})
    write_csv(DATA / "input_sha256_manifest.csv", source_hashes)

    orig_n = len(orig_inv)
    cand_n = len(cand_inv)
    anomalies = [r for r in anomaly_rows if r["TS_below_any_endpoint"]]
    dom_counts = Counter(r["dominant_state_cautious_label"] for r in dominant_rows)
    dominant_table = "\n".join(
        "| {surface} | {state} | {coverage:.6g} | {ooh:.6g} | {vac:.6g} | {tof:.6g} |".format(
            surface=r["surface_name"], state=r["dominant_state_cautious_label"],
            coverage=float(r["dominant_state_coverage"]), ooh=float(r["theta_OOH"]),
            vac=float(r["theta_vacancy"]), tof=float(r["TOF_C6H12O_s-1"]))
        for r in dominant_rows
    )
    report = f"""# 八步失败位点与完整模型物理映射诊断

## 摘要

- 批次内候选模型：CatMAP 官方解析为 {len(cand_rxns)} 个基元反应；原始显式网络解析为 {len(orig_rxns)} 个基元反应、八个化学阶段。
- 表面计量矩阵的守恒零空间维数：原始网络 {orig_n}，九反应候选 {cand_n}。原始网络的两个库存不变量与候选网络的单个 `Ha−OOH` 不变量不等价；因此候选网络不是已证明的严格动力学粗粒化。
- 原始 TS 低于至少一个端点的 NEB 路径：{len(anomalies)} 条，涉及 {len({r['surface_name'] for r in anomalies})} 个位点。每条记录均保留原始 TS−IS、TS−FS 和 CatMAP 端点最高态正逆势垒。
- 10 个已收敛候选中主导状态分布：{'; '.join(f'{k}: {v}' for k, v in dom_counts.items())}。

## 失败位点诊断与恢复尝试

`failed_site_root_cause_diagnostics.csv` 保留能量输入跨度、已有尝试残差、求解异常证据及可证范围。结论仅限数值：

- Ti–Zr：既有批次在求解诊断前中断；本轮恢复流程只完成 `uniform_0` 一个初猜，尚未运行 donor 初猜，也未完成完整验收；不能把它归类为已证明的 CatMAP 数值失败根因。
- Ti–Ti、Ti–Ni：多个初猜在 120/240 位精度下有很大的相对表面残差，属于未找到稳态根；当前证据不证明不存在根。变量仍采用正覆盖度 log 参数化、形式位点守恒和由计量矩阵推导的库存约束。
- Ti–W：尝试日志出现数值奇异矩阵且 240 位任务超时。其 NEB1 原始 TS−IS = −0.566651 eV；端点最高态法因此给候选正向有效势垒 0 eV。该物理异常标记与求解奇异是两类问题，不能互相替代解释。

本轮 Ti–Zr donor-seed 恢复脚本首个 uniform 初猜达到缩放残差 `4.5881e-86`，但旧相对残差门槛仍为 1，未通过。该恢复进程达到本轮时间预算后停止；没有生成可接受根或 TOF。SCIPY 初始 Jacobian 条件数约 `3.58e5`，该数值只针对这个初猜，不能外推为全局条件数。日志/脚本保留了该次有限尝试；未降低门槛。

## NEB 异常与势垒语义

7 个要求核对的位点均由自动规则得到明确标记：Ti 和 Ti–Zr 为 `TS_BELOW_FS`；Ti–W、Ti–Ta、Ti–Mo、Ti–V、Ti–Cr 为 `TS_BELOW_IS`。这是原始有符号比较，不将负值钳制或重写为真实零势垒。CatMAP 候选有效势垒按 `max(G_IS,G_TS,G_FS)` 计算，分别列出正、逆方向；端点值不是新 DFT TS。

## 主导覆盖度与候选机制解释

候选中 NEB1 终态复合态覆盖接近 1 的位点确实包括 Re、Mo、Mn、Cr、V、Ta；Hf 的最大覆盖是水标记态；Ti、Co、Fe 的最大覆盖是烯烃标记态。此分布与各自能量表和候选网络平衡相容，但复合态是粗粒化的一有效位点表达，不证明真实几何占位。

| 位点 | 候选模型最大覆盖状态 | 最大覆盖度 | OOH* 覆盖度 | 空位覆盖度 | 产物 TOF (s⁻¹) |
|---|---|---:|---:|---:|---:|
{dominant_table}

Re 与 Mo 的 NEB1 终态高覆盖同时伴随较高候选 TOF；所以“NEB1终态覆盖高必然堵塞空位并降低 TOF”不成立于这组模型数据。Ta 的空位覆盖极低，可能在候选模型内限制需要空位的连接，但仅凭覆盖度不能定为速率控制。跨位点 TOF 差异同时依赖各能量输入、端点最高态势垒、气相储库和候选连接步骤。

**连接动力学敏感性边界：**现存 `MULTISITE_LINK_ENDPOINT_SENSITIVITY.csv` 中，Ti–Co 的连接端点/假设敏感性来自历史条件性计算，Ti–Fe 仅有少量复用的端点能量敏感性；Re、Mo 对应记录均为 `NOT_RUN`。本轮没有为 Re/Mo 引入任意连接 TS 或借用 Ti–Co 势垒，因此不能据现有数据断言 Re/Mo 的 TOF 排序对连接动力学稳健。表中 `conditional_sensitivity` 也不是 CatMAP 官方 DRC。

## 原始八步与九反应候选的数学关系

实际 CatMAP parser 解析所得计量矩阵和 nullspace 系数见 `data/`。原始网络库存基底等价于：

```text
I1 = θ(Ha*) − θ(OOH*) − θ(OH*) − θ(O*)
I2 = θ(Hb*) − θ(O*)
```

九反应候选的库存只留下：

```text
I' = θ(Ha*) − θ(OOH*)
```

候选将 Hb*、O*、OH*、独立环氧吸附态及其库存变化并入 NEB 复合态/连接，不能从原始网络严格推导出此减少。两网络按所声明分子式和每个吸附物种 `n_sites=1` 做逐反应原子与形式位点计数均通过；这只验证代数平衡，不验证复合态结构身份、真实几何占位、能量连续性或连接路径的 DFT 势垒。完整矩阵、逐反应检查与状态映射见 `data/original_8step_surface_stoichiometric_matrix.csv`、`data/candidate_9reaction_surface_stoichiometric_matrix.csv`、`data/network_balance_checks.csv`、`data/original_to_candidate_reaction_mapping.csv`。

## 结论与门控

现有数据可以构成原子/形式位点守恒的九反应**条件性候选**，但没有足够依据将它认定为原始八步网络的严格粗粒化。NEB1/NEB2 原始三态可在候选中保留；缺失的连接 TS、Hb/O/OH 独立库存、复合态真实 `n_sites` 以及从洁净表面初态的可达性仍阻止物理升级。现在不能据此开展正式 14 位点 TOF/DRC 或正式构效拟合。

`MULTISITE_CALCULATION_STATUS.csv` 的 `reused_existing_result` 对部分新算位点为 True，确有口径歧义。旧表未覆盖；本次数据整理应区分“本批 CatMAP 重求解”“复用历史稳态”和“仅复用输入/预处理产物”，对应于 `calculation_provenance_audit.csv`。

本次另生成 `data/multisite_candidate_master.csv`，按 `surface_name` 将14个目标位点状态、已有覆盖度/TOF、NEB1/NEB2有符号端点差与候选有效势垒、旧的条件性TS敏感性及连接端点敏感性状态合并。空字段表示来源结果未提供该项，不以0替代；该表仍属于候选模型数据。
"""
    (REPORTS / "FAILURE_AND_PHYSICAL_MAPPING_DIAGNOSTICS.md").write_text(report, encoding="utf-8")
    print(json.dumps({"original_reactions": len(orig_rxns), "candidate_reactions": len(cand_rxns),
                      "original_inventory_nullity": orig_n, "candidate_inventory_nullity": cand_n,
                      "flagged_neb_paths": len(anomalies), "coverage_sites": len(dominant_rows),
                      "report": str(REPORTS / "FAILURE_AND_PHYSICAL_MAPPING_DIAGNOSTICS.md")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
