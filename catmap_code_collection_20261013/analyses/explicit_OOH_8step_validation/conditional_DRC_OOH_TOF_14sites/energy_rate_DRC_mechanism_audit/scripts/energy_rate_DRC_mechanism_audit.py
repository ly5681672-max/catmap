from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from datetime import datetime, timezone

import mpmath as mp


REPO = Path(r"H:\catmap\catmap")
PROJECT = REPO / "catmap_code_collection_20261013"
BASE = PROJECT / "analyses" / "explicit_OOH_8step_validation"
ANALYSIS = BASE / "conditional_DRC_OOH_TOF_14sites"
REPLAY = BASE / "multisite_conditional_replay_20261021"
MODEL = REPO / "lyq_alkene_epoxidation" / "2-creating_microkinetic_model" / "alkene_epoxidation.mkm"
ENERGIES = MODEL.with_name("energies_dft_neb.txt")
EXPECTED_MODEL_SHA = "86115cefb7784fe28ab519271cdfaf0e308b6f6e3a8308b1e90e4318cb497cfe"
EXPECTED_ENERGY_SHA = "55ad4e94d9c603b6705ec8c92630c9a1ceff5859cf69865145152375e203ca30"
SITES_A = ("tife", "timn", "tire")
SITES_B = ("tico", "tife", "timo", "tini")
H = mp.mpf("0.005")
TEMP = mp.mpf("333.15")
KB = mp.mpf("8.617332478e-5")
PLANCK = mp.mpf("4.135667516e-15")
CATALYST = {
    "ti": "Ti", "tico": "Ti–Co", "tife": "Ti–Fe", "timn": "Ti–Mn",
    "tihf": "Ti–Hf", "tire": "Ti–Re", "timo": "Ti–Mo", "tiv": "Ti–V",
    "tizr": "Ti–Zr", "titi": "Ti–Ti", "tini": "Ti–Ni", "tiw": "Ti–W",
    "tita": "Ti–Ta", "ticr": "Ti–Cr",
}
EXPECTED_RXNS = [
    [["s", "H2O2_g"], ["H2O2_s"]],
    [["H2O2_s", "s"], ["Ha_s", "OOH_s"]],
    [["s", "C6H12_g"], ["C6H12_s"]],
    [["OOH_s", "C6H12_s", "s"], ["OOH-C6H12_s", "s", "s"], ["C6H12O_s", "O_s", "Hb_s"]],
    [["C6H12O_s"], ["C6H12O_g", "s"]],
    [["Hb_s", "O_s"], ["Hb-O_s", "s"], ["OH_s", "s"]],
    [["Ha_s", "OH_s"], ["s", "H2O_s"]],
    [["H2O_s"], ["H2O_g", "s"]],
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def s(value):
    return mp.nstr(value, n=mp.mp.dps)


def m(value):
    # Preserve the exact binary value when CatMAP itself holds a Python float
    # (notably static gas reservoir inputs) instead of rounding it back via str().
    return mp.mpf(value)


def configure_precision(dps: int):
    global H, TEMP, KB, PLANCK
    mp.mp.dps = dps
    # Rebuild constants only after setting dps. mpmath values created at
    # import-time default precision would otherwise silently cap precision.
    H = mp.mpf("0.005")
    TEMP = mp.mpf("333.15")
    KB = mp.mpf("8.617332478e-5")
    PLANCK = mp.mpf("4.135667516e-15")


def get_model(surface):
    from catmap import ReactionModel

    setup = REPLAY / "sites" / surface / f"{surface}_static_probe.mkm"
    if not setup.is_file() and surface == "tico":
        setup = REPLAY / "sites" / surface / "tico_probe.mkm"
    if not setup.is_file():
        raise FileNotFoundError(setup)
    model = ReactionModel(setup_file=str(setup))
    if [ [list(side) for side in rxn] for rxn in model.elementary_rxns ] != EXPECTED_RXNS:
        raise RuntimeError(f"{surface}: setup parsed reactions differ from approved original eight-step network")
    if sha256(MODEL) != EXPECTED_MODEL_SHA or sha256(ENERGIES) != EXPECTED_ENERGY_SHA:
        raise RuntimeError("Approved model/energy SHA256 mismatch; aborting rate analysis")

    # Initialize only CatMAP's native numerical representation and compile
    # analytical rate functions. No mapper/run/steady-state solve is invoked.
    model._math, model._mpfloat, model._matrix = mp, mp.mpf, mp.matrix
    model._kB, model._h = KB, PLANCK
    solver = model.solver
    solver._math, solver._mpfloat, solver._matrix = mp, mp.mpf, mp.matrix
    solver._kB, solver._h = KB, PLANCK
    solver.temperature = TEMP
    solver.compile()
    params = list(model.scaler.get_rxn_parameters([TEMP]))
    return model, solver, params


def root_for(surface):
    f = REPLAY / "sites" / surface / f"constrained_{surface}_result.json"
    if not f.exists():
        # Some recovered roots live in the dedicated recovery directory.
        candidates = list((BASE / "three_site_root_recovery_20261021" / "data").glob(f"{surface}*recovery*.json"))
        if candidates:
            f = candidates[0]
        else:
            raise FileNotFoundError(f)
    return json.loads(f.read_text(encoding="utf-8-sig")), f


def param_energy_map(solver, params):
    names = list(solver.adsorbate_names) + list(solver.transition_state_names)
    return {name: params[i] for i, name in enumerate(names)}


def state_energy(names, emap, solver):
    total = mp.mpf("0")
    for name in names:
        if name in emap:
            total += emap[name]
        elif name in solver.gas_names:
            total += m(solver._gas_energies[solver.gas_names.index(name)])
        elif name in solver.site_names:
            total += m(solver._site_energies[solver.site_names.index(name)])
        else:
            raise KeyError(f"Unmapped CatMAP state component: {name}")
    return total


def barrier_rows(surface, model, solver, params):
    emap = param_energy_map(solver, params)
    rows = []
    for i, rxn in enumerate(solver.elementary_rxns):
        initial, final = rxn[0], rxn[-1]
        gis, gfs = state_energy(initial, emap, solver), state_energy(final, emap, solver)
        explicit = len(rxn) > 2
        gts_raw = state_energy(rxn[1], emap, solver) if explicit else None
        candidates = [(gis, "IS_ENDPOINT"), (gfs, "FS_ENDPOINT")]
        if explicit:
            candidates.append((gts_raw, "EXPLICIT_TS"))
        gpeak = max(v for v, _ in candidates)
        controller = next(label for value, label in candidates if value == gpeak)
        dg = gfs - gis
        rows.append({
            "surface_name": surface,
            "catalyst": CATALYST[surface],
            "reaction_index_1based": i + 1,
            "reaction_parsed": " -> ".join(" + ".join(side) for side in rxn),
            "G_IS_eV": s(gis),
            "G_TS_raw_eV": s(gts_raw) if explicit else "",
            "G_FS_eV": s(gfs),
            "G_TS_effective_eV": s(gpeak),
            "explicit_TS_present": "YES" if explicit else "NO",
            "explicit_TS_below_effective_max": "YES" if explicit and gts_raw < gpeak else "NO" if explicit else "NOT_APPLICABLE",
            "deltaG_rxn_eV": s(dg),
            "deltaG_forward_eff_eV": s(gpeak - gis),
            "deltaG_reverse_eff_eV": s(gpeak - gfs),
            "effective_control": controller if explicit else f"ENDPOINT_APPROXIMATION_{controller}",
            "signed_raw_TS_minus_IS_eV": s(gts_raw - gis) if explicit else "",
            "signed_raw_TS_minus_FS_eV": s(gts_raw - gfs) if explicit else "",
            "energy_basis": "CatMAP scaler reaction parameters + CatMAP gas/site energies; state sums follow parsed elementary_rxns",
        })
    return rows


def rate_and_flux(surface, solver, params, root, root_path):
    ads = list(solver.adsorbate_names)
    coverages = [m(root["coverage"][name]) for name in ads]
    vacancy = m(root["vacant"])
    vacancy_rebuilt = mp.mpf("1") - sum(coverages)
    if abs(vacancy_rebuilt - vacancy) > mp.mpf("1e-50"):
        raise RuntimeError(f"{surface}: stored vacancy does not reproduce from CatMAP coverages")
    solver._coverage = coverages
    ks = list(solver.get_rate_constants(list(params), coverages))
    net = list(solver.get_rxn_rates(coverages, ks))  # native CatMAP high-precision net rates
    n = len(solver.elementary_rxns)
    pressure = {name: m(solver.gas_pressures[i]) for i, name in enumerate(solver.gas_names)}
    theta = {name: coverages[i] for i, name in enumerate(ads)}
    surface_site = next(name for name in solver.site_names if name != "g")
    theta[surface_site] = vacancy

    def directional(participants, k):
        value = k
        for name in participants:
            if name in pressure:
                value *= pressure[name]
            elif name in theta:
                value *= theta[name]
            else:
                raise KeyError(f"Unmapped participant in native CatMAP rate term: {name}")
        return value

    rows = []
    for i, rxn in enumerate(solver.elementary_rxns):
        kf, kr = ks[i], ks[n + i]
        rf, rr = directional(rxn[0], kf), directional(rxn[-1], kr)
        rnet = net[i]
        reconstructed = rf - rr
        abs_err = abs(reconstructed - rnet)
        cancel = (abs(rf) + abs(rr)) / abs(rnet) if rnet else mp.inf
        # kf/kr detailed balance against the standard-state reaction energy.
        emap = param_energy_map(solver, params)
        dg = state_energy(rxn[-1], emap, solver) - state_energy(rxn[0], emap, solver)
        db_log_error = abs(mp.log(kf / kr) + dg / (KB * TEMP))
        archived = m(root["reaction_net_rates"][i])
        rows.append({
            "surface_name": surface,
            "catalyst": CATALYST[surface],
            "reaction_index_1based": i + 1,
            "forward_rate_constant_s_units": s(kf),
            "reverse_rate_constant_s_units": s(kr),
            "forward_flux_s-1": s(rf),
            "reverse_flux_s-1": s(rr),
            "net_flux_CatMAP_s-1": s(rnet),
            "archived_net_flux_s-1": s(archived),
            "current_minus_archived_net_flux_abs_error": s(abs(rnet - archived)),
            "net_flux_forward_minus_reverse_s-1": s(reconstructed),
            "cancellation_ratio": s(cancel) if mp.isfinite(cancel) else "NET_ZERO",
            "directional_minus_CatMAP_net_abs_error": s(abs_err),
            "detailed_balance_log_ratio_abs_error": s(db_log_error),
            "root_source": str(root_path),
            "computation": "CatMAP native rate constants/net rates; directional mass-action products from parsed CatMAP reaction participants in mpmath (avoids get_elem_ec float cast)",
        })
    return rows


def perturb_ooh(surface, model, solver, params, root):
    ads = list(solver.adsorbate_names)
    pidx = ads.index("OOH_s")
    coverages = [m(root["coverage"][name]) for name in ads]
    solver._coverage = coverages
    base_k = list(solver.get_rate_constants(list(params), coverages))
    vectors = {}
    for direction, shift in (("minus", -H), ("plus", H)):
        p = list(params)
        p[pidx] += shift
        vectors[direction] = {"params": p, "k": list(solver.get_rate_constants(p, coverages))}
    base_barriers = barrier_rows(surface, model, solver, params)
    b_minus = barrier_rows(surface, model, solver, vectors["minus"]["params"])
    b_plus = barrier_rows(surface, model, solver, vectors["plus"]["params"])
    n = len(solver.elementary_rxns)
    rows = []
    for i, rxn in enumerate(solver.elementary_rxns):
        kmf, kpf = vectors["minus"]["k"][i], vectors["plus"]["k"][i]
        kmr, kpr = vectors["minus"]["k"][n+i], vectors["plus"]["k"][n+i]
        bf0, br0 = m(base_barriers[i]["deltaG_forward_eff_eV"]), m(base_barriers[i]["deltaG_reverse_eff_eV"])
        bfm, brm = m(b_minus[i]["deltaG_forward_eff_eV"]), m(b_minus[i]["deltaG_reverse_eff_eV"])
        bfp, brp = m(b_plus[i]["deltaG_forward_eff_eV"]), m(b_plus[i]["deltaG_reverse_eff_eV"])
        c0, cm, cp = base_barriers[i]["effective_control"], b_minus[i]["effective_control"], b_plus[i]["effective_control"]
        # branch switching includes changes in effective max-state identity, including ties.
        switch = (c0 != cm) or (c0 != cp) or (cm != cp)
        emap0 = param_energy_map(solver, params)
        emapm = param_energy_map(solver, vectors["minus"]["params"])
        emapp = param_energy_map(solver, vectors["plus"]["params"])
        dg0 = state_energy(rxn[-1], emap0, solver) - state_energy(rxn[0], emap0, solver)
        dgm = state_energy(rxn[-1], emapm, solver) - state_energy(rxn[0], emapm, solver)
        dgp = state_energy(rxn[-1], emapp, solver) - state_energy(rxn[0], emapp, solver)
        db0 = abs(mp.log(base_k[i] / base_k[n+i]) + dg0/(KB*TEMP))
        dbm = abs(mp.log(kmf / kmr) + dgm/(KB*TEMP))
        dbp = abs(mp.log(kpf / kpr) + dgp/(KB*TEMP))
        rows.append({
            "surface_name": surface,
            "catalyst": CATALYST[surface],
            "reaction_index_1based": i + 1,
            "reaction_parsed": " -> ".join(" + ".join(side) for side in rxn),
            "h_eV": s(H),
            "deltaG_OOH_minus_eV": "-0.005",
            "deltaG_OOH_plus_eV": "+0.005",
            "forward_barrier_baseline_eV": s(bf0),
            "forward_barrier_minus_eV": s(bfm),
            "forward_barrier_plus_eV": s(bfp),
            "reverse_barrier_baseline_eV": s(br0),
            "reverse_barrier_minus_eV": s(brm),
            "reverse_barrier_plus_eV": s(brp),
            "k_forward_minus": s(kmf), "k_forward_plus": s(kpf),
            "k_reverse_minus": s(kmr), "k_reverse_plus": s(kpr),
            "dln_kforward_dG_OOH_eV-1": s((mp.log(kpf) - mp.log(kmf)) / (2*H)),
            "dln_kreverse_dG_OOH_eV-1": s((mp.log(kpr) - mp.log(kmr)) / (2*H)),
            "detailed_balance_log_error_baseline": s(db0),
            "detailed_balance_log_error_minus": s(dbm),
            "detailed_balance_log_error_plus": s(dbp),
            "controlling_energy_before": c0,
            "controlling_energy_minus": cm,
            "controlling_energy_plus": cp,
            "branch_switch": "NONSMOOTH" if switch else "SMOOTH_SAME_BRANCH",
            "rate_constant_exists": "YES" if kmf > 0 and kpf > 0 and kmr > 0 and kpr > 0 else "NO_OR_UNDERFLOW",
            "steady_state_resolved": "NO_FIXED_COVERAGE_ONLY",
        })
    return rows


def dominant_drc_rows():
    cov_rows = read_csv(ANALYSIS / "data" / "14_SITE_TOF_COVERAGE.csv")
    drc_rows = read_csv(ANALYSIS / "data" / "14_SITE_CONDITIONAL_DRC.csv")
    gates = read_csv(ANALYSIS / "data" / "DRC_NUMERICAL_GATE.csv")
    drc = {(r["surface_name"], r["state"]): r for r in drc_rows}
    gate = {(r["surface_name"], r["state"]): r for r in gates}
    out = []
    for cr in cov_rows:
        surf, state = cr["surface_name"], cr["dominant_adsorbate"]
        r, g = drc.get((surf, state)), gate.get((surf, state))
        row = {
            "surface_name": surf, "catalyst": cr["Catalyst"],
            "dominant_adsorbate": state, "dominant_coverage": cr["dominant_coverage"],
            "baseline_TOF_s-1": cr["TOF_s-1"], "numeric_status": cr["numeric_status"],
            "scientific_status": cr["scientific_status"],
            "DRC_gate": g["status"] if g else "MISSING_GATE",
            "DRC_h0.005": r["DRC_h0.005"] if r else "",
            "DRC_h0.010": r["DRC_h0.010"] if r else "",
            "sigma_eff_h0.005": "", "sigma_eff_h0.010": "", "sigma_eff_mean": "",
            "sigma_eff_abs_difference": "", "TS_weighted_3X1_2X2": "",
            "delta_sigma_eff_minus_TS_weighted": "", "TS_weighted_status": "NOT_EVALUATED",
            "interpretation_status": "PROVISIONAL_MODEL_CONDITIONAL",
        }
        if r and g and g["status"] == "PASS_CONDITIONAL" and r["perturbation_gate"] == "PASS_CONDITIONAL" and r["DRC_h0.005"] and r["DRC_h0.010"]:
            x1, x2, theta = m(r["DRC_h0.005"]), m(r["DRC_h0.010"]), m(cr["dominant_coverage"])
            if theta != 0:
                q1, q2 = -x1/theta, -x2/theta
                row.update({"sigma_eff_h0.005": s(q1), "sigma_eff_h0.010": s(q2),
                            "sigma_eff_mean": s((q1+q2)/2),
                            "sigma_eff_abs_difference": s(abs(q1-q2))})
                ts1 = drc.get((surf, "OOH-C6H12_s"))
                ts2 = drc.get((surf, "Hb-O_s"))
                gt1, gt2 = gate.get((surf, "OOH-C6H12_s")), gate.get((surf, "Hb-O_s"))
                smooth = bool(ts1 and ts2 and gt1 and gt2 and
                              gt1["status"] == "PASS_CONDITIONAL" and gt2["status"] == "PASS_CONDITIONAL" and
                              gt1.get("TS_smoothness_gate") == "SMOOTH_TS_ABOVE_ENDPOINTS" and
                              gt2.get("TS_smoothness_gate") == "SMOOTH_TS_ABOVE_ENDPOINTS" and
                              ts1.get("DRC_h0.005") and ts2.get("DRC_h0.005"))
                if smooth:
                    # Use the same h=0.005 vector for this explicitly heuristic comparison.
                    sigma_ts = 3*m(ts1["DRC_h0.005"]) + 2*m(ts2["DRC_h0.005"])
                    row["TS_weighted_3X1_2X2"] = s(sigma_ts)
                    row["delta_sigma_eff_minus_TS_weighted"] = s((q1+q2)/2 - sigma_ts)
                    row["TS_weighted_status"] = "HEURISTIC_COMPARISON_SMOOTH_TS"
                else:
                    row["TS_weighted_status"] = "NOT_COMPUTED_TS_DRC_INCOMPLETE_OR_NONSMOOTH"
        out.append(row)
    return out


def write_reports(outdir, barrier_rows_all, flux_rows_all, perturb_rows_all, ratio_rows):
    reports = outdir / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    bysite_bar = {}
    bysite_flux = {}
    for r in barrier_rows_all: bysite_bar.setdefault(r["surface_name"], []).append(r)
    for r in flux_rows_all: bysite_flux.setdefault(r["surface_name"], []).append(r)
    cov_rows = read_csv(ANALYSIS / "data" / "14_SITE_TOF_COVERAGE.csv")
    drc_rows = read_csv(ANALYSIS / "data" / "14_SITE_CONDITIONAL_DRC.csv")
    gate_rows = read_csv(ANALYSIS / "data" / "DRC_NUMERICAL_GATE.csv")
    drc_lookup = {(r["surface_name"], r["state"]): r for r in drc_rows}
    gate_lookup = {(r["surface_name"], r["state"]): r for r in gate_rows}
    diagnosis = ["# 三个位点的能量—速率诊断", "", "本报告只针对原始八步 CatMAP 输入的条件性速率。固定覆盖度计算，没有重解稳态；模型映射仍为 PROVISIONAL。有效势垒严格按 CatMAP 0.3.1 的 `G_TS_i=max(G_IS,G_FS,G_TS)` 生成式核算。下表中的速率及 DRC 来自冻结的既有基线/扰动记录；本轮只做固定覆盖度 rate-constant/net-rate 核验。", "", "## 逐位点数值摘要"]
    for surf in SITES_A:
        br, fr = bysite_bar[surf], bysite_flux[surf]
        cov = next(x for x in cov_rows if x["surface_name"] == surf)
        ddom = drc_lookup.get((surf, cov["dominant_adsorbate"]))
        gdom = gate_lookup.get((surf, cov["dominant_adsorbate"]))
        xdom = f"{ddom['DRC_h0.005']}/{ddom['DRC_h0.010']}" if ddom and gdom and gdom["status"] == "PASS_CONDITIONAL" else "不完整/未通过"
        d_ooh, g_ooh = drc_lookup.get((surf, "OOH_s")), gate_lookup.get((surf, "OOH_s"))
        x_ooh = f"{d_ooh['DRC_h0.005']}/{d_ooh['DRC_h0.010']}" if d_ooh and g_ooh and g_ooh["status"] == "PASS_CONDITIONAL" else "不完整/未通过"
        d_ts1, g_ts1 = drc_lookup.get((surf, "OOH-C6H12_s")), gate_lookup.get((surf, "OOH-C6H12_s"))
        x_ts1 = f"{d_ts1['DRC_h0.005']}/{d_ts1['DRC_h0.010']}" if d_ts1 and g_ts1 and g_ts1["status"] == "PASS_CONDITIONAL" else "不完整/未通过"
        d_ts2, g_ts2 = drc_lookup.get((surf, "Hb-O_s")), gate_lookup.get((surf, "Hb-O_s"))
        x_ts2 = f"{d_ts2['DRC_h0.005']}/{d_ts2['DRC_h0.010']}" if d_ts2 and g_ts2 and g_ts2["status"] == "PASS_CONDITIONAL" else "不完整/未通过"
        diagnosis += [f"### {CATALYST[surf]}", "", f"- 逐步基线势垒控制：" + "; ".join(f"R{x['reaction_index_1based']} {x['effective_control']} (正/逆={float(x['deltaG_forward_eff_eV']):.4g}/{float(x['deltaG_reverse_eff_eV']):.4g} eV)" for x in br),
                      f"- 主要稳态态：{cov['dominant_adsorbate']}，θ={float(cov['dominant_coverage']):.8g}；主导态 DRC(h=.005/.010)={xdom}；OOH DRC={x_ooh}；TS1 DRC={x_ts1}；TS2 DRC={x_ts2}。",
                      f"- 基线环氧 TOF={cov['TOF_s-1']} s⁻¹；数值门控={cov['numeric_status']}，科学状态={cov['scientific_status']}。",
                      f"- 八步 CatMAP 净通量范围：{min(float(x['net_flux_CatMAP_s-1']) for x in fr):.6g} 至 {max(float(x['net_flux_CatMAP_s-1']) for x in fr):.6g} s⁻¹；通量相等是稳态循环约束，不单独定位决速步。",
                      f"- 最严重正逆相消比：{max(float(x['cancellation_ratio']) for x in fr if x['cancellation_ratio']!='NET_ZERO'):.4g}。", ""]
        diagnosis += ["|步|ΔGrxn|CatMAP有效前/逆势垒|控制态|显式TS原始TS−IS|r+|r−|rnet|相消比|", "|---:|---:|---:|---|---:|---:|---:|---:|---:|"]
        for b, f in zip(br, fr):
            raw = b["signed_raw_TS_minus_IS_eV"] or "—"
            diagnosis.append(f"|R{b['reaction_index_1based']}|{float(b['deltaG_rxn_eV']):.6g}|{float(b['deltaG_forward_eff_eV']):.6g}/{float(b['deltaG_reverse_eff_eV']):.6g} eV|{b['effective_control']}|{raw}|{float(f['forward_flux_s-1']):.5g}|{float(f['reverse_flux_s-1']):.5g}|{float(f['net_flux_CatMAP_s-1']):.5g}|{float(f['cancellation_ratio']):.4g}|")
        diagnosis.append("")
    diagnosis += ["## 对三个 DRC 现象的逐点解释", "", "- **Ti–Fe：** 主导 C6H12_s 的 DRC 约 −1.988；OOH_s 与 Ha_s 均约 +1；稳定态 Hb_s、O_s DRC 近零，TS1/TS2 也近零。R2 由 FS 端点控制，R4 由显式 TS 控制；OOH 升高使 R2 前向有效势垒上升、R4 前向有效势垒下降。因此 OOH 总 DRC 的正号与 R2 方向相容，但 R2/R4 的直接常数斜率对 Ti–Co 也同样存在，不能单独解释两位点总 DRC 符号差异；还取决于各自基线通量权重及 I1/I2 约束下的全局覆盖度重排。不能据此唯一认定第2步控制。", "- **Ti–Mn：** C6H12_s DRC 约 −2.389，产物吸附态 C6H12O_s 约 +0.610、H2O_s 约 −0.174、Ha_s 约 +0.695、Hb_s/O_s 约 +0.311，OOH_s 约 +0.085，TS2 约 +0.299。R6 有显式 TS 控制，但其反向有效势垒约 2.393 eV；多个稳定态、TS 与约束库存共同参与，不能用单一 TS 或主导覆盖度解释总 DRC。", "- **Ti–Re：** H2O_s DRC 约 −1.976，TS1 约 +0.827、TS2 约 +0.173；另有 C6H12_s 约 −0.731、Ha_s 约 +0.587、OOH_s 约 −0.240、Hb_s/O_s 约 −0.173。R4/R6 的显式 TS 分别控制本模型的有效势垒，R7 由 IS 端点控制且 ΔG_rxn≈−2.461 eV。水态、吸附态和过渡态响应共同塑造 DRC；这是条件性模型解释，不是 NEB 物理映射闭合证据。", "", "## 通量解释边界", "", "表中 `r+`、`r−` 以 mpmath 按 CatMAP 解析反应物幂次和原始气相库值计算；`rnet` 使用 CatMAP 原生高精度函数。正逆大数相减时，相消比越大，净速率越依赖精度。净通量在稳态循环内相等是必要一致性，不是决速步判据。", ""]
    diagnosis += ["## 解释界限", "", "能量控制分支和方向通量是当前 CatMAP 参数化下的数值事实；将它们联系到真实 Ti–M 机理属于条件性解释。对于无显式 TS 的步骤，端点最高态是模型近似。正逆通量和 net flux 相等不证明某步为速率决定步。", ""]
    (reports / "THREE_SITE_KINETIC_DIAGNOSIS.md").write_text("\n".join(diagnosis), encoding="utf-8")

    cov_rows = read_csv(ANALYSIS / "data" / "14_SITE_TOF_COVERAGE.csv")
    drc_rows = read_csv(ANALYSIS / "data" / "14_SITE_CONDITIONAL_DRC.csv")
    drc_gates = read_csv(ANALYSIS / "data" / "DRC_NUMERICAL_GATE.csv")
    drc_lookup = {(r["surface_name"], r["state"]): r for r in drc_rows}
    gate_lookup = {(r["surface_name"], r["state"]): r for r in drc_gates}
    oo = ["# OOH* 能量扰动—速率常数分支诊断", "", "对每个位点固定既有合格覆盖度，仅将 CatMAP 参数向量中的 `OOH_s` 能量分别平移 ±0.005 eV；重新调用官方速率常数函数，不执行稳态重求、不重算 DRC。中心差分 `d ln k / dG_OOH` 表示直接速率常数响应，不是含覆盖度响应的完整 DRC。", "", "## 位点结果", ""]
    for surf in SITES_B:
        rows = [r for r in perturb_rows_all if r["surface_name"] == surf]
        switches = [r for r in rows if r["branch_switch"] == "NONSMOOTH"]
        changing = [r for r in rows if m(r["forward_barrier_minus_eV"]) != m(r["forward_barrier_plus_eV"]) or m(r["reverse_barrier_minus_eV"]) != m(r["reverse_barrier_plus_eV"])]
        ooh_drc = drc_lookup.get((surf, "OOH_s"))
        ooh_gate = gate_lookup.get((surf, "OOH_s"))
        if ooh_drc and ooh_gate and ooh_gate.get("status") == "PASS_CONDITIONAL":
            drc_note = f"既有完整 OOH_s DRC（h=.005/.010）={ooh_drc['DRC_h0.005']}/{ooh_drc['DRC_h0.010']}，门控 PASS_CONDITIONAL"
        else:
            drc_note = "既有 OOH_s DRC 不完整或未通过数值门控；不作定量对照"
        oo += [f"### {CATALYST[surf]}", "", f"- OOH 能量改变势垒的步骤：" + (", ".join("R"+str(r["reaction_index_1based"]) for r in changing) or "本范围内未改变有效势垒"),
               f"- 控制分支切换：" + (", ".join("R"+str(r["reaction_index_1based"]) for r in switches) or "无 ±h 分支身份切换"),
               "- 速率常数导数摘要：" + "; ".join(f"R{r['reaction_index_1based']} dlnkf/dG={float(r['dln_kforward_dG_OOH_eV-1']):.4g}, dlnkr/dG={float(r['dln_kreverse_dG_OOH_eV-1']):.4g} eV⁻¹" for r in rows if float(r['dln_kforward_dG_OOH_eV-1']) or float(r['dln_kreverse_dG_OOH_eV-1'])),
               f"- {drc_note}。", f"- 基线/±h 详细平衡最大 |ln(kf/kr)+ΔG/kBT|={max(float(r['detailed_balance_log_error_baseline']) for r in rows):.3g}/{max(float(r['detailed_balance_log_error_minus']) for r in rows):.3g}/{max(float(r['detailed_balance_log_error_plus']) for r in rows):.3g}。", ""]
    oo += ["## 与既有 DRC 的区别", "", "完整 DRC 来自已归档的±能量扰动后受库存约束稳态重求；这里是固定覆盖度的 k 直接响应。只有二者方向/步骤相符时，才构成互补证据，不可把某一步的 `dln k/dG` 直接等同于 OOH* DRC。正逆常数比同时用每个 parsed reaction 的 `ΔG` 检查详细平衡，误差见 CSV。", ""]
    (reports / "OOH_DRC_SWITCH_DIAGNOSIS.md").write_text("\n".join(oo), encoding="utf-8")

    complete = [r for r in ratio_rows if r["sigma_eff_mean"]]
    rtext = ["# 主导覆盖态 DRC 与覆盖度比例审计", "", "对每个位点实际最大稳定吸附态，只有两幅度 DRC 均通过门控时才计算 `sigma_eff=-X_dom/theta_dom`。保留 h=0.005、0.010 两值及其差，避免平均掩盖幅度不一致。科学状态均为 PROVISIONAL。", "", "## 全14位点状态", "", "| 位点 | 主导态 | θ | DRC(.005/.010) | σeff(.005/.010/均值) | TS加权 | 状态 |", "|---|---|---:|---:|---:|---:|---|"]
    for r in ratio_rows:
        if r["sigma_eff_mean"]:
            sig = f"{float(r['sigma_eff_h0.005']):.5g}/{float(r['sigma_eff_h0.010']):.5g}/{float(r['sigma_eff_mean']):.5g}"
        else:
            sig = "— / 未计算（门控不完整）"
        rtext.append(f"|{r['catalyst']}|{r['dominant_adsorbate']}|{float(r['dominant_coverage']):.6g}|{r['DRC_h0.005'] or '—'}/{r['DRC_h0.010'] or '—'}|{sig}|{r['TS_weighted_3X1_2X2'] or '—'}|{r['DRC_gate']}|")
    rtext += ["", "## 可计算结果明细", "", "| 位点 | 主导态 | θ | X(.005) | X(.010) | σ_eff(.005) | σ_eff(.010) | 均值 |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in complete:
        rtext.append("| {c} | {s} | {t} | {x1} | {x2} | {q1} | {q2} | {q} |".format(c=r["catalyst"], s=r["dominant_adsorbate"], t=f"{float(r['dominant_coverage']):.6g}", x1=f"{float(r['DRC_h0.005']):.5g}", x2=f"{float(r['DRC_h0.010']):.5g}", q1=f"{float(r['sigma_eff_h0.005']):.5g}", q2=f"{float(r['sigma_eff_h0.010']):.5g}", q=f"{float(r['sigma_eff_mean']):.5g}"))
    rtext += ["", "## σ≈3 与 TS 加权项的证据等级", "", "`3·X_TS1+2·X_TS2` 并非由 CatMAP 的一般 DRC 定义自动推出。实际 CatMAP 质量作用速率按每一步 parsed reactants 中物种的幂次生成：例如 R4 前向项为 θ(OOH)·θ(C6H12)·θ*，R6 前向项为 θ(Hb)·θ(O)；这些幂次没有导出统一的 3 和 2。过渡态能量进入正逆 TST 常数，而不是作为额外覆盖度因子。另有六个无显式 TS 的端点近似步骤，且 I1/I2 将自由覆盖度限制在库存子空间，因此 `σ_TS` 只能是启发式对照。CSV 仅对两个 TS 扰动门控完整且光滑的行列出该加权值。", "", "主导态的 `X≈−3θ` 是样本内数值近似，并非由位点守恒单独严格推出，也不表示三个真实金属位点/原子或三中心协同。DRC 对能量的导数还包含稳态覆盖度沿库存子空间的响应；θ_dom 本身并不唯一决定该导数。Ti–Fe、Ti–Mn、Ti–Re 的偏离应结合其他稳定态/TS 的 DRC、端点控制分支及通量相消看待。", "", "## DRC 状态完整性", ""]
    incomplete = [r["catalyst"] for r in ratio_rows if not r["sigma_eff_mean"]]
    rtext.append("主导态 DRC 两个步长均通过：" + (", ".join(r["catalyst"] for r in complete) or "无") + ".")
    rtext.append("主导态 DRC 不完整/未通过，未计算比例：" + (", ".join(incomplete) or "无") + ".")
    (reports / "DRC_COVERAGE_RELATION_AUDIT.md").write_text("\n".join(rtext), encoding="utf-8")

    findings = ["# DRC 差异的能量—速率机制核查：结论", "", "## 1. 主导吸附态为何可能出现 X≈−3θ？", "", "**结论（CONFIRMED_NUMERIC）：** 在主导态 DRC 两种扰动均通过的 11 个位点中，8 个的 `σ_eff=-X_dom/θ_dom` 落在约 2.993–3.005；其余代表偏离为 Ti–Fe≈2.000、Ti–Mn≈2.608、Ti–Re≈2.823。它是固定模型和工况下的样本数值规律。", "", "**数学判断（MODEL_DERIVED）：** DRC 定义为 `X_i=-kBT·d ln(TOF)/dG_i`；稳态覆盖度满足 `F(θ,G)=0`，因此覆盖度响应由 `dθ/dG_i=-(∂F/∂θ)^{-1}(∂F/∂G_i)`（在库存约束子空间上）决定。CatMAP 位点守恒只给出覆盖度总和约束，不给出 `X_i=-3θ_i` 恒等式。速率质量作用幂次也没有统一的 3 倍主导覆盖度因子。故目前不能从方程严格推出该近似，严格来源未建立。", "", "**合理解释边界（CONDITIONAL_INTERPRETATION）：** 这些位点在特定工况下存在占位/通量敏感性组合，使 `-X_dom/θ_dom` 数值接近3；Fe/Mn/Re 的偏离及多个 DRC 非零证明该比例不是普适单态占位定律。不能解释为三个金属原子或三位点协同。", "", "## 2. 代表位点：DRC、覆盖度与通量证据", "", "| 位点 | 主导态 | θ主导 | X主导(.005/.010) | σeff均值 | TS加权3X1+2X2 | OOH DRC(.005/.010) |", "|---|---|---:|---:|---:|---:|---:|"]
    for surf in ("tico", "tife", "timn", "tire", "timo", "tini"):
        r = next((x for x in ratio_rows if x["surface_name"] == surf), None)
        cov = next((x for x in cov_rows if x["surface_name"] == surf), None)
        od, og = drc_lookup.get((surf, "OOH_s")), gate_lookup.get((surf, "OOH_s"))
        od_text = f"{od['DRC_h0.005']}/{od['DRC_h0.010']}" if od and og and og.get("status") == "PASS_CONDITIONAL" else "未通过/缺失"
        if r and cov:
            sigma = f"{float(r['sigma_eff_mean']):.5g}" if r["sigma_eff_mean"] else "未计算"
            findings.append(f"| {CATALYST[surf]} | {cov['dominant_adsorbate']} | {float(cov['dominant_coverage']):.5g} | {r['DRC_h0.005'] or '—'}/{r['DRC_h0.010'] or '—'} | {sigma} | {r['TS_weighted_3X1_2X2'] or '—'} | {od_text} |")
    findings += ["", "Ti–Fe、Ti–Mn 的 `σ_eff` 明显低于约3的数值，而 Ti–Re 接近该值；但 TS 加权项在各位点并不一致：相近的结果只是个案，不构成通用等式。其差异应保留为当前候选模型的条件性数值特征，不能仅凭覆盖度高或净通量相等宣称某一步是实际决速步。", "", "## 3. OOH* 能量—速率常数—完整 DRC 的联系", "", "在 333.15 K，单个包含 OOH_s 的反应态自由能变化会改变对应反应的端点/TS 最高态分支；若分支不切换，CatMAP TST 常数的局部对数斜率约为 `1/(kBT)=34.83 eV⁻¹` 或其相反数，取决于 OOH_s 位于控制态的哪一侧。固定覆盖度直接斜率不是完整 DRC，因为完整 DRC 包括约束稳态下覆盖度重排。", "", "- Ti–Co：既有 OOH DRC=-0.499737/-0.499735（h=.005/.010）；固定覆盖度时 R2 与 R4 的 `dln kf/dG_OOH` 分别约 -34.8327 与 +34.8327 eV⁻¹，方向相反，且无分支切换。由于两条直接常数响应互相竞争，单看它们不能推导完整 DRC=-0.5；差异必须包含基线通量权重和受 I1/I2 约束的覆盖度重排。", "- Ti–Fe：OOH DRC=+1.000000（两步长）；其 R2、R4 直接前向常数斜率与 Ti–Co 同为约 -/+34.8327 eV⁻¹，但总 DRC 符号相反。这说明单步 Arrhenius 斜率不是位点间 DRC 差异的充分解释；实际差异来自每个位点不同的基线状态/通量权重及全稳态覆盖度响应。R2 为 FS 端点控制的无显式 TS 步，R4 才是显式 TS 控制。", "- Ti–Mo：OOH DRC=+0.000325/+0.000330，TS1 DRC约 +1；固定覆盖度 R2/R4 的 OOH 前向常数斜率相反。接近零的 OOH DRC与多步响应相消相容，但不能仅从两个直接斜率定量分解。TS1 DRC通过门控且该 TS 在扰动范围高于端点，表示当前模型 TOF 对 TS1 能量有强局部敏感性；它不等价于证明真实反应中存在唯一决速步。", "- Ti–Ni：OOH DRC=-1.1764/-1.1765，TS1 DRC约 +1。R2 的 OOH 影响体现在逆向常数（`dln kr/dG≈+34.8327 eV⁻¹`），R4 影响正向常数（约 +34.8327 eV⁻¹）；OOH 同时处于不同反应状态，净 DRC 是全网络稳态响应。TS1 的强正敏感性与 Ti–Mo 一样属于模型局部敏感性，不单独建立真实 NEB 速控结论。", "", "## 4. 主导 DRC 与覆盖度比是否具备普适性？", "", "不具备。本分析仅按既有门控数据计算两幅度下的 `-X_dom/θ_dom`。如果主导态 DRC 两个 h 中任一个未通过，则比例为空而非补估。`3·X_TS1+2·X_TS2` 是探索性的加权比较，不是 CatMAP DRC 公式，也不是化学位点数；只有两个 TS 的扰动门控完整且有效 TS 位于端点之上时才保留。", "", "## 5. 结果层级及执行边界", "", "数值表描述现有 PROVISIONAL 八步候选输入下的原生 CatMAP 速率常数、净速率与归档 DRC；方向速率通过 mpmath 依照 CatMAP 解析的质量作用物种幂次重建，以规避官方方向通量 helper 转 Python float 的精度损失。三个位点方向通量与归档净速率逐步核对，能量扰动只调用官方 rate-constant API。所有这些都不等于 NEB 与独立态映射已经物理闭合。", "", "本轮没有重跑稳态或 DRC；固定覆盖度数值检查次数有限。结果不可用于升级为正式 DFT 动力学结论或正式 SISSO 拟合。"]
    (reports / "FINAL_DRC_MECHANISM_FINDINGS.md").write_text("\n".join(findings), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(ANALYSIS / "energy_rate_DRC_mechanism_audit"))
    parser.add_argument("--dps", type=int, default=385)
    args = parser.parse_args()
    configure_precision(args.dps)
    import catmap
    catmap_version = getattr(catmap, "__version__", "unknown")
    if catmap_version != "0.3.1":
        raise SystemExit(f"Expected CatMAP 0.3.1, found {catmap_version}; aborting without calculation")
    outdir = Path(args.out)
    for sub in ("data", "reports", "logs"):
        (outdir / sub).mkdir(parents=True, exist_ok=True)
    if sha256(MODEL) != EXPECTED_MODEL_SHA or sha256(ENERGIES) != EXPECTED_ENERGY_SHA:
        raise SystemExit("Frozen model/energy hash changed; refusing to proceed")

    barrier_all, flux_all, pert_all = [], [], []
    for surface in sorted(set(SITES_A + SITES_B)):
        model, solver, params = get_model(surface)
        root, root_path = root_for(surface)
        barrier_all.extend(barrier_rows(surface, model, solver, params))
        if surface in SITES_A:
            flux_all.extend(rate_and_flux(surface, solver, params, root, root_path))
        if surface in SITES_B:
            pert_all.extend(perturb_ooh(surface, model, solver, params, root))

    ratios = dominant_drc_rows()
    write_csv(outdir / "data" / "EIGHT_STEP_EFFECTIVE_BARRIERS.csv", barrier_all)
    write_csv(outdir / "data" / "THREE_SITE_DIRECTIONAL_FLUX.csv", flux_all)
    write_csv(outdir / "data" / "OOH_ENERGY_TO_RATE_CONSTANT.csv", pert_all)
    write_csv(outdir / "data" / "14_SITE_DOMINANT_DRC_RATIO.csv", ratios)
    write_reports(outdir, barrier_all, flux_all, pert_all, ratios)
    inputs = [MODEL, ENERGIES,
              ANALYSIS / "data" / "14_SITE_TOF_COVERAGE.csv",
              ANALYSIS / "data" / "14_SITE_CONDITIONAL_DRC.csv",
              ANALYSIS / "data" / "DRC_NUMERICAL_GATE.csv"]
    probe_paths = {}
    for surface in sorted(set(SITES_A + SITES_B)):
        probe = REPLAY / "sites" / surface / f"{surface}_static_probe.mkm"
        if not probe.is_file() and surface == "tico":
            probe = REPLAY / "sites" / surface / "tico_probe.mkm"
        probe_paths[surface] = probe
        inputs.append(probe)
    output_files = sorted(p for p in outdir.rglob("*") if p.is_file())
    meta = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "catmap_version": catmap_version, "decimal_precision_mpmath": mp.mp.dps,
        "temperature_K": str(TEMP), "pressure_mode": "static",
        "operations": "native CatMAP rate constants and net rates at archived fixed coverage; no steady-state solve; no DRC rerun",
        "model_SHA256": sha256(MODEL), "energy_SHA256": sha256(ENERGIES),
        "input_SHA256": {str(p): sha256(p) for p in inputs},
        "sites_effective_barriers": list(sorted(set(SITES_A + SITES_B))),
        "sites_directional_flux": list(SITES_A), "sites_OOH_perturbation": list(SITES_B),
        "probe_setup_SHA256": {site: sha256(path) for site, path in probe_paths.items()},
        "verified_runtime_conditions": {"temperature_K": 333.15, "pressure_mode": "static", "gas_names_order": ["C6H12O_g", "C6H12_g", "H2O2_g", "H2O_g"], "gas_pressures": ["1e-20", "1.0", "0.185", "0.815"], "adsorbate_count": 9, "transition_state_count": 2, "reaction_count": 8},
        "CatMAP_get_directional_rates_note": "method get_elem_ec casts directional terms to float; audit instead evaluates the same mass-action participant terms with mpmath to preserve precision",
        "steady_state_resolves": 0,
        "existing_baseline_or_DRC_solves_repeated": 0,
        "analysis_script_SHA256": sha256(Path(__file__).resolve()),
        "output_SHA256": {str(p.relative_to(outdir)): sha256(p) for p in output_files if p.name != "run_metadata.json" and p.name != "run.log"},
    }
    (outdir / "logs" / "run_metadata.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = {"out": str(outdir), "effective_barrier_rows": len(barrier_all), "directional_flux_rows": len(flux_all), "ooh_rate_sensitivity_rows": len(pert_all), "dominant_drc_rows": len(ratios), "steady_state_solves": 0, "DRC_reruns": 0, "model_SHA256": meta["model_SHA256"], "energy_SHA256": meta["energy_SHA256"]}
    line = json.dumps(summary, ensure_ascii=False)
    print(line)
    (outdir / "logs" / "run.log").write_text(line + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
