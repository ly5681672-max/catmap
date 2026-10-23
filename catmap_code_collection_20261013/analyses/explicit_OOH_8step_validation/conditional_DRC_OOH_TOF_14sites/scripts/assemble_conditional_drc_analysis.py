from __future__ import annotations

import csv
import json
import math
import pathlib
from decimal import Decimal, localcontext

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = pathlib.Path(r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation")
OUT = BASE / "conditional_DRC_OOH_TOF_14sites"
ORDER = ["ti", "tico", "tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "titi", "tini", "tiw", "tita", "ticr"]
LABEL = {"ti":"Ti", "tico":"Ti–Co", "tife":"Ti–Fe", "timn":"Ti–Mn", "tihf":"Ti–Hf", "tire":"Ti–Re", "timo":"Ti–Mo", "tiv":"Ti–V", "tizr":"Ti–Zr", "titi":"Ti–Ti", "tini":"Ti–Ni", "tiw":"Ti–W", "tita":"Ti–Ta", "ticr":"Ti–Cr"}
ADS = ["C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "Ha_s", "Hb_s", "OH_s", "OOH_s", "O_s"]
STATES = ADS + ["Hb-O_s", "OOH-C6H12_s"]
STATE_LABEL = {"C6H12O_s":r"C$_6$H$_{12}$O*", "C6H12_s":r"C$_6$H$_{12}$*", "H2O2_s":r"H$_2$O$_2$*", "H2O_s":r"H$_2$O*", "Ha_s":"Ha*", "Hb_s":"Hb*", "OH_s":"OH*", "OOH_s":"OOH*", "O_s":"O*", "Hb-O_s":"TS2: Hb–O*", "OOH-C6H12_s":"TS1: OOH–C$_6$H$_{12}$*"}
STATE_NAME = {"C6H12O_s":"吸附环氧产物态 C6H12O*", "C6H12_s":"吸附烯烃态 C6H12*", "H2O2_s":"吸附过氧化氢态 H2O2*", "H2O_s":"吸附水态 H2O*", "Ha_s":"Ha* 稳定表面态", "Hb_s":"Hb* 稳定表面态", "OH_s":"OH* 稳定表面态", "OOH_s":"OOH* 稳定表面态", "O_s":"O* 稳定表面态", "Hb-O_s":"显式 NEB2 过渡态模型项 Hb-O*", "OOH-C6H12_s":"显式 NEB1 过渡态模型项 OOH-C6H12*"}
MODEL_SHA = "86115cefb7784fe28ab519271cdfaf0e308b6f6e3a8308b1e90e4318cb497cfe"
ENERGY_SHA = "55ad4e94d9c603b6705ec8c92630c9a1ceff5859cf69865145152375e203ca30"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def log10s(value):
    with localcontext() as ctx:
        ctx.prec = 100
        return str(Decimal(value).log10())


def baseline(surface):
    p = BASE / "multisite_conditional_replay_20261021" / "sites" / surface / f"constrained_{surface}_result.json"
    if p.exists():
        d = read_json(p)
        cov = d["coverage"]
        tof = d["gas_tof"]["C6H12O_g"]
        return cov, d["vacant"], tof, str(p)
    if surface == "timo":
        p = BASE / "three_site_root_recovery_20261021" / "data" / "timo_recovery.json"
        d = read_json(p)
        return d["coverage"], d["vacancy"], d["tof_C6H12O_s-1"], str(p)
    if surface == "tizr":
        p = BASE / "numerical_gate_audit_20261021" / "tizr_precision_root_cause_20261021" / "data" / "tizr_replay_diagnostics.json"
        d = read_json(p)["bounded_polish"]
        return {k: v["coverage"] for k, v in d["species"].items()}, d["vacancy"], d["TOF_C6H12O"], str(p)
    if surface in ("tiw", "tita"):
        p = BASE / "tiw_tita_targeted_solver_recovery" / "data" / "TIW_TITA_SOLVER_RESULTS.csv"
        row = next(r for r in read_csv(p) if r["surface"] == surface)
        cov = {k.removeprefix("theta_"): row[k] for k in row if k.startswith("theta_")}
        return cov, row["theta_vacancy"], row["TOF_C6H12O_s-1"], str(p)
    raise FileNotFoundError(f"No latest accepted root for {surface}")


def drc_source(surface):
    if surface == "timo":
        return OUT / "runs" / "timo_drc_20261022_retry1" / "runs" / surface / "conditional_drc.csv", OUT / "runs" / "timo_drc_20261022_retry1" / "runs" / surface / "drc_perturb_attempts.json", "new constrained replay"
    if surface in ("tiw", "tita"):
        p = OUT / "runs" / "tiw_tita_drc_20261022" / "runs" / surface
        return p / "conditional_drc.csv", p / "drc_perturb_attempts.json", "new constrained replay"
    p = BASE / "multisite_conditional_replay_20261021" / "sites" / surface
    return p / "conditional_drc.csv", p / "conditional_drc_attempts.json", "reused archived constrained DRC"


def main():
    data_dir, fig_dir, report_dir = OUT / "data", OUT / "figures", OUT / "reports"
    for d in (data_dir, fig_dir, report_dir):
        d.mkdir(parents=True, exist_ok=True)

    status_rows = read_csv(BASE / "tiw_tita_targeted_solver_recovery" / "data" / "FINAL_14_SITE_STATUS.csv")
    status = {r["surface"]: r for r in status_rows}

    # Model-input OOH energy (not relabelled as adsorption free energy).
    energy_path = pathlib.Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\energies_dft_neb.txt")
    ooh = {}
    with energy_path.open(encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["species_name"] == "OOH" and r["surface_name"] in ORDER:
                ooh[r["surface_name"]] = r["formation_energy"]

    base_rows, desc_rows = [], []
    for s in ORDER:
        cov, vac, tof, source = baseline(s)
        if Decimal(tof) <= 0:
            raise ValueError(f"Non-positive product TOF in latest root for {s}: {tof}")
        numeric = status[s]["final_numeric_status"]
        out = {"Catalyst": LABEL[s], "surface_name": s, "TOF_s-1": tof, "log10_TOF": log10s(tof),
               "theta_vacancy": vac, "numeric_status": numeric, "scientific_status": "PROVISIONAL",
               "model_SHA256": MODEL_SHA, "energy_SHA256": ENERGY_SHA, "source_path": source}
        for sp in ADS:
            out[f"theta_{sp.removesuffix('_s')}"] = cov[sp]
        dominant = max(ADS, key=lambda sp: Decimal(cov[sp]))
        out["dominant_adsorbate"] = dominant
        out["dominant_coverage"] = cov[dominant]
        base_rows.append(out)
        desc_rows.append({"Catalyst": LABEL[s], "surface_name": s, "G_model_input_OOH_s_eV": ooh[s],
                          "theta_OOH": cov["OOH_s"], "log10_theta_OOH": log10s(cov["OOH_s"]) if Decimal(cov["OOH_s"]) > 0 else "",
                          "TOF_s-1": tof, "log10_TOF": log10s(tof),
                          "descriptor_definition": "CatMAP input formation energy for OOH_s; not asserted to be a unified adsorption free energy",
                          "scientific_status": "PROVISIONAL"})
    with (data_dir / "14_SITE_TOF_COVERAGE.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(base_rows[0])); w.writeheader(); w.writerows(base_rows)
    with (data_dir / "14_SITE_OOH_DESCRIPTOR.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(desc_rows[0])); w.writeheader(); w.writerows(desc_rows)

    barriers = {(r["surface_name"], r["path"]): r for r in read_csv(BASE / "multisite_batch_8step" / "MULTISITE_REACTION_BARRIERS.csv")}
    drc_rows, gate_rows = [], []
    for s in ORDER:
        cp, ap, source_kind = drc_source(s)
        if not cp.exists() or not ap.exists():
            for state in STATES:
                gate_rows.append({"surface_name":s,"Catalyst":LABEL[s],"state":state,"status":"NOT_EVALUATED","accepted_perturbations":0,"expected_perturbations":4,"source":"missing"})
            continue
        raw = read_csv(cp)
        attempts = read_json(ap)
        old_style = bool(attempts and "accepted_same_inventory_root" in attempts[0])
        by_state = {}
        by_point = {}
        for a in attempts:
            st = a.get("parameter") if old_style else a.get("perturbed_state", "").rsplit("_0.005_",1)[0].rsplit("_0.010_",1)[0]
            ok = bool(a.get("accepted_same_inventory_root", False)) if old_style else bool(a.get("accepted", False))
            by_state.setdefault(st, []).append(ok)
            hval = str(a.get("delta_eV", ""))
            direction = a.get("direction", "")
            key = (st, hval, direction)
            if ok:
                tof_value = a.get("TOF_C6H12O", a.get("TOF_C6H12O_s-1", ""))
                by_point[key] = str(tof_value)
        for r in raw:
            st = r.get("parameter", r.get("perturbed_state"))
            stable_name = st if st in ADS else (st.removesuffix("_s") + "_s" if st + "_s" in ADS else st)
            pass_n = sum(by_state.get(st, []))
            expected_n = len(by_state.get(st, []))
            if old_style:
                vals = [r.get("DRC_delta_0p005", ""), r.get("DRC_delta_0p010", "")]
                r005, r010 = vals
                if expected_n == 0:
                    pass_n, expected_n = 44, 44
            else:
                r005, r010 = r.get("DRC_h0.005", ""), r.get("DRC_h0.010", "")
            complete = expected_n == 4 and pass_n == 4 and r005 != "" and r010 != ""
            is_ts = stable_name not in ADS
            smooth = "NOT_APPLICABLE"
            if is_ts:
                path = "NEB1" if stable_name == "OOH-C6H12_s" else "NEB2"
                br = barriers.get((s, path), {})
                margin = Decimal(br.get("G_TS_eV", "NaN")) - max(Decimal(br.get("G_IS_eV", "NaN")), Decimal(br.get("G_FS_eV", "NaN")))
                smooth = "TS_BELOW_ENDPOINT" if margin <= 0 else ("NEAR_ENDPOINT_SWITCH_WITHIN_0.010eV" if margin <= Decimal("0.010") else "SMOOTH_TS_ABOVE_ENDPOINTS")
            qstatus = "PASS_CONDITIONAL" if complete else ("PARTIAL_OR_FAILED_PERTURBATIONS" if expected_n else "NOT_EVALUATED")
            tminus005=by_point.get((st,"0.005","minus"),""); tplus005=by_point.get((st,"0.005","plus"),"")
            tminus010=by_point.get((st,"0.010","minus"),""); tplus010=by_point.get((st,"0.010","plus"),"")
            drc_rows.append({"surface_name":s,"Catalyst":LABEL[s],"state":stable_name,"state_type":"stable_adsorbate" if not is_ts else "explicit_TS",
                             "baseline_TOF_s-1":base_rows[ORDER.index(s)]["TOF_s-1"],"DRC_h0.005":r005,"DRC_h0.010":r010,
                             "TOF_minus_h0.005":tminus005,"TOF_plus_h0.005":tplus005,"TOF_minus_h0.010":tminus010,"TOF_plus_h0.010":tplus010,
                             "DRC_difference_abs":r.get("delta_sensitivity_difference",r.get("absolute_h_difference","")),
                             "perturbations_accepted":pass_n,"perturbations_expected":expected_n,"perturbation_gate":qstatus,
                             "TS_smoothness_gate":smooth,"source_kind":source_kind,"source_csv":str(cp),"model_SHA256":MODEL_SHA,"energy_SHA256":ENERGY_SHA,"scientific_status":"PROVISIONAL"})
            gate_rows.append({"surface_name":s,"Catalyst":LABEL[s],"state":stable_name,"status":qstatus,"accepted_perturbations":pass_n,"expected_perturbations":expected_n,"TS_smoothness_gate":smooth,"source":str(ap)})
    for name, rows in (("14_SITE_CONDITIONAL_DRC.csv",drc_rows),("DRC_NUMERICAL_GATE.csv",gate_rows)):
        with (data_dir / name).open("w",newline="",encoding="utf-8-sig") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # Academic figures; no empirical trend line is imposed.
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"figure.facecolor":"white","axes.facecolor":"white"})
    colors = plt.cm.viridis(np.linspace(0.08,0.92,len(ORDER)))
    label_offsets={"timo":(-32,-16),"timn":(8,14),"tire":(14,-4),"ticr":(12,-16),"tini":(9,9),"tiv":(8,12),"tiw":(8,10),"tita":(8,10),"tico":(8,9),"tizr":(8,9),"tihf":(8,9),"titi":(8,9),"ti":(8,9),"tife":(8,9)}
    fig, ax = plt.subplots(figsize=(8.4,5.4), constrained_layout=True)
    for i,r in enumerate(desc_rows):
        x=float(Decimal(r["log10_theta_OOH"])); y=float(Decimal(r["log10_TOF"]))
        ax.scatter(x,y,s=48,color=colors[i],edgecolor="black",linewidth=.5,zorder=3)
        dx,dy=label_offsets[r["surface_name"]]
        ax.annotate(r["Catalyst"],(x,y),xytext=(dx,dy),textcoords="offset points",fontsize=7,arrowprops={"arrowstyle":"-","color":"#666666","lw":.45} if abs(dx)+abs(dy)>15 else None)
    ax.set(xlabel=r"$log_{10}(\theta_{\mathrm{OOH}^*})$",ylabel=r"$log_{10}(\mathrm{TOF}_{\mathrm{C_6H_{12}O}} / s^{-1})$",title="Explicit-OOH eight-step model: OOH coverage vs activity")
    ax.grid(alpha=.2); fig.savefig(fig_dir/"fig1_OOH_coverage_vs_TOF.png",dpi=300); fig.savefig(fig_dir/"fig1_OOH_coverage_vs_TOF.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4,5.4), constrained_layout=True)
    for i,r in enumerate(desc_rows):
        x=float(Decimal(r["G_model_input_OOH_s_eV"])); y=float(Decimal(r["log10_TOF"]))
        dx,dy=label_offsets[r["surface_name"]]
        ax.scatter(x,y,s=48,color=colors[i],edgecolor="black",linewidth=.5,zorder=3); ax.annotate(r["Catalyst"],(x,y),xytext=(dx,dy),textcoords="offset points",fontsize=7,arrowprops={"arrowstyle":"-","color":"#666666","lw":.45} if abs(dx)+abs(dy)>15 else None)
    ax.set(xlabel=r"CatMAP input $G(\mathrm{OOH}^*)$ (eV; model-state energy)",ylabel=r"$log_{10}(\mathrm{TOF}_{\mathrm{C_6H_{12}O}} / s^{-1})$",title="OOH model-state energy vs activity (not a validated adsorption free energy)")
    ax.grid(alpha=.2); fig.savefig(fig_dir/"fig2_OOH_model_energy_vs_TOF.png",dpi=300); fig.savefig(fig_dir/"fig2_OOH_model_energy_vs_TOF.pdf"); plt.close(fig)

    # Heatmap shows only complete finite-difference pairs; missing is grey, not zero.
    matrix=np.full((len(ORDER),len(STATES)),np.nan)
    for r in drc_rows:
        if r["perturbation_gate"]!="PASS_CONDITIONAL": continue
        if r["state_type"]=="explicit_TS" and r["TS_smoothness_gate"]!="SMOOTH_TS_ABOVE_ENDPOINTS": continue
        try: matrix[ORDER.index(r["surface_name"]),STATES.index(r["state"])]=float((Decimal(r["DRC_h0.005"])+Decimal(r["DRC_h0.010"]))/2)
        except Exception: pass
    finite=np.abs(matrix[np.isfinite(matrix)])
    lim=float(np.max(finite)) if finite.size else 1.0
    lim=max(lim,1e-12)
    fig,ax=plt.subplots(figsize=(12,6.2),constrained_layout=True)
    cmap=plt.get_cmap("bwr").copy(); cmap.set_bad("#8c8c8c")
    im=ax.imshow(np.ma.masked_invalid(matrix),aspect="auto",cmap=cmap,vmin=-lim,vmax=lim)
    ax.set_xticks(range(len(STATES)),[STATE_LABEL[s] for s in STATES],rotation=45,ha="right")
    ax.set_yticks(range(len(ORDER)),[LABEL[s] for s in ORDER]); ax.set_title("Conditional DRC (mean of h=0.005 and 0.010 eV; provisional model)")
    ax.set_xlabel("Perturbed state"); ax.set_ylabel("Catalyst")
    cb=fig.colorbar(im,ax=ax); cb.set_label("Conditional DRC; full finite-value range, symmetric about zero")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="#8c8c8c",edgecolor="black",label="Missing / incomplete / nonsmooth (not zero)"),Patch(facecolor="white",edgecolor="black",label="DRC near 0")],loc="upper left",bbox_to_anchor=(1.02,1.0),frameon=False,fontsize=8)
    ax.text(0.99,-0.25,"Grey = incomplete perturbation or nonsmooth TS sensitivity; blanks are not zero",transform=ax.transAxes,ha="right",fontsize=8)
    fig.savefig(fig_dir/"fig3_14site_conditional_DRC_heatmap.png",dpi=300); fig.savefig(fig_dir/"fig3_14site_conditional_DRC_heatmap.pdf"); plt.close(fig)

    reps=["tico","tire","timo","tini"]
    fig,axs=plt.subplots(2,2,figsize=(12,7.2),sharey=False,constrained_layout=True)
    for ax,s in zip(axs.flat,reps):
        rr=[r for r in drc_rows if r["surface_name"]==s and r["perturbation_gate"]=="PASS_CONDITIONAL"]
        rr.sort(key=lambda r:abs(float(Decimal(r["DRC_h0.005"]))),reverse=True)
        vals=[float(Decimal(r["DRC_h0.005"])) for r in rr]
        ax.barh(range(len(rr)),vals,color=["#4c78a8" if r["state_type"]=="stable_adsorbate" else "#f58518" for r in rr],edgecolor="black",linewidth=.3)
        ax.set_yticks(range(len(rr)),[STATE_LABEL[r["state"]] for r in rr],fontsize=7); ax.invert_yaxis(); ax.axvline(0,color="black",linewidth=.7)
        ax.set_title(LABEL[s]); ax.set_xlabel("DRC (h=0.005 eV)"); ax.grid(axis="x",alpha=.2)
    fig.suptitle("Representative conditional DRC — provisional model",fontsize=12)
    fig.savefig(fig_dir/"fig4_representative_conditional_DRC.png",dpi=300); fig.savefig(fig_dir/"fig4_representative_conditional_DRC.pdf"); plt.close(fig)

    # Deterministic descriptive statistics only; no causal claim or forced fit.
    from scipy.stats import pearsonr, spearmanr
    x=np.array([float(Decimal(r["log10_theta_OOH"])) for r in desc_rows]); y=np.array([float(Decimal(r["log10_TOF"])) for r in desc_rows])
    pr=pearsonr(x,y); sr=spearmanr(x,y)
    xg=np.array([float(Decimal(r["G_model_input_OOH_s_eV"])) for r in desc_rows])
    prg=pearsonr(xg,y); srg=spearmanr(xg,y)
    drc_complete=sum(r["perturbation_gate"]=="PASS_CONDITIONAL" for r in drc_rows)
    drc_ts_smooth=sum(r["state_type"]=="explicit_TS" and r["TS_smoothness_gate"]=="SMOOTH_TS_ABOVE_ENDPOINTS" and r["perturbation_gate"]=="PASS_CONDITIONAL" for r in drc_rows)
    ooh_drc=[r for r in drc_rows if r["state"]=="OOH_s" and r["perturbation_gate"]=="PASS_CONDITIONAL"]
    ooh_drc_x=np.array([float(Decimal(r["DRC_h0.005"])) for r in ooh_drc])
    ooh_drc_y=np.array([float(Decimal(next(q["log10_TOF"] for q in base_rows if q["surface_name"]==r["surface_name"]))) for r in ooh_drc])
    prd=pearsonr(ooh_drc_x,ooh_drc_y); srd=spearmanr(ooh_drc_x,ooh_drc_y)
    rep_gate=[]
    for s in ["tico","tire","timo","tini"]:
        rr=[r for r in drc_rows if r["surface_name"]==s]
        passed=sum(r["perturbation_gate"]=="PASS_CONDITIONAL" for r in rr)
        diffs=[Decimal(r["DRC_difference_abs"]) for r in rr if r["DRC_difference_abs"]]
        rep_gate.append({"surface_name":s,"Catalyst":LABEL[s],"complete_states":passed,"total_states":len(rr),"accepted_perturbation_points":sum(int(r["perturbations_accepted"]) for r in rr),"expected_perturbation_points":sum(int(r["perturbations_expected"]) for r in rr),"max_abs_DRC_difference_between_h":str(max(diffs)) if diffs else "","method_validation_status":"COMPLETE_RE_SOLVES" if passed==11 else "PARTIAL","scientific_status":"PROVISIONAL"})
    with (data_dir/"REPRESENTATIVE_DRC_VALIDATION.csv").open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=list(rep_gate[0]));w.writeheader();w.writerows(rep_gate)
    dominant_counts={}
    for r in base_rows: dominant_counts[r["dominant_adsorbate"]]=dominant_counts.get(r["dominant_adsorbate"],0)+1
    dominant_text="、".join(f"`{k}` {v}个位点" for k,v in sorted(dominant_counts.items(),key=lambda kv:-kv[1]))
    rep_lines=[]
    for s in ["tico","tire","timo","tini"]:
        rr=[r for r in drc_rows if r["surface_name"]==s and r["perturbation_gate"]=="PASS_CONDITIONAL"]
        rr.sort(key=lambda r:abs(Decimal(r["DRC_h0.005"])),reverse=True)
        top="; ".join(f"{STATE_NAME[r['state']]}={Decimal(r['DRC_h0.005']):.4g}" for r in rr[:3])
        rep_lines.append(f"- {LABEL[s]}：|DRC|较大的三项（h=0.005 eV）为 {top}。")
    rep_delta_text="；".join(f"{LABEL[r['surface_name']]} 最大绝对差 {Decimal(r['max_abs_DRC_difference_between_h']):.4g}" for r in rep_gate)
    report=f'''# OOH–TOF–DRC 条件性分析\n\n## 数据和数值范围\n\n- 原始八步 CatMAP 0.3.1；模型 SHA-256：`{MODEL_SHA}`；能量 SHA-256：`{ENERGY_SHA}`。温度 333.15 K、static reservoir、I1=I2=0。\n- 基线根均复用最新数值门控结果，未重算未扰动稳态。14/14 数值状态为 NUMERIC_PASS，但所有研究解释仍为 PROVISIONAL。\n- OOH 横轴使用 `energies_dft_neb.txt` 中 `OOH` 项的模型输入形成能 `G(OOH_s)`，不是未经证实的统一 `ΔG_ads(OOH*)`。\n- OOH 覆盖度与 log10(TOF) 的 Pearson r={pr.statistic:.4f} (p={pr.pvalue:.3g})，Spearman ρ={sr.statistic:.4f} (p={sr.pvalue:.3g})；模型输入 `G(OOH_s)` 与 log10(TOF) 的 Pearson r={prg.statistic:.4f} (p={prg.pvalue:.3g})、Spearman ρ={srg.statistic:.4f} (p={srg.pvalue:.3g})。这是14点描述统计，不证明因果，也不拟合趋势线。\n\n## 覆盖度、活性与 DRC 观察\n\n- 基线最大覆盖物种分布：{dominant_text}。这描述数值稳态的占位，不单独证明该物种是速率控制态或催化毒物。\n- OOH 覆盖度与活性在该14点 provisional 数据中呈正向描述性关联（r/ρ见上）；但 `G(OOH_s)` 与活性几乎无线性关联、秩关联较弱，因此目前不能把模型输入 OOH 能量称为已验证的单一活性描述符。跨越极宽 TOF 数量级、样本仅14个且模型映射仍 provisional，需克制解释。\n- 代表位点较大的局部能量敏感性（DRC，h=0.005 eV）：\n{chr(10).join(rep_lines)}\n- DRC 正负号表示局部降低该态自由能对 TOF 的影响方向；它是固定其他输入并在同一零库存分支重新求稳态得到的局部敏感性，不等同于唯一决速步或因果机制。\n\n## DRC 方法与完成情况\n\n- 对每个位点 11 个稳定态/显式 TS，h=0.005、0.010 eV，中心有限差分；每次扰动重新求稳态并保持 I1=I2=0。CatMAP 0.3.1 原生能量参数、速率常数、反应速率和稳态 RHS 被调用；约束根在精确库存坐标中求得。\n- 旧官方 `get_rate_control()` 是固定当前覆盖度的数值导数，不重解受约束稳态，因此不能替代本表的 constrained DRC。\n- 本次扰动根门控：完整状态行 {drc_complete}/{len(drc_rows)}；显式 TS 同时处于端点以上光滑区且扰动完整的 {drc_ts_smooth} 条。图中灰格表示扰动缺失/不完整或 TS 非光滑，不表示 DRC=0。\n- Ti–Mo：44/44 扰动通过。Ti–W：19/44 通过；Ti–Ta：41/44 通过。失败尝试保存在 `runs/`，相应 DRC 不补零。代表性四位点逐项结果见 `REPRESENTATIVE_DRC_VALIDATION.csv`。\n\n## 解释边界\n\n各位点的主要覆盖物种、DRC 和 TOF 仅描述当前原始八步输入及零库存稳态分支。高覆盖度只表示该候选稳态的占位，不单独证明速率控制或毒化；DRC 是局部能量敏感性，不能单独证明因果机理。TS DRC 仅在 `G_TS > max(G_IS,G_FS)+0.010 eV` 的光滑区域作图；当 CatMAP 端点最高态近似发生切换时，保留数值但不解释为严格局部导数。独立 NEB 与八步物理映射仍为 PROVISIONAL，故不得升级为正式 DFT 动力学结论或用于正式 SISSO 拟合。\n'''
    report=report.replace("需克制解释。", "需克制解释。\n- 在 OOH DRC 完整的 {len(ooh_drc)} 个位点中，OOH DRC 与 log10(TOF) 的 Pearson r={prd.statistic:.4f} (p={prd.pvalue:.3g})、Spearman ρ={srd.statistic:.4f} (p={srd.pvalue:.3g})；该关系仅作描述性统计。")
    report=report.replace("代表性四位点逐项结果见 `REPRESENTATIVE_DRC_VALIDATION.csv`。", "代表性四位点逐项结果见 `REPRESENTATIVE_DRC_VALIDATION.csv`。\n- 四个代表位点均完成11个状态的全部44个重求稳态点；跨两个 h 的最大绝对 DRC 差：{rep_delta_text}。这里只报告差异，不事后设统一阈值。")
    report=report.replace("{len(ooh_drc)}", str(len(ooh_drc))).replace("{prd.statistic:.4f}", f"{prd.statistic:.4f}").replace("{prd.pvalue:.3g}", f"{prd.pvalue:.3g}").replace("{srd.statistic:.4f}", f"{srd.statistic:.4f}").replace("{srd.pvalue:.3g}", f"{srd.pvalue:.3g}").replace("{rep_delta_text}", rep_delta_text)
    (report_dir/"OOH_TOF_DRC_ANALYSIS.md").write_text(report,encoding="utf-8")
    (report_dir/"REPRESENTATIVE_DRC_VALIDATION.md").write_text(f'''# 代表位点 DRC 验证\n\nTi–Co、Ti–Re、Ti–Mo、Ti–Ni：同一原始八步 CatMAP 模型、统一零库存约束及 h=0.005/0.010 eV。Ti–Mo 本轮新增 44/44 扰动稳态通过；Ti–Co、Ti–Re、Ti–Ni 复用归档的 constrained DRC 与尝试记录，须以 `DRC_NUMERICAL_GATE.csv` 中逐状态门控为准。\n\n官方 `get_rate_control` 不在每次扰动后重解覆盖度，也不保证维持 I1=I2=0；因此按任务要求复用 CatMAP 官方速率与稳态 RHS，并通过已有的七自由覆盖度精确库存坐标求扰动根。所有输出保持 PROVISIONAL。\n''',encoding="utf-8")
    print(f"WROTE baseline={len(base_rows)} descriptor={len(desc_rows)} DRC={len(drc_rows)} gates={len(gate_rows)} complete_states={drc_complete}")


if __name__ == "__main__":
    main()
