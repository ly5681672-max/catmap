from __future__ import annotations

import csv, json, math, os
from pathlib import Path
from collections import Counter
import mpmath as mp
from three_site_recovery import setup_model, INDEPENDENT, ROOT, TMP

OUT = Path(r"C:\Users\Administrator\Documents\codex探索\three_site_recovery_deliverables")
TARGETS = ("timo", "tiw", "tita")
SITE_MAP = {"timo":"Ti-Mo", "tiw":"Ti-W", "tita":"Ti-Ta"}
OLD_STATUS = Path(r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\numerical_gate_audit_20261021\data\OLD_VS_NEW_QUALITY_STATUS.csv")


def write_csv(path, rows):
    keys=[]
    for row in rows:
        for k in row:
            if k not in keys: keys.append(k)
    with path.open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)


def evaluate(site):
    mp.mp.dps=385
    wd=TMP/site
    res=json.loads((wd/f"{site}_recovery.json").read_text(encoding="utf-8"))
    model,solver,params,old=setup_model(site,wd)
    try:
        # reset all numerical hooks at the precision used for root validation
        mp.mp.dps=385
        model._math=mp; model._mpfloat=mp.mpf; model._matrix=mp.matrix; model._Axb_solver=mp.lu_solve
        model._kB=mp.mpf("8.617332478e-5"); model._h=mp.mpf("4.135667516e-15")
        solver.compile()
        params=list(solver.scaler.get_rxn_parameters([mp.mpf("333.15")]))
        names=list(solver.adsorbate_names); idx={n:i for i,n in enumerate(names)}
        free=[idx[n] for n in INDEPENDENT]
        derived={"Ha_s":("OOH_s","OH_s","O_s"),"Hb_s":("O_s",)}
        weights=(1,1,1,1,2,2,3)
        def state(z):
            odds=[mp.exp(v) for v in z]; den=1+sum(a*b for a,b in zip(weights,odds))
            cv=[mp.mpf(0)]*len(names)
            for i,q in zip(free,odds): cv[i]=q/den
            for target,srcs in derived.items(): cv[idx[target]]=sum(cv[idx[s]] for s in srcs)
            return cv,1/den
        root=None
        if res.get("accepted"):
            root=[mp.mpf(x) for x in res["root_logits"]]
        else:
            sci=next((a for a in reversed(res.get("attempts",[])) if a.get("polished_root_logits")),None)
            candidate=sci.get("polished_root_logits") if sci else None
            if candidate is None:
                candidates=[a for a in res.get("attempts",[]) if a.get("root_logits")]
                if candidates: candidate=min(candidates,key=lambda a:mp.mpf(a["scaled_residual_inf"])).get("root_logits")
            root=[mp.mpf(x) for x in candidate] if candidate else [mp.mpf("-4")]*len(free)
        cv,vac=state(root)
        ks=solver.get_rate_constants(params,list(cv))
        native=list(solver.get_rxn_rates(list(cv),ks))
        steady=list(solver.steady_state_function(list(cv)))
        n_rxn=len(solver.elementary_rxns)
        pressures=[mp.mpf(v) for v in solver.gas_pressures]
        empty=[]
        for s in solver.site_names:
            total=mp.mpf(solver.species_definitions[s]["total"])
            empty.append(total-sum(cv[i]*mp.mpf(solver.species_definitions[n]["n_sites"])
                                   for i,n in enumerate(names) if solver.species_definitions[n]["site"]==s))
        env={"kf":ks[:n_rxn],"kr":ks[n_rxn:],"theta":cv,"p":pressures,"s":empty}
        fwd=[]; rev=[]
        for j,rxn in enumerate(solver.elementary_rxns):
            ef=solver.rate_equation_term(rxn[0],f"kf[{j}]")
            er=solver.rate_equation_term(rxn[-1],f"kr[{j}]")
            fwd.append(eval(ef,{"__builtins__":{}},env)); rev.append(eval(er,{"__builtins__":{}},env))
        net=[a-b for a,b in zip(fwd,rev)]
        diff_native=max(abs(a-b) for a,b in zip(net,native))
        stoich=[[0]*n_rxn for _ in names]
        for j,rxn in enumerate(solver.elementary_rxns):
            ci,cf=Counter(rxn[0]),Counter(rxn[-1])
            for n in names: stoich[idx[n]][j]=cf.get(n,0)-ci.get(n,0)
        local={}
        for n in INDEPENDENT:
            i=idx[n]; gross=sum(abs(stoich[i][j])*(fwd[j]+rev[j]) for j in range(n_rxn))
            local[n]={"abs_residual":abs(steady[i]),"gross":gross,
                      "relative":abs(steady[i])/gross if gross else mp.inf}
        max_local=max(v["relative"] for v in local.values())
        i1=cv[idx["Ha_s"]]-cv[idx["OOH_s"]]-cv[idx["OH_s"]]-cv[idx["O_s"]]
        i2=cv[idx["Hb_s"]]-cv[idx["O_s"]]
        tof=list(solver.get_turnover_frequency(params,rates=native,verify_coverages=False))[solver.gas_names.index("C6H12O_g")]
        # Jacobian conditioning diagnostic in the seven-coordinate chart.
        local_scales=[local[n]["gross"] for n in INDEPENDENT]
        def F(z):
            cc,_=state(z); solver.get_rate_constants(params,list(cc))
            ff=list(solver.steady_state_function(list(cc)))
            return tuple(ff[idx[n]]/local_scales[k] for k,n in enumerate(INDEPENDENT))
        try:
            J=mp.matrix([[mp.diff(lambda q:F([q if kk==col else root[kk] for kk in range(len(root))])[rr],root[col])
                          for col in range(len(root))] for rr in range(len(root))])
            sv=mp.svd(J,compute_uv=False)
            svals=sorted([abs(x) for x in sv],reverse=True)
            cond=svals[0]/svals[-1] if svals[-1] else mp.inf
            jac={"rank_estimate":sum(x > svals[0]*mp.mpf("1e-30") for x in svals),
                 "sigma_max":mp.nstr(svals[0],25),"sigma_min":mp.nstr(svals[-1],25),
                 "condition_number":mp.nstr(cond,25)}
        except Exception as exc:
            jac={"diagnostic_error":f"{type(exc).__name__}: {exc}"}
        numeric_status=res.get("numeric_status","NUMERIC_FAIL")
        # Only Ti-Mo has a complete independent-start cross-precision replay.
        summary={"surface":site,"catalyst":SITE_MAP[site],"numeric_status":numeric_status,
                 "scientific_status":"PROVISIONAL","tof_s-1":mp.nstr(tof,50) if res.get("accepted") else "",
                 "log10_tof":mp.nstr(mp.log10(tof),35) if res.get("accepted") and tof>0 else "",
                 "theta_OOH":mp.nstr(cv[idx["OOH_s"]],35),"theta_vacancy":mp.nstr(vac,35),
                 "theta_C6H12":mp.nstr(cv[idx["C6H12_s"]],35),"theta_H2O":mp.nstr(cv[idx["H2O_s"]],35),
                 "theta_H2O2":mp.nstr(cv[idx["H2O2_s"]],35),"theta_C6H12O":mp.nstr(cv[idx["C6H12O_s"]],35),
                 "max_abs_species_residual":mp.nstr(max(abs(steady[idx[n]]) for n in INDEPENDENT),30),
                 "max_local_normalized_residual":mp.nstr(max_local,30),
                 "site_sum_error":mp.nstr(abs(sum(cv)+vac-1),25),"I1":mp.nstr(i1,25),"I2":mp.nstr(i2,25),
                 "max_abs_directional_reconstruction_error_vs_CatMAP_net":mp.nstr(diff_native,25),
                 "Jacobian":jac,"candidate_only":not bool(res.get("accepted"))}
        flux=[]
        for j,rxn in enumerate(solver.elementary_rxns):
            C=(abs(fwd[j])+abs(rev[j]))/abs(net[j]) if net[j] else mp.inf
            flux.append({"surface":site,"catalyst":SITE_MAP[site],"candidate_validity":"ACCEPTED_ROOT" if res.get("accepted") else "FAILED_SEARCH_CANDIDATE_ONLY",
                         "reaction_index":j+1,"reaction":" -> ".join([" + ".join(side) for side in (rxn[0],rxn[-1])]),
                         "forward_rate_s-1":mp.nstr(fwd[j],40),"reverse_rate_s-1":mp.nstr(rev[j],40),
                         "net_rate_s-1":mp.nstr(net[j],40),"cancellation_index":mp.nstr(C,25),
                         "CatMAP_net_rate_s-1":mp.nstr(native[j],40),"net_reconstruction_abs_error":mp.nstr(abs(net[j]-native[j]),20)})
        return summary,flux,{"site":site,"attempts":res.get("attempts",[]),"root_candidate_logits": [mp.nstr(x,80) for x in root],
                             "coverage":{n:mp.nstr(cv[i],60) for i,n in enumerate(names)},"vacancy":mp.nstr(vac,60),
                             "local_residuals":{n:{k:mp.nstr(v,35) if hasattr(v,'_mpf_') else v for k,v in vals.items()} for n,vals in local.items()},
                             "jacobian":jac,"tof_official":mp.nstr(tof,60) if res.get("accepted") else "",
                             "net_rates": [mp.nstr(x,40) for x in native],"forward": [mp.nstr(x,40) for x in fwd],"reverse": [mp.nstr(x,40) for x in rev],
                             "high_precision_directional_reconstruction_error":mp.nstr(diff_native,25)}
    finally:
        os.chdir(old)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    summaries=[]; fluxes=[]; details={}
    for s in TARGETS:
        a,b,c=evaluate(s); summaries.append(a); fluxes.extend(b); details[s]=c
    write_csv(OUT/"THREE_SITE_SOLVER_RESULTS.csv",summaries)
    write_csv(OUT/"THREE_SITE_FLUX_VALIDATION.csv",fluxes)
    (OUT/"candidate_diagnostics.json").write_text(json.dumps(details,indent=2,ensure_ascii=False),encoding="utf-8")
    gate={"model":"original explicit OOH 8-step CatMAP model","catmap_version":"0.3.1","temperature_K":333.15,
          "scientific_status":"PROVISIONAL","conditions_unchanged":True,"original_data_modified":False,
          "numeric_gate": {"local_species_relative_residual_max": "<=1e-12", "site_balance": "chart construction plus independent check", "zero_inventory": "I1=I2=0 chart", "cross_precision_TOF_relative": "<=5 ppm", "major_and_OOH_coverage_relative": "<=1e-4", "physical_result_not_upgraded":True},
          "sites":[]}
    for s,summary in zip(TARGETS,summaries):
        rr=json.loads((TMP/s/f"{s}_recovery.json").read_text(encoding="utf-8"))
        item={"surface":s,"catalyst":SITE_MAP[s],"status":summary["numeric_status"],"scientific_status":"PROVISIONAL",
              "accepted_root":bool(rr.get("accepted")),"attempt_count":len(rr.get("attempts",[])),
              "root_seed_method":rr.get("accepted_seed","none"),"cross_precision_replay":rr.get("cross_precision_replay",None),
              "limitation":"No accepted physical steady state" if not rr.get("accepted") else "Numerical validation passed on original candidate model; energy mapping remains PROVISIONAL"}
        gate["sites"].append(item)
    (OUT/"THREE_SITE_NUMERICAL_GATE.json").write_text(json.dumps(gate,indent=2,ensure_ascii=False),encoding="utf-8")
    # Update a copy of prior status, never overwrite the historical audit table.
    old=list(csv.DictReader(OLD_STATUS.open(encoding="utf-8-sig",newline="")))
    mapped={"timo":"NUMERIC_PASS","tiw":"NUMERIC_FAIL","tita":"NUMERIC_FAIL"}
    for r in old:
        s=r["surface"]
        r["current_numeric_status"]=mapped.get(s,r["new_status"])
        r["scientific_status"]="PROVISIONAL"
        if s=="timo": r["current_note"]="Recovered with bounded numerical seed search and CatMAP NewtonRoot polish; independent 250 dps replay passes 5 ppm."
        elif s in ("tiw","tita"): r["current_note"]="Bounded recovery attempted; no accepted root. See three-site root cause report."
        else: r["current_note"]="Historical numerical gate retained; no new solve this round."
    write_csv(OUT/"UPDATED_14_SITE_NUMERICAL_STATUS.csv",old)
    print(OUT)


if __name__=="__main__": main()
