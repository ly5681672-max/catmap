"""Validate the newly generated conditional Ti-Co run set and CSV outputs."""
from __future__ import annotations
import csv, json
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parent
CASES = [ROOT / "runs" / "strict_baseline_replay"]
CASES += sorted((ROOT / "runs" / "explicit_neb_ts").glob("NEB*"))
CASES += sorted((ROOT / "runs" / "endpoint_link_barriers").glob("step*"))
failures=[]; accepted=0; max_rel=mp.mpf("0"); max_scaled=mp.mpf("0"); max_site=mp.mpf("0"); max_inventory=mp.mpf("0"); max_flux_rel=mp.mpf("0")
for case in CASES:
    ap=case/"candidate_solver_attempts.json"
    if not ap.exists(): failures.append(f"missing attempts: {case.name}"); continue
    rec=json.loads(ap.read_text(encoding="utf-8"))
    if rec.get("catmap_version")!="0.3.1" or rec.get("status")!="ACCEPTED": failures.append(f"not accepted CatMAP 0.3.1: {case.name}"); continue
    accepted+=1
    for a in rec["attempts"]:
        if a.get("accepted"):
            max_rel=max(max_rel,mp.mpf(a["relative_residual"])); max_scaled=max(max_scaled,mp.mpf(a["scaled_residual_inf"]))
    with (case/"candidate_steady_state_coverages.csv").open(encoding="utf-8-sig",newline="") as f: cov={r["species"]:mp.mpf(r["coverage"]) for r in csv.DictReader(f)}
    site=abs(sum(v for k,v in cov.items() if k!="*")+cov["*"]-1); inv=abs(cov["Ha_s"]-cov["OOH_s"])
    max_site=max(max_site,site); max_inventory=max(max_inventory,inv)
    if site>mp.mpf("1e-60") or inv>mp.mpf("1e-60") or min(cov.values())<0: failures.append(f"coverage/site/inventory quality: {case.name}")
    with (case/"candidate_reaction_rates.csv").open(encoding="utf-8-sig",newline="") as f: rates=list(csv.DictReader(f))
    if len(rates)!=9: failures.append(f"reaction-rate rows != 9: {case.name}"); continue
    vals=[mp.mpf(r["net_rate"]) for r in rates]
    prod=next(mp.mpf(r["turnover_frequency"]) for r in csv.DictReader((case/"candidate_tof.csv").open(encoding="utf-8-sig",newline="")) if r["gas_species"]=="C6H12O_g")
    denom=max(abs(prod),max(map(abs,vals)),mp.mpf("1e-300"))
    fluxrel=max(abs(x-prod) for x in vals)/denom; max_flux_rel=max(max_flux_rel,fluxrel)
    if fluxrel>mp.mpf("1e-20"): failures.append(f"flux/TOF mismatch: {case.name} {mp.nstr(fluxrel,8)}")

with (ROOT/"FORWARD_REVERSE_RATE_AUDIT.csv").open(encoding="utf-8-sig",newline="") as f: rate_rows=list(csv.DictReader(f))
with (ROOT/"TICO_PROVISIONAL_DRC.csv").open(encoding="utf-8-sig",newline="") as f: drc=list(csv.DictReader(f))
with (ROOT/"EXPLICIT_TS_PERTURBATION_CASES.csv").open(encoding="utf-8-sig",newline="") as f: exp_cases=list(csv.DictReader(f))
with (ROOT/"ENDPOINT_BARRIER_CONTROL.csv").open(encoding="utf-8-sig",newline="") as f: endpoint=list(csv.DictReader(f))
if len(rate_rows)!=9: failures.append("directional audit must contain 9 reactions")
if len(drc)!=4 or len(exp_cases)!=8: failures.append("expected 4 DRC summaries plus 8 TS perturbation cases")
if len(endpoint)!=21: failures.append("expected 7 endpoint steps x (baseline, .005, .010)")
result={"catmap_version":"0.3.1","candidate_solve_cases_expected":len(CASES),"candidate_solve_cases_accepted":accepted,
        "explicit_TS_perturbation_cases":len(exp_cases),"endpoint_barrier_cases_including_baseline":len(endpoint),
        "directional_rate_rows":len(rate_rows),"DRC_summary_rows":len(drc),
        "max_accepted_relative_surface_residual":mp.nstr(max_rel,12),"max_scaled_residual_inf":mp.nstr(max_scaled,12),
        "max_site_balance_error":mp.nstr(max_site,8),"max_Ha_minus_OOH_inventory_error":mp.nstr(max_inventory,8),
        "max_reaction_net_rate_vs_product_TOF_relative_difference":mp.nstr(max_flux_rel,12),"failures":failures,
        "validation_status":"PASS" if not failures and accepted==len(CASES) else "FAIL"}
(ROOT/"validation_summary.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
print(json.dumps(result,ensure_ascii=False,indent=2))
if failures or accepted!=len(CASES): raise SystemExit(1)
