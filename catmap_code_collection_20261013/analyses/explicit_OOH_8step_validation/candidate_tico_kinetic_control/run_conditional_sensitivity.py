"""Run conditional CatMAP TS DRC and endpoint-barrier sensitivities."""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import mpmath as mp

LINKS = {
  1: ("* + H2O2_g <-> H2O2*", {"H":2,"O":2}, 1),
  2: ("H2O2* + * -> Ha* + OOH*", {"H":2,"O":2}, 2),
  3: ("* + C6H12_g <-> C6H12*", {"C":6,"H":12}, 1),
  4: ("Ha* + OOH* + C6H12* <-> NEB1IS* + 2*", {"C":6,"H":14,"O":2}, 3),
  6: ("NEB1FS* <-> NEB2IS* + C6H12O_g", {"C":6,"H":14,"O":2}, 1),
  8: ("NEB2FS* <-> H2O*", {"H":2,"O":1}, 1),
  9: ("H2O* <-> H2O_g + *", {"H":2,"O":1}, 1),
}
H_STEPS = (0.005, 0.010)


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f: return list(csv.DictReader(f))


def coverage_map(path): return {r["species"]: mp.mpf(r["coverage"]) for r in read_csv(path)}


def run_case(case, model_src, energy_text, solver_script, seed_cov, precision=180):
    case.mkdir(parents=True, exist_ok=True)
    (case / "candidate_energies.txt").write_text(energy_text, encoding="utf-8")
    (case / "candidate_tico.mkm").write_text(model_src, encoding="utf-8")
    cmd = [sys.executable, str(solver_script), "--model", str(case/"candidate_tico.mkm"),
           "--precision", str(precision), "--seed-coverages", str(seed_cov)]
    p = subprocess.run(cmd, cwd=case, capture_output=True, text=True, encoding="utf-8", errors="replace")
    (case/"solver_stdout.log").write_text(p.stdout, encoding="utf-8")
    (case/"solver_stderr.log").write_text(p.stderr, encoding="utf-8")
    attempts = case/"candidate_solver_attempts.json"
    state = {"exit_code": p.returncode, "status": "FAILED", "accepted_root_count": 0}
    if p.returncode == 0 and attempts.exists():
        a = json.loads(attempts.read_text(encoding="utf-8")); state["status"] = a["status"]
        roots = a.get("accepted_roots", []); state["accepted_root_count"] = len(roots)
        if roots: state["tof"] = roots[0]["product_TOF"]
    if state["status"] == "ACCEPTED":
        state["coverages"] = {k: mp.nstr(v, 70) for k,v in coverage_map(case/"candidate_steady_state_coverages.csv").items()}
    state["case"] = str(case); state["log"] = str(case/"solver_stdout.log")
    return state


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--base-run",type=Path,required=True); ap.add_argument("--solver-script",type=Path,required=True); ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args(); mp.mp.dps=120; base=a.base_run.resolve(); out=a.output.resolve(); runs=out/"runs"; runs.mkdir(parents=True,exist_ok=True)
    model_text=(base/"candidate_tico.mkm").read_text(encoding="utf-8")
    energy_source=(base/"candidate_energies.txt").read_text(encoding="utf-8-sig")
    seed_cov=base/"candidate_steady_state_coverages.csv"
    strict_base_case=runs/"strict_baseline_replay"
    strict_base=run_case(strict_base_case,model_text,energy_source,a.solver_script.resolve(),seed_cov,precision=220)
    if strict_base["status"]!="ACCEPTED":
        raise RuntimeError(f"strict baseline replay failed: {strict_base}")
    seed_cov=strict_base_case/"candidate_steady_state_coverages.csv"
    base_cov=coverage_map(seed_cov)
    base_tof=mp.mpf(strict_base["tof"])
    # Preserve the model root as the first numerical seed; do not copy old outputs.
    baseline={"status":"STRICT_BASELINE_REPLAY_ACCEPTED","tof":mp.nstr(base_tof,100),"coverages":{k:mp.nstr(v,100) for k,v in base_cov.items()},"accepted_root_count":strict_base["accepted_root_count"],"case":strict_base["case"]}
    base_energies={}
    for row in csv.DictReader(energy_source.splitlines(),delimiter="\t"):
        base_energies[row["species_name"]+"_"+row["site_name"]]=mp.mpf(row["formation_energy"])
    explicit=[]; endpoint=[]
    for name, species in (("NEB1TS_s","NEB1TS"),("NEB2TS_s","NEB2TS")):
        case_rows=[]
        for h in H_STEPS:
            for sign in (-1,1):
                delta=mp.mpf(str(sign*h)); energy=base_energies[species+"_111"]+delta
                lines=energy_source.splitlines(); rebuilt=[]
                for line in lines:
                    cells=line.split("\t")
                    if len(cells)>=4 and cells[2]==species: cells[3]=mp.nstr(energy,30)
                    rebuilt.append("\t".join(cells))
                label=f"{species}_{'plus' if sign>0 else 'minus'}_{h:.3f}eV"
                case=runs/"explicit_neb_ts"/label
                res=run_case(case,model_text,"\n".join(rebuilt)+"\n",a.solver_script.resolve(),seed_cov,precision=220)
                res.update({"ts":species,"delta_eV":mp.nstr(delta,10),"h_eV":h,"sign":sign,"energy_eV":mp.nstr(energy,25)})
                explicit.append(res)
        for h in H_STEPS:
            minus=next(r for r in explicit if r.get("ts")==species and r.get("h_eV")==h and r.get("sign")==-1)
            plus=next(r for r in explicit if r.get("ts")==species and r.get("h_eV")==h and r.get("sign")==1)
            if minus.get("status")!="ACCEPTED" or plus.get("status")!="ACCEPTED": continue
            lminus,lplus=mp.log(mp.mpf(minus["tof"])),mp.log(mp.mpf(plus["tof"]))
            kb=mp.mpf("8.617333262145e-5"); T=mp.mpf("333.15")
            drc=-kb*T*(lplus-lminus)/(2*mp.mpf(str(h)))
            branch=branch_distance(base_cov,coverage_map(Path(minus["case"])/"candidate_steady_state_coverages.csv"),coverage_map(Path(plus["case"])/"candidate_steady_state_coverages.csv"))
            case_rows.append({"ts":species,"h_eV":h,"drc_central":mp.nstr(drc,40),"tof_minus":minus["tof"],"tof_plus":plus["tof"],"branch_log_coverage_distance_max":mp.nstr(branch,15),"status":"PROVISIONAL_DRC_CANDIDATE" if branch < mp.mpf("0.5") else "INVALID_BRANCH_SHIFT"})
        explicit.extend(case_rows)
    # Endpoint steps: delta=0 is the verified provisional baseline; delta>0 adds
    # the same height to forward and reverse barriers by placing scenario TS at
    # max(endpoint energies)+delta (so reaction delta-G remains invariant).
    for step,(old_expression,composition,n_sites) in LINKS.items():
        # CatMAP-derived, thermodynamically consistent endpoint state energies.
        import os
        from catmap import ReactionModel
        cwd=Path.cwd(); os.chdir(base)
        try:
            m=ReactionModel(setup_file="candidate_tico.mkm"); m.run(); ens=m.scaler.get_free_energies([mp.mpf("333.15"),mp.mpf("1.0")]); rxn=m.elementary_rxns[step-1]
            gi=m.get_state_energy(rxn[0],ens); gf=m.get_state_energy(rxn[-1],ens); top=max(gi,gf)
        finally: os.chdir(cwd)
        endpoint.append({"step":step,"delta_eV":"0","scenario":"original CatMAP endpoint max approximation","status":"BASELINE_REUSED","tof":mp.nstr(base_tof,100),"case":str(strict_base_case),"branch_log_coverage_distance_max":"0"})
        for delta in (0.005,0.010):
            ts=f"Link{step}TS"; energy=top+mp.mpf(str(delta)); ts_species=ts+"_s"
            expr=old_expression.replace("<->",f"<-> {ts}* ->",1) if "<->" in old_expression else old_expression.replace("->",f"<-> {ts}* ->",1)
            # Step 2 is irreversible in the accepted candidate baseline. Keep
            # that directionality; its extra-height test is forward-only and
            # is not a symmetric detailed-balance-preserving barrier test.
            if step==2: expr=f"H2O2* + * -> {ts}* -> Ha* + OOH*"
            config=model_text+f"\n# hypothetical additional symmetric barrier; not DFT\nrxn_expressions[{step-1}] = {expr!r}\nspecies_definitions[{ts_species!r}] = {{'composition': {composition!r}, 'n_sites': {n_sites}}}\n"
            lines=energy_source.splitlines(); lines.append("\t".join(["tico","111",ts,f"{mp.nstr(energy,25)}","[]","HYPOTHETICAL_NOT_DFT_ENDPOINT_BARRIER_SENSITIVITY"]))
            label=f"step{step}_extra_{delta:.3f}eV"; case=runs/"endpoint_link_barriers"/label
            res=run_case(case,config,"\n".join(lines)+"\n",a.solver_script.resolve(),seed_cov,precision=220)
            res.update({"step":step,"delta_eV":f"{delta:.3f}","hypothetical_TS_eV":mp.nstr(energy,25),"scenario":"HYPOTHETICAL_NOT_DFT","endpoint_max_eV":mp.nstr(top,25)})
            endpoint.append(res)
        for delta in (0.005,0.010):
            row=next(r for r in endpoint if r.get("step")==step and r.get("delta_eV")==f"{delta:.3f}" and r.get("status") not in ("BASELINE_REUSED","PROVISIONAL_ENDPOINT_SENSITIVITY"))
            if row.get("status")!="ACCEPTED":
                row["sensitivity_status"]="INVALID_SOLVE"; continue
            cov=coverage_map(Path(row["case"])/"candidate_steady_state_coverages.csv")
            branch=branch_distance(base_cov,cov,cov)
            row["branch_log_coverage_distance_max"]=mp.nstr(branch,15)
            row["one_sided_dlnTOF_ddelta_eV-1"]=mp.nstr((mp.log(mp.mpf(row["tof"]))-mp.log(base_tof))/mp.mpf(str(delta)),40)
            row["conditional_barrier_sensitivity_X"] = mp.nstr(-mp.mpf("8.617333262145e-5")*mp.mpf("333.15")*(mp.log(mp.mpf(row["tof"]))-mp.log(base_tof))/mp.mpf(str(delta)),40)
            row["sensitivity_status"]="PROVISIONAL_ONE_SIDED" if branch < mp.mpf("0.5") else "INVALID_BRANCH_SHIFT"
            if step==2:
                row["sensitivity_status"]="PROVISIONAL_FORWARD_ONLY_IRREVERSIBLE_STEP"
    # Data tables, retaining failed cases and explicit branch evidence.
    write_rows(out/"TICO_PROVISIONAL_DRC.csv",explicit)
    write_rows(out/"ENDPOINT_BARRIER_CONTROL.csv",endpoint)
    (out/"baseline_reference.json").write_text(json.dumps(baseline,indent=2),encoding="utf-8")
    print(json.dumps({"explicit_cases":len(explicit),"endpoint_cases":len(endpoint),"out":str(out)},ensure_ascii=False))


def branch_distance(base,minus,plus):
    keys=set(base)&set(minus)&set(plus); vals=[]
    for k in keys:
        if base[k]>0 and minus[k]>0 and plus[k]>0:
            vals.append(max(abs(mp.log(minus[k]/base[k])),abs(mp.log(plus[k]/base[k]))))
    return max(vals) if vals else mp.inf


def write_rows(path,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


if __name__=="__main__": main()
