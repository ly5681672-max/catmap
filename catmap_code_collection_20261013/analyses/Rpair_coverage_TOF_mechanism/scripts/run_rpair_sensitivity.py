"""Optional CatMAP 0.3.1 Rpair* energy and concentration sensitivity scans.

No hidden run: this script must be explicitly executed on a machine with
CatMAP 0.3.1. Originals under official_strict_tolerance remain untouched.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
ANALYSIS=Path(__file__).resolve().parents[1]
BASE=ROOT/"catmap_official_run"/"results"/"official_strict_tolerance"
BASELINE_CSV=ROOT/"analyses"/"dG3_coverage_DRC_mechanism"/"data"/"combined_dG3_coverage_DRC_TOF.csv"
KB=8.617333262145e-5 # eV K-1
TEMP=333.15
DEFAULT_SITES=("ti","tico","tife")
DELTAS=(-0.05,-0.025,0.,0.025,0.05)
C_H2O2=(0.0925,0.185,0.37)
C_ALKENE=(0.5,1.,2.)

DRIVER=r'''import json
from catmap import ReactionModel, __version__
if __version__ != "0.3.1":
    raise RuntimeError("CatMAP 0.3.1 is required; found "+str(__version__))
model=ReactionModel(setup_file="alkene_epoxidation.mkm")
model.run(recalculate=True)
tof_map=model.turnover_frequency_map
cover_map=model.coverage_map
assert len(tof_map)==1 and len(cover_map)==1
gas=list(model.gas_names)
names=list(model.adsorbate_names)+[s for s in model.site_names if s!="g"]
p=tof_map[0][1]
c=cover_map[0][1]
assert "C6H12O_g" in gas and "Rpair_s" in names
tof=p[gas.index("C6H12O_g")]
theta=c[names.index("Rpair_s")]
assert tof>0 and theta>=0
res=model.solver.get_residual(c)
assert res<=model.solver.tolerance
with open("result.json","w",encoding="utf-8") as f:
    json.dump(dict(tof=str(tof),theta_Rpair=str(theta),residual=str(res),
                   version=__version__,surface=model.surface_names[0]),f,indent=2)
'''

def set_concentrations(txt, changes):
    for key,val in changes.items():
        pattern=r"(species_definitions\['"+re.escape(key)+r"'\]\s*=\s*\{\s*'concentration':\s*)([-+0-9.eE]+)(\s*\})"
        txt,n=re.subn(pattern,lambda m:m.group(1)+repr(float(val))+m.group(3),txt)
        if n!=1:
            raise RuntimeError(f"Expected one concentration definition for {key}: found {n}")
    return txt

def perturb_rpair(src,dest,surface,delta):
    with src.open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        fieldnames=reader.fieldnames
        items=list(reader)
    changed=0
    for row in items:
        if row["surface_name"]==surface and row["species_name"]=="Rpair" and row["site_name"]=="active":
            row["formation_energy"]=repr(float(row["formation_energy"])+delta)
            changed+=1
    if changed!=1:
        raise RuntimeError(f"Expected one Rpair energy row for {surface}, found {changed}")
    with dest.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames,delimiter="\t",lineterminator="\n")
        writer.writeheader();writer.writerows(items)

def baseline_values():
    with BASELINE_CSV.open(encoding="utf-8-sig",newline="") as f:
        d=list(csv.DictReader(f))
    return {r["surface"]:float(r["TOF_s-1"]) for r in d if r["status"]=="CONVERGED"}

def run_one(surface,tag,delta,changes,out,python):
    src=BASE/surface
    if not (src/"alkene_epoxidation.mkm").is_file() or not (src/"energies.txt").is_file():
        raise FileNotFoundError(src)
    dst=out/"runs"/surface/tag
    dst.mkdir(parents=True,exist_ok=True)
    setup=(src/"alkene_epoxidation.mkm").read_text(encoding="utf-8")
    if "pressure_mode = 'concentration'" not in setup:
        raise RuntimeError("Expected strict concentration baseline; refuse non-matching input")
    if "surface_names = ['"+surface+"']" not in setup:
        raise RuntimeError("Strict setup does not match expected site")
    setup=set_concentrations(setup,changes)
    setup+="\ndata_file = 'rpair_sensitivity.pkl'\nlog_file = 'rpair_sensitivity.log'\n"
    (dst/"alkene_epoxidation.mkm").write_text(setup,encoding="utf-8")
    if abs(delta)<1e-15:
        shutil.copy2(src/"energies.txt",dst/"energies.txt")
    else:
        perturb_rpair(src/"energies.txt",dst/"energies.txt",surface,delta)
    (dst/"run_one.py").write_text(DRIVER,encoding="utf-8")
    p=subprocess.run([str(python),"-X","utf8","run_one.py"],cwd=dst,
                     capture_output=True,text=True,encoding="utf-8",
                     errors="replace",check=False)
    (dst/"console.log").write_text(p.stdout+"\n"+p.stderr,encoding="utf-8")
    if p.returncode!=0 or not (dst/"result.json").exists():
        raise RuntimeError(f"CatMAP run failed ({p.returncode}); see {dst/'console.log'}")
    d=json.loads((dst/"result.json").read_text(encoding="utf-8"))
    return float(d["tof"]),float(d["theta_Rpair"]),d["residual"]

def write_csv(file,fieldnames,rows):
    file.parent.mkdir(parents=True,exist_ok=True)
    with file.open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=fieldnames)
        w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--surfaces",nargs="+",default=list(DEFAULT_SITES))
    ap.add_argument("--python",default=sys.executable,help="Python executable with CatMAP 0.3.1")
    ap.add_argument("--out",type=Path,default=ANALYSIS/"results")
    args=ap.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    baseline=baseline_values()
    jobs=[]
    for surface in args.surfaces:
        if surface not in baseline:
            raise ValueError(f"{surface} has no CONVERGED archived baseline")
        jobs.append((surface,"baseline",0.,{}))
        jobs += [(surface,f"E_Rpair_{v:+.3f}",v,{}) for v in DELTAS if v!=0]
        jobs += [(surface,f"H2O2_{v:.5g}",0.,{"H2O2_g":v}) for v in C_H2O2 if v!=0.185]
        jobs += [(surface,f"alkene_{v:.5g}",0.,{"C6H12_g":v}) for v in C_ALKENE if v!=1]
    records=[]
    for k,(surface,tag,delta,changes) in enumerate(jobs,1):
        record=dict(surface=surface,tag=tag,delta_E_Rpair_eV=delta,
                    concentration_H2O2=changes.get("H2O2_g",0.185),
                    concentration_C6H12=changes.get("C6H12_g",1.),
                    TOF="",theta_Rpair="",residual="",status="",error="")
        try:
            tof,theta,res=run_one(surface,tag,delta,changes,out,args.python)
            if tag=="baseline":
                distance=abs(math.log10(tof/baseline[surface]))
                if distance>0.05:
                    raise RuntimeError(f"Baseline log10 TOF mismatch: {distance:.5f} dex")
            record.update(TOF=repr(tof),theta_Rpair=repr(theta),residual=res,status="CONVERGED")
            print(f"{k}/{len(jobs)} {surface} {tag}: TOF={tof:.6g}")
        except Exception as err:
            record.update(status="FAILED",error=str(err))
            print(f"{k}/{len(jobs)} FAILED {surface} {tag}: {err}",file=sys.stderr)
            records.append(record)
            write_csv(out/"scan_rows.csv",list(record),records)
            if tag=="baseline":raise RuntimeError("Baseline validation failed; stop perturbations") from err
            continue
        records.append(record)
        write_csv(out/"scan_rows.csv",list(record),records)
    sensitivity=[]
    for surface in args.surfaces:
        points={r["tag"]:r for r in records if r["surface"]==surface and r["status"]=="CONVERGED"}
        minus=points.get("E_Rpair_-0.025");plus=points.get("E_Rpair_+0.025")
        if not minus or not plus:continue
        dq=0.025
        drc=-KB*TEMP*(math.log(float(plus["TOF"]))-math.log(float(minus["TOF"])))/(2*dq)
        sensitivity.append(dict(surface=surface,numerical_Rpair_DRC=drc,energy_delta_eV=dq))
    write_csv(out/"energy_sensitivity.csv",["surface","numerical_Rpair_DRC","energy_delta_eV"],sensitivity)
    print("Scan records:",out/"scan_rows.csv")
    print("Rpair finite-difference energy sensitivity:",out/"energy_sensitivity.csv")
    print("These are new numerical results ONLY if the above CatMAP runs have actually completed.")

if __name__=="__main__":
    main()
