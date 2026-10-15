"""Independent official CatMAP 0.3.1 spot checks of new OOH-BEP grid.
No alteration to archived inputs. Executes only selected points in isolated dirs.
"""
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"catmap_official_run"/"results"/"official_strict_tolerance"/"tico"
DATA=ROOT/"analyses"/"BEP_energy_correction_optimization"/"data"/"source_carbon_O_H_TS1.csv"
DRIVER="""from catmap import ReactionModel,__version__
import json
assert __version__=="0.3.1",__version__
m=ReactionModel(setup_file="alkene_epoxidation.mkm")
m.run(recalculate=True)
labs=m.output_labels.get("turnover_frequency") or list(m.gas_names)
idx=labs.index("C6H12O_g")
tof=m.turnover_frequency_map[0][1][idx]
assert tof>0
with open("official_result.json","w",encoding="utf-8") as f:
 json.dump({"tof":str(tof),"version":__version__},f)
"""
def main():
    if not SOURCE.exists():raise RuntimeError("Strict Ti-Co CatMAP baseline missing")
    with DATA.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==13
    ref=next(r for r in rows if r["surface_name"]=="tico")
    xr=float(ref["NEB1_reaction_energy_EpoxideHHO_minus_Rpair"])
    zr=float(ref["OOH_star_legacy"])
    rp=float(ref["Rpair_strict"])
    coeff=(.41212309733267594,.3703773005809142,.20021749611103448,-.30977739698033324)
    points=[("baseline",xr,zr,True),
            ("new_BEP_anchor",xr,zr,False),
            ("middle", -1.2,-4.31,False)]
    for name,x,z,baseline in points:
        target=HERE/"native_catmap_checks"/name
        target.mkdir(parents=True,exist_ok=True)
        lines=(SOURCE/"alkene_epoxidation.mkm").read_text(encoding="utf-8")
        lines+="\ndata_file = 'native_check.pkl'\nlog_file = 'native_check.log'\n"
        (target/"alkene_epoxidation.mkm").write_text(lines,encoding="utf-8")
        if baseline:
            shutil.copy2(SOURCE/"energies.txt",target/"energies.txt")
        else:
            y=coeff[0]+coeff[1]*(x+2.0922766923076916)+coeff[2]*(z+4.309874384615384)+coeff[3]*(z+4.309874384615384)**2
            with (SOURCE/"energies.txt").open(encoding="utf-8",newline="") as f:
                reader=csv.DictReader(f,delimiter="\t");headers=reader.fieldnames;records=list(reader)
            counts={"TS1":0,"EpoxideHHO":0}
            for row in records:
                if row["surface_name"]=="tico" and row["site_name"]=="active" and row["species_name"] in counts:
                    sp=row["species_name"]
                    row["formation_energy"]=str(rp+(y if sp=="TS1" else x))
                    counts[sp]+=1
            assert counts=={"TS1":1,"EpoxideHHO":1},counts
            with (target/"energies.txt").open("w",encoding="utf-8",newline="") as f:
                w=csv.DictWriter(f,fieldnames=headers,delimiter="\t");w.writeheader();w.writerows(records)
        (target/"run_check.py").write_text(DRIVER,encoding="utf-8")
        r=subprocess.run([sys.executable,"-X","utf8","run_check.py"],cwd=target,encoding="utf-8",errors="replace",capture_output=True)
        (target/"console.log").write_text(r.stdout+"\n"+r.stderr,encoding="utf-8")
        if r.returncode:
            raise RuntimeError(f"CatMAP failed at {name}; see {target/'console.log'}")
        result=json.loads((target/"official_result.json").read_text(encoding="utf-8"))
        print(name,x,z,"CatMAP 0.3.1 TOF",result["tof"])
if __name__=="__main__":main()
