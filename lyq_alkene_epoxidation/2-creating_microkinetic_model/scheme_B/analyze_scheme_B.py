"""Scheme B diagnostics: 1% flux audit and independent catalyst-point solving.

Does not alter DFT energies or the eight reaction steps. Only unpickle data
produced by your own trusted CatMAP installation.
"""
import argparse
import ast
import csv
import math
import pickle
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STATES = ("S1","S2","S3","S4","S5","S6","S7")

def settings():
    tree=ast.parse((ROOT/"alkene_epoxidation.mkm").read_text(encoding="utf-8"))
    need={"descriptor_ranges","resolution","surface_names"}
    result={}
    for node in tree.body:
        if isinstance(node,ast.Assign) and len(node.targets)==1:
            target=node.targets[0]
            if isinstance(target,ast.Name) and target.id in need:
                result[target.id]=ast.literal_eval(node.value)
    if set(result)!=need:
        raise ValueError("Could not read descriptor grid and surfaces from .mkm")
    return result

def read_input():
    with (ROOT/"energy_audit.csv").open(encoding="utf-8-sig",newline="") as file:
        return {r["surface"]:r for r in csv.DictReader(file)}

def f(x):
    try: return float(x)
    except (ValueError,TypeError,OverflowError): return float("nan")

def xy(pt):
    return tuple(round(float(x),8) for x in pt)

def get_grid(ranges,resolution):
    if isinstance(resolution,int): resolution=[resolution,resolution]
    axes=[[lo+(hi-lo)*i/(n-1) for i in range(n)] for (lo,hi),n in zip(ranges,resolution)]
    return {xy((x,y)) for x in axes[0] for y in axes[1]}

def saved_maps(data):
    out={}
    for name in ("coverage_map","rate_map","production_rate_map","numbers_map"):
        if not data.get(name): raise ValueError("Missing "+name+"; rerun CatMAP.")
        rows={}
        for pt,item in data[name]:
            key=xy(pt)
            if key in rows: raise ValueError("Repeated descriptor point in "+name)
            rows[key]=item
        out[name]=rows
    return out

def writecsv(path,rows):
    if not rows: return
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="",encoding="utf-8-sig") as file:
        writer=csv.DictWriter(file,fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

def energy_stats(E):
    out={}
    for ts,initial,final in (("TS1","S3","S4"),("TS2","S5","S6")):
        g0,gts,g1=(f(E[k]) for k in (initial,ts,final))
        top=max(g0,gts,g1)
        out[ts+"_minus_IS_eV"]=gts-g0
        out[ts+"_minus_FS_eV"]=gts-g1
        out[ts+"_effective_Eaf_eV"]=top-g0
        out[ts+"_effective_Ear_eV"]=top-g1
        out[ts+"_effective_top"]=("TS" if gts>=max(g0,g1) else
                                    "IS" if g0>=g1 else "FS")
    return out

def scaler_values(model,pt):
    params=model.scaler.get_rxn_parameters(list(pt))
    names=list(model.adsorbate_names)+list(model.transition_state_names)
    if len(params)<len(names): raise ValueError("Invalid scaler parameter size")
    return {name.split("_")[0]:f(params[i]) for i,name in enumerate(names)},params

def diagnostics(pt,coverages,rates,product=None,tolerance=.01,min_signal=0):
    cover=[f(x) for x in coverages]
    rate=[f(x) for x in rates]
    if len(cover)!=8 or len(rate)!=8: raise ValueError("Expected eight states and eight reactions")
    ok=all(math.isfinite(x) for x in cover+rate)
    maxabs=max(abs(x) for x in rate) if ok else float("nan")
    absdiff=max(abs(x-rate[4]) for x in rate) if ok else float("nan")
    rel=absdiff/max(maxabs,1e-300) if ok else float("nan")
    cs=sum(cover)
    cover_ok=ok and min(cover)>=-1e-8 and max(cover)<=1+1e-8 and abs(cs-1)<1e-6
    flux_ok=ok and rel<=tolerance
    gas_rate=f(product[0]) if product is not None else float("nan")
    product_ok=rate[4]>0 and (not math.isfinite(gas_rate) or gas_rate>0)
    signal_ok=maxabs>=min_signal
    qualified=cover_ok and flux_ok and product_ok and signal_ok
    if not ok: status="NONFINITE"
    elif not cover_ok: status="COVERAGE_INVALID"
    elif not flux_ok: status="FLUX_MISMATCH"
    elif not product_ok: status="NONPOSITIVE_PRODUCT"
    elif not signal_ok: status="LOW_SIGNAL"
    else: status="PASS"
    out={"descriptor_S3_eV":float(pt[0]),"descriptor_S5_eV":float(pt[1]),
         "status":status,"volcano_eligible":int(qualified),
         "flux_mismatch_relative":rel,"flux_mismatch_absolute":absdiff,
         "max_abs_step_rate":maxabs,"epoxide_net_step5":rate[4],
         "water_net_step8":rate[7],"epoxide_production_map":gas_rate,
         "coverage_sum":cs,"theta_empty":cover[7],
         "dominant_state":(STATES+("empty",))[max(range(8),key=cover.__getitem__)],
         "dominant_coverage":max(cover)}
    out.update({"net_r"+str(i+1):r for i,r in enumerate(rate)})
    out.update({"theta_"+s:cover[i] for i,s in enumerate(STATES)})
    return out

def solve_catalysts(model,data,tolerance,min_signal):
    """Solve the 14 measured catalyst coordinates, never a nearest-map proxy."""
    source=read_input()
    conf=settings()
    ranges=conf["descriptor_ranges"]
    seeds=data["numbers_map"] if model.use_numbers_solver else data["coverage_map"]
    model.solver.compile()
    result=[]
    for surface in conf["surface_names"]:
        src=source[surface]
        pt=(float(src["S3"]),float(src["S5"]))
        energy,params=scaler_values(model,pt)
        row={"surface":surface,"exact_point_solver":1}
        row.update(energy_stats(energy))
        for state in STATES+("TS1","TS2"):
            row["fit_minus_DFT_"+state+"_eV"]=energy[state]-float(src[state])
        options=[]
        for oldpt,seed in seeds:
            distance=math.sqrt(sum(((float(a)-float(b))/(hi-lo))**2
                        for a,b,(lo,hi) in zip(oldpt,pt,ranges)))
            options.append((distance,seed))
        options.sort(key=lambda v:v[0])
        errors=[]
        for distance,seed in options[:8]:
            try:
                # Explicit initial guess: exact nonlinear solve, never rounded-map retrieval.
                cov=list(model.solver.get_coverage(list(params),c0=seed))
                rates=list(model.solver.get_rate(list(params),coverages=cov,
                                                 verify_coverages=False))
                row.update(diagnostics(pt,cov,rates,tolerance=tolerance,min_signal=min_signal))
                row["solve_status"]="SOLVED"
                row["seed_distance_normalized"]=distance
                break
            except (ArithmeticError,ValueError,TypeError,OverflowError) as exc:
                errors.append(type(exc).__name__+": "+str(exc)[:100])
        else:
            row.update({"solve_status":"FAILED","status":"EXACT_SOLVE_FAILED",
                        "volcano_eligible":0,"error":" | ".join(errors[:3])})
        result.append(row)
    return result

def analyze(model=None,tolerance=.01,min_signal=0):
    dest=ROOT/"diagnostics"
    dest.mkdir(parents=True,exist_ok=True)
    cfg=settings()
    regular=get_grid(cfg["descriptor_ranges"],cfg["resolution"])
    with (ROOT/"scheme_B.pkl").open("rb") as file:
        data=pickle.load(file)  # Only use the locally generated, trusted CatMAP pickle.
    maps=saved_maps(data)
    points=set(maps["rate_map"])
    for name,series in maps.items():
        if regular-set(series): raise ValueError("Missing regular grid points from "+name)
    grid=[];aux=[]
    for pt in sorted(points):
        row=diagnostics(pt,maps["coverage_map"][pt],maps["rate_map"][pt],
                        maps["production_rate_map"][pt],tolerance,min_signal)
        if model is not None:
            energies,_=scaler_values(model,pt)
            row.update(energy_stats(energies))
        (grid if pt in regular else aux).append(row)
    writecsv(dest/"grid400_diagnostics.csv",grid)
    ranked=sorted((r for r in grid if r["volcano_eligible"]),
                  key=lambda r:r["epoxide_net_step5"],reverse=True)
    writecsv(dest/"grid400_volcano_eligible.csv",ranked)
    writecsv(dest/"bisection_extra_points.csv",aux)
    source=read_input()
    dft=[]
    for surface in cfg["surface_names"]:
        item=source[surface]
        x,y,z,w=(float(item[k]) for k in
                   ("TS1_minus_IS","TS1_minus_FS","TS2_minus_IS","TS2_minus_FS"))
        dft.append({"surface":surface,"S3_eV":item["S3"],"S5_eV":item["S5"],
                    "DFT_TS1_minus_IS_eV":x,"DFT_TS1_minus_FS_eV":y,
                    "DFT_TS2_minus_IS_eV":z,"DFT_TS2_minus_FS_eV":w,
                    "DFT_TS1_effective_Eaf_eV":max(0,x,x-y),
                    "DFT_TS2_effective_Eaf_eV":max(0,z,z-w)})
    writecsv(dest/"catalysts_input_barriers.csv",dft)
    summary=[{"metric":"regular_grid_points","value":len(grid)},
             {"metric":"extra_bisection_points","value":len(aux)},
             {"metric":"eligible_grid_points","value":len(ranked)},
             {"metric":"excluded_grid_points","value":len(grid)-len(ranked)},
             {"metric":"flux_tolerance","value":tolerance},
             {"metric":"min_signal_floor","value":min_signal}]
    if model is not None:
        summary += [
            {"metric":"grid_TS1_below_FS","value":sum(r["TS1_minus_FS_eV"]<0 for r in grid)},
            {"metric":"grid_TS2_below_FS","value":sum(r["TS2_minus_FS_eV"]<0 for r in grid)}]
        exact=solve_catalysts(model,data,tolerance,min_signal)
        writecsv(dest/"catalysts_exact.csv",exact)
        summary += [
            {"metric":"exact_catalysts_solved",
             "value":sum(r["solve_status"]=="SOLVED" for r in exact)},
            {"metric":"exact_catalysts_eligible",
             "value":sum(r["volcano_eligible"] for r in exact)}]
    else:
        summary.append({"metric":"exact_catalysts_solved","value":"grid_only"})
    writecsv(dest/"summary.csv",summary)
    print("Scheme B grid:",len(grid),"regular nodes;",
          len(ranked),"eligible;",len(grid)-len(ranked),"excluded.")
    if model is not None:
        print("Exact catalyst points solved:",summary[-2]["value"],"/",
              len(cfg["surface_names"]))
    print("Outputs:",dest)
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--grid-only",action="store_true",
                   help="Optional: grid-only is now the default for direct execution")
    p.add_argument("--flux-tol",type=float,default=.01)
    p.add_argument("--min-signal",type=float,default=0)
    args=p.parse_args()
    if not (0<args.flux_tol<1) or args.min_signal<0:
        p.error("flux-tol must be 0..1 and min-signal >=0")
    # Direct execution from PyCharm defaults to saved-grid diagnostics.
    # The --grid-only option is retained for backward compatibility.
    print("Grid-only mode: examining saved scheme_B.pkl; exact catalyst solves are NOT performed.")
    print("For 14 independent catalyst-point solves, execute run_scheme_B.py.")
    analyze(tolerance=args.flux_tol,min_signal=args.min_signal)
