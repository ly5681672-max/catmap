#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Select C/O/H-derived descriptors against NEB epoxidation barrier.

Independent analysis stage feeding the ORIGINAL 8-step CatMAP model.
This never edits DFT energies, energies.txt, or the .mkm configuration.
Run: python 1-generating_input_file/bep_descriptor_analysis.py
"""
from __future__ import annotations
import ast
import csv
from datetime import datetime
from itertools import combinations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent / "2-creating_microkinetic_model"
AUDIT = HERE / "energy_audit.csv"
ENERGIES = MODEL_DIR / "energies.txt"
MKM = MODEL_DIR / "alkene_epoxidation.mkm"

# Distinct choices within the requested C/O/H family.
# 'C' = 1-hexene adsorption relative to its gas-phase reference, not fixed mu_C.
# 'O' = O* formation energy, 'H' = H* formation energy (Ha/Hb positions).
CANDIDATES = ("C6H12_ads", "O_form", "Ha_form", "Hb_form")
EXPECTED_MODEL_DESCRIPTORS = ("O_s", "Ha_s")

def read_model_surfaces():
    exprs = ast.parse(MKM.read_text(encoding="utf-8"))
    for node in exprs.body:
        if isinstance(node, ast.Assign):
            if any(isinstance(t, ast.Name) and t.id == "surface_names" for t in node.targets):
                return tuple(ast.literal_eval(node.value))
    raise ValueError("surface_names missing from original-framework .mkm")

def read_rows():
    with ENERGIES.open(encoding="utf-8", newline="") as fh:
        erows = list(csv.DictReader(fh, delimiter="\t"))
    e, gas = {}, {}
    for row in erows:
        if row["surface_name"] == "None":
            gas[row["species_name"]] = float(row["formation_energy"])
        else:
            e.setdefault(row["surface_name"], {})[row["species_name"]] = float(row["formation_energy"])
    with AUDIT.open(encoding="utf-8", newline="") as fh:
        audit = {r["surface"]: r for r in csv.DictReader(fh)}
    result = []
    for surface in read_model_surfaces():
        item = e[surface]
        g = audit[surface]
        raw = float(g["marker_minus_IS"])
        dFS = float(g["FS_minus_IS"])
        # The three-point net-NEB curve does not contain a verified saddle
        # for all surfaces. Keep both signed marker and provisional effective Ea.
        Ea = float(g["effective_barrier"])
        if abs(Ea - max(0.0, raw, dFS)) > 5e-6:
            raise ValueError("Unexpected activation-barrier definition: " + surface)
        result.append(dict(
            surface=surface,
            C6H12_ads=item["C6H12"] - gas["C6H12"],
            O_form=item["O"],
            Ha_form=item["Ha"],
            Hb_form=item["Hb"],
            OH_form=item["OH"],
            OOH_form=item["OOH"],
            NEB_marker_minus_IS=raw,
            reaction_G=dFS,
            Ea_provisional=Ea,
            downward_marker=int(raw < 0),
        ))
    return result

def fit_ols(rows, names):
    x = np.asarray([[1.0] + [r[f] for f in names] for r in rows], dtype=float)
    y = np.asarray([r["Ea_provisional"] for r in rows], dtype=float)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    pred = x @ beta
    ss = np.sum((y - pred) ** 2)
    denom = np.sum((y - y.mean()) ** 2)
    r2 = float(1 - ss / denom) if denom > 1e-15 else float("nan")
    return beta, r2, pred

def loocv(rows, names):
    pred = np.asarray([
        float(np.dot([1.0]+[r[f] for f in names],
                     fit_ols([v for j,v in enumerate(rows) if j != i], names)[0]))
        for i,r in enumerate(rows)
    ])
    actual = np.asarray([r["Ea_provisional"] for r in rows])
    return float(np.sqrt(np.mean((pred - actual)**2))), pred

def save_csv(filename, records, columns):
    with (MODEL_DIR / filename).open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(records)

def save_plot(fig, filename):
    target = MODEL_DIR / filename
    try:
        fig.savefig(target)
    except PermissionError:
        folder = MODEL_DIR / "plots" / ("rerun_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
        folder.mkdir(parents=True, exist_ok=True)
        target = folder / filename
        fig.savefig(target)
    print("Saved", target)

def run(make_plots=True):
    rows = read_rows()
    reports=[]
    # Screening includes single-descriptor linear fits and all pairs.
    sets = [(x,) for x in CANDIDATES] + list(combinations(CANDIDATES, 2))
    for names in sets:
        beta, r2, _ = fit_ols(rows, names)
        rmse, _ = loocv(rows, names)
        reports.append(dict(
            descriptors="+".join(names), count=len(rows),
            R2=r2, LOOCV_RMSE_eV=rmse,
            intercept_eV=beta[0], coef_1=beta[1],
            coef_2=beta[2] if len(beta)>2 else "",
            note="linear association of provisional net-NEB barrier; not elementary-step BEP",
        ))
    reports.sort(key=lambda r:r["LOOCV_RMSE_eV"])
    save_csv("bep_screening.csv", reports, list(reports[0]))
    pairs=[r for r in reports if "+" in r["descriptors"]]
    best = pairs[0]
    names=tuple(best["descriptors"].split("+"))
    beta, r2, pred=fit_ols(rows,names)
    cv, cvpred=loocv(rows,names)
    data=[]
    for i,r in enumerate(rows):
        data.append({**r,"E_barrier_fit_eV":float(pred[i]),
                     "E_barrier_LOOCV_eV":float(cvpred[i]),
                     "fit_residual_eV":float(r["Ea_provisional"]-pred[i])})
    save_csv("bep_predictions.csv",data,list(data[0]))
    # Inspect strongest unrestricted pair separately; never silently broaden C/O/H choice.
    extnames=("OH_form","Ha_form")
    _, extR2,_ =fit_ols(rows,extnames)
    extcv,_=loocv(rows,extnames)
    print("BEP-like linear screening: 14 surfaces matching original .mkm")
    print("C/O/H-family best pair:", names, "R2 =",round(r2,4),"LOOCV RMSE =",round(cv,4),"eV")
    print("Coefficients [intercept, x, y] =", [round(float(v),6) for v in beta])
    print("Extended (OH*,Ha*) diagnostic: R2 =",round(extR2,4),
          "LOOCV RMSE =",round(extcv,4),"eV")
    print("This is NOT classical elementary-step BEP and does NOT validate full-NEB saddles.")
    if names != ("O_form","Ha_form"):
        print("WARNING: selected pair changed. Revisit .mkm descriptor_names before plotting volcano.")
    if make_plots:
        import matplotlib.pyplot as plt
        import matplotlib
        actual=np.asarray([r["Ea_provisional"] for r in rows])
        fig,ax=plt.subplots(figsize=(6.4,5.4))
        ax.scatter(actual, pred, label="training fit")
        lo=float(min(actual.min(),pred.min(),cvpred.min())-.05)
        hi=float(max(actual.max(),pred.max(),cvpred.max())+.05)
        ax.plot([lo,hi],[lo,hi],linestyle="--",color="gray",label="y=x")
        for r,x,y in zip(rows,actual,pred):
            ax.annotate(r["surface"],(x,y),fontsize=7,xytext=(4,3),textcoords="offset points")
        ax.set(xlabel="NEB-derived provisional forward barrier (eV)",
               ylabel="Linear-model prediction (eV)",title=f"C/O/H screening: R²={r2:.3f}; LOO RMSE={cv:.3f} eV")
        ax.legend()
        fig.tight_layout()
        save_plot(fig, "bep_fit.pdf")
        plt.close(fig)
        # Plane of activation barriers, NOT a turnover-frequency volcano.
        x=np.array([r[names[0]] for r in rows]);y=np.array([r[names[1]] for r in rows])
        xx,yy=np.meshgrid(np.linspace(x.min()-.15,x.max()+.15,85),
                          np.linspace(y.min()-.15,y.max()+.15,85))
        z=np.maximum(0., beta[0]+beta[1]*xx+beta[2]*yy)
        fig,ax=plt.subplots(figsize=(7,5.6))
        cf=ax.contourf(xx,yy,z,levels=20)
        fig.colorbar(cf,ax=ax,label="fitted effective forward barrier (eV)")
        ax.scatter(x,y,edgecolors="black",facecolors="white",s=35)
        for r,a,b in zip(rows,x,y):
            ax.annotate(r["surface"],(a,b),fontsize=8,xytext=(4,3),textcoords="offset points")
        ax.set(xlabel=names[0]+" (eV)",ylabel=names[1]+" (eV)",
               title="Linear BEP-like descriptor plane (NOT TOF volcano)")
        fig.tight_layout()
        save_plot(fig, "bep_barrier_map.pdf")
        plt.close(fig)
        print("Saved bep_fit.pdf, bep_barrier_map.pdf")
    return rows,best

if __name__ == "__main__":
    run()
