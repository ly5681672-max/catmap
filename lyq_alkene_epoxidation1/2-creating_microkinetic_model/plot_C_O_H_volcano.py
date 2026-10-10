#!/usr/bin/env python3
"""COH descriptor activity projections (not a reparameterized CatMAP C/O/H model).

Requires numerical E_C,E_O,E_H values per surface in descriptor_C_O_H.csv.
Activity at R/P is interpolated on the *existing* 2D R/P CatMAP map, then
projected onto C/O, C/H, and O/H descriptor scatter plots.  This MUST NOT be
misinterpreted as a genuine 3-descriptor microkinetic fit.
"""
from __future__ import annotations
import csv
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent
CSV = ROOT / "descriptor_C_O_H.csv"
RATE = ROOT / "production_rate_table.txt"
AUDIT = ROOT / "ts_scaling_audit.csv"
PLOT = ROOT / "volcano_C_O_H.pdf"

def rows(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))

def main():
    if not CSV.is_file():
        raise FileNotFoundError(f"Missing {CSV}. Values E_C/E_O/E_H must be supplied.")
    descriptors = rows(CSV)
    provided = {}
    for r in descriptors:
        s = r["surface"].lower().strip()
        try:
            vals = tuple(float(r[key]) for key in ("E_C", "E_O", "E_H"))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{s}: missing E_C/E_O/E_H; supply the validated "
                             "surface adsorption/descriptor data first.") from exc
        if not all(np.isfinite(vals)):
            raise ValueError(f"{s}: descriptor contains nonfinite values")
        if s in provided:
            raise ValueError(f"Duplicate surface {s}")
        provided[s] = vals
    records = rows(AUDIT)
    missing = [r["surface"] for r in records if r["surface"] not in provided]
    if missing:
        raise ValueError("Missing descriptors for: " + ", ".join(missing))
    production = rows(RATE, delimiter="\t")
    if not production or "C6H12O_g" not in production[0]:
        raise ValueError("production_rate_table.txt must contain C6H12O_g. "
                         "Run notebook's model and export cells first.")
    grid_xy = np.asarray([[float(r["descriptor-R_s"]),float(r["descriptor-P_s"])]
                          for r in production], dtype=float)
    grid_log_rate = np.asarray([
        np.log10(max(float(r["C6H12O_g"]), 1e-30)) for r in production
    ], dtype=float)
    names = []
    vals = []
    rate_values = []
    for record in records:
        name = record["surface"]
        xy = [float(record["R_input"]),float(record["P_input"])]
        predicted = griddata(grid_xy, grid_log_rate, [xy], method="linear")[0]
        if not np.isfinite(predicted):
            raise ValueError(f"{name}: R/P outside numerical descriptor grid; "
                             "expand descriptor_ranges and recompute.")
        names.append(name)
        vals.append(provided[name])
        rate_values.append(predicted)
    xyz = np.asarray(vals,dtype=float)
    rate_values = np.asarray(rate_values,dtype=float)
    # With 3 input descriptors, report pairwise 2D maps rather than invent
    # a fourth dimension or silently hold the third at arbitrary values.
    pages = []
    for xidx,yidx in [(0,1),(0,2),(1,2)]:
        x,y=xyz[:,xidx],xyz[:,yidx]
        fig,ax=plt.subplots(figsize=(7.4,5.6))
        if (np.ptp(x)>1e-9 and np.ptp(y)>1e-9):
            try:
                cn=ax.tricontourf(x,y,rate_values,levels=14,cmap="viridis")
                fig.colorbar(cn,ax=ax,label="log10(epoxide production, R/P-interpolated)")
            except (ValueError,RuntimeError):
                sc=ax.scatter(x,y,c=rate_values,cmap="viridis",s=65)
                fig.colorbar(sc,ax=ax,label="log10(epoxide production, R/P-interpolated)")
        else:
            sc=ax.scatter(x,y,c=rate_values,cmap="viridis",s=65)
            fig.colorbar(sc,ax=ax,label="log10(epoxide production, R/P-interpolated)")
        ax.scatter(x,y,c="white",s=14,edgecolor="black",linewidth=.4,zorder=4)
        for name,xx,yy in zip(names,x,y):
            ax.annotate(name,(xx,yy),xytext=(3,3),textcoords="offset points",fontsize=7)
        ax.set(xlabel=["E_C","E_O","E_H"][xidx]+" [eV]",
               ylabel=["E_C","E_O","E_H"][yidx]+" [eV]",
               title="C/O/H descriptor activity projection (exploratory)")
        fig.tight_layout()
        pages.append(fig)
    target = PLOT
    try:
        with PdfPages(target) as pdf:
            for fig in pages:
                pdf.savefig(fig,bbox_inches="tight")
    except PermissionError:
        out=ROOT/"plots"/("rerun_"+datetime.now().strftime("%Y%m%d_%H%M%S"))
        out.mkdir(parents=True,exist_ok=True)
        target=out/PLOT.name
        with PdfPages(target) as pdf:
            for fig in pages:
                pdf.savefig(fig,bbox_inches="tight")
    finally:
        for fig in pages:
            plt.close(fig)
    print("Saved:",target)
    print("NOTE: plotted activity is interpolated from the R/P CatMAP map.")
    print("This is NOT a true C/O/H-descriptor reaction scaling model.")

if __name__=="__main__":
    main()
