#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot CatMAP 8-step activity map with C/O/H-screened O* and Ha* descriptors.

Run AFTER the updated test.ipynb exports production_rate_table.txt:
    python plot_C_O_H_volcano.py
Output: volcano_C_O_H.pdf. No CatMAP import is needed for this postprocessor.

E_C is C6H12 adsorption energy, E_O is O* formation energy,
E_H is Ha* formation energy. NONE is a fixed elemental reference mu.
"""
from __future__ import annotations
import ast
import csv
from datetime import datetime
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
TABLE = ROOT / "production_rate_table.txt"
DESC = ROOT / "descriptor_C_O_H.csv"
MKM = ROOT / "alkene_epoxidation.mkm"

def rows(p, delimiter=","):
    with p.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))

def model_surfaces():
    for node in ast.parse(MKM.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "surface_names" for t in node.targets
        ):
            return list(ast.literal_eval(node.value))
    raise ValueError("No surface_names in .mkm")

def main():
    desc = {r["surface"]: r for r in rows(DESC)}
    records = rows(TABLE, "\t")
    if not records:
        raise ValueError("Empty production_rate_table.txt")
    axes = ["descriptor-O_s", "descriptor-Ha_s"]
    if not all(x in records[0] for x in [*axes, "C6H12O_g"]):
        raise ValueError(
            "Stale three-step R/P activity table. Re-run the eight-step "
            "lyq_alkene_epoxidation1/test.ipynb first. "
            "Expected O_s and Ha_s descriptor columns."
        )
    pts = np.array([[float(r[col]) for col in axes] for r in records], dtype=float)
    rates = np.array([float(r["C6H12O_g"]) for r in records], dtype=float)
    if not np.all(np.isfinite(rates)):
        raise ValueError("Nonfinite CatMAP production rates; validate solver first")
    if np.any(rates < 0):
        factor = max(float(np.percentile(np.abs(rates), 75)), 1e-30)
        z = np.arcsinh(rates/factor)
        bar = "asinh(net epoxide production / scale)"
    else:
        z = np.log10(np.maximum(rates, 1e-40))
        bar = "log10(epoxide production rate, model units)"
    fig, ax = plt.subplots(figsize=(8.3, 6.4))
    filled = ax.tricontourf(pts[:,0], pts[:,1], z, levels=22)
    fig.colorbar(filled, ax=ax, label=bar)
    for surf in model_surfaces():
        if surf not in desc:
            raise ValueError("Missing C/O/H descriptors: "+surf)
        r=desc[surf]
        # E_O/E_H have been independently derived from original energies.txt
        x, y = float(r["E_O"]), float(r["E_H"])
        ax.scatter(x,y,s=32,c="white",edgecolors="black",zorder=3)
        ax.annotate(surf,(x,y),fontsize=8,xytext=(4,3),
                    textcoords="offset points",zorder=4)
    ax.set(xlabel="O* formation energy (eV)",
           ylabel="Ha* formation energy (eV)",
           title="Eight-step CatMAP epoxidation activity: C/O/H-screened descriptors")
    fig.tight_layout()
    path = ROOT/"volcano_C_O_H.pdf"
    try:
        fig.savefig(path)
    except PermissionError:
        dest=ROOT/"plots"/("rerun_"+datetime.now().strftime("%Y%m%d_%H%M%S"))
        dest.mkdir(parents=True,exist_ok=True)
        path=dest/path.name
        fig.savefig(path)
    finally:
        plt.close(fig)
    print("Saved",path)
    print("Model: 8-step CatMAP. Axes: O_s / Ha_s. Color: epoxide production rate.")
    print("Note: exploratory descriptor fit, not experimentally validated TOF.")

if __name__ == "__main__":
    main()
