#!/usr/bin/env python3
"""Standalone entry point for the exploratory CatMAP coadsorption model.

Run within this folder or from anywhere: python run_model.py
Options: --validate-only, --rebuild (requires artifact_tool), --allow-provisional
This deliberately does not report physically validated TOFs.
"""
from __future__ import annotations
import argparse
import csv
import os
from pathlib import Path
import runpy
import sys

ROOT = Path(__file__).resolve().parent
SURFACES = ("tife","timn","tihf","tire","tinb","timo","tiv","tizr",
            "tico","titi","tiw","tita","ticr","ti","tini")

def validate():
    from collections import defaultdict
    e=defaultdict(dict)
    p=ROOT/"energies.txt"
    with p.open(encoding="utf-8",newline="") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            e[row["surface_name"]][row["species_name"]]=float(row["formation_energy"])
    a=list(csv.DictReader((ROOT/"energy_audit.csv").open(encoding="utf-8",newline="")))
    assert len(a)==len(SURFACES),f"Expected {len(SURFACES)} surfaces"
    gases=e["None"]
    assert set(gases)>= {"H2O2","C6H12","H2O","C6H12O"}
    for row in a:
        s=row["surface"]
        assert s in SURFACES
        assert {"R","P","RP"}<=e[s].keys(),s
        d=e[s]
        delta=d["P"]-d["R"]
        actual=float(row["FS_minus_IS"])
        assert abs(delta-actual)<2e-6,(s,delta,actual)
        expected_ts=max(0.0,float(row["marker_minus_IS"]),actual)
        actual_ts=d["RP"]-d["R"]
        assert abs(actual_ts-expected_ts)<2e-6,(s,actual_ts,expected_ts)
    assert abs((e["tiw"]["P"]-e["tiw"]["R"])+2.675784)<2e-6
    print(f"PASS: {len(a)} surfaces, {len(gases)} gas references, formation-G mapping, provisional barriers.")
    return a

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild",action="store_true",help="Regenerate from bundled xlsx (requires artifact_tool)")
    parser.add_argument("--validate-only",action="store_true",help="Validate input and exit; CatMAP not required")
    parser.add_argument("--allow-provisional",action="store_true",help="Acknowledge kinetic barriers not validated against full NEB images")
    args=parser.parse_args()
    os.chdir(ROOT)
    if args.rebuild:
        runpy.run_path(str(ROOT/"build_coadsorption.py"),run_name="__main__")
    validate()
    if args.validate_only:
        return
    if not args.allow_provisional:
        parser.error("Kinetic barriers are provisional. Re-run with --allow-provisional to explicitly acknowledge this.")
    try:
        from catmap import ReactionModel
    except ImportError as exc:
        raise SystemExit("CatMAP is not importable. Activate your CatMAP environment and try again. Details: "+str(exc))
    model=ReactionModel(setup_file=str(ROOT/"coadsorption_model.mkm"))
    model.run()
    print("CatMAP run completed. Results are exploratory, NOT experimentally validated TOFs.")
    print("Output files remain in "+str(ROOT))

if __name__=="__main__":
    main()
