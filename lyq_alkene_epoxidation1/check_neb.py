#!/usr/bin/env python3
"""Audit source three-state free energies against original CatMAP step 4.
Does not modify any source files. Python standard library only.
"""
from __future__ import annotations
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "lyq_alkene_epoxidation" / "2-creating_microkinetic_model" / "energies.txt"


def load_energies(path: Path) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            surface = row["surface_name"].lower()
            if surface == "none":
                continue
            out.setdefault(surface, {})[row["species_name"]] = float(row["formation_energy"])
    return out


def main() -> None:
    with (HERE / "neb_three_state_summary.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as fh:
        source = list(csv.DictReader(fh))
    old = load_energies(OLD)
    print(
        "surface | Excel TS-IS | step4 TS-IS | Excel FS-IS | step4 FS-IS | "
        "difference in FS definitions | 3-point trend"
    )
    print("-" * 118)
    for row in source:
        surface = row["surface"]
        g_marker = float(row["deltaG_marker_eV"])
        g_fs = float(row["deltaG_FS_eV"])
        e = old[surface]
        old_is = e["OOH"] + e["C6H12"]
        old_ts = e["OOH-C6H12"] - old_is
        old_fs = e["C6H12O"] + e["O"] + e["Hb"] - old_is
        trend = "downhill 3 points" if 0 >= g_marker >= g_fs else "not downhill"
        print(
            f"{surface:6} | {g_marker:+11.3f} | {old_ts:+11.3f} | "
            f"{g_fs:+11.3f} | {old_fs:+11.3f} | "
            f"{old_fs-g_fs:+16.3f} | {trend}"
        )
        if abs(g_marker - old_ts) > 0.005:
            raise ValueError(f"{surface}: TS reference mismatch >5 meV")
    assert len(source) == 15
    print("\nPASS: all 15 nominal TS relative energies match old CatMAP within 5 meV.")
    print("NOTE: Excel FS is the FULL coadsorbed product state; step4 FS is not.")
    print("NOTE: a three-point downhill order is NOT a check of every NEB image.")


if __name__ == "__main__":
    main()
