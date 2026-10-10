#!/usr/bin/env python3
"""Build coadsorption CatMAP energies from the ORIGINAL Excel workbook.
Requires: artifact_tool. No writes to the original workbook or old model.
Run: python generate_input1.py (from 1-generating_input_file).
"""
import csv
from pathlib import Path
from artifact_tool import Blob, SpreadsheetFile

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "己烯环氧化.xlsx"
OUTPUT = HERE.parent / "2-creating_microkinetic_model" / "energies.txt"
AUDIT = HERE / "energy_audit.csv"
# 与原版 data_only_not_executable.py 相同的原子参考能，单位 eV。
# 固定元素参考值不能直接作为火山图的横纵坐标。
ref_dict = {
    'O': -4.37,
    'C': -9.28,
    'H': -1.11,
}
# 15 个不同催化剂必须使用各自的裸板自由能作为位点参考。
STOICH = {
    "C6H12": {"C":6,"H":12}, "H2O2": {"H":2,"O":2},
    "C6H12O": {"C":6,"H":12,"O":1}, "H2O": {"H":2,"O":1},
    "R": {"C":6,"H":14,"O":2}, "P": {"C":6,"H":14,"O":2},
    "RP": {"C":6,"H":14,"O":2}
}

def atom_ref(formula):
    return sum(ref_dict[k] * v for k,v in STOICH[formula].items())

def n(value, what):
    if not isinstance(value,(float,int)):
        raise ValueError(f"Missing numeric cell: {what}: {value!r}")
    return float(value)

def main():
    book = SpreadsheetFile.import_xlsx(Blob.load(str(SOURCE)))
    ws = book.worksheets.get_item("吉布斯自由能汇总")
    values = ws.get_range("A1:DE81").values
    def cell(row, col): return values[row-1][col-1]
    surfaces = [(1,"tife"),(9,"timn"),(17,"tihf"),(24,"tire"),(32,"tinb"),
                (39,"timo"),(46,"tiv"),(53,"tizr"),(60,"tico"),(67,"titi"),
                (74,"tiw"),(81,"tita"),(89,"ticr"),(96,"ti"),(103,"tini")]
    gases = {name:n(cell(r,6), f"gas {name} G") for name,r in
             [("H2O2",57),("C6H12",58),("H2O",59),("C6H12O",60)]}
    output = []
    for name,g in gases.items():
        output.append(("None","gas",name,g-atom_ref(name)))
    rows=[]
    for col,surface in surfaces:
        gcol = col+6 if surface == 'tita' else col+5
        # Raw corrected free energies, from the SAME sheet:
        # slab row 4, coadsorbed reactant row 12, pathway marker row 16,
        # coadsorbed product row 20 (TiMo composite total: row 21).
        slab=n(cell(4,gcol), f"{surface} slab G")
        site_ref_dict = {**ref_dict, "111": slab}
        initial=n(cell(12,gcol), f"{surface} IS G")
        marker=n(cell(16,gcol), f"{surface} NEB-marker G")
        final=n(cell(21 if surface == "timo" else 20,gcol), f"{surface} FS G")
        delta_marker=marker-initial
        delta_final=final-initial
        # Rate-model effective barrier, not a claim about actual NEB saddle.
        # For downhill marker choose zero forward activation; for others
        # marker is used as a provisional candidate until full NEB validation.
        effective=max(initial, marker, final)
        gfs={name:g-site_ref_dict["111"]-atom_ref(name) for name,g in
             [("R",initial),("P",final),("RP",effective)]}
        for name,gf in gfs.items():
            output.append((surface,"111",name,gf))
        rows.append((surface,slab,initial,marker,final,delta_marker,
                     delta_final,effective-initial,
                     "UNVERIFIED_FULL_NEB",
                     "EFFECTIVE_TS_CLAMPED" if effective!=marker else "NEB_MARKER_CANDIDATE"))
        if abs((gfs["P"]-gfs["R"])-delta_final)>1e-7:
            raise AssertionError(surface+" inconsistent free-energy reference")
    with OUTPUT.open("w",encoding="utf-8",newline="") as fh:
        wr=csv.writer(fh,delimiter="\t",lineterminator="\n")
        wr.writerow(["surface_name","site_name","species_name","formation_energy","frequencies","reference"])
        for surface,site,name,gf in output:
            wr.writerow([surface,site,name,f"{gf:.9f}","[]",
                         "Original workbook Gibbs summary; 333.15 K"])
    with AUDIT.open("w",encoding="utf-8",newline="") as fh:
        wr=csv.writer(fh)
        wr.writerow(["surface","G_slab","G_IS","G_path_marker","G_FS",
                     "marker_minus_IS","FS_minus_IS","effective_barrier",
                     "full_NEB_status","kinetic_TS_status"])
        wr.writerows(rows)
    print(f"Wrote {len(output)} energy entries and {len(rows)} audit rows.")
    print("NOTE: effective barriers are MODEL ASSUMPTIONS until full NEB validation.")

if __name__=="__main__":
    main()
