"""Prepare one isolated provisional site model from project source tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROJECT = ROOT.parents[2]
SNAPSHOT = ROOT / "inputs" / "source_snapshot"
RUNS = ROOT / "runs"
SOURCE_DIR = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model")
SOURCE_PROJECT = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation")
TICO_TEMPLATE = PROJECT / "analyses" / "explicit_OOH_8step_validation" / "candidate_tico_kinetic_control" / "candidate_tico.mkm"
FILES = {
    "alkene_epoxidation.mkm": SOURCE_DIR / "alkene_epoxidation.mkm",
    "energies_dft_neb.txt": SOURCE_DIR / "energies_dft_neb.txt",
    "DFT_formation_catmap.txt": SOURCE_DIR / "DFT_formation_catmap.txt",
    "DFT_TS_relative.csv": SOURCE_DIR / "DFT_TS_relative.csv",
    "DFT_adsorption.csv": SOURCE_DIR / "DFT_adsorption.csv",
    "convert_dft_to_catmap.py": SOURCE_PROJECT / "1-generating_input_file" / "convert_dft_to_catmap.py",
    "source_workbook.xlsx": SOURCE_PROJECT / "己烯环氧化.xlsx",
    "candidate_network_template.mkm": TICO_TEMPLATE,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def snapshot_sources() -> dict[str, str]:
    SNAPSHOT.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name, source in FILES.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        dest = SNAPSHOT / name
        if dest.exists() and sha256(dest) != sha256(source):
            raise RuntimeError(f"refusing to overwrite a different snapshot: {dest}")
        if not dest.exists():
            shutil.copy2(source, dest)
        hashes[str(source)] = sha256(source)
    (SNAPSHOT / "source_sha256.json").write_text(json.dumps(hashes, indent=2), encoding="utf-8")
    return hashes


def load_converter():
    path = SNAPSHOT / "convert_dft_to_catmap.py"
    spec = importlib.util.spec_from_file_location("project_conversion_rules", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load conversion rules: {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
        return list(csv.DictReader(f, delimiter=delimiter))


def prepare(surface: str) -> Path:
    hashes = snapshot_sources()
    converter = load_converter()
    source_data = converter.load_sheets(SNAPSHOT / "source_workbook.xlsx")
    single, summary = source_data["catmap数据集"], source_data["吉布斯自由能汇总"]
    single_cols = converter.surfaces(single)
    summary_cols = converter.surfaces(summary)
    if surface not in single_cols or surface not in summary_cols:
        raise ValueError(f"surface absent from source workbook: {surface}")
    single_g_col = single_cols[surface][1]
    summary_g_col = summary_cols[surface][1]
    clean_slab = float(single[converter.cell(single_g_col, 51)])

    neb_row = next((r for r in read_table(SNAPSHOT / "DFT_TS_relative.csv") if r["surface"] == surface), None)
    if neb_row is None:
        raise ValueError(f"no complete NEB row for {surface}")
    ref1 = converter.ref_energy({"C": 6, "H": 14, "O": 2})
    ref2 = converter.ref_energy({"H": 2, "O": 1})
    neb_states = {
        "NEB1IS": float(neb_row["G_IS1"]) - clean_slab - ref1,
        "NEB1TS": float(neb_row["G_TS1"]) - clean_slab - ref1,
        "NEB1FS": float(neb_row["G_FS1_total"]) - clean_slab - ref1,
        "NEB2IS": float(neb_row["G_IS2"]) - clean_slab - ref2,
        "NEB2TS": float(neb_row["G_TS2"]) - clean_slab - ref2,
        "NEB2FS": float(neb_row["G_FS2_adjusted"]) - clean_slab - ref2,
    }

    source_rows = read_table(SNAPSHOT / "DFT_formation_catmap.txt")
    gas_rows = [r for r in source_rows if r["site_name"] == "gas"]
    site_rows = [r for r in source_rows if r["surface_name"] == surface and r["site_name"] != "gas"]
    required = {"H2O2", "Ha", "OOH", "C6H12", "H2O", "C6H12O"}
    available = {r["species_name"] for r in site_rows}
    if not required.issubset(available):
        raise ValueError(f"missing independent states for {surface}: {sorted(required-available)}")
    selected = gas_rows + site_rows + [
        {"surface_name": surface, "site_name": "111", "species_name": species,
         "formation_energy": f"{energy:.12f}", "frequencies": "[]",
         "reference": "FULL_NEB_STATE_FORMATION_FROM_SOURCE_G_MINUS_SURFACE_G_MINUS_ELEMENTAL_REFERENCE"}
        for species, energy in neb_states.items()
    ]

    run_dir = RUNS / f"{surface}_candidate"
    run_dir.mkdir(parents=True, exist_ok=True)
    energy_file = run_dir / f"candidate_{surface}_energies.txt"
    fields = ["surface_name", "site_name", "species_name", "formation_energy", "frequencies", "reference"]
    with energy_file.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows({k: r[k] for k in fields} for r in selected)

    text = (SNAPSHOT / "candidate_network_template.mkm").read_text(encoding="utf-8")
    text, n_surface = re.subn(r"surface_names\s*=\s*\['tico'\]", f"surface_names = ['{surface}']", text)
    text, n_data = re.subn(r"input_file\s*=\s*'candidate_energies\.txt'", f"input_file = '{energy_file.name}'", text)
    text, n_pickle = re.subn(r"data_file\s*=\s*'candidate_tico_endpoint_model\.pkl'", f"data_file = 'candidate_{surface}.pkl'", text)
    if (n_surface, n_data, n_pickle) != (1, 1, 1):
        raise RuntimeError(f"template anchors not unique: surface={n_surface}, input={n_data}, data={n_pickle}")
    model_file = run_dir / f"candidate_{surface}.mkm"
    model_file.write_text(text, encoding="utf-8")

    neb_deltas = {
        "NEB1_TS_minus_IS_eV": neb_states["NEB1TS"] - neb_states["NEB1IS"],
        "NEB1_TS_minus_FS_eV": neb_states["NEB1TS"] - neb_states["NEB1FS"],
        "NEB2_TS_minus_IS_eV": neb_states["NEB2TS"] - neb_states["NEB2IS"],
        "NEB2_TS_minus_FS_eV": neb_states["NEB2TS"] - neb_states["NEB2FS"],
    }
    metadata = {
        "surface": surface, "candidate_status": "PROVISIONAL",
        "source_hashes": hashes,
        "clean_slab_G_source_sheet": "catmap数据集", "clean_slab_G_source_row": 51,
        "clean_slab_G_eV": clean_slab,
        "NEB_state_source": "DFT_TS_relative.csv; raw Gibbs totals; same surface-specific clean slab and converter atomic references",
        "NEB_formation_eV": neb_states, "NEB_signed_differences_eV": neb_deltas,
        "conditions": {"temperature_K": 333.15, "concentrations": {"C6H12_g": 1.0, "H2O2_g": 0.185, "H2O_g": 0.815, "C6H12O_g": 1e-20}},
        "site_model": "one effective site per NEB composite state; coarse-graining, not structure-verified",
        "missing_link_TS": "CatMAP endpoint-barrier approximation; not DFT TS",
        "model_mapping_warning": "reduced 9-reaction candidate, not a strict lumping proof of the original 8-reaction model",
        "candidate_model": str(model_file), "candidate_energies": str(energy_file),
    }
    (run_dir / "preparation_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    with (run_dir / "reaction_energy_barriers.csv").open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction", "value_eV", "meaning"])
        w.writeheader()
        for key, value in neb_deltas.items():
            w.writerow({"reaction": key, "value_eV": f"{value:.12g}", "meaning": "signed NEB TS minus endpoint; not replaced by max(0, value)"})
    print(json.dumps(metadata, ensure_ascii=False, indent=2))
    return model_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--surface", required=True, help="lowercase surface_name")
    args = parser.parse_args()
    prepare(args.surface)
