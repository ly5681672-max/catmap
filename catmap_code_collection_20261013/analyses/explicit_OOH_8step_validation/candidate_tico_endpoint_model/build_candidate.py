"""Build an isolated, explicitly conditional Ti-Co CatMAP candidate input.

All numerical energies are copied from the archived Ti-Co CatMAP table and
NEB state-alignment audit. No DFT values are recalculated or edited.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


MODEL_SOURCE = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\alkene_epoxidation.mkm")
ENERGY_SOURCE = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\energies_dft_neb.txt")
NEB_ALIGNMENT = Path(r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\data\source_audit\tico_neb_state_alignment.csv")

NETWORK = """rxn_expressions = [
    '* + H2O2_g <-> H2O2*',
    'H2O2* + * -> Ha* + OOH*',
    '* + C6H12_g <-> C6H12*',
    'Ha* + OOH* + C6H12* <-> NEB1IS* + 2*',
    'NEB1IS* <-> NEB1TS* -> NEB1FS*',
    'NEB1FS* <-> NEB2IS* + C6H12O_g',
    'NEB2IS* <-> NEB2TS* -> NEB2FS*',
    'NEB2FS* <-> H2O*',
    'H2O* <-> H2O_g + *',
]

surface_names = ['tico']
temperature = 333.15
pressure_mode = 'concentration'
species_definitions = {
    's': {'site_names': ['111'], 'total': 1},
    'Ha_s': {'composition': {'H': 1}, 'n_sites': 1},
    'NEB1IS_s': {'composition': {'C': 6, 'H': 14, 'O': 2}, 'n_sites': 1},
    'NEB1TS_s': {'composition': {'C': 6, 'H': 14, 'O': 2}, 'n_sites': 1},
    'NEB1FS_s': {'composition': {'C': 6, 'H': 14, 'O': 2}, 'n_sites': 1},
    'NEB2IS_s': {'composition': {'H': 2, 'O': 1}, 'n_sites': 1},
    'NEB2TS_s': {'composition': {'H': 2, 'O': 1}, 'n_sites': 1},
    'NEB2FS_s': {'composition': {'H': 2, 'O': 1}, 'n_sites': 1},
    'C6H12_g': {'concentration': 1.0},
    'H2O2_g': {'concentration': 0.185},
    'H2O_g': {'concentration': 0.815},
    'C6H12O_g': {'concentration': 1e-20},
}

scaler = 'ThermodynamicScaler'
descriptor_names = ['temperature', 'pressure']
descriptor_ranges = None
resolution = 1
input_file = 'candidate_energies.txt'
data_file = 'candidate_tico_endpoint_model.pkl'
gas_thermo_mode = 'frozen_gas'
adsorbate_thermo_mode = 'frozen_adsorbate'
avoid_scaling = True
output_variables = ['turnover_frequency', 'forward_rate', 'reverse_rate', 'coverage']
decimal_precision = 100
tolerance = 1e-50
max_rootfinding_iterations = 300
max_bisections = 5
use_numbers_solver = False
"""


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    with ENERGY_SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
        fields = rows[0].keys()
    selected = [r for r in rows if r["surface_name"] in {"None", "tico"}]
    gas_names = {r["species_name"] for r in selected if r["surface_name"] == "None"}
    tico_names = {r["species_name"] for r in selected if r["surface_name"] == "tico"}
    required_gases = {"H2O2", "C6H12", "C6H12O", "H2O"}
    required_ads = {"H2O2", "Ha", "OOH", "C6H12", "C6H12O", "O", "OH", "H2O"}
    missing = (required_gases - gas_names) | (required_ads - tico_names)
    if missing:
        raise ValueError(f"required source rows missing: {sorted(missing)}")

    with NEB_ALIGNMENT.open(encoding="utf-8-sig", newline="") as stream:
        aligned = list(csv.DictReader(stream))
    neb_values: dict[tuple[str, str], str] = {}
    for row in aligned:
        state = row["state"]
        if row["path"] in {"NEB1", "NEB2"} and state in {"IS", "TS", "FS"}:
            neb_values[(row["path"], state)] = row["G_NEB_formation_eV"]
    if len(neb_values) != 6:
        raise ValueError(f"expected six NEB formation energies; found {len(neb_values)}")

    new_rows = [r.copy() for r in selected]
    for path, prefix, composition in (
        ("NEB1", "NEB1", "C6H14O2"),
        ("NEB2", "NEB2", "H2O"),
    ):
        for state in ("IS", "TS", "FS"):
            internal = f"{prefix}{state}"
            new_rows.append({
                "surface_name": "tico",
                "site_name": "111",
                "species_name": internal,
                "formation_energy": neb_values[(path, state)],
                "frequencies": "[]",
                "reference": f"FULL_{path}_COADSORPTION_{state}_GIBBS_FROM_ALIGNMENT_AUDIT",
            })

    energies = out / "candidate_energies.txt"
    with energies.open("w", encoding="utf-8", newline="") as stream:
        stream.write("\t".join(fields) + "\n")
        writer = csv.DictWriter(stream, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writerows(new_rows)
    model = out / "candidate_tico.mkm"
    model.write_text(
        "# EXPLORATORY / PROVISIONAL: full NEB states plus CatMAP endpoint barriers for links.\n"
        "# One effective site per complete coadsorbed state is an explicit coarse-graining assumption.\n"
        "# Not a validated geometric site count or a fully DFT-parameterized kinetic model.\n"
        + NETWORK,
        encoding="utf-8",
    )

    manifest = out / "candidate_input_manifest.sha256"
    manifest.write_text(
        "\n".join([
            f"{file_hash(MODEL_SOURCE)}  {MODEL_SOURCE}",
            f"{file_hash(ENERGY_SOURCE)}  {ENERGY_SOURCE}",
            f"{file_hash(NEB_ALIGNMENT)}  {NEB_ALIGNMENT}",
            f"{file_hash(model)}  {model}",
            f"{file_hash(energies)}  {energies}",
            "",
        ]), encoding="utf-8",
    )
    print(f"model={model}")
    print(f"energies={energies}")
    print(f"manifest={manifest}")
    print(f"source_rows={len(selected)}; appended_full_NEB_states={len(neb_values)}")


if __name__ == "__main__":
    main()
