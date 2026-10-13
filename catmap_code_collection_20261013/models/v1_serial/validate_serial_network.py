"""Parse and audit the six-step serial network with local CatMAP 0.3.1."""

from __future__ import annotations

import json
import runpy
from fractions import Fraction
from pathlib import Path
from typing import Any

from catmap.model import ReactionModel
from catmap.parsers.parser_base import ParserBase


HERE = Path(__file__).resolve().parent
MODEL_FILE = HERE / "ti_re_serial_neb.mkm"
DATA_DIR = HERE.parents[1] / "data"
UNACTIVATED_BARRIER_STATUS = "SOFTWARE_ASSUMPTION_NOT_DFT"

STOICH_SPECIES = (
    "Rpair_s", "EpoxideHHO_s", "HHO_s", "OHHa_s", "H2O_s",
    "H2O2_g", "C6H12_g", "C6H12O_g", "H2O_g", "*",
)

COMPOSITIONS = {
    "Rpair_s": {"C": 6, "H": 14, "O": 2},
    "TS1_s": {"C": 6, "H": 14, "O": 2},
    "EpoxideHHO_s": {"C": 6, "H": 14, "O": 2},
    "HHO_s": {"H": 2, "O": 1},
    "TS2_s": {"H": 2, "O": 1},
    "OHHa_s": {"H": 2, "O": 1},
    "H2O_s": {"H": 2, "O": 1},
    "H2O2_g": {"H": 2, "O": 2},
    "C6H12_g": {"C": 6, "H": 12},
    "C6H12O_g": {"C": 6, "H": 12, "O": 1},
    "H2O_g": {"H": 2, "O": 1},
    "*": {},
}

REACTION_STOICHIOMETRY = (
    {"Rpair_s": 1, "H2O2_g": -1, "C6H12_g": -1, "*": -1},
    {"Rpair_s": -1, "EpoxideHHO_s": 1},
    {"EpoxideHHO_s": -1, "HHO_s": 1, "C6H12O_g": 1},
    {"HHO_s": -1, "OHHa_s": 1},
    {"OHHa_s": -1, "H2O_s": 1},
    {"H2O_s": -1, "H2O_g": 1, "*": 1},
)


def build_model() -> ReactionModel:
    config = runpy.run_path(str(MODEL_FILE))
    model = ReactionModel()
    model.surface_names = config["surface_names"]
    model.species_definitions = {key: dict(value) for key, value in config["species_definitions"].items()}
    model.parse_elementary_rxns(config["rxn_expressions"])
    ParserBase(model)._baseparse()
    model.prefactor_list = None
    return model


def _rank(matrix: list[list[int]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next((r for r in range(pivot_row, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        divisor = work[pivot_row][column]
        work[pivot_row] = [value / divisor for value in work[pivot_row]]
        for r in range(len(work)):
            if r != pivot_row and work[r][column]:
                scale = work[r][column]
                work[r] = [a - scale * b for a, b in zip(work[r], work[pivot_row])]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def _matrix() -> list[list[int]]:
    return [[reaction.get(species, 0) for reaction in REACTION_STOICHIOMETRY]
            for species in STOICH_SPECIES]


def _conserved_totals() -> list[str]:
    totals = {
        element: [COMPOSITIONS[s].get(element, 0) for s in STOICH_SPECIES]
        for element in ("C", "H", "O")
    }
    totals["site"] = [1 if species.endswith("_s") or species == "*" else 0 for species in STOICH_SPECIES]
    matrix = _matrix()
    conserved = []
    for name, vector in totals.items():
        if all(sum(vector[row] * matrix[row][column] for row in range(len(STOICH_SPECIES))
                   ) == 0 for column in range(6)):
            conserved.append(name)
    return conserved


def _conserved_total_rank() -> int:
    vectors = [
        [COMPOSITIONS[s].get(element, 0) for s in STOICH_SPECIES]
        for element in ("C", "H", "O")
    ]
    vectors.append([1 if species.endswith("_s") or species == "*" else 0
                    for species in STOICH_SPECIES])
    return _rank([[vectors[column][row] for column in range(len(vectors))]
                  for row in range(len(STOICH_SPECIES))])


def stoichiometry_report() -> dict[str, Any]:
    matrix = _matrix()
    rank = _rank(matrix)
    overall = {species: sum(reaction.get(species, 0) for reaction in REACTION_STOICHIOMETRY)
               for species in STOICH_SPECIES}
    overall = {species: coefficient for species, coefficient in overall.items() if coefficient}
    expected = {"C6H12_g": -1, "H2O2_g": -1, "C6H12O_g": 1, "H2O_g": 1}
    return {
        "species_order": list(STOICH_SPECIES),
        "reaction_order": [f"step_{n}" for n in range(1, 7)],
        "stoichiometric_matrix": matrix,
        "matrix_rank": rank,
        "right_nullity": len(REACTION_STOICHIOMETRY) - rank,
        "left_nullity": len(STOICH_SPECIES) - rank,
        "independent_conserved_totals": _conserved_totals(),
        "conserved_total_rank": _conserved_total_rank(),
        "overall_net": overall,
        "overall_matches_target": overall == expected,
        "ts_species_excluded_from_stable_state_matrix": ["TS1_s", "TS2_s"],
    }


def reaction_balance_rows() -> list[dict[str, Any]]:
    reactions = [
        ("step_1_feed_capture", ["*", "H2O2_g", "C6H12_g"], ["Rpair_s"], []),
        ("step_2_NEB1", ["Rpair_s"], ["EpoxideHHO_s"], ["TS1_s"]),
        ("step_3_epoxide_desorption", ["EpoxideHHO_s"], ["HHO_s", "C6H12O_g"], []),
        ("step_4_NEB2", ["HHO_s"], ["OHHa_s"], ["TS2_s"]),
        ("step_5_water_formation", ["OHHa_s"], ["H2O_s"], []),
        ("step_6_water_desorption", ["H2O_s"], ["H2O_g", "*"], []),
    ]
    site_counts = {species: (1 if species.endswith("_s") or species == "*" else 0)
                   for species in COMPOSITIONS}
    result = []
    for name, initial, final, transition in reactions:
        def atoms(side: list[str]) -> dict[str, int]:
            total: dict[str, int] = {}
            for species in side:
                for element, count in COMPOSITIONS[species].items():
                    total[element] = total.get(element, 0) + count
            return total
        left_atoms, right_atoms = atoms(initial), atoms(final)
        left_sites = sum(site_counts[species] for species in initial)
        right_sites = sum(site_counts[species] for species in final)
        ts_atoms = atoms(transition) if transition else {}
        result.append({
            "reaction": name,
            "initial_state": "+".join(initial),
            "transition_state": "+".join(transition),
            "final_state": "+".join(final),
            "initial_atoms": json.dumps(left_atoms, sort_keys=True),
            "final_atoms": json.dumps(right_atoms, sort_keys=True),
            "atom_balance": "PASS" if left_atoms == right_atoms else "FAIL",
            "transition_state_atoms": json.dumps(ts_atoms, sort_keys=True),
            "initial_site_count": left_sites,
            "final_site_count": right_sites,
            "site_balance": "PASS" if left_sites == right_sites else "FAIL",
            "n_sites_evidence": "ASSUMED_SINGLE_EFFECTIVE_CENTER",
        })
    return result


def model_metadata() -> dict[str, Any]:
    config = runpy.run_path(str(MODEL_FILE))
    return {
        "catmap_version": __import__("catmap").__version__,
        "temperature_K": config["temperature"],
        "reservoir_pressures": {name: config["species_definitions"][name]["pressure"]
                                for name in ("C6H12_g", "H2O2_g", "H2O_g", "C6H12O_g")},
        "site_occupancy_evidence": "ASSUMED_SINGLE_EFFECTIVE_CENTER",
        "kinetics_ready": False,
        "tof_calculation_authorized": config["tof_calculation_authorized"],
        "unactivated_step_barrier_status": UNACTIVATED_BARRIER_STATUS,
    }


def audit() -> dict[str, Any]:
    model = build_model()
    model.verify()
    return {
        "metadata": model_metadata(),
        "reaction_count": len(model.elementary_rxns),
        "stable_adsorbates": list(model.adsorbate_names),
        "transition_states": list(model.transition_state_names),
        "stoichiometry": stoichiometry_report(),
        "catmap_balance_check": "PASS",
    }


def write_balance_csv(path: Path) -> None:
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = reaction_balance_rows()
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    write_balance_csv(DATA_DIR / "SERIAL_REACTION_BALANCE.csv")
    print(json.dumps(audit(), ensure_ascii=False, indent=2))
