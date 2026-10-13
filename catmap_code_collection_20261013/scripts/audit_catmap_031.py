from __future__ import annotations

import ast
import csv
import json
from pathlib import Path
from typing import Any

from catmap.model import ReactionModel
from catmap.parsers.parser_base import ParserBase


def load_literal_assignment(path: Path, name: str) -> Any:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == name for target in node.targets
        ):
            return ast.literal_eval(node.value)
    raise KeyError(f"No literal assignment {name!r} in {path}")


def make_model(rxn_expressions: list[str], extra_definitions: dict[str, dict[str, Any]] | None = None) -> ReactionModel:
    model = ReactionModel()
    model.species_definitions = {
        "Ha_s": {"composition": {"H": 1}},
        "Hb_s": {"composition": {"H": 1}},
        "Hb-O_s": {"composition": {"H": 1, "O": 1}},
    }
    if extra_definitions:
        model.species_definitions.update(extra_definitions)
    model.parse_elementary_rxns(rxn_expressions)
    ParserBase(model)._baseparse()
    model.prefactor_list = None
    return model


def audit_original_model(model_path: Path) -> dict[str, Any]:
    expressions = load_literal_assignment(model_path, "rxn_expressions")
    model = make_model(expressions)
    model.verify()
    parsed = []
    for index, rxn in enumerate(model.elementary_rxns, 1):
        parsed.append({
            "reaction_index": index,
            "initial_state": rxn[0],
            "transition_state": rxn[1] if len(rxn) == 3 else None,
            "final_state": rxn[-1],
        })
    return {
        "catmap_version": __import__("catmap").__version__,
        "model_path": str(model_path),
        "reaction_count": len(model.elementary_rxns),
        "adsorbate_names": list(model.adsorbate_names),
        "transition_state_names": list(model.transition_state_names),
        "gas_names": list(model.gas_names),
        "parsed_reactions": parsed,
        "OOH-C6H12_is_stable_adsorbate": "OOH-C6H12_s" in model.adsorbate_names,
        "OOH-C6H12_is_transition_state": "OOH-C6H12_s" in model.transition_state_names,
        "Hb-O_is_stable_adsorbate": "Hb-O_s" in model.adsorbate_names,
        "Hb-O_is_transition_state": "Hb-O_s" in model.transition_state_names,
        "default_adsorbate_n_sites": model.species_definitions["OOH_s"]["n_sites"],
        "n_sites_override_supported": "n_sites" in model.species_definitions["OOH_s"],
        "full_element_and_site_balance_verification": "PASS",
    }


def barrier_values(model: ReactionModel, e_is: float, e_ts: float, e_fs: float) -> dict[str, float]:
    forward = model.get_rxn_energy([["IS_s"], ["TS-s"], ["FS_s"]],
                                   {"IS_s": e_is, "TS-s": e_ts, "FS_s": e_fs})
    reverse = model.get_rxn_energy([["FS_s"], ["TS-s"], ["IS_s"]],
                                   {"IS_s": e_is, "TS-s": e_ts, "FS_s": e_fs})
    return {
        "deltaG_forward": forward[0],
        "raw_TS_minus_IS": e_ts - e_is,
        "catmap_Ea_forward": forward[1],
        "raw_TS_minus_FS": e_ts - e_fs,
        "catmap_Ea_reverse": reverse[1],
    }


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    audit = audit_original_model(args.model)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
