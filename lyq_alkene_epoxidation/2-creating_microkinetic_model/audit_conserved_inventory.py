"""Check exact CatMAP network invariants and 333.15 K root inventories."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from catmap import ReactionModel
from mpmath import mp
from sympy import Matrix

import branch_audit as audit


INVARIANTS = {
    "I1_Ha_minus_OOH_OH_O": {"Ha_s": 1, "OOH_s": -1, "OH_s": -1, "O_s": -1},
    "I2_Hb_minus_O": {"Hb_s": 1, "O_s": -1},
}
ZERO_TOLERANCE = mp.mpf("1e-10")


def verify_network(config: Path) -> dict:
    model = ReactionModel(setup_file=str(config))
    names = list(model.adsorbate_names)
    stoich = audit._surface_stoichiometric_rows(model.elementary_rxns, names)
    matrix = Matrix(stoich)
    if matrix.rank() != len(names) - len(INVARIANTS):
        raise ValueError("Unexpected CatMAP surface stoichiometric rank")
    for label, weights in INVARIANTS.items():
        vector = Matrix([weights.get(name, 0) for name in names])
        if matrix.T * vector != Matrix.zeros(len(model.elementary_rxns), 1):
            raise ValueError(f"{label} is not conserved by parsed CatMAP reactions")
    return {
        "adsorbate_order": names,
        "reaction_count": len(model.elementary_rxns),
        "stoichiometric_rank": matrix.rank(),
        "independent_conserved_inventories": len(INVARIANTS),
        "invariants": INVARIANTS,
        "parsed_elementary_rxns": model.elementary_rxns,
        "source_config": str(config),
    }


def main() -> None:
    source = audit.AUDIT_ROOT / "combined_validated_roots.csv"
    roots = pd.read_csv(source, dtype=str, keep_default_na=False)
    first = roots.iloc[0]
    surface = first["source_surface"]
    config = Path(first["replay_root"]) / surface / f"baseline_{surface}.mkm"
    network = verify_network(config)
    (audit.AUDIT_ROOT / "stoichiometric_invariants.json").write_text(
        json.dumps(network, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    names = network["adsorbate_order"]
    mp.dps = 100
    output = []
    for row in roots.to_dict("records"):
        if row["coverage_order"].split(";")[:-1] != names:
            raise ValueError(f"Coverage order mismatch: {row['root_id']}")
        values = audit._mp_vector(row["full_coverage_vector_high_precision"])
        if len(values) != len(names) + 1:
            raise ValueError(f"Coverage length mismatch: {row['root_id']}")
        by_name = dict(zip(names, values))
        inventories = {
            label: sum(mp.mpf(coefficient) * by_name[name] for name, coefficient in weights.items())
            for label, weights in INVARIANTS.items()
        }
        output.append({
            "root_id": row["root_id"],
            "surface": row["source_surface"],
            "source_tag": row["source_tag"],
            "seed_id": row["source_seed_id"],
            **{label: mp.nstr(value, 90) for label, value in inventories.items()},
            "zero_inventory_compatible": all(abs(value) <= ZERO_TOLERANCE for value in inventories.values()),
            "invariant_tolerance": mp.nstr(ZERO_TOLERANCE, 20),
        })
    result = pd.DataFrame(output)
    result.to_csv(audit.AUDIT_ROOT / "conserved_inventory.csv", index=False)
    print(f"roots={len(result)} zero_inventory_compatible={result.zero_inventory_compatible.sum()}")
    print(result.groupby("surface")["zero_inventory_compatible"].agg(["count", "sum"]).to_string())


if __name__ == "__main__":
    main()
