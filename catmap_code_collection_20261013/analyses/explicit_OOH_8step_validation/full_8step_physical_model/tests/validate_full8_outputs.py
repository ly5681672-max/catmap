from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def close(actual: float, expected: float, tol: float = 1e-9) -> None:
    assert math.isclose(actual, expected, rel_tol=0, abs_tol=tol), (actual, expected)


def main() -> None:
    parser = json.loads((DATA / "PARSER_VALIDATION.json").read_text(encoding="utf-8"))
    assert parser["candidate_network_elementary_reactions"] == 12
    assert parser["old_candidate_elementary_reactions"] == 9
    assert parser["old_candidate_inventory_nullity"] == 1
    assert parser["original_inventory_nullity"] == 2
    assert parser["candidate_inventory_nullity"] == 2
    assert parser["candidate_atom_balance"] == "PASS"
    assert parser["candidate_formal_site_balance"] == "PASS"

    candidate_balances = rows("FULL_8STEP_BALANCE_CHECKS.csv")
    assert len(candidate_balances) == 12
    assert all(r["atom_balance"] == r["formal_site_balance"] == "PASS" for r in candidate_balances)
    assert len(rows("ORIGINAL_8STEP_STOICHIOMETRY.csv")) == 9
    assert len(rows("OLD_NINE_REACTION_STOICHIOMETRY.csv")) == 9
    assert len(rows("OLD_NINE_REACTION_INVENTORY_BASIS.csv")) == 1
    assert len(rows("ORIGINAL_8STEP_INVENTORY_BASIS.csv")) == 2
    assert len(rows("FULL_8STEP_INVENTORY_BASIS.csv")) == 2

    rxns = rows("REACTION_ENERGY_AND_BARRIERS.csv")
    assert len(rxns) == 12
    close(sum(float(r["deltaG_forward_eV"]) for r in rxns), -2.088914)
    neb1, neb2 = rxns[4], rxns[8]
    assert neb1["barrier_basis"] == neb2["barrier_basis"] == "EXPLICIT_DFT_NEB_TS"
    close(float(neb1["forward_barrier_eV"]), 0.556344)
    close(float(neb1["reverse_barrier_eV"]), 0.796030)
    close(float(neb2["forward_barrier_eV"]), 0.225994)
    close(float(neb2["reverse_barrier_eV"]), 0.706574)

    conn = {r["connection"]: float(r["deltaG_eV"]) for r in rows("NEB_CONNECTION_ENERGY_MAP.csv")}
    close(conn["NEB1-FS -> NEB2-IS + C6H12O_g"], 0.811237)
    close(conn["full NEB2-FS -> H2O*"], -0.601109)
    mkm = (ROOT / "inputs" / "FULL_8STEP_MODEL.mkm").read_text(encoding="utf-8")
    assert "OOH*" in mkm and "NEB1TS*" in mkm and "NEB2TS*" in mkm
    assert "CONDITIONAL" in mkm and "use_numbers_solver = False" in mkm
    db = rows("../runs/tico_provisional_zero_inventory/thermodynamic_detailed_balance.csv")
    assert len(db) == 12 and all(r["status"] == "PASS" for r in db)
    print("PASS: CatMAP parser artifacts, old-vs-full inventory dimensions, atom/formal-site balances, detailed balance, explicit NEB barriers, and net reaction energy closure")


if __name__ == "__main__":
    main()
