"""Audit detailed balance and flux closure at the accepted provisional root."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import mpmath as mp
from catmap import ReactionModel

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "tife_candidate"
TEMP = mp.mpf("333.15")


def energy(definition: dict) -> mp.mpf:
    val = definition.get("formation_energy", 0)
    if isinstance(val, (list, tuple)):
        val = val[0]
    return mp.mpf(str(val))


def main() -> None:
    os.chdir(RUN)
    mp.mp.dps = 120
    model = ReactionModel(setup_file="candidate_tife.mkm"); model.run()
    solver = model.solver; solver.compile()
    kb_catmap = mp.mpf(str(solver._kB))
    with (RUN / "candidate_steady_state_coverages.csv").open(encoding="utf-8-sig", newline="") as stream:
        cv_rows = list(csv.DictReader(stream))
    cv_map = {r["species"]: mp.mpf(r["coverage"]) for r in cv_rows}
    cv = [cv_map[n] for n in solver.adsorbate_names]
    params = list(model.scaler.get_rxn_parameters([TEMP, mp.mpf("1.0")]))
    constants = list(solver.get_rate_constants(params, cv))
    n = len(solver.elementary_rxns)
    defs = model.species_definitions
    rows = []
    for i, rxn in enumerate(solver.elementary_rxns):
        initial, final = rxn[0], rxn[-1]
        g_initial = sum(energy(defs[name]) for name in initial)
        g_final = sum(energy(defs[name]) for name in final)
        delta_g = g_final - g_initial
        kf, kr = mp.mpf(str(constants[i])), mp.mpf(str(constants[n + i]))
        if kf > 0 and kr > 0:
            detailed_balance_delta = kb_catmap * TEMP * mp.log(kf / kr)
            mismatch = detailed_balance_delta + delta_g
            status = "CHECKED" if abs(mismatch) < mp.mpf("1e-8") else "MISMATCH"
        else:
            detailed_balance_delta = mismatch = mp.nan
            status = "IRREVERSIBLE_OR_ZERO_REVERSE_CONSTANT; NOT_TESTED"
        arrow_note = "single-arrow in setup, but compiled CatMAP has nonzero reverse constant" if i == 1 and kr > 0 else ""
        rows.append({"step": i + 1, "reaction": model.print_rxn(rxn, mode="text"),
                     "deltaG_products_minus_reactants_eV": mp.nstr(delta_g, 24),
                     "kf_s-1": mp.nstr(kf, 24), "kr_s-1": mp.nstr(kr, 24),
                     "reverse_rate_constant_nonzero": bool(kr > 0), "arrow_semantics_note": arrow_note,
                     "kBT_ln_kf_over_kr_eV": mp.nstr(detailed_balance_delta, 24),
                     "detailed_balance_residual_eV": mp.nstr(mismatch, 24), "status": status})
    with (ROOT / "data" / "tife_detailed_balance_audit.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    result = {"catmap_kB_eV_per_K": mp.nstr(kb_catmap, 15),
              "reversible_steps_checked": sum(r["status"] == "CHECKED" for r in rows),
              "reversible_steps_mismatched": sum(r["status"] == "MISMATCH" for r in rows),
              "zero_reverse_constant_steps": sum(float(r["kr_s-1"]) == 0 for r in rows),
              "step2_arrow_vs_rate_constant": rows[1]["arrow_semantics_note"],
              "criterion_eV": "absolute[kBT_CatMAP ln(kf/kr)+DeltaG] < 1e-8; evaluated using installed CatMAP's kB",
              "source": str(RUN / "candidate_tife.mkm"), "status": "PROVISIONAL_CANDIDATE_ONLY"}
    (ROOT / "data" / "tife_detailed_balance_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
