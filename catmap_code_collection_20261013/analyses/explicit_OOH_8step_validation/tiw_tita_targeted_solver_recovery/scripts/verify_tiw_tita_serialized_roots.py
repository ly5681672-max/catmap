from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import csv
import mpmath as mp

PROJECT = Path(r"H:\catmap\catmap\catmap_code_collection_20261013")
BASE = Path(os.environ.get("TIW_TITA_OUTPUT_DIR", str(PROJECT / "analyses" / "explicit_OOH_8step_validation" / "tiw_tita_targeted_solver_recovery")))
RECOVERY = BASE / "combined_targeted_recovery.json"
(BASE / "data").mkdir(parents=True, exist_ok=True)
SCRIPT = Path(__file__).resolve().with_name("tiw_tita_damped_recovery.py")
SPEC = importlib.util.spec_from_file_location("targeted", SCRIPT)
targeted = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(targeted)

all_data = json.loads(RECOVERY.read_text(encoding="utf-8"))
verified = {}
for site, dps in (("tiw", 300), ("tita", 180)):
    candidates = all_data.get(site, {}).get("catmap_newton_results", [])
    record = next((r for r in reversed(candidates)
                   if r.get("dps") == dps and r.get("accepted_native_equation")), None)
    model, solver, params, oldcwd = targeted.setup(site, dps, tag="serialized_verify_20261022")
    try:
        ctx = targeted.make_context(solver, params)
        names, idx, free_idx, rows, nu, state, residual, _, get_fluxes = ctx
        if record is not None:
            # Preferred path: replay the full-precision NewtonRoot coordinates.
            z = [mp.mpf(v) for v in record["root_logits"]]
            source = "full_precision_root_logits"
            saved_coverage = record["coverage"]
        else:
            # Recovery path: the final solver JSON may be absent after an
            # interrupted rerun. Reconstruct the exact-inventory chart from
            # the already archived high-significance coverage CSV. Vacancy is
            # recomputed from the chart identity rather than its rounded CSV
            # field, preserving the tiny Ti-W vacancy.
            table = BASE / "deliverables" / "data" / "TIW_TITA_SOLVER_RESULTS.csv"
            with table.open(encoding="utf-8-sig", newline="") as handle:
                row = next(r for r in csv.DictReader(handle) if r["surface"] == site)
            saved_coverage = {n: row[f"theta_{n}"] for n in names}
            free_coverage = [mp.mpf(row[f"theta_{n}"]) for n in targeted.FREE]
            vacancy = mp.mpf(1) - sum(mp.mpf(w) * c
                                      for w, c in zip(targeted.WEIGHTS, free_coverage))
            z = [mp.log(c / vacancy) for c in free_coverage]
            source = "reconstructed_from_archived_high_precision_coverage_csv"
        cv, vac = state(z)
        parts = residual(z, return_parts=True)
        cv, vac, ks, fwd, rev, net, f, prod, cons, gross, local, logratio = parts
        tof = net[4]
        site_sum = sum(cv) + (mp.mpf(1) - sum(cv))
        i1 = cv[idx["Ha_s"]] - cv[idx["OOH_s"]] - cv[idx["OH_s"]] - cv[idx["O_s"]]
        i2 = cv[idx["Hb_s"]] - cv[idx["O_s"]]
        verified[site] = {
            "dps": dps,
            "source_precision": dps,
            "replay_source": source,
            "tof_s-1": mp.nstr(tof, 90),
            "vacancy": mp.nstr(vac, 70),
            "OOH_s": mp.nstr(cv[idx["OOH_s"]], 70),
            "C6H12_s": mp.nstr(cv[idx["C6H12_s"]], 100),
            "reconstructed_vs_saved_C6H12_abs": mp.nstr(abs(cv[idx["C6H12_s"]] - mp.mpf(saved_coverage["C6H12_s"])), 50),
            "reconstructed_vs_saved_coverage_abs": {n: mp.nstr(abs(cv[i] - mp.mpf(saved_coverage[n])), 20)
                                                    for i, n in enumerate(names)},
            "rate_constants_first_pair": [mp.nstr(ks[0], 50), mp.nstr(ks[len(net)], 50)],
            "fwd_rxn1": mp.nstr(fwd[0], 80), "rev_rxn1": mp.nstr(rev[0], 80),
            "root_logit_C6H12": mp.nstr(z[1], 100),
            "max_local_residual": mp.nstr(max(abs(v) for v in local), 40),
            "I1": mp.nstr(i1, 50), "I2": mp.nstr(i2, 50),
            "site_sum": mp.nstr(site_sum, 50),
            "net_rates": [mp.nstr(v, 70) for v in net],
            "serialized_state_replayed": True,
            "replay_coordinates": source,
            "weighted_site_balance_reconstructed_vacancy": mp.nstr(
                mp.mpf(1) - sum(mp.mpf(w) * cv[idx[n]]
                                for n, w in zip(targeted.FREE, targeted.WEIGHTS)), 70),
        }
    finally:
        import os
        os.chdir(oldcwd)

(BASE / "data" / "serialized_root_replay.json").write_text(
    json.dumps(verified, indent=2, ensure_ascii=False), encoding="utf-8"
)
print(json.dumps(verified, indent=2, ensure_ascii=False))
