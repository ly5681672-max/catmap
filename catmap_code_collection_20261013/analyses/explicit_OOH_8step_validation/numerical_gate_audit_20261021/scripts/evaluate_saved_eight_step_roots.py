from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from collections import Counter
import mpmath as mp

RESULTS = Path(r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\multisite_conditional_replay_20261021")
OUT_DIR = Path(r"C:\Users\Administrator\Documents\codex探索\numerical_gate_audit_20261021")
SITES = ["tico", "ti", "tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "titi", "tiw", "tita", "ticr", "tini"]
INDEPENDENT = ["C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "OH_s", "OOH_s", "O_s"]


def result_path(site: str) -> Path:
    return RESULTS / "sites" / site / f"constrained_{site}_result.json"


def main() -> None:
    mp.mp.dps = 110
    rows = []
    for site in SITES:
        p = result_path(site)
        if not p.exists():
            rows.append({"surface": site, "status": "NO_ACCEPTED_ROOT_RESULT"})
            continue
        result = json.loads(p.read_text(encoding="utf-8"))
        work = p.parent
        old = Path.cwd()
        os.chdir(work)
        try:
            from catmap import ReactionModel
            setup = f"{site}_static_probe.mkm"
            if not (work / setup).exists():
                setup = f"{site}_probe.mkm"
            model = ReactionModel(setup_file=setup)
            model.parse()  # parse/build only; deliberately no model.run()
            solver, scaler = model.solver, model.scaler
            mp.mp.dps = 385
            # Mirror only CatMAP's documented numerical-context setup from
            # ReactionModel.run(); avoid invoking any descriptor/steady solve.
            model._math = mp
            model._mpfloat = mp.mpf
            model._matrix = mp.matrix
            model._Axb_solver = mp.lu_solve
            model._kB = mp.mpf("8.617332478e-5")
            model._h = mp.mpf("4.135667516e-15")
            model._math.infnorm = lambda x: mp.norm(x, "inf")
            model._math.lsqnorm = lambda x: mp.norm(x, 2) ** 2
            solver.compile()
            params = list(scaler.get_rxn_parameters([mp.mpf("333.15")]))
            names = list(solver.adsorbate_names)
            idx = {name: i for i, name in enumerate(names)}
            cv = [mp.mpf(result["coverage"][name]) for name in names]
            solver._coverage = list(cv)
            ks = solver.get_rate_constants(params, list(cv))
            directions = list(solver.get_directional_rates(params))
            n_rxn = len(solver.elementary_rxns)
            net = list(solver.steady_state_function(list(cv)))
            stoich = [[0] * n_rxn for _ in names]
            for j, rxn in enumerate(solver.elementary_rxns):
                ci, cf = Counter(rxn[0]), Counter(rxn[-1])
                for sp in names:
                    stoich[idx[sp]][j] = cf.get(sp, 0) - ci.get(sp, 0)
            species = {}
            for name in INDEPENDENT:
                i = idx[name]
                gross = sum(abs(stoich[i][j]) * (abs(directions[j]) + abs(directions[n_rxn + j])) for j in range(n_rxn))
                residual = abs(net[i])
                species[name] = {
                    "abs_residual": mp.nstr(residual, 35),
                    "gross_directional_scale": mp.nstr(gross, 35),
                    "normalized_residual": mp.nstr(residual / gross, 35) if gross else "UNDEFINED_ZERO_GROSS_FLOW",
                    "coverage": result["coverage"][name],
                }
            vacancy = mp.mpf(result["vacant"])
            inv1 = cv[idx["Ha_s"]] - cv[idx["OOH_s"]] - cv[idx["OH_s"]] - cv[idx["O_s"]]
            inv2 = cv[idx["Hb_s"]] - cv[idx["O_s"]]
            tof = list(solver.get_turnover_frequency(params, rates=solver.get_rxn_rates(list(cv), ks), verify_coverages=False))
            product_idx = solver.gas_names.index("C6H12O_g")
            # CatMAP 0.3.1 get_directional_rates calls get_elem_ec(), which
            # casts each directional rate to Python float. Do not use it to
            # subtract near-equilibrium fluxes at ultra-low TOF. Instead,
            # reconstruct turnover from CatMAP's high-precision net rates,
            # matching get_turnover_frequency's signed reaction membership.
            net_rates = list(solver.get_rxn_rates(list(cv), ks))
            product_rxn_indices = [j for j, rxn in enumerate(solver.elementary_rxns)
                                   if "C6H12O_g" in rxn[0] or "C6H12O_g" in rxn[-1]]
            product_net_from_catmap_rates = sum(
                (-net_rates[j] if "C6H12O_g" in solver.elementary_rxns[j][0] else net_rates[j])
                for j in product_rxn_indices
            )
            product_tof = tof[product_idx]
            product_net_error = abs(product_tof - product_net_from_catmap_rates)
            product_net_relative_error = product_net_error / max(abs(product_tof), abs(product_net_from_catmap_rates)) if product_tof or product_net_from_catmap_rates else mp.mpf(0)
            rows.append({
                "surface": site,
                "status": "EVALUATED_FROM_SAVED_ROOT_NO_RESOLVE",
                "n_rxn": n_rxn,
                "site_sum_error": mp.nstr(abs(sum(cv) + vacancy - 1), 20),
                "I1": mp.nstr(abs(inv1), 20),
                "I2": mp.nstr(abs(inv2), 20),
                "max_abs_residual_independent": mp.nstr(max(abs(net[idx[n]]) for n in INDEPENDENT), 35),
                "max_species_normalized_residual": mp.nstr(max(mp.mpf(v["normalized_residual"]) for v in species.values() if v["normalized_residual"] != "UNDEFINED_ZERO_GROSS_FLOW"), 35),
                "TOF_C6H12O": mp.nstr(product_tof, 40),
                "product_net_from_catmap_net_rates": mp.nstr(product_net_from_catmap_rates, 40),
                "product_TOF_vs_net_rate_abs_error": mp.nstr(product_net_error, 20),
                "product_TOF_vs_net_rate_relative_error": mp.nstr(product_net_relative_error, 20),
                "net_reaction_rates": json.dumps([mp.nstr(x, 35) for x in net_rates]),
                "species": json.dumps(species),
            })
        finally:
            os.chdir(old)
    out = OUT_DIR / "saved_root_local_flux_audit.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8-sig", newline="") as f:
        fieldnames = list(dict.fromkeys(k for r in rows for k in r))
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(out)
    for row in rows:
        print(row["surface"], row.get("status"), row.get("max_abs_residual_independent"), row.get("max_species_normalized_residual"), row.get("TOF_C6H12O"))


if __name__ == "__main__":
    main()
