from __future__ import annotations

import json
import os
from pathlib import Path

import mpmath as mp


OUT = Path(r"C:\Users\Administrator\Documents\codex探索\tizr_precision_root_cause_20261021")
ROOT = Path(r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\multisite_conditional_replay_20261021\sites\tizr\constrained_tizr_result.json")


def main() -> None:
    data = json.loads(ROOT.read_text(encoding="utf-8"))
    old = Path.cwd()
    os.chdir(OUT / "inputs")
    try:
        from catmap import ReactionModel
        model = ReactionModel(setup_file="tizr_precision_replay.mkm")
        model.parse()
        mp.mp.dps = 385
        model._math, model._mpfloat, model._matrix, model._Axb_solver = mp, mp.mpf, mp.matrix, mp.lu_solve
        model._kB, model._h = mp.mpf("8.617332478e-5"), mp.mpf("4.135667516e-15")
        model._math.infnorm = lambda x: mp.norm(x, "inf")
        model._math.lsqnorm = lambda x: mp.norm(x, 2) ** 2
        solver = model.solver
        solver.compile()
        names = list(solver.adsorbate_names)
        coverage = [mp.mpf(data["coverage"][name]) for name in names]
        vacancy = mp.mpf(data["vacant"])
        params = list(model.scaler.get_rxn_parameters([mp.mpf("333.15")]))
        constants = solver.get_rate_constants(params, coverage)
        net = list(solver.get_rxn_rates(coverage, constants))
        tof = list(solver.get_turnover_frequency(params, rates=net, verify_coverages=False))
        product_i = solver.gas_names.index("C6H12O_g")
        n_rxn = len(solver.elementary_rxns)
        forward_exprs = [solver.rate_equation_term(rxn[0], f"kf[{j}]")
                         for j, rxn in enumerate(solver.elementary_rxns)]
        reverse_exprs = [solver.rate_equation_term(rxn[-1], f"kr[{j}]")
                         for j, rxn in enumerate(solver.elementary_rxns)]
        env = {"kf": constants[:n_rxn], "kr": constants[n_rxn:],
               "p": list(solver.gas_pressures), "theta": coverage,
               "s": [vacancy], "mp": mp}
        forward = [eval(expr, {"__builtins__": {}}, env) for expr in forward_exprs]
        reverse = [eval(expr, {"__builtins__": {}}, env) for expr in reverse_exprs]
        product_reactions = [(i + 1, rate) for i, (reaction, rate) in enumerate(zip(solver.elementary_rxns, net))
                             if "C6H12O_g" in reaction[0] or "C6H12O_g" in reaction[-1]]
        out = {
            "coverage_snapshot_re_evaluation_dps": mp.mp.dps,
            "stored_tof": data["C6H12O_net_tof"],
            "recomputed_catmap_tof_from_stored_coverage": mp.nstr(tof[product_i], 80),
            "stored_product_net_rate": data["reaction_net_rates"][4],
            "recomputed_product_reaction_rates": [(i, mp.nstr(v, 80)) for i, v in product_reactions],
            "recomputed_reaction_net_rates": [mp.nstr(v, 80) for v in net],
            "directional_rates_mpmath": [
                {"reaction_index": j + 1, "forward": mp.nstr(forward[j], 80),
                 "reverse": mp.nstr(reverse[j], 80), "net": mp.nstr(net[j], 80),
                 "forward_reverse_identity_error": mp.nstr(abs(forward[j] - reverse[j] - net[j]), 40),
                 "cancellation_factor": mp.nstr((abs(forward[j]) + abs(reverse[j])) / abs(net[j]), 40)
                    if net[j] != 0 else "NET_ZERO"}
                for j in range(n_rxn)],
            "stored_vs_recomputed_tof_relative_error": mp.nstr(
                abs(tof[product_i] - mp.mpf(data["C6H12O_net_tof"])) / abs(mp.mpf(data["C6H12O_net_tof"])), 50),
            "stored_vs_recomputed_tof_delta_log10": mp.nstr(
                abs(mp.log10(tof[product_i]) - mp.log10(mp.mpf(data["C6H12O_net_tof"]))), 50),
            "vacancy_plus_adsorbates": mp.nstr(vacancy + sum(coverage), 50),
            "surface_rhs": [mp.nstr(v, 50) for v in solver.steady_state_function(coverage)],
            "catmap_reaction_order": [str(x) for x in solver.elementary_rxns],
        }
        (OUT / "data/archived_coverage_re_evaluation.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(json.dumps({k: out[k] for k in ("stored_tof", "recomputed_catmap_tof_from_stored_coverage",
                                                "stored_vs_recomputed_tof_relative_error",
                                                "stored_vs_recomputed_tof_delta_log10",
                                                "recomputed_product_reaction_rates", "surface_rhs",
                                                "directional_rates_mpmath")}, indent=2))
    finally:
        os.chdir(old)


if __name__ == "__main__":
    main()
