"""Conditional Ti-Co steady-state solve using CatMAP 0.3.1 rate functions.

The reaction rates, rate constants, and site-balance model are CatMAP-native.
An outer log-coverage chart enforces the candidate network's exact Ha=OOH
inventory invariant. Results remain PROVISIONAL because five links use the
CatMAP endpoint-barrier approximation and the effective site count of full
coadsorbed states is not structurally verified.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path

import mpmath as mp


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--precision", type=int, default=160, help="working mpmath decimal precision")
    p.add_argument("--seed-coverages", type=Path, help="optional positive coverage seed from a compatible candidate root")
    args = p.parse_args()
    model_path = args.model.resolve()
    run_dir = model_path.parent
    previous_cwd = Path.cwd()
    os.chdir(run_dir)
    try:
        from catmap import ReactionModel
        from catmap.solvers.solver_base import NewtonRoot
        import catmap

        model = ReactionModel(setup_file=model_path.name)
        model.run()
        solver, scaler = model.solver, model.scaler
        solver.compile()
        mp.mp.dps = args.precision
        names = list(solver.adsorbate_names)
        required = ["Ha_s", "OOH_s", "C6H12O_g"]
        if any(name not in names + list(solver.gas_names) for name in required):
            raise RuntimeError(f"candidate species missing: {required}; ads={names}; gas={solver.gas_names}")
        ha_i, ooh_i = names.index("Ha_s"), names.index("OOH_s")
        independent_names = [name for name in names if name != "Ha_s"]
        independent_indices = [names.index(name) for name in independent_names]
        if "OOH_s" not in independent_names:
            raise RuntimeError("OOH must be an independent coordinate")
        descriptors = [mp.mpf("333.15"), mp.mpf("1.0")]
        rxn_params = list(scaler.get_rxn_parameters(descriptors))
        n_rxn = len(solver.elementary_rxns)
        if n_rxn != 9:
            raise RuntimeError(f"expected 9 reactions, found {n_rxn}")

        # q_i = exp(z_i); Ha is reconstructed equal to OOH, preserving the
        # zero-inventory class. All adsorbates have n_sites=1 in this
        # explicitly coarse-grained candidate; vacancy follows exact balance.
        def state_from_logits(z):
            weights = [mp.exp(v) for v in z]
            ooh_weight = weights[independent_names.index("OOH_s")]
            denominator = 1 + sum(weights) + ooh_weight
            cv = [mp.mpf("0") for _ in names]
            for i, w in zip(independent_indices, weights):
                cv[i] = w / denominator
            cv[ha_i] = ooh_weight / denominator
            vacancy = 1 / denominator
            return cv, vacancy

        def evaluate(z, include_rates=False):
            cv, vacancy = state_from_logits(z)
            ks = list(solver.get_rate_constants(rxn_params, cv))
            net = list(solver.steady_state_function(cv))
            if not include_rates:
                return cv, vacancy, net
            rates = list(solver.get_rxn_rates(cv, ks))
            tof = list(solver.get_turnover_frequency(rxn_params, rates=rates, verify_coverages=False))
            return cv, vacancy, net, rates, tof

        # Stoichiometric rows for stable adsorbates; Ha equation is redundant
        # with OOH in this exact invariant subspace, so solve the other 8 rows.
        stoich = [[0] * n_rxn for _ in names]
        for j, rxn in enumerate(solver.elementary_rxns):
            initial, final = Counter(rxn[0]), Counter(rxn[-1])
            for i, name in enumerate(names):
                stoich[i][j] = final.get(name, 0) - initial.get(name, 0)

        seeds = []
        # Generic positive seeds plus a numerical seed projected from the
        # archived explicit-OOH provisional root, if available.
        if args.seed_coverages:
            seed_rows = {}
            with args.seed_coverages.resolve().open(encoding="utf-8-sig", newline="") as f:
                for row in csv.DictReader(f):
                    seed_rows[row["species"]] = mp.mpf(row["coverage"])
            if set(independent_names).issubset(seed_rows) and "*" in seed_rows:
                vacancy0 = seed_rows["*"]
                if vacancy0 > 0 and all(seed_rows[name] > 0 for name in independent_names):
                    seeds.append(("base_candidate_root_numerical_seed_only",
                                  [mp.log(seed_rows[name] / vacancy0) for name in independent_names]))
        for value in ("-30", "-15", "-5", "0", "3"):
            seeds.append((f"uniform_{value}", [mp.mpf(value)] * len(independent_names)))
        archived = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\analysis_ooh\zero_inventory_tico_diagnostic\data\tico_zero_inventory_root.json")
        seed_source = "uniform multistart only"
        if archived.is_file():
            root_data = json.loads(archived.read_text(encoding="utf-8"))
            cov_map = root_data.get("coverages", root_data.get("coverage", {}))
            if isinstance(cov_map, dict) and all(n in cov_map for n in independent_names):
                vac0 = mp.mpf(str(root_data.get("vacancy", 0)))
                if vac0 > 0 and all(mp.mpf(str(cov_map[n])) > 0 for n in independent_names):
                    z = [mp.log(mp.mpf(str(cov_map[n])) / vac0) for n in independent_names]
                    seeds.insert(0, ("archived_root_numerical_seed_only", z))
                    seed_source = "old root used as numerical seed only; no old output reused"

        attempts = []
        accepted_roots = []
        for seed_label, start in seeds:
            cv0, _, _ = evaluate(start)
            solver._coverage = list(cv0)
            directional = list(solver.get_directional_rates(rxn_params))
            row_scales = []
            for i in independent_indices:
                gross = sum(abs(stoich[i][j]) * (abs(directional[j]) + abs(directional[n_rxn + j]))
                            for j in range(n_rxn))
                row_scales.append(max(gross, mp.mpf("1e-300")))
            global_scale = max(row_scales)
            row_scales = [max(s, global_scale * mp.mpf("1e-250")) for s in row_scales]

            def residual(z):
                _, _, net = evaluate(z)
                return tuple(net[i] / scale for i, scale in zip(independent_indices, row_scales))

            def jacobian(z):
                z = list(z)
                return mp.matrix([[mp.diff(lambda q: residual([q if k == col else z[k]
                                                               for k in range(len(z))])[row], z[col])
                                   for col in range(len(z))] for row in range(len(z))])

            root = list(start)
            scaled = mp.inf
            error = None
            iterations = 0
            try:
                iterator = NewtonRoot(residual, start, mp.matrix, mp.mpf, mp.lu_solve,
                                      J=jacobian, norm=lambda x: mp.norm(x, "inf"), verbose=0)
                for iterations, (root_matrix, scaled) in enumerate(iterator, 1):
                    root = list(root_matrix)
                    if scaled < mp.mpf("1e-45") or iterations >= 200:
                        break
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
            cv, vacancy, net, rates, tof = evaluate(root, True)
            max_abs = max(abs(net[i]) for i in independent_indices)
            max_rate = max(abs(v) for v in rates)
            relative = max_abs / max(max_rate, mp.mpf("1e-300"))
            inventory = cv[ha_i] - cv[ooh_i]
            site_error = abs(sum(cv) + vacancy - 1)
            accepted_here = (scaled < mp.mpf("1e-40") and relative < mp.mpf("1e-7")
                             and site_error < mp.mpf("1e-100") and abs(inventory) < mp.mpf("1e-100")
                             and vacancy >= 0 and all(v >= 0 for v in cv))
            attempt = {"seed": seed_label, "iterations": iterations,
                       "scaled_residual_inf": mp.nstr(scaled, 24),
                       "max_abs_surface_residual": mp.nstr(max_abs, 24),
                       "relative_residual": mp.nstr(relative, 24),
                       "site_balance_error": mp.nstr(site_error, 8),
                       "Ha_minus_OOH": mp.nstr(inventory, 8),
                       "accepted": bool(accepted_here), "error": error}
            attempts.append(attempt)
            print(json.dumps(attempt, ensure_ascii=False))
            if accepted_here:
                accepted_roots.append((root, cv, vacancy, net, rates, tof, attempt))

        (run_dir / "candidate_solver_attempts.json").write_text(
            json.dumps({"status": "ACCEPTED" if accepted_roots else "NO_ACCEPTED_ROOT",
                        "catmap_version": getattr(catmap, "__version__", "unknown"),
                        "conditions": {"temperature_K": "333.15", "concentrations": {
                            "C6H12_g": "1.0", "H2O2_g": "0.185", "H2O_g": "0.815", "C6H12O_g": "1e-20"}},
                        "inventory_constraint": "theta(Ha*) = theta(OOH*) = 0-inventory class",
                        "seed_note": seed_source, "attempts": attempts,
                        "accepted_roots": [{"seed": r[6]["seed"],
                                            "product_TOF": mp.nstr(r[5][solver.gas_names.index("C6H12O_g")], 40),
                                            "vacancy": mp.nstr(r[2], 30)} for r in accepted_roots]}, indent=2), encoding="utf-8")
        if not accepted_roots:
            print("NO_ACCEPTED_ROOT; no TOF or coverage accepted as candidate result")
            return

        root, cv, vacancy, net, rates, tof, accepted_attempt = accepted_roots[0]
        coverage_path = run_dir / "candidate_steady_state_coverages.csv"
        with coverage_path.open("w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=["species", "coverage", "n_sites", "site_occupancy"])
            w.writeheader()
            w.writerow({"species": "*", "coverage": mp.nstr(vacancy, 70), "n_sites": 1,
                        "site_occupancy": mp.nstr(vacancy, 70)})
            for name, value in zip(names, cv):
                n_sites = int(model.species_definitions[name].get("n_sites", 1))
                w.writerow({"species": name, "coverage": mp.nstr(value, 70), "n_sites": n_sites,
                            "site_occupancy": mp.nstr(value * n_sites, 70)})

        rate_path = run_dir / "candidate_reaction_rates.csv"
        with rate_path.open("w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=["reaction_index", "reaction_internal", "net_rate",
                                               "max_abs_surface_residual"])
            w.writeheader()
            for j, (reaction, rate) in enumerate(zip(solver.elementary_rxns, rates), 1):
                w.writerow({"reaction_index": j, "reaction_internal": model.print_rxn(reaction, mode="text"),
                            "net_rate": mp.nstr(rate, 70),
                            "max_abs_surface_residual": mp.nstr(max(abs(net[i]) for i in independent_indices), 24)})
        tof_path = run_dir / "candidate_tof.csv"
        with tof_path.open("w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=["gas_species", "turnover_frequency", "status"])
            w.writeheader()
            for gas, value in zip(solver.gas_names, tof):
                w.writerow({"gas_species": gas, "turnover_frequency": mp.nstr(value, 70),
                            "status": "PROVISIONAL_ENDPOINT_BARRIER_CANDIDATE"})
        print("accepted_root=" + json.dumps(accepted_attempt, ensure_ascii=False))
        print(f"coverages={coverage_path}")
        print(f"rates={rate_path}")
        print(f"tof={tof_path}")
    finally:
        os.chdir(previous_cwd)


if __name__ == "__main__":
    main()
