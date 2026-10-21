"""One bounded provisional Ti-Co solve for the explicit full-state network.

CatMAP 0.3.1 supplies the parsed mechanism, thermodynamic rate constants,
surface-rate function and NewtonRoot iteration. A positive log-coverage chart
enforces the two zero-inventory relations and the formal multi-site balance.
The composite-state n_sites values and endpoint-link kinetics remain explicit
conditional assumptions; this is not a physically validated prediction.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
from collections import Counter
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import least_squares


HERE = Path(__file__).resolve().parents[1]
MODEL = HERE / "inputs" / "FULL_8STEP_MODEL.mkm"
RUN = HERE / "runs" / "tico_provisional_zero_inventory"
NSITES = {
    "C6H12O_s": 1, "C6H12_s": 1, "H2O2_s": 1, "H2O_s": 1,
    "Ha_s": 1, "Hb_s": 1, "NEB1FS_s": 4, "NEB1IS_s": 3,
    "NEB2FS_s": 2, "NEB2IS_s": 3, "OH_s": 1, "OOH_s": 1, "O_s": 1,
}
OMIT = {"Ha_s", "Hb_s"}


def save_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--precision", type=int, default=100)
    ap.add_argument("--max-iterations", type=int, default=250)
    args = ap.parse_args()
    RUN.mkdir(parents=True, exist_ok=True)
    log_path = RUN / "solver_console.jsonl"
    started = time.time()
    old_cwd = Path.cwd()
    os.chdir(MODEL.parent)
    try:
        from catmap import ReactionModel
        from catmap.solvers.solver_base import NewtonRoot
        import catmap

        model = ReactionModel(setup_file=MODEL.name)
        # Native run initializes CatMAP's parser/scaler/solver. With no
        # descriptor_values/ranges in this setup it does not map a TOF grid.
        model.run()
        solver, scaler = model.solver, model.scaler
        solver.compile()
        mp.mp.dps = args.precision
        names = list(solver.adsorbate_names)
        if set(names) != set(NSITES):
            raise RuntimeError(f"stable species mismatch: parsed={names}; expected={sorted(NSITES)}")
        ha_i, hb_i, o_i = names.index("Ha_s"), names.index("Hb_s"), names.index("O_s")
        oh_i, ooh_i = names.index("OH_s"), names.index("OOH_s")
        independent = [n for n in names if n not in OMIT]
        indices = [names.index(n) for n in independent]
        params = list(scaler.get_rxn_parameters([mp.mpf("333.15"), mp.mpf("1.0")]))
        n_rxn = len(solver.elementary_rxns)
        if n_rxn != 12:
            raise RuntimeError(f"expected 12 CatMAP elementary reactions, got {n_rxn}")

        # The parsed invariant basis is Ha-Hb-OH-OOH=0 and Hb-O=0.
        # Thus Hb=O and Ha=O+OH+OOH, preserving the requested original basis.
        def coverage_from_logweights(z):
            weights = {n: mp.exp(v) for n, v in zip(independent, z)}
            w_o, w_oh, w_ooh = weights["O_s"], weights["OH_s"], weights["OOH_s"]
            denominator = 1 + sum(NSITES[n] * weights[n] for n in independent)
            denominator += NSITES["Hb_s"] * w_o
            denominator += NSITES["Ha_s"] * (w_o + w_oh + w_ooh)
            vacancy = 1 / denominator
            cv = [mp.mpf("0")] * len(names)
            for n, w in weights.items():
                cv[names.index(n)] = w * vacancy
            cv[hb_i] = w_o * vacancy
            cv[ha_i] = (w_o + w_oh + w_ooh) * vacancy
            return cv, vacancy

        def calc(z, with_rates=False):
            cv, vacancy = coverage_from_logweights(z)
            k = list(solver.get_rate_constants(params, cv))
            residual = list(solver.steady_state_function(cv))
            if not with_rates:
                return cv, vacancy, residual
            net_rates = list(solver.get_rxn_rates(cv, k))
            tof = list(solver.get_turnover_frequency(params, rates=net_rates, verify_coverages=False))
            return cv, vacancy, residual, net_rates, tof

        # Check thermodynamic detailed balance from CatMAP's native rate
        # constants against the same CatMAP-parsed standard reaction energies.
        probe, _ = coverage_from_logweights([mp.mpf("-10")] * len(independent))
        solver._coverage = list(probe)
        k_all = list(solver.get_rate_constants(params, probe))
        kf, kr = k_all[:n_rxn], k_all[n_rxn:]
        def energy(state):
            total = mp.mpf("0")
            for ident in state:
                if ident in model.site_names:
                    continue
                e = model.species_definitions[ident]["formation_energy"]
                total += mp.mpf(str(e[0] if isinstance(e, (list, tuple)) else e))
            return total
        thermo_rows = []
        kbt = mp.mpf("8.617332478e-5") * mp.mpf("333.15")
        for j, reaction in enumerate(solver.elementary_rxns):
            dg = energy(list(reaction[-1])) - energy(list(reaction[0]))
            k_ratio = kf[j] / kr[j]
            expected = mp.exp(-dg / kbt)
            rel_err = abs(k_ratio - expected) / max(abs(expected), mp.mpf("1e-100"))
            thermo_rows.append({"reaction_index": j + 1,
                                "reaction": model.print_rxn(reaction, mode="text"),
                                "deltaG_eV": mp.nstr(dg, 35),
                                "kf_over_kr": mp.nstr(k_ratio, 35),
                                "exp_minus_deltaG_over_kBT": mp.nstr(expected, 35),
                                "relative_error": mp.nstr(rel_err, 20),
                                "status": "PASS" if rel_err < mp.mpf("1e-12") else "FAIL"})
        save_csv(RUN / "thermodynamic_detailed_balance.csv", list(thermo_rows[0]), thermo_rows)

        stoich = [[0] * n_rxn for _ in names]
        for j, rxn in enumerate(solver.elementary_rxns):
            left, right = Counter(rxn[0]), Counter(rxn[-1])
            for i, n in enumerate(names):
                stoich[i][j] = right.get(n, 0) - left.get(n, 0)

        results, accepted = [], []
        seeds: list[tuple[str, list[mp.mpf]]] = []
        legacy_seed = HERE.parents[0] / "candidate_tico_kinetic_control" / "candidate_steady_state_coverages.csv"
        if legacy_seed.is_file():
            old = {}
            with legacy_seed.open(encoding="utf-8-sig", newline="") as f:
                for row in csv.DictReader(f):
                    old[row["species"]] = mp.mpf(row["coverage"])
            seed_values = {n: max(old.get(n, mp.mpf("1e-40")), mp.mpf("1e-40")) for n in independent}
            # New O/OH/Hb coordinates are numerical seed values only. Project
            # Ha and Hb from the two exact zero-inventory relations below.
            seed_values["O_s"] = mp.mpf("1e-30")
            seed_values["OH_s"] = mp.mpf("1e-30")
            weighted = sum(NSITES[n] * v for n, v in seed_values.items())
            weighted += NSITES["Hb_s"] * seed_values["O_s"]
            weighted += NSITES["Ha_s"] * (3 * seed_values["O_s"] + seed_values["OH_s"] + seed_values["OOH_s"])
            v0 = max(mp.mpf("1") - weighted, mp.mpf("1e-40"))
            seeds.append(("legacy_candidate_root_projected_numerical_seed_only",
                          [mp.log(seed_values[n] / v0) for n in independent]))
        seeds.extend((f"uniform_log_{level}", [mp.mpf(level)] * len(independent))
                     for level in ("-20", "-8", "-2"))
        with log_path.open("w", encoding="utf-8") as logf:
            for label, start in seeds:
                cv0, _, _ = calc(start)
                solver._coverage = list(cv0)
                directional = list(solver.get_directional_rates(params))
                scales = []
                for i in indices:
                    gross = sum(abs(stoich[i][j]) *
                                (abs(directional[j]) + abs(directional[n_rxn + j]))
                                for j in range(n_rxn))
                    scales.append(max(gross, mp.mpf("1e-300")))
                gmax = max(scales)
                scales = [max(s, gmax * mp.mpf("1e-250")) for s in scales]

                def f(z):
                    cv, _, net = calc(z)
                    return tuple(net[i] / scale for i, scale in zip(indices, scales))

                def jac(z):
                    z = list(z)
                    return mp.matrix([[mp.diff(lambda q: f([q if k == col else z[k]
                                                           for k in range(len(z))])[row], z[col])
                                       for col in range(len(z))] for row in range(len(z))])

                # CatMAP's NewtonRoot documentation recommends a SciPy
                # low-precision seed followed by high-precision root polishing.
                # The kinetic equations and all rates remain CatMAP-native.
                try:
                    seed_fit = least_squares(
                        lambda x: np.asarray([float(v) for v in f([mp.mpf(str(v)) for v in x])]),
                        np.asarray([float(v) for v in start]),
                        bounds=(-200.0, 20.0), max_nfev=500, xtol=1e-12, ftol=1e-12, gtol=1e-12,
                    )
                    start_refined = [mp.mpf(str(v)) for v in seed_fit.x]
                except Exception:
                    seed_fit = None
                    start_refined = list(start)
                root, norm, error, iteration = list(start_refined), mp.inf, None, 0
                try:
                    iterator = NewtonRoot(f, start_refined, mp.matrix, mp.mpf, mp.lu_solve,
                                          J=jac, norm=lambda x: mp.norm(x, "inf"), verbose=0)
                    for iteration, (candidate, norm) in enumerate(iterator, 1):
                        root = list(candidate)
                        if norm < mp.mpf("1e-60") or iteration >= args.max_iterations:
                            break
                except Exception as exc:
                    error = f"{type(exc).__name__}: {exc}"

                cv, vacancy, net, rates, tof = calc(root, True)
                max_abs = max(abs(net[i]) for i in range(len(net)))
                solver._coverage = list(cv)
                directional = list(solver.get_directional_rates(params))
                gross_max = max((abs(r) for r in directional), default=mp.mpf("0"))
                rel = max_abs / max(gross_max, mp.mpf("1e-300"))
                site_sum = vacancy + sum(NSITES[n] * cv[i] for i, n in enumerate(names))
                i1 = cv[ha_i] - cv[ooh_i] - cv[oh_i] - cv[o_i]
                i2 = cv[hb_i] - cv[o_i]
                ok = (norm < mp.mpf("1e-50") and rel < mp.mpf("1e-25")
                      and abs(site_sum - 1) < mp.mpf("1e-70")
                      and abs(i1) < mp.mpf("1e-70") and abs(i2) < mp.mpf("1e-70")
                      and vacancy > 0 and all(v >= 0 for v in cv))
                row = {"seed": label, "iterations": iteration,
                       "scipy_seed_success": bool(seed_fit.success) if seed_fit is not None else False,
                       "scipy_seed_cost": float(seed_fit.cost) if seed_fit is not None else "",
                       "scaled_residual_inf": mp.nstr(norm, 25),
                       "max_abs_surface_residual": mp.nstr(max_abs, 25),
                       "relative_residual": mp.nstr(rel, 25),
                       "weighted_site_sum": mp.nstr(site_sum, 25),
                       "I1_Ha_OOH_OH_O": mp.nstr(i1, 12), "I2_Hb_O": mp.nstr(i2, 12),
                       "accepted": str(bool(ok)), "error": error or ""}
                results.append(row)
                logf.write(json.dumps(row) + "\n")
                logf.flush()
                print(json.dumps(row, ensure_ascii=False))
                if ok:
                    accepted.append((cv, vacancy, net, rates, tof, row))
        save_csv(RUN / "solver_attempts.csv", list(results[0]), results)
        summary = {"status": "ACCEPTED_PROVISIONAL" if accepted else "NO_ACCEPTED_ROOT",
                   "model_class": "CONDITIONAL_ENDPOINT_BARRIER_CANDIDATE",
                   "catmap_version": getattr(catmap, "__version__", "unknown"),
                   "precision_digits": args.precision, "wall_seconds": round(time.time()-started, 3),
                   "conditions": {"temperature_K": 333.15, "concentrations":
                       {"C6H12_g": 1.0, "H2O2_g": 0.185, "H2O_g": 0.815, "C6H12O_g": 1e-20}},
                   "inventory_constraints": {"I1": "theta(Ha)-theta(OOH)-theta(OH)-theta(O)=0",
                                             "I2": "theta(Hb)-theta(O)=0"},
                   "site_counts": NSITES, "attempts": results}
        summary["detailed_balance"] = {
            "checked_reactions": len(thermo_rows),
            "status": "PASS" if all(r["status"] == "PASS" for r in thermo_rows) else "FAIL",
            "max_relative_error": max(float(r["relative_error"]) for r in thermo_rows),
            "note": "native CatMAP kf/kr compared with exp(-standard_deltaG/kBT); gas concentrations remain in rate law, not standard K"
        }
        (RUN / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        if not accepted:
            return 2
        cv, vacancy, net, rates, tof, chosen = accepted[0]
        save_csv(RUN / "steady_state_coverages.csv",
                 ["species", "coverage", "n_sites", "site_occupancy"],
                 [{"species": "*", "coverage": mp.nstr(vacancy, 80), "n_sites": 1,
                   "site_occupancy": mp.nstr(vacancy, 80)}] +
                 [{"species": n, "coverage": mp.nstr(cv[i], 80), "n_sites": NSITES[n],
                   "site_occupancy": mp.nstr(cv[i] * NSITES[n], 80)} for i, n in enumerate(names)])
        save_csv(RUN / "reaction_net_rates.csv",
                 ["reaction_index", "reaction", "net_rate"],
                 [{"reaction_index": j, "reaction": model.print_rxn(r, mode="text"),
                   "net_rate": mp.nstr(rates[j-1], 80)} for j, r in enumerate(solver.elementary_rxns, 1)])
        save_csv(RUN / "gas_tof.csv", ["gas_species", "turnover_frequency", "status"],
                 [{"gas_species": n, "turnover_frequency": mp.nstr(v, 80),
                   "status": "PROVISIONAL_ENDPOINT_BARRIER_CANDIDATE"}
                  for n, v in zip(solver.gas_names, tof)])
        print("Accepted one provisional root; see", RUN)
        return 0
    finally:
        os.chdir(old_cwd)


if __name__ == "__main__":
    raise SystemExit(main())
