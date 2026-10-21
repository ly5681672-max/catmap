"""Run one provisional site with CatMAP 0.3.1 and stoichiometry-derived inventory."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import Counter
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from pathlib import Path

import mpmath as mp
from scipy.optimize import least_squares


def nullspace(A: list[list[int]]) -> list[list[int]]:
    """Exact integer basis of null(A), returned using rational row reduction."""
    M = [[Fraction(v) for v in row] for row in A]
    nrow, ncol = len(M), len(M[0]) if M else 0
    pivots: list[int] = []
    r = 0
    for c in range(ncol):
        pivot = next((i for i in range(r, nrow) if M[i][c]), None)
        if pivot is None:
            continue
        M[r], M[pivot] = M[pivot], M[r]
        scale = M[r][c]
        M[r] = [v / scale for v in M[r]]
        for i in range(nrow):
            if i != r and M[i][c]:
                factor = M[i][c]
                M[i] = [a - factor * b for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == nrow:
            break
    free = [c for c in range(ncol) if c not in pivots]
    basis = []
    for f in free:
        v = [Fraction(0) for _ in range(ncol)]
        v[f] = Fraction(1)
        for row, pivot_col in enumerate(pivots):
            v[pivot_col] = -M[row][f]
        den = reduce(lcm, (x.denominator for x in v), 1)
        ints = [int(x * den) for x in v]
        common = reduce(gcd, (abs(x) for x in ints if x), 0) or 1
        ints = [x // common for x in ints]
        first = next((x for x in ints if x), 1)
        if first < 0:
            ints = [-x for x in ints]
        basis.append(ints)
    return basis


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def solve(model_path: Path, precision: int) -> dict:
    model_path = model_path.resolve()
    out = model_path.parent
    old_cwd = Path.cwd()
    os.chdir(out)
    try:
        from catmap import ReactionModel
        from catmap.solvers.solver_base import NewtonRoot
        import catmap

        model = ReactionModel(setup_file=model_path.name)
        model.run()
        solver = model.solver
        scaler = model.scaler
        solver.compile()
        mp.mp.dps = precision
        names = list(solver.adsorbate_names)
        reactions = list(solver.elementary_rxns)
        n_rxn = len(reactions)
        if n_rxn != 9:
            raise RuntimeError(f"candidate reaction count expected 9, observed {n_rxn}")

        stoich = [[0] * n_rxn for _ in names]
        for j, rxn in enumerate(reactions):
            initial, final = Counter(rxn[0]), Counter(rxn[-1])
            for i, name in enumerate(names):
                stoich[i][j] = final.get(name, 0) - initial.get(name, 0)
        invariants = nullspace([list(row) for row in zip(*stoich)])
        invariant_rows = []
        for vec in invariants:
            invariant_rows.append({"invariant": " + ".join(f"({coef})theta({name})" for coef, name in zip(vec, names) if coef),
                                   **{name: coef for coef, name in zip(vec, names)}})
        write_csv(out / "inventory_invariants.csv", invariant_rows, ["invariant", *names])
        if len(invariants) != 1:
            raise RuntimeError(f"candidate zero-inventory constraint unsupported: nullity={len(invariants)}")
        inv = invariants[0]
        support = [(i, c) for i, c in enumerate(inv) if c]
        if len(support) != 2 or support[0][1] != -support[1][1] or abs(support[0][1]) != 1:
            raise RuntimeError(f"candidate inventory is not a simple pair equality: {inv}")
        dep_i, dep_c = support[0]
        partner_i, partner_c = support[1]
        if dep_c < 0:
            dep_i, partner_i, dep_c, partner_c = partner_i, dep_i, partner_c, dep_c
        dependent, partner = names[dep_i], names[partner_i]
        independent = [n for n in names if n != dependent]
        independent_idx = [names.index(n) for n in independent]

        site_n = {name: int(model.species_definitions.get(name, {}).get("n_sites", 1)) for name in names}
        if any(n != 1 for n in site_n.values()):
            raise RuntimeError(f"current positive-log parameterization requires n_sites=1; observed={site_n}")
        descriptors = [mp.mpf("333.15"), mp.mpf("1.0")]
        rxn_params = list(scaler.get_rxn_parameters(descriptors))

        def state(z):
            weights = [mp.exp(v) for v in z]
            weight_map = dict(zip(independent, weights))
            dep_weight = weight_map[partner]
            denom = 1 + sum(weights) + dep_weight
            cv = [mp.mpf("0") for _ in names]
            for name, weight in weight_map.items():
                cv[names.index(name)] = weight / denom
            cv[dep_i] = dep_weight / denom
            vacancy = 1 / denom
            return cv, vacancy

        def evaluate(z, full=False):
            cv, vacancy = state(z)
            k = list(solver.get_rate_constants(rxn_params, cv))
            net_species = list(solver.steady_state_function(cv))
            if not full:
                return cv, vacancy, net_species
            rates = list(solver.get_rxn_rates(cv, k))
            tof = list(solver.get_turnover_frequency(rxn_params, rates=rates, verify_coverages=False))
            return cv, vacancy, net_species, rates, tof

        seeds = [[mp.mpf(v)] * len(independent) for v in ("0", "-5", "3")]
        labels = [f"uniform_{v}" for v in ("0", "-5", "3")]
        attempts, accepted = [], []
        for label, start in zip(labels, seeds):
            cv0, _, directional = evaluate(start)
            solver._coverage = list(cv0)
            d_rates = list(solver.get_directional_rates(rxn_params))
            row_scales = []
            for i in independent_idx:
                gross = sum(abs(stoich[i][j]) * (abs(d_rates[j]) + abs(d_rates[n_rxn+j])) for j in range(n_rxn))
                row_scales.append(max(mp.mpf(str(gross)), mp.mpf("1e-300")))
            global_scale = max(row_scales)
            row_scales = [max(x, global_scale * mp.mpf("1e-250")) for x in row_scales]

            def residual(z):
                _, _, n = evaluate(z)
                return tuple(n[i] / s for i, s in zip(independent_idx, row_scales))

            def jac(z):
                z = list(z)
                # Differentiate the exact CatMAP steady-state function after
                # applying the inventory/site-balance parameterization. This
                # avoids assuming ordering or semantics for internal Jacobian
                # rows in the installed CatMAP build.
                return mp.matrix([[mp.diff(lambda q: residual([q if k == col else z[k]
                                                               for k in range(len(z))])[row], z[col])
                                   for col in range(len(z))] for row in range(len(z))])

            # CatMAP's NewtonRoot documentation recommends a low-cost initial
            # guess followed by high-precision polishing. SciPy is used only
            # to seed the official CatMAP NewtonRoot; acceptance is evaluated
            # exclusively after high-precision polishing below.
            scipy_message = ""
            precision_digits = mp.mp.dps
            try:
                mp.mp.dps = min(precision_digits, 50)
                seed_result = least_squares(
                    lambda x: [float(v) for v in residual([mp.mpf(str(a)) for a in x])],
                    [float(v) for v in start], method="trf", bounds=(-700.0, 700.0),
                    x_scale="jac", max_nfev=100, ftol=1e-10, xtol=1e-10, gtol=1e-10)
                if all(mp.isfinite(mp.mpf(str(v))) for v in seed_result.x):
                    polish_start = [mp.mpf(str(v)) for v in seed_result.x]
                    scipy_message = f"{seed_result.success}: {seed_result.message}"
                else:
                    polish_start = list(start)
                    scipy_message = f"nonfinite seed; fallback to {label}"
            except Exception as exc:
                polish_start = list(start)
                scipy_message = f"seed failed; fallback to {label}: {type(exc).__name__}: {exc}"
            finally:
                mp.mp.dps = precision_digits
            root, scaled, error, iterations = list(polish_start), mp.inf, None, 0
            try:
                iterator = NewtonRoot(residual, polish_start, mp.matrix, mp.mpf, mp.lu_solve, J=jac,
                                      norm=lambda x: mp.norm(x, "inf"), verbose=0)
                for iterations, (root_matrix, scaled) in enumerate(iterator, 1):
                    root = list(root_matrix)
                    if scaled < mp.mpf("1e-85") or iterations >= 300:
                        break
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
            cv, vacancy, net, rates, tof = evaluate(root, True)
            max_res = max(abs(net[i]) for i in independent_idx)
            max_rate = max(abs(x) for x in rates)
            relative = max_res / max(max_rate, mp.mpf("1e-300"))
            inventory_exact = sum(Fraction(c) * Fraction(str(cv[i])) for i, c in enumerate(inv))
            inventory_err = abs(mp.mpf(inventory_exact.numerator) / inventory_exact.denominator)
            site_err = abs(sum(cv) + vacancy - 1)
            accepted_here = (scaled < mp.mpf("1e-80") and relative < mp.mpf("1e-25") and
                             site_err < mp.mpf("1e-100") and inventory_err < mp.mpf("1e-100") and
                             vacancy >= 0 and all(x >= 0 for x in cv))
            attempt = {"seed": label, "accepted": bool(accepted_here), "iterations": iterations,
                       "scaled_residual_inf": mp.nstr(scaled, 24), "relative_surface_residual": mp.nstr(relative, 24),
                       "site_balance_error": mp.nstr(site_err, 10), "inventory_error": mp.nstr(inventory_err, 10),
                       "scipy_initial_guess": scipy_message, "solver_error": error}
            attempts.append(attempt)
            print(json.dumps(attempt, ensure_ascii=False), flush=True)
            if accepted_here:
                accepted.append((root, cv, vacancy, net, rates, tof, attempt))

        status = "ACCEPTED" if accepted else "NO_ACCEPTED_ROOT"
        inv_terms = [f"({coef})*theta({name})" for coef, name in zip(inv, names) if coef]
        inventory_equation = " + ".join(inv_terms) + " = 0"
        (out / "candidate_solver_attempts.json").write_text(json.dumps({
            "status": status, "catmap_version": getattr(catmap, "__version__", "unknown"),
            "constraint_derived_from_stoichiometric_nullspace": {"coefficients": dict(zip(names, inv)), "dependent_species": dependent, "paired_species": partner,
                                                                   "zero_inventory_equation": inventory_equation},
            "seed_policy": "fixed positive numerical seeds only; seeds are not physical initial states",
            "attempts": attempts,
            "accepted_roots": [{"seed": r[6]["seed"], "product_TOF": mp.nstr(r[5][solver.gas_names.index("C6H12O_g")], 90)} for r in accepted]}, indent=2), encoding="utf-8")
        if not accepted:
            raise RuntimeError(f"no candidate steady root accepted for {model_path}")

        # Different converged roots indicate branch ambiguity, even when each is
        # numerically valid. Do not silently select a branch as a unique result.
        cover_vectors = [[r[1][i] for i in range(len(names))] + [r[2]] for r in accepted]
        branch_spread = mp.mpf("0")
        for i in range(len(cover_vectors[0])):
            vals = [v[i] for v in cover_vectors]
            branch_spread = max(branch_spread, max(vals) - min(vals))
        if branch_spread > mp.mpf("1e-20"):
            raise RuntimeError(
                "multiple numerically accepted roots are not equivalent; "
                f"max absolute coverage spread={mp.nstr(branch_spread, 12)}"
            )
        chosen = accepted[0]
        root, cv, vacancy, net, rates, tof, attempt = chosen
        # Export CatMAP's own forward/reverse fluxes at this exact accepted
        # constrained steady state; do not reconstruct directions from net rates.
        solver._coverage = list(cv)
        solver.get_rate_constants(rxn_params, list(cv))
        directional_rates = list(solver.get_directional_rates(rxn_params))
        forward_rates = directional_rates[:n_rxn]
        reverse_rates = directional_rates[n_rxn:]
        cover_rows = [{"species": "*", "coverage": mp.nstr(vacancy, 80), "n_sites": 1, "site_occupancy": mp.nstr(vacancy, 80)}]
        cover_rows += [{"species": n, "coverage": mp.nstr(v, 80), "n_sites": site_n[n], "site_occupancy": mp.nstr(v * site_n[n], 80)} for n, v in zip(names, cv)]
        write_csv(out / "candidate_steady_state_coverages.csv", cover_rows)
        rate_rows = []
        for j, (rxn, rate) in enumerate(zip(reactions, rates), 1):
            rate_rows.append({"reaction_index": j, "reaction_internal": model.print_rxn(rxn, mode="text"),
                              "forward_rate": mp.nstr(forward_rates[j-1], 80),
                              "reverse_rate": mp.nstr(reverse_rates[j-1], 80),
                              "net_rate": mp.nstr(rate, 80), "status": "PROVISIONAL_CANDIDATE"})
        write_csv(out / "candidate_reaction_rates.csv", rate_rows)
        tof_rows = [{"gas_species": g, "turnover_frequency": mp.nstr(v, 80), "status": "PROVISIONAL_CANDIDATE"} for g, v in zip(solver.gas_names, tof)]
        write_csv(out / "candidate_tof.csv", tof_rows)

        stoich_rows = []
        for name, row in zip(names, stoich):
            stoich_rows.append({"stable_species": name, **{f"reaction_{j+1}": x for j, x in enumerate(row)}})
        write_csv(out / "candidate_stoichiometry.csv", stoich_rows)
        all_rates_match = max(abs(x - tof[solver.gas_names.index("C6H12O_g")]) for x in rates) / max(abs(tof[solver.gas_names.index("C6H12O_g")]), mp.mpf("1e-300"))
        summary = {"surface": model.surface_names[0], "status": "NUMERICAL_PASS_PROVISIONAL_PHYSICS",
                   "catmap_version": getattr(catmap, "__version__", "unknown"), "temperature_K": 333.15,
                   "candidate_T_formula_s-1": mp.nstr(tof[solver.gas_names.index("C6H12O_g")], 80),
                   "inventory_invariant": f"{dependent} - {partner} = 0 (derived from this candidate's CatMAP stoichiometry)",
                   "accepted_roots": len(accepted), "chosen_seed": attempt["seed"],
                   "accepted_root_max_absolute_coverage_spread": mp.nstr(branch_spread, 12),
                   "max_relative_surface_residual": attempt["relative_surface_residual"],
                   "site_balance_error": attempt["site_balance_error"], "inventory_error": attempt["inventory_error"],
                   "max_net_reaction_rate_vs_product_TOF_relative_difference": mp.nstr(all_rates_match, 12),
                   "physical_status": "PROVISIONAL; endpoint approximations, one-effective-site representation, and reduced-state mapping are unverified"}
        (out / "numerical_gate.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return summary
    finally:
        os.chdir(old_cwd)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--precision", type=int, default=220)
    args = ap.parse_args()
    solve(args.model, args.precision)
