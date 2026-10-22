from __future__ import annotations

"""Bounded seed location + CatMAP NewtonRoot validation for Ti-W/Ti-Ta.

Only numerical coordinates/seeds are altered. Official CatMAP 0.3.1 rates,
steady-state equations, original eight-step model, energies, reservoirs, and
zero-inventory chart are retained.
"""

import importlib.util
import json
import math
import os
import csv
from collections import Counter
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import least_squares

PROJECT = Path(r"H:\catmap\catmap\catmap_code_collection_20261013")
BASE = PROJECT / "analyses" / "explicit_OOH_8step_validation"
REPLAY = BASE / "multisite_conditional_replay_20261021"
OLD_RECOVERY = BASE / "three_site_root_recovery_20261021"
OUT = Path(os.environ.get("TIW_TITA_OUTPUT_DIR", str(BASE / "tiw_tita_targeted_solver_recovery")))
OUT.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location(
    "three_site_recovery", OLD_RECOVERY / "scripts" / "three_site_recovery.py"
)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)

FREE = ("C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "OH_s", "OOH_s", "O_s")
WEIGHTS = (1, 1, 1, 1, 2, 2, 3)
SEEDS = ("tico", "tihf", "tini", "tire", "timn", "tiv", "tife", "ticr", "titi", "ti", "timo", "tizr")
MAX_SEEDS = 6
SCIPY_NFEV = 160
NEWTON_STEPS = 30


def dump_json(path: Path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def setup(site, dps, tag="solve"):
    work = OUT / "runs" / f"{site}_{tag}"
    work.mkdir(parents=True, exist_ok=True)
    model, solver, params, oldcwd = old.setup_model(site, work)
    mp.mp.dps = dps
    # Preserve the exact binary-float reservoir values while preventing
    # NumPy/Python float arithmetic from collapsing high-precision net fluxes.
    # This is a dtype promotion only: each value is converted from the actual
    # runtime float, not re-entered from a different decimal concentration.
    solver.gas_pressures = [mp.mpf(float(v)) for v in solver.gas_pressures]
    return model, solver, params, oldcwd


def make_context(solver, params):
    names = list(solver.adsorbate_names)
    idx = {name: i for i, name in enumerate(names)}
    free_idx = [idx[n] for n in FREE]
    row_names = FREE
    rows = [idx[n] for n in row_names]
    n_rxn = len(solver.elementary_rxns)
    nu = []
    for rxn in solver.elementary_rxns:
        ci, cf = Counter(rxn[0]), Counter(rxn[-1])
        nu.append([cf.get(n, 0) - ci.get(n, 0) for n in names])

    # Stable log-sum-exp chart for exact site balance and I1=I2=0.
    logw = [mp.log(mp.mpf(w)) for w in WEIGHTS]

    def state(z):
        terms = [mp.mpf(0)] + [z[i] + logw[i] for i in range(7)]
        m = max(terms)
        den = sum(mp.exp(t - m) for t in terms)
        vac = mp.exp(-m) / den
        cv = [mp.mpf(0)] * len(names)
        for i, j in enumerate(free_idx):
            cv[j] = mp.exp(z[i] - m) / den
        cv[idx["Ha_s"]] = cv[idx["OOH_s"]] + cv[idx["OH_s"]] + cv[idx["O_s"]]
        cv[idx["Hb_s"]] = cv[idx["O_s"]]
        return cv, vac

    def get_fluxes(cv):
        ks = list(solver.get_rate_constants(params, list(cv)))
        p = [mp.mpf(v) for v in solver.gas_pressures]
        # Mirror CatMAP 0.3.1's generated mean-field expression exactly:
        # s[0] = 1 - theta[0] - theta[1] - ... .  In particular, do not
        # substitute the algebraically equivalent softmax vacancy here: at
        # near-zero vacancy, the order/precision of subtraction is material.
        empty = [mp.mpf("1.0")]
        for value in cv:
            empty[0] -= value
        empty.append(mp.mpf("0.0"))
        env = {"kf": ks[:n_rxn], "kr": ks[n_rxn:], "theta": cv, "p": p, "s": empty}
        fwd, rev = [], []
        for j, rxn in enumerate(solver.elementary_rxns):
            fexpr = solver.rate_equation_term(rxn[0], f"kf[{j}]")
            rexpr = solver.rate_equation_term(rxn[-1], f"kr[{j}]")
            fwd.append(eval(fexpr, {"__builtins__": {}}, env))
            rev.append(eval(rexpr, {"__builtins__": {}}, env))
        return ks, fwd, rev

    def residual_mp(z, return_parts=False):
        cv, vac = state(z)
        ks, fwd, rev = get_fluxes(cv)
        net = [fwd[j] - rev[j] for j in range(n_rxn)]
        # CatMAP's official steady_state_function is the governing equation.
        f = list(solver.steady_state_function(list(cv)))
        production, consumption, gross, local, logratio = [], [], [], [], []
        for row in rows:
            prod = mp.mpf(0)
            cons = mp.mpf(0)
            for j in range(n_rxn):
                coeff = nu[j][row]
                if coeff > 0:
                    prod += coeff * abs(fwd[j])
                    cons += coeff * abs(rev[j])
                elif coeff < 0:
                    prod += -coeff * abs(rev[j])
                    cons += -coeff * abs(fwd[j])
            g = prod + cons
            production.append(prod)
            consumption.append(cons)
            gross.append(g)
            local.append((prod - cons) / g if g else mp.mpf(0))
            # Log production/consumption ratio avoids the +/-1 saturation of
            # normalized net residuals far from a root; it has exactly the same
            # zero as CatMAP's species balance equation.
            logratio.append(mp.log(prod / cons) if prod > 0 and cons > 0 else mp.nan)
        if return_parts:
            return cv, vac, ks, fwd, rev, net, f, production, consumption, gross, local, logratio
        return logratio

    def residual_float(zf):
        z = [mp.mpf(float(v)) for v in zf]
        try:
            vals = residual_mp(z)
            out = np.asarray([float(v) for v in vals], dtype=float)
            if not np.all(np.isfinite(out)):
                return np.full(len(rows), 1e100)
            return out
        except (ValueError, OverflowError, ZeroDivisionError):
            return np.full(len(rows), 1e100)

    return names, idx, free_idx, rows, nu, state, residual_mp, residual_float, get_fluxes


def seed_vectors(site, idx, state):
    seeds = []
    for src in SEEDS:
        path = REPLAY / "sites" / src / f"constrained_{src}_result.json"
        if not path.exists():
            continue
        d = json.loads(path.read_text(encoding="utf-8"))
        try:
            # Archived successful-root schema names the empty-site field
            # ``vacant`` (not ``vac``); preserve exact decimal text.
            vac = mp.mpf(d.get("vacant", d.get("vac")))
            c = d["coverage"]
            z = [mp.log(mp.mpf(c[n]) / vac) for n in FREE]
            seeds.append((src, z))
        except Exception:
            continue
    # Use documented failed-site coordinates only as extra initial guesses, never as roots.
    prev = json.loads((OLD_RECOVERY / "data" / f"{site}_recovery.json").read_text(encoding="utf-8"))
    for i, a in enumerate(prev.get("attempts", [])):
        z = a.get("root_logits")
        if isinstance(z, list) and len(z) == 7:
            try:
                zz = [mp.mpf(q) for q in z]
                if max(abs(q) for q in zz) < 100:
                    seeds.append((f"old_candidate_{i}_{a.get('seed','x')}", zz))
            except Exception:
                pass
    # Perturb only coordinates of initial guesses; exact inventory/site chart retained.
    base = list(seeds)
    patterns = [
        [0, 0, 0, 0, 0, 0, 0],
        [2, -2, 1, -1, 2, -2, 1],
        [-2, 2, -1, 1, -2, 2, -1],
        [1, 1, -2, 2, -1, 1, -2],
    ]
    for src, z in base[:4]:
        for k, pat in enumerate(patterns[1:3]):
            seeds.append((f"{src}_perturb{k}", [q + mp.mpf(pat[j]) for j, q in enumerate(z)]))
    for value in (-8, -4, -2, 0):
        seeds.append((f"uniform_{value}", [mp.mpf(value)] * 7))
    # Rank candidates at model-specific rates and preserve distinct seed profiles.
    unique, seen = [], set()
    for name, z in seeds:
        key = tuple(round(float(v), 6) for v in z)
        if key not in seen:
            unique.append((name, z))
            seen.add(key)
    return unique


def scipy_locate(site, solver, params, context):
    names, idx, free_idx, rows, nu, state, rmp, rf, get_fluxes = context
    seeds = seed_vectors(site, idx, state)
    # Rank all existing seeds by the native local-flow residual; process best diverse seeds.
    ranked = []
    for name, z in seeds:
        try:
            score = max(abs(float(v)) for v in rmp(z))
            if math.isfinite(score):
                ranked.append((score, name, z))
        except Exception:
            pass
    ranked.sort(key=lambda x: x[0])
    # Include seeds with chemically distinct dominant coverage profiles, then fill by score.
    picked, doms = [], set()
    for score, name, z in ranked:
        cv, _ = state(z)
        dom = max(FREE, key=lambda n: cv[idx[n]])
        if dom not in doms:
            picked.append((score, name, z))
            doms.add(dom)
        if len(picked) >= MAX_SEEDS:
            break
    for item in ranked:
        if len(picked) >= MAX_SEEDS:
            break
        if all(item[1] != q[1] for q in picked):
            picked.append(item)

    attempts = []
    candidates = []
    for score, name, z in picked:
        z0 = np.asarray([float(v) for v in z], dtype=float)
        # Numerical search box only bounds this seed locator; no boundary candidate passes.
        res = least_squares(rf, z0, method="trf", bounds=(-300.0, 300.0),
                            x_scale="jac", max_nfev=SCIPY_NFEV,
                            ftol=1e-12, xtol=1e-12, gtol=1e-12)
        zf = [mp.mpf(float(v)) for v in res.x]
        val = rmp(zf, return_parts=True)
        local_max = max(abs(float(v)) for v in val[-2])
        logratio_max = max(abs(float(v)) for v in val[-1])
        record = {"seed": name, "initial_local_residual": score, "success": bool(res.success),
                  "message": str(res.message), "nfev": int(res.nfev), "cost": float(res.cost),
                  "optimality": float(res.optimality), "final_local_residual_max": local_max,
                  "final_log_production_consumption_ratio_max": logratio_max,
                  "logits": [mp.nstr(v, 30) for v in zf],
                  "vacancy": mp.nstr(val[1], 25),
                  "dominant": max(FREE, key=lambda n: val[0][idx[n]])}
        attempts.append(record)
        candidates.append((local_max, name, zf))
    candidates.sort(key=lambda x: x[0])
    return attempts, candidates


def native_newton(site, solver, params, context, candidates):
    names, idx, free_idx, rows, nu, state, rmp, rf, get_fluxes = context
    from catmap.solvers.solver_base import NewtonRoot
    n_rxn = len(solver.elementary_rxns)

    results = []
    for score, seed_name, z0 in candidates[:4]:
        precision_levels = (260, 300) if site == "tiw" else (120, 180)
        for dps in precision_levels:
            mp.mp.dps = dps
            z0 = [mp.mpf(v) for v in z0]
            def F(z):
                return tuple(rmp(z))
            def J(z):
                z = list(z)
                return mp.matrix([[mp.diff(lambda q: F([q if k == col else z[k] for k in range(7)])[row], z[col])
                                   for col in range(7)] for row in range(7)])
            iterator = NewtonRoot(F, z0, mp.matrix, mp.mpf, mp.lu_solve,
                                  J=J, norm=lambda v: mp.norm(v, "inf"), verbose=0)
            z = z0
            err = None
            norm = mp.inf
            steps = 0
            polish_target = mp.mpf("1e-220") if site == "tiw" else mp.mpf("1e-100")
            try:
                for steps, (x, norm) in enumerate(iterator, start=1):
                    z = list(x)
                    if norm < polish_target or steps >= NEWTON_STEPS:
                        break
            except Exception as e:
                err = f"{type(e).__name__}: {e}"
            vals = rmp(z, return_parts=True)
            cv, vac, ks, fwd, rev, net, f, production, consumption, gross, local, logratio = vals
            official_net = list(solver.get_rxn_rates(list(cv), ks))
            native_elementary = list(solver.elementary_rates(ks, list(cv), solver.gas_pressures,
                                                              solver._mpfloat, solver._matrix))
            stoich_from_official = [sum(nu[j][row] * official_net[j] for j in range(n_rxn)) for row in rows]
            catmap_empty = mp.mpf(1) - sum(cv)
            manual_empty = mp.mpf(1) - sum(cv[i] * mp.mpf(str(solver.species_definitions[names[i]]["n_sites"]))
                                      for i in range(len(names)))
            direct_net0 = ks[0] * mp.mpf(str(solver.gas_pressures[2])) * catmap_empty - ks[n_rxn] * cv[idx["H2O2_s"]]
            rate_reconstruction_rel = max(
                (abs(official_net[j] - net[j]) / max(abs(official_net[j]), abs(net[j]), mp.mpf("1e-500"))
                 for j in range(n_rxn)), default=mp.mpf(0))
            balance_reconstruction_rel = max(
                (abs(f[row] - stoich_from_official[k]) /
                 max(abs(f[row]), abs(stoich_from_official[k]), mp.mpf("1e-500"))
                 for k, row in enumerate(rows)), default=mp.mpf(0))
            loc = [abs(v) for v in local]
            # Independent diagnostics: the seven chart coordinates plus the
            # two dependent species reconstructed from I1=I2=0.
            i1 = cv[idx["Ha_s"]] - cv[idx["OOH_s"]] - cv[idx["OH_s"]] - cv[idx["O_s"]]
            i2 = cv[idx["Hb_s"]] - cv[idx["O_s"]]
            full_site_sum = sum(cv) + catmap_empty
            site_error = full_site_sum - 1
            cycle_rel = max(
                (abs(official_net[j] - official_net[4]) /
                 max(abs(official_net[4]), abs(official_net[j]), mp.mpf("1e-1000"))
                 for j in range(n_rxn)), default=mp.mpf(0))
            cancel = [((abs(fwd[j]) + abs(rev[j])) / abs(net[j]) if net[j] else mp.inf)
                      for j in range(n_rxn)]
            try:
                singular_values = mp.svd(J(z), compute_uv=False)
                jac_smax = max(singular_values)
                jac_smin = min(singular_values)
                jac_cond = jac_smax / jac_smin if jac_smin else mp.inf
            except Exception:
                jac_smax = jac_smin = jac_cond = mp.nan
            record = {"seed": seed_name, "dps": dps, "steps": steps,
                      "polish_target": mp.nstr(polish_target, 10),
                      "newton_norm": mp.nstr(norm, 25), "error": err,
                      "max_local_residual": mp.nstr(max(loc), 25),
                      "max_log_production_consumption_ratio": mp.nstr(max(abs(v) for v in logratio), 25),
                      "vacancy": mp.nstr(vac, 35),
                      # Preserve enough significant digits to reconstruct
                      # cancellation-sensitive roots after JSON serialization.
                      "root_logits": [mp.nstr(v, max(100, dps - 10)) for v in z],
                      "coverage": {n: mp.nstr(cv[i], max(100, dps - 10)) for i, n in enumerate(names)},
                      "net_rates": [mp.nstr(v, max(100, dps - 10)) for v in net],
                      "official_net_rates": [mp.nstr(v, max(100, dps - 10)) for v in official_net],
                      "native_elementary_rates": [mp.nstr(v, max(100, dps - 10)) for v in native_elementary],
                      "manual_forward_rates": [mp.nstr(v, max(100, dps - 10)) for v in fwd],
                      "manual_reverse_rates": [mp.nstr(v, max(100, dps - 10)) for v in rev],
                      "empty_catmap_formula": mp.nstr(catmap_empty, max(100, dps - 10)),
                      "empty_manual_formula": mp.nstr(manual_empty, max(100, dps - 10)),
                      "direct_formula_step1": mp.nstr(direct_net0, max(100, dps - 10)),
                      "k_first_forward_return": mp.nstr(ks[0], 40),
                      "k_first_reverse_return": mp.nstr(ks[n_rxn], 40),
                      "k_first_forward_internal": mp.nstr(solver._kf[0], 40),
                      "k_first_reverse_internal": mp.nstr(solver._kr[0], 40),
                      "pressure_H2O2_value": mp.nstr(mp.mpf(solver.gas_pressures[2]), 40),
                      "theta_H2O2_value": mp.nstr(cv[idx["H2O2_s"]], 40),
                      "official_balance_from_net": [mp.nstr(v, 40) for v in stoich_from_official],
                      "catmap_steady_state_function": [mp.nstr(v, 50) for v in f],
                      "local_residual_species": {names[row]: mp.nstr(local[k], 40)
                                                  for k, row in enumerate(rows)},
                      "gross_local_flux_species": {names[row]: mp.nstr(gross[k], 40)
                                                    for k, row in enumerate(rows)},
                      "production_local_flux_species": {names[row]: mp.nstr(production[k], 40)
                                                         for k, row in enumerate(rows)},
                      "consumption_local_flux_species": {names[row]: mp.nstr(consumption[k], 40)
                                                         for k, row in enumerate(rows)},
                      "inventory_I1": mp.nstr(i1, 50), "inventory_I2": mp.nstr(i2, 50),
                      "full_site_sum": mp.nstr(full_site_sum, 50),
                      "site_balance_error": mp.nstr(site_error, 50),
                      "max_cycle_flux_relative_deviation_vs_product_TOF": mp.nstr(cycle_rel, 30),
                      "forward_reverse_cancellation_ratio": [mp.nstr(v, 40) for v in cancel],
                      "logratio_jacobian_sigma_max": mp.nstr(jac_smax, 30),
                      "logratio_jacobian_sigma_min": mp.nstr(jac_smin, 30),
                      "logratio_jacobian_condition_number": mp.nstr(jac_cond, 30),
                      "max_directional_reconstruction_relative_error": mp.nstr(rate_reconstruction_rel, 25),
                      "max_balance_reconstruction_relative_error": mp.nstr(balance_reconstruction_rel, 25)}
            if (max(loc) < mp.mpf("1e-12") and vac > 0 and all(v >= 0 for v in cv)
                    and abs(site_error) < mp.mpf("1e-100") and abs(i1) < mp.mpf("1e-100")
                    and abs(i2) < mp.mpf("1e-100") and cycle_rel < mp.mpf("1e-12")
                    and rate_reconstruction_rel < mp.mpf("1e-12")):
                tof = net[4]  # C6H12O* desorption is elementary step 5
                record["tof_s-1"] = mp.nstr(tof, 90)
                record["log10_tof"] = mp.nstr(mp.log10(tof), 60) if tof > 0 else None
                record["accepted_native_equation"] = True
                results.append((record, z))
            else:
                record["accepted_native_equation"] = False
                results.append((record, z))
    return results


def run(site):
    # A fresh CatMAP cache is essential after precision/math-context changes.
    # Keep earlier diagnostic runs untouched instead of reusing their pickle.
    work = OUT / "runs" / f"{site}_final_20261022"
    # Build log(2), log(3), and the exact-inventory chart at the target
    # precision. Raising mpmath precision after these constants are created
    # cannot restore their discarded digits.
    context_dps = 300 if site == "tiw" else 180
    model, solver, params, oldcwd = setup(site, context_dps, tag="final_root_20261022")
    try:
        context = make_context(solver, params)
        attempts, candidates = scipy_locate(site, solver, params, context)
        if not candidates:
            # If the bounded float seed locator rejects every historical seed,
            # recover the already archived high-precision root from the prior
            # result table. Reconstruct vacancy from the exact weighted site
            # identity (not the rounded vacancy text column) and replay this
            # state through CatMAP's native equations before Newton polishing.
            table = OUT / "deliverables" / "data" / "TIW_TITA_SOLVER_RESULTS.csv"
            if table.exists():
                with table.open(encoding="utf-8-sig", newline="") as handle:
                    row = next(r for r in csv.DictReader(handle) if r["surface"] == site)
                _, idx, _, _, _, _, residual, _, _ = context
                coverages = [mp.mpf(row[f"theta_{name}"]) for name in FREE]
                vacancy = mp.mpf(1) - sum(mp.mpf(w) * c
                                          for w, c in zip(WEIGHTS, coverages))
                z = [mp.log(c / vacancy) for c in coverages]
                replay_score = max(abs(float(v)) for v in residual(z))
                if math.isfinite(replay_score) and replay_score < 1e-12:
                    candidates = [(replay_score, "archived_high_precision_coverage_replay", z)]
        result = native_newton(site, solver, params, context, candidates)
        payload = {"site": site, "catmap_version": "0.3.1", "model_sha256": "86115cefb7784fe28ab519271cdfaf0e308b6f6e3a8308b1e90e4318cb497cfe",
                   "energy_sha256": "55ad4e94d9c603b6705ec8c92630c9a1ceff5859cf69865145152375e203ca30",
                   "search_bounds_logits_only": [-300, 300], "scipy_nfev_per_seed": SCIPY_NFEV,
                   "scipy_attempts": attempts, "catmap_newton_results": [r for r, _ in result]}
        dump_json(OUT / f"{site}_targeted_recovery.json", payload)
        return payload
    finally:
        os.chdir(oldcwd)


if __name__ == "__main__":
    all_results = {}
    for site in ("tiw", "tita"):
        try:
            all_results[site] = run(site)
            print(site, "SCIPY_CANDIDATES", len(all_results[site]["scipy_attempts"]),
                  "CATMAP_POLISH", len(all_results[site]["catmap_newton_results"]))
            print(json.dumps(all_results[site]["catmap_newton_results"][-1] if all_results[site]["catmap_newton_results"] else {}, ensure_ascii=False))
        except Exception as e:
            print(site, "ERROR", type(e).__name__, str(e))
            all_results[site] = {"error": f"{type(e).__name__}: {e}"}
    dump_json(OUT / "combined_targeted_recovery.json", all_results)
