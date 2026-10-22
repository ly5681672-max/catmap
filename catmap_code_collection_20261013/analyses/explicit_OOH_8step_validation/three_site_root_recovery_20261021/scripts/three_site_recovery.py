from __future__ import annotations

"""Finite, same-model CatMAP root-seed recovery for Ti-Mo/Ti-W/Ti-Ta.

This uses CatMAP 0.3.1's parsed rate constants, steady-state function and
NewtonRoot. It changes only numerical initial guesses; model energies,
reaction network, temperature, reservoirs and inventory constraints are fixed.
"""

import json
import os
from pathlib import Path
from collections import Counter

import mpmath as mp

ROOT = Path(r"H:\catmap\catmap\catmap_code_collection_20261013\analyses\explicit_OOH_8step_validation\multisite_conditional_replay_20261021")
TMP = Path(r"C:\Users\Administrator\Documents\codex探索\three_site_recovery_work")
TARGETS = ("timo", "tiw", "tita")
SEED_SURFACES = ("tico", "tihf", "tini", "tire", "timn", "tiv", "tife", "ticr", "titi", "ti", "tizr")
INDEPENDENT = ("C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "OH_s", "OOH_s", "O_s")
MAX_STEPS_PER_SEED = 35
MAX_SEEDS_PER_SITE = 6
ACCEPT_LOCAL_REL = mp.mpf("1e-12")


def json_dump(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def setup_model(site: str, work: Path):
    src = ROOT / "sites" / site / f"{site}_static_probe.mkm"
    source = src.read_text(encoding="utf-8")
    # Keep all source model definitions as-is; only isolate CatMAP cache/output.
    setup = work / f"{site}_recovery.mkm"
    setup.write_text(source + f"\n\ndata_file = {str(work / (site + '_recovery.pkl'))!r}\n", encoding="utf-8")
    from catmap import ReactionModel

    old = Path.cwd()
    os.chdir(work)
    try:
        model = ReactionModel(setup_file=setup.name)
        model.parse()
        solver, scaler = model.solver, model.scaler
        mp.mp.dps = 385
        model._math = mp
        model._mpfloat = mp.mpf
        model._matrix = mp.matrix
        model._Axb_solver = mp.lu_solve
        model._kB = mp.mpf("8.617332478e-5")
        model._h = mp.mpf("4.135667516e-15")
        mp.infnorm = lambda x: mp.norm(x, "inf")
        mp.lsqnorm = lambda x: mp.norm(x, 2) ** 2
        solver.compile()
        params = list(scaler.get_rxn_parameters([mp.mpf("333.15")]))
        return model, solver, params, old
    except Exception:
        os.chdir(old)
        raise


def run_site(site: str) -> dict:
    work = TMP / site
    work.mkdir(parents=True, exist_ok=True)
    model, solver, params, old_cwd = setup_model(site, work)
    try:
        names = list(solver.adsorbate_names)
        idx = {name: i for i, name in enumerate(names)}
        free_idx = [idx[n] for n in INDEPENDENT]
        derived = {"Ha_s": ("OOH_s", "OH_s", "O_s"), "Hb_s": ("O_s",)}
        weights = (1, 1, 1, 1, 2, 2, 3)
        rows = [idx[n] for n in INDEPENDENT]
        n_rxn = len(solver.elementary_rxns)
        stoich = [[0] * n_rxn for _ in names]
        for j, rxn in enumerate(solver.elementary_rxns):
            ci, cf = Counter(rxn[0]), Counter(rxn[-1])
            for sp in names:
                stoich[idx[sp]][j] = cf.get(sp, 0) - ci.get(sp, 0)

        def state(z):
            odds = [mp.exp(v) for v in z]
            den = 1 + sum(w * q for w, q in zip(weights, odds))
            cv = [mp.mpf(0)] * len(names)
            for i, q in zip(free_idx, odds):
                cv[i] = q / den
            for target, sources in derived.items():
                cv[idx[target]] = sum(cv[idx[s]] for s in sources)
            # Compute vacancy directly from the chart denominator instead of
            # subtracting a nearly unit site sum (catastrophic cancellation).
            vacancy = 1 / den
            return cv, vacancy

        def residual(z):
            cv, _ = state(z)
            solver.get_rate_constants(params, list(cv))
            return [solver.steady_state_function(list(cv))[i] for i in rows]

        # Seed from validated roots on other sites. These are numerical initial
        # guesses only; every candidate is solved against the target's own input.
        seeds = []
        seed_errors = []
        for source in SEED_SURFACES:
            p = ROOT / "sites" / source / f"constrained_{source}_result.json"
            if not p.exists():
                continue
            dat = json.loads(p.read_text(encoding="utf-8"))
            try:
                cv = [mp.mpf(dat["coverage"][n]) for n in names]
                vac = mp.mpf(dat["vacant"])
                if vac <= 0 or any(cv[i] <= 0 for i in free_idx):
                    continue
                logits = [mp.log(cv[i] / vac) for i in free_idx]
                # Root profile / dominant species provides a reproducible seed rank.
                dominant = max(dat["coverage"], key=lambda n: mp.mpf(dat["coverage"][n]))
                seeds.append((source, logits, dominant))
            except Exception as exc:
                seed_errors.append({"source": source, "error": f"{type(exc).__name__}: {exc}"})

        # Retain structurally distinct starting branches, then fill remaining
        # slots in fixed order. No seed is accepted without target residual tests.
        picked = []
        seen_dominants = set()
        for item in seeds:
            if item[2] not in seen_dominants:
                picked.append(item)
                seen_dominants.add(item[2])
            if len(picked) >= MAX_SEEDS_PER_SITE:
                break
        for item in seeds:
            if len(picked) >= MAX_SEEDS_PER_SITE:
                break
            if item not in picked:
                picked.append(item)
        for k in ("-8", "-4", "-2"):
            picked.append((f"uniform_logit_{k}", [mp.mpf(k)] * len(INDEPENDENT), "uniform"))

        def seed_metric(item):
            _, z, _ = item
            cv, _ = state(z)
            solver._coverage = list(cv)
            d = list(solver.get_directional_rates(params))
            scales = []
            for row in rows:
                gross = sum(abs(stoich[row][j]) * (abs(d[j]) + abs(d[n_rxn + j]))
                            for j in range(n_rxn))
                scales.append(mp.mpf(str(gross)) if gross > 0 else mp.mpf("1e-300"))
            f = residual(z)
            score = max(abs(f[k] / scales[k]) for k in range(len(rows)))
            return score

        picked.sort(key=seed_metric)

        attempts = []
        accepted = None
        from catmap.solvers.solver_base import NewtonRoot

        for seed_name, start, dominant in picked:
            cv0, vac0 = state(start)
            solver._coverage = list(cv0)
            # Native directional rates are used only to form fixed row scales,
            # never to subtract rates or calculate TOF.
            directions = list(solver.get_directional_rates(params))
            scales = []
            for row in rows:
                gross = sum(abs(stoich[row][j]) * (abs(directions[j]) + abs(directions[n_rxn + j]))
                            for j in range(n_rxn))
                scales.append(mp.mpf(str(gross)) if gross > 0 else mp.mpf("1e-300"))

            def scaled(z):
                f = residual(z)
                return tuple(f[k] / scales[k] for k in range(len(rows)))

            def jac(z):
                z = list(z)
                return mp.matrix([[mp.diff(lambda q: scaled([q if k == col else z[k]
                                                              for k in range(len(z))])[row], z[col])
                                   for col in range(len(z))] for row in range(len(rows))])

            iterator = NewtonRoot(scaled, list(start), mp.matrix, mp.mpf, mp.lu_solve,
                                  J=jac, norm=lambda x: mp.norm(x, "inf"), verbose=0)
            z = list(start)
            err = None
            norm = mp.inf
            steps = 0
            try:
                for steps, (x, norm) in enumerate(iterator, start=1):
                    z = list(x)
                    if norm < mp.mpf("1e-30") or steps >= MAX_STEPS_PER_SEED:
                        break
            except Exception as exc:
                err = f"{type(exc).__name__}: {exc}"
                norm = mp.norm(mp.matrix(scaled(z)), "inf")
            cv, vac = state(z)
            f = residual(z)
            absmax = max(abs(v) for v in f)
            # Local-flow residual is a gate only; get_directional_rates is
            # floating-point in this CatMAP version, so record its limitation.
            solver._coverage = list(cv)
            dirs = list(solver.get_directional_rates(params))
            local = []
            for row in rows:
                gross = sum(abs(stoich[row][j]) * (abs(dirs[j]) + abs(dirs[n_rxn + j]))
                            for j in range(n_rxn))
                local.append(float(abs(f[rows.index(row)]) / mp.mpf(str(gross))) if gross else None)
            inv1 = cv[idx["Ha_s"]] - cv[idx["OOH_s"]] - cv[idx["OH_s"]] - cv[idx["O_s"]]
            inv2 = cv[idx["Hb_s"]] - cv[idx["O_s"]]
            row = {"seed": seed_name, "seed_dominant_species": dominant, "steps": steps,
                   "scaled_residual_inf": mp.nstr(norm, 20), "max_abs_residual": mp.nstr(absmax, 20),
                   "max_local_normalized_residual_float_scale": max((x for x in local if x is not None), default=None),
                   "vacancy": mp.nstr(vac, 30), "I1": mp.nstr(inv1, 12), "I2": mp.nstr(inv2, 12),
                   "site_sum_error": mp.nstr(abs(sum(cv) + vac - 1), 15),
                   "root_logits": [mp.nstr(x, 30) for x in z],
                   "coverage": {name: mp.nstr(cv[i], 25) for i, name in enumerate(names)},
                   "solver_error": err}
            attempts.append(row)
            if norm < mp.mpf("1e-30") and max((x for x in local if x is not None), default=1) < float(ACCEPT_LOCAL_REL):
                # Re-evaluate all outputs with CatMAP high-precision net-rate API.
                ks = solver.get_rate_constants(params, list(cv))
                netrates = list(solver.get_rxn_rates(list(cv), ks))
                tofs = list(solver.get_turnover_frequency(params, rates=netrates, verify_coverages=False))
                prod_i = solver.gas_names.index("C6H12O_g")
                tof = tofs[prod_i]
                if vac > 0 and all(x >= 0 for x in cv) and tof > 0:
                    accepted = {"root": z, "coverage": cv, "vacancy": vac, "tof": tof,
                                "netrates": netrates, "local": local, "seed": seed_name,
                                "steps": steps, "norm": norm}
                    break

        # CatMAP's own NewtonRoot is local and may stop at a stationary point.
        # As documented in its implementation, use a mature low-precision
        # least-squares call only to locate a seed, then require CatMAP
        # NewtonRoot to polish and verify that seed with native equations.
        if accepted is None:
            import numpy as np
            from scipy.optimize import least_squares
            seed_name, start, dominant = picked[0]
            cv0, _ = state(start)
            solver._coverage = list(cv0)
            d0 = list(solver.get_directional_rates(params))
            scales0 = []
            for row in rows:
                gross = sum(abs(stoich[row][j]) * (abs(d0[j]) + abs(d0[n_rxn + j]))
                            for j in range(n_rxn))
                scales0.append(mp.mpf(str(gross)) if gross > 0 else mp.mpf("1e-300"))

            def scipy_resid(zf):
                zmp = [mp.mpf(float(x)) for x in zf]
                f = residual(zmp)
                vals = tuple(f[k] / scales0[k] for k in range(len(rows)))
                return np.asarray([float(v) for v in vals], dtype=float)

            scipy_attempt = {"seed": seed_name, "method": "SciPy least_squares seed only; CatMAP NewtonRoot required for acceptance"}
            try:
                ls = least_squares(scipy_resid, np.asarray([float(x) for x in start]),
                                   bounds=(-500.0, 500.0), max_nfev=220,
                                   xtol=1e-12, ftol=1e-12, gtol=1e-12)
                z0 = [mp.mpf(float(x)) for x in ls.x]
                scipy_attempt.update({"success": bool(ls.success), "status": int(ls.status),
                                      "message": str(ls.message), "nfev": int(ls.nfev),
                                      "cost": float(ls.cost), "optimality": float(ls.optimality),
                                      "max_scaled_residual_seed": float(np.max(np.abs(scipy_resid(ls.x))))})
                # Recompute local scales at the candidate and use the official
                # CatMAP NewtonRoot at current high precision for final polish.
                cv0, _ = state(z0)
                solver._coverage = list(cv0)
                d0 = list(solver.get_directional_rates(params))
                scales0 = []
                for row in rows:
                    gross = sum(abs(stoich[row][j]) * (abs(d0[j]) + abs(d0[n_rxn + j]))
                                for j in range(n_rxn))
                    scales0.append(mp.mpf(str(gross)) if gross > 0 else mp.mpf("1e-300"))

                def scipy_scaled(z):
                    f = residual(z)
                    return tuple(f[k] / scales0[k] for k in range(len(rows)))

                def scipy_jac(z):
                    z = list(z)
                    return mp.matrix([[mp.diff(lambda q: scipy_scaled([q if k == col else z[k]
                                                                         for k in range(len(z))])[row], z[col])
                                       for col in range(len(z))] for row in range(len(rows))])

                polish = NewtonRoot(scipy_scaled, z0, mp.matrix, mp.mpf, mp.lu_solve,
                                    J=scipy_jac, norm=lambda x: mp.norm(x, "inf"), verbose=0)
                zp = list(z0)
                pnorm = mp.inf
                psteps = 0
                perror = None
                try:
                    for psteps, (x, pnorm) in enumerate(polish, start=1):
                        zp = list(x)
                        if pnorm < mp.mpf("1e-30") or psteps >= MAX_STEPS_PER_SEED:
                            break
                except Exception as exc:
                    perror = f"{type(exc).__name__}: {exc}"
                    pnorm = mp.norm(mp.matrix(scipy_scaled(zp)), "inf")
                cvp, vp = state(zp)
                fp = residual(zp)
                solver._coverage = list(cvp)
                dp = list(solver.get_directional_rates(params))
                lp = []
                for row in rows:
                    gross = sum(abs(stoich[row][j]) * (abs(dp[j]) + abs(dp[n_rxn + j]))
                                for j in range(n_rxn))
                    lp.append(float(abs(fp[rows.index(row)] / mp.mpf(str(gross)))) if gross else None)
                scipy_attempt.update({"polish_steps": psteps, "polish_error": perror,
                                      "polished_scaled_residual_inf": mp.nstr(pnorm, 30),
                                      "polished_max_abs_residual": mp.nstr(max(abs(v) for v in fp), 30),
                                      "polished_max_local_normalized_residual_float_scale": max(x for x in lp if x is not None),
                                      "polished_vacancy": mp.nstr(vp, 30),
                                      "polished_root_logits": [mp.nstr(x, 60) for x in zp]})
                attempts.append(scipy_attempt)
                if pnorm < mp.mpf("1e-30") and max(x for x in lp if x is not None) < float(ACCEPT_LOCAL_REL) and vp > 0:
                    ksp = solver.get_rate_constants(params, list(cvp))
                    nrp = list(solver.get_rxn_rates(list(cvp), ksp))
                    tofp = list(solver.get_turnover_frequency(params, rates=nrp, verify_coverages=False))
                    product_i = solver.gas_names.index("C6H12O_g")
                    if tofp[product_i] > 0 and all(v >= 0 for v in cvp):
                        accepted = {"root": zp, "coverage": cvp, "vacancy": vp, "tof": tofp[product_i],
                                    "netrates": nrp, "local": lp, "seed": "scipy_seed_then_CatMAP_NewtonRoot",
                                    "steps": psteps, "norm": pnorm}
            except Exception as exc:
                scipy_attempt.update({"success": False, "error": f"{type(exc).__name__}: {exc}"})
                attempts.append(scipy_attempt)

        out = {"surface": site, "catmap_version": "0.3.1", "method": "CatMAP native steady-state equations + NewtonRoot; target-specific energies; exact zero-inventory chart",
               "dps": mp.mp.dps, "attempt_limit_per_seed": MAX_STEPS_PER_SEED,
               "seed_limit": MAX_SEEDS_PER_SITE, "attempts": attempts,
               "seed_parse_errors": seed_errors,
               "accepted": accepted is not None,
               "numeric_status": "NUMERIC_PARTIAL" if accepted else "NUMERIC_FAIL"}
        if accepted:
            out["accepted_seed"] = accepted["seed"]
            out["tof_C6H12O_s-1"] = mp.nstr(accepted["tof"], 80)
            out["log10_tof"] = mp.nstr(mp.log10(accepted["tof"]), 60)
            out["coverage"] = {name: mp.nstr(accepted["coverage"][i], 80) for i, name in enumerate(names)}
            out["vacancy"] = mp.nstr(accepted["vacancy"], 80)
            out["net_reaction_rates"] = [mp.nstr(v, 50) for v in accepted["netrates"]]
            out["max_local_normalized_residual_float_scale"] = max(x for x in accepted["local"] if x is not None)
            out["root_logits"] = [mp.nstr(x, 100) for x in accepted["root"]]
            out["important_note"] = "Candidate only; independent cross-precision re-solve and official flux audit still required before NUMERIC_PASS. Directional API float scales were used only for local normalization."
            if site == "timo":
                # Exact arbitrary-precision directional mass-action diagnostic
                # using CatMAP's parsed elementary reactions and rate constants.
                def directional_exact(cv, ks):
                    # Mirror CatMAP's actual Python-float reservoir inputs
                    # exactly; decimal-string conversion would subtly alter
                    # the inputs and be magnified by large forward/reverse
                    # rates when reconstructing a small net flux.
                    pressures = [mp.mpf(v) for v in solver.gas_pressures]
                    empty = [mp.mpf(str(total)) - sum(cv[i] * mp.mpf(str(solver.species_definitions[n]["n_sites"]))
                                                       for i,n in enumerate(solver.adsorbate_names)
                                                       if solver.species_definitions[n]["site"] == site)
                             for site,total in ((s,solver.species_definitions[s]["total"]) for s in solver.site_names)]
                    env = {"kf": ks[:n_rxn], "kr": ks[n_rxn:], "theta": cv,
                           "p": pressures, "s": empty}
                    vals_f, vals_r = [], []
                    for j, rxn in enumerate(solver.elementary_rxns):
                        # Build exactly the same mass-action term strings as
                        # CatMAP's mean-field solver, then evaluate with mpf.
                        fexpr = solver.rate_equation_term(rxn[0], f"kf[{j}]")
                        rexpr = solver.rate_equation_term(rxn[-1], f"kr[{j}]")
                        vals_f.append(eval(fexpr, {"__builtins__": {}}, env))
                        vals_r.append(eval(rexpr, {"__builtins__": {}}, env))
                    return vals_f, vals_r

                root = list(accepted["root"])
                cvb, vb = state(root)
                ksb = solver.get_rate_constants(params, list(cvb))
                fwd, rev = directional_exact(cvb, ksb)
                net_exact = [a - b for a, b in zip(fwd, rev)]
                net_official = list(solver.get_rxn_rates(list(cvb), ksb))
                residual_b = list(solver.steady_state_function(list(cvb)))
                exact_stoich = [[0] * n_rxn for _ in names]
                for j, rxn in enumerate(solver.elementary_rxns):
                    ci, cf = Counter(rxn[0]), Counter(rxn[-1])
                    for sp in names:
                        exact_stoich[idx[sp]][j] = cf.get(sp, 0) - ci.get(sp, 0)
                local_flux = {}
                for n in INDEPENDENT:
                    i = idx[n]
                    gross = sum(abs(exact_stoich[i][j]) * (fwd[j] + rev[j]) for j in range(n_rxn))
                    local_flux[n] = {"abs_dtheta_dt": mp.nstr(abs(residual_b[i]), 35),
                                     "gross_forward_reverse_scale": mp.nstr(gross, 35),
                                     "normalized_residual": mp.nstr(abs(residual_b[i]) / gross, 35) if gross else "UNDEFINED_ZERO_GROSS_FLOW"}
                product_i = solver.gas_names.index("C6H12O_g")
                tof_b = list(solver.get_turnover_frequency(params, rates=net_official, verify_coverages=False))[product_i]
                local_net_tof = net_official[4]
                site_error_b = abs(sum(cvb) + vb - 1)
                inv1_b = cvb[idx["Ha_s"]] - cvb[idx["OOH_s"]] - cvb[idx["OH_s"]] - cvb[idx["O_s"]]
                inv2_b = cvb[idx["Hb_s"]] - cvb[idx["O_s"]]
                out["high_precision_flux_audit"] = {
                    "directional_rates_method": "Arbitrary-precision mass-action evaluation of CatMAP parsed elementary_rxns, CatMAP rate constants and gas reservoir values; not get_directional_rates (which casts to float).",
                    "directional_forward": [mp.nstr(x, 40) for x in fwd],
                    "directional_reverse": [mp.nstr(x, 40) for x in rev],
                    "net_rates_exact_forward_minus_reverse": [mp.nstr(x, 40) for x in net_exact],
                    "max_abs_difference_vs_CatMAP_get_rxn_rates": mp.nstr(max(abs(a-b) for a,b in zip(net_exact,net_official)), 25),
                    "local_species_flux": local_flux,
                    "max_local_normalized_residual": mp.nstr(max(mp.mpf(v["normalized_residual"]) for v in local_flux.values()), 30),
                    "site_sum_error": mp.nstr(site_error_b, 25), "I1": mp.nstr(inv1_b, 25), "I2": mp.nstr(inv2_b, 25),
                    "official_product_TOF": mp.nstr(tof_b, 40),
                    "product_desorption_net_rate": mp.nstr(local_net_tof, 40),
                    "TOF_minus_product_step_rate": mp.nstr(tof_b-local_net_tof, 25),
                    "net_rate_cycle_spread_relative": mp.nstr(max(abs(x-net_exact[0]) for x in net_exact)/max(abs(x) for x in net_exact), 25),
                }

                # Independent nearby-start re-solve at 250 decimal digits.
                mp.mp.dps = 250
                for memo_name in ("_rate_constant_memoize", "_steady_state_memoize"):
                    memo = getattr(solver, memo_name, None)
                    if hasattr(memo, "clear"):
                        memo.clear()
                params250 = [mp.mpf(str(x)) for x in params]
                perturb = [mp.mpf("0.01") if j % 2 == 0 else mp.mpf("-0.01") for j in range(len(root))]
                zstart = [x+d for x,d in zip(root, perturb)]
                cvstart, _ = state(zstart)
                ksstart = solver.get_rate_constants(params250, list(cvstart))
                fstart, rstart = directional_exact(cvstart, ksstart)
                scales = []
                for n in INDEPENDENT:
                    i = idx[n]
                    scales.append(sum(abs(exact_stoich[i][j])*(fstart[j]+rstart[j]) for j in range(n_rxn)))
                def cross_resid(z):
                    cv, _ = state(z)
                    solver.get_rate_constants(params250, list(cv))
                    ff = list(solver.steady_state_function(list(cv)))
                    return tuple(ff[idx[n]]/scales[k] for k,n in enumerate(INDEPENDENT))
                def cross_jac(z):
                    z=list(z)
                    return mp.matrix([[mp.diff(lambda q: cross_resid([q if kk==col else z[kk] for kk in range(len(z))])[rr],z[col])
                                       for col in range(len(z))] for rr in range(len(rows))])
                from catmap.solvers.solver_base import NewtonRoot as _NewtonRoot
                cross_iter = _NewtonRoot(cross_resid,zstart,mp.matrix,mp.mpf,mp.lu_solve,
                                         J=cross_jac,norm=lambda x:mp.norm(x,"inf"),verbose=0)
                zc=list(zstart); cnorm=mp.inf; csteps=0; cerr=None
                try:
                    for csteps,(xx,nn) in enumerate(cross_iter,start=1):
                        zc=list(xx); cnorm=nn
                        if cnorm < mp.mpf("1e-35") or csteps >= 80:
                            break
                except Exception as exc:
                    cerr=f"{type(exc).__name__}: {exc}"
                    cnorm=mp.norm(mp.matrix(cross_resid(zc)),"inf")
                cvc, vc=state(zc)
                ksc=solver.get_rate_constants(params250,list(cvc))
                nc=list(solver.get_rxn_rates(list(cvc),ksc))
                tof_c=list(solver.get_turnover_frequency(params250,rates=nc,verify_coverages=False))[product_i]
                dlog=abs(mp.log10(tof_c)-mp.log10(tof_b)) if tof_c>0 and tof_b>0 else mp.inf
                cov_rel={n: mp.nstr(abs(cvc[idx[n]]-cvb[idx[n]])/max(abs(cvc[idx[n]]),abs(cvb[idx[n]])),25)
                                    if max(abs(cvc[idx[n]]),abs(cvb[idx[n]]))>mp.mpf("1e-250") else "UNRESOLVED_BELOW_PRECISION"
                         for n in ("C6H12_s","H2O2_s","H2O_s","OOH_s","C6H12O_s")}
                out["cross_precision_replay"] = {
                    "method": "CatMAP native steady_state_function + NewtonRoot; same exact zero-inventory chart; alternating 0.01 logit perturbation",
                    "dps": 250, "iterations": csteps, "solver_error": cerr,
                    "scaled_residual_inf": mp.nstr(cnorm,30),
                    "local_coverage_relative_differences": cov_rel,
                    "baseline_TOF": mp.nstr(tof_b,40), "replay_TOF": mp.nstr(tof_c,40),
                    "delta_log10_TOF": mp.nstr(dlog,30),
                    "relative_TOF_difference_ppm": mp.nstr(abs(tof_c-tof_b)/max(abs(tof_c),abs(tof_b))*mp.mpf("1e6"),25),
                    "max_logit_difference_same_branch": mp.nstr(max(abs(a-b) for a,b in zip(root,zc)),25),
                    "replay_vacancy": mp.nstr(vc,35), "replay_site_sum_error": mp.nstr(abs(sum(cvc)+vc-1),20),
                    "replay_I1": mp.nstr(cvc[idx["Ha_s"]]-cvc[idx["OOH_s"]]-cvc[idx["OH_s"]]-cvc[idx["O_s"]],20),
                    "replay_I2": mp.nstr(cvc[idx["Hb_s"]]-cvc[idx["O_s"]],20),
                }
                numeric_pass = (mp.mpf(out["high_precision_flux_audit"]["max_local_normalized_residual"]) < mp.mpf("1e-12")
                                and dlog <= mp.mpf("2.17147e-6")
                                and max(mp.mpf(v) for v in cov_rel.values() if v != "UNRESOLVED_BELOW_PRECISION") <= mp.mpf("1e-4")
                                and vc > 0 and tof_c > 0)
                out["numeric_status"] = "NUMERIC_PASS" if numeric_pass else "NUMERIC_PARTIAL"
                out["cross_precision_replay"]["gate_pass"] = bool(numeric_pass)
                mp.mp.dps = 385
        json_dump(work / f"{site}_recovery.json", out)
        return out
    finally:
        os.chdir(old_cwd)


def main() -> None:
    TMP.mkdir(parents=True, exist_ok=True)
    all_results = []
    for site in TARGETS:
        try:
            res = run_site(site)
        except Exception as exc:
            res = {"surface": site, "accepted": False, "numeric_status": "NUMERIC_FAIL",
                   "fatal_error": f"{type(exc).__name__}: {exc}"}
            json_dump(TMP / site / f"{site}_recovery.json", res)
        all_results.append(res)
        print(json.dumps({"surface": site, "accepted": res.get("accepted"),
                          "status": res.get("numeric_status"), "error": res.get("fatal_error"),
                          "attempts": len(res.get("attempts", []))}, ensure_ascii=False))
    json_dump(TMP / "THREE_SITE_RECOVERY_ALL.json", all_results)


if __name__ == "__main__":
    main()
