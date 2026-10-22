from __future__ import annotations

import os
import csv
import json
from collections import Counter
from pathlib import Path

SOURCE = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\alkene_epoxidation.mkm")
ENERGIES = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\energies_dft_neb.txt")
SURFACE = os.environ.get("CATMAP_8STEP_SURFACE", "tico").strip().lower()
WORK = Path(os.environ.get(
    "CATMAP_8STEP_WORK",
    str(Path(r"C:\Users\Administrator\Documents\codex探索\catmap_original8_zero_inventory_multisite") / SURFACE),
))
ROOTS = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\analysis_ooh\steady_state_branch_audit\combined_validated_roots.csv")


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    setup = WORK / f"{SURFACE}_static_probe.mkm"
    source = SOURCE.read_text(encoding="utf-8")
    overrides = f"""
# Isolated replay of the archived zero-inventory root with the approved
# original static reservoir configuration. Scientific inputs are not altered.
surface_names = [{SURFACE!r}]
scaler = 'ThermodynamicScaler'
descriptor_names = ['temperature']
descriptor_ranges = None
resolution = 1
output_variables = []
pressure_mode = 'static'
gas_thermo_mode = 'frozen_gas'
adsorbate_thermo_mode = 'frozen_adsorbate'
input_file = {str(ENERGIES)!r}
data_file = {f'{SURFACE}_static_probe.pkl'!r}
decimal_precision = 360
from mpmath import mp
tolerance = mp.mpf('1e-260')
max_rootfinding_iterations = 250
use_numbers_solver = False
"""
    setup.write_text(source + "\n\n" + overrides, encoding="utf-8")

    from catmap import ReactionModel

    previous = Path.cwd()
    os.chdir(WORK)
    try:
        model = ReactionModel(setup_file=setup.name)
        model.run()
        solver = model.solver
        scaler = model.scaler
        print("model surfaces:", model.surface_names)
        print("adsorbates:", solver.adsorbate_names)
        print("TS:", solver.transition_state_names)
        print("gas:", solver.gas_names)
        print("solver type:", type(solver).__name__)
        print("scaler type:", type(scaler).__name__)
        solver.compile()
        descriptors = [model._mpfloat("333.15")]
        params = scaler.get_rxn_parameters(descriptors)
        print("parameter count:", len(params))
        print("gas pressures:", scaler.thermodynamics.gas_pressures)
        print("solver pressure:", solver.gas_pressures)
        print("free energies loaded:", len(scaler._electronic_energy_dict))
        print("compiled:", solver._compiled)
        print("steady-state function:", solver.steady_state_function.__name__)

        # A seven-coordinate chart of the zero-inventory class. Ha and Hb
        # are reconstructed from the two exact invariants; the vacant site
        # is reconstructed from CatMAP's site balance.
        import mpmath as mp
        mp.mp.dps = 385
        names = list(solver.adsorbate_names)
        indep = ["C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "OH_s", "OOH_s", "O_s"]
        idx = {name: names.index(name) for name in names}
        free_idx = [idx[name] for name in indep]
        derived = {"Ha_s": ("OOH_s", "OH_s", "O_s"), "Hb_s": ("O_s",)}
        desc = [mp.mpf("333.15")]
        rxn_params = list(scaler.get_rxn_parameters(desc))

        # The softmax-like site-fraction coordinates ensure positive species
        # and vacancy coverages while retaining the exact site balance.
        site_weights = [1, 1, 1, 1, 2, 2, 3]
        def state_from_logits(z):
            odds = [mp.exp(v) for v in z]
            den = 1 + sum(w * q for w, q in zip(site_weights, odds))
            cv = [mp.mpf("0")] * len(names)
            for i, q in zip(free_idx, odds):
                cv[i] = q / den
            for target, sources in derived.items():
                cv[idx[target]] = sum(cv[idx[s]] for s in sources)
            vacancy = 1 - sum(cv)
            return cv, vacancy

        def evaluate(z, want_rates=False):
            cv, vacancy = state_from_logits(z)
            ks = solver.get_rate_constants(rxn_params, list(cv))
            net = list(solver.steady_state_function(list(cv)))
            if want_rates:
                rates = list(solver.get_rxn_rates(list(cv), ks))
                tof = list(solver.get_turnover_frequency(rxn_params, rates=rates, verify_coverages=False))
                return cv, vacancy, net, rates, tof
            return cv, vacancy, net

        def csv_root_seed():
            candidates = []
            with ROOTS.open(encoding="utf-8-sig", newline="") as fh:
                for row in csv.DictReader(fh):
                    if row.get("surface") == SURFACE and row.get("quality_pass", "").lower() == "true":
                        candidates.append(row)
            # Prefer the audited Hb/O branch, but project it onto I1=I2=0;
            # this is a numerical initial guess only, not a selected root.
            row = next((r for r in candidates if r.get("dominant_species") == "Hb_s"), candidates[0])
            raw = [mp.mpf(x) for x in row["full_coverage_vector_high_precision"].split(";")]
            current = {name: raw[i] for i, name in enumerate(names)}
            free = [current[name] for name in indep]
            ha = free[5] + free[4] + free[6]
            hb = free[6]
            vac = 1 - sum(free) - ha - hb
            if vac <= 0 or any(x <= 0 for x in free):
                raise ValueError("projected historical root is not in the positive zero-inventory domain")
            logits = [mp.log(x / vac) for x in free]
            seedinfo = {"source_root_id": row.get("root_id"), "source_dominant_species": row.get("dominant_species"),
                        "projected_I1": str((ha-free[5]-free[4]-free[6])), "projected_I2": str(hb-free[6]),
                        "projected_vacancy": mp.nstr(vac, 30)}
            return logits, seedinfo

        try:
            z0, seedinfo = csv_root_seed()
            cv0, vac0, net0 = evaluate(z0)
        except (ValueError, IndexError, KeyError) as exc:
            # A historical root is only a numerical seed. If its projection
            # lies outside the positive zero-inventory chart, do not treat
            # that as evidence that the physical model has no root: retain
            # the exact model/inventory and try the bounded uniform-logit
            # initial guesses below.
            z0 = None
            seedinfo = {
                "source": "uniform_logit_fallback_only",
                "historical_seed_rejected": str(exc),
                "scientific_inputs_changed": False,
            }
        independent_rows = [idx[n] for n in indep]
        # Per-species gross-flow scales are derived from CatMAP's own
        # directional-rate expressions, solely to condition the residual.
        n_rxn = len(solver.elementary_rxns)
        stoich = [[0] * n_rxn for _ in names]
        for j, rxn in enumerate(solver.elementary_rxns):
            initial, final = rxn[0], rxn[-1]
            ci, cf = Counter(initial), Counter(final)
            for sp in names:
                stoich[idx[sp]][j] = cf.get(sp, 0) - ci.get(sp, 0)
        from catmap.solvers.solver_base import NewtonRoot
        seed_candidates = []
        if z0 is not None:
            seed_candidates.append(("projected_hb_branch", list(z0)))
        seed_candidates.extend((f"uniform_logit_{v}", [mp.mpf(v)] * len(indep)) for v in ("-20", "-10", "-3", "-1"))
        attempts = []
        accepted = None
        for seed_name, start in seed_candidates:
            start_cv, start_vac, start_net = evaluate(start)
            solver._coverage = list(start_cv)
            directional = list(solver.get_directional_rates(rxn_params))
            row_scales = []
            for row in independent_rows:
                gross = sum(abs(stoich[row][j]) * (abs(directional[j]) + abs(directional[n_rxn + j]))
                            for j in range(n_rxn))
                row_scales.append(mp.mpf(str(gross)) if gross > 0 else mp.mpf("1e-300"))
            global_scale = max(row_scales)
            row_scales = [max(s, global_scale * mp.mpf("1e-250")) for s in row_scales]

            def residual(z):
                z = list(z)
                _, _, net = evaluate(z)
                return tuple(net[i] / row_scales[j] for j, i in enumerate(independent_rows))

            def jacobian(z):
                z = list(z)
                return mp.matrix([[mp.diff(
                    lambda q: residual([q if kk == col else z[kk] for kk in range(len(indep))])[row],
                    z[col],
                ) for col in range(len(indep))] for row in range(len(indep))])

            root_iterator = NewtonRoot(
                residual, start, mp.matrix, mp.mpf, mp.lu_solve,
                J=jacobian, norm=lambda x: mp.norm(x, "inf"), verbose=0,
            )
            root_try = list(start)
            scaled_try = mp.inf
            count = 0
            root_error = None
            try:
                for count, (root_matrix, scaled_try) in enumerate(root_iterator, start=1):
                    root_try = list(root_matrix)
                    if scaled_try < mp.mpf("1e-50") or count >= 250:
                        break
            except Exception as exc:
                root_error = f"{type(exc).__name__}: {exc}"
                scaled_try = mp.norm(mp.matrix(residual(root_try)), "inf")
            cv_try, vac_try, net_try, rates_try, tof_try = evaluate(root_try, want_rates=True)
            abs_resid = max(abs(net_try[i]) for i in independent_rows)
            rel_resid = abs_resid / max(abs(x) for x in rates_try)
            inv1_try = cv_try[idx["Ha_s"]] - cv_try[idx["OOH_s"]] - cv_try[idx["OH_s"]] - cv_try[idx["O_s"]]
            inv2_try = cv_try[idx["Hb_s"]] - cv_try[idx["O_s"]]
            ok = (scaled_try < mp.mpf("1e-50") and rel_resid < mp.mpf("1e-4")
                  and abs_resid < mp.mpf("1e-45") and abs(inv1_try) < mp.mpf("1e-300")
                  and abs(inv2_try) < mp.mpf("1e-300") and vac_try >= 0
                  and all(value >= 0 for value in cv_try) and sum(cv_try) + vac_try == 1)
            attempts.append({"seed": seed_name, "iterations": count, "scaled_residual_inf": mp.nstr(scaled_try, 30),
                             "max_abs_surface_residual": mp.nstr(abs_resid, 30), "relative_residual": mp.nstr(rel_resid, 30),
                             "I1": mp.nstr(inv1_try, 20), "I2": mp.nstr(inv2_try, 20), "vacancy": mp.nstr(vac_try, 20),
                             "accepted": bool(ok), "rootfinder_error": root_error})
            print("attempt:", json.dumps(attempts[-1]))
            if ok:
                accepted = (root_try, cv_try, vac_try, net_try, rates_try, tof_try, scaled_try, count, seed_name)
                break
            if scaled_try < mp.mpf("1e-40") and not seed_name.startswith("high_precision_polish"):
                # One precision-escalation polish is allowed for the closest
                # numerical candidate; this does not alter the model inputs.
                mp.mp.dps = 550
                seed_candidates.append((f"high_precision_polish_{seed_name}", root_try))

        (WORK / "constrained_solver_attempts.json").write_text(json.dumps(attempts, indent=2), encoding="utf-8")
        if accepted is None:
            print("NO ACCEPTED ZERO-INVENTORY ROOT; attempts retained. No TOF accepted.")
            return
        root, cv, vacancy, net, reaction_rates, tof, scaled_resid, iteration_count, accepted_seed = accepted
        accepted_root_high = list(root)
        # Independent nearby numerical initial guess at the lower working
        # precision; this is not a physical transient or attainability test.
        mp.mp.dps = 385
        for memo_name in ("_rate_constant_memoize", "_steady_state_memoize"):
            memo = getattr(solver, memo_name, None)
            if hasattr(memo, "clear"):
                memo.clear()
        perturb = [mp.mpf("0.01") if j % 2 == 0 else mp.mpf("-0.01") for j in range(len(indep))]
        cross_seed = [v + d for v, d in zip(accepted_root_high, perturb)]
        cross_cv, _, _ = evaluate(cross_seed)
        solver._coverage = list(cross_cv)
        cross_directional = list(solver.get_directional_rates(rxn_params))
        cross_scales = []
        for row in independent_rows:
            gross = sum(abs(stoich[row][j]) * (abs(cross_directional[j]) + abs(cross_directional[n_rxn + j]))
                        for j in range(n_rxn))
            cross_scales.append(mp.mpf(str(gross)) if gross > 0 else mp.mpf("1e-300"))
        cross_global = max(cross_scales)
        row_scales = [max(s, cross_global * mp.mpf("1e-250")) for s in cross_scales]
        cross_iterator = NewtonRoot(
            residual, cross_seed, mp.matrix, mp.mpf, mp.lu_solve,
            J=jacobian, norm=lambda x: mp.norm(x, "inf"), verbose=0,
        )
        cross_root = list(cross_seed)
        cross_scaled_resid = mp.inf
        cross_iterations = 0
        cross_error = None
        try:
            for cross_iterations, (x, norm_value) in enumerate(cross_iterator, start=1):
                cross_root, cross_scaled_resid = list(x), norm_value
                if norm_value < mp.mpf("1e-50") or cross_iterations >= 250:
                    break
        except Exception as exc:
            cross_error = f"{type(exc).__name__}: {exc}"
            cross_scaled_resid = mp.norm(mp.matrix(residual(cross_root)), "inf")
        cv2, vacancy2, net2, rates2, tof2 = evaluate(cross_root, want_rates=True)
        cross_abs = max(abs(net2[i]) for i in independent_rows)
        cross_rel = cross_abs / max(abs(x) for x in rates2)
        cross_i1 = cv2[idx["Ha_s"]] - cv2[idx["OOH_s"]] - cv2[idx["OH_s"]] - cv2[idx["O_s"]]
        cross_i2 = cv2[idx["Hb_s"]] - cv2[idx["O_s"]]
        cross_ok = (cross_scaled_resid < mp.mpf("1e-50") and cross_rel < mp.mpf("1e-4")
                    and cross_abs < mp.mpf("1e-45") and abs(cross_i1) < mp.mpf("1e-300")
                    and abs(cross_i2) < mp.mpf("1e-300") and vacancy2 >= 0)
        cross_log_delta = abs(mp.log10(tof2[0]) - mp.log10(tof[0])) if tof2[0] > 0 and tof[0] > 0 else mp.inf
        cross_cov_delta = max(abs(a - b) for a, b in zip(cv2, cv))
        cross_precision = {
            "seed": "accepted_root_with_alternating_0.01_logit_perturbation",
            "working_dps": 385, "iterations": cross_iterations, "rootfinder_error": cross_error,
            "scaled_residual_inf": mp.nstr(cross_scaled_resid, 30),
            "max_abs_surface_residual_independent": mp.nstr(cross_abs, 30),
            "relative_residual_to_max_net_rate": mp.nstr(cross_rel, 30),
            "I1": mp.nstr(cross_i1, 30), "I2": mp.nstr(cross_i2, 30),
            "site_sum": mp.nstr(sum(cv2) + vacancy2, 30),
            "log10_TOF": mp.nstr(mp.log10(tof2[0]), 40) if tof2[0] > 0 else "nonpositive",
            "delta_log10_TOF_vs_high_precision": mp.nstr(cross_log_delta, 30),
            "max_absolute_coverage_difference": mp.nstr(cross_cov_delta, 30),
            "passed_existing_steady_quality_gate": bool(cross_ok),
        }
        (WORK / "cross_precision_replay.json").write_text(json.dumps(cross_precision, indent=2), encoding="utf-8")
        print("cross precision replay:", json.dumps(cross_precision))
        mp.mp.dps = 550
        for memo_name in ("_rate_constant_memoize", "_steady_state_memoize"):
            memo = getattr(solver, memo_name, None)
            if hasattr(memo, "clear"):
                memo.clear()
        cv, vacancy, net, reaction_rates, tof = evaluate(root, want_rates=True)
        i1 = cv[idx["Ha_s"]] - cv[idx["OOH_s"]] - cv[idx["OH_s"]] - cv[idx["O_s"]]
        i2 = cv[idx["Hb_s"]] - cv[idx["O_s"]]
        relative_resid = max(abs(net[i]) for i in independent_rows) / max(abs(x) for x in reaction_rates)
        solver_success = True

        # Independent parsed-stoichiometry reconstruction of surface and gas
        # fluxes from CatMAP's eight official net reaction rates.
        gas_stoich = [[0] * n_rxn for _ in solver.gas_names]
        for j, rxn in enumerate(solver.elementary_rxns):
            ci, cf = Counter(rxn[0]), Counter(rxn[-1])
            for gi, gas_name in enumerate(solver.gas_names):
                gas_stoich[gi][j] = cf.get(gas_name, 0) - ci.get(gas_name, 0)
        surface_from_rates = [sum(stoich[i][j] * reaction_rates[j] for j in range(n_rxn)) for i in range(len(names))]
        gas_from_rates = [sum(gas_stoich[i][j] * reaction_rates[j] for j in range(n_rxn)) for i in range(len(solver.gas_names))]
        max_rate = max(abs(x) for x in reaction_rates)
        cycle_abs = max(abs(x - reaction_rates[3]) for x in reaction_rates)
        cycle_rel = cycle_abs / max_rate
        surface_balance_abs = max(abs(x) for x in surface_from_rates)
        gas_balance_abs = max(abs(a - b) for a, b in zip(gas_from_rates, tof))
        stoich_flux_rel = max(surface_balance_abs, gas_balance_abs) / max_rate

        # Local stability restricted to the same zero-inventory manifold.
        # Use a numerical Jacobian of CatMAP's generated steady-state function
        # in transformed coordinates, then convert zdot = D^{-1} F to its
        # local Jacobian at the fixed point.
        def selected_rates_at(z):
            state, _ = state_from_logits(z)
            solver.get_rate_constants(rxn_params, list(state))
            return list(solver.steady_state_function(list(state)))
        jlog = mp.matrix([[mp.diff(
            lambda q: selected_rates_at([q if kk == col else root[kk] for kk in range(len(indep))])[row], root[col]
        ) for col in range(len(indep))] for row in independent_rows])
        odds_root = [mp.exp(v) for v in root]
        den_root = 1 + sum(w * q for w, q in zip(site_weights, odds_root))
        p_root = [w * q / den_root for w, q in zip(site_weights, odds_root)]
        dmat = mp.matrix([[cv[idx[indep[i]]] * ((1 if i == j else 0) - p_root[j])
                           for j in range(len(indep))] for i in range(len(indep))])
        local_jacobian = mp.inverse(dmat) * jlog
        stability_error = None
        eigenvalues = []
        try:
            eigs = mp.eig(local_jacobian, left=False, right=False)
            eigenvalues = [mp.nstr(e, 30) for e in eigs]
            stable_local = all(mp.re(e) < 0 for e in eigs)
        except Exception as exc:
            stability_error = f"{type(exc).__name__}: {exc}"
            stable_local = None
        result = {
            "surface": SURFACE, "catmap_version": "0.3.1", "temperature_K": "333.15",
            "pressure_mode": "static",
            "effective_gas_reservoirs": {"C6H12_g": "1.0", "H2O2_g": "0.185", "H2O_g": "0.815", "C6H12O_g": "1e-20"},
            "gas_thermo_mode": "frozen_gas", "adsorbate_thermo_mode": "frozen_adsorbate",
            "decimal_precision_configured": 360, "tolerance_configured": "1e-260",
            "method": "CatMAP official rates and steady-state function; seven independent coverage coordinates; exact zero-inventory reconstruction",
            "seed": seedinfo, "adsorbate_order": names,
            "coverage": {name: mp.nstr(cv[i], 80) for i, name in enumerate(names)},
            "vacant": mp.nstr(vacancy, 80), "I1": mp.nstr(i1, 40), "I2": mp.nstr(i2, 40),
            "surface_residuals": {name: mp.nstr(net[i], 40) for i, name in enumerate(names)},
            "reaction_net_rates": [mp.nstr(x, 40) for x in reaction_rates],
            "reaction_labels": [str(r) for r in solver.elementary_rxns],
            "gas_tof": {name: mp.nstr(value, 40) for name, value in zip(solver.gas_names, tof)},
            "C6H12O_net_tof": mp.nstr(tof[solver.gas_names.index("C6H12O_g")], 80),
            "log10_C6H12O_net_tof": mp.nstr(mp.log10(tof[solver.gas_names.index("C6H12O_g")]), 60),
            "tof_unit_convention": "s^-1 as CatMAP turnover_frequency under the configured rate prefactors",
            "max_independent_residual_abs": mp.nstr(max(abs(net[i]) for i in independent_rows), 40),
            "relative_residual_to_max_reaction_rate": mp.nstr(relative_resid, 40),
            "site_sum": mp.nstr(sum(cv) + vacancy, 40),
            "max_rate_abs": mp.nstr(max_rate, 40), "cycle_max_relative_error": mp.nstr(cycle_rel, 40),
            "surface_stoichiometry_balance_abs": mp.nstr(surface_balance_abs, 40),
            "gas_stoichiometry_reconstruction_abs": mp.nstr(gas_balance_abs, 40),
            "stoichiometry_flux_relative_error": mp.nstr(stoich_flux_rel, 40),
            "surface_flux_reconstruction": {name: mp.nstr(value, 30) for name, value in zip(names, surface_from_rates)},
            "gas_flux_reconstruction": {name: mp.nstr(value, 30) for name, value in zip(solver.gas_names, gas_from_rates)},
            "zero_inventory_local_stability_eigenvalues_per_s": eigenvalues,
            "locally_stable_in_zero_inventory_subspace": stable_local,
            "stability_diagnostic_error": stability_error,
            "scaled_residual_inf": mp.nstr(scaled_resid, 40), "newton_iterations": iteration_count,
            "accepted_numeric_seed": accepted_seed, "solver_converged_and_validated": bool(solver_success),
            "root_logits": [mp.nstr(x, 40) for x in root],
        }
        (WORK / f"constrained_{SURFACE}_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        result["cross_precision_replay"] = cross_precision
        print("result:", json.dumps({k: result[k] for k in ("coverage", "vacant", "I1", "I2", "gas_tof", "max_independent_residual_abs", "relative_residual_to_max_reaction_rate", "site_sum", "scaled_residual_inf", "newton_iterations", "solver_converged_and_validated")}, indent=2))
        (WORK / f"constrained_{SURFACE}_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

        # Conditional DRC finite differences: perturb CatMAP's parsed energy
        # parameter vector, while re-solving from the accepted root inside the
        # same exact inventory chart. This intentionally does not call the
        # stock unconstrained get_rate_control() path.
        baseline_params = list(rxn_params)
        baseline_root = list(root)
        target_gas = solver.gas_names.index("C6H12O_g")
        kbt = mp.mpf("8.617332478e-5") * mp.mpf("333.15")
        drc_rows = []
        attempt_rows = []
        for parameter_index, parameter_name in enumerate(solver.adsorbate_names + solver.transition_state_names):
            tof_by_delta = {}
            continuation_root = {"minus": list(baseline_root), "plus": list(baseline_root)}
            for delta_text in ("0.005", "0.010"):
                delta = mp.mpf(delta_text)
                perturbed_log_tof = {}
                pair_accepted = True
                for sign in (-1, 1):
                    label = "minus" if sign < 0 else "plus"
                    rxn_params = list(baseline_params)
                    rxn_params[parameter_index] += sign * delta
                    for memo_name in ("_rate_constant_memoize", "_steady_state_memoize"):
                        memo = getattr(solver, memo_name, None)
                        if hasattr(memo, "clear"):
                            memo.clear()

                    def pert_residual(z):
                        _, _, net_state = evaluate(z)
                        return tuple(net_state[i] / row_scales[j] for j, i in enumerate(independent_rows))

                    def pert_jacobian(z):
                        z = list(z)
                        return mp.matrix([[mp.diff(
                            lambda q: pert_residual([q if kk == col else z[kk] for kk in range(len(z))])[row], z[col]
                        ) for col in range(len(indep))] for row in range(len(independent_rows))])

                    iter_pert = NewtonRoot(
                        pert_residual, continuation_root[label], mp.matrix, mp.mpf, mp.lu_solve,
                        J=pert_jacobian, norm=lambda x: mp.norm(x, "inf"), verbose=0,
                    )
                    root_pert = list(baseline_root)
                    norm_pert = mp.inf
                    iterations_pert = 0
                    err_pert = None
                    try:
                        for iterations_pert, (x_pert, norm_pert) in enumerate(iter_pert, start=1):
                            root_pert = list(x_pert)
                            if norm_pert < mp.mpf("1e-50") or iterations_pert >= 100:
                                break
                    except Exception as exc:
                        err_pert = f"{type(exc).__name__}: {exc}"
                        norm_pert = mp.norm(mp.matrix(pert_residual(root_pert)), "inf")
                    cv_pert, vac_pert, net_pert, rates_pert, tof_pert = evaluate(root_pert, want_rates=True)
                    abs_pert = max(abs(net_pert[i]) for i in independent_rows)
                    rel_pert = abs_pert / max((abs(x) for x in rates_pert), default=mp.mpf("1e-300"))
                    inv1_pert = cv_pert[idx["Ha_s"]] - cv_pert[idx["OOH_s"]] - cv_pert[idx["OH_s"]] - cv_pert[idx["O_s"]]
                    inv2_pert = cv_pert[idx["Hb_s"]] - cv_pert[idx["O_s"]]
                    tof_value = tof_pert[target_gas]
                    accepted_pert = (norm_pert < mp.mpf("1e-45") and rel_pert < mp.mpf("1e-20")
                                     and abs(inv1_pert) < mp.mpf("1e-300") and abs(inv2_pert) < mp.mpf("1e-300")
                                     and vac_pert > 0 and all(v >= 0 for v in cv_pert) and tof_value > 0)
                    if accepted_pert:
                        continuation_root[label] = list(root_pert)
                    pair_accepted = pair_accepted and accepted_pert
                    log_value = mp.log(tof_value) if tof_value > 0 else mp.nan
                    perturbed_log_tof[label] = log_value
                    attempt_rows.append({
                        "parameter": parameter_name, "delta_eV": delta_text, "direction": label,
                        "accepted_same_inventory_root": bool(accepted_pert), "iterations": iterations_pert,
                        "scaled_residual_inf": mp.nstr(norm_pert, 20),
                        "relative_residual": mp.nstr(rel_pert, 20),
                        "I1": mp.nstr(inv1_pert, 12), "I2": mp.nstr(inv2_pert, 12),
                        "branch_logit_distance_from_baseline": mp.nstr(max(abs(a-b) for a, b in zip(root_pert, baseline_root)), 16),
                        "TOF_C6H12O": mp.nstr(tof_value, 30), "error": err_pert or "",
                    })
                if pair_accepted:
                    tof_by_delta[delta_text] = perturbed_log_tof
            valid_deltas = [d for d, values in tof_by_delta.items() if set(values) == {"minus", "plus"}]
            derivs = {}
            for d in valid_deltas:
                delta = mp.mpf(d)
                derivs[d] = -(kbt) * (tof_by_delta[d]["plus"] - tof_by_delta[d]["minus"]) / (2 * delta)
            drc_rows.append({
                "parameter": parameter_name,
                "DRC_delta_0p005": mp.nstr(derivs["0.005"], 30) if "0.005" in derivs else "",
                "DRC_delta_0p010": mp.nstr(derivs["0.010"], 30) if "0.010" in derivs else "",
                "delta_sensitivity_difference": mp.nstr(abs(derivs["0.005"]-derivs["0.010"]), 15)
                    if set(derivs) == {"0.005", "0.010"} else "",
                "status": "PROVISIONAL_SAME_INVENTORY_DRC" if set(derivs) == {"0.005", "0.010"} else "FAILED_OR_INCOMPLETE",
            })
        with (WORK / "conditional_drc.csv").open("w", encoding="utf-8-sig", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(drc_rows[0]))
            writer.writeheader()
            writer.writerows(drc_rows)
        (WORK / "conditional_drc_attempts.json").write_text(json.dumps(attempt_rows, indent=2), encoding="utf-8")
        (WORK / "conditional_drc_method.txt").write_text(
            "Conditional finite-difference DRC from CatMAP 0.3.1 native rates, with each +/- energy perturbation "
            "re-solved from the accepted root in the exact I1=I2=0 chart. This is not the stock "
            "MeanFieldSolver.get_rate_control output; stock DRC invokes unconstrained get_turnover_frequency/coverage. "
            "All values remain PROVISIONAL because NEB-to-independent-state physical mapping and global reachability "
            "are not established. Both perturbation sizes 0.005 and 0.010 eV are reported for stability checking.\n",
            encoding="utf-8",
        )
    finally:
        os.chdir(previous)


if __name__ == "__main__":
    main()
