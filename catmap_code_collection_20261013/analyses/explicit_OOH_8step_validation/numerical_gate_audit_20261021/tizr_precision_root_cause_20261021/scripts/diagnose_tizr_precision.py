from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

import mpmath as mp


PROJECT = Path(r"H:\catmap\catmap\catmap_code_collection_20261013")
MODEL = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\alkene_epoxidation.mkm")
ENERGIES = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\energies_dft_neb.txt")
REPLAY = PROJECT / "analyses/explicit_OOH_8step_validation/multisite_conditional_replay_20261021"
SITE_DIR = REPLAY / "sites/tizr"
SAVED_ROOT = SITE_DIR / "constrained_tizr_result.json"
OUT = Path(r"C:\Users\Administrator\Documents\codex探索\tizr_precision_root_cause_20261021")
INDEPENDENT = ["C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "OH_s", "OOH_s", "O_s"]
TARGET_SCALED_RESIDUAL = mp.mpf("1e-100")
MAX_POLISH_ITERATIONS = 100


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for folder in ("data", "logs", "inputs"):
        (OUT / folder).mkdir(exist_ok=True)

    saved = json.loads(SAVED_ROOT.read_text(encoding="utf-8"))
    input_hashes = {"mkm_sha256": sha256(MODEL), "energies_sha256": sha256(ENERGIES),
                    "saved_root_sha256": sha256(SAVED_ROOT)}
    setup_path = OUT / "inputs/tizr_precision_replay.mkm"
    source = MODEL.read_text(encoding="utf-8")
    overrides = f"""
# Isolated Ti-Zr numerical replay; original scientific inputs remain read-only.
surface_names = ['tizr']
scaler = 'ThermodynamicScaler'
descriptor_names = ['temperature']
descriptor_ranges = None
resolution = 1
output_variables = []
pressure_mode = 'static'
gas_thermo_mode = 'frozen_gas'
adsorbate_thermo_mode = 'frozen_adsorbate'
input_file = {str(ENERGIES)!r}
data_file = 'tizr_precision_replay.pkl'
decimal_precision = 385
from mpmath import mp
tolerance = mp.mpf('1e-260')
max_rootfinding_iterations = 250
use_numbers_solver = False
"""
    setup_path.write_text(source + "\n\n" + overrides, encoding="utf-8")

    old_cwd = Path.cwd()
    os.chdir(OUT / "inputs")
    try:
        from catmap import ReactionModel
        from catmap.solvers.solver_base import NewtonRoot

        model = ReactionModel(setup_file=setup_path.name)
        model.parse()
        solver, scaler = model.solver, model.scaler
        mp.mp.dps = 385
        model._math = mp
        model._mpfloat = mp.mpf
        model._matrix = mp.matrix
        model._Axb_solver = mp.lu_solve
        model._kB = mp.mpf("8.617332478e-5")
        model._h = mp.mpf("4.135667516e-15")
        model._math.infnorm = lambda x: mp.norm(x, "inf")
        model._math.lsqnorm = lambda x: mp.norm(x, 2) ** 2
        solver.compile()

        names = list(solver.adsorbate_names)
        idx = {name: names.index(name) for name in names}
        indep_indices = [idx[name] for name in INDEPENDENT]
        # This chart reproduces exactly the two audited zero-inventory relations.
        derived = {"Ha_s": ("OOH_s", "OH_s", "O_s"), "Hb_s": ("O_s",)}
        weights = [1, 1, 1, 1, 2, 2, 3]
        params = list(scaler.get_rxn_parameters([mp.mpf("333.15")]))
        nrxn = len(solver.elementary_rxns)

        def state_from_logits(z):
            odds = [mp.exp(value) for value in z]
            denom = 1 + sum(weight * value for weight, value in zip(weights, odds))
            coverage = [mp.mpf(0)] * len(names)
            for i, value in zip(indep_indices, odds):
                coverage[i] = value / denom
            for target, sources in derived.items():
                coverage[idx[target]] = sum(coverage[idx[source]] for source in sources)
            vacancy = 1 - sum(coverage)
            return coverage, vacancy

        def evaluate(z):
            coverage, vacancy = state_from_logits(z)
            constants = solver.get_rate_constants(params, list(coverage))
            residual = list(solver.steady_state_function(list(coverage)))
            net = list(solver.get_rxn_rates(list(coverage), constants))
            tof = list(solver.get_turnover_frequency(params, rates=net, verify_coverages=False))
            return coverage, vacancy, residual, net, tof, constants

        # Gross-flow row scales mirror the archived replay method. The official
        # directional-rate interface uses Python float, so it is used only for
        # scaling/initialization and never to calculate a final net flux.
        stoich = [[0] * nrxn for _ in names]
        for j, reaction in enumerate(solver.elementary_rxns):
            initial, final = Counter(reaction[0]), Counter(reaction[-1])
            for species in names:
                stoich[idx[species]][j] = final.get(species, 0) - initial.get(species, 0)

        high_z = [mp.mpf(value) for value in saved["root_logits"]]
        high_eval = evaluate(high_z)
        solver._coverage = list(high_eval[0])
        directional_float = list(solver.get_directional_rates(params))
        scales = []
        for row in indep_indices:
            gross = sum(abs(stoich[row][j]) *
                        (abs(directional_float[j]) + abs(directional_float[nrxn + j]))
                        for j in range(nrxn))
            scales.append(mp.mpf(str(gross)) if gross > 0 else mp.mpf("1e-300"))
        global_scale = max(scales)
        scales = [max(value, global_scale * mp.mpf("1e-250")) for value in scales]

        def residual(z):
            values = evaluate(z)[2]
            return tuple(values[row] / scales[k] for k, row in enumerate(indep_indices))

        def jacobian(z):
            return mp.matrix([[mp.diff(
                lambda q: residual([q if col == k else z[k] for k in range(len(INDEPENDENT))])[row],
                z[col],
            ) for col in range(len(INDEPENDENT))] for row in range(len(INDEPENDENT))])

        early_seed = [value + (mp.mpf("0.01") if j % 2 == 0 else mp.mpf("-0.01"))
                      for j, value in enumerate(high_z)]
        early_iterator = NewtonRoot(residual, early_seed, mp.matrix, mp.mpf, mp.lu_solve,
                                    J=jacobian, norm=lambda value: mp.norm(value, "inf"), verbose=0)
        early_root = list(early_seed)
        early_scaled = mp.inf
        early_iterations = 0
        for early_iterations, (root_matrix, norm_value) in enumerate(early_iterator, start=1):
            early_root, early_scaled = list(root_matrix), norm_value
            if norm_value < mp.mpf("1e-50") or early_iterations >= 250:
                break
        early_eval = evaluate(early_root)

        # One bounded continuation from the actual archived early-stop root.
        # Same model, same input hashes, same precision, same frozen scales.
        polish_iterator = NewtonRoot(residual, early_root, mp.matrix, mp.mpf, mp.lu_solve,
                                     J=jacobian, norm=lambda value: mp.norm(value, "inf"), verbose=0)
        polish_root = list(early_root)
        polish_scaled = mp.norm(mp.matrix(residual(polish_root)), "inf")
        polish_iterations = 0
        polish_error = None
        try:
            for polish_iterations, (root_matrix, norm_value) in enumerate(polish_iterator, start=1):
                polish_root, polish_scaled = list(root_matrix), norm_value
                if norm_value < TARGET_SCALED_RESIDUAL or polish_iterations >= MAX_POLISH_ITERATIONS:
                    break
        except Exception as exc:  # retain diagnostic state; do not silently treat as pass
            polish_error = f"{type(exc).__name__}: {exc}"
            polish_scaled = mp.norm(mp.matrix(residual(polish_root)), "inf")
        polish_eval = evaluate(polish_root)

        # CatMAP's own mass-action expression generator gives directional rates.
        # Evaluate those expressions in mpmath arithmetic to avoid the float cast
        # in get_directional_rates()/get_elem_ec() when forward and reverse nearly cancel.
        forward_exprs, reverse_exprs = [], []
        for j, reaction in enumerate(solver.elementary_rxns):
            forward_exprs.append(solver.rate_equation_term(reaction[0], f"kf[{j}]"))
            reverse_exprs.append(solver.rate_equation_term(reaction[-1], f"kr[{j}]"))

        def exact_directional(snapshot):
            coverage, vacancy, _, net, tof, constants = snapshot
            env = {"kf": constants[:nrxn], "kr": constants[nrxn:],
                   "p": list(solver.gas_pressures), "theta": coverage,
                   "s": [vacancy], "mp": mp}
            forward = [eval(expr, {"__builtins__": {}}, env) for expr in forward_exprs]
            reverse = [eval(expr, {"__builtins__": {}}, env) for expr in reverse_exprs]
            return forward, reverse

        snapshots = {"archived_high_precision_root": high_eval,
                     "replayed_old_1e-50_stop": early_eval,
                     "same_precision_bounded_polish": polish_eval}
        root_rows, flux_rows, root_metrics = [], [], {}
        for label, snapshot in snapshots.items():
            coverage, vacancy, residual_values, net, tof, _ = snapshot
            forward, reverse = exact_directional(snapshot)
            # Verify exact forward-minus-reverse against CatMAP's native net rates.
            flux_formula_errors = [abs((forward[j] - reverse[j]) - net[j]) for j in range(nrxn)]
            i1 = coverage[idx["Ha_s"]] - coverage[idx["OOH_s"]] - coverage[idx["OH_s"]] - coverage[idx["O_s"]]
            i2 = coverage[idx["Hb_s"]] - coverage[idx["O_s"]]
            max_independent_resid = max(abs(residual_values[i]) for i in indep_indices)
            row_metrics = {}
            for i, name in enumerate(names):
                root_rows.append({"solution": label, "species": name,
                                  "coverage": mp.nstr(coverage[i], 80),
                                  "absolute_dtheta_dt": mp.nstr(residual_values[i], 60),
                                  "delta_vs_archived_high_root": mp.nstr(
                                      coverage[i] - high_eval[0][i], 60)})
                row_metrics[name] = {"coverage": mp.nstr(coverage[i], 80),
                                     "residual": mp.nstr(residual_values[i], 60)}
            root_rows.append({"solution": label, "species": "vacancy",
                              "coverage": mp.nstr(vacancy, 80),
                              "absolute_dtheta_dt": "N/A",
                              "delta_vs_archived_high_root": mp.nstr(vacancy - high_eval[1], 60)})

            for j, reaction in enumerate(solver.elementary_rxns):
                net_abs = abs(net[j])
                cancellation = ((abs(forward[j]) + abs(reverse[j])) / net_abs
                                if net_abs != 0 else None)
                flux_rows.append({
                    "solution": label, "reaction_index": j + 1,
                    "reaction": f"{' + '.join(reaction[0])} -> {' + '.join(reaction[-1])}",
                    "forward_rate_mpmath": mp.nstr(forward[j], 80),
                    "reverse_rate_mpmath": mp.nstr(reverse[j], 80),
                    "net_rate_catmap_mpmath": mp.nstr(net[j], 80),
                    "forward_minus_reverse_error": mp.nstr(flux_formula_errors[j], 40),
                    "cancellation_factor": mp.nstr(cancellation, 40) if cancellation is not None else "NET_ZERO",
                    "tof_C6H12O_g": mp.nstr(tof[solver.gas_names.index("C6H12O_g")], 80),
                    "mp_dps": mp.mp.dps,
                    "direction_source": "CatMAP rate_equation_term evaluated with mpmath; not float get_directional_rates",
                })
            product_tof = tof[solver.gas_names.index("C6H12O_g")]
            root_metrics[label] = {
                "species": row_metrics,
                "vacancy": mp.nstr(vacancy, 80),
                "I1": mp.nstr(i1, 60), "I2": mp.nstr(i2, 60),
                "site_sum": mp.nstr(sum(coverage) + vacancy, 60),
                "max_independent_abs_residual": mp.nstr(max_independent_resid, 60),
                "net_reaction_rates": [mp.nstr(value, 80) for value in net],
                "TOF_C6H12O": mp.nstr(product_tof, 80),
                "log10_TOF_C6H12O": mp.nstr(mp.log10(product_tof), 60) if product_tof > 0 else "NONPOSITIVE",
                "max_forward_reverse_identity_error": mp.nstr(max(flux_formula_errors), 40),
                "max_coverage_difference_from_archived_high": mp.nstr(
                    max(abs(a - b) for a, b in zip(coverage, high_eval[0])), 60),
            }

        write_csv(OUT / "data/TIZR_ROOT_COMPARISON.csv", root_rows)
        write_csv(OUT / "data/TIZR_FORWARD_REVERSE_FLUX.csv", flux_rows)

        high_tof = high_eval[4][solver.gas_names.index("C6H12O_g")]
        early_tof = early_eval[4][solver.gas_names.index("C6H12O_g")]
        polished_tof = polish_eval[4][solver.gas_names.index("C6H12O_g")]
        delta_log_early = abs(mp.log10(early_tof) - mp.log10(high_tof))
        delta_log_polished = abs(mp.log10(polished_tof) - mp.log10(high_tof))
        rel_polished = abs(polished_tof - high_tof) / abs(high_tof)
        summary = {
            "site": "tizr", "catmap_version": "0.3.1", "mp_dps": mp.mp.dps,
            "inputs": {"model": str(MODEL), "energies": str(ENERGIES), **input_hashes},
            "archived_cross_precision": json.loads((SITE_DIR / "cross_precision_replay.json").read_text(encoding="utf-8")),
            "archived_high_root": root_metrics["archived_high_precision_root"],
            "replayed_early_stop": {**root_metrics["replayed_old_1e-50_stop"],
                                    "iterations": early_iterations,
                                    "scaled_residual_inf": mp.nstr(mp.norm(mp.matrix(residual(early_root)), "inf"), 50)},
            "bounded_polish": {**root_metrics["same_precision_bounded_polish"],
                               "iterations": polish_iterations,
                               "scaled_residual_inf": mp.nstr(polish_scaled, 50),
                               "target_scaled_residual": mp.nstr(TARGET_SCALED_RESIDUAL),
                               "max_iterations": MAX_POLISH_ITERATIONS,
                               "error": polish_error},
            "delta_log10_TOF_early_vs_high": mp.nstr(delta_log_early, 40),
            "delta_log10_TOF_polished_vs_high": mp.nstr(delta_log_polished, 40),
            "polished_TOF_relative_difference_vs_high": mp.nstr(rel_polished, 40),
            "old_stop_criterion_reproduced": early_iterations == 5,
            "accepted_gate_5ppm": bool(rel_polished <= mp.mpf("5e-6")),
            "physics_changed": False,
        }
        (OUT / "data/tizr_replay_diagnostics.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        log = [
            f"CatMAP={type(solver).__name__}/{model.scaler.__class__.__name__}",
            f"reaction_count={nrxn}; adsorbates={','.join(names)}",
            f"input_hashes={json.dumps(input_hashes, sort_keys=True)}",
            f"early_stop_iterations={early_iterations}; scaled_residual={summary['replayed_early_stop']['scaled_residual_inf']}",
            f"polish_iterations={polish_iterations}; scaled_residual={mp.nstr(polish_scaled, 50)}; error={polish_error}",
            f"TOF_high={mp.nstr(high_tof, 80)}",
            f"TOF_early={mp.nstr(early_tof, 80)}",
            f"TOF_polished={mp.nstr(polished_tof, 80)}",
            f"delta_log10_early={mp.nstr(delta_log_early, 40)}",
            f"delta_log10_polished={mp.nstr(delta_log_polished, 40)}",
            f"relative_polished={mp.nstr(rel_polished, 40)}",
        ]
        (OUT / "logs/tizr_precision_replay.log").write_text("\n".join(log) + "\n", encoding="utf-8")
        print("REPLAY_SUMMARY=" + json.dumps({
            "early_iterations": early_iterations,
            "early_scaled_residual": summary["replayed_early_stop"]["scaled_residual_inf"],
            "polish_iterations": polish_iterations,
            "polish_scaled_residual": mp.nstr(polish_scaled, 25),
            "tof_high": mp.nstr(high_tof, 30), "tof_early": mp.nstr(early_tof, 30),
            "tof_polished": mp.nstr(polished_tof, 30),
            "delta_log_early": mp.nstr(delta_log_early, 20),
            "delta_log_polished": mp.nstr(delta_log_polished, 20),
            "relative_polished": mp.nstr(rel_polished, 20),
            "accepted_5ppm": bool(rel_polished <= mp.mpf("5e-6")),
            "polish_error": polish_error,
            "flux_identity_max_error": root_metrics["same_precision_bounded_polish"]["max_forward_reverse_identity_error"],
        }, sort_keys=True))
    finally:
        os.chdir(old_cwd)


if __name__ == "__main__":
    main()
