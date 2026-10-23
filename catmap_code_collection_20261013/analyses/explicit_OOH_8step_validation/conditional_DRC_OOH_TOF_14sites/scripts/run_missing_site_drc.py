from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import pathlib
import sys
import time
from collections import Counter

ROOT = pathlib.Path(r"H:\catmap\catmap\catmap_code_collection_20261013")
ANALYSIS = ROOT / "analyses" / "explicit_OOH_8step_validation"
REPLAY = ANALYSIS / "multisite_conditional_replay_20261021"
MODEL = pathlib.Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\alkene_epoxidation.mkm")
ENERGIES = pathlib.Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\energies_dft_neb.txt")
EXPECTED_MODEL_SHA = "86115cefb7784fe28ab519271cdfaf0e308b6f6e3a8308b1e90e4318cb497cfe"
EXPECTED_ENERGY_SHA = "55ad4e94d9c603b6705ec8c92630c9a1ceff5859cf69865145152375e203ca30"
ADS_ORDER = ["C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "Ha_s", "Hb_s", "OH_s", "OOH_s", "O_s"]
FREE = ["C6H12O_s", "C6H12_s", "H2O2_s", "H2O_s", "OH_s", "OOH_s", "O_s"]
DERIVED = {"Ha_s": ("OOH_s", "OH_s", "O_s"), "Hb_s": ("O_s",)}
H_STEPS = ("0.005", "0.010")


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def root_data(surface: str):
    p = REPLAY / "sites" / surface / f"constrained_{surface}_result.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8-sig"))
        return {k: str(v) for k, v in d["coverage"].items()}, d.get("root_logits"), d.get("C6H12O_net_tof"), p
    if surface == "timo":
        p = ANALYSIS / "three_site_root_recovery_20261021" / "data" / "timo_recovery.json"
        d = json.loads(p.read_text(encoding="utf-8-sig"))
        return {k: str(v) for k, v in d["coverage"].items()}, d.get("root_logits"), d.get("tof_C6H12O_s-1"), p
    if surface in ("tiw", "tita"):
        p = ANALYSIS / "tiw_tita_targeted_solver_recovery" / "data" / "TIW_TITA_SOLVER_RESULTS.csv"
        with p.open(encoding="utf-8-sig", newline="") as f:
            row = next(r for r in csv.DictReader(f) if r["surface"] == surface)
        cov = {name.removeprefix("theta_"): row[name] for name in row if name.startswith("theta_")}
        return cov, None, row["TOF_C6H12O_s-1"], p
    raise FileNotFoundError(f"No approved latest root package for {surface}")


def setup_catmap(surface: str, outdir: pathlib.Path):
    # Existing probe contains the original eight-step model settings for this surface.
    setup = REPLAY / "sites" / surface / f"{surface}_static_probe.mkm"
    if not setup.exists():
        raise FileNotFoundError(setup)
    import mpmath as mp
    import catmap
    from catmap import ReactionModel

    if getattr(catmap, "__version__", None) != "0.3.1":
        raise RuntimeError(f"CatMAP version mismatch: {getattr(catmap, '__version__', None)}")
    previous = pathlib.Path.cwd()
    work = outdir / "runs" / surface
    work.mkdir(parents=True, exist_ok=True)
    os.chdir(work)
    model = ReactionModel(setup_file=str(setup))
    mp.mp.dps = max(180, int(model.decimal_precision) + 25)
    model._math = mp.mp
    model._mpfloat = mp.mpf
    model._matrix = mp.matrix
    model._kB = mp.mpf("8.617332478e-5")
    model._h = mp.mpf("4.135667516e-15")
    model._math.infnorm = lambda x: mp.norm(x, "inf")
    model._math.lsqnorm = lambda x: mp.norm(x, 2) ** 2
    model.solver.compile()
    # Contract checks against the parsed model; no mapper run and no base root solve.
    solver = model.solver
    if len(model.rxn_expressions) != 8 or len(solver.adsorbate_names) != 9 or len(solver.transition_state_names) != 2:
        raise RuntimeError("MODEL_CONTRACT=FAIL: parsed reaction/state counts differ from approved original model")
    if getattr(model, "pressure_mode", None) != "static":
        raise RuntimeError("MODEL_CONTRACT=FAIL: pressure_mode is not static")
    if abs(float(solver.temperature) - 333.15) > 1e-10:
        raise RuntimeError(f"MODEL_CONTRACT=FAIL: temperature {solver.temperature}")
    if list(solver.adsorbate_names) != ADS_ORDER:
        raise RuntimeError(f"MODEL_CONTRACT=FAIL: adsorbate order changed: {solver.adsorbate_names}")
    params = list(model.scaler.get_rxn_parameters([mp.mpf("333.15")]))
    os.chdir(previous)
    return mp, model, solver, params, previous


def run_site(surface: str, outdir: pathlib.Path):
    start_time = time.time()
    mp, model, solver, base_params, previous = setup_catmap(surface, outdir)
    try:
        mp.mp.dps = max(180, int(model.decimal_precision) + 25)
        idx = {n: i for i, n in enumerate(solver.adsorbate_names)}
        if solver.adsorbate_names != tuple(ADS_ORDER):
            raise RuntimeError("Unexpected CatMAP species order")
        cov_dict, saved_logits, saved_tof, root_path = root_data(surface)
        # Reconstruct a precise root chart from stored coverages. Stored JSON/CSV
        # values are finite precision; check their projection, then enforce the
        # approved inventory identities algebraically in the solver coordinates.
        cv_saved = [mp.mpf(cov_dict[n]) for n in ADS_ORDER]
        serialized_i1 = cv_saved[idx["Ha_s"]] - cv_saved[idx["OOH_s"]] - cv_saved[idx["OH_s"]] - cv_saved[idx["O_s"]]
        serialized_i2 = cv_saved[idx["Hb_s"]] - cv_saved[idx["O_s"]]
        if abs(serialized_i1) > mp.mpf("1e-60") or abs(serialized_i2) > mp.mpf("1e-60"):
            raise ValueError(f"stored root projection inconsistent with zero-inventory class: I1={serialized_i1}, I2={serialized_i2}")
        weights = [1, 1, 1, 1, 2, 2, 3]
        free_values = [mp.mpf(cov_dict[n]) for n in FREE]
        vacancy = 1 - sum(free_values) - free_values[5] - free_values[4] - 2 * free_values[6]
        if vacancy <= 0 or any(v <= 0 for v in free_values):
            raise ValueError(f"saved root projection outside positive zero-inventory chart: vacancy={vacancy}")
        z0 = [mp.log(v / vacancy) for v in free_values]
        free_idx = [idx[n] for n in FREE]

        def state(z):
            odds = [mp.exp(q) for q in z]
            den = 1 + sum(w * q for w, q in zip(weights, odds))
            vals = [mp.mpf("0")] * len(ADS_ORDER)
            for j, q in zip(free_idx, odds):
                vals[j] = q / den
            for target, sources in DERIVED.items():
                vals[idx[target]] = sum(vals[idx[s]] for s in sources)
            vac = 1 - sum(vals)
            return vals, vac

        cv, vacancy_chart = state(z0)
        i1 = cv[idx["Ha_s"]] - cv[idx["OOH_s"]] - cv[idx["OH_s"]] - cv[idx["O_s"]]
        i2 = cv[idx["Hb_s"]] - cv[idx["O_s"]]
        free_idx = [idx[n] for n in FREE]

        n_rxn = len(solver.elementary_rxns)
        stoich = [[0] * n_rxn for _ in ADS_ORDER]
        for j, rxn in enumerate(solver.elementary_rxns):
            ci, cf = Counter(rxn[0]), Counter(rxn[-1])
            for sp in ADS_ORDER:
                stoich[idx[sp]][j] = cf.get(sp, 0) - ci.get(sp, 0)
        independent = [idx[n] for n in FREE]

        def eval_at(z, params, rates=False):
            vals, vac = state(z)
            ks = solver.get_rate_constants(list(params), list(vals))
            net = list(solver.steady_state_function(list(vals)))
            if not rates:
                return vals, vac, net
            rxrates = list(solver.get_rxn_rates(list(vals), ks))
            tof = list(solver.get_turnover_frequency(list(params), rates=rxrates, verify_coverages=False))
            return vals, vac, net, rxrates, tof

        # CatMAP native directional rates set per-species gross-flow residual scales;
        # subtraction of tiny net rates is still done by mpmath on the native RHS.
        solver._coverage = list(cv)
        directional = list(solver.get_directional_rates(list(base_params)))
        scales = []
        for row in independent:
            gross = sum(abs(stoich[row][j]) * (abs(mp.mpf(str(directional[j]))) + abs(mp.mpf(str(directional[n_rxn + j])))) for j in range(n_rxn))
            scales.append(gross if gross > 0 else mp.mpf("1e-300"))
        global_scale = max(scales)
        scales = [max(q, global_scale * mp.mpf("1e-250")) for q in scales]

        def solve_perturbed(params, seed, label):
            for memo_name in ("_rate_constant_memoize", "_steady_state_memoize"):
                memo = getattr(solver, memo_name, None)
                if hasattr(memo, "clear"):
                    memo.clear()

            def residual(z):
                _, _, net = eval_at(z, params)
                return tuple(net[i] / scales[k] for k, i in enumerate(independent))

            def jac(z):
                return mp.matrix([[mp.diff(lambda x: residual([x if c == col else z[c] for c in range(len(z))])[row], z[col])
                                   for col in range(len(z))] for row in range(len(independent))])

            from catmap.solvers.solver_base import NewtonRoot
            iterator = NewtonRoot(residual, list(seed), mp.matrix, mp.mpf, mp.lu_solve,
                                  J=jac, norm=lambda x: mp.norm(x, "inf"), verbose=0)
            z = list(seed)
            norm = mp.inf
            iterations = 0
            error = ""
            try:
                for iterations, (x, norm) in enumerate(iterator, start=1):
                    z = list(x)
                    if norm < mp.mpf("1e-50") or iterations >= 120:
                        break
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
                norm = mp.norm(mp.matrix(residual(z)), "inf")
            vals, vac, net, rxrates, tof = eval_at(z, params, rates=True)
            inv1 = vals[idx["Ha_s"]] - vals[idx["OOH_s"]] - vals[idx["OH_s"]] - vals[idx["O_s"]]
            inv2 = vals[idx["Hb_s"]] - vals[idx["O_s"]]
            loc = max(abs(net[i]) / scales[k] for k, i in enumerate(independent))
            maxrate = max((abs(x) for x in rxrates), default=mp.mpf("0"))
            cycle = max(abs(x - rxrates[3]) for x in rxrates) / maxrate if maxrate else mp.inf
            product = tof[solver.gas_names.index("C6H12O_g")]
            ok = (loc < mp.mpf("1e-12") and vac > 0 and all(x >= 0 for x in vals)
                  and abs(inv1) < mp.mpf("1e-100") and abs(inv2) < mp.mpf("1e-100")
                  and cycle < mp.mpf("1e-12") and product > 0)
            return {
                "surface": surface, "state": label, "accepted": bool(ok), "iterations": iterations,
                "scaled_local_residual_max": mp.nstr(loc, 35), "site_balance": mp.nstr(sum(vals) + vac, 35),
                "I1": mp.nstr(inv1, 25), "I2": mp.nstr(inv2, 25), "vacancy": mp.nstr(vac, 35),
                "TOF_C6H12O_s-1": mp.nstr(product, 70), "log10_TOF": mp.nstr(mp.log10(product), 50),
                "net_rates": [mp.nstr(x, 45) for x in rxrates], "cycle_flux_relative_deviation": mp.nstr(cycle, 25),
                "coverages": {n: mp.nstr(vals[i], 45) for i, n in enumerate(ADS_ORDER)},
                "root_chart": [mp.nstr(x, 45) for x in z], "chart_distance_from_baseline_inf": mp.nstr(max(abs(a-b) for a,b in zip(z,z0)), 20),
                "error": error,
            }

        # Baseline is not solved here; saved accepted root is only evaluated for provenance.
        base_vals, base_vac, base_net, base_rx, base_tof = eval_at(z0, base_params, rates=True)
        base_product = base_tof[solver.gas_names.index("C6H12O_g")]
        if saved_tof and abs(mp.log10(base_product) - mp.log10(mp.mpf(saved_tof))) > mp.mpf("5e-6"):
            raise RuntimeError(f"baseline replay from stored root differs >5 ppm in log10 TOF: {base_product} vs {saved_tof}")

        rows, attempts = [], []
        params0 = list(base_params)
        state_names = list(solver.adsorbate_names) + list(solver.transition_state_names)
        kbt = mp.mpf("8.617332478e-5") * mp.mpf("333.15")
        for pidx, name in enumerate(state_names):
            tof_logs = {}
            for htxt in H_STEPS:
                h = mp.mpf(htxt)
                pair = {}
                for sign, direction in ((-1, "minus"), (1, "plus")):
                    pars = list(params0)
                    pars[pidx] += sign * h
                    result = solve_perturbed(pars, z0, f"{name}_{htxt}_{direction}")
                    result.update({"perturbed_state": name, "delta_eV": htxt, "direction": direction})
                    attempts.append(result)
                    if result["accepted"]:
                        pair[direction] = mp.mpf(result["log10_TOF"]) * mp.log(10)
                if len(pair) == 2:
                    tof_logs[htxt] = pair
            drcs = {}
            for htxt in H_STEPS:
                if htxt in tof_logs:
                    h = mp.mpf(htxt)
                    drcs[htxt] = -kbt * (tof_logs[htxt]["plus"] - tof_logs[htxt]["minus"]) / (2*h)
            rows.append({
                "surface_name": surface, "perturbed_state": name,
                "state_type": "stable_adsorbate" if name in solver.adsorbate_names else "explicit_transition_state",
                "baseline_TOF_s-1": mp.nstr(base_product, 70),
                "DRC_h0.005": mp.nstr(drcs["0.005"], 40) if "0.005" in drcs else "",
                "DRC_h0.010": mp.nstr(drcs["0.010"], 40) if "0.010" in drcs else "",
                "absolute_h_difference": mp.nstr(abs(drcs["0.005"] - drcs["0.010"]), 30) if len(drcs)==2 else "",
                "accepted_perturbations": sum(1 for a in attempts if a["perturbed_state"] == name and a["accepted"]),
                "expected_perturbations": 4, "DRC_status": "PASS_SAME_INVENTORY_CONDITIONAL" if len(drcs)==2 else "INCOMPLETE",
                "root_source": str(root_path), "source_model_SHA256": sha256(MODEL), "source_energy_SHA256": sha256(ENERGIES),
                "scientific_status": "PROVISIONAL",
            })
            print(f"{surface} {name}: {rows[-1]['DRC_status']} {rows[-1]['DRC_h0.005']} {rows[-1]['DRC_h0.010']}", flush=True)
        site_dir = outdir / "runs" / surface
        (site_dir / "drc_perturb_attempts.json").write_text(json.dumps(attempts, indent=2), encoding="utf-8")
        with (site_dir / "conditional_drc.csv").open("w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
        (site_dir / "site_replay_metadata.json").write_text(json.dumps({
            "surface": surface, "CatMAP_version": "0.3.1", "base_root_was_reused_not_solved": True,
            "root_source": str(root_path), "baseline_TOF_replayed_from_saved_coverage": mp.nstr(base_product,70),
            "saved_baseline_TOF": saved_tof, "model_sha256": sha256(MODEL), "energies_sha256": sha256(ENERGIES),
            "n_adsorbates": len(solver.adsorbate_names), "n_transition_states": len(solver.transition_state_names),
            "n_parsed_elementary_reactions": len(solver.elementary_rxns), "I1": mp.nstr(i1,20), "I2": mp.nstr(i2,20),
            "baseline_vacancy": mp.nstr(base_vac,45), "baseline_net_rates": [mp.nstr(x,40) for x in base_rx],
            "serialized_root_I1_projection": mp.nstr(serialized_i1, 30), "serialized_root_I2_projection": mp.nstr(serialized_i2, 30),
            "elapsed_seconds": time.time()-start_time, "all_44_perturbations_accepted": len(attempts)==44 and all(x["accepted"] for x in attempts),
            "stock_get_rate_control_used": False,
            "reason_stock_method_not_used": "Stock method differentiates get_turnover_frequency at current solver state and does not re-solve on the exact I1=I2=0 manifold.",
        }, indent=2), encoding="utf-8")
        return rows, attempts
    finally:
        os.chdir(previous)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--surfaces", nargs="+", choices=["timo", "tiw", "tita"], required=True)
    args = ap.parse_args()
    out = pathlib.Path(args.outdir)
    if not out.exists():
        raise FileNotFoundError(f"Output directory must be created explicitly and must be new: {out}")
    mh, eh = sha256(MODEL), sha256(ENERGIES)
    if mh != EXPECTED_MODEL_SHA or eh != EXPECTED_ENERGY_SHA:
        raise RuntimeError(f"MODEL_CONTRACT=FAIL: source hashes differ {mh} {eh}")
    all_rows, all_attempts = [], []
    for surface in args.surfaces:
        rows, attempts = run_site(surface, out)
        all_rows.extend(rows); all_attempts.extend(attempts)
        (out / "logs" / f"{surface}_RUN_COMPLETE.txt").write_text(
            f"surface={surface}\nperturbations={len(attempts)}\naccepted={sum(bool(x['accepted']) for x in attempts)}\n", encoding="utf-8")
    (out / "data" / "NEW_SITE_DRC.json").write_text(json.dumps({"model_contract": "PASS", "new_sites": args.surfaces, "results": all_rows}, indent=2), encoding="utf-8")
    print("MODEL_CONTRACT=PASS; new-site DRC rows:", len(all_rows), "perturbations:", len(all_attempts), "accepted:", sum(bool(x["accepted"]) for x in all_attempts), flush=True)


if __name__ == "__main__":
    main()
