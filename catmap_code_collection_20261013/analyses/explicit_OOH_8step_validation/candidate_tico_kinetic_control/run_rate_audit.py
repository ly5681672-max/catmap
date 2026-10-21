"""Audit CatMAP 0.3.1 directional rates at the archived provisional root.

Forward/reverse mass-action factors follow CatMAP's compiled elementary
reaction lists. CatMAP's public get_directional_rates currently casts each
directional flux to Python float; this script records that API output and
reconstructs the same terms with the official high-precision rate constants
and mpmath before subtracting them.
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--coverages", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--barriers", type=Path, required=True, help="candidate's audited reaction free-energy/barrier table")
    ap.add_argument("--precision", type=int, default=180)
    a = ap.parse_args()
    model_path, cov_path, out = a.model.resolve(), a.coverages.resolve(), a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    old = Path.cwd(); os.chdir(model_path.parent)
    try:
        import catmap
        from catmap import ReactionModel
        model = ReactionModel(setup_file=model_path.name)
        model.run(); solver = model.solver; solver.compile(); mp.mp.dps = a.precision
        with cov_path.open(encoding="utf-8-sig", newline="") as f:
            coverage_rows = list(csv.DictReader(f))
        cvmap = {r["species"]: mp.mpf(r["coverage"]) for r in coverage_rows}
        ads = [cvmap[n] for n in solver.adsorbate_names]
        vacancy = cvmap["*"]
        if abs(sum(ads) + vacancy - 1) > mp.mpf("1e-60"):
            raise RuntimeError("input coverage violates site balance")
        if abs(cvmap["Ha_s"] - cvmap["OOH_s"]) > mp.mpf("1e-60"):
            raise RuntimeError("input coverage violates zero-inventory class")
        solver._coverage = ads
        rxn_parameters = list(model.scaler.get_rxn_parameters([mp.mpf("333.15"), mp.mpf("1.0")]))
        constants = list(solver.get_rate_constants(rxn_parameters, ads))
        n = len(solver.elementary_rxns)
        if n != 9 or len(constants) != 2 * n:
            raise RuntimeError(f"unexpected reaction/constants counts: {n}, {len(constants)}")
        official_dirs = list(solver.get_directional_rates(rxn_parameters))
        native_nets = list(solver.get_rxn_rates(ads, constants))
        gas_p = {name: mp.mpf(str(value)) for name, value in zip(solver.gas_names, solver.gas_pressures)}
        ads_map = dict(zip(solver.adsorbate_names, ads))
        ads_map.update({"s": vacancy, "*_s": vacancy})
        barrier_rows = {int(r["step"]): r for r in csv.DictReader(a.barriers.resolve().open(encoding="utf-8-sig", newline=""))}
        rows = []
        for i, rxn in enumerate(solver.elementary_rxns):
            sides = []
            for direction, species in (("forward", rxn[0]), ("reverse", rxn[-1])):
                factors = []
                for sp, count in Counter(species).items():
                    if sp in gas_p:
                        val = gas_p[sp]
                    elif sp in ads_map:
                        val = ads_map[sp]
                    else:
                        raise RuntimeError(f"unmapped mass-action species {sp}")
                    factors.append(val ** count)
                k = constants[i] if direction == "forward" else constants[n + i]
                sides.append((k * mp.fprod(factors), k))
            fwd, kf = sides[0]; rev, kr = sides[1]; net = fwd - rev
            native = native_nets[i]
            # Official directional API returns floats; compare after converting
            # those exact binary values to mp, but use high-precision values below.
            api_f, api_r = mp.mpf(str(official_dirs[i])), mp.mpf(str(official_dirs[n + i]))
            api_net = api_f - api_r
            rxn_text = model.print_rxn(rxn, mode="text")
            rows.append({
                "reaction_index": i + 1, "reaction_internal": rxn_text,
                "deltaG_eV": barrier_rows[i+1]["deltaG_eV"],
                "effective_forward_barrier_eV": barrier_rows[i+1]["forward_barrier_eV"],
                "effective_reverse_barrier_eV": barrier_rows[i+1]["reverse_barrier_eV"],
                "barrier_source": barrier_rows[i+1]["barrier_status"],
                "kf_catmap": mp.nstr(kf, 70), "kr_catmap": mp.nstr(kr, 70),
                "forward_flux_high_precision": mp.nstr(fwd, 70),
                "reverse_flux_high_precision": mp.nstr(rev, 70),
                "net_flux_high_precision": mp.nstr(net, 70),
                "catmap_native_net_rate": mp.nstr(native, 70),
                "native_net_relative_difference": mp.nstr(abs(net-native)/max(abs(net), abs(native), mp.mpf("1e-300")), 12),
                "official_directional_api_forward_float": mp.nstr(api_f, 25),
                "official_directional_api_reverse_float": mp.nstr(api_r, 25),
                "official_directional_api_subtracted_net_float": mp.nstr(api_net, 25),
                "theta_Ha": mp.nstr(cvmap["Ha_s"], 35), "theta_OOH": mp.nstr(cvmap["OOH_s"], 35),
                "theta_C6H12": mp.nstr(cvmap["C6H12_s"], 35), "theta_vacancy": mp.nstr(vacancy, 35),
                "theta_H2O2": mp.nstr(cvmap["H2O2_s"], 35),
                "theta_H2O": mp.nstr(cvmap["H2O_s"], 35),
            })
        for row in rows:
            row["theta_NEB1IS"] = mp.nstr(cvmap["NEB1IS_s"], 35)
            row["theta_NEB1FS"] = mp.nstr(cvmap["NEB1FS_s"], 35)
            row["theta_NEB2IS"] = mp.nstr(cvmap["NEB2IS_s"], 35)
            row["theta_NEB2FS"] = mp.nstr(cvmap["NEB2FS_s"], 35)
        fields = list(rows[0])
        with (out / "FORWARD_REVERSE_RATE_AUDIT.csv").open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
        tof_file = cov_path.parent / "candidate_tof.csv"
        with tof_file.open(encoding="utf-8-sig", newline="") as f:
            tof_rows = list(csv.DictReader(f))
        tof = next(mp.mpf(r["turnover_frequency"]) for r in tof_rows if r["gas_species"] == "C6H12O_g")
        step4 = rows[3]
        formula = constants[3] * cvmap["Ha_s"] * cvmap["OOH_s"] * cvmap["C6H12_s"]
        report = {
            "catmap_version": catmap.__version__, "precision_digits": a.precision,
            "source_model": str(model_path), "coverage_source": str(cov_path),
            "coverage_site_sum": mp.nstr(sum(ads) + vacancy, 70),
            "inventory_Ha_minus_OOH": mp.nstr(cvmap["Ha_s"] - cvmap["OOH_s"], 15),
            "step4_forward_formula": "kf4 * theta(Ha*) * theta(OOH*) * theta(C6H12*)",
            "step4_forward_formula_value": mp.nstr(formula, 70),
            "step4_forward_flux": step4["forward_flux_high_precision"],
            "product_TOF_reference": mp.nstr(tof, 70),
            "step4_forward_over_TOF": mp.nstr(formula/tof, 30),
            "official_directional_rates_precision_note": "CatMAP 0.3.1 get_directional_rates/get_elem_ec converts each directional flux to Python float; high-precision reconstruction from CatMAP rate constants and reaction stoichiometry is reported separately.",
            "native_net_rates_all_steps": [r["catmap_native_net_rate"] for r in rows],
        }
        (out / "rate_audit_metadata.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, ensure_ascii=False, indent=2))
    finally:
        os.chdir(old)


if __name__ == "__main__":
    main()
