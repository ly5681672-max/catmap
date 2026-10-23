"""Replay two Ti-Co full-NEB single steps with CatMAP 0.3.1 only.

This deliberately initializes/compiles the official model without invoking
mapper, steady-state solve, coverage, TOF, or DRC calculations.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
CASES = {
    "NEB1": {
        "model": ROOT / "NEB1_single_step.mkm",
        "input": ROOT / "NEB1_energies.tsv",
        "expected_dG": mp.mpf("-0.239686"),
        "expected_fwd": mp.mpf("0.556344"),
        "expected_rev": mp.mpf("0.796030"),
    },
    "NEB2": {
        "model": ROOT / "NEB2_single_step.mkm",
        "input": ROOT / "NEB2_energies.tsv",
        "expected_dG": mp.mpf("-0.480580"),
        "expected_fwd": mp.mpf("0.225994"),
        "expected_rev": mp.mpf("0.706574"),
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def s(value) -> str:
    return mp.nstr(value, 45)


def run_case(name: str, cfg: dict) -> dict:
    import catmap
    from catmap import ReactionModel

    previous = Path.cwd()
    os.chdir(ROOT)
    try:
        model = ReactionModel(setup_file=cfg["model"].name)
        model.verify()
        if len(model.elementary_rxns) != 1:
            raise RuntimeError(f"{name}: expected exactly one CatMAP elementary reaction")
        rxn = model.elementary_rxns[0]
        if len(rxn) != 3:
            raise RuntimeError(f"{name}: explicit IS-TS-FS reaction not parsed: {rxn}")

        # Run only the model's initialization path with no descriptor values or
        # ranges; this installs CatMAP's official numerical functions but does
        # not run descriptor mapping or any steady-state calculation.
        model.output_variables = []
        model.descriptor_ranges = None
        model.descriptors = None
        model.descriptor_values = None
        model.run()
        solver = model.solver
        solver.compile()
        mp.mp.dps = 80
        descriptors = [mp.mpf("333.15")]
        rxn_params = list(model.scaler.get_rxn_parameters(descriptors))
        energy_dict = {
            key: val for key, val in zip(
                list(model.adsorbate_names) + list(model.transition_state_names),
                rxn_params,
            )
        }
        energy_dict["s"] = mp.mpf("0")
        delta_g, cat_fwd = model.get_rxn_energy(rxn, energy_dict)
        _, cat_rev = model.get_rxn_energy(model.reverse_rxn(rxn), energy_dict)

        # Call the official compiled CatMAP rate-constant function through
        # SteadyStateSolver.get_rate_constants. A unit formal coverage is only
        # needed by the API; ideal single-step rate constants do not depend on it.
        solver.get_rate_constants(
            list(rxn_params),
            [mp.mpf("1") if i == 0 else mp.mpf("0")
             for i in range(len(solver.adsorbate_names))],
        )
        kf, kr = solver._kf[0], solver._kr[0]
        kT = model._kB * mp.mpf("333.15")
        expected_ratio = mp.exp(-delta_g / kT)
        ratio = kf / kr

        close = lambda x, y: abs(x-y) <= mp.mpf("1e-35")
        checks = {
            "catmap_031": catmap.__version__ == "0.3.1",
            "one_reaction": len(model.elementary_rxns) == 1,
            "explicit_ts_parsed": len(rxn) == 3,
            "mass_and_site_balance": True,  # model.verify() returned normally
            "deltaG_matches": close(delta_g, cfg["expected_dG"]),
            "forward_barrier_matches": close(cat_fwd, cfg["expected_fwd"]),
            "reverse_barrier_matches": close(cat_rev, cfg["expected_rev"]),
            # CatMAP's generated function embeds its constants as Python
            # literals; comparison to the high-precision scalar path is thus
            # checked at a conservative near-machine precision.
            "detailed_balance_rate_ratio": abs(ratio / expected_ratio - 1) <= mp.mpf("1e-12"),
        }
        return {
            "step": name,
            "catmap_version": catmap.__version__,
            "reaction_expression": str(model.rxn_expressions[0]),
            "parsed_reaction": repr(rxn),
            "adsorbates": list(model.adsorbate_names),
            "transition_states": list(model.transition_state_names),
            "state_compositions": {
                key: model.species_definitions[key].get("composition")
                for key in model.adsorbate_names + model.transition_state_names
            },
            "formal_n_sites_bookkeeping_only": {
                key: model.species_definitions[key].get("n_sites", 1)
                for key in model.adsorbate_names + model.transition_state_names
            },
            "deltaG_eV": s(delta_g),
            "forward_barrier_CatMAP_eV": s(cat_fwd),
            "reverse_barrier_CatMAP_eV": s(cat_rev),
            "forward_rate_constant_s-1": s(kf),
            "reverse_rate_constant_s-1": s(kr),
            "rate_constant_ratio_kf_over_kr": s(ratio),
            "detailed_balance_exp_minus_dG_over_kT": s(expected_ratio),
            "detailed_balance_relative_error": s(abs(ratio / expected_ratio - 1)),
            "checks": checks,
            "model_sha256": sha256(cfg["model"]),
            "energy_input_sha256": sha256(cfg["input"]),
            "scope": "single-step barriers/rate constants only; no steady state, coverage, TOF, or DRC",
        }
    finally:
        os.chdir(previous)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    RESULTS.mkdir(exist_ok=True)
    rows = []
    for name, cfg in CASES.items():
        rows.append(run_case(name, cfg))
    payload = {"status": "PASS" if all(all(r["checks"].values()) for r in rows) else "FAIL",
               "results": rows}
    (RESULTS / "native_replay_results.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    fields = ["step", "catmap_version", "reaction_expression", "deltaG_eV",
              "forward_barrier_CatMAP_eV", "reverse_barrier_CatMAP_eV",
              "forward_rate_constant_s-1", "reverse_rate_constant_s-1",
              "rate_constant_ratio_kf_over_kr", "detailed_balance_exp_minus_dG_over_kT",
              "detailed_balance_relative_error", "checks"]
    with (RESULTS / "native_replay_results.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fields if key != "checks"} |
                            {"checks": json.dumps(row["checks"], ensure_ascii=False)})
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if payload["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
