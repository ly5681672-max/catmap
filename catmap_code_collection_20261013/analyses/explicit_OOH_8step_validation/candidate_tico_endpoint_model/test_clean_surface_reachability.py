"""Try CatMAP 0.3.1's native ODE integrator from a clean surface state."""
from __future__ import annotations
import argparse, json, os
from pathlib import Path
import mpmath as mp

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--model", type=Path, required=True)
    args = p.parse_args()
    model_path = args.model.resolve()
    cwd = Path.cwd(); os.chdir(model_path.parent)
    out = model_path.parent / "clean_surface_ode_reachability.json"
    record = {"method": "CatMAP 0.3.1 SteadyStateSolver.get_ode_solution_coverage",
              "initial_state": "all adsorbate coverages zero; clean vacant surface",
              "status": "FAILED_OR_NOT_REACHED"}
    try:
        from catmap import ReactionModel
        model = ReactionModel(setup_file=model_path.name); model.run()
        solver = model.solver; solver.compile()
        mp.mp.dps = 100
        pars = list(model.scaler.get_rxn_parameters([mp.mpf("333.15"), mp.mpf("1.0")]))
        solver._descriptors = [mp.mpf("333.15"), mp.mpf("1.0")]
        c0 = [mp.mpf("0") for _ in solver.adsorbate_names]
        solver.get_rate_constants(pars, c0)
        # Keep official configured iteration cap and tolerance.
        cv = solver.get_ode_solution_coverage(pars, solver.ideal_steady_state_function, c0)
        vacancy = 1 - sum(cv)
        record.update({"status": "ODE_RETURNED", "coverages": {n: mp.nstr(v, 30)
                       for n, v in zip(solver.adsorbate_names, cv)},
                       "vacancy": mp.nstr(vacancy, 30),
                       "Ha_minus_OOH": mp.nstr(cv[solver.adsorbate_names.index("Ha_s")] -
                                                  cv[solver.adsorbate_names.index("OOH_s")], 12),
                       "note": "ODE return is not by itself proof of physically validated structural mapping"})
    except Exception as exc:
        record["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        out.write_text(json.dumps(record, indent=2), encoding="utf-8")
        os.chdir(cwd)
    print(json.dumps(record, ensure_ascii=False))

if __name__ == "__main__": main()
