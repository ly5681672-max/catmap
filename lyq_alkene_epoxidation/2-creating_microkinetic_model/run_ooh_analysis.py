# -*- coding: utf-8 -*-
"""CatMAP OOH-focused analysis for lyq_alkene_epoxidation.

Use this file in:
    H:\\catmap\\catmap\\lyq_alkene_epoxidation\\2-creating_microkinetic_model

Based on CatMAP official tutorials:
    https://catmap.readthedocs.io/en/latest/topics/thermodynamic_descriptors.html
    https://catmap.readthedocs.io/en/latest/topics/output_variables.html

This script never edits the existing alkene_epoxidation.mkm or test.ipynb.
Baseline TOF is exploratory until transition-state energetics are validated.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SETUP = ROOT / "alkene_epoxidation.mkm"
ENERGIES = ROOT / "energies_dft_neb.txt"
TS_CSV = ROOT / "DFT_TS_relative.csv"
OUTPUT = ROOT / "analysis_ooh"
SURFACES = [
    "tife", "timn", "tihf", "tire", "timo", "tiv", "tizr",
    "tico", "titi", "tiw", "tita", "ticr", "ti", "tini",
]
# CatMAP numbers_solver compares the SUM OF SQUARED residuals to tolerance.
# Old tolerance=1e-50 permits raw d(theta)/dt around 1e-25, making ~1e-26 TOFs unreliable.
# Ref: https://catmap.readthedocs.io/en/latest/tutorials/refining_a_microkinetic_model.html
# Source: https://github.com/SUNCAT-Center/catmap/blob/master/catmap/solvers/numbers_solver.py
DEFAULT_PRECISION = 180
DEFAULT_TOLERANCE = "1e-120"
DEFAULT_MIN_TOF = 0.0  # compatibility argument only; never use as a physical/numerical validity filter
MAX_RELATIVE_COVERAGE_RESIDUAL = 1e-4
TOF_STABILITY_LOG10_TOL = 0.02  # consecutive-stage TOFs agree within about 4.7% in relative terms
REFINEMENT_STAGES = (
    (180, "1e-120"),
    (260, "1e-180"),
    (360, "1e-260"),
    (460, "1e-340"),
)
MAX_RELATIVE_FLUX_ERROR = 1e-4
MAX_ABSOLUTE_COVERAGE_RESIDUAL = 1e-45
DIAGNOSTIC_SURFACES = ("tife", "titi", "tiw")

GAS_EFFECTIVE_ACTIVITIES = {
    "C6H12_g": 1.0,
    "H2O2_g": 0.185,
    "H2O_g": 0.815,
    "C6H12O_g": 1e-20,
}


def serialize_mp_vector(values, digits: int = 80) -> str:
    """Serialize an mpmath/vector-like sequence without float conversion."""
    from mpmath import mp

    return ";".join(mp.nstr(mp.mpf(value), n=digits) for value in values)


def read_data() -> pd.DataFrame:
    if not SETUP.exists() or not ENERGIES.exists():
        raise FileNotFoundError("把此脚本放在 .mkm 与 energies_dft_neb.txt 所在目录运行")
    table = pd.read_csv(ENERGIES, sep="\t")
    table["formation_energy"] = pd.to_numeric(table["formation_energy"], errors="raise")
    ts = pd.read_csv(TS_CSV).set_index("surface") if TS_CSV.exists() else None
    rows = []
    for surface in SURFACES:
        sp = table.loc[table["surface_name"].eq(surface)].set_index("species_name")
        missing = {"OOH", "OH", "O"} - set(sp.index)
        if missing:
            raise ValueError(f"{surface} 缺少关键DFT物种: {missing}")
        row = {
            "surface": surface,
            "G_OOH_eV": float(sp.at["OOH", "formation_energy"]),
            "G_OH_eV": float(sp.at["OH", "formation_energy"]),
            "G_O_eV": float(sp.at["O", "formation_energy"]),
        }
        if ts is not None and surface in ts.index:
            for key in (
                "dG_TS1_minus_IS1", "dG_TS1_minus_FS1",
                "dG_TS2_minus_IS2", "dG_TS2_minus_FS2",
            ):
                row[key] = float(ts.at[surface, key])
            row["TS_requires_check"] = any(row[key] < 0 for key in (
                "dG_TS1_minus_IS1", "dG_TS1_minus_FS1",
                "dG_TS2_minus_IS2", "dG_TS2_minus_FS2",
            ))
        rows.append(row)
    return pd.DataFrame(rows)


def fit_linear(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    a, b = np.polyfit(x, y, 1)
    pred = a * x + b
    denom = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - float(np.sum((y - pred) ** 2)) / denom if denom else float("nan")
    mae = float(np.mean(np.abs(y - pred)))
    return float(a), float(b), float(r2), mae


def remove_stale_plot(filename: str) -> None:
    """Remove *only* obsolete OOH-analysis plots; never touch figures_latest/."""
    for path in (OUTPUT / filename, (OUTPUT / filename).with_suffix(".pdf")):
        if path.exists():
            path.unlink()
            print(f"[已删除不再有效的历史图] {path}")


def plot_linear(df: pd.DataFrame, x_col: str, y_col: str, filename: str,
                xlabel: str, ylabel: str, log_tof: bool = False,
                min_points: int = 3) -> None:
    data = df[["surface", x_col, y_col]].copy()
    data[y_col] = pd.to_numeric(data[y_col], errors="coerce")
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if len(data) < min_points:
        remove_stale_plot(filename)
        print(f"[跳过] {filename}: 有效点数不足{min_points}，当前不能做拟合。")
        return
    x = data[x_col].to_numpy(dtype=float)
    y = data[y_col].to_numpy(dtype=float)
    if log_tof:
        good = y > 0
        x, y = x[good], np.log10(y[good])
        data = data.loc[good]
    if len(x) < min_points:
        remove_stale_plot(filename)
        print(f"[跳过] {filename}: 有效TOF点数不足{min_points}。")
        return
    a, b, r2, mae = fit_linear(x, y)
    xx = np.linspace(float(min(x)), float(max(x)), 150)
    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    ax.scatter(x, y, s=42)
    ax.plot(xx, a * xx + b, linestyle="--", linewidth=1.5)
    for _, row in data.iterrows():
        yy = np.log10(float(row[y_col])) if log_tof else float(row[y_col])
        ax.annotate(str(row["surface"]), (float(row[x_col]), yy),
                    xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"n={len(x)}  slope={a:.4f}  intercept={b:.4f}  R²={r2:.4f}  MAE={mae:.4f}")
    ax.grid(alpha=0.18)
    fig.tight_layout()
    path = OUTPUT / filename
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[已保存] {path}  n={len(x)} R²={r2:.4f} MAE={mae:.4f}")


def create_single_surface_setup(
    surface: str, work: Path, precision: int = DEFAULT_PRECISION,
    tolerance: str = DEFAULT_TOLERANCE, numbers_solver: bool = True,
) -> Path:
    """Build an isolated single-surface DFT model; retain the user's source .mkm."""
    source = SETUP.read_text(encoding="utf-8")
    extra = [
        "# One surface, original DFT energetics; no cross-surface scaling.",
        "scaler = 'ThermodynamicScaler'",
        f"surface_names = {[surface]!r}",
        "descriptor_names = ['temperature', 'pressure']",
        "descriptor_ranges = [[333.15, 333.15], [1.0, 1.0]]",
        "resolution = 1",
        "pressure_mode = 'concentration'",
        f"input_file = {str(ENERGIES)!r}",
        f"data_file = {'baseline_' + surface + '.pkl'!r}",
        f"decimal_precision = {precision}",
        # Use a true mpmath value: 1e-340 as a Python float becomes ZERO.
        "from mpmath import mp",
        f"tolerance = mp.mpf({tolerance!r})",
        "max_rootfinding_iterations = 250",
        f"use_numbers_solver = {numbers_solver}",
    ]
    for gas, value in GAS_EFFECTIVE_ACTIVITIES.items():
        extra.append(f"species_definitions[{gas!r}] = {{'concentration': {value!r}}}")
    path = work / f"baseline_{surface}.mkm"
    path.write_text(source + "\n\n" + "\n".join(extra) + "\n", encoding="utf-8")
    return path


def build_stoichiometric_matrices(elementary_rxns, gas_species=GAS_EFFECTIVE_ACTIVITIES):
    """Build gas/surface net-stoichiometry from CatMAP's parsed reactions.

    A CatMAP elementary reaction can contain transition-state/intermediate
    arrows. Its net balance is therefore taken from the first and last
    states, rather than inferred from a hard-coded step number. Bare ``s``
    and ``g`` site markers are omitted from the species-balance rows; the
    site balance is checked separately from the returned coverages.
    """
    from collections import Counter

    gas_names = tuple(gas_species)
    gas_set = set(gas_names)
    surface_names = []
    for reaction in elementary_rxns:
        if len(reaction) < 2:
            raise ValueError(f"Malformed CatMAP elementary reaction: {reaction!r}")
        for state in (reaction[0], reaction[-1]):
            for species in state:
                if species in gas_set or species in {"s", "g"}:
                    continue
                if species not in surface_names:
                    surface_names.append(species)

    def delta(reaction, species):
        initial = Counter(reaction[0])
        final = Counter(reaction[-1])
        return final.get(species, 0) - initial.get(species, 0)

    return {
        "gas_species": gas_names,
        "surface_species": tuple(surface_names),
        "gas_matrix": [
            [delta(reaction, species) for reaction in elementary_rxns]
            for species in gas_names
        ],
        "surface_matrix": [
            [delta(reaction, species) for reaction in elementary_rxns]
            for species in surface_names
        ],
    }


def check_steady_state(
    rates, gases, coverages, solver_residual, numbers_solver: bool,
    min_tof: float = DEFAULT_MIN_TOF, stoich_matrix: dict | None = None,
) -> dict:
    """Independent steady-state test in native mpmath precision.

    The CatMAP numbers solver compares an L2-squared surface residual with
    tolerance; the usual coverage solver uses an unsquared rate norm.
    Never treat a fixed TOF plot threshold as a scientific activity cutoff.
    """
    from mpmath import mp

    rr = [mp.mpf(str(v)) for v in rates]
    if len(rr) == 0:
        raise ValueError("At least one elementary rate is required")
    gas = {k: mp.mpf(str(v)) for k, v in gases.items()}
    keys = ("C6H12O_g", "C6H12_g", "H2O2_g", "H2O_g")
    if not all(k in gas for k in keys):
        raise ValueError("Gas turnover frequency lacks one or more net stoichiometric fluxes")

    rmax = max(abs(r) for r in rr)
    # Retain this legacy metric for this linear eight-step network. The
    # parsed stoichiometric balance below is the model-derived gate.
    cycle_reference = rr[3] if len(rr) > 3 else rr[0]
    cycle_abs = max(abs(r - cycle_reference) for r in rr)
    cycle_rel = cycle_abs / rmax if rmax else mp.inf

    net_terms = (gas["C6H12O_g"], -gas["C6H12_g"],
                 -gas["H2O2_g"], gas["H2O_g"])
    gas_max = max(abs(v) for v in net_terms)
    gas_abs = max(abs(v - net_terms[0]) for v in net_terms[1:])
    gas_rel = gas_abs / gas_max if gas_max else mp.inf

    residual = mp.mpf(str(solver_residual))
    if numbers_solver:
        residual = mp.sqrt(max(residual, mp.mpf("0")))
    residual_rel = residual / rmax if rmax else mp.inf

    cover = [mp.mpf(str(v)) for v in coverages]
    coverage_error = abs(sum(cover) - 1)
    nonnegative = all(v >= -mp.mpf("1e-25") for v in cover)
    net_tof = gas["C6H12O_g"]
    correct_signs = (net_tof > 0 and gas["C6H12_g"] < 0 and
                     gas["H2O2_g"] < 0 and gas["H2O_g"] > 0)
    stoich_surface_abs = mp.mpf("0")
    stoich_surface_rel = mp.mpf("0")
    gas_reconstruction_abs = mp.mpf("0")
    gas_reconstruction_rel = mp.mpf("0")
    if stoich_matrix is not None:
        if not stoich_matrix["gas_matrix"] or len(stoich_matrix["gas_matrix"][0]) != len(rr):
            raise ValueError("Stoichiometric matrix/rate vector length mismatch")
        gas_expected = [
            sum(mp.mpf(str(coeff)) * rate for coeff, rate in zip(row, rr))
            for row in stoich_matrix["gas_matrix"]
        ]
        gas_reconstruction_abs = max(
            abs(gas_expected[i] - gas[name])
            for i, name in enumerate(stoich_matrix["gas_species"])
        )
        gas_reconstruction_rel = (
            gas_reconstruction_abs / rmax if rmax else mp.inf
        )
        surface_expected = [
            sum(mp.mpf(str(coeff)) * rate for coeff, rate in zip(row, rr))
            for row in stoich_matrix["surface_matrix"]
        ]
        if surface_expected:
            stoich_surface_abs = max(abs(value) for value in surface_expected)
            stoich_surface_rel = (
                stoich_surface_abs / rmax if rmax else mp.inf
            )
    parsed_stoich_ok = (
        stoich_matrix is None
        or (
            stoich_surface_rel < mp.mpf(str(MAX_RELATIVE_FLUX_ERROR))
            and gas_reconstruction_rel < mp.mpf(str(MAX_RELATIVE_FLUX_ERROR))
        )
    )
    finite = all(mp.isfinite(v) for v in
                 (rmax, cycle_rel, gas_rel, residual, residual_rel, coverage_error))
    consistent = bool(
        finite and nonnegative and correct_signs
        and parsed_stoich_ok
        and coverage_error < mp.mpf("1e-8")
        and cycle_rel < mp.mpf(str(MAX_RELATIVE_FLUX_ERROR))
        and gas_rel < mp.mpf(str(MAX_RELATIVE_FLUX_ERROR))
        and residual_rel < mp.mpf(str(MAX_RELATIVE_COVERAGE_RESIDUAL))
        and residual < mp.mpf(str(MAX_ABSOLUTE_COVERAGE_RESIDUAL))
    )
    quality = ("validated_single_precision" if consistent else
               "failed_steady_state_checks")
    return {
        "cycle_max_relative_error": float(cycle_rel),
        "cycle_max_absolute_error": float(cycle_abs),
        "balance_max_relative_error": float(gas_rel),
        "balance_max_absolute_error": float(gas_abs),
        "steady_state_residual_max_abs": float(residual),
        "steady_state_residual_relative_to_max_rate": float(residual_rel),
        "coverage_sum_error": float(coverage_error),
        "cycle_max_rate_abs": float(rmax),
        "parsed_surface_balance_max_abs": float(stoich_surface_abs),
        "parsed_surface_balance_relative_to_max_rate": float(stoich_surface_rel),
        "parsed_gas_reconstruction_max_abs": float(gas_reconstruction_abs),
        "parsed_gas_reconstruction_relative_to_max_rate": float(gas_reconstruction_rel),
        "parsed_stoichiometry_checked": stoich_matrix is not None,
        "parsed_stoichiometry_pass": parsed_stoich_ok,
        "quality_status": quality,
        "quality_pass": consistent,
        "fit_rate_threshold": 0.0,  # deprecated; zero denotes no arbitrary TOF floor
        **{f"step_rate_{i + 1}": float(v) for i, v in enumerate(rr)},
    }


def compare_refinement_runs(previous: dict, current: dict,
                            max_log_delta: float = TOF_STABILITY_LOG10_TOL) -> dict:
    """Test stability across independently recomputed CatMAP precision/tolerance.

    A small absolute squared residual alone can misclassify ultra-low fluxes;
    two converged, stoichiometrically consistent roots with matching log(TOF)
    are required before an activity point can be used for regression.
    """
    import math
    old_log = float(previous.get("log10_net_C6H12O", float("nan")))
    new_log = float(current.get("log10_net_C6H12O", float("nan")))
    delta = abs(new_log - old_log)
    valid = bool(previous.get("quality_pass") and current.get("quality_pass"))
    stable = valid and math.isfinite(delta) and delta <= max_log_delta
    return {
        "refinement_stable": bool(stable),
        "refinement_delta_log10_TOF": float(delta),
    }




def run_single_surface(
    surface: str, precision: int = DEFAULT_PRECISION,
    tolerance: str = DEFAULT_TOLERANCE, numbers_solver: bool = True,
    output_root: Path = OUTPUT, min_tof: float = DEFAULT_MIN_TOF,
    seed_file: Path | None = None, max_iterations: int = 250,
    max_bisections: int = 3, capture_state: bool = False,
) -> dict:
    from catmap import ReactionModel

    work = output_root / surface
    work.mkdir(parents=True, exist_ok=True)
    config = create_single_surface_setup(
        surface, work, precision=precision, tolerance=tolerance,
        numbers_solver=numbers_solver
    )
    # Keep every attempt isolated. Only explicitly provided seed files may
    # populate CatMAP's numbers_map/coverage_map, as in the official tutorial.
    if max_iterations != 250 or max_bisections != 3:
        with config.open("a", encoding="utf-8") as fh:
            fh.write(f"\nmax_rootfinding_iterations = {int(max_iterations)}\n")
            fh.write(f"max_bisections = {int(max_bisections)}\n")
    for suffix in (".log", ".pkl"):
        cache = work / (config.stem + suffix)
        if cache.is_file():
            cache.unlink()
    if seed_file is not None:
        import shutil
        seed = Path(seed_file).resolve()
        if not seed.is_file():
            raise FileNotFoundError(f"CatMAP warm-start seed missing: {seed}")
        shutil.copy2(seed, work / (config.stem + ".pkl"))

    old_cwd = Path.cwd()
    try:
        os.chdir(work)
        model = ReactionModel(setup_file=config.name)
        model.output_variables = ["coverage", "rate", "turnover_frequency"]
        model.run(recalculate=(seed_file is not None))

        # CatMAP only creates output labels/maps after a point has been
        # accepted by the mapper.  When no point converges (the legacy
        # coverage solver commonly reports a singular Jacobian), indexing
        # output_labels directly would hide the numerical cause as a bare
        # KeyError.  Fail closed with the official CatMAP state exposed.
        labels = getattr(model, "output_labels", {}) or {}
        required = ("coverage", "rate", "turnover_frequency")
        missing = [name for name in required if name not in labels]
        if missing:
            available = ", ".join(sorted(labels)) or "<none>"
            map_state = ", ".join(
                f"{name}_map={'present' if getattr(model, name + '_map', None) else 'empty'}"
                for name in required
            )
            log_path = config.with_suffix(".log")
            raise RuntimeError(
                "CatMAP did not produce the requested output labels "
                f"{missing!r}; available={available}; {map_state}. "
                "This means the mapper did not accept a usable steady-state "
                f"point. Inspect the CatMAP log: {log_path}"
            )

        def target_values(variable):
            entries = getattr(model, variable + "_map", None)
            if not entries:
                raise RuntimeError(
                    f"CatMAP did not converge at 333.15 K: missing {variable}_map;"
                    f" inspect {config.with_suffix('.log')}"
                )
            matches = [
                values for coords, values in entries
                if abs(float(coords[0]) - 333.15) < 1e-5
                and abs(float(coords[1]) - 1.0) < 1e-7
            ]
            if not matches:
                raise RuntimeError(
                    f"CatMAP map contains no result at the required 333.15 K, 1.0 "
                    f"coordinate ({variable})."
                )
            return matches[-1]

        gas_labels = list(labels["turnover_frequency"])
        gas_raw = dict(zip(gas_labels, target_values("turnover_frequency")))
        cover_names = list(model.output_labels["coverage"])
        covers_raw = list(target_values("coverage"))
        rates_raw = list(target_values("rate"))
        numbers_raw = list(target_values("numbers")) if capture_state else None

        # CatMAP supports an additional last coverage for vacant sites
        # when the numbers-based solver is in use.
        if len(covers_raw) not in (len(cover_names), len(cover_names) + 1):
            raise ValueError("Unexpected coverage array/label length")
        covers = dict(zip(cover_names, covers_raw))
        vacant = (covers_raw[-1] if len(covers_raw) == len(cover_names) + 1
                  else 1 - sum(covers_raw))
        from mpmath import mp
        full_coverages = (
            list(covers_raw) if len(covers_raw) == len(cover_names) + 1
            else list(covers_raw) + [vacant]
        )

        # Independently recompute the CatMAP coverage residual.
        # Numbers solver returns a *squared L2 norm*; coverage solver
        # instead returns maximum absolute d(theta)/dt.
        # Coverage-based solver expects adsorbate coverages only; the
        # numbers solver uses adsorbates PLUS the vacant-site coverage.
        # Warm-start runs have TWO temperature grid points. Ensure the
        # independent residual is evaluated at the final target (333.15 K),
        # not at the most recently visited high-temperature seed point.
        model._descriptors = [333.15, 1.0]
        model.solver._descriptors = [333.15, 1.0]
        rxn_parameters = model.scaler.get_rxn_parameters([333.15, 1.0])
        residual_input = full_coverages if numbers_solver else list(covers_raw)
        residual = model.solver.get_residual(
            residual_input, validate_coverages=False, refresh_rate_constants=True
        )
        stoich_matrix = build_stoichiometric_matrices(
            model.elementary_rxns, gas_species=gas_raw.keys()
        )
        quality = check_steady_state(
            rates_raw, gas_raw, full_coverages, residual,
            numbers_solver=numbers_solver, min_tof=min_tof,
            stoich_matrix=stoich_matrix,
        )
        dominant = sorted(
            zip(cover_names + ["vacant"], [float(v) for v in full_coverages]),
            key=lambda x: x[1], reverse=True
        )[0]
        log_tof = (float(mp.log10(gas_raw["C6H12O_g"]))
                   if gas_raw["C6H12O_g"] > 0 else float("nan"))
        result = {
            "surface": surface,
            "status": "solved",  # solver accepted a root; NOT a physics quality gate
            "run_precision": precision,
            "run_tolerance": tolerance,
            "temperature_K": float(model.temperature),
            "pressure_mode": str(model.pressure_mode),
            "gas_thermo_mode": str(model.gas_thermo_mode),
            "adsorbate_thermo_mode": str(model.adsorbate_thermo_mode),
            "gas_effective_activities": str(GAS_EFFECTIVE_ACTIVITIES),
            "kinetic_prefactor_assumption": "CatMAP setup default (no explicit prefactor_list in base file)",
            "initial_seed_temperature_K": None,
            "initial_seed_source": str(seed_file) if seed_file else "",
            "solver_mode": "numbers" if numbers_solver else "coverages",
            "net_C6H12O": float(gas_raw["C6H12O_g"]),
            "net_C6H12O_high_precision": mp.nstr(gas_raw["C6H12O_g"], 40),
            "log10_net_C6H12O": log_tof,
            "net_C6H12": float(gas_raw["C6H12_g"]),
            "net_H2O2": float(gas_raw["H2O2_g"]),
            "net_H2O": float(gas_raw["H2O_g"]),
            "theta_OOH": float(covers.get("OOH_s", mp.nan)),
            "theta_OH": float(covers.get("OH_s", mp.nan)),
            "theta_O": float(covers.get("O_s", mp.nan)),
            "theta_vacant": float(vacant),
            "dominant_species": dominant[0],
            "dominant_coverage": dominant[1],
        }
        result.update(quality)
        if capture_state:
            # Recompute rates from the captured physical coverage on a fresh
            # solver state.  This is deliberately separate from CatMAP's map
            # extraction so stale rate constants cannot silently pass.
            for memo_name in ("_rate_constant_memoize", "_steady_state_memoize"):
                memo = getattr(model.solver, memo_name, None)
                if hasattr(memo, "clear"):
                    memo.clear()
            recomputed_rates = list(model.solver.get_rate(
                rxn_parameters,
                coverages=full_coverages,
                verify_coverages=False,
            ))
            recomputed_gas_values = list(model.solver.get_turnover_frequency(
                rxn_parameters,
                rates=recomputed_rates,
                verify_coverages=False,
            ))
            recomputed_gas = dict(zip(gas_labels, recomputed_gas_values))
            for memo_name in ("_rate_constant_memoize", "_steady_state_memoize"):
                memo = getattr(model.solver, memo_name, None)
                if hasattr(memo, "clear"):
                    memo.clear()
            recomputed_residual = model.solver.get_residual(
                full_coverages,
                validate_coverages=False,
                refresh_rate_constants=True,
            )
            recomputed_quality = check_steady_state(
                recomputed_rates,
                recomputed_gas,
                full_coverages,
                recomputed_residual,
                numbers_solver=numbers_solver,
                min_tof=min_tof,
                stoich_matrix=stoich_matrix,
            )
            from mpmath import mp
            rate_abs = max(
                abs(mp.mpf(a) - mp.mpf(b))
                for a, b in zip(rates_raw, recomputed_rates)
            )
            rate_scale = max(abs(mp.mpf(value)) for value in rates_raw)
            tof_abs = max(
                abs(mp.mpf(a) - mp.mpf(b))
                for a, b in zip(gas_raw.values(), recomputed_gas_values)
            )
            tof_scale = max(abs(mp.mpf(value)) for value in gas_raw.values())
            result.update({
                "coverage_order": ";".join(cover_names + ["vacant"]),
                "numbers_order": ";".join(
                    list(model.solver.adsorbate_names)
                    + [name for name in model.solver.site_names if name != "g"]
                ),
                "full_coverage_vector_high_precision": serialize_mp_vector(
                    full_coverages
                ),
                "full_numbers_vector_high_precision": serialize_mp_vector(
                    numbers_raw
                ),
                "net_rate_vector_high_precision": serialize_mp_vector(rates_raw),
                "gas_tof_vector_high_precision": serialize_mp_vector(
                    gas_raw.values()
                ),
                "recomputed_rate_vector_high_precision": serialize_mp_vector(
                    recomputed_rates
                ),
                "recomputed_gas_tof_vector_high_precision": serialize_mp_vector(
                    recomputed_gas_values
                ),
                "rate_recompute_max_abs": float(rate_abs),
                "rate_recompute_relative_error": float(
                    rate_abs / rate_scale if rate_scale else mp.inf
                ),
                "tof_recompute_max_abs": float(tof_abs),
                "tof_recompute_relative_error": float(
                    tof_abs / tof_scale if tof_scale else mp.inf
                ),
                "recompute_cache_cleared": True,
                "recomputed_quality_status": recomputed_quality["quality_status"],
                "recomputed_quality_pass": recomputed_quality["quality_pass"],
                "recomputed_steady_state_residual_relative_to_max_rate": (
                    recomputed_quality[
                        "steady_state_residual_relative_to_max_rate"
                    ]
                ),
                "recomputed_parsed_surface_balance_relative_to_max_rate": (
                    recomputed_quality[
                        "parsed_surface_balance_relative_to_max_rate"
                    ]
                ),
            })
        return result
    finally:
        os.chdir(old_cwd)


def run_adaptive_surface(surface: str, stages: int = 3) -> tuple[dict, list[dict]]:
    """Refine only where needed, retaining separate files for every attempt."""
    snapshots = []
    previous = None
    final = None
    for stage_idx, (prec, tol) in enumerate(REFINEMENT_STAGES[:stages], 1):
        label = f"stage_{stage_idx:02d}"
        try:
            row = run_single_surface(
                surface, precision=prec, tolerance=tol,
                output_root=OUTPUT / "refinement" / label
            )
        except Exception as exc:
            row = {
                "surface": surface, "status": f"FAILED: {type(exc).__name__}: {exc}",
                "quality_status": "solver_failed", "quality_pass": False,
            }
        row["refinement_stage"] = stage_idx
        row["refinement_stable"] = False
        row["refinement_delta_log10_TOF"] = float("nan")
        if previous is not None:
            row.update(compare_refinement_runs(previous, row))
        snapshots.append(row)
        final = dict(row)
        print(
            f"[{surface}] {label}: {row['quality_status']}"
            f" | TOF={row.get('net_C6H12O_high_precision', 'none')}"
            f" | relative_flux_error={row.get('cycle_max_relative_error', 'none')}"
            f" | relative_residual={row.get('steady_state_residual_relative_to_max_rate', 'none')}"
            f" | refinement_stable={row['refinement_stable']}"
        )
        if row["refinement_stable"]:
            final["quality_status"] = "validated_refinement_stable"
            final["quality_pass"] = True
            break
        previous = row
    if not final["refinement_stable"]:
        final["quality_pass"] = False
        final["quality_status"] = "unresolved_precision_or_flux"
    return final, snapshots


def run_diagnostics(min_tof: float, compare_coverage_solver: bool = False) -> None:
    """Official CatMAP tolerance/precision refinement, isolated for TiFe/TiTi/TiW.

    Legacy vs high-precision numbers solver are diagnostic comparisons.
    The optional coverage solver provides an independent convergence formulation.
    Output never overwrites dft_baseline.csv or the user's figures_latest/.
    """
    diag_root = OUTPUT / "diagnostics"
    settings = [
        ("legacy_numbers", 100, "1e-50", True),
        ("strict_numbers", DEFAULT_PRECISION, DEFAULT_TOLERANCE, True),
        ("refined_numbers", 260, "1e-180", True),
        ("deep_numbers", 360, "1e-260", True),
    ]
    if compare_coverage_solver:
        settings.append(("strict_coverages", 180, "1e-70", False))
    rows = []
    for tag, prec, tol, numbers in settings:
        for surface in DIAGNOSTIC_SURFACES:
            print(f"\n[诊断] {surface} | {tag} | digits={prec} tol={tol}")
            try:
                record = run_single_surface(
                    surface, precision=prec, tolerance=tol,
                    numbers_solver=numbers, output_root=diag_root / tag,
                    min_tof=min_tof
                )
            except Exception as exc:
                record = {
                    "surface": surface, "status": f"FAILED: {type(exc).__name__}: {exc}",
                    "quality_status": "solver_failed", "quality_pass": False,
                }
            record["diagnostic_stage"] = tag
            rows.append(record)
            print("[诊断结果]", record["status"], record["quality_status"],
                  "cycle_err=", record.get("cycle_max_relative_error"),
                  "dtheta=", record.get("steady_state_residual_max_abs"))
    result = pd.DataFrame(rows)
    target = diag_root / "diagnostics_comparison.csv"
    target.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(target, index=False)
    print(f"\n已导出诊断表：{target}")
    print("请先比较严格求解与旧设置的step_rate_1...8以及残差，再判断TOF。")


def main() -> None:
    parser = argparse.ArgumentParser(description="OOH-centered CatMAP DFT diagnostics")
    parser.add_argument("--recover", action="store_true",
                        help="用CatMAP官方MinResidMapper温度路径+缓存种子尝试恢复未验证的表面")
    parser.add_argument("--recover-surfaces", nargs="*", default=None,
                        help="指定恢复表面，如: --recover-surfaces titi timn ticr")
    parser.add_argument("--fit-only", action="store_true",
                        help="只读取已有CSV，不运行CatMAP；未验证的旧结果不会自动拟合")
    parser.add_argument("--diagnose", action="store_true",
                        help="专门对TiFe/TiTi/TiW比较旧求解与高精度求解")
    parser.add_argument("--coverage-crosscheck", action="store_true",
                        help="在--diagnose中额外对照传统coverage solver")
    # Precision/tolerance are selected as matched refinement stages to prevent
    # choosing an apparently converged residual that exceeds the actual TOF.
    parser.add_argument("--min-tof", type=float, default=DEFAULT_MIN_TOF,
                        help="旧命令兼容参数：不再作为TOF物理/数值有效性门槛")
    parser.add_argument("--max-stage", type=int, choices=(2, 3, 4), default=3,
                        help="独立求解精度级数：2/3/4，默认3")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    df = read_data()
    df.to_csv(OUTPUT / "dft_formation_descriptors.csv", index=False)
    plot_linear(df, "G_OOH_eV", "G_OH_eV", "OOH_OH_linear.png",
                "OOH* formation energy (eV)", "OH* formation energy (eV)")

    if args.recover:
        from recover_ooh_solver import recover_unvalidated_surfaces
        recover_unvalidated_surfaces(
            surfaces=args.recover_surfaces,
        )
        return

    if args.diagnose:
        run_diagnostics(args.min_tof, args.coverage_crosscheck)
        print("诊断模式不重写旧DFT基准表，也不生成未经验证的OOH-TOF拟合。")
        return

    if args.fit_only:
        target = OUTPUT / "dft_baseline.csv"
        if not target.exists():
            raise FileNotFoundError("没有dft_baseline.csv，请先正常运行脚本")
        baseline = pd.read_csv(target)
    else:
        baseline_rows = []
        stage_rows = []
        for surface in SURFACES:
            print(f"\n========== {surface}: strict DFT single-point baseline ==========")
            try:
                row, snapshots = run_adaptive_surface(surface, stages=args.max_stage)
                baseline_rows.append(row)
                stage_rows.extend(snapshots)
                print("[复算完成]", surface, row["quality_status"],
                      "cycle_err=", row.get("cycle_max_relative_error"))
            except Exception as exc:
                print(f"[求解失败] {surface}: {type(exc).__name__}: {exc}")
                baseline_rows.append({
                    "surface": surface,
                    "status": f"FAILED: {type(exc).__name__}: {exc}",
                    "quality_status": "solver_failed",
                    "quality_pass": False,
                })
        baseline = pd.DataFrame(baseline_rows)
        pd.DataFrame(stage_rows).to_csv(
            OUTPUT / "refinement_convergence.csv", index=False
        )
        baseline.to_csv(OUTPUT / "dft_baseline.csv", index=False)

    result = df.merge(baseline, on="surface", how="left")
    if "net_C6H12O" not in result.columns:
        result["net_C6H12O"] = np.nan

    # Fail closed: pre-v2 'solved' rows without independent numerical diagnostics
    # must NOT be used as scientifically validated TOF data.
    if "quality_pass" not in result.columns:
        result["quality_pass"] = False
        result["quality_status"] = "legacy_unverified"
    if "refinement_stable" not in result.columns:
        result["refinement_stable"] = False
    if "log10_net_C6H12O" not in result.columns:
        result["log10_net_C6H12O"] = np.nan
    result["quality_pass"] = result["quality_pass"].astype(str).str.lower().eq("true")
    result["refinement_stable"] = result["refinement_stable"].astype(str).str.lower().eq("true")
    # Do NOT reject any original DFT row based on the sign of an input
    # TS difference. CatMAP's effective barrier is non-negative by design;
    # the original DFT is kept unchanged and all TS metadata is descriptive.
    result["eligible_for_tof_fit"] = (
        result["quality_pass"] & result["refinement_stable"]
        & np.isfinite(pd.to_numeric(result.get("log10_net_C6H12O"), errors="coerce"))
    )
    result.to_csv(OUTPUT / "OOH_TOF_results.csv", index=False)

    selected = result.loc[result["eligible_for_tof_fit"]]
    print("\n通过两次独立精度收敛、八步通量、物料守恒的数据：", len(selected))
    print("催化剂：", ", ".join(selected["surface"]) or "无")
    plot_linear(selected, "G_OOH_eV", "log10_net_C6H12O",
                "OOH_logTOF_linear.png", "OOH* formation energy (eV)",
                "log10(net epoxide TOF)", min_points=5)
    print("\n输出目录：", OUTPUT)
    print("OOH-OH拟合直接使用原始DFT能量；OOH-TOF拟合要求两次独立求解的稳定性。")
    print("DFT能量及相对TS数据保持原样；不以负相对能量排除催化剂。")
    print("新的收敛轨迹在 analysis_ooh/refinement_convergence.csv。")


if __name__ == "__main__":
    main()
