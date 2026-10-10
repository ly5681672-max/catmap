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
# CatMAP 0.2.79: numbers_solver compares the SUM OF SQUARED residuals to tolerance.
# Old tolerance=1e-50 permits raw d(theta)/dt around 1e-25, making ~1e-26 TOFs unreliable.
# Ref: https://catmap.readthedocs.io/en/latest/tutorials/refining_a_microkinetic_model.html
# Source: https://github.com/SUNCAT-Center/catmap/blob/master/catmap/solvers/numbers_solver.py
DEFAULT_PRECISION = 180
DEFAULT_TOLERANCE = "1e-120"
DEFAULT_MIN_TOF = 1e-40  # diagnostic significance floor, not a physical TOF limit
MAX_RELATIVE_FLUX_ERROR = 1e-4
MAX_ABSOLUTE_COVERAGE_RESIDUAL = 1e-45
DIAGNOSTIC_SURFACES = ("tife", "titi", "tiw")

GAS_EFFECTIVE_ACTIVITIES = {
    "C6H12_g": 1.0,
    "H2O2_g": 0.185,
    "H2O_g": 0.815,
    "C6H12O_g": 1e-20,
}


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
                xlabel: str, ylabel: str, log_tof: bool = False) -> None:
    data = df[["surface", x_col, y_col]].copy()
    data[y_col] = pd.to_numeric(data[y_col], errors="coerce")
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if len(data) < 3:
        remove_stale_plot(filename)
        print(f"[跳过] {filename}: 有效点数不足3，当前不能做拟合。")
        return
    x = data[x_col].to_numpy(dtype=float)
    y = data[y_col].to_numpy(dtype=float)
    if log_tof:
        good = y > 0
        x, y = x[good], np.log10(y[good])
        data = data.loc[good]
    if len(x) < 3:
        remove_stale_plot(filename)
        print(f"[跳过] {filename}: 正TOF点数不足3。")
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
        f"tolerance = {tolerance}",
        "max_rootfinding_iterations = 250",
        f"use_numbers_solver = {numbers_solver}",
    ]
    for gas, value in GAS_EFFECTIVE_ACTIVITIES.items():
        extra.append(f"species_definitions[{gas!r}] = {{'concentration': {value!r}}}")
    path = work / f"baseline_{surface}.mkm"
    path.write_text(source + "\n\n" + "\n".join(extra) + "\n", encoding="utf-8")
    return path


def check_steady_state(
    rates, gases, coverages, solver_residual, numbers_solver: bool,
    min_tof: float = DEFAULT_MIN_TOF
) -> dict:
    """Independent diagnostics of eight-step flux and gas stoichiometry.

    Important: in CatMAP numbers_solver the convergence residual is L2**2.
    Convert to ordinary L2 using sqrt before comparing to physical rates.
    Preserve mpmath precision until writing summary floats.
    """
    from mpmath import mp

    rr = [mp.mpf(str(v)) for v in rates]
    if len(rr) != 8:
        raise ValueError(f"Expected 8 elementary rates, got {len(rr)}")
    gas = {k: mp.mpf(str(v)) for k, v in gases.items()}
    needed = ("C6H12O_g", "C6H12_g", "H2O2_g", "H2O_g")
    if any(k not in gas for k in needed):
        raise ValueError("Missing net gas flux labels, cannot check stoichiometry")

    rmax = max(abs(r) for r in rr)
    cycle_abs = max(abs(r - rr[3]) for r in rr)
    cycle_rel = cycle_abs / rmax if rmax else mp.nan

    gas_terms = (gas["C6H12O_g"], -gas["C6H12_g"],
                 -gas["H2O2_g"], gas["H2O_g"])
    gas_scale = max(abs(v) for v in gas_terms)
    gas_abs = max(abs(v - gas_terms[0]) for v in gas_terms[1:])
    gas_rel = gas_abs / gas_scale if gas_scale else mp.nan

    raw_resid = mp.mpf(str(solver_residual))
    if numbers_solver:
        raw_resid = mp.sqrt(max(raw_resid, mp.mpf("0")))

    cover = [mp.mpf(str(v)) for v in coverages]
    coverage_sum_err = abs(sum(cover) - 1)
    nonnegative = all(v >= -mp.mpf("1e-25") for v in cover)
    tof = gas["C6H12O_g"]
    signs_ok = (tof > 0 and gas["C6H12_g"] < 0 and
                gas["H2O2_g"] < 0 and gas["H2O_g"] > 0)

    finite = all(mp.isfinite(v) for v in
                 (rmax, cycle_rel, gas_rel, raw_resid, coverage_sum_err))
    numeric_ok = bool(
        finite and nonnegative and coverage_sum_err < mp.mpf("1e-8")
        and cycle_rel < mp.mpf(str(MAX_RELATIVE_FLUX_ERROR))
        and gas_rel < mp.mpf(str(MAX_RELATIVE_FLUX_ERROR))
        and raw_resid < mp.mpf(str(MAX_ABSOLUTE_COVERAGE_RESIDUAL))
        and signs_ok
    )
    significant = bool(tof >= mp.mpf(str(min_tof))) if mp.isfinite(tof) else False
    if not finite:
        quality = "invalid_nonfinite"
    elif not numeric_ok:
        quality = "failed_steady_state_checks"
    elif not significant:
        quality = "below_diagnostic_rate_floor"
    else:
        quality = "validated_numerically"
    return {
        "cycle_max_relative_error": float(cycle_rel),
        "cycle_max_absolute_error": float(cycle_abs),
        "balance_max_relative_error": float(gas_rel),
        "balance_max_absolute_error": float(gas_abs),
        "steady_state_residual_max_abs": float(raw_resid),
        "coverage_sum_error": float(coverage_sum_err),
        "cycle_max_rate_abs": float(rmax),
        "quality_status": quality,
        "quality_pass": quality == "validated_numerically",
        "fit_rate_threshold": min_tof,
        **{f"step_rate_{i + 1}": float(v) for i, v in enumerate(rr)}
    }


def run_single_surface(
    surface: str, precision: int = DEFAULT_PRECISION,
    tolerance: str = DEFAULT_TOLERANCE, numbers_solver: bool = True,
    output_root: Path = OUTPUT, min_tof: float = DEFAULT_MIN_TOF
) -> dict:
    from catmap import ReactionModel

    work = output_root / surface
    work.mkdir(parents=True, exist_ok=True)
    config = create_single_surface_setup(
        surface, work, precision=precision, tolerance=tolerance,
        numbers_solver=numbers_solver
    )
    # Do not carry over solutions from a previous, possibly invalid configuration.
    for suffix in (".log", ".pkl"):
        cache = work / (config.stem + suffix)
        if cache.is_file():
            cache.unlink()

    old_cwd = Path.cwd()
    try:
        os.chdir(work)
        model = ReactionModel(setup_file=config.name)
        model.output_variables = ["coverage", "rate", "turnover_frequency"]
        model.run()

        gas_labels = list(model.output_labels["turnover_frequency"])
        gas_raw = dict(zip(gas_labels, model.turnover_frequency_map[0][1]))
        cover_names = list(model.output_labels["coverage"])
        covers_raw = list(model.coverage_map[0][1])
        rates_raw = list(model.rate_map[0][1])

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
        residual = model.solver.get_residual(
            full_coverages, validate_coverages=False, refresh_rate_constants=True
        )
        quality = check_steady_state(
            rates_raw, gas_raw, full_coverages, residual,
            numbers_solver=numbers_solver, min_tof=min_tof
        )
        dominant = sorted(
            zip(cover_names + ["vacant"], [float(v) for v in full_coverages]),
            key=lambda x: x[1], reverse=True
        )[0]
        result = {
            "surface": surface,
            "status": "solved",  # solver accepted a root; NOT a physics quality gate
            "run_precision": precision,
            "run_tolerance": tolerance,
            "solver_mode": "numbers" if numbers_solver else "coverages",
            "net_C6H12O": float(gas_raw["C6H12O_g"]),
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
        return result
    finally:
        os.chdir(old_cwd)


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
    parser.add_argument("--fit-only", action="store_true",
                        help="只读取已有CSV，不运行CatMAP；未验证的旧结果不会自动拟合")
    parser.add_argument("--diagnose", action="store_true",
                        help="专门对TiFe/TiTi/TiW比较旧求解与高精度求解")
    parser.add_argument("--coverage-crosscheck", action="store_true",
                        help="在--diagnose中额外对照传统coverage solver")
    parser.add_argument("--precision", type=int, default=DEFAULT_PRECISION)
    parser.add_argument("--tolerance", default=DEFAULT_TOLERANCE,
                        help="numbers solver使用残差平方和；默认1e-120")
    parser.add_argument("--min-tof", type=float, default=DEFAULT_MIN_TOF,
                        help="进入拟合所需的最低可分辨TOF，仅为诊断数值门槛")
    args = parser.parse_args()
    if args.precision < 80:
        parser.error("precision须至少80位；建议180位")
    try:
        from decimal import Decimal
        tol = Decimal(args.tolerance)
        if tol <= 0 or tol < Decimal(10) ** (-args.precision):
            parser.error("tolerance必须大于0且不小于10^(-decimal_precision)")
    except (ValueError, ArithmeticError):
        parser.error("tolerance必须是正数科学计数法，如1e-120")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    df = read_data()
    df.to_csv(OUTPUT / "dft_formation_descriptors.csv", index=False)
    plot_linear(df, "G_OOH_eV", "G_OH_eV", "OOH_OH_linear.png",
                "OOH* formation energy (eV)", "OH* formation energy (eV)")

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
        for surface in SURFACES:
            print(f"\n========== {surface}: strict DFT single-point baseline ==========")
            try:
                row = run_single_surface(
                    surface, precision=args.precision,
                    tolerance=args.tolerance, min_tof=args.min_tof
                )
                baseline_rows.append(row)
                print("[计算完成]", surface, row["quality_status"],
                      "cycle_err=", row["cycle_max_relative_error"])
            except Exception as exc:
                print(f"[求解失败] {surface}: {type(exc).__name__}: {exc}")
                baseline_rows.append({
                    "surface": surface,
                    "status": f"FAILED: {type(exc).__name__}: {exc}",
                    "quality_status": "solver_failed",
                    "quality_pass": False,
                })
        baseline = pd.DataFrame(baseline_rows)
        baseline.to_csv(OUTPUT / "dft_baseline.csv", index=False)

    result = df.merge(baseline, on="surface", how="left")
    if "net_C6H12O" not in result.columns:
        result["net_C6H12O"] = np.nan

    # Fail closed: pre-v2 'solved' rows without independent numerical diagnostics
    # must NOT be used as scientifically validated TOF data.
    if "quality_pass" not in result.columns:
        result["quality_pass"] = False
        result["quality_status"] = "legacy_unverified"
    result["quality_pass"] = result["quality_pass"].astype(str).str.lower().eq("true")
    if "TS_requires_check" not in result.columns:
        result["TS_requires_check"] = True
    result["TS_requires_check"] = (
        result["TS_requires_check"].fillna(True).astype(str).str.lower().eq("true")
    )
    result["eligible_for_tof_fit"] = (
        result["quality_pass"] & ~result["TS_requires_check"]
        & (pd.to_numeric(result["net_C6H12O"], errors="coerce") >= args.min_tof)
    )
    result.to_csv(OUTPUT / "OOH_TOF_results.csv", index=False)

    selected = result.loc[result["eligible_for_tof_fit"]]
    print("\n通过稳态、化学计量、数值速率阈值且无已知TS异常的数据：", len(selected))
    print("催化剂：", ", ".join(selected["surface"]) or "无")
    plot_linear(selected, "G_OOH_eV", "net_C6H12O",
                "OOH_logTOF_linear.png", "OOH* formation energy (eV)",
                "log10(net epoxide TOF)", log_tof=True)
    print("\n输出目录：", OUTPUT)
    print("OOH-OH拟合是DFT能量关系；OOH-TOF拟合仅纳入明确通过质量门槛的行。")
    print("注意：TS异常标记仅基于现有NEB映射表，其他TS/溶剂效应问题仍需DFT复核。")


if __name__ == "__main__":
    main()