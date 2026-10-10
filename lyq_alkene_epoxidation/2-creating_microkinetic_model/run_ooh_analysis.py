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


def plot_linear(df: pd.DataFrame, x_col: str, y_col: str, filename: str,
                xlabel: str, ylabel: str, log_tof: bool = False) -> None:
    data = df[["surface", x_col, y_col]].copy()
    data[y_col] = pd.to_numeric(data[y_col], errors="coerce")
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if len(data) < 3:
        print(f"[跳过] {filename}: 有效点数不足3，当前不能做拟合。")
        return
    x = data[x_col].to_numpy(dtype=float)
    y = data[y_col].to_numpy(dtype=float)
    if log_tof:
        good = y > 0
        x, y = x[good], np.log10(y[good])
        data = data.loc[good]
    if len(x) < 3:
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


def create_single_surface_setup(surface: str, work: Path) -> Path:
    """Produce a separate .mkm from the unchanged 8-step setup."""
    source = SETUP.read_text(encoding="utf-8")
    extra = [
        "# Baseline overrides: one DFT surface; no linear scaling.",
        "scaler = 'ThermodynamicScaler'",
        f"surface_names = {[surface]!r}",
        "descriptor_names = ['temperature', 'pressure']",
        "descriptor_ranges = [[333.15, 333.15], [1.0, 1.0]]",
        "resolution = 1",
        "pressure_mode = 'concentration'",
        f"input_file = {str(ENERGIES)!r}",
        f"data_file = {'baseline_' + surface + '.pkl'!r}",
    ]
    for gas, value in GAS_EFFECTIVE_ACTIVITIES.items():
        extra.append(f"species_definitions[{gas!r}] = {{'concentration': {value!r}}}")
    path = work / f"baseline_{surface}.mkm"
    path.write_text(source + "\n\n" + "\n".join(extra) + "\n", encoding="utf-8")
    return path


def run_single_surface(surface: str) -> dict:
    from catmap import ReactionModel

    work = OUTPUT / surface
    work.mkdir(parents=True, exist_ok=True)
    config = create_single_surface_setup(surface, work)
    # Only remove our own analysis caches; no user inputs or figures are touched.
    for suffix in (".log", ".pkl"):
        candidate = work / (config.stem + suffix)
        if candidate.is_file():
            candidate.unlink()
    data_cache = work / f"baseline_{surface}.pkl"
    if data_cache.is_file():
        data_cache.unlink()

    old_cwd = Path.cwd()
    try:
        os.chdir(work)
        model = ReactionModel(setup_file=config.name)
        model.output_variables = ["coverage", "rate", "turnover_frequency"]
        model.run()
        gases = list(model.output_labels["turnover_frequency"])
        net = [float(v) for v in model.turnover_frequency_map[0][1]]
        gas_fluxes = dict(zip(gases, net))
        cover_names = list(model.output_labels["coverage"])
        cover_values = [float(v) for v in model.coverage_map[0][1]]
        covers = dict(zip(cover_names, cover_values))
        result = {
            "surface": surface,
            "status": "solved",
            "net_C6H12O": gas_fluxes.get("C6H12O_g", np.nan),
            "net_C6H12": gas_fluxes.get("C6H12_g", np.nan),
            "net_H2O2": gas_fluxes.get("H2O2_g", np.nan),
            "net_H2O": gas_fluxes.get("H2O_g", np.nan),
            "theta_OOH": covers.get("OOH_s", np.nan),
            "theta_OH": covers.get("OH_s", np.nan),
            "theta_O": covers.get("O_s", np.nan),
        }
        # OOH/OH label mismatches should never silently change columns.
        if len(cover_values) == len(cover_names) + 1:
            result["theta_vacant"] = cover_values[-1]
        else:
            result["theta_vacant"] = float(1 - sum(cover_values[:len(cover_names)]))
        terms = [
            result["net_C6H12O"], -result["net_C6H12"],
            -result["net_H2O2"], result["net_H2O"],
        ]
        scale = max(abs(t) for t in terms)
        if np.isfinite(scale) and scale > 0:
            result["balance_max_relative_error"] = max(
                abs(x - terms[0]) for x in terms[1:]
            ) / scale
        else:
            result["balance_max_relative_error"] = np.nan
        return result
    finally:
        os.chdir(old_cwd)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fit-only", action="store_true",
                        help="不运行CatMAP，直接用已有analysis_ooh/dft_baseline.csv作拟合")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    df = read_data()
    df.to_csv(OUTPUT / "dft_formation_descriptors.csv", index=False)
    plot_linear(df, "G_OOH_eV", "G_OH_eV", "OOH_OH_linear.png",
                "OOH* formation energy (eV)", "OH* formation energy (eV)")

    if args.fit_only:
        if not (OUTPUT / "dft_baseline.csv").exists():
            raise FileNotFoundError("未找到analysis_ooh/dft_baseline.csv；请先正常运行脚本")
        baseline = pd.read_csv(OUTPUT / "dft_baseline.csv")
    else:
        baseline_rows = []
        for surface in SURFACES:
            print(f"\n========== {surface}: single-surface DFT baseline ==========")
            try:
                baseline_rows.append(run_single_surface(surface))
            except Exception as e:
                print(f"[需排查] {surface}: {type(e).__name__}: {e}")
                baseline_rows.append({"surface": surface, "status": f"FAILED: {e}"})
        baseline = pd.DataFrame(baseline_rows)
        baseline.to_csv(OUTPUT / "dft_baseline.csv", index=False)

    result = df.merge(baseline, on="surface", how="left")
    if "net_C6H12O" not in result.columns:
        result["net_C6H12O"] = np.nan
    result.to_csv(OUTPUT / "OOH_TOF_results.csv", index=False)
    plot_linear(result, "G_OOH_eV", "net_C6H12O", "OOH_logTOF_linear.png",
                "OOH* formation energy (eV)", "log10(net epoxide TOF)", log_tof=True)
    print("\n输出目录：", OUTPUT)
    print("请先检查TS_requires_check以及balance_max_relative_error，再解释TOF拟合。")
    print("注意：当前负势垒/TS映射未经严格验证，速率结果只能视为诊断数据。")


if __name__ == "__main__":
    main()