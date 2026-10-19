from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr


ARCHIVE = Path(r"H:\catmap\catmap\catmap_code_collection_20261013")
SOURCE = Path(r"H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model")
DESCRIPTOR_FILE = ARCHIVE / "analyses/BEP_energy_correction_optimization/data/source_carbon_O_H_TS1.csv"
TOF_FILE = ARCHIVE / "analyses/microkinetic_mechanism_13sites/data/gas_turnover_frequencies_13sites.csv"
OUT = SOURCE / "analysis_ooh/strict6step_ooh_association"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    descriptor_rows = read_csv(DESCRIPTOR_FILE)
    all_tof_rows = read_csv(TOF_FILE)
    product_rows = [row for row in all_tof_rows if row["gas_species"] == "C6H12O_g"]

    descriptors = {row["surface_name"]: row for row in descriptor_rows}
    tofs = {row["surface"]: row for row in product_rows}
    if len(descriptors) != len(descriptor_rows) or len(tofs) != len(product_rows):
        raise ValueError("surface 键存在重复，拒绝不确定连接")
    if set(descriptors) != set(tofs) or len(tofs) != 13:
        raise ValueError("描述符与产品 TOF 未形成相同的 13 位点集合")

    points: list[dict[str, object]] = []
    for surface in sorted(tofs):
        desc, result = descriptors[surface], tofs[surface]
        if desc["status"] != "CONVERGED" or float(desc["temperature_K"]) != 333.15:
            raise ValueError(f"{surface}: 六步严格 TOF 状态或温度不符合筛选条件")
        if result["sign"] != "1":
            raise ValueError(f"{surface}: C6H12O_g TOF 不是正向")
        tof = float(result["turnover_frequency"])
        if tof <= 0:
            raise ValueError(f"{surface}: 非正 TOF，不能取 log10")
        x = float(desc["OOH_star_legacy"])
        y = float(result["log10_abs_turnover_frequency"])
        if abs(y - math.log10(tof)) > 2e-9:
            raise ValueError(f"{surface}: TOF 与导出 log10 值不一致")
        points.append({
            "surface_name": surface,
            "catalyst_name": desc["catalyst"],
            "OOH_star_legacy_eV": x,
            "OOH_energy_reference": desc["OOH_star_reference"],
            "strict6step_C6H12O_TOF_s-1": tof,
            "log10_strict6step_C6H12O_TOF": y,
            "strict6step_status": desc["status"],
            "temperature_K": desc["temperature_K"],
            "OOH_source": str(DESCRIPTOR_FILE),
            "TOF_source": str(TOF_FILE),
        })

    xs = np.asarray([float(row["OOH_star_legacy_eV"]) for row in points])
    ys = np.asarray([float(row["log10_strict6step_C6H12O_TOF"]) for row in points])
    pearson = pearsonr(xs, ys)
    spearman = spearmanr(xs, ys)
    slope, intercept = np.polyfit(xs, ys, 1)
    fitted = slope * xs + intercept
    residuals = ys - fitted
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((ys - np.mean(ys))**2))
    r2 = 1.0 - ss_res / ss_tot
    rmse = float(np.sqrt(np.mean(residuals**2)))
    mae = float(np.mean(np.abs(residuals)))

    data_dir, figure_dir, script_dir = OUT / "data", OUT / "figures", OUT / "scripts"
    for path in (data_dir, figure_dir, script_dir):
        path.mkdir(parents=True, exist_ok=True)
    output_csv = data_dir / "OOH_legacy_Eform_vs_strict6step_logTOF_13sites.csv"
    png = figure_dir / "OOH_legacy_Eform_vs_strict6step_logTOF_13sites.png"
    pdf = figure_dir / "OOH_legacy_Eform_vs_strict6step_logTOF_13sites.pdf"
    notes = OUT / "OOH_legacy_Eform_vs_strict6step_logTOF_notes.md"
    # This directory was created specifically for this run; refresh only its
    # own generated artifacts when correcting figure layout.

    with output_csv.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(points[0]))
        writer.writeheader()
        writer.writerows(points)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "legend.fontsize": 8,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "axes.unicode_minus": False,
    })
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    other = [row for row in points if row["surface_name"] not in {"ti", "tico"}]
    ti = [row for row in points if row["surface_name"] == "ti"]
    tico = [row for row in points if row["surface_name"] == "tico"]
    ax.scatter([r["OOH_star_legacy_eV"] for r in other],
               [r["log10_strict6step_C6H12O_TOF"] for r in other],
               s=55, marker="o", color="#0072B2", edgecolor="white", linewidth=0.7,
               label="Other validated sites", zorder=3)
    ax.scatter([r["OOH_star_legacy_eV"] for r in ti],
               [r["log10_strict6step_C6H12O_TOF"] for r in ti],
               s=86, marker="D", color="#E69F00", edgecolor="black", linewidth=0.6,
               label="Ti", zorder=4)
    ax.scatter([r["OOH_star_legacy_eV"] for r in tico],
               [r["log10_strict6step_C6H12O_TOF"] for r in tico],
               s=95, marker="*", color="#D55E00", edgecolor="black", linewidth=0.5,
               label="Ti–Co", zorder=5)
    line_x = np.linspace(float(xs.min()), float(xs.max()), 200)
    ax.plot(line_x, slope * line_x + intercept, color="#555555", lw=1.5,
            linestyle="--", label="Descriptive linear fit", zorder=2)
    offsets = {
        "ti": (7, 7), "tico": (7, 7), "ticr": (-34, -13), "tife": (7, 8),
        "tihf": (7, 13), "timn": (7, -13), "timo": (-38, 8), "tini": (7, 7),
        "tire": (7, 7), "tita": (7, -13), "titi": (7, 7), "tiv": (7, -13),
        "tiw": (7, 7),
    }
    for row in points:
        ax.annotate(str(row["catalyst_name"]),
                    (row["OOH_star_legacy_eV"], row["log10_strict6step_C6H12O_TOF"]),
                    xytext=offsets[row["surface_name"]], textcoords="offset points",
                    fontsize=7.3, color="#202020",
                    arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.5},
                    zorder=6)
    ax.set_xlabel(r"Historical isolated-state $E_{form}(OOH^*)$ (eV)")
    ax.set_ylabel(r"Strict six-step CatMAP $\log_{10}(TOF_{C_6H_{12}O})$ (s$^{-1}$)")
    ax.set_title("OOH descriptor vs epoxide TOF — cross-model descriptive association")
    ax.text(0.02, 0.98,
            f"N=13  Pearson r={pearson.statistic:.3f} (p={pearson.pvalue:.3f})  "
            f"Spearman ρ={spearman.statistic:.3f} (p={spearman.pvalue:.3f})  "
            f"linear R²={r2:.3f}", transform=ax.transAxes, fontsize=8,
            va="top", bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#cccccc", "alpha": 0.92})
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    ranked = sorted(zip(points, residuals), key=lambda pair: abs(float(pair[1])), reverse=True)
    note_text = f"""# 历史 OOH* 形成能—正式六步 TOF 关联（探索性）

- N = {len(points)}；仅按 `surface_name` 一对一匹配，13/13 位点一致。
- 描述符：`OOH_star_legacy`，来源为 `{DESCRIPTOR_FILE}`，参考字段均为 `DFT_Gibbs_single`。
- 目标：正式六步 CatMAP 结果中 `C6H12O_g` 的 turnover frequency；输入表 13/13 状态为 CONVERGED、333.15 K、正向。TOF 与输入表内 log10 值最大差异小于 5e-16。
- Pearson r = {pearson.statistic:.6f}, p = {pearson.pvalue:.6f}; Spearman rho = {spearman.statistic:.6f}, p = {spearman.pvalue:.6f}.
- 样本内直线描述：R2 = {r2:.6f}, RMSE = {rmse:.6f} log10 units, MAE = {mae:.6f} log10 units；这只是简要描述统计，不是机理拟合或因果验证。
- 绝对线性残差较大的样本（仅作诊断线索）：{', '.join(str(x[0]['catalyst_name']) for x in ranked[:4])}。

## 解释边界

OOH* 是历史独立吸附态能量描述符；严格六步网络不含 OOH* 物种，因此该图不提供六步模型中的 OOH* 覆盖度。两边能量/TOF来自不同模型与参考数据路径，尚未证明slab与热修正完全统一。图中相关性只能用于探索，不能解释为八步 OOH* 稳态覆盖度导致 TOF 的因果证据。

## 文件

- 配对数据：`{output_csv}`
- 图：`{png}`、`{pdf}`
- 复现脚本：`{script_dir / Path(__file__).name}`
"""
    notes.write_text(note_text, encoding="utf-8")
    print(f"N={len(points)} Pearson r={pearson.statistic:.6f} p={pearson.pvalue:.6f}; "
          f"Spearman rho={spearman.statistic:.6f} p={spearman.pvalue:.6f}; "
          f"R2={r2:.6f} RMSE={rmse:.6f}")
    print("OUTPUT", output_csv)
    print("OUTPUT", png)
    print("OUTPUT", pdf)
    print("OUTPUT", notes)


if __name__ == "__main__":
    main()
