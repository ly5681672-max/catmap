"""Plot existing descriptor data and CatMAP's built-in DRC outputs."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress, pearsonr, spearmanr


HERE = Path(__file__).resolve().parent
FOCUS = {"ti", "tico", "tini", "tire", "tife"}
COLORS = {
    "ti": "#0072B2", "tico": "#E69F00", "tini": "#009E73",
    "tire": "#CC79A7", "tife": "#D55E00", "other": "#858B93",
}
MARKERS = {"ti": "o", "tico": "s", "tini": "^", "tire": "D", "tife": "P"}
LABEL_OFFSETS = {
    "tife": (7, 9), "tire": (8, 10), "titi": (8, 0), "timo": (8, 10),
    "timn": (8, 10), "ticr": (8, 13), "tiw": (9, -14), "tiv": (-23, 8),
    "tita": (8, 7), "tihf": (8, 8),
}


def rows(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def main():
    data = [r for r in rows(HERE / "activity_analysis_data.csv") if r["status"] == "CONVERGED"]
    if len(data) != 13:
        raise RuntimeError(f"Expected 13 existing converged sites, found {len(data)}")
    x = np.array([float(r["step3_standard_delta_G_eV"]) for r in data])
    y = np.array([float(r["log10_TOF"]) for r in data])
    rho = spearmanr(x, y).statistic
    pearson = pearsonr(x, y).statistic
    fit = linregress(x, y)
    effective_shifts = np.array([
        float(r["step3_effective_delta_G_eV"]) - float(r["step3_standard_delta_G_eV"])
        for r in data
    ])
    shift = float(effective_shifts.mean())
    if float(effective_shifts.max() - effective_shifts.min()) > 1e-10:
        raise RuntimeError("Step-3 reservoir correction is not a common x-axis shift")

    plt.rcParams.update({
        "font.family": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(8.4, 6.0), layout="constrained")
    for r, xi, yi in zip(data, x, y):
        key = r["surface"]
        highlight = key in FOCUS
        ax.scatter(xi, yi, s=86 if highlight else 52,
                   color=COLORS.get(key, COLORS["other"]),
                   marker=MARKERS.get(key, "o"),
                   edgecolor="#222222" if highlight else "white",
                   linewidth=0.7, zorder=3 if highlight else 2)
        ax.annotate(r["label"], (xi, yi), xytext=LABEL_OFFSETS.get(key, (5, 5)),
                    textcoords="offset points", fontsize=8)
    xx = np.linspace(float(x.min()), float(x.max()), 200)
    ax.plot(xx, fit.slope * xx + fit.intercept, color="#444444",
            linestyle="--", linewidth=1.2, zorder=1,
            label="线性拟合（仅描述趋势）")
    ax.set_xlabel(r"$\Delta G_3$ standard formation-energy reference (eV)")
    ax.set_ylabel(r"$\log_{10}(\mathrm{TOF}_{C_6H_{12}O})$")
    ax.set_title("Step-3 free energy versus product TOF\n13 converged sites; Ti–Zr unresolved")
    ax.grid(color="#E1E4E8", linewidth=0.65)
    ax.set_axisbelow(True)
    ax.legend(loc="lower left", frameon=False)
    ax.text(0.02, 0.98,
            f"Spearman ρ = {rho:.3f}\nPearson r = {pearson:.3f}\nOLS R² = {fit.rvalue**2:.4f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#C8CDD2", "alpha": 0.95})
    sec = ax.secondary_xaxis("top", functions=(lambda v: v + shift, lambda v: v - shift))
    sec.set_xlabel(r"$\Delta G_3$ effective at fixed reservoir (eV)")
    ax.text(0.01, -0.16,
            f"Top axis: effective = standard {shift:+.6f} eV. Ti–Zr unresolved; no value imputed. n=13 exploratory sites.",
            transform=ax.transAxes, fontsize=8.5, color="#444444")
    fig.savefig(HERE / "descriptor_dg3_logtof.png", dpi=300, facecolor="white")
    fig.savefig(HERE / "descriptor_dg3_logtof.pdf", facecolor="white")
    plt.close(fig)

    drc_rows = rows(HERE / "official_drc" / "five_site_product_drc.csv")
    surfaces = ["ti", "tico", "tini", "tire", "tife"]
    parameters = ["Rpair_s", "TS1_s", "EpoxideHHO_s", "HHO_s", "TS2_s", "OHHa_s", "H2O_s"]
    drc = np.zeros((len(surfaces), len(parameters)))
    for i, surface in enumerate(surfaces):
        current = {r["energy_parameter"]: float(r["degree_of_rate_control"])
                   for r in drc_rows if r["surface"] == surface}
        if set(current) != set(parameters):
            raise RuntimeError(f"DRC species mismatch for {surface}: {sorted(current)}")
        drc[i] = [current[p] for p in parameters]

    fig, ax = plt.subplots(figsize=(10.0, 4.5), layout="constrained")
    lim = max(1.1, float(np.max(np.abs(drc))))
    image = ax.imshow(drc, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(parameters)), [
        "Rpair*", "TS1*\n(NEB1)", "EpoxideHHO*", "HHO*", "TS2*\n(NEB2)", "OHHa*", "H2O*"
    ])
    ax.set_yticks(range(len(surfaces)), ["Ti", "Ti–Co", "Ti–Ni", "Ti–Re", "Ti–Fe"])
    ax.set_title("CatMAP 0.3.1 built-in product-TOF degree of rate control\n"
                 r"$X_i=-k_BT\,\partial\ln(TOF)/\partial G_i$; fixed conditions, 333.15 K")
    for i in range(drc.shape[0]):
        for j in range(drc.shape[1]):
            val = drc[i, j]
            ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(val) > lim * 0.58 else "#202020")
    ax.set_xticks(np.arange(-0.5, len(parameters), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(surfaces), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False)
    cb = fig.colorbar(image, ax=ax, shrink=0.88, pad=0.025)
    cb.set_label("CatMAP DRC (dimensionless)")
    fig.savefig(HERE / "official_drc_product_tof.png", dpi=300, facecolor="white")
    fig.savefig(HERE / "official_drc_product_tof.pdf", facecolor="white")
    plt.close(fig)

    with (HERE / "descriptor_drc_statistics.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["n_sites", "spearman_rho", "pearson_r", "linear_R2",
                         "step3_effective_minus_standard_mean_eV",
                         "step3_shift_range_eV", "unresolved_sites"])
        writer.writerow([len(data), repr(rho), repr(pearson), repr(fit.rvalue**2),
                         repr(shift), repr(float(effective_shifts.max()-effective_shifts.min())), "tizr"])
    print(f"descriptor Spearman={rho:.6f}; Pearson={pearson:.6f}; R2={fit.rvalue**2:.6f}")
    print(f"step3 effective-standard common shift={shift:.12f} eV")
    print("Wrote descriptor_dg3_logtof.png/.pdf and official_drc_product_tof.png/.pdf")


if __name__ == "__main__":
    main()
