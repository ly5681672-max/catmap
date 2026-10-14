"""Summarize and plot existing official CatMAP outputs; does not run CatMAP."""

import csv
import math
from decimal import Decimal
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
RESULTS_ROOT = HERE.parent
ANALYSIS_DATA = PROJECT_ROOT / "analyses" / "microkinetic_mechanism_13sites" / "data"
ENERGY_FILE = RESULTS_ROOT / "official_strict_tolerance" / "tire" / "energies.txt"
OUT = HERE
PLOT_OUT = HERE / "rerun_from_archived_tables"
SURFACES = [
    "tife", "timn", "tihf", "tire", "timo", "tiv", "tizr", "tico",
    "titi", "tiw", "tita", "ticr", "ti", "tini",
]
LABELS = {
    "tife": "Ti–Fe", "timn": "Ti–Mn", "tihf": "Ti–Hf", "tire": "Ti–Re",
    "timo": "Ti–Mo", "tiv": "Ti–V", "tizr": "Ti–Zr", "tico": "Ti–Co",
    "titi": "Ti–Ti", "tiw": "Ti–W", "tita": "Ti–Ta", "ticr": "Ti–Cr",
    "ti": "Ti", "tini": "Ti–Ni",
}
FOCUS = {"ti", "tico", "tini", "tire", "tife"}
PALETTE = {
    "ti": "#0072B2", "tico": "#E69F00", "tini": "#009E73",
    "tire": "#CC79A7", "tife": "#D55E00",
}
MARKERS = {"ti": "o", "tico": "s", "tini": "^", "tire": "D", "tife": "P"}


def read_csv(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def index_surface_rows(rows, source, expected_surfaces):
    indexed = {}
    for row in rows:
        surface = row["surface"]
        if surface in indexed:
            raise RuntimeError(f"Duplicate surface {surface!r} in {source}")
        indexed[surface] = row
    if set(indexed) != expected_surfaces:
        missing = sorted(expected_surfaces - set(indexed))
        extra = sorted(set(indexed) - expected_surfaces)
        raise RuntimeError(f"Surface mismatch in {source}: missing={missing}, extra={extra}")
    return indexed


def index_step_rows(rows, source, expected_surfaces):
    indexed = {}
    for row in rows:
        key = (row["surface"], int(row["step"]))
        if key in indexed:
            raise RuntimeError(f"Duplicate surface/step {key!r} in {source}")
        indexed[key] = row
    expected = {(surface, step) for surface in expected_surfaces for step in range(1, 7)}
    if set(indexed) != expected:
        missing = sorted(expected - set(indexed))
        extra = sorted(set(indexed) - expected)
        raise RuntimeError(f"Surface/step mismatch in {source}: missing={missing[:12]}, extra={extra[:12]}")
    return indexed


def main():
    mp.mp.dps = 75
    expected_valid = set(SURFACES) - {"tizr"}
    summary_file = ANALYSIS_DATA / "convergence_13sites.csv"
    coverage_file = ANALYSIS_DATA / "coverage_13sites.csv"
    rate_file = ANALYSIS_DATA / "directional_reaction_rates_13sites.csv"
    barrier_file = ANALYSIS_DATA / "effective_barriers_13sites.csv"
    required = (summary_file, coverage_file, rate_file, barrier_file, ENERGY_FILE)
    missing_files = [str(path) for path in required if not path.is_file()]
    if missing_files:
        raise FileNotFoundError("Required archived CatMAP inputs are missing:\n" + "\n".join(missing_files))

    summary = index_surface_rows(read_csv(summary_file), summary_file.name, expected_valid)
    if any(row["status"] != "CONVERGED" for row in summary.values()):
        raise RuntimeError("The archived 13-site convergence table contains a non-CONVERGED site")
    # Ti-Zr has no validated numerical result: retain it explicitly without assigning TOF=0.
    summary["tizr"] = {"surface": "tizr", "status": "UNRESOLVED", "tof_C6H12O": ""}
    coverage_by_surface = index_surface_rows(read_csv(coverage_file), coverage_file.name, expected_valid)
    rate_by_surface_step = index_step_rows(read_csv(rate_file), rate_file.name, expected_valid)
    barrier_by_surface_step = index_step_rows(read_csv(barrier_file), barrier_file.name, expected_valid)
    energy_rows = read_csv(ENERGY_FILE, delimiter="\t")
    energies = {}
    gas_energies = {}
    for row in energy_rows:
        value = Decimal(row["formation_energy"])
        if row["surface_name"] == "None":
            gas_energies[row["species_name"]] = value
        else:
            energies[(row["surface_name"], row["species_name"])] = value

    data = []
    for surface in SURFACES:
        s = summary[surface]
        row = {
            "surface": surface,
            "label": LABELS[surface],
            "status": s["status"],
            "tof_C6H12O_g": s["tof_C6H12O"],
            "coverage_vacancy": "",
            "coverage_EpoxideHHO": "",
            "log10_TOF": "",
            "log10_vacancy_coverage": "",
            "log10_EpoxideHHO_coverage": "",
        }
        if s["status"] != "CONVERGED":
            data.append(row)
            continue

        coverages = coverage_by_surface[surface]
        e = lambda species: energies[(surface, species)]
        standard_dg = {
            2: e("EpoxideHHO") - e("Rpair"),
            3: e("HHO") + gas_energies["C6H12O"] - e("EpoxideHHO"),
            4: e("OHHa") - e("HHO"),
        }
        row.update(
            coverage_vacancy=coverages["s"],
            coverage_EpoxideHHO=coverages["EpoxideHHO_s"],
            log10_TOF=mp.nstr(mp.log10(mp.mpf(row["tof_C6H12O_g"])), 18),
            log10_vacancy_coverage=mp.nstr(mp.log10(mp.mpf(coverages["s"])), 18),
            log10_EpoxideHHO_coverage=mp.nstr(mp.log10(mp.mpf(coverages["EpoxideHHO_s"])), 18),
        )
        for step in (2, 3, 4):
            barrier = barrier_by_surface_step[(surface, step)]
            rate = rate_by_surface_step[(surface, step)]
            row[f"step{step}_standard_delta_G_eV"] = str(standard_dg[step])
            row[f"step{step}_effective_delta_G_eV"] = barrier["delta_G_eV"]
            row[f"step{step}_forward_barrier_eV"] = barrier["forward_effective_barrier_eV"]
            row[f"step{step}_reverse_barrier_eV"] = barrier["reverse_effective_barrier_eV"]
            # Keep CatMAP's directly exported high-precision net rate; do not reconstruct it
            # from lower-precision forward/reverse columns.
            row[f"step{step}_net_rate"] = rate["net_rate_high_precision"]
        data.append(row)

    fields = [
        "surface", "label", "status", "tof_C6H12O_g", "log10_TOF",
        "coverage_vacancy", "log10_vacancy_coverage", "coverage_EpoxideHHO",
        "log10_EpoxideHHO_coverage",
    ] + [
        f"step{step}_{metric}"
        for step in (2, 3, 4)
        for metric in (
            "standard_delta_G_eV", "effective_delta_G_eV", "forward_barrier_eV",
            "reverse_barrier_eV", "net_rate",
        )
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    PLOT_OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "activity_analysis_data.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)

    valid = [r for r in data if r["status"] == "CONVERGED"]
    ordered = sorted(valid, key=lambda r: mp.mpf(r["tof_C6H12O_g"]), reverse=True)
    plt.rcParams.update({"font.family": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42})

    fig, ax = plt.subplots(figsize=(8.2, 8), layout="constrained")
    y = list(range(len(ordered)))
    for yi, row in zip(y, ordered):
        key = row["surface"]
        highlight = key in FOCUS
        ax.scatter(float(row["log10_TOF"]), yi, s=78 if highlight else 43,
                   marker=MARKERS.get(key, "o"),
                   color=PALETTE.get(key, "#8A8F98"),
                   edgecolor="#222222" if highlight else "white", linewidth=0.7,
                   zorder=3 if highlight else 2)
    ax.set_yticks(y, [r["label"] for r in ordered])
    ax.invert_yaxis()
    ax.set_xlabel(r"$\log_{10}(\mathrm{TOF}_{C_6H_{12}O})$ (CatMAP output; no zero imputation)")
    ax.set_title("333.15 K exploratory activity ranking\n13 converged sites; Ti–Zr unresolved")
    ax.grid(axis="x", color="#D9DDE2", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.text(0.01, -0.105,
            "Ti–Zr: unresolved, excluded from numerical ranking (not assigned TOF = 0).",
            transform=ax.transAxes, fontsize=9, color="#444444")
    fig.savefig(PLOT_OUT / "tof_rank_log10.png", dpi=300, facecolor="white")
    fig.savefig(PLOT_OUT / "tof_rank_log10.pdf", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.7), sharey=True, layout="constrained")
    xfields = ["log10_vacancy_coverage", "log10_EpoxideHHO_coverage"]
    xlabels = [r"$\log_{10}(\theta_{vacancy})$", r"$\log_{10}(\theta_{EpoxideHHO})$"]
    for ax, xfield, xlabel in zip(axes, xfields, xlabels):
        for row in valid:
            key = row["surface"]
            highlight = key in FOCUS
            ax.scatter(float(row[xfield]), float(row["log10_TOF"]),
                       s=78 if highlight else 43,
                       marker=MARKERS.get(key, "o"),
                       color=PALETTE.get(key, "#8A8F98"),
                       edgecolor="#222222" if highlight else "white", linewidth=0.7,
                       zorder=3 if highlight else 2)
            if highlight:
                ax.annotate(LABELS[key], (float(row[xfield]), float(row["log10_TOF"])),
                            xytext=(5, 5), textcoords="offset points", fontsize=8)
        ax.set_xlabel(xlabel + " (fractional coverage)")
        ax.grid(color="#E1E4E8", linewidth=0.65)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(r"$\log_{10}(\mathrm{TOF}_{C_6H_{12}O})$")
    fig.suptitle("Coverage–activity comparison at 333.15 K\nFocus sites labeled; Ti–Zr omitted as unresolved")
    fig.savefig(PLOT_OUT / "coverage_vs_tof.png", dpi=300, facecolor="white")
    fig.savefig(PLOT_OUT / "coverage_vs_tof.pdf", facecolor="white")
    plt.close(fig)

    print(f"Wrote {OUT / 'activity_analysis_data.csv'}")
    print(f"Wrote {PLOT_OUT / 'tof_rank_log10.png'} and PDF")
    print(f"Wrote {PLOT_OUT / 'coverage_vs_tof.png'} and PDF")
    print(f"Plotted {len(valid)} converged sites; unresolved sites remain missing.")


if __name__ == "__main__":
    main()
