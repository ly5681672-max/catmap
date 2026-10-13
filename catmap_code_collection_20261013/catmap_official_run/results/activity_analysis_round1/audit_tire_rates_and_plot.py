"""Audit saved Ti-Re directional rates and plot existing 13-site data only."""

import csv
from decimal import Decimal
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "official_strict_tolerance" / "tire"
LABELS = {
    "ti": "Ti", "tico": "Ti–Co", "tini": "Ti–Ni", "tire": "Ti–Re",
    "tife": "Ti–Fe", "titi": "Ti–Ti", "timo": "Ti–Mo", "timn": "Ti–Mn",
    "ticr": "Ti–Cr", "tiw": "Ti–W", "tiv": "Ti–V", "tita": "Ti–Ta",
    "tihf": "Ti–Hf",
}
LABEL_OFFSETS = {
    "tife": (7, 8), "tire": (8, 9), "titi": (8, -10), "timo": (8, 7),
    "timn": (8, 8), "ticr": (8, 11), "tiw": (8, -13), "tiv": (-24, 8),
    "tita": (8, 6), "tihf": (8, 7),
}


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def audit_directional_rates():
    with (RESULTS / "directional_reaction_rates.csv").open(
            encoding="utf-8-sig", newline="") as handle:
        raw_directional = list(csv.reader(handle))
    if not raw_directional or raw_directional[0] != ["direction", "reaction", "rate"]:
        raise RuntimeError("Unexpected original directional-rate CSV header")
    if any(len(row) != 2 for row in raw_directional[1:]):
        raise RuntimeError("Expected the known two-field CatMAP data rows")
    net = read_rows(RESULTS / "net_reaction_rates.csv")
    if len(raw_directional) != 13 or len(net) != 6:
        raise RuntimeError("Expected 12 directional and 6 net Ti-Re rate records")

    parsed = []
    for row in raw_directional[1:]:
        label = row[0]
        reaction, direction = label.rsplit(" ", 1)
        if direction not in {"forward", "reverse"}:
            raise RuntimeError(f"Unrecognized CatMAP directional label: {label}")
        parsed.append((reaction, direction, Decimal(row[1])))

    net_by_step = {int(row["step"]): row for row in net}
    out = []
    for step in range(1, 7):
        forward = parsed[step - 1]
        reverse = parsed[step + 5]
        net_row = net_by_step[step]
        if forward[0] != reverse[0] or forward[1] != "forward" or reverse[1] != "reverse":
            raise RuntimeError(f"Forward/reverse labels are not paired for step {step}")
        if forward[0] != net_row["reaction"]:
            raise RuntimeError(f"Reaction label mismatch at step {step}")
        difference = forward[2] - reverse[2]
        official_net = Decimal(net_row["net_rate"])
        delta = difference - official_net
        scale = max(abs(official_net), Decimal("1e-100"))
        matched = abs(delta) <= scale * Decimal("1e-12")
        out.append({
            "step": step,
            "reaction": forward[0],
            "forward_rate_from_existing_csv": str(forward[2]),
            "reverse_rate_from_existing_csv": str(reverse[2]),
            "forward_minus_reverse": str(difference),
            "CatMAP_net_rate_existing_output": str(official_net),
            "difference_directional_minus_net": str(delta),
            "relative_difference": str(delta / scale),
            "comparison": "MATCH" if matched else "MISMATCH",
            "same_exact_coverage_proven": "not recorded; CatMAP re-verifies coverage between output calculations",
        })

    out_path = HERE / "tire_directional_rate_audit.csv"
    with out_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out[0]))
        writer.writeheader()
        writer.writerows(out)
    return out


def plot_existing_data():
    data = [r for r in read_rows(HERE / "activity_analysis_data.csv")
            if r["status"] == "CONVERGED"]
    if len(data) != 13 or any(r["surface"] == "tizr" for r in data):
        raise RuntimeError("Expected exactly 13 converged sites, with Ti-Zr unresolved")

    plt.rcParams.update({
        "font.family": ["Microsoft YaHei", "Arial", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(9.2, 6.4), layout="constrained")
    xs = [float(r["step3_standard_delta_G_eV"]) for r in data]
    ys = [float(r["log10_TOF"]) for r in data]
    cvals = [float(r["step2_forward_barrier_eV"]) for r in data]
    sizes_raw = [float(r["step4_forward_barrier_eV"]) for r in data]
    low, high = min(sizes_raw), max(sizes_raw)
    sizes = [48 + 115 * (v - low) / (high - low) if high > low else 80 for v in sizes_raw]
    focus = {"ti", "tico", "tini", "tire", "tife"}
    image = None
    for r, x, y, color_value, size in zip(data, xs, ys, cvals, sizes):
        key = r["surface"]
        image = ax.scatter(
            x, y, c=[color_value], cmap="viridis", vmin=min(cvals), vmax=max(cvals),
            s=size, marker="o", edgecolor="#111111" if key in focus else "white",
            linewidth=0.9 if key in focus else 0.65, zorder=3 if key in focus else 2,
        )
        ax.annotate(LABELS[key], (x, y), xytext=LABEL_OFFSETS.get(key, (5, 5)),
                    textcoords="offset points", fontsize=8.5)

    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("NEB1 forward barrier, TS1 (eV)")
    size_values = [low, (low + high) / 2, high]
    size_handles = [Line2D([], [], marker="o", linestyle="", color="#777777",
                           markersize=(48 + 115 * (v - low) / (high - low)) ** 0.5,
                           label=f"{v:.2f} eV") for v in size_values]
    ax.legend(handles=size_handles, title="Marker area: NEB2 forward barrier",
              loc="lower left", frameon=True, fontsize=8, title_fontsize=8)
    ax.set_xlabel(r"Step-3 standard endpoint $\Delta G_3$ (eV)")
    ax.set_ylabel(r"$\log_{10}(\mathrm{TOF}_{C_6H_{12}O})$")
    ax.set_title("Existing CatMAP results: energy–activity relation\n"
                 "Color = NEB1 barrier; marker area = NEB2 barrier; Ti–Zr unresolved")
    ax.grid(color="#E1E4E8", linewidth=0.65)
    ax.set_axisbelow(True)
    ax.text(0.01, -0.15,
            "13 existing converged sites; no curve fitting/interpolation. "
            "NEB1/2 are explicit TS inputs; step 3 uses endpoint approximation.",
            transform=ax.transAxes, fontsize=8.5, color="#444444")
    fig.savefig(HERE / "joint_energy_activity_neb.png", dpi=300, facecolor="white")
    fig.savefig(HERE / "joint_energy_activity_neb.pdf", facecolor="white")
    plt.close(fig)


def main():
    audit = audit_directional_rates()
    plot_existing_data()
    for row in audit:
        if row["comparison"] == "MISMATCH":
            print("rate_mismatch step={} directional_minus_net={} relative={}".format(
                row["step"], row["difference_directional_minus_net"], row["relative_difference"]))
    print("Wrote Ti-Re aligned audit CSV and joint energy-activity figure from existing outputs only.")


if __name__ == "__main__":
    main()
