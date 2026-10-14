"""Run the official CatMAP 0.3.1 mapper/scaler/solver on a TS1-TS2 grid."""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from catmap import ReactionModel, __version__
from catmap.analyze import VectorMap


CATALYST_LABELS = {
    "ti": "Ti", "tico": "Ti-Co", "tini": "Ti-Ni", "tife": "Ti-Fe",
    "tire": "Ti-Re", "timo": "Ti-Mo", "timn": "Ti-Mn", "ticr": "Ti-Cr",
    "tiw": "Ti-W", "tiv": "Ti-V", "tita": "Ti-Ta", "tihf": "Ti-Hf", "titi": "Ti-Ti",
}


def load_validated_points():
    energy = {}
    with Path("energies.txt").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["surface_name"] in CATALYST_LABELS and row["species_name"] in {"TS1", "TS2"}:
                key = (row["surface_name"], row["species_name"])
                if key in energy:
                    raise ValueError(f"Duplicate descriptor energy: {key}")
                energy[key] = float(row["formation_energy"])
    points = []
    with Path("../data/matched_data.csv").open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            if row["surface_name"] not in CATALYST_LABELS:
                continue
            if row["status"] != "CONVERGED" or float(row["temperature_K"]) != 333.15:
                raise ValueError(f"Unvalidated point: {row}")
            points.append({
                "surface_name": row["surface_name"], "catalyst": CATALYST_LABELS[row["surface_name"]],
                "TS1_s_eV": energy[(row["surface_name"], "TS1")],
                "TS2_s_eV": energy[(row["surface_name"], "TS2")],
                "actual_TOF_C6H12O": float(row["TOF"]),
                "actual_log10_TOF_C6H12O": float(row["log10_TOF"]),
            })
    if len(points) != 13 or {p["surface_name"] for p in points} != set(CATALYST_LABELS):
        raise ValueError(f"Expected 13 unique validated catalysts, got {len(points)}")
    return points


def fmt(value):
    return str(value)


def write_grid(nodes, output_dir):
    xvals = sorted({float(row["TS1_s_eV"]) for row in nodes})
    yvals = sorted({float(row["TS2_s_eV"]) for row in nodes})
    if len(nodes) != len(xvals) * len(yvals):
        raise RuntimeError("CatMAP nodes do not form a complete Cartesian descriptor grid")
    z = np.full((len(xvals), len(yvals)), np.nan)
    indexes_x = {v: i for i, v in enumerate(xvals)}
    indexes_y = {v: i for i, v in enumerate(yvals)}
    for row in nodes:
        z[indexes_x[float(row["TS1_s_eV"])], indexes_y[float(row["TS2_s_eV"])]] = float(row["log10_TOF_C6H12O"])
    if not np.isfinite(z).all():
        raise RuntimeError("Official CatMAP grid includes missing/non-finite nodes")
    xyz_path = output_dir / f"catmap_TS1_TS2_{len(xvals)}x{len(yvals)}_xyz.csv"
    with xyz_path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["TS1_s_eV", "TS2_s_eV", "CatMAP_log10_TOF_C6H12O"])
        for row in nodes:
            w.writerow([fmt(row["TS1_s_eV"]), fmt(row["TS2_s_eV"]), fmt(row["log10_TOF_C6H12O"])])
    matrix_path = output_dir / f"catmap_TS1_TS2_{len(xvals)}x{len(yvals)}_matrix.csv"
    with matrix_path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["TS1\\TS2", *[format(v, ".17g") for v in yvals]])
        for i, xv in enumerate(xvals):
            w.writerow([format(xv, ".17g"), *[format(float(v), ".17g") for v in z[i]]])
    return np.asarray(xvals), np.asarray(yvals), z


def label_real_points(ax, points):
    # Place every validated catalyst label without altering its coordinates.
    offsets = [(8, 8), (8, 19), (8, -18), (-46, 8), (-48, -18),
               (-46, 20), (0, 25), (0, -26), (26, 0), (-55, 0),
               (24, 18), (24, -19), (-54, 20), (-54, -20),
               (12, 30), (-12, -30)]
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer)
    used = []
    for p in points:
        xy = (p["TS1_s_eV"], p["TS2_s_eV"])
        candidates = []
        for offset in offsets:
            trial = ax.annotate(p["catalyst"], xy, xytext=offset, textcoords="offset points",
                                fontsize=8, arrowprops={"arrowstyle": "-", "lw": .35, "color": "#666"})
            fig.canvas.draw()
            box = trial.get_window_extent(renderer).expanded(1.04, 1.05)
            outside = box.x0 < axes_box.x0 or box.x1 > axes_box.x1 or box.y0 < axes_box.y0 or box.y1 > axes_box.y1
            overlap = any(box.overlaps(previous) for previous in used)
            score = 1000 * outside + 100 * overlap + abs(offset[0]) + abs(offset[1])
            candidates.append((score, offset, box))
            trial.remove()
            if not outside and not overlap:
                break
        _, offset, _ = min(candidates, key=lambda item: item[0])
        annotation = ax.annotate(p["catalyst"], xy, xytext=offset, textcoords="offset points",
                                 fontsize=8, color="#162b4d", zorder=8,
                                 arrowprops={"arrowstyle": "-", "lw": .35, "color": "#666"})
        fig.canvas.draw()
        used.append(annotation.get_window_extent(renderer).expanded(1.04, 1.05))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, nargs=2, default=[5, 5])
    parser.add_argument("--tag", default=None)
    parser.add_argument("--reuse-cache", action="store_true",
                        help="Rebuild exports/plots from the existing CatMAP pickle without re-solving")
    args = parser.parse_args()
    if __version__ != "0.3.1":
        raise RuntimeError(f"Expected CatMAP 0.3.1; found {__version__}")
    nx, ny = args.resolution
    tag = args.tag or f"{nx}x{ny}"
    out = Path("runs") / tag
    out.mkdir(parents=True, exist_ok=True)
    model = ReactionModel(setup_file="ts1_ts2_descriptor_grid.mkm")
    model.resolution = [nx, ny]
    model.data_file = str(out / f"catmap_ts1_ts2_{tag}.pkl")
    model.log_file = str(out / f"catmap_ts1_ts2_{tag}.log")
    if args.reuse_cache:
        model.load_data_file(overwrite=True)
        if model.turnover_frequency_map is None or model.coverage_map is None:
            raise RuntimeError("Requested cached CatMAP maps are missing from the pickle")
        cached_summary_path = out / "grid_summary.json"
        if not cached_summary_path.exists():
            raise RuntimeError("Cached run summary is required to reuse CatMAP maps")
        cached_summary = json.loads(cached_summary_path.read_text(encoding="utf-8"))
        cached_max_residual = mp.mpf(cached_summary["maximum_recomputed_steady_state_residual"])
        if cached_max_residual > mp.mpf(str(model.tolerance)):
            raise RuntimeError("Cached CatMAP nodes did not pass the strict residual check")
    else:
        model.run(recalculate=True)

    labels = model.output_labels.get("turnover_frequency") or list(model.gas_names)
    model.output_labels = dict(model.output_labels)
    model.output_labels["turnover_frequency"] = labels
    product_idx = labels.index("C6H12O_g")
    tof_map = model.turnover_frequency_map
    if len(tof_map) != nx * ny:
        raise RuntimeError(f"Expected {nx*ny} TOF nodes; got {len(tof_map)}")
    coverage_map = model.coverage_map
    if len(coverage_map) != nx * ny:
        raise RuntimeError("CatMAP did not return coverages at every grid node")
    coverage_by_point = {tuple(float(v) for v in point): values for point, values in coverage_map}
    if set(coverage_by_point) != {tuple(float(v) for v in point) for point, _ in tof_map}:
        raise RuntimeError("CatMAP coverage and TOF maps do not share the same complete descriptor coordinates")
    nodes = []
    residuals = []
    for descriptor, values in tof_map:
        point_key = tuple(float(v) for v in descriptor)
        if not args.reuse_cache:
            model.solver._descriptors = list(point_key)
            residual = model.solver.get_residual(
                coverage_by_point[point_key], validate_coverages=False,
                refresh_rate_constants=True)
            residuals.append(residual)
            if residual > model.solver.tolerance:
                raise RuntimeError(f"Official residual exceeds tolerance at {descriptor}: {residual}")
        coverages = coverage_by_point[point_key]
        if any(value < 0 for value in coverages):
            raise RuntimeError(f"Negative coverage at {descriptor}: {coverages}")
        tof = values[product_idx]
        if tof <= 0:
            raise RuntimeError(f"Non-positive product TOF at {descriptor}: {tof}")
        log_tof = mp.log10(mp.mpf(str(tof)))
        nodes.append({"TS1_s_eV": float(descriptor[0]), "TS2_s_eV": float(descriptor[1]),
                      "TOF_C6H12O": str(tof), "log10_TOF_C6H12O": float(log_tof)})
    actual_points = load_validated_points()
    with (out / "validated_catalyst_points.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(actual_points[0]))
        w.writeheader()
        w.writerows(actual_points)
    xgrid, ygrid, z = write_grid(nodes, out)

    max_idx = np.unravel_index(np.argmax(z), z.shape)
    peak_x, peak_y, peak_z = float(xgrid[max_idx[0]]), float(ygrid[max_idx[1]]), float(z[max_idx])
    boundary = max_idx[0] in {0, len(xgrid)-1} or max_idx[1] in {0, len(ygrid)-1}
    peak = {"TS1_s_eV": peak_x, "TS2_s_eV": peak_y, "CatMAP_log10_TOF_C6H12O": peak_z,
            "CatMAP_TOF_C6H12O": float(10**peak_z), "grid_resolution": [nx, ny],
            "boundary_maximum": bool(boundary),
            "classification": "maximum within sampled range (boundary)" if boundary else "interior grid maximum; local neighborhood requires inspection"}
    with (out / "catmap_grid_maximum.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(peak))
        w.writeheader(); w.writerow(peak)

    # Native CatMAP analysis output. The validated catalyst markers/labels are
    # a minimal overlay; the underlying map is CatMAP VectorMap data.
    vector = VectorMap(model)
    vector.plot_variable = "turnover_frequency"
    vector.include_labels = ["C6H12O_g"]
    vector.plot_mode = "separate"
    # VectorMap 0.3.1's map_to_array interpolates all values cubically; applying
    # log_scale after that can create log10 warnings from small cubic undershoot.
    # Pass only the selected product channel transformed to log10 before the
    # official plotter; the original ReactionModel outputs and exported TOFs
    # remain untouched.
    vector.turnover_frequency_map = [
        [point, [float(mp.log10(mp.mpf(str(values[product_idx]))))]]
        for point, values in model.turnover_frequency_map
    ]
    vector.output_labels = dict(model.output_labels)
    vector.output_labels["turnover_frequency"] = ["C6H12O_g"]
    vector.log_scale = False
    # We add validated-catalyst labels from the baseline table below; suppress
    # duplicate lowercase anchor labels that VectorMap would otherwise add.
    vector.descriptor_dict = {}
    vector.plot_size = 8
    # CatMAP 0.3.1 VectorMap's auto-range scan includes unselected gas labels,
    # whose net TOFs can be negative. Bound the log colour scale to the selected
    # positive product channel to avoid log10(negative/zero) NaNs.
    product_logs = [float(row["log10_TOF_C6H12O"]) for row in nodes]
    vector.min = min(product_logs)
    vector.max = max(product_logs)
    # CatMAP's native VectorMap uses ReactionModel.map_to_array to interpolate
    # its solver nodes for display. Request 400x400 here for the Origin-ready
    # visualization matrix; these are official-map interpolates, not 160,000
    # independent steady-state solves. The original solver-node table is saved
    # separately above and remains the source for the discrete grid maximum.
    vector.resolution = [400, 400]
    vector_array = model.map_to_array(
        vector.turnover_frequency_map, model.descriptor_ranges, [400, 400],
        log_interpolate=False)
    vector_z_yx = vector_array[:, :, 0]
    if not np.isfinite(vector_z_yx).all():
        raise RuntimeError("CatMAP VectorMap 400x400 interpolation contains NaN/Inf")
    vector_x = np.linspace(model.descriptor_ranges[0][0], model.descriptor_ranges[0][1], 400)
    vector_y = np.linspace(model.descriptor_ranges[1][0], model.descriptor_ranges[1][1], 400)
    interp_max_idx = np.unravel_index(np.argmax(vector_z_yx), vector_z_yx.shape)
    interp_peak_x = float(vector_x[interp_max_idx[1]])
    interp_peak_y = float(vector_y[interp_max_idx[0]])
    interp_peak_z = float(vector_z_yx[interp_max_idx])
    interp_peak = {
        "TS1_s_eV": interp_peak_x, "TS2_s_eV": interp_peak_y,
        "CatMAP_VectorMap_interpolated_log10_TOF_C6H12O": interp_peak_z,
        "interpolated_TOF_C6H12O": float(10**interp_peak_z),
        "source": "CatMAP 0.3.1 ReactionModel.map_to_array via VectorMap; not an independent solver node",
    }
    with (out / "catmap_VectorMap_interpolated_400x400_xyz.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["TS1_s_eV", "TS2_s_eV", "CatMAP_VectorMap_interpolated_log10_TOF_C6H12O"])
        for i, xv in enumerate(vector_x):
            for j, yv in enumerate(vector_y):
                w.writerow([format(float(xv), ".17g"), format(float(yv), ".17g"),
                            format(float(vector_z_yx[j, i]), ".17g")])
    with (out / "catmap_VectorMap_interpolated_400x400_matrix.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["TS1\\TS2", *[format(float(v), ".17g") for v in vector_y]])
        for i, xv in enumerate(vector_x):
            w.writerow([format(float(xv), ".17g"), *[format(float(v), ".17g") for v in vector_z_yx[:, i]]])
    with (out / "catmap_VectorMap_interpolated_maximum.csv").open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(interp_peak))
        w.writeheader(); w.writerow(interp_peak)
    vector.descriptor_labels = ["TS1_s formation energy (eV)", "TS2_s formation energy (eV)"]
    vector.output_file = str(out / f"catmap_VectorMap_{tag}.pdf")
    fig = vector.plot(save=False)
    ax = fig.axes[0]
    pad_x = max((model.descriptor_ranges[0][1] - model.descriptor_ranges[0][0]) * .018, .005)
    pad_y = max((model.descriptor_ranges[1][1] - model.descriptor_ranges[1][0]) * .018, .005)
    ax.set_xlim(model.descriptor_ranges[0][0] - pad_x, model.descriptor_ranges[0][1] + pad_x)
    ax.set_ylim(model.descriptor_ranges[1][0] - pad_y, model.descriptor_ranges[1][1] + pad_y)
    real_x = [p["TS1_s_eV"] for p in actual_points]
    real_y = [p["TS2_s_eV"] for p in actual_points]
    ax.scatter(real_x, real_y, marker="D", s=55, facecolor="white", edgecolor="#111827",
               zorder=7, label="13 validated strict-baseline catalysts")
    label_real_points(ax, actual_points)
    ax.scatter([peak_x], [peak_y], marker="*", s=280, color="#ff3434", edgecolor="black",
               zorder=10, label="Actual CatMAP solver-node maximum")
    ax.scatter([interp_peak_x], [interp_peak_y], marker="o", s=95, facecolor="none",
               edgecolor="#ff8c00", linewidth=2, zorder=9, label="VectorMap 400x400 interpolated maximum")
    ax.legend(loc="best", fontsize=8)
    ax.set_title(f"CatMAP 0.3.1 VectorMap: TS1_s + TS2_s\n400×400 official map interpolation from {nx}×{ny} converged solver nodes")
    fig.savefig(out / f"catmap_VectorMap_400x400_from_{tag}.pdf")
    fig.savefig(out / f"catmap_VectorMap_400x400_from_{tag}.png", dpi=300)
    plt.close(fig)

    # Contours are rendered directly from all official CatMAP grid nodes,
    # without adding interpolated or fitted activity values.
    X, Y = np.meshgrid(xgrid, ygrid, indexing="ij")
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    cf = ax.contourf(X, Y, z, levels=40, cmap="viridis")
    fig.colorbar(cf, ax=ax, label="CatMAP log10(TOF of C6H12O)")
    ax.scatter(real_x, real_y, marker="D", s=58, facecolor="white", edgecolor="#111827",
               zorder=7, label="13 validated strict-baseline catalysts")
    label_real_points(ax, actual_points)
    ax.scatter([peak_x], [peak_y], marker="*", s=280, color="#ff3434", edgecolor="black",
               zorder=10, label="CatMAP grid maximum")
    ax.annotate(f"{'Boundary maximum' if boundary else 'Grid maximum'}\nTS1={peak_x:.6f} eV\nTS2={peak_y:.6f} eV\nlog10TOF={peak_z:.6f}",
                (peak_x, peak_y), xytext=(.03, .96), textcoords="axes fraction", ha="left", va="top",
                bbox={"boxstyle": "round,pad=.35", "fc": "white", "alpha": .94},
                arrowprops={"arrowstyle": "->", "lw": .8})
    pad_x = max((xgrid[-1]-xgrid[0])*.018, 0.005)
    pad_y = max((ygrid[-1]-ygrid[0])*.018, 0.005)
    ax.set_xlim(xgrid[0]-pad_x, xgrid[-1]+pad_x)
    ax.set_ylim(ygrid[0]-pad_y, ygrid[-1]+pad_y)
    ax.set_xlabel("TS1_s formation energy (eV)")
    ax.set_ylabel("TS2_s formation energy (eV)")
    ax.set_title(f"CatMAP 0.3.1 official TS1_s + TS2_s grid ({nx}×{ny})\nGeneralizedLinearScaler + SteadyStateSolver; grid nodes only")
    ax.legend(loc="best", fontsize=8)
    fig.savefig(out / f"catmap_TS1_TS2_activity_{tag}.png", dpi=300)
    fig.savefig(out / f"catmap_TS1_TS2_activity_{tag}.pdf")
    plt.close(fig)

    max_residual = max(residuals) if residuals else cached_max_residual
    summary = {
        "catmap_version": __version__, "resolution": [nx, ny], "computed_nodes": len(nodes),
        "surface_names": list(model.surface_names), "descriptor_names": list(model.descriptor_names),
        "temperature_K": model.temperature, "reservoir_concentrations": {
            "C6H12_g": 1.0, "H2O2_g": 0.185, "H2O_g": 0.815, "C6H12O_g": 1e-20},
        "scaler": type(model.scaler).__name__, "solver": type(model.solver).__name__,
        "mapper": type(model.mapper).__name__, "tolerance": str(model.solver.tolerance),
        "maximum_recomputed_steady_state_residual": str(max_residual),
        "decimal_precision": model.decimal_precision, "avoid_scaling": model.scaler.avoid_scaling,
        "descriptor_decimal_precision": model.descriptor_decimal_precision,
        "max_bisections": model.mapper.max_bisections, "max_initial_guesses": model.mapper.max_initial_guesses,
        "all_grid_nodes_positive_and_finite": True, "maximum": peak,
        "official_VectorMap_400x400_interpolated_maximum": interp_peak,
        "vector_map_interpolation": "CatMAP 0.3.1 ReactionModel.map_to_array cubic interpolation used by VectorMap",
        "vector_map_400x400_is_independent_solver_run": False,
        "grid_descriptor_ranges": model.descriptor_ranges,
        "actual_catalysts_overlaid": [p["catalyst"] for p in actual_points],
        "outputs": [p.name for p in sorted(out.iterdir()) if p.suffix in {".csv", ".pdf", ".png"}],
    }
    (out / "grid_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
