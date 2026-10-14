"""Reproduce and plot the empirical TS1/TS2 quadratic surface from frozen inputs."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "empirical"
LABELS = {
    "ti": "Ti", "tico": "Ti-Co", "tini": "Ti-Ni", "tife": "Ti-Fe",
    "tire": "Ti-Re", "timo": "Ti-Mo", "timn": "Ti-Mn", "ticr": "Ti-Cr",
    "tiw": "Ti-W", "tiv": "Ti-V", "tita": "Ti-Ta", "tihf": "Ti-Hf", "titi": "Ti-Ti",
}


def load_points():
    with (DATA / "energies.txt").open(encoding="utf-8-sig", newline="") as stream:
        energy_rows = list(csv.DictReader(stream, delimiter="\t"))
    energies = {}
    for row in energy_rows:
        surface, state = row["surface_name"], row["species_name"]
        if surface in LABELS and state in {"TS1", "TS2"}:
            key = surface, state
            if key in energies:
                raise ValueError(f"Duplicate energy record: {key}")
            energies[key] = float(row["formation_energy"])
    with (DATA / "matched_data.csv").open(encoding="utf-8-sig", newline="") as stream:
        tof_rows = list(csv.DictReader(stream))
    if len(tof_rows) != 13 or {r["surface_name"] for r in tof_rows} != set(LABELS):
        raise ValueError("Expected exactly the 13 validated sites; Ti-Zr must be excluded")
    if len(energies) != 26:
        raise ValueError(f"Expected 26 TS1/TS2 values for 13 surfaces; found {len(energies)}")
    points = []
    for row in tof_rows:
        surface = row["surface_name"]
        if row["status"] != "CONVERGED" or float(row["temperature_K"]) != 333.15:
            raise ValueError(f"TOF is not a validated 333.15 K result: {surface}")
        tof = float(row["TOF"])
        log_tof = float(row["log10_TOF"])
        if tof <= 0 or abs(log_tof - math.log10(tof)) > 1e-10:
            raise ValueError(f"TOF/log10(TOF) mismatch: {surface}")
        points.append({
            "surface_name": surface,
            "catalyst": LABELS[surface],
            "TS1_s_eV": energies[surface, "TS1"],
            "TS2_s_eV": energies[surface, "TS2"],
            "TOF_C6H12O": tof,
            "log10_TOF_C6H12O": log_tof,
        })
    return points


def design_matrix(x, y, scales):
    xm, xs, ym, ys = scales
    u, v = (x-xm)/xs, (y-ym)/ys
    return np.column_stack([np.ones_like(u), u, v, u*u, u*v, v*v])


def in_polygon(vertices, point):
    p = np.asarray(point, float)
    edge = np.roll(vertices, -1, axis=0) - vertices
    cross = edge[:, 0]*(p[1]-vertices[:, 1]) - edge[:, 1]*(p[0]-vertices[:, 0])
    return bool(np.all(cross >= -1e-10) or np.all(cross <= 1e-10))


def label_points(ax, x, y, names):
    offsets = [(7, 8), (8, 18), (8, -17), (-42, 8), (-42, -18),
               (-42, 18), (0, 23), (0, -25), (24, 0), (-52, 0),
               (22, 16), (22, -17), (-52, 17), (-52, -17)]
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer)
    used = []
    for xx, yy, name in zip(x, y, names):
        opts = []
        for off in offsets:
            ann = ax.annotate(name, (xx, yy), xytext=off, textcoords="offset points",
                              fontsize=8, arrowprops={"arrowstyle": "-", "lw": .35, "color": "#666"})
            fig.canvas.draw()
            box = ann.get_window_extent(renderer).expanded(1.04, 1.05)
            outside = box.x0 < axes_box.x0 or box.x1 > axes_box.x1 or box.y0 < axes_box.y0 or box.y1 > axes_box.y1
            overlap = any(box.overlaps(old) for old in used)
            opts.append((1000*outside + 100*overlap + abs(off[0]) + abs(off[1]), off, box))
            ann.remove()
            if not outside and not overlap:
                break
        _, off, _ = min(opts, key=lambda z: z[0])
        ann = ax.annotate(name, (xx, yy), xytext=off, textcoords="offset points", fontsize=8,
                          color="#162b4d", zorder=7,
                          arrowprops={"arrowstyle": "-", "lw": .35, "color": "#666"})
        fig.canvas.draw()
        used.append(ann.get_window_extent(renderer).expanded(1.04, 1.05))


def write_csv(path, fields, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    points = load_points()
    x = np.array([r["TS1_s_eV"] for r in points], float)
    y = np.array([r["TS2_s_eV"] for r in points], float)
    z = np.array([r["log10_TOF_C6H12O"] for r in points], float)
    scales = (float(x.mean()), float(x.std()), float(y.mean()), float(y.std()))
    matrix = design_matrix(x, y, scales)
    coef, _, rank, singular = np.linalg.lstsq(matrix, z, rcond=None)
    if rank != 6:
        raise ValueError(f"Quadratic design is rank deficient: {rank}/6")
    fitted = matrix @ coef
    r2 = 1 - float(np.sum((fitted-z)**2)/np.sum((z-z.mean())**2))
    residual_rmse = float(np.sqrt(np.mean((fitted-z)**2)))
    hessian = np.array([[2*coef[3], coef[4]], [coef[4], 2*coef[5]]])
    stat_uv = np.linalg.solve(hessian, -coef[1:3])
    peak = np.array([scales[0]+scales[1]*stat_uv[0], scales[2]+scales[3]*stat_uv[1]])
    peak_log = float(np.array([1, *stat_uv, stat_uv[0]**2, stat_uv.prod(), stat_uv[1]**2]) @ coef)
    eig = np.linalg.eigvalsh(hessian)
    if not np.all(eig < 0):
        raise ValueError(f"Fitted stationary point is not a two-direction concave peak: Hessian eigenvalues={eig}")
    points_xy = np.column_stack([x, y])
    hull_indices = ConvexHull(points_xy).vertices
    hull = points_xy[hull_indices]
    in_rect = bool(x.min() < peak[0] < x.max() and y.min() < peak[1] < y.max())
    in_hull = in_polygon(hull, peak)
    if not in_rect:
        raise ValueError(f"Reported empirical peak is not inside observed-coordinate rectangle: {peak}")

    nx = ny = 400
    xgrid = np.linspace(float(x.min()), float(x.max()), nx)
    ygrid = np.linspace(float(y.min()), float(y.max()), ny)
    XX, YY = np.meshgrid(xgrid, ygrid, indexing="ij")
    ZZ = (design_matrix(XX.ravel(), YY.ravel(), scales) @ coef).reshape(nx, ny)
    if not np.isfinite(ZZ).all():
        raise ValueError("Empirical grid contains non-finite values")

    # Origin XYZ and matrix; matrix first column is TS1 and top row after corner is TS2.
    write_csv(OUT / "empirical_TS1_TS2_fullgrid_400x400_xyz.csv",
              ["TS1_s_eV", "TS2_s_eV", "predicted_log10_TOF"],
              [{"TS1_s_eV": format(float(a), ".17g"), "TS2_s_eV": format(float(b), ".17g"),
                "predicted_log10_TOF": format(float(c), ".17g")}
               for a, b, c in zip(XX.ravel(), YY.ravel(), ZZ.ravel())])
    with (OUT / "empirical_TS1_TS2_fullgrid_400x400_matrix.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["TS1\\TS2", *[format(float(v), ".17g") for v in ygrid]])
        for i, xv in enumerate(xgrid):
            writer.writerow([format(float(xv), ".17g"), *[format(float(v), ".17g") for v in ZZ[i, :]]])
    write_csv(OUT / "validated_catalyst_points.csv", list(points[0].keys()), points)
    peak_row = {
        "TS1_s_eV": float(peak[0]), "TS2_s_eV": float(peak[1]),
        "predicted_log10_TOF": peak_log, "predicted_TOF": 10**peak_log,
        "inside_observed_rectangle": in_rect, "inside_catalyst_convex_hull": in_hull,
        "peak_type": "empirical quadratic stationary maximum; outside catalyst convex hull" if not in_hull else "empirical quadratic stationary maximum inside catalyst convex hull",
    }
    write_csv(OUT / "empirical_peak_summary.csv", list(peak_row.keys()), [peak_row])
    with (OUT / "empirical_quadratic_coefficients.json").open("w", encoding="utf-8") as stream:
        json.dump({"model": "z=c0+c1*u+c2*v+c3*u^2+c4*u*v+c5*v^2; u=(TS1-mean)/std; v=(TS2-mean)/std",
                   "coefficients_standardized": coef.tolist(), "scales": dict(zip(["TS1_mean", "TS1_std", "TS2_mean", "TS2_std"], scales)),
                   "n": len(points), "in_sample_R2": r2, "in_sample_RMSE_log10TOF": residual_rmse,
                   "hessian_eigenvalues": eig.tolist(), "peak": peak_row,
                   "descriptor_ranges_eV": {"TS1_min": float(x.min()), "TS1_max": float(x.max()), "TS2_min": float(y.min()), "TS2_max": float(y.max())}},
                  stream, ensure_ascii=False, indent=2)

    # Empirical two-dimensional contour; full rectangular domain, no hull mask.
    fig, ax = plt.subplots(figsize=(10.0, 8.0), constrained_layout=True)
    contour = ax.contourf(XX, YY, ZZ, levels=60, cmap="viridis")
    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label("Empirical predicted log10(TOF of C6H12O)")
    ax.scatter(x, y, s=62, marker="D", facecolor="white", edgecolor="#111827", zorder=6,
               label="13 validated catalyst points")
    label_points(ax, x, y, [r["catalyst"] for r in points])
    ax.scatter([peak[0]], [peak[1]], marker="*", s=310, color="#ff3434", edgecolor="black", zorder=10,
               label="Empirical quadratic peak")
    ax.annotate(f"Empirical peak\nTS1={peak[0]:.6f} eV\nTS2={peak[1]:.6f} eV\nlog10TOF={peak_log:.6f}\nOutside catalyst hull",
                peak, xytext=(.03, .96), textcoords="axes fraction", ha="left", va="top", fontsize=9,
                bbox={"boxstyle": "round,pad=.35", "fc": "white", "alpha": .94},
                arrowprops={"arrowstyle": "->", "lw": .8})
    # Grid evaluation remains exactly on the observed rectangle; expand only
    # the visible axes slightly so markers at extrema are not clipped.
    pad_x = (x.max()-x.min()) * 0.015
    pad_y = (y.max()-y.min()) * 0.015
    ax.set(xlabel="TS1_s formation energy (eV)", ylabel="TS2_s formation energy (eV)",
           xlim=(x.min()-pad_x, x.max()+pad_x), ylim=(y.min()-pad_y, y.max()+pad_y),
           title="Empirical quadratic response surface: TS1_s + TS2_s\nNot a CatMAP continuous-grid calculation")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=.15)
    fig.savefig(OUT / "empirical_TS1_TS2_contour.png", dpi=300)
    fig.savefig(OUT / "empirical_TS1_TS2_contour.pdf")
    plt.close(fig)

    # 3D view is the same empirical surface; sparse mesh rendering preserves clarity.
    fig = plt.figure(figsize=(11.0, 8.5), constrained_layout=True)
    ax3 = fig.add_subplot(111, projection="3d", computed_zorder=False)
    surf = ax3.plot_surface(XX, YY, ZZ, rstride=6, cstride=6, cmap="viridis", linewidth=0, antialiased=True, alpha=.88, zorder=1)
    ax3.scatter(x, y, z, marker="D", color="white", edgecolor="black", s=35, depthshade=False, zorder=20, label="13 validated CatMAP TOF points")
    for row, xv, yv, zv in zip(points, x, y, z):
        ax3.text(xv, yv, zv, row["catalyst"], fontsize=7)
    ax3.scatter([peak[0]], [peak[1]], [peak_log], marker="*", color="#ff3434", edgecolor="black", s=420, depthshade=False, zorder=1000, label="Empirical quadratic peak")
    ax3.text(peak[0], peak[1], peak_log, "  fitted peak", color="#b30000", fontsize=9, zorder=1001)
    ax3.set(xlabel="TS1_s (eV)", ylabel="TS2_s (eV)", zlabel="Empirical log10(TOF)",
            title="TS1_s + TS2_s empirical quadratic response surface\nNot CatMAP grid; peak outside catalyst convex hull")
    fig.colorbar(surf, ax=ax3, shrink=.62, pad=.12, label="Predicted log10(TOF of C6H12O)")
    ax3.view_init(elev=28, azim=-55)
    ax3.legend(loc="upper left", fontsize=8)
    fig.savefig(OUT / "empirical_TS1_TS2_surface3d.png", dpi=300)
    fig.savefig(OUT / "empirical_TS1_TS2_surface3d.pdf")
    plt.close(fig)

    # Two slices through the fitted stationary point; preserve the fixed coordinate.
    profile_records = {}
    for axis, variable, fixed, vals in (("TS1", "TS1_s_eV", "TS2_s_eV", xgrid),
                                        ("TS2", "TS2_s_eV", "TS1_s_eV", ygrid)):
        if axis == "TS1":
            xp, yp = vals, np.full_like(vals, peak[1])
        else:
            xp, yp = np.full_like(vals, peak[0]), vals
        zp = design_matrix(xp, yp, scales) @ coef
        recs = [{variable: float(a), fixed: float(b), "predicted_log10_TOF": float(c)} for a, b, c in zip(xp, yp, zp)]
        profile_records[axis] = recs
        write_csv(OUT / f"profile_{axis}_through_peak.csv", [variable, fixed, "predicted_log10_TOF"], recs)
        fig, ax = plt.subplots(figsize=(8.2, 5.8), constrained_layout=True)
        ax.plot(vals, zp, color="#225ea8", lw=2.4, label="Empirical quadratic slice")
        peak_value = float(design_matrix(np.array([peak[0]]), np.array([peak[1]]), scales)[0] @ coef)
        ax.scatter([peak[0] if axis == "TS1" else peak[1]], [peak_value], marker="*", color="#ff3434", edgecolor="black", s=180, zorder=5, label="Fitted peak")
        ax.set(xlabel=f"{axis}_s formation energy (eV)", ylabel="Predicted log10(TOF of C6H12O)",
               title=f"Empirical profile through fitted peak; fixed {fixed}={peak[1] if axis == 'TS1' else peak[0]:.6f} eV")
        ax.grid(alpha=.2)
        ax.legend()
        fig.savefig(OUT / f"empirical_{axis}_profile.png", dpi=300)
        fig.savefig(OUT / f"empirical_{axis}_profile.pdf")
        plt.close(fig)

    summary = {
        "n_catalysts": 13, "temperature_K": 333.15, "in_sample_R2": r2,
        "in_sample_RMSE_log10TOF": residual_rmse, "descriptor_rectangle_eV": {
            "TS1_min": float(x.min()), "TS1_max": float(x.max()), "TS2_min": float(y.min()), "TS2_max": float(y.max())},
        "hessian_eigenvalues": eig.tolist(), "peak": peak_row,
        "peak_nearest_catalysts_standardized_distance": sorted((float(np.hypot((x[i]-peak[0])/scales[1], (y[i]-peak[1])/scales[3])), points[i]["catalyst"]) for i in range(13))[:3],
    }
    (OUT / "empirical_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "EMPIRICAL_NOTES.md").write_text(
        "# TS1_s + TS2_s 经验二次响应面\n\n"
        "本结果使用归档 `energies.txt` 的 TS1/TS2 formation_energy 与按 surface_name 匹配的 13 个严格收敛 C6H12O TOF；排除 Ti–Zr。模型为普通六参数二次响应面，基于 13 点最小二乘拟合，不是 CatMAP 网格。图与 400×400 Origin 网格均在全部 13 点坐标构成的矩形内直接评估，不做凸包遮罩。\n\n"
        f"样本内 R²={r2:.8f}；RMSE={residual_rmse:.8f} log10(TOF)。Hessian 特征值={eig[0]:.8f}, {eig[1]:.8f}（两者为负才支持局部下凹）。拟合驻点：TS1={peak[0]:.8f} eV，TS2={peak[1]:.8f} eV，log10TOF={peak_log:.8f}，预测 TOF={10**peak_log:.6g}；矩形内={in_rect}，催化剂凸包内={in_hull}。凸包外峰仅是经验响应面的数学峰，不能解读成已实现催化剂活性。\n\n"
        "输入的 O/Ha/Hb 历史形成能不参与本模型；TS1/TS2 来自 CatMAP energies.txt。经验图不代表 CatMAP 官方连续网格。实际 13 点及拟合系数均单独导出。\n",
        encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
