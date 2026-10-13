"""Export the previously fitted O+Ha and O+Hb quadratic surfaces on full grids.

This script reads the frozen matched data and coefficients in fit_metrics.json.
It does not refit, run CatMAP, or modify any source/result files already present.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull


OUT = Path(__file__).resolve().parent
METRICS_PATH = OUT / "fit_metrics.json"
DATA_PATH = OUT / "matched_data.csv"
EXPECTED_CATALYSTS = {
    "Ti", "Ti-Co", "Ti-Ni", "Ti-Fe", "Ti-Re", "Ti-Mo", "Ti-Mn",
    "Ti-Cr", "Ti-W", "Ti-V", "Ti-Ta", "Ti-Hf", "Ti-Ti",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_frozen_inputs():
    metrics = json.loads(METRICS_PATH.read_text(encoding="utf-8"))
    with DATA_PATH.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 13 or len({r["surface_name"] for r in rows}) != 13:
        raise ValueError("matched_data.csv must contain 13 unique surface_name records")
    if {r["catalyst"] for r in rows} != EXPECTED_CATALYSTS:
        raise ValueError("Catalyst names do not match the 13 validated sites")
    if any(r["status"] != "CONVERGED" or float(r["temperature_K"]) != 333.15 for r in rows):
        raise ValueError("Input records must be the existing converged 333.15 K results")
    return rows, metrics


def quadratic(u, v, coef):
    c0, c1, c2, c3, c4, c5 = coef
    return c0 + c1*u + c2*v + c3*u*u + c4*u*v + c5*v*v


def evaluate(x, y, coef, scales):
    xm, xs, ym, ys = scales
    return quadratic((np.asarray(x)-xm)/xs, (np.asarray(y)-ym)/ys, coef)


def exact_maximum(coef, scales, polygon):
    """Maximize the same quadratic on a rectangle or a convex polygon."""
    xm, xs, ym, ys = scales
    points = np.asarray(polygon, dtype=float)
    candidates = [(float(x), float(y)) for x, y in points]
    c0, c1, c2, c3, c4, c5 = coef
    hessian = np.array([[2*c3, c4], [c4, 2*c5]], dtype=float)
    gradient = np.array([c1, c2], dtype=float)

    # A stationary point can be a maximum only for a negative-definite Hessian.
    if np.all(np.linalg.eigvalsh(hessian) < 0) and abs(np.linalg.det(hessian)) > 1e-14:
        uv = np.linalg.solve(hessian, -gradient)
        xy = np.array([xm + xs*uv[0], ym + ys*uv[1]])
        if polygon_contains(points, xy):
            candidates.append((float(xy[0]), float(xy[1])))

    # Every boundary segment is a one-dimensional quadratic; include its exact
    # interior maximum where the curvature is concave.
    for p0, p1 in zip(points, np.roll(points, -1, axis=0)):
        midpoint = (p0 + p1) / 2
        q0 = float(evaluate(p0[0], p0[1], coef, scales))
        q1 = float(evaluate(p1[0], p1[1], coef, scales))
        qm = float(evaluate(midpoint[0], midpoint[1], coef, scales))
        a = 2 * (q1 + q0 - 2*qm)
        b = q1 - q0 - a
        if a < -1e-14:
            t = -b / (2*a)
            if 0 < t < 1:
                p = p0 + t*(p1-p0)
                candidates.append((float(p[0]), float(p[1])))
    scored = [(float(evaluate(x, y, coef, scales)), x, y) for x, y in candidates]
    z, x, y = max(scored, key=lambda item: item[0])
    return {"x_peak": x, "y_peak": y, "predicted_log10_TOF_peak": z,
            "predicted_TOF_peak": 10.0**z}


def polygon_contains(vertices, point):
    # Convex polygon containment, including its boundary, by consistent cross
    # product signs. This avoids treating a boundary stationary point as outside.
    p = np.asarray(point, dtype=float)
    edges = np.roll(vertices, -1, axis=0) - vertices
    cross = edges[:, 0]*(p[1]-vertices[:, 1]) - edges[:, 1]*(p[0]-vertices[:, 0])
    return bool(np.all(cross >= -1e-10) or np.all(cross <= 1e-10))


def point_on_boundary(vertices, point, tolerance=1e-8):
    p = np.asarray(point, dtype=float)
    for a, b in zip(vertices, np.roll(vertices, -1, axis=0)):
        segment = b-a
        length2 = float(segment @ segment)
        t = min(1.0, max(0.0, float((p-a) @ segment / length2)))
        if np.linalg.norm(p-(a+t*segment)) <= tolerance:
            return True
    return False


def label_points(ax, x, y, rows):
    offsets = [(7, 8), (8, 18), (8, -18), (-42, 8), (-42, -18),
               (-42, 18), (0, 22), (0, -24), (24, 0), (-52, 0),
               (22, 16), (22, -17), (-52, 17), (-52, -17)]
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer)
    occupied = []
    for xx, yy, row in zip(x, y, rows):
        choices = []
        for offset in offsets:
            ann = ax.annotate(row["catalyst"], (xx, yy), xytext=offset,
                              textcoords="offset points", fontsize=8,
                              color="#172b4d", zorder=7,
                              arrowprops={"arrowstyle": "-", "lw": .4, "color": "#777"})
            ax.figure.canvas.draw()
            box = ann.get_window_extent(renderer).expanded(1.04, 1.06)
            outside = box.x0 < axes_box.x0 or box.x1 > axes_box.x1 or box.y0 < axes_box.y0 or box.y1 > axes_box.y1
            overlaps = any(box.overlaps(prev) for prev in occupied)
            choices.append((1000*outside + 100*overlaps + abs(offset[0]) + abs(offset[1]), offset, box))
            ann.remove()
            if not outside and not overlaps:
                break
        _, offset, _ = min(choices, key=lambda choice: choice[0])
        ann = ax.annotate(row["catalyst"], (xx, yy), xytext=offset,
                          textcoords="offset points", fontsize=8,
                          color="#172b4d", zorder=7,
                          arrowprops={"arrowstyle": "-", "lw": .4, "color": "#777"})
        ax.figure.canvas.draw()
        occupied.append(ann.get_window_extent(renderer).expanded(1.04, 1.06))


def write_csv(path, fieldnames, records):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def validate_exports():
    for prefix in ("O_Ha", "O_Hb"):
        for n in (200, 400):
            xyz_path = OUT / f"{prefix}_empirical_fullgrid_{n}x{n}_xyz.csv"
            with xyz_path.open(encoding="utf-8-sig", newline="") as stream:
                xyz_rows = list(csv.DictReader(stream))
            if len(xyz_rows) != n*n or any(not all(math.isfinite(float(row[k])) for k in ("x", "y", "predicted_log10_TOF")) for row in xyz_rows):
                raise ValueError(f"XYZ validation failed: {xyz_path.name}")
            matrix_path = OUT / f"{prefix}_empirical_fullgrid_{n}x{n}_matrix.csv"
            with matrix_path.open(encoding="utf-8-sig", newline="") as stream:
                matrix = list(csv.reader(stream))
            if len(matrix) != n+1 or any(len(row) != n+1 for row in matrix):
                raise ValueError(f"Matrix dimensions invalid: {matrix_path.name}")
        with (OUT / f"{prefix}_catalyst_points.csv").open(encoding="utf-8-sig", newline="") as stream:
            points = list(csv.DictReader(stream))
        with (OUT / f"{prefix}_peak_summary.csv").open(encoding="utf-8-sig", newline="") as stream:
            peaks = list(csv.DictReader(stream))
        if len(points) != 13 or len({r["catalyst_name"] for r in points}) != 13 or len(peaks) != 2:
            raise ValueError(f"Point/peak table validation failed: {prefix}")
        for ext in ("png", "pdf"):
            image_path = OUT / f"{prefix}_TOF_empirical_map_full.{ext}"
            if not image_path.exists() or image_path.stat().st_size < 10_000:
                raise ValueError(f"Figure missing or suspiciously small: {image_path.name}")


def export_grid(prefix, descriptor, rows, coef, scales, bounds, resolution):
    xmin, xmax, ymin, ymax = bounds
    xvals = np.linspace(xmin, xmax, resolution)
    yvals = np.linspace(ymin, ymax, resolution)
    # x is the row coordinate and y the column coordinate, matching the
    # requested Origin matrix convention: first column x, top row y.
    XX, YY = np.meshgrid(xvals, yvals, indexing="ij")
    ZZ = evaluate(XX, YY, coef, scales)
    if ZZ.shape != (resolution, resolution) or not np.isfinite(ZZ).all():
        raise ValueError(f"Invalid full grid for {prefix}: shape={ZZ.shape}")

    xyz = [{"x": format(float(x), ".17g"), "y": format(float(y), ".17g"),
            "predicted_log10_TOF": format(float(z), ".17g")}
           for x, y, z in zip(XX.ravel(), YY.ravel(), ZZ.ravel())]
    write_csv(OUT / f"{prefix}_empirical_fullgrid_{resolution}x{resolution}_xyz.csv",
              ["x", "y", "predicted_log10_TOF"], xyz)

    matrix_path = OUT / f"{prefix}_empirical_fullgrid_{resolution}x{resolution}_matrix.csv"
    with matrix_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["x\\y", *[format(float(v), ".17g") for v in yvals]])
        for i, x in enumerate(xvals):
            writer.writerow([format(float(x), ".17g"), *[format(float(v), ".17g") for v in ZZ[i, :]]])
    return xvals, yvals, ZZ


def make_figure(prefix, descriptor, rows, coef, scales, bounds, peaks):
    xmin, xmax, ymin, ymax = bounds
    gx = np.linspace(xmin, xmax, 220)
    gy = np.linspace(ymin, ymax, 220)
    XX, YY = np.meshgrid(gx, gy, indexing="xy")
    ZZ = evaluate(XX, YY, coef, scales)
    fig, ax = plt.subplots(figsize=(10.2, 8.0), constrained_layout=True)
    surface = ax.contourf(XX, YY, ZZ, levels=60, cmap="viridis")
    colorbar = fig.colorbar(surface, ax=ax)
    colorbar.set_label("Predicted log10(TOF of C6H12O)")
    x = np.array([float(r["O"]) for r in rows])
    y = np.array([float(r[descriptor]) for r in rows])
    ax.scatter(x, y, s=60, marker="D", facecolor="white", edgecolor="#111827",
               linewidth=1.0, zorder=6, label="validated catalyst points (13)")
    label_points(ax, x, y, rows)
    ax.scatter([peaks[0]["x_peak"]], [peaks[0]["y_peak"]], marker="*", s=270,
               color="#ff3434", edgecolor="black", linewidth=1.0, zorder=10,
               label="full-range maximum")
    ax.scatter([peaks[1]["x_peak"]], [peaks[1]["y_peak"]], marker="X", s=130,
               color="#ffd166", edgecolor="black", linewidth=1.0, zorder=11,
               label="maximum within catalyst convex hull")
    full, hull = peaks
    if prefix == "O_Ha":
        full_text_pos, hull_text_pos = (0.02, 0.025), (0.30, 0.72)
    else:
        full_text_pos, hull_text_pos = (0.55, 0.055), (0.48, 0.88)
    ax.annotate(f"Full-range maximum (boundary)\nO={full['x_peak']:.4f} eV; {descriptor}={full['y_peak']:.4f} eV\n"
                f"pred. log10(TOF)={full['predicted_log10_TOF_peak']:.4f}",
                (full["x_peak"], full["y_peak"]), xytext=full_text_pos, textcoords="axes fraction",
                ha="left", va="bottom", fontsize=8,
                bbox={"boxstyle": "round,pad=.3", "fc": "white", "alpha": .94},
                arrowprops={"arrowstyle": "->", "lw": .7})
    ax.annotate(f"Convex-hull maximum\nO={hull['x_peak']:.4f} eV; {descriptor}={hull['y_peak']:.4f} eV\n"
                f"pred. log10(TOF)={hull['predicted_log10_TOF_peak']:.4f}",
                (hull["x_peak"], hull["y_peak"]), xytext=hull_text_pos, textcoords="axes fraction",
                ha="left", va="bottom", fontsize=8,
                bbox={"boxstyle": "round,pad=.3", "fc": "white", "alpha": .94},
                arrowprops={"arrowstyle": "->", "lw": .7})
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("O* historical formation energy (eV)")
    ax.set_ylabel(f"{descriptor}* historical formation energy (eV)")
    ax.set_title(f"O* + {descriptor}*: empirical quadratic response surface\n"
                 "Full rectangular evaluation; NOT a CatMAP continuous grid")
    ax.grid(alpha=.15)
    ax.legend(loc="best", fontsize=8)
    fig.savefig(OUT / f"{prefix}_TOF_empirical_map_full.png", dpi=300)
    fig.savefig(OUT / f"{prefix}_TOF_empirical_map_full.pdf")
    plt.close(fig)


def main():
    rows, metrics = read_frozen_inputs()
    point_xy = np.array([[float(r["O"]), float(r["Ha"])] for r in rows])
    hull_indices = ConvexHull(point_xy).vertices
    results = {}
    for descriptor, prefix in (("Ha", "O_Ha"), ("Hb", "O_Hb")):
        model = metrics["two_dimensional"][descriptor]
        coef = np.asarray(model["coefficients_standardized"], dtype=float)
        if coef.shape != (6,):
            raise ValueError(f"Unexpected saved quadratic coefficients for {descriptor}")
        xs = np.array([float(r["O"]) for r in rows])
        ys = np.array([float(r[descriptor]) for r in rows])
        scales = (float(xs.mean()), float(xs.std()), float(ys.mean()), float(ys.std()))
        design = np.column_stack([np.ones(13), (xs-scales[0])/scales[1], (ys-scales[2])/scales[3],
                                  ((xs-scales[0])/scales[1])**2,
                                  ((xs-scales[0])/scales[1])*((ys-scales[2])/scales[3]),
                                  ((ys-scales[2])/scales[3])**2])
        fitted = design @ coef
        actual = np.array([float(r["log10_TOF"]) for r in rows])
        r2 = 1 - float(np.sum((fitted-actual)**2) / np.sum((actual-actual.mean())**2))
        if abs(r2-float(model["in_sample_r2"])) > 1e-8:
            raise ValueError(f"Saved coefficients do not reproduce previous {descriptor} fit (R2 {r2})")

        def padded_range(values):
            span = float(values.max()-values.min())
            if span <= 0:
                raise ValueError("Descriptor range must be nonzero")
            return float(values.min()-.05*span), float(values.max()+.05*span)

        xmin, xmax = padded_range(xs)
        ymin, ymax = padded_range(ys)
        rectangle = np.array([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])
        catalyst_hull = point_xy[hull_indices] if descriptor == "Ha" else np.array([[float(r["O"]), float(r["Hb"])] for r in rows])[ConvexHull(np.array([[float(r["O"]), float(r["Hb"])] for r in rows])).vertices]
        full_peak = exact_maximum(coef, scales, rectangle)
        hull_peak = exact_maximum(coef, scales, catalyst_hull)
        peaks = [("global maximum in full rectangular plotting range", full_peak),
                 ("maximum within catalyst convex hull", hull_peak)]
        peak_rows = [{**peak, "peak_type": kind} for kind, peak in peaks]
        write_csv(OUT / f"{prefix}_peak_summary.csv",
                  ["x_peak", "y_peak", "predicted_log10_TOF_peak", "predicted_TOF_peak", "peak_type"], peak_rows)
        point_rows = [{"catalyst_name": r["catalyst"], "O_energy": r["O"],
                       f"{descriptor}_energy": r[descriptor], "actual_log10_TOF": r["log10_TOF"],
                       "actual_TOF": r["TOF"], "point_type": "validated catalyst"} for r in rows]
        write_csv(OUT / f"{prefix}_catalyst_points.csv", list(point_rows[0].keys()), point_rows)
        for resolution in (200, 400):
            export_grid(prefix, descriptor, rows, coef, scales, (xmin, xmax, ymin, ymax), resolution)
        make_figure(prefix, descriptor, rows, coef, scales, (xmin, xmax, ymin, ymax), [full_peak, hull_peak])
        results[descriptor] = {
            "bounds_eV": {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax},
            "resolution_points": [200, 400],
            "quadratic_in_sample_R2_reproduced": r2,
            "full_range_maximum": full_peak,
            "convex_hull_maximum": hull_peak,
            "convex_hull_maximum_is_boundary": point_on_boundary(catalyst_hull, (hull_peak["x_peak"], hull_peak["y_peak"])),
            "peaks_coincide": bool(np.allclose([full_peak["x_peak"], full_peak["y_peak"]],
                                                [hull_peak["x_peak"], hull_peak["y_peak"]], rtol=0, atol=1e-10)),
            "nearest_catalysts_to_full_range_maximum": sorted(
                ((math.hypot((float(r["O"])-full_peak["x_peak"])/np.std(xs),
                             (float(r[descriptor])-full_peak["y_peak"])/np.std(ys)), r["catalyst"]) for r in rows))[:3],
        }

    validate_exports()
    notes = [
        "# O*+H*经验二次响应面的全矩形网格导出",
        "",
        "## 方法与数据边界",
        "- 输入：本目录 `matched_data.csv`（13 个按 surface_name 匹配的已验证 CatMAP TOF 与历史 O/Ha/Hb 形成能）及 `fit_metrics.json` 中保存的原二次响应面系数。",
        "- 沿用原六参数标准化二次响应面；未重新拟合。脚本以原系数复算13点样本内R²并与既有指标核对。",
        "- 未运行 CatMAP，未更改 TOF、DFT/NEB、反应网络或旧图。",
        "- 补全方式：仅移除原图的凸包遮罩，在完整矩形上直接评估同一个二次多项式。矩形外推区只是经验模型数学延伸，不代表新催化剂结果。特别是若最大值落在扩展矩形边界且远离催化剂凸包，不得解释为物理活性峰。",
        "- 形成能为历史参考下的 O、Ha、Hb 描述符；Ha 与 Hb 分开。历史单物种形成能与当前 CatMAP 能量参考尚未完成物理一致性验证。",
        "- 新图标题明确说明为经验二次响应面、不是 CatMAP 连续网格。点为13个真实已验证催化剂；模型表面不是新的 CatMAP 计算。",
        "- 网格 XYZ 每种分辨率均为 N×N 行且无 NaN。矩阵 CSV 左上角为 `x\\y`，第一列为 x，首行其余列为 y，内部为预测 log10(TOF)。",
        "",
        "## 坐标范围与最大值",
        "",
        "| 图 | xmin (eV) | xmax (eV) | ymin (eV) | ymax (eV) | full-range maximum (x, y, log10TOF) | convex-hull maximum (x, y, log10TOF) |",
        "|---|---:|---:|---:|---:|---|---|",
    ]
    for species, label in (("Ha", "O*+Ha*"), ("Hb", "O*+Hb*")):
        item = results[species]
        b, f, h = item["bounds_eV"], item["full_range_maximum"], item["convex_hull_maximum"]
        notes.append(f"| {label} | {b['xmin']:.10f} | {b['xmax']:.10f} | {b['ymin']:.10f} | {b['ymax']:.10f} | ({f['x_peak']:.10f}, {f['y_peak']:.10f}, {f['predicted_log10_TOF_peak']:.10f}) | ({h['x_peak']:.10f}, {h['y_peak']:.10f}, {h['predicted_log10_TOF_peak']:.10f}) |")
        full_class = "boundary maximum" if any(abs(f[k]-b[lim]) < 1e-9 for k, lim in (("x_peak", "xmin"), ("x_peak", "xmax"), ("y_peak", "ymin"), ("y_peak", "ymax"))) else "interior maximum"
        hull_class = "boundary maximum" if item["convex_hull_maximum_is_boundary"] else "interior maximum"
        notes.append(f"- {label} full-range maximum: {full_class}; convex-hull maximum: {hull_class}. 两者{'重合' if item['peaks_coincide'] else '不重合'}。全域最大值附近催化剂（标准化描述符距离最近）：{', '.join(name for _, name in item['nearest_catalysts_to_full_range_maximum'])}。")
    notes.extend([
        "",
        "## Origin 文件",
        "- 精细绘图优先导入 `*_empirical_fullgrid_400x400_xyz.csv`：XYZ 长表最直接，三列为 x、y、predicted_log10_TOF。",
        "- 若使用 Origin Matrix/Contour，可导入 `*_empirical_fullgrid_400x400_matrix.csv`；200×200版本用于较轻量绘图。",
        "- 催化剂叠加点与峰值分别见 `*_catalyst_points.csv` 和 `*_peak_summary.csv`。峰值在连续解析二次面/凸包上的精确优化位置，不一定恰好落在离散网格节点。",
        "",
        "## 输出文件",
    ])
    for name in sorted(p.name for p in OUT.iterdir() if p.name.startswith(("O_Ha_", "O_Hb_")) and p.name not in {"O_Ha_TOF_empirical_map.png", "O_Hb_TOF_empirical_map.png", "O_Ha_TOF_empirical_map.pdf", "O_Hb_TOF_empirical_map.pdf"}):
        notes.append(f"- `{name}`")
    notes.extend(["", f"输入 matched_data.csv SHA256：`{sha256(DATA_PATH)}`", f"输入 fit_metrics.json SHA256：`{sha256(METRICS_PATH)}`", ""])
    (OUT / "FULLGRID_NOTES.md").write_text("\n".join(notes), encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
