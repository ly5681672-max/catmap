"""Render default and readability-enhanced CatMAP 0.3.1 maps from cache.

No ReactionModel.run() call is made: the existing 15x15 grid and kinetic
outputs are loaded from the round2 CatMAP data file and only VectorMap display
settings are changed. Original PDFs and solver inputs remain untouched.
"""

from __future__ import annotations

import ast
import os
from pathlib import Path

import numpy as np
from catmap import ReactionModel, __version__, plt
from catmap.analyze import VectorMap


OUT_DIR = Path(__file__).resolve().parent
RUN_DIR = OUT_DIR.parents[1]
SOURCE_DIR = RUN_DIR / "results" / "official_plots_round2"
SETUP = SOURCE_DIR / "alkene_epoxidation_descriptor_map.mkm"
LOG = SOURCE_DIR / "catmap_descriptor_grid.log"
PRODUCT = "C6H12O_g"
EXPECTED_POINTS = 15 * 15


def labels_from_catmap_log() -> dict:
    # CatMAP's text log stores this small label map as a Python literal; parse
    # it without executing the log or loading any solver code from it.
    for line in LOG.read_text(encoding="utf-8").splitlines():
        if line.startswith("output_labels = "):
            return ast.literal_eval(line.split("=", 1)[1].strip())
    raise RuntimeError("CatMAP output_labels not found in the saved official log")


def values_for(model, variable: str, *, selected_label: str | None = None) -> np.ndarray:
    labels = model.output_labels[variable]
    if selected_label is None:
        return np.asarray(
            [float(value) for _, row in getattr(model, variable + "_map") for value in row],
            dtype=float,
        )
    index = labels.index(selected_label)
    return np.asarray(
        [float(row[index]) for _, row in getattr(model, variable + "_map")],
        dtype=float,
    )


def render(model, variable: str, out_path: Path, *, labels=None,
           minimum=None, maximum=None, log_scale=True, n_ticks=6) -> None:
    plotter = VectorMap(model)
    plotter.plot_variable = variable
    plotter.plot_mode = "separate"
    plotter.include_labels = labels
    plotter.include_descriptors = True
    plotter.descriptor_labels = [
        "TS1_s formation energy (eV)",
        "EpoxideHHO_s formation energy (eV)",
    ]
    plotter.default_descriptor_pt_args = {
        "marker": "o", "linestyle": "None", "markersize": 3.5, "color": "black",
    }
    plotter.default_descriptor_label_args = {
        "xytext": (3, 3), "textcoords": "offset points", "fontsize": 6,
    }
    plotter.plot_function = "contourf"
    plotter.log_scale = log_scale
    plotter.min = minimum
    plotter.max = maximum
    plotter.n_ticks = n_ticks
    plotter.output_file = str(out_path)
    figure = plotter.plot(save=True)
    plt.close(figure)


def main() -> None:
    if __version__ != "0.3.1":
        raise RuntimeError(f"Expected CatMAP 0.3.1; found {__version__}")
    if not SETUP.is_file() or not LOG.is_file():
        raise FileNotFoundError("Existing round2 CatMAP setup or log is missing")

    old_cwd = Path.cwd()
    try:
        os.chdir(SOURCE_DIR)
        model = ReactionModel(setup_file=str(SETUP))
    finally:
        os.chdir(old_cwd)

    model.output_labels = labels_from_catmap_log()
    if len(model.elementary_rxns) != 6:
        raise RuntimeError("Cached model does not contain the unchanged six-step mechanism")
    for map_name in ("turnover_frequency_map", "production_rate_map", "rate_map", "coverage_map"):
        if len(getattr(model, map_name, []) or []) != EXPECTED_POINTS:
            raise RuntimeError(f"Expected {EXPECTED_POINTS} cached points in {map_name}")

    default_dir = OUT_DIR / "default_official"
    enhanced_dir = OUT_DIR / "enhanced_official"
    default_dir.mkdir(exist_ok=True)
    enhanced_dir.mkdir(exist_ok=True)

    # Reproduce the existing round2 official VectorMap settings as the default
    # variant; only the output destination changes.
    specs = [
        ("turnover_frequency", "product_tof", [PRODUCT], PRODUCT, 1e-300),
        ("production_rate", "product_production_rate", [PRODUCT], PRODUCT, 1e-300),
        ("coverage", "coverage", None, None, None),
        ("rate", "six_step_net_rate", None, None, None),
    ]
    limits: dict[str, tuple[float, float]] = {}
    for variable, stem, selected_labels, quantile_label, default_minimum in specs:
        render(
            model, variable, default_dir / f"{stem}.pdf", labels=selected_labels,
            minimum=default_minimum, maximum=None, log_scale=True, n_ticks=6,
        )

        vals = values_for(model, variable, selected_label=quantile_label)
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if vals.size == 0:
            raise RuntimeError(f"No positive cached values available for {variable}")
        lo, hi = np.percentile(vals, [2, 98])
        if not lo < hi:
            raise RuntimeError(f"Degenerate percentile color limits for {variable}: {lo}, {hi}")
        limits[variable] = (float(lo), float(hi))
        render(
            model, variable, enhanced_dir / f"{stem}.pdf", labels=selected_labels,
            minimum=float(lo), maximum=float(hi), log_scale=True, n_ticks=10,
        )

    print(f"CatMAP {__version__}: rendered cached maps only; solver was not run.")
    print(f"Default maps: {default_dir}")
    print(f"Enhanced maps: {enhanced_dir}")
    for variable, (lo, hi) in limits.items():
        print(f"{variable} enhanced color limits (positive-value p2-p98): {lo:.8g} .. {hi:.8g}")
    print("Ti-Zr remains unverified and is excluded from this 13-site cached grid.")


if __name__ == "__main__":
    main()
