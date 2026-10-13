"""Re-render CatMAP 0.3.1 product maps with all cached catalyst points labeled.

The official 15x15 map cache is reused without ReactionModel.run(). The only
display changes are restoring the cached descriptor point dictionary, adding
official VectorMap point/annotation arguments, and padding the axes after the
native map is drawn so boundary markers and labels are not clipped.
"""

from __future__ import annotations

import ast
import os
from pathlib import Path

import numpy as np
from catmap import ReactionModel, __version__, plt
from catmap.analyze import VectorMap


OUT = Path(__file__).resolve().parent
RUN = OUT.parents[1]
SOURCE = RUN / "results" / "official_plots_round2"
SETUP = SOURCE / "alkene_epoxidation_descriptor_map.mkm"
LOG = SOURCE / "catmap_descriptor_grid.log"
SURFACES = [
    "tife", "timn", "tihf", "tire", "timo", "tiv", "tico",
    "titi", "tiw", "tita", "ticr", "ti", "tini",
]
DISPLAY_NAMES = {
    "tife": "Ti–Fe", "timn": "Ti–Mn", "tihf": "Ti–Hf",
    "tire": "Ti–Re", "timo": "Ti–Mo", "tiv": "Ti–V",
    "tico": "Ti–Co", "titi": "Ti–Ti", "tiw": "Ti–W",
    "tita": "Ti–Ta", "ticr": "Ti–Cr", "ti": "Ti",
    "tini": "Ti–Ni",
}
LABEL_OFFSETS = {
    "Ti–Fe": (7, 9, "left"),
    "Ti–Mn": (8, 12, "left"),
    "Ti–Hf": (8, 13, "left"),
    "Ti–Re": (-9, -15, "right"),
    "Ti–Mo": (-10, 12, "right"),
    "Ti–V": (10, 10, "left"),
    "Ti–Co": (0, 12, "center"),
    "Ti–Ti": (-11, -13, "right"),
    "Ti–W": (-20, 2, "right"),
    "Ti–Ta": (8, -13, "left"),
    "Ti–Cr": (10, -14, "left"),
    "Ti": (-8, 11, "right"),
    "Ti–Ni": (-2, -15, "center"),
}
PRODUCT = "C6H12O_g"


def log_literal(variable: str):
    prefix = variable + " = "
    for line in LOG.read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix):
            return ast.literal_eval(line[len(prefix):])
    raise RuntimeError(f"CatMAP log lacks {variable}")


def load_cached_model():
    old_cwd = Path.cwd()
    try:
        os.chdir(SOURCE)
        model = ReactionModel(setup_file=str(SETUP))
    finally:
        os.chdir(old_cwd)
    if __version__ != "0.3.1":
        raise RuntimeError(f"Expected CatMAP 0.3.1, got {__version__}")
    if list(model.surface_names) != SURFACES:
        raise RuntimeError("Cached catalyst list differs from the validated 13-site grid")
    if len(model.elementary_rxns) != 6:
        raise RuntimeError("Expected the unchanged six-step mechanism")
    if len(model.turnover_frequency_map) != 225 or len(model.production_rate_map) != 225:
        raise RuntimeError("Expected 225 cached descriptor-grid points")
    model.output_labels = log_literal("output_labels")
    raw_coords = log_literal("descriptor_dict")
    if list(raw_coords) != SURFACES:
        raise RuntimeError("Cached descriptor points do not match the validated site list")
    # Use display names as keys only in this in-memory plotting copy; each
    # coordinate is carried over exactly from CatMAP's saved run log.
    model.descriptor_dict = {
        DISPLAY_NAMES[surface]: raw_coords[surface] for surface in SURFACES
    }
    return model, raw_coords


def selected_values(model, variable: str) -> np.ndarray:
    label_index = model.output_labels[variable].index(PRODUCT)
    return np.asarray(
        [float(row[label_index]) for _, row in getattr(model, variable + "_map")],
        dtype=float,
    )


def render(model, variable: str, output_name: str) -> dict:
    plotter = VectorMap(model)
    plotter.plot_variable = variable
    plotter.plot_mode = "separate"
    plotter.include_labels = [PRODUCT]
    plotter.include_descriptors = True
    plotter.descriptor_labels = [
        "TS1_s formation energy (eV)",
        "EpoxideHHO_s formation energy (eV)",
    ]
    plotter.plot_function = "contourf"
    plotter.log_scale = True
    plotter.n_ticks = 10
    plotter.default_descriptor_pt_args = {
        "marker": "o", "linestyle": "None", "markersize": 4.2,
        "markerfacecolor": "black", "markeredgecolor": "white",
        "markeredgewidth": 0.55,
    }
    plotter.descriptor_pt_args = {
        label: dict(plotter.default_descriptor_pt_args)
        for label in DISPLAY_NAMES.values()
    }
    plotter.descriptor_label_args = {}
    for label, (dx, dy, align) in LABEL_OFFSETS.items():
        plotter.descriptor_label_args[label] = {
            "xytext": (dx, dy), "textcoords": "offset points",
            "fontsize": 7, "ha": align, "va": "center",
            "bbox": {
                "facecolor": "white", "edgecolor": "none",
                "alpha": 0.78, "pad": 0.55,
            },
        }

    values = selected_values(model, variable)
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size == 0:
        raise RuntimeError(f"No positive values in cached {variable} map")
    plotter.min, plotter.max = [float(v) for v in np.percentile(positive, [2, 98])]
    plotter.output_file = str(OUT / output_name)
    fig = plotter.plot(save=False)
    ax = fig.axes[0]

    observed = {text.get_text() for text in ax.texts}
    expected = set(DISPLAY_NAMES.values())
    if observed != expected:
        plt.close(fig)
        raise RuntimeError(
            f"Official descriptor annotations mismatch: missing={sorted(expected-observed)}, "
            f"unexpected={sorted(observed-expected)}"
        )

    # Expand only the displayed viewport after CatMAP has rendered the map at
    # its original ranges; this leaves the cached map values and contour domain
    # unchanged while preventing edge markers/labels from being cropped.
    x0, x1 = model.descriptor_ranges[0]
    y0, y1 = model.descriptor_ranges[1]
    padx, pady = (x1 - x0) * 0.06, (y1 - y0) * 0.06
    ax.set_xlim(x0 - padx, x1 + padx)
    ax.set_ylim(y0 - pady, y1 + pady)
    fig.savefig(plotter.output_file)
    plt.close(fig)
    return {
        "file": output_name,
        "labels": len(observed),
        "points": len(raw_coords := model.descriptor_dict),
        "color_limits": [plotter.min, plotter.max],
        "original_descriptor_ranges": model.descriptor_ranges,
        "display_ranges": [[x0 - padx, x1 + padx], [y0 - pady, y1 + pady]],
    }


if __name__ == "__main__":
    model, _ = load_cached_model()
    for variable, name in (
        ("turnover_frequency", "product_tof.pdf"),
        ("production_rate", "product_production_rate.pdf"),
    ):
        print(render(model, variable, name))
    print("CatMAP ReactionModel.run() was not called; only cached maps were plotted.")
