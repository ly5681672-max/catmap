"""Try every relevant CatMAP 0.3.1 native plot on the cached round2 model.

This script never calls ReactionModel.run(); it only loads the existing mapper
outputs, invokes CatMAP analyze classes, and copies previously generated
CatMAP-native PDFs into this review directory for a single complete inventory.
"""

from __future__ import annotations

import ast
import os
import shutil
import traceback
from pathlib import Path

from catmap import ReactionModel, __version__, plt
from catmap.analyze import MatrixMap, MechanismAnalysis, ScalingAnalysis, VectorMap


HERE = Path(__file__).resolve().parent
RUN = HERE.parents[1]
RESULTS = RUN / "results"
SOURCE = RESULTS / "official_plots_round2"
SETUP = SOURCE / "alkene_epoxidation_descriptor_map.mkm"
LOG = SOURCE / "catmap_descriptor_grid.log"
PRODUCT = "C6H12O_g"


def saved_labels() -> dict:
    for line in LOG.read_text(encoding="utf-8").splitlines():
        if line.startswith("output_labels = "):
            return ast.literal_eval(line.split("=", 1)[1].strip())
    raise RuntimeError("CatMAP saved output_labels are missing")


def report_attempt(log, name: str, fn) -> None:
    try:
        result = fn()
        log.write(f"PASS\t{name}\t{result or 'generated'}\n")
    except Exception as exc:  # retain exact native-interface failure for audit
        log.write(f"UNAVAILABLE\t{name}\t{type(exc).__name__}: {exc}\n")
        log.write(traceback.format_exc() + "\n")


def vector_plot(model, variable: str, target: Path, *, log_scale=True) -> str:
    plotter = VectorMap(model)
    plotter.plot_variable = variable
    if variable == "selectivity":
        # Matplotlib's contour routines do not accept CatMAP's cached mpmath
        # scalar objects. Convert only the temporary plotting copy to float;
        # the saved CatMAP cache and the underlying values remain untouched.
        plotter._rxm.selectivity_map = [
            (point, [float(value) for value in row])
            for point, row in model.selectivity_map
        ]
    plotter.plot_mode = "separate"
    plotter.include_labels = [PRODUCT] if variable in {
        "turnover_frequency", "production_rate", "selectivity",
        "consumption_rate",
    } else None
    plotter.include_descriptors = True
    plotter.descriptor_labels = [
        "TS1_s formation energy (eV)",
        "EpoxideHHO_s formation energy (eV)",
    ]
    plotter.plot_function = "contourf"
    plotter.log_scale = log_scale
    plotter.n_ticks = 10
    plotter.output_file = str(target)
    fig = plotter.plot(save=True)
    plt.close(fig)
    return str(target)


def main() -> None:
    if __version__ != "0.3.1":
        raise RuntimeError(f"Expected CatMAP 0.3.1, found {__version__}")
    if not SETUP.is_file() or not LOG.is_file():
        raise FileNotFoundError("The existing CatMAP descriptor setup/cache is missing")

    copied = HERE / "previously_generated_official"
    copied.mkdir(exist_ok=True)
    # Package earlier official-native figures unchanged; preserve source layout.
    for name in ("official_plots_round2", "official_plots_round3"):
        src = RESULTS / name
        if not src.exists():
            continue
        dst = copied / name
        if name.endswith("round2"):
            files = list(src.glob("*.pdf"))
            dst.mkdir(exist_ok=True)
            for file in files:
                shutil.copy2(file, dst / file.name)
        else:
            for sub in ("default_official", "enhanced_official"):
                sub_src = src / sub
                if sub_src.exists():
                    sub_dst = dst / sub
                    sub_dst.mkdir(parents=True, exist_ok=True)
                    for file in sub_src.glob("*.pdf"):
                        shutil.copy2(file, sub_dst / file.name)

    old_cwd = Path.cwd()
    try:
        os.chdir(SOURCE)
        model = ReactionModel(setup_file=str(SETUP))
    finally:
        os.chdir(old_cwd)
    model.output_labels = saved_labels()

    log_path = HERE / "official_plot_attempts.tsv"
    with log_path.open("w", encoding="utf-8", newline="") as log:
        log.write(f"CatMAP\t{__version__}\n")
        log.write("RUN\tReactionModel.run was not called; cached maps only\n")
        log.write("CACHED_MAPS\t" + ",".join(
            name.removesuffix("_map") for name in (
                "coverage_map", "rate_map", "turnover_frequency_map",
                "production_rate_map", "selectivity_map",
            ) if getattr(model, name, None)
        ) + "\n")

        # This output exists in the saved grid but its earlier figure was
        # skipped before invoking VectorMap. Try both native log and linear modes.
        report_attempt(log, "selectivity VectorMap/log", lambda: vector_plot(
            model, "selectivity", HERE / "selectivity_log_official.pdf",
            log_scale=True,
        ))
        report_attempt(log, "selectivity VectorMap/linear", lambda: vector_plot(
            model, "selectivity", HERE / "selectivity_linear_official.pdf",
            log_scale=False,
        ))

        # Native solver outputs not requested in the cached model: invoke the
        # plotting class to record the actual missing-map result, without solving.
        for variable in ("forward_rate", "reverse_rate", "consumption_rate"):
            report_attempt(log, f"VectorMap/{variable}", lambda v=variable: vector_plot(
                model, v, HERE / f"{v}_official.pdf", log_scale=True,
            ))

        for variable in ("rate_control", "selectivity_control", "rxn_order"):
            def matrix_attempt(v=variable):
                plotter = MatrixMap(model)
                plotter.plot_variable = v
                plotter.plot_mode = "separate"
                plotter.output_file = str(HERE / f"matrix_{v}_official.pdf")
                fig = plotter.plot(save=True)
                plt.close(fig)
                return plotter.output_file
            report_attempt(log, f"MatrixMap/{variable}", matrix_attempt)

        def mechanism_graph():
            analysis = MechanismAnalysis(model)
            analysis.rxn_mechanisms = {"forward_six_step": [0, 1, 2, 3, 4, 5]}
            out = HERE / "six_step_mechanism_graph"
            analysis.create_graph(
                mechanism="forward_six_step", filename=str(out),
                exclude_sites=True, exclude_ts=False,
            )
            return str(out)
        report_attempt(log, "MechanismAnalysis.create_graph", mechanism_graph)

        def scaling_plot():
            analysis = ScalingAnalysis(model)
            out = HERE / "official_scaling_analysis.pdf"
            fig = analysis.plot(save=str(out))
            plt.close(fig)
            return str(out)
        report_attempt(log, "ScalingAnalysis.plot", scaling_plot)

    print(f"CatMAP {__version__}; no kinetic calculation run")
    print(f"Attempt log: {log_path}")
    print(f"Output directory: {HERE}")


if __name__ == "__main__":
    main()
