"""Generate isolated CatMAP 0.3.1 native mechanism and descriptor-map plots."""

from pathlib import Path
from math import log10
import sys

from catmap import ReactionModel, __version__, plt
from catmap.analyze import MechanismAnalysis, VectorMap


ROOT = Path(__file__).resolve().parent
SETUP = ROOT / 'alkene_epoxidation_descriptor_map.mkm'
SURFACES = [
    'tife', 'timn', 'tihf', 'tire', 'timo', 'tiv', 'tico',
    'titi', 'tiw', 'tita', 'ticr', 'ti', 'tini',
]
PRODUCT = 'C6H12O_g'
REACTIONS = [1, 2, 3, 4, 5, 6]


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    if __version__ != '0.3.1':
        raise RuntimeError('Expected CatMAP 0.3.1; found {}'.format(__version__))

    model = ReactionModel(setup_file=str(SETUP))
    if list(model.surface_names) != SURFACES:
        raise RuntimeError('Parsed surface list differs from the isolated setup')
    if len(model.elementary_rxns) != 6:
        raise RuntimeError('Expected the unchanged six-step reaction network')

    # This is a new, isolated official descriptor-grid calculation. It never
    # writes to the original single-point model or its prior result directories.
    model.run(recalculate=True)
    expected_points = 15 * 15
    actual_points = len(model.rate_map)
    print('CatMAP {}; descriptor points: {} / {}'.format(
        __version__, actual_points, expected_points
    ))

    # Free-energy diagrams use the same official MechanismAnalysis interface;
    # avoid_scaling=True restores the tabulated energies at each catalyst point.
    for surface in SURFACES:
        analysis = MechanismAnalysis(model)
        analysis.rxn_mechanisms = {'forward_six_step': REACTIONS}
        analysis.output_file = str(ROOT / 'free_energy_{}.pdf'.format(surface))
        analysis.plot(plot_variants=[surface], save=True)
        plt.close('all')

    def map_plot(variable, filename, *, labels=None, log_scale=False,
                 minimum=None, maximum=None):
        plotter = VectorMap(model)
        plotter.plot_variable = variable
        plotter.plot_mode = 'separate'
        plotter.include_labels = labels
        plotter.include_descriptors = True
        plotter.descriptor_labels = [
            'TS1_s formation energy (eV)',
            'EpoxideHHO_s formation energy (eV)',
        ]
        plotter.default_descriptor_pt_args = {
            'marker': 'o', 'linestyle': 'None', 'markersize': 3.5,
            'color': 'black',
        }
        plotter.default_descriptor_label_args = {
            'xytext': (3, 3), 'textcoords': 'offset points', 'fontsize': 6,
        }
        plotter.log_scale = log_scale
        plotter.min = minimum
        plotter.max = maximum
        plotter.output_file = str(ROOT / filename)
        figure = plotter.plot(save=True)
        plt.close(figure)
        print('official {} plot: {}'.format(variable, filename))

    map_plot('turnover_frequency', 'product_tof_descriptor_map.pdf',
             labels=[PRODUCT], log_scale=True, minimum=1e-300)
    map_plot('production_rate', 'product_production_rate_descriptor_map.pdf',
             labels=[PRODUCT], log_scale=True, minimum=1e-300)
    map_plot('coverage', 'coverage_descriptor_maps.pdf', log_scale=True)
    map_plot('rate', 'six_step_net_rate_descriptor_maps.pdf', log_scale=True)

    selectivity_index = model.output_labels['selectivity'].index(PRODUCT)
    selectivities = [float(values[selectivity_index])
                     for _, values in model.selectivity_map]
    if min(selectivities) > 0 and int(log10(min(selectivities))) == int(
            log10(max(selectivities))):
        print('official selectivity map skipped: CatMAP 0.3.1 log contour '
              'levels degenerate over this narrow range ({:.9g}..{:.9g})'.format(
                  min(selectivities), max(selectivities)
              ))
    else:
        map_plot('selectivity', 'selectivity_descriptor_map.pdf',
                 labels=[PRODUCT], log_scale=True)

    print('Map labels: {}'.format(model.output_labels))
    print('Map point counts: coverage={}, rate={}, TOF={}, production={}, selectivity={}'.format(
        len(model.coverage_map), len(model.rate_map),
        len(model.turnover_frequency_map), len(model.production_rate_map),
        len(model.selectivity_map),
    ))


if __name__ == '__main__':
    main()
