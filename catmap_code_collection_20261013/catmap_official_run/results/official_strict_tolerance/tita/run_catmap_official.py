"""Invoke the unmodified CatMAP 0.3.1 ReactionModel on this isolated input."""

import csv
from catmap import ReactionModel, __version__
from catmap.analyze import MechanismAnalysis


if __version__ != '0.3.1':
    raise RuntimeError('Expected CatMAP 0.3.1; found {}'.format(__version__))

print('CatMAP {}'.format(__version__))
model = ReactionModel(setup_file='alkene_epoxidation.mkm')
print('surface_names={!r}'.format(model.surface_names))
print('temperature={!r}'.format(model.temperature))
print('numerical_representation={!r}; decimal_precision={!r}; use_numbers_solver={!r}'.format(
    model.numerical_representation, model.decimal_precision, model.use_numbers_solver
))
model.run(recalculate=True)
print('ReactionModel.run() completed')
print('output_variables={!r}'.format(model.output_variables))
print('output_labels={!r}'.format(model.output_labels))


def single_output(name):
    points = getattr(model, name + '_map')
    if len(points) != 1:
        raise RuntimeError('{} contains {} descriptor points, expected one'.format(name, len(points)))
    return points[0][1]


coverages = single_output('coverage')
net_rates = single_output('rate')
tofs = single_output('turnover_frequency')
directional_rates = single_output('directional_rates')
residual = model.solver.get_residual(coverages)
site_names = [name for name in model.site_names if name != 'g']
coverage_names = list(model.adsorbate_names) + site_names
coverage_sum = model._math.fsum(coverages)
if len(coverages) != len(coverage_names):
    raise RuntimeError('CatMAP coverage vector and official labels have different lengths')
if residual > model.solver.tolerance:
    raise RuntimeError('Official CatMAP residual exceeds tolerance: {}'.format(residual))
if any(value < 0 for value in coverages):
    raise RuntimeError('CatMAP returned a negative coverage')
if abs(coverage_sum - 1) > model._mpfloat('1e-30'):
    raise RuntimeError('CatMAP coverages do not sum to one: {}'.format(coverage_sum))

with open('coverage.csv', 'w', newline='', encoding='utf-8') as handle:
    writer = csv.writer(handle)
    writer.writerow(['species', 'coverage'])
    writer.writerows(zip(coverage_names, map(str, coverages)))
with open('net_reaction_rates.csv', 'w', newline='', encoding='utf-8') as handle:
    writer = csv.writer(handle)
    writer.writerow(['step', 'reaction', 'net_rate'])
    writer.writerows((i, repr(rxn), str(rate)) for i, (rxn, rate) in enumerate(
        zip(model.elementary_rxns, net_rates), start=1
    ))
with open('directional_reaction_rates.csv', 'w', newline='', encoding='utf-8') as handle:
    writer = csv.writer(handle)
    writer.writerow(['direction', 'reaction', 'rate'])
    labels = model.output_labels['directional_rates']
    writer.writerows((label, str(value)) for label, value in zip(labels, directional_rates))
with open('gas_turnover_frequencies.csv', 'w', newline='', encoding='utf-8') as handle:
    writer = csv.writer(handle)
    writer.writerow(['gas_species', 'turnover_frequency'])
    writer.writerows(zip(model.gas_names, map(str, tofs)))
forward_analysis = MechanismAnalysis(model)
forward_analysis.rxn_mechanisms = {'forward_six_step': [1, 2, 3, 4, 5, 6]}
forward_analysis.output_file = 'six_step_free_energy.pdf'
forward_analysis.plot(
    plot_variants=[[333.15, 1.0]],
    save=(model.surface_names[0] == 'tire'),
)
forward_deltas, forward_barriers = forward_analysis.data_dict['forward_six_step']
reverse_analysis = MechanismAnalysis(model)
reverse_analysis.rxn_mechanisms = {'reverse_six_step': [-6, -5, -4, -3, -2, -1]}
reverse_analysis.plot(plot_variants=[[333.15, 1.0]], save=False)
_, reverse_path_barriers = reverse_analysis.data_dict['reverse_six_step']
with open('effective_barriers.csv', 'w', newline='', encoding='utf-8') as handle:
    writer = csv.writer(handle)
    writer.writerow(['step', 'reaction', 'barrier_model', 'delta_G',
                     'forward_effective_barrier', 'reverse_effective_barrier'])
    for i, rxn in enumerate(model.elementary_rxns, start=1):
        writer.writerow([i, repr(rxn), 'explicit_TS' if len(rxn) == 3 else 'CatMAP_endpoint',
                         str(forward_deltas[i]), str(forward_barriers[i - 1]),
                         str(reverse_path_barriers[6 - i])])
with open('convergence.csv', 'w', newline='', encoding='utf-8') as handle:
    writer = csv.writer(handle)
    writer.writerow(['catmap_version', 'surface', 'temperature_K', 'pressure_mode',
                     'numerical_representation', 'decimal_precision', 'residual',
                     'tolerance', 'coverage_sum', 'status'])
    writer.writerow([__version__, model.surface_names[0], model.temperature,
                     model.pressure_mode, model.numerical_representation,
                     model.decimal_precision, str(residual), str(model.solver.tolerance),
                     str(coverage_sum), 'CONVERGED'])
print('official_steady_state_residual={!r}'.format(residual))
print('official_solver_tolerance={!r}'.format(model.solver.tolerance))
print('coverage_names={!r}'.format(coverage_names))
print('coverage_sum={!r}'.format(coverage_sum))
print('net_rates={!r}'.format(net_rates))
print('turnover_frequency_labels={!r}'.format(model.gas_names))
print('turnover_frequencies={!r}'.format(tofs))
print('pressure_corrected_CatMAP_MechanismAnalysis_barriers=effective_barriers.csv')

if model.surface_names and model.surface_names[0] == 'tire':
    from catmap.analyze import VectorMap

    print('official_MechanismAnalysis=six_step_free_energy.pdf')

    for variable, filename in (
        ('coverage', 'coverage_map.pdf'),
        ('rate', 'net_rate_map.pdf'),
    ):
        vector = VectorMap(model)
        vector.plot_variable = variable
        vector.output_file = filename
        try:
            vector.plot(save=True)
            print('official_VectorMap_{}={}'.format(variable, filename))
        except Exception as error:
            print('official_VectorMap_{}_unsupported={}: {}'.format(
                variable, type(error).__name__, error
            ))
