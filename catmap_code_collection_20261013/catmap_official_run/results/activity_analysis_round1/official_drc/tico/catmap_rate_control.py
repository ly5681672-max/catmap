import csv
from catmap import ReactionModel, __version__

if __version__ != '0.3.1':
    raise RuntimeError('Expected CatMAP 0.3.1')
model = ReactionModel(setup_file='alkene_epoxidation.mkm')
model.run(recalculate=True)
labels = model.output_labels['rate_control']
points = model.rate_control_map
if len(points) != 1:
    raise RuntimeError('Expected one fixed-condition DRC point')
matrix = points[0][1]
if len(matrix) != len(labels[0]):
    raise RuntimeError('DRC gas dimension mismatch')
if any(len(row) != len(labels[1]) for row in matrix):
    raise RuntimeError('DRC parameter dimension mismatch')
with open('rate_control.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gas_species', 'energy_parameter', 'degree_of_rate_control'])
    for gas, row in zip(labels[0], matrix):
        for parameter, value in zip(labels[1], row):
            w.writerow([gas, parameter, repr(value)])
print('CatMAP', __version__)
print('surface', model.surface_names[0], 'T', model.temperature,
      'pressure_mode', model.pressure_mode)
print('output_labels', labels)
print('product_TOF_DRC', matrix[labels[0].index('C6H12O_g')])
print('official_rate_control_csv=rate_control.csv')
