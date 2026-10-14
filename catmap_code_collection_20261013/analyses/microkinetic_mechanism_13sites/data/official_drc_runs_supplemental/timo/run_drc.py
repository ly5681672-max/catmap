import csv
import json
from catmap import ReactionModel, __version__

if __version__ != "0.3.1":
    raise RuntimeError("Expected CatMAP 0.3.1, got " + str(__version__))
model = ReactionModel(setup_file="alkene_epoxidation.mkm")
model.run(recalculate=True)
labels = model.output_labels["rate_control"]
points = model.rate_control_map
if len(points) != 1:
    raise RuntimeError("Expected exactly one fixed-condition descriptor point")
matrix = points[0][1]
if len(matrix) != len(labels[0]) or any(len(row) != len(labels[1]) for row in matrix):
    raise RuntimeError("Official DRC dimensions do not match output_labels")
with open("rate_control.csv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(["gas_species", "energy_parameter", "degree_of_rate_control"])
    for gas, row in zip(labels[0], matrix):
        for parameter, value in zip(labels[1], row):
            writer.writerow([gas, parameter, repr(value)])
metadata = {
    "catmap_version": __version__,
    "surface": model.surface_names[0],
    "temperature_K": float(model.temperature),
    "pressure_mode": model.pressure_mode,
    "solver_tolerance": str(model.solver.tolerance),
    "output_labels": labels,
    "product_gas": "C6H12O_g",
    "product_DRC": matrix[labels[0].index("C6H12O_g")],
    "method": "CatMAP official get_rate_control; numerical Jacobian of turnover_frequency",
}
with open("drc_run_metadata.json", "w", encoding="utf-8") as handle:
    json.dump(metadata, handle, indent=2, ensure_ascii=False, default=str)
print(json.dumps(metadata, ensure_ascii=False, default=str))
