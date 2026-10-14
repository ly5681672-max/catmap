"""Read-only compatibility and avoid_scaling probe for CatMAP 0.3.1."""
import json
from pathlib import Path

from catmap import ReactionModel, __version__


if __version__ != "0.3.1":
    raise RuntimeError(f"Expected CatMAP 0.3.1, got {__version__}")
model = ReactionModel(setup_file="ts1_ts2_descriptor_grid.mkm")
model.parse()
model.verify()
scaler = model.scaler
coefficient_matrix = scaler.get_coefficient_matrix()
state_names = list(scaler.adsorbate_names) + list(scaler.transition_state_names)
if not state_names:
    raise RuntimeError("No adsorbate or transition-state species available for anchor checks")
descriptor_pair = scaler.descriptor_names
anchor_results = []
for surface in model.surface_names:
    coords = [float(x) for x in scaler.descriptor_dict[surface]]
    values = scaler.get_electronic_energies(coords)
    surface_id = model.surface_names.index(surface)
    per_state = {}
    for species in state_names:
        raw_value = scaler.parameter_dict[species][surface_id]
        if raw_value is None:
            continue
        raw = float(raw_value)
        returned = float(values[species])
        if abs(returned - raw) > 1e-10:
            raise RuntimeError(
                f"avoid_scaling did not return raw anchor {surface}/{species}: "
                f"{returned} != {raw}"
            )
        per_state[species] = {"raw_eV": raw, "returned_eV": returned,
                              "absolute_error_eV": abs(returned - raw)}
    if not per_state:
        raise RuntimeError(f"No raw surface-state energies available for {surface}")
    anchor_results.append({"surface_name": surface, "descriptor_1": coords[0],
                           "descriptor_2": coords[1], "state_checks": per_state,
                           "is_raw_anchor": True})

mid = [sum(bounds) / 2 for bounds in model.descriptor_ranges]
near = [mid[0] + 0.0001234567, mid[1] + 0.0002345678]
near_values = scaler.get_electronic_energies(near)
non_anchor_checks = {}
for species in state_names:
    # Descriptor states are set directly to the requested coordinate at the
    # end of get_electronic_energies(); they are not dependent scaler outputs.
    if species not in scaler.parameter_names or species in descriptor_pair:
        continue
    row = scaler.parameter_names.index(species)
    predicted = float(sum(float(coefficient_matrix[row, j]) * (near + [1])[j]
                          for j in range(len(near) + 1)))
    near_returned = float(near_values[species])
    if abs(near_returned - predicted) > 1e-8:
        raise RuntimeError(
            f"non-anchor state is not scaled as expected for {species}: "
            f"{near_returned} != {predicted}"
        )
    non_anchor_checks[species] = {"scaled_prediction_eV": predicted,
                                  "returned_eV": near_returned,
                                  "absolute_error_eV": abs(near_returned - predicted)}
if not non_anchor_checks:
    raise RuntimeError("No non-anchor scaling checks were completed")

report = {
    "catmap_version": __version__, "surface_count": len(model.surface_names),
    "surface_names": list(model.surface_names), "descriptor_names": list(descriptor_pair),
    "descriptor_ranges": model.descriptor_ranges, "resolution": model.resolution,
    "scaler": type(scaler).__name__, "solver": type(model.solver).__name__,
    "mapper": type(model.mapper).__name__, "coefficient_matrix_shape": list(coefficient_matrix.shape),
    "avoid_scaling": scaler.avoid_scaling,
    "descriptor_decimal_precision": scaler.descriptor_decimal_precision,
    "all_raw_anchor_checks_passed": True,
    "raw_anchor_state_check_count": sum(len(row["state_checks"]) for row in anchor_results),
    "non_anchor_scaled_check_passed": True,
    "non_anchor_state_checks": non_anchor_checks,
    "anchor_checks": anchor_results,
}
Path("avoid_scaling_probe.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
