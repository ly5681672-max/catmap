"""CatMAP 0.3.1 parse and thermodynamic validation for the isolated candidate."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, required=True)
    args = parser.parse_args()
    model_path = args.model.resolve()
    run_dir = model_path.parent
    old_cwd = Path.cwd()
    os.chdir(run_dir)
    try:
        from catmap import ReactionModel

        model = ReactionModel(setup_file=model_path.name)
        model.run()
        if len(model.elementary_rxns) != 9:
            raise AssertionError(f"expected 9 elementary reactions, got {len(model.elementary_rxns)}")
        if "OOH_s" not in model.adsorbate_names:
            raise AssertionError("explicit stable OOH* is missing")

        # Two non-grid descriptor coordinates are used only to request the
        # single Ti-Co thermodynamic state; all energies come from TableParser.
        descriptor_values = [model.temperature, 1.0]
        energies = model.scaler.get_free_energies(descriptor_values)
        rows = []
        for index, reaction in enumerate(model.elementary_rxns, start=1):
            delta_g, forward_barrier = model.get_rxn_energy(reaction, energies)
            initial = model.get_state_energy(reaction[0], energies)
            final = model.get_state_energy(reaction[-1], energies)
            if len(reaction) == 3:
                ts_energy = model.get_state_energy(reaction[1], energies)
                reverse_barrier = max(0, ts_energy - final, -delta_g)
                ts_label = "+".join(reaction[1])
                barrier_type = "explicit NEB TS" if ts_label in {"NEB1TS_s", "NEB2TS_s"} else "explicit candidate TS"
            else:
                ts_energy = max(initial, final)
                reverse_barrier = max(0, ts_energy - final, -delta_g)
                ts_label = "CatMAP endpoint max"
                barrier_type = "CatMAP endpoint approximation; no DFT TS"
            detailed_balance_error = (forward_barrier - reverse_barrier) - delta_g
            if abs(detailed_balance_error) > 1e-10:
                raise AssertionError(f"detailed-balance barrier identity failed in reaction {index}: {detailed_balance_error}")
            rows.append({
                "step": str(index),
                "reaction_internal": model.print_rxn(reaction, mode="text"),
                "TS_source": ts_label,
                "barrier_status": barrier_type,
                "deltaG_eV": f"{float(delta_g):.12g}",
                "forward_barrier_eV": f"{float(forward_barrier):.12g}",
                "reverse_barrier_eV": f"{float(reverse_barrier):.12g}",
                "detailed_balance_error_eV": f"{float(detailed_balance_error):.4g}",
            })

        out_csv = run_dir / "candidate_reaction_energy_barriers.csv"
        with out_csv.open("w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

        species_rows = []
        for species in model.adsorbate_names + model.transition_state_names:
            definition = model.species_definitions[species]
            species_rows.append({
                "species": species,
                "class": definition.get("type", ""),
                "composition": str(definition.get("composition", {})),
                "n_sites": str(definition.get("n_sites", "")),
                "formation_energy_eV": str(energies.get(species)),
            })
        with (run_dir / "candidate_species_audit.csv").open("w", newline="", encoding="utf-8-sig") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(species_rows[0]))
            writer.writeheader()
            writer.writerows(species_rows)

        import catmap
        print(f"CatMAP={getattr(catmap, '__version__', 'version unavailable')}")
        print(f"elementary_reactions={len(model.elementary_rxns)}")
        print(f"stable_adsorbates={len(model.adsorbate_names)}: {','.join(model.adsorbate_names)}")
        print(f"transition_states={len(model.transition_state_names)}: {','.join(model.transition_state_names)}")
        print(f"atom_and_site_balance=PASS (ReactionModel setup checks)")
        print(f"detailed_balance=PASS for {len(rows)} reactions")
        print(f"reaction_table={out_csv}")
        print(f"species_table={run_dir / 'candidate_species_audit.csv'}")
    finally:
        os.chdir(old_cwd)


if __name__ == "__main__":
    main()
