# lyq_alkene_epoxidation1 — NEB-consistent epoxidation reference

**Isolated experiment, not a replacement for `lyq_alkene_epoxidation`.**
Created on branch `feature/lyq-alkene-epoxidation1-neb`; do **not** merge until the mechanistic assumptions are reviewed.

## Goal and provenance
- Source: `../lyq_alkene_epoxidation/己烯环氧化.xlsx`, worksheet **吉布斯自由能汇总**.
- `neb_three_state_summary.csv`: 15 surfaces, **G** relative to the Excel coadsorption initial state. TiW has full precision; most other entries are rounded 0.001 eV from the previously extracted worksheet summary. This is **not** the complete CI-NEB image series.
- Never alter the source Excel electronic energies (`E0`), free-energy corrections (`ΔG`), or corrected `G`.
- Do **not** mix the source's corrected free energies with the old model's formation energies unless the reference convention and correction levels are reconciled.

## State mapping
| Workbook state | Formula | Mapping |
|---|---|---|
| IS | coadsorbed (H2O2)(C6H12)* | `R*` (full C6H14O2) |
| labelled TS | complete C6H14O2 intermediate NEB image | pathway marker, **not automatically a CatMAP saddle point** |
| FS | coadsorbed (C6H12O)(H2O)* | `P*` (full C6H14O2) |

The old step `OOH* + C6H12* -> C6H12O* + O* + Hb*` is **not** this net coadsorbed reaction. In particular it omits the `Ha*` atom (handled elsewhere in the old network), so substituting its pseudo-TS species with this full-composition NEB marker is physically inconsistent.

## Files
- `neb_three_state_summary.csv`: read-only snapshot of the three reported points.
- `check_neb.py`: standard-library validation of state ordering and the old model's mismatched step-4 final state. It **never writes energies.txt**.
- `coadsorption_model.mkm.template`: a separate **CatMAP model template**, not executable as-is until coadsorption and gas chemical potential input is provided.
- `README.md`: assumptions, exact needed next inputs and usage.

## Generated Gibbs energies and provisional CatMAP model
The branch includes `energies.txt`, `energy_audit.csv`, `build_coadsorption.py`, and `coadsorption_model.mkm`. They use gas G (rows 57–60), catalyst slab G (row 4), coadsorbed IS G (row 12), path marker G (row 16), and coadsorbed FS G (row 20; TiMo row 21 combined total). TiTa has an extra address column, so its G column is 87 rather than 86. Element reference convention C=-9.28, H=-1.11, O=-4.37 eV, identical to the earlier CatMAP input convention. For surfaces with a downhill path marker, `RP` is set to `max(G_IS,G_marker,G_FS)` for a nonnegative model barrier; **this is a provisional kinetic assumption, not a measured NEB saddle**. The gas reservoir activities continue the old model approximation and require sensitivity/solvation validation. Running CatMAP to convergence has NOT been verified in this environment.

## Why no immediately validated CatMAP kinetics?
CatMAP needs thermodynamically consistent **absolute formation free energies** for `R*`, `P*`, gases and empty sites on *the same reference scale*. The three relative energies alone determine `G_P-G_R`, **not** `G_R-(G_H2O2,g+G_C6H12,g)` or product desorption. Choosing an arbitrary adsorption energy would silently invent kinetics. The workbook gives free-energy-corrected values but the complete reference-state crosswalk for every surface and gas correction must be verified first.

In a one-site minimal network, after reference alignment:
```
* + H2O2_g + C6H12_g <-> R*
R* -> P*                  # downhill barrierless elementary conversion only where NEB validates it
P* <-> C6H12O_g + H2O_g + *
```
For barriered surfaces, the true NEB saddle energy is needed; **do not treat a single labelled midpoint as TS without checking all images**. An irreversible first-order step alone also cannot capture liquid-phase mass transfer.

## Run checks
```powershell
python lyq_alkene_epoxidation1/check_neb.py
```
This reports all surfaces and compares the old model's step-4 final state, explicitly marking different chemical states. No files in the original directory are changed.

## What data is still required for defensible TOF?
1. Full NEB image sequence for each surface, with energy correction convention specified.
2. Clean-slab and gas free energies in a consistent convention (same temperature, corrections and reference).
3. For every surface, verified absolute formation free energies of coadsorbed `R*` and `P*`.
4. Confirm adsorbed site occupancy, adsorption/desorption mechanisms and kinetic prefactors.
