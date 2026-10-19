# Eight-step explicit OOH* network — conversion note (2026-10-19)

**Status: a CatMAP reaction-expression/configuration draft, not a converged kinetic result.** Do not compare or replace the existing six-step baseline/TOF until all energy inputs and the physical site model are validated.

## Model location and objective

- [explicit_ooh_8step.mkm](explicit_ooh_8step.mkm): independent candidate setup.
- Intended analyses: the **steady-state OOH* coverage vs TOF** and **common-reference OOH* adsorption thermodynamics vs TOF**.
- The previous six-step scheme retains `Rpair*` as a lumped effective state; the eight-step scheme follows the researcher's original elementary steps instead. They are *different reaction models* and must not be expected to reproduce the same TOF without further physical parameterization.

## Stoichiometry and formal CatMAP site accounting

| Step | Substrates / intermediates | Atom balance | Effective site balance |
|---|---|---|---|
| 1 | `H2O2_g + * <-> H2O2*` | H2O2 | 1 → 1 |
| 2 | `H2O2* + * -> Ha* + OOH*` | H2O2 | 2 → 2 |
| 3 | `C6H12_g + * <-> C6H12*` | C6H12 | 1 → 1 |
| 4 | `OOH* + C6H12* + * <-> TS1* -> C6H12O* + O* + Hb*` | C6H13O2 | 3 → 3 → 3 |
| 5 | `C6H12O* <-> C6H12O_g + *` | C6H12O | 1 → 1 |
| 6 | `Hb* + O* -> OH* + *` | HO | 2 → 2 |
| 7 | `Ha* + OH* -> H2O* + *` | H2O | 2 → 2 |
| 8 | `H2O* <-> H2O_g + *` | H2O | 1 → 1 |

Net stoichiometry: `H2O2_g + C6H12_g -> C6H12O_g + H2O_g`. All eight steps pass an independent elemental + effective site count check.

## Required energy provenance before running CatMAP

1. Provide consistent `surface_name/site_name/species_name/formation_energy/frequencies/reference` rows for `H2O2_s`, `Ha_s`, `OOH_s`, `C6H12_s`, `TS1_s`, `C6H12O_s`, `O_s`, `Hb_s`, `OH_s`, `H2O_s`, plus the four gas species; don't fabricate missing values.
2. **Do not reuse the six-step Rpair / TS1 energies unmodified**. In the six-step network, `TS1_s` formally contains `C6H14O2` in **one** coarse-grained effective center; this draft step 4 requires `C6H13O2` and **three** simultaneous positions. The original TS1 energy is not directly transferable; separate, reference-consistent TS energies or a justified alternative model are necessary. Likewise, historical standalone adsorption energies need clean-slab/thermal reference reconciliation before insertion.
3. The legacy `E(OOH*)` already used as a BEP descriptor is **not** itself a dynamically solved OOH* coverage. Fitting this descriptor vs TOF does *not* require rebuilding CatMAP; coverage requires the explicit network.
4. In CatMAP, two-state expressions without explicit TS specify the endpoint-barrier form; the `->` marker alone must **not** be treated as proof of strictly irreversible or instantaneous kinetics. If negligible reverse flux is desired, quantify or impose appropriate consistent kinetic assumptions in a separate step.
5. The three-site mean-field reaction in step 4 is a **hypothesis** about independent surface positions and their coverage/product dependence. It is not automatically equivalent to the old single-effective-center Ti–O–M coarse graining. Validate physical occupancy/geometry before interpreting TOF.
6. Once reference-consistent energies are available, run only a single Ti-Co CatMAP 0.3.1 test first, verify parse and mass/site conservation, convergence, flux consistency, and OOH* coverage; only then consider 13-site extension. Do not rerun the original six-step model or BEP grid for this task.

## Existing protected baseline

[Six-step Ti-Co setup](../../../catmap_official_run/results/official_strict_tolerance/tico/alkene_epoxidation.mkm) and [existing BEP-TOF outputs](../../BEP_OOH_activity_map/catmap_official_run/data/13site_BEP_fitted_TOF.csv) remain unchanged.
