# Scheme B: coadsorbed eight-step microkinetic model

This is a **separate, non-destructive alternative** within lyq_alkene_epoxidation. The original `../alkene_epoxidation.mkm`, `../energies.txt` and `../test.ipynb` are left unchanged, so original Scheme A remains reproducible.

## Origin and definitions

Single source of truth: `../../己烯环氧化.xlsx` (`吉布斯自由能汇总` plus `catmap数据集`). All energies use the original μ(C)=-9.28, μ(H)=-1.11, μ(O)=-4.37 eV and **one slab per entire coadsorbed state**.

- S1 = H2O2 (summary row 8); S2 = Ha+OOH, **proxy sum of separately adsorbed species** (single-species sheet rows 8,16).
- S3 = Ha+OOH+C6H12 (summary row 12), TS1 = epoxidation (row 16), S4 = Ha+Hb+O+C6H12O (row 20, TiMo uses row 21 total including gas product).
- S5 = Ha+Hb+O (summary row 24), TS2 = Hb transfers H to O (row 28), S6 = Ha+OH (row 32; TiMo's combined row 32 has an included gas epoxide removed to isolate this state).
- S7 = adsorbed H2O (single-species sheet row 39).

CatMAP keeps eight net chemical steps (S1 to S7 coarse-grained), treating each S-state as **one complete catalytic unit**. It cannot be compared one-to-one with Scheme A's individual adsorbate coverages. The epoxide on TiMo was already gas phase; its summed total is used as an **effective lumped state**, not presented as a new adsorbed DFT geometry.

## Running and regeneration (PowerShell)

From `lyq_alkene_epoxidation/2-creating_microkinetic_model/scheme_B`:

```powershell
python build_scheme_B.py
python run_scheme_B.py
# Regenerate from Excel only if desired:
python build_scheme_B.py --write
```

This uses **Python standard library** for XLSX extraction (no openpyxl/pandas needed). The output `energies.txt` has 15 catalysts; the .mkm fits the 14 original catalysts and still excludes TiNb. The new descriptors are S3_s and S5_s (old OOH_s and O_s are no longer explicit species).

## Data checks and limitations

- All 15 original Excel TS1 values are above the corresponding S4 final states. The original Scheme A showed 9/15 TS1-below-FS cases caused by mixing independent-state and coadsorbed-state energies.
- Each catalyst's eight step ΔG values sum to **−2.089 eV**, identical to gas net ΔG (absolute error ≲ 1e-13 eV).
- The **original** nonphysical TS-position flags remain: TS1 below IS for TiNb, TiW, TiTa; TS2 below IS for TiMo, TiV, TiCr; TS2 below FS for Ti, TiZr. No TS values were raised/clamped in the input. Effective kinetic barriers can differ when the CatMAP solver uses max(IS, TS, FS).
- The XLSX G values represent free energies at 333.15 K; `frozen_gas` / `frozen_adsorbate` avoid double applying a thermal correction. Gas pressure parameters represent **effective liquid-phase activities**, not literal pressures.
- S2 and S7 are **approximate states** from the existing single-adsorption data. B's full-complex state model is **not a site-resolved coadsorption interaction model**.
- Scaling over new descriptors is an additional approximation; first validate all actual catalyst input states/flags before interpreting interpolated volcano maps.

**Validation boundary:** Static energy, stoichiometry, and closure checks passed before commit. Full CatMAP kinetic convergence was **not** run in the remote-editing environment (CatMAP module unavailable); run `python run_scheme_B.py` in your local CatMAP environment.

## Windows file encoding

Some CatMAP installations use Windows' GBK default when reading `energies.txt`. The `reference` column therefore uses **ASCII-only labels** (`Excel:GibbsSummary/...`, `Excel:SingleAds/...`). These are provenance descriptions only and do **not** change any numerical energies. The original Excel worksheet names inside `build_scheme_B.py` remain in Chinese.

If `UnicodeDecodeError: 'gbk'` occurred on an older checkout, pull the newest main branch and rerun `python build_scheme_B.py` followed by `python run_scheme_B.py`.

## Net-flux reliability and 14 exact catalyst results

After the CatMAP run, run_scheme_B.py automatically audits the saved map and solves 14 measured catalyst points at their own descriptor coordinates, not nearest grid points. Input energies, the Excel file, and all eight reaction steps are untouched.

Output directory: diagnostics/ (UTF-8 BOM CSV, compatible with Excel):

- grid400_diagnostics.csv: 400 regular points, all eight net step rates, steady-state flux consistency, S1-S7/empty-site coverages, the two fitted TS energy positions, and forward/reverse effective kinetic barriers.
- grid400_volcano_eligible.csv: only qualified points, sorted by positive net epoxide rate (step 5).
- bisection_extra_points.csv: additional solver bisection points, not counted as grid nodes.
- catalysts_input_barriers.csv: measured catalyst TS energy gaps, without scaling.
- catalysts_exact.csv: newly solved actual catalyst points (14) and fitted-vs-measured energy gaps, with explicit failed-solve labels.
- summary.csv: aggregate pass/fail and energy-position counts.

Flux quality criterion: max(abs(r_i-r_5)) / max_i(abs(r_i)) <= 0.01 for the eight net rates, plus physical total coverage and positive product rate. The optional min-signal floor defaults to 0; no arbitrary kinetic cutoff is imposed. Near-zero forward-minus-reverse cancellation can still produce numerical flux discrepancies, so excluded points must not be ranked in the volcano plot.

Run from the scheme_B folder in the same CatMAP environment:

    python run_scheme_B.py

For already saved scheme_B.pkl, without re-running CatMAP and without exact catalyst solves:

    python analyze_scheme_B.py

The optional legacy form `python analyze_scheme_B.py --grid-only` is also supported.
Running this file directly in PyCharm analyzes the saved 400-point map **only**;
it never performs the 14 exact catalyst solves. To calculate those, run
`python run_scheme_B.py` instead.

The uploaded earlier saved map gave 240/400 flux-consistent grid points and 160/400 excluded; exact catalyst solves require the user's installed CatMAP. A solved grid is not proof that all transition-state positions are physically validated.
