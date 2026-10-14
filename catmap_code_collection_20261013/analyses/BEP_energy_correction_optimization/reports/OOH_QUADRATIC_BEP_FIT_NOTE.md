# ΔG₂ + OOH* + (OOH*)²: four-parameter quadratic NEB1 BEP fitting

Source: `../data/source_carbon_O_H_TS1.csv` (13 complete Ti/M sites, Ti–Zr unresolved and excluded). Read-only fit; no new DFT, NEB or CatMAP run. All differences are signed G(TS1*)−G(Rpair*) and **should not be called positive activation barriers**.

Centered predictors: `u = ΔG₂ + 2.092276692 eV`; `z = E(OOH*) + 4.309874385 eV` (historical single-adsorbate formation energy).

Model: `Y = 0.412123097 + 0.370377301*u + 0.200217496*z - 0.309777397*z²`.

- n = 13; number of coefficients = 4 (intercept + linear ΔG₂ + linear OOH + quadratic OOH).
- Training R² = **0.8866652**; training RMSE = **0.1577516 eV**.
- LOOCV RMSE (refitting 12 sites in every fold) = **0.3814319 eV**.
- Existing traditional BEP training R² = 0.2054; do not compare fitting and predictive metrics as identical.

Files:
- `../figures/08_OOH_quadratic_BEP_fit.svg`: the committed **vector** parity/residual figure, readable on GitHub.
- `../data/08_OOH_quadratic_BEP_predictions.csv`: all 13 sites' original, fitted, LOOCV, residuals.
- `../scripts/plot_ooh_quadratic_bep.py`: regenerate original image as `08_OOH_quadratic_BEP_fit.png/.pdf/.svg` using archived CSV in-place.

Run:

```powershell
cd H:\catmap\catmap\catmap_code_collection_20261013\analyses\BEP_energy_correction_optimization
python .\scripts\plot_ooh_quadratic_bep.py
```

**Important:** the OOH descriptor comes from legacy adsorption energies; single-adsorbate reference/slab and strict six-step network have not been reconciled. This is a high-R² exploratory fit, not a separately validated mechanism. No input/TOF files are changed. The SVG is a separate vector illustration of the exact same numerical fit; the Python script generates the PNG/PDF publication versions.
