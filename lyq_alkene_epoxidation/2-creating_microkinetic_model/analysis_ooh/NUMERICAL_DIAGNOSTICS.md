# OOH numerical diagnostics — status before rerun

The existing CSV and OOH_logTOF_linear.png/.pdf in this directory were produced using the **previous** single-point model (100-digit precision; `tolerance=1e-50`). They are **historical exploratory output, not validated TOF results**.

## Why the original results were not valid

In the CatMAP numbers-based solver, the convergence norm is the **sum of squared** changes in surface coverages. Hence, a configured `tolerance=1e-50` can accept individual residuals around `1e-25`. The previous Ti–Fe result had an elementary epoxidation rate around `2.24e-56` but an epoxide desorption rate around `3.77e-26`, which is incompatible with a genuine steady-state single-cycle flux.

Source: [CatMAP numerical-accuracy tutorial](https://catmap.readthedocs.io/en/latest/tutorials/refining_a_microkinetic_model.html) and [numbers_solver.py](https://github.com/SUNCAT-Center/catmap/blob/master/catmap/solvers/numbers_solver.py).

## Rerun in the CatMAP environment

```powershell
python -m unittest -v test_ooh_numerics
python run_ooh_analysis.py --diagnose
python run_ooh_analysis.py --diagnose --coverage-crosscheck
python run_ooh_analysis.py
```

- `--diagnose` analyzes only Ti–Fe, Ti–Ti, Ti–W under both old and refined number-solver settings. Results go to `analysis_ooh/diagnostics/diagnostics_comparison.csv`. It does **not** overwrite `dft_baseline.csv`.
- Refined settings default to `decimal_precision=180`, `tolerance=1e-120`, `max_rootfinding_iterations=250`, with fresh per-surface caches.
- The optional coverage-solver crosscheck uses a non-squared maximum residual criterion.
- The 14-surface rerun exports eight elementary net rates, an independent physical flux check, gas balance, raw coverage residual, and quality status.
- Only `quality_pass=True` **and** no known abnormal TS flags can enter the new OOH–TOF regression. If fewer than three valid points remain, the local obsolete TOF plot is deleted rather than depicting invalid fit results.
- `min_tof=1e-40` is a **diagnostic significance cutoff**, not a physical TOF limit; adjust only with a numerical error study.

**Scientific caveat:** numerical convergence does not validate the DFT-derived activation energies, solvation/activity approximations, or catalyst rankings. The OOH–OH DFT-energy correlation is independent of numerical TOF convergence.

This patch does not change `alkene_epoxidation.mkm`, `test.ipynb`, the DFT input tables, or the seven legacy CatMAP figure recipes.
