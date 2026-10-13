"""Run CatMAP 0.3.1's built-in rate_control output on five selected sites."""

import csv
import shutil
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "results" / "official_strict_tolerance"
PYTHON = Path(r"G:\anaconda\python.exe")
SURFACES = ["ti", "tico", "tini", "tire", "tife"]
EXPECTED = "output_variables = ['coverage', 'rate', 'turnover_frequency', 'directional_rates']"


def build_run(surface):
    run_dir = HERE / "official_drc" / surface
    run_dir.mkdir(parents=True, exist_ok=True)
    source_dir = SOURCE / surface
    setup = (source_dir / "alkene_epoxidation.mkm").read_text(encoding="utf-8")
    if setup.count(EXPECTED) != 1:
        raise RuntimeError(f"Expected one output_variables anchor in {surface} setup")
    setup = setup.replace(EXPECTED, "output_variables = ['rate_control']")
    (run_dir / "alkene_epoxidation.mkm").write_text(setup, encoding="utf-8")
    shutil.copy2(source_dir / "energies.txt", run_dir / "energies.txt")
    script = run_dir / "catmap_rate_control.py"
    script.write_text(
        """import csv
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
""",
        encoding="utf-8",
    )
    result = subprocess.run(
        [str(PYTHON), "-X", "utf8", script.name],
        cwd=run_dir,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    log = result.stdout + result.stderr
    (run_dir / "catmap_drc_console.log").write_text(log, encoding="utf-8")
    if result.returncode != 0 or not (run_dir / "rate_control.csv").is_file():
        raise RuntimeError(f"Official DRC failed for {surface}; inspect {run_dir / 'catmap_drc_console.log'}")
    if "jacobian_fail" in log.lower():
        raise RuntimeError(f"CatMAP logged a DRC Jacobian failure for {surface}")
    with (run_dir / "rate_control.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    product = [r for r in rows if r["gas_species"] == "C6H12O_g"]
    if not product or all(float(r["degree_of_rate_control"]) == 0 for r in product):
        raise RuntimeError(f"Product DRC missing or all zero for {surface}")
    print(f"{surface}: official DRC exported ({len(product)} parameters)")


def main():
    if not PYTHON.is_file():
        raise FileNotFoundError(PYTHON)
    for surface in SURFACES:
        build_run(surface)
    combined = HERE / "official_drc" / "five_site_product_drc.csv"
    with combined.open("w", encoding="utf-8-sig", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["surface", "gas_species", "energy_parameter", "degree_of_rate_control"])
        for surface in SURFACES:
            with (HERE / "official_drc" / surface / "rate_control.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                for row in csv.DictReader(handle):
                    if row["gas_species"] == "C6H12O_g":
                        writer.writerow([surface, row["gas_species"], row["energy_parameter"], row["degree_of_rate_control"]])
    print(f"combined_product_DRC={combined}")


if __name__ == "__main__":
    main()
