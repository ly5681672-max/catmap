"""Run one independent official CatMAP 0.3.1 model per approved surface."""

import csv
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE_DIR = ROOT / 'results' / 'tire_official'
RUN_ROOT = Path(__file__).resolve().parent
PYTHON = Path(r'G:\anaconda\python.exe')
SURFACES = [
    'tife', 'timn', 'tihf', 'tire', 'timo', 'tiv', 'tizr', 'tico',
    'titi', 'tiw', 'tita', 'ticr', 'ti', 'tini',
]


def read_single_csv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def main():
    if not PYTHON.is_file():
        raise FileNotFoundError(PYTHON)
    template = (TEMPLATE_DIR / 'alkene_epoxidation.mkm').read_text(encoding='utf-8')
    runner = TEMPLATE_DIR / 'run_catmap_official.py'
    energies = TEMPLATE_DIR / 'energies.txt'
    aggregate = []

    for surface in SURFACES:
        run_dir = RUN_ROOT / surface
        run_dir.mkdir(exist_ok=True)
        setup = template.replace("surface_names = ['tire']", "surface_names = ['{}']".format(surface))
        if setup == template and surface != 'tire':
            raise RuntimeError('Could not uniquely substitute surface in the CatMAP setup')
        (run_dir / 'alkene_epoxidation.mkm').write_text(setup, encoding='utf-8')
        shutil.copy2(energies, run_dir / 'energies.txt')
        shutil.copy2(runner, run_dir / 'run_catmap_official.py')

        result = subprocess.run(
            [str(PYTHON), '-X', 'utf8', 'run_catmap_official.py'],
            cwd=str(run_dir), capture_output=True, text=True,
            encoding='utf-8', errors='replace', check=False,
        )
        console = result.stdout + result.stderr
        (run_dir / 'catmap_console.log').write_text(console, encoding='utf-8')

        row = {'surface': surface, 'status': 'FAILED', 'returncode': result.returncode,
               'tof_C6H12O_g': '', 'residual': '', 'tolerance': '', 'coverage_sum': ''}
        convergence_file = run_dir / 'convergence.csv'
        if result.returncode == 0 and convergence_file.is_file():
            conv = read_single_csv(convergence_file)[0]
            row.update(status=conv['status'], residual=conv['residual'],
                       tolerance=conv['tolerance'], coverage_sum=conv['coverage_sum'])
            for item in read_single_csv(run_dir / 'gas_turnover_frequencies.csv'):
                if item['gas_species'] == 'C6H12O_g':
                    row['tof_C6H12O_g'] = item['turnover_frequency']
                    break
        else:
            row['failure_log'] = 'catmap_console.log'
        aggregate.append(row)
        print('{}: {} (exit={})'.format(surface, row['status'], result.returncode))

    fields = ['surface', 'status', 'returncode', 'tof_C6H12O_g', 'residual',
              'tolerance', 'coverage_sum', 'failure_log']
    with (RUN_ROOT / 'official_14_summary.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(aggregate)

    failed = [row['surface'] for row in aggregate if row['status'] != 'CONVERGED']
    print('summary={}'.format(RUN_ROOT / 'official_14_summary.csv'))
    print('failed_surfaces={!r}'.format(failed))
    if failed:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
