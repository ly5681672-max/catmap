"""Retry only first-pass failures with a larger official damping-search limit."""

import csv
import shutil
import subprocess
from pathlib import Path
from decimal import Decimal


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE_DIR = ROOT / 'results' / 'tire_official'
RUN_ROOT = Path(__file__).resolve().parent
PYTHON = Path(r'G:\anaconda\python.exe')
SURFACES = ['tife', 'titi', 'ti', 'tini']


def read_csv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def main():
    template = (TEMPLATE_DIR / 'alkene_epoxidation.mkm').read_text(encoding='utf-8')
    runner = TEMPLATE_DIR / 'run_catmap_official.py'
    aggregate = []
    for surface in SURFACES:
        run_dir = RUN_ROOT / surface
        run_dir.mkdir(exist_ok=True)
        setup = template.replace("surface_names = ['tire']", "surface_names = ['{}']".format(surface))
        setup = setup.replace('resolution = 1', 'resolution = 1\nmax_damping_iterations = 50')
        if 'max_damping_iterations = 50' not in setup:
            raise RuntimeError('Could not add the uniquely located official damping option')
        (run_dir / 'alkene_epoxidation.mkm').write_text(setup, encoding='utf-8')
        shutil.copy2(TEMPLATE_DIR / 'energies.txt', run_dir / 'energies.txt')
        shutil.copy2(runner, run_dir / 'run_catmap_official.py')
        result = subprocess.run(
            [str(PYTHON), '-X', 'utf8', 'run_catmap_official.py'], cwd=str(run_dir),
            capture_output=True, text=True, encoding='utf-8', errors='replace', check=False,
        )
        (run_dir / 'catmap_console.log').write_text(result.stdout + result.stderr, encoding='utf-8')
        row = {'surface': surface, 'status': 'FAILED', 'returncode': result.returncode,
               'tof_C6H12O_g': '', 'residual': '', 'tolerance': '', 'coverage_sum': ''}
        conv_path = run_dir / 'convergence.csv'
        if result.returncode == 0 and conv_path.is_file():
            row.update(read_csv(conv_path)[0])
            row['returncode'] = result.returncode
            for item in read_csv(run_dir / 'gas_turnover_frequencies.csv'):
                if item['gas_species'] == 'C6H12O_g':
                    row['tof_C6H12O_g'] = item['turnover_frequency']
                    break
        aggregate.append(row)
        print('{}: {} (exit={})'.format(surface, row['status'], result.returncode))

    fields = ['surface', 'status', 'returncode', 'tof_C6H12O_g', 'residual',
              'tolerance', 'coverage_sum']
    with (RUN_ROOT / 'retry_summary.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(aggregate)

    first_pass = read_csv(ROOT / 'results' / 'official_14' / 'official_14_summary.csv')
    retry_by_surface = {row['surface']: row for row in aggregate}
    final_rows = []
    for row in first_pass:
        chosen = retry_by_surface.get(row['surface'], row)
        chosen = dict(chosen)
        chosen['selected_run'] = (
            'official_14_retry_damping50/' + row['surface']
            if row['surface'] in retry_by_surface
            else 'official_14/' + row['surface']
        )
        chosen['max_damping_iterations'] = 50 if row['surface'] in retry_by_surface else 10
        selected_dir = RUN_ROOT.parent / chosen['selected_run']
        coverage_rows = read_csv(selected_dir / 'coverage.csv')
        adsorbate_rows = [item for item in coverage_rows if item['species'] != 's']
        dominant = max(adsorbate_rows, key=lambda item: Decimal(item['coverage']))
        chosen['dominant_adsorbate'] = dominant['species']
        chosen['dominant_coverage'] = dominant['coverage']
        barriers = {item['step']: item for item in read_csv(selected_dir / 'effective_barriers.csv')}
        chosen['NEB1_forward_barrier'] = barriers['2']['forward_effective_barrier']
        chosen['NEB2_forward_barrier'] = barriers['4']['forward_effective_barrier']
        chosen['epoxide_desorption_forward_barrier'] = barriers['3']['forward_effective_barrier']
        final_rows.append(chosen)
    final_fields = ['surface', 'status', 'returncode', 'tof_C6H12O_g', 'residual',
                    'tolerance', 'coverage_sum', 'selected_run',
                    'max_damping_iterations', 'dominant_adsorbate',
                    'dominant_coverage', 'NEB1_forward_barrier',
                    'NEB2_forward_barrier', 'epoxide_desorption_forward_barrier']
    with (RUN_ROOT / 'final_14_summary.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=final_fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(final_rows)
    print('final_summary={}'.format(RUN_ROOT / 'final_14_summary.csv'))


if __name__ == '__main__':
    main()
