"""Read-only original-eight-stage site and NEB mapping diagnostic; no kinetics."""
from pathlib import Path
from collections import Counter
from decimal import Decimal as D
import csv
import hashlib
import inspect
import json
import os
import shutil

ROOT = Path(__file__).resolve().parents[1]
PROJECT = Path(r'H:\catmap\catmap\catmap_code_collection_20261013')
SOURCE = PROJECT / 'analyses/explicit_OOH_8step_validation/candidate_tico_kinetic_control/candidate_energies.txt'
ORIGINAL = Path(r'H:\catmap\catmap\lyq_alkene_epoxidation\2-creating_microkinetic_model\alkene_epoxidation.mkm')


def csv_out(name, rows):
    with (ROOT / 'data' / name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    from catmap import ReactionModel
    import catmap
    import mpmath as mp
    from catmap.solvers.mean_field_solver import MeanFieldSolver
    for name in ('data', 'inputs', 'reports'):
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    inputs = [ORIGINAL, ORIGINAL.with_name('energies_dft_neb.txt'), SOURCE]
    before = {str(p): digest(p) for p in inputs}
    # Copy only into this isolated diagnostic directory; original setup is unchanged.
    for p in inputs[:2]:
        shutil.copy2(p, ROOT / 'inputs' / p.name)
    os.chdir(ROOT / 'inputs')
    m = ReactionModel(setup_file=ORIGINAL.name)
    m.parse()
    m.verify()
    names = list(m.adsorbate_names)
    rxns = list(m.elementary_rxns)
    assert len(names) == 9 and len(rxns) == 8
    assert set(m.transition_state_names) == {'OOH-C6H12_s', 'Hb-O_s'}
    assert all(m.species_definitions[n]['n_sites'] == 1 for n in names)
    matrix = []
    for n in names:
        row = {'species': n, 'n_sites': m.species_definitions[n]['n_sites']}
        row.update({f'r{j+1}': Counter(r[-1])[n]-Counter(r[0])[n] for j, r in enumerate(rxns)})
        matrix.append(row)
    csv_out('ORIGINAL8_STOICHIOMETRY.csv', matrix)
    inventories = [dict(Ha_s=1, OOH_s=-1, OH_s=-1, O_s=-1), dict(Hb_s=1, O_s=-1)]
    for inv in inventories:
        assert all(sum(inv.get(n, 0) * row[f'r{j+1}'] for n, row in zip(names, matrix)) == 0 for j in range(8))
    m._mpfloat = mp.mpf
    baseline = m.solver.site_string_list()
    # In-memory probe only: changing one metadata n_sites must not be interpreted
    # as a supported physical model. No verify/solve is performed on this probe.
    old_n = m.species_definitions['OOH_s']['n_sites']
    m.species_definitions['OOH_s']['n_sites'] = 3
    changed = m.solver.site_string_list()
    m.species_definitions['OOH_s']['n_sites'] = old_n
    assert baseline == changed
    cv = [D('0.01')] * len(names)
    unweighted = D(1)-sum(cv)
    weighted_probe = unweighted-D('0.02')
    rows = list(csv.DictReader(SOURCE.open(encoding='utf-8-sig'), delimiter='\t'))
    energies = {r['species_name']: D(r['formation_energy']) for r in rows if r['surface_name'] == 'tico'}
    original_rows = list(csv.DictReader(inputs[1].open(encoding='utf-8-sig'), delimiter='\t'))
    original_energies = {r['species_name']: D(r['formation_energy']) for r in original_rows if r['surface_name'] == 'tico'}
    assert all(abs(original_energies[n.removesuffix('_s')]-energies[n.removesuffix('_s')]) < D('1e-9')
               for n in names + list(m.transition_state_names))
    edict = {n+'_s': float(v) for n, v in energies.items()}
    edict.update({n: 0.0 for n in m.site_names})
    mapping = []
    for label, index in [('NEB1', 3), ('NEB2', 5)]:
        gis, gts, gfs = [energies[label+s] for s in ('IS', 'TS', 'FS')]
        delta, activation = m.get_rxn_energy(rxns[index], edict)
        independent = D(str(round(delta, 10)))
        forward, reverse = gts-gis, gts-gfs
        # A single TS shift cannot change the forward-minus-reverse identity.
        assert abs((forward-reverse)-(gfs-gis)) < D('1e-12')
        mismatch = independent-(gfs-gis)
        assert abs(mismatch) > D('1')
        mapping.append(dict(path=label, original_step=index+1, full_IS=gis, full_TS=gts,
                            full_FS=gfs, full_dG=gfs-gis, signed_TS_IS=forward,
                            signed_TS_FS=reverse, independent_dG=independent,
                            independent_minus_full_dG=mismatch,
                            legacy_effective_forward=activation,
                            legacy_effective_reverse=activation-delta,
                            status='FAIL_DIRECT_ENDPOINT_EQUIVALENCE'))
    csv_out('NEB_MAPPING_IDENTITY.csv', mapping)
    assert all(digest(p) == before[str(p)] for p in inputs)
    source_path = Path(inspect.getfile(MeanFieldSolver))
    result = dict(catmap_version=catmap.__version__, stable_species=names,
                  transition_states=list(m.transition_state_names), chemical_stages=8,
                  elementary_expressions=len(rxns), original_n_sites=1,
                  geometric_site_count='UNVERIFIED_NOT_FRAGMENT_COUNT',
                  site_strings=baseline, metadata_probe_site_strings=changed,
                  n_sites_metadata_probe_changes_vacancy=False,
                  test_coverage_per_species='0.01', native_vacancy=str(unweighted),
                  hypothetical_weighted_probe_vacancy=str(weighted_probe),
                  inventory_checks='PASS', parser_and_balance='PASS',
                  direct_NEB_mapping='FAIL', steady_state_run=False,
                  original_energy_rows_match_reused_table=True,
                  original_hashes_unchanged=True,
                  source_hashes=before,
                  official_mean_field_source=dict(path=str(source_path), sha256=digest(source_path)))
    (ROOT / 'data' / 'VALIDATION.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
