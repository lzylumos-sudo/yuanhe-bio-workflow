# -*- coding: utf-8 -*-
"""Smoke tests for the pepwedge pipeline (no pytest required).

    python tests/test_smoke.py
"""
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np

from pepwedge.core import COMPLEXES, contacts, load_complex, membrane_frame, residue_sasa
from pepwedge.design import HEADS, LINKS, PHENYLENE, TAILS, compose, enumerate_scaffolds
from pepwedge.pharmacophore import consensus_features
from pepwedge.pocket import SUBPOCKETS, build_pocket

FAILED = []


def check(name, cond, detail=''):
    if cond:
        print(f'  ok   {name} {detail}')
    else:
        print(f'  FAIL {name} {detail}')
        FAILED.append(name)


def test_structures():
    print('[structures]')
    cxs = {t: load_complex(os.path.join(ROOT, f), t) for t, f in COMPLEXES.items()}
    for tag, cx in cxs.items():
        n_rec, n_lig = len(cx.receptor), len(cx.ligand)
        check(f'{tag} chains', n_rec > 300 and 10 < n_lig < 60, f'{n_rec} MurJ / {n_lig} peptide')
        ct = contacts(cx)
        check(f'{tag} contacts', len(ct) > 10, f'{len(ct)} residues')
    n, c0 = membrane_frame(cxs['SglCJ3'])
    check('membrane normal', abs(np.linalg.norm(n) - 1.0) < 1e-6, f'n={np.round(n, 2)}')


def test_assembly():
    print('[assembly]')
    sc = compose(HEADS['COOH'], [PHENYLENE, PHENYLENE],
                 [LINKS['bond'], LINKS['alkyne'], LINKS['bond']], TAILS['Ph'])
    check('span = 14.02', abs(sc.span - 14.02) < 0.02, f'{sc.span:.2f} A')
    check('SMILES balanced', sc.smiles.count('(') == sc.smiles.count(')'), sc.smiles)
    check('elements', len(sc.elements) == len(sc.atoms))
    check('no clashing atom pairs', min(
        float(np.linalg.norm(sc.atoms[i] - sc.atoms[j]))
        for i in range(len(sc.atoms)) for j in range(i + 1, len(sc.atoms))) > 1.0)
    check('props in range', 250 < sc.mw < 350 and 0 < sc.hbd <= 2, f'MW {sc.mw:.0f}')


def test_enumeration():
    print('[enumeration]')
    a = enumerate_scaffolds(14.1, tol=0.55, ring_counts=(2, 3))
    b = enumerate_scaffolds(13.3, tol=0.55, ring_counts=(2, 3))
    check('D=14.1 non-empty', len(a) > 200, f'{len(a)} skeletons')
    check('D=13.3 non-empty', len(b) > 100, f'{len(b)} skeletons')
    check('span within tolerance', all(abs(s.span - 14.1) <= 0.56 for s in a))
    check('unique SMILES', len({s.smiles for s in a}) == len(a),
          f'{len({s.smiles for s in a})}/{len(a)}')
    check('deterministic', [s.smiles for s in a] == [s.smiles for s in enumerate_scaffolds(
        14.1, tol=0.55, ring_counts=(2, 3))])


def test_pharmacophore_and_pocket():
    print('[pharmacophore / pocket]')
    cxs = {t: load_complex(os.path.join(ROOT, f), t) for t, f in COMPLEXES.items()}
    feats = consensus_features(cxs, min_structures=1, max_features=12)
    check('feature count', 8 <= len(feats) <= 12, f'{len(feats)}')
    anions = [f for f in feats if f.ftype == 'anion']
    check('has anionic anchors', len(anions) >= 2, f'{len(anions)}')
    pk = build_pocket(cxs['SglCJ3'])
    keep = pk.buriedness >= 0.5
    vol = float(keep.sum() * pk.step ** 3)
    check('buried volume', vol > 2000, f'{vol:.0f} A^3')
    check('subpockets partitioned', set(pk.subpocket) == set(SUBPOCKETS),
          f'{sorted(set(pk.subpocket))}')


def test_end_to_end():
    print('[end-to-end]')
    from pepwedge.pipeline import run
    out = os.path.join(ROOT, 'results', '_smoke')
    shutil.rmtree(out, ignore_errors=True)
    res = run(ROOT, out, verbose=False, figdir=os.path.join(out, 'figures'))
    check('scaffolds scored', len(res['results']) > 500, f'{len(res["results"])}')
    check('csv written', os.path.exists(os.path.join(out, 'scaffolds_ranked.csv')))
    check('sdf written', os.path.exists(os.path.join(out, 'top_scaffold_poses.sdf')))
    top = res['rows'][0]
    check('top has salt bridges', int(top['salt']) >= 1, top['smiles'])
    check('top contacts anchor residue',
          any(k in top['murj_contacts'] for k in ('K46', 'R53', 'R270', 'K368')),
          top['murj_contacts'])
    shutil.rmtree(out, ignore_errors=True)




def test_priors():
    print('[priors]')
    from pepwedge import priors as P
    d = P.load_priors(strict=True)
    check('schema', d['schema_version'] == '1.0', d['schema_version'])
    check('residue range', P.residue_range() == (5, 500), str(P.residue_range()))
    check('structures', len(P.structure_map()) == 3, str(P.structure_map()))
    check('essential anchors', set(P.essential_map()) == {46, 53, 270, 368},
          str(sorted(P.essential_map())))
    check('resistance alleles', set(P.resistance_map()) == {50, 244},
          str(sorted(P.resistance_map())))
    check('peptide priors', sum(len(v) for v in P.peptide_functional_map().values()) == 18)
    # a malformed prior must be rejected, not silently accepted
    import copy
    bad = copy.deepcopy(d)
    bad['murj_residues'][0]['verify'] = {'type': 'not_a_rule'}
    rejected = False
    try:
        P.validate(bad)
    except P.PriorError:
        rejected = True
    check('invalid verify type rejected', rejected)
    # provenance is mandatory: a prior without its quote must be refused
    no_quote = copy.deepcopy(d)
    no_quote['peptide_residues'][0].pop('quote', None)
    try:
        P.validate(no_quote)
        check('prior without quote rejected', False)
    except P.PriorError:
        check('prior without quote rejected', True)
    # a peptide rule attached to a MurJ residue must be refused, not silently
    # dispatched to the wrong verifier
    wrong_sec = copy.deepcopy(d)
    wrong_sec['murj_residues'][0]['verify'] = {'type': 'peptide_contact'}
    try:
        P.validate(wrong_sec)
        check('cross-section verify type rejected', False)
    except P.PriorError:
        check('cross-section verify type rejected', True)


def test_construct_filter():
    print('[construct filter]')
    cxs = {t: load_complex(os.path.join(ROOT, f), t) for t, f in COMPLEXES.items()}
    for tag, cx in cxs.items():
        nums = sorted(cx.receptor)
        check(f'{tag} MurJ range', nums[0] >= 5 and nums[-1] <= 500, f'{nums[0]}-{nums[-1]}')
        check(f'{tag} 496 residues', len(cx.receptor) == 496, str(len(cx.receptor)))
    check('9NU4 construct residues', len(cxs['SglM'].extra) == 106, str(len(cxs['SglM'].extra)))
    check('9NU5 construct residues', len(cxs['SglPP7'].extra) == 7)
    check('9NU8 construct residues', len(cxs['SglCJ3'].extra) == 14)


def test_evidence_chain():
    print('[evidence chain]')
    from pepwedge.evidence import (atom_budget, anchor_table, census, dsasa_both,
                                   interface_census)
    from pepwedge.verify import run_checks, summarize
    from pepwedge import priors as P

    cxs = {t: load_complex(os.path.join(ROOT, f), t) for t, f in COMPLEXES.items()}
    dsasa, ic = {}, {}
    for tag, cx in cxs.items():
        dsasa[tag] = dsasa_both(cx)
        ic[tag] = interface_census(cx, dsasa[tag])
    for tag, v in ic.items():
        check(f'{tag} interface size', 20 <= v['murj_residues_in_interface'] <= 45,
              f"{v['murj_residues_in_interface']} MurJ / {v['peptide_residues_in_interface']} peptide")
        check(f'{tag} two-sided area agrees', v['asymmetry'] < 0.15, str(v['asymmetry']))
    check('peptide residues in interface', ic['SglCJ3']['peptide_residues_in_interface'] == 28)

    cen = {t: census(cx) for t, cx in cxs.items()}
    check('census atoms identical across structures',
          len({c['murj_atoms'] for c in cen.values()}) == 1,
          str([c['murj_atoms'] for c in cen.values()]))

    at = anchor_table(cxs['SglCJ3'])
    check('anchor pairs computed', len(at['anchor_pairs']) == 6, str(len(at['anchor_pairs'])))

    b = atom_budget(cxs, dsasa_by_tag=dsasa)
    check('small-molecule heavy atoms', 15 <= b['small_molecule']['heavy_atoms_median'] <= 35,
          str(b['small_molecule']))

    checks = run_checks(P.load_priors(), cxs, dict(dsasa=dsasa, contacts={
        t: contacts(cx) for t, cx in cxs.items()}))
    summ = summarize(checks)
    check('checks count', summ['total'] == 24, str(summ['total']))
    check('no unexplained failure', summ['failed'] == 0 and summ['mismatch'] == 0,
          str(summ['by_status']))
    check('pass rate', summ['passed'] >= 20, f"{summ['passed']}/{summ['total']}")




def test_docs_match_data():
    """The docs quote concrete numbers; they must not drift from the outputs."""
    print('[docs vs data]')
    import csv
    import json
    import re
    doc_path = os.path.join(ROOT, 'docs', '01_分子骨架设计算法.md')
    res = os.path.join(ROOT, 'results')
    doc = open(doc_path, encoding='utf-8').read()
    ph = json.load(open(os.path.join(res, 'pharmacophore.json'), encoding='utf-8'))
    pk = json.load(open(os.path.join(res, 'pocket.json'), encoding='utf-8'))
    su = json.load(open(os.path.join(res, 'summary.json'), encoding='utf-8'))
    rows = list(csv.DictReader(open(os.path.join(res, 'interface_hotspots.csv'), encoding='utf-8')))

    bad = []
    for f in ph['features']:                       # table lists the top 5 only
        m = re.search(r'\*\*' + f['fid'] + r'\*\* \| [^|]+ \| \(([-\d.]+), ([-\d.]+), ([-\d.]+)\) '
                      r'\| ([\d.]+) \| ([\d.]+)', doc)
        if not m:
            continue
        pos = [float(m.group(i)) for i in (1, 2, 3)]
        if max(abs(a - b) for a, b in zip(pos, f['pos'])) > 0.1:
            bad.append(f['fid'] + ' pos')
        if abs(float(m.group(4)) - f['convergence']) > 1e-9 or abs(float(m.group(5)) - f['weight']) > 0.01:
            bad.append(f['fid'] + ' weight')

    short = {'S1_basic_bowl': 'S1', 'S2_carboxylate_gate': 'S2',
             'S3_TM2_TM7_groove': 'S3', 'S4_TM7_floor': 'S4'}
    for sp in pk['subpockets']:
        key = short[sp['subpocket']]
        m = re.search(r'\*\*' + key + r'[^|]*\|[^|]*\|[^|]*\|[^|]*\| ([^|]+) \|', doc)
        got = [x.strip() for x in m.group(1).split(',')] if m else None
        if got != sp['residues']:
            bad.append(key + ' lining residues')

    for label, key in (('enumerated', 'enumerated'), ('scored', 'property_pass'),
                       ('no_clash', 'no_clash'), ('mimic', 'peptide_mimic'),
                       ('drug_like', 'drug_like')):
        if f"{su['funnel'][key]:,}" not in doc:
            bad.append('funnel ' + label)

    vol_m = re.search(r'\*\*([\d,]+) Å³\*\* 的埋藏腔体', doc)
    vol_doc = float(vol_m.group(1).replace(',', '')) if vol_m else -1.0
    if abs(vol_doc - pk['total_volume_A3']) > 1.0:
        bad.append(f'pocket volume {vol_doc} vs {pk["total_volume_A3"]:.1f}')
    check('docs match current outputs', not bad, ', '.join(bad))
    check('summary best has contacts', len(su['best'][0]['contacts']) > 0)
    check('summary best has heavy_atoms', su['best'][0]['props']['heavy_atoms'] > 10)




def test_filters():
    print('[filters]')
    import json
    from pepwedge import filters as F
    cxs = {t: load_complex(os.path.join(ROOT, f), t) for t, f in COMPLEXES.items()}
    pocket = json.load(open(os.path.join(ROOT, 'results', 'pocket.json'), encoding='utf-8'))
    cal = F.calibration_report(cxs, pocket)
    check('calibration has 3 entries', len(cal) == 3, str(len(cal)))
    # the peptides must satisfy the thresholds the filters use -- that is the
    # whole point of calibrating on them instead of guessing
    ok = [c for c in cal if 'anchor' in c]
    check('peptide satisfies anchor cutoff',
          all(c['dist'] <= F.ANCHOR_CUT for c in ok), str(ok))
    check('peptide shows no like-charge contact',
          all(c['repulsion'] > F.REPULSION_CUT for c in ok), str(ok))
    polar, hydro = F._pocket_residues(pocket)
    check('pocket residues split', len(polar) > 0 and len(hydro) > 0,
          f'{len(polar)}/{len(hydro)}')


if __name__ == '__main__':
    test_structures()
    test_assembly()
    test_enumeration()
    test_pharmacophore_and_pocket()
    test_priors()
    test_construct_filter()
    test_evidence_chain()
    test_end_to_end()
    test_docs_match_data()
    test_filters()
    print()
    if FAILED:
        print(f'FAILED: {len(FAILED)} -> {FAILED}')
        sys.exit(1)
    print('all smoke tests passed')
