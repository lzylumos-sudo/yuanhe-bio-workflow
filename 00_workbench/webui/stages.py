# -*- coding: utf-8 -*-
"""webui.stages -- data contract of the 7-step workflow UI.

Stages 1-2 are backed by files the existing pipeline already wrote.  Stages 3-7
(generation, docking, screening, MD, results) are *slots*: the UI lights up as
soon as the files listed in ``CONTRACT`` appear under ``results/``.  Nothing is
generated here and no model is called.

File contract (all optional, UTF-8, relative to ``results/``)
-------------------------------------------------------------
<key>.csv           one row per molecule.  Required columns: ``id``, ``smiles``.
                    Any further column is shown in the table; numeric columns
                    can be used for sorting.  An optional ``keep`` column
                    (true/false) marks molecules that pass this stage.
<key>_summary.json  ``{"funnel": [{"label": "...", "count": 123}, ...]}``
                    -- the screening-count change shown on that page.
<key>_poses.sdf     optional 3D poses in the SglCJ3 (9NU8) frame; the molecule
                    is matched by SDF title line == ``id`` (or == ``smiles``).

``key`` is one of: generated, docking, screening, md, final_hits.
"""
from __future__ import annotations

import json
import os

from . import data

ELEMENT = data.ELEMENT_COLORS

STAGES = [
    dict(n=1, key='complex', title='Target complex', short='Import',
         subtitle='MurJ complexes with three phage lysis peptides',
         unit='Complex', mode='builtin'),
    dict(n=2, key='scaffold', title='Scaffold design', short='Scaffold',
         subtitle='Assemble and filter scaffolds using key-residue pharmacophores',
         unit='Scaffold', mode='builtin'),
    dict(n=3, key='generated', title='Batch generation', short='Generate',
         subtitle='Scaffold-conditioned molecular generation',
         unit='Molecule', mode='slot',
         funnel=['Requested molecules', 'Valid and deduplicated', 'Properties QC Passed', 'Sent to docking']),
    dict(n=4, key='docking', title='Molecular docking', short='Docking',
         subtitle='Ensemble docking and scoring across three complex conformations',
         unit='Molecule', mode='slot',
         funnel=['Input molecules', 'Docking succeeded', 'Score threshold passed', 'Three-conformation agreement']),
    dict(n=5, key='screening', title='Virtual screening', short='Screening',
         subtitle='Binding mode, drug-likeness and diversity screening',
         unit='Molecule', mode='slot',
         funnel=['Input molecules', 'Key salt bridges retained', 'Drug-likeness / PAINS', 'Diversity clustering']),
    dict(n=6, key='md', title='Molecular dynamics', short='MD',
         subtitle='Complex stability and binding free energy',
         unit='Molecule', mode='slot',
         funnel=['Entered MD', 'Trajectory complete', 'Ligand RMSD Stable', 'Salt bridge retained']),
    dict(n=7, key='final_hits', title='Results overview', short='Result',
         subtitle='Candidates and workflow funnel',
         unit='Molecule', mode='slot',
         funnel=['Final candidates']),
]
BY_KEY = {s['key']: s for s in STAGES}

CONTRACT = {
    'generated': dict(
        csv='results/generated.csv', summary='results/generated_summary.json',
        poses='results/generated_poses.sdf',
        desc='Generative model output:id, smiles, parent_scaffold, …'),
    'docking': dict(
        csv='results/docking.csv', summary='results/docking_summary.json',
        poses='results/docking_poses.sdf',
        desc='Docking results:id, smiles, vina_score, cnn_score, keep, …'),
    'screening': dict(
        csv='results/screening.csv', summary='results/screening_summary.json',
        poses='results/screening_poses.sdf',
        desc='Screening results:id, smiles, keep, reason, …'),
    'md': dict(
        csv='results/md.csv', summary='results/md_summary.json',
        poses='results/md_poses.sdf',
        desc='MD Summary:id, smiles, ns, rmsd_lig, saltbridge_frac, dG, keep, …'),
    'final_hits': dict(
        csv='results/final_hits.csv', summary='results/final_hits_summary.json',
        poses='results/final_hits_poses.sdf',
        desc='Final candidates:id, smiles, rank, …'),
}


# --------------------------------------------------------------------------
def _drug_like_poses():
    return data.read_json('results/drug_like_poses.json', {}) or {}


def _stage2_rows():
    """All 538 drug-like scaffolds with their structure-filter verdict."""
    def produce():
        rows = data.filter_rows()
        poses = _drug_like_poses()
        out = []
        for i, r in enumerate(rows):
            key = f"{r.get('smiles')}|{r.get('axis')}"
            out.append(dict(
                id=f'S2-{i + 1:03d}', smiles=r.get('smiles'), axis=r.get('axis'),
                score=r.get('total'), heavy=r.get('heavy_atoms'),
                anchors=r.get('anchor_ids'), bridge=r.get('bridge'),
                spread=r.get('ens_spread'), keep=bool(r.get('keep')),
                has_pose=key in poses, _key=key))
        out.sort(key=lambda x: (not x['keep'], -(x['score'] or 0)))
        for i, r in enumerate(out):
            r['id'] = f'S2-{i + 1:03d}'
        return out
    return data._memo(('stage2',), (data._mtime(data.rel('results/structure_filters.csv')),
                                    data._mtime(data.rel('results/drug_like_poses.json'))), produce)


def _slot_rows(key):
    c = CONTRACT[key]
    rows = data.read_csv(c['csv'])
    out = []
    for i, r in enumerate(rows):
        row = {}
        for k, v in r.items():
            if v is None or v == '':
                row[k] = None
                continue
            row[k] = data._num(v)
        row.setdefault('id', f'{key[:3].upper()}-{i + 1:04d}')
        if 'keep' in row:
            row['keep'] = bool(row['keep']) if not isinstance(row['keep'], str) else row['keep'].lower() == 'true'
        out.append(row)
    return out


def _slot_summary(key):
    return data.read_json(CONTRACT[key]['summary'], None)


def _slot_poses(key):
    path = CONTRACT[key]['poses']
    return data._memo(('slotposes', key), (data._mtime(data.rel(path)),),
                      lambda: data._parse_sdf(path))


# --------------------------------------------------------------------------
# funnels
# --------------------------------------------------------------------------
def _f(label, count, note=None, tone=None):
    return dict(label=label, count=count, note=note, tone=tone)


def _funnel_complex():
    hot = data.hotpsots_raw()
    ph = data.pharmacophore() or {}
    pk = data.pocket() or {}
    common = sum(1 for r in hot if float(r.get('convergence') or 0) >= 1.0)
    return [
        _f('Interface residues', len(hot), 'Peptide contact interfaces across three complexes'),
        _f('Convergent hotspots across three structures', common, 'convergence ≥ 1.0'),
        _f('Pharmacophore features', len(ph.get('features', [])), 'Anchor / Aromatic / Hydrophobic / Donor'),
        _f('Available pharmacophore axes', len(data.axes()), 'Pairwise feature axes'),
        _f('Pocket subregions', len(pk.get('subpockets', [])), 'S1–S4'),
    ]


def _funnel_scaffold():
    s = data.summary() or {}
    f = s.get('funnel', {})
    fs = data.filter_summary() or {}
    return [
        _f('Enumerated scaffolds', f.get('enumerated'), 'Head / Ring / Linker / Tail combinations'),
        _f('Property checks passed', f.get('property_pass'), 'MW, cLogP, TPSA etc.'),
        _f('Three drug-likeness checks passed', f.get('drug_like'), 'No clashes + Peptide mimicry + Properties'),
        _f('Structure filters passed', fs.get('survivors'), 'Salt bridge / Repulsion / Bridging / Ensemble'),
    ]


def _extras_scaffold():
    s = data.summary() or {}
    f = s.get('funnel', {})
    fs = data.filter_summary() or {}
    return dict(
        independent=[
            dict(label='No clashes', count=f.get('no_clash')),
            dict(label='Peptide overlap', count=f.get('peptide_mimic')),
            dict(label='Anchor salt bridge', count=(fs.get('by_stage') or {}).get('anchor')),
            dict(label='Subregion bridging', count=(fs.get('by_stage') or {}).get('bridge')),
            dict(label='Ensemble consistency', count=(fs.get('by_stage') or {}).get('ensemble')),
        ])


def _demo_counts():
    fs = data.filter_summary() or {}
    by = fs.get('by_stage') or {}
    n_in = fs.get('n_input')
    return dict(n_in=n_in, bridge=by.get('bridge'), ens=by.get('ensemble'))


def _demo_rows(key):
    rows = _stage2_rows()
    if key in ('generated', 'docking'):
        return rows
    if key == 'screening':
        return [r for r in rows if r['bridge']]
    return [r for r in rows if r['keep']]


def _demo_funnel(key):
    d = _demo_counts()
    labels = BY_KEY[key]['funnel']
    n_all, n_br, n_k = d['n_in'], d['bridge'], d['ens']
    table = {
        'generated': [n_all, n_all, n_all, n_all],
        'docking': [n_all, n_all, n_br, n_k],
        'screening': [n_br, n_br, n_k, n_k],
        'md': [n_k, n_k, n_k, n_k],
        'final_hits': [n_k],
    }[key]
    return [_f(lb, c) for lb, c in zip(labels, table)]


# --------------------------------------------------------------------------
def stage_state(key, demo=False):
    """Return dict(status, rows, funnel, files) for one stage."""
    st = BY_KEY[key]
    if st['mode'] == 'builtin':
        rows = _stage2_rows() if key == 'scaffold' else []
        available = bool(rows) if key == 'scaffold' else any(data.exists(s['pdb']) for s in data.STRUCTURES.values())
        return dict(status='historical' if available else 'pending', rows=rows,
                    funnel=_funnel_scaffold() if key == 'scaffold' else _funnel_complex(),
                    files=[])
    c = CONTRACT[key]
    files = [dict(path=c['csv'], role='Molecule table', exists=data.exists(c['csv'])),
             dict(path=c['summary'], role='Funnel counts', exists=data.exists(c['summary'])),
             dict(path=c['poses'], role='3D Poses (optional)', exists=data.exists(c['poses']))]
    real = _slot_rows(key) if data.exists(c['csv']) else []
    summ = _slot_summary(key)
    if real or summ:
        if summ and summ.get('funnel'):
            funnel = [_f(x.get('label'), x.get('count'), x.get('note')) for x in summ['funnel']]
        else:
            funnel = [_f('Molecule count', len(real))]
        return dict(status='ready', rows=real, funnel=funnel, files=files)
    if demo:
        return dict(status='demo', rows=_demo_rows(key), funnel=_demo_funnel(key), files=files)
    return dict(status='pending', rows=[],
                funnel=[_f(lb, None) for lb in st['funnel']], files=files)


def workflow(demo=False):
    out = []
    for st in STAGES:
        state = stage_state(st['key'], demo)
        funnel = state['funnel']
        last = funnel[-1]['count'] if funnel else None
        if st['key'] == 'complex':
            last = sum(data.exists(s['pdb']) for s in data.STRUCTURES.values())
        out.append(dict(
            n=st['n'], key=st['key'], title=st['title'], short=st['short'],
            subtitle=st['subtitle'], unit=st['unit'], status=state['status'],
            funnel=funnel, out=last, n_rows=len(state['rows']),
            files=state['files'], contract=CONTRACT.get(st['key'], {}).get('desc'),
            extras=_extras_scaffold() if st['key'] == 'scaffold' else None))
    return dict(stages=out, demo=demo)


_PUBLIC = ('id', 'smiles', 'axis', 'score', 'heavy', 'anchors', 'bridge', 'spread', 'keep', 'has_pose')


def molecules(key, demo=False):
    state = stage_state(key, demo)
    rows = state['rows']
    if key in ('scaffold',) or state['status'] == 'demo':
        rows = [{k: r.get(k) for k in _PUBLIC} for r in rows]
        cols = [('id', 'ID'), ('smiles', 'SMILES'), ('axis', 'Pharmacophore axis'),
                ('score', 'Pose score'), ('spread', 'Ensemble dispersion Å'), ('keep', 'Passed')]
    else:
        cols = []
        if rows:
            keys = [k for k in rows[0].keys() if not k.startswith('_')]
            lead = [k for k in ('id', 'smiles') if k in keys]
            rest = [k for k in keys if k not in lead][:6]
            cols = [(k, k) for k in lead + rest]
    return dict(status=state['status'], rows=rows, columns=[dict(key=k, label=lb) for k, lb in cols],
                n=len(rows))


def _lig(atoms, bonds, smiles, legend):
    return dict(smiles=smiles, atoms=atoms, bonds=bonds, legend=legend)


def pose(key, mol_id, demo=False):
    """3D ligand for one molecule, in the viewer's scene format (or None)."""
    state = stage_state(key, demo)
    row = next((r for r in state['rows'] if str(r.get('id')) == str(mol_id)), None)
    if row is None:
        return dict(ligand=None, reason='Molecule not found')
    if key == 'scaffold' or state['status'] == 'demo':
        s2 = next((r for r in _stage2_rows() if r['id'] == mol_id), None)
        p = _drug_like_poses().get(s2['_key']) if s2 else None
        if not p:
            return dict(ligand=None, reason='None 3D Pose')
        atoms = [dict(e=e, p=xyz, rgb=ELEMENT.get(e.upper(), '#8fa3bf'))
                 for e, xyz in zip(p['el'], p['xyz'])]
        return dict(ligand=_lig(atoms, p['bonds'], p['smiles'], f"{mol_id} · {p['axis']}"))
    for m in _slot_poses(key):
        if m['smiles'] in (str(mol_id), str(row.get('smiles'))):
            return dict(ligand=m)
    return dict(ligand=None, reason='This molecule has no 3D pose yet (see file contract)')


def mol2d(smiles):
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except Exception:                                       # noqa: BLE001
        return dict(svg=None, reason='Not installed RDKit, Cannot draw 2D Structure')
    m = Chem.MolFromSmiles(smiles or '')
    if m is None:
        return dict(svg=None, reason='SMILES could not be parsed')
    d = rdMolDraw2D.MolDraw2DSVG(360, 240)
    d.drawOptions().clearBackground = False
    d.drawOptions().bondLineWidth = 2
    d.DrawMolecule(m)
    d.FinishDrawing()
    return dict(svg=d.GetDrawingText())
