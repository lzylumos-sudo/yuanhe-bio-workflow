# -*- coding: utf-8 -*-
"""webui.data -- read-only access layer for the WEDGE-Design web console.

Every number served by the API is read from a file produced by the pipeline
(`results/`, `evidence/`) or recomputed from the deposited coordinates.  Nothing
is hard-coded in the UI layer, so re-running the pipeline is reflected after a
page refresh (files are cached by modification time).

Only the Python standard library is used here; `pepwedge` is imported lazily so
that the server still starts if numpy is missing.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEBUI = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(WEBUI, 'cache')

STRUCTURES = {
    'SglM': dict(pdb='9nu4.pdb', pdbid='9NU4', title='MurJ in complex with SglM',
                 resolution=3.60, peptide='SglM', color='#ecc94b'),
    'SglPP7': dict(pdb='9NU5.pdb', pdbid='9NU5', title='MurJ in complex with SglPP7',
                   resolution=3.70, peptide='SglPP7', color='#dd6b20'),
    'SglCJ3': dict(pdb='9NU8.pdb', pdbid='9NU8', title='MurJ in complex with SglCJ3',
                   resolution=3.60, peptide='SglCJ3', color='#38a169'),
}
REF_TAG = 'SglCJ3'          # frame in which poses / pharmacophore / pocket live
MURJ_RANGE = (5, 500)       # EcMurJ numbering; everything else is fusion construct

ELEMENT_COLORS = {'C': '#8fa3bf', 'N': '#4f8ef7', 'O': '#e55353', 'S': '#e0b64a',
                  'P': '#d97757', 'F': '#7ad1a0', 'CL': '#7ad1a0', 'BR': '#b06a4a',
                  'I': '#8a6fd1', 'H': '#c8d2df'}
CHEM_COLORS = {'pos': '#4f8ef7', 'neg': '#e55353', 'polar': '#38a169',
               'aromatic': '#c99a3a', 'hydrophobic': '#8b97a8'}
SUBPOCKET_COLORS = {'S1_basic_bowl': '#4f8ef7', 'S2_carboxylate_gate': '#e55353',
                    'S3_TM2_TM7_groove': '#38a169', 'S4_TM7_floor': '#c99a3a'}

# --------------------------------------------------------------------------
# small file helpers with mtime-keyed memoisation
# --------------------------------------------------------------------------
_lock = threading.RLock()
_cache: dict = {}


def _mtime(path):
    try:
        return os.path.getmtime(path)
    except OSError:
        return 0.0


def _memo(key, mtimes, producer):
    stamp = tuple(round(m, 3) for m in mtimes)
    with _lock:
        hit = _cache.get(key)
        if hit and hit[0] == stamp:
            return hit[1]
    value = producer()
    with _lock:
        _cache[key] = (stamp, value)
    return value


def rel(path):
    local = os.path.join(ROOT, path)
    if path in {s['pdb'] for s in STRUCTURES.values()} and not os.path.isfile(local):
        return os.path.join(os.path.dirname(ROOT), 'data', 'structures', path.lower())
    return local


def exists(path):
    return os.path.exists(rel(path))


def read_json(path, default=None):
    p = rel(path)

    def produce():
        if not os.path.exists(p):
            return default
        with open(p, encoding='utf-8') as fh:
            return json.load(fh)

    return _memo(('json', path), (_mtime(p),), produce)


def read_text(path, default=''):
    p = rel(path)

    def produce():
        if not os.path.exists(p):
            return default
        with open(p, encoding='utf-8', errors='replace') as fh:
            return fh.read()

    return _memo(('text', path), (_mtime(p),), produce)


def _num(value):
    if value is None or value == '':
        return None
    if isinstance(value, (int, float)):
        return value
    txt = str(value).strip()
    low = txt.lower()
    if low in ('true', 'false'):
        return low == 'true'
    try:
        if txt.lstrip('+-').isdigit():
            return int(txt)
        return float(txt)
    except ValueError:
        return value


def read_csv(path, numeric=None):
    p = rel(path)

    def produce():
        if not os.path.exists(p):
            return []
        numeric_set = set(numeric or ())
        with open(p, encoding='utf-8', newline='') as fh:
            rows = []
            for row in csv.DictReader(fh):
                out = {}
                for k, v in row.items():
                    out[k] = _num(v) if k in numeric_set else v
                rows.append(out)
            return rows

    return _memo(('csv', path), (_mtime(p),), produce)


# --------------------------------------------------------------------------
# results
# --------------------------------------------------------------------------
NUM_SCAFFOLD = {'rank', 'axis_length_A', 'MW', 'cLogP', 'TPSA', 'HBD', 'HBA', 'RotB',
                'Rings', 'heavy_atoms', 'pose_score', 'total', 'clash', 'salt', 'hb',
                'hyd', 'peptide_overlap', 'buried', 'property_penalty', 'rigidity'}
NUM_FILTER = {'total', 'heavy_atoms', 'n_anchor', 'n_repulsion', 'n_polar', 'n_hydro',
              'ens_rmsd', 'ens_spread', 'ens_min_score', 'ens_max_clash', 'score_SglM',
              'score_SglPP7', 'score_SglCJ3', 'keep', 'bridge'}


def summary():
    return read_json('results/summary.json')


def hotpsots_raw():
    return read_csv('results/interface_hotspots.csv',
                    numeric={'murj_resid', 'score', 'convergence', 'dSASA_A2',
                             'n_contacts', 'n_peptide_residues', 'dmin_A'})


def pharmacophore():
    return read_json('results/pharmacophore.json', {'features': [], 'pairs': []})


def pocket():
    return read_json('results/pocket.json')


def axes():
    return read_json('results/axes.json', [])


def filter_summary():
    return read_json('results/structure_filters_summary.json')


def filter_rows():
    return read_csv('results/structure_filters.csv', numeric=NUM_FILTER)


def scaffold_rows():
    """Ranked scaffolds joined with the structure-filter verdict (by SMILES)."""
    def produce():
        rows = read_csv('results/scaffolds_ranked.csv', numeric=NUM_SCAFFOLD)
        filt = {r.get('smiles'): r for r in filter_rows()}
        poses = pose_index()
        for r in rows:
            f = filt.get(r.get('smiles'))
            if f:
                r['f_anchor'] = f.get('n_anchor')
                r['f_anchor_ids'] = f.get('anchor_ids')
                r['f_repulsion'] = f.get('n_repulsion')
                r['f_bridge'] = f.get('bridge')
                r['f_ens_rmsd'] = f.get('ens_rmsd')
                r['f_ens_spread'] = f.get('ens_spread')
                r['f_keep'] = f.get('keep')
                r['f_scores'] = {k[6:]: f[k] for k in f if k.startswith('score_')}
            r['pose_index'] = poses.get(r.get('smiles'))
        return rows

    return _memo(('scaffolds',), (_mtime(rel('results/scaffolds_ranked.csv')),
                                  _mtime(rel('results/structure_filters.csv'))), produce)


# --------------------------------------------------------------------------
# poses (SDF)
# --------------------------------------------------------------------------
def _parse_sdf(path):
    p = rel(path)
    if not os.path.exists(p):
        return []
    mols = []
    with open(p, encoding='utf-8', errors='replace') as fh:
        text = fh.read()
    for block in text.split('$$$$'):
        lines = block.splitlines()
        while lines and not lines[0].strip():
            lines.pop(0)
        if len(lines) < 4:
            continue
        try:
            n_atoms = int(lines[2][0:3])
        except (ValueError, IndexError):
            continue
        header = lines[0].strip()
        atoms, bonds = [], []
        for line in lines[3:3 + n_atoms]:
            if len(line) < 31:
                continue
            try:
                x, y, z = float(line[0:10]), float(line[10:20]), float(line[20:30])
            except ValueError:
                continue
            el = line[31:34].strip().upper() or 'C'
            atoms.append(dict(e=el.capitalize() if el != 'CL' else 'Cl',
                              p=[round(x, 3), round(y, 3), round(z, 3)],
                              rgb=ELEMENT_COLORS.get(el, '#8fa3bf')))
        for line in lines[3 + n_atoms:]:
            if line.startswith('M  END'):
                break
            if len(line) < 9:
                continue
            try:
                a, b = int(line[0:3]) - 1, int(line[3:6]) - 1
            except ValueError:
                continue
            if 0 <= a < len(atoms) and 0 <= b < len(atoms):
                bonds.append([a, b])
        if atoms:
            mols.append(dict(smiles=header, atoms=atoms, bonds=bonds,
                             legend=lines[1].strip() if len(lines) > 1 else ''))
    return mols


def poses_top():
    return _memo(('sdf', 'results/top_scaffold_poses.sdf'),
                 (_mtime(rel('results/top_scaffold_poses.sdf')),),
                 lambda: _parse_sdf('results/top_scaffold_poses.sdf'))


def poses_diverse():
    return _memo(('sdf', 'results/diverse_scaffolds.sdf'),
                 (_mtime(rel('results/diverse_scaffolds.sdf')),),
                 lambda: _parse_sdf('results/diverse_scaffolds.sdf'))


def pose_index():
    return _memo(('poseidx',), (_mtime(rel('results/top_scaffold_poses.sdf')),),
                 lambda: {m['smiles']: i for i, m in enumerate(poses_top())})


# --------------------------------------------------------------------------
# evidence
# --------------------------------------------------------------------------
def evidence_report():
    return read_text('evidence/report.md')


def evidence_file(path):
    return read_text(path)


def priors_check():
    return read_csv('evidence/04_priors_check.csv')


def census():
    return read_csv('evidence/00_census.csv', numeric={'construct_residues', 'construct_atoms',
                                                       'murj_residues', 'murj_atoms',
                                                       'peptide_residues', 'peptide_atoms',
                                                       'contacts_heavy_atoms'})


def interface_murj():
    return read_csv('evidence/01_interface_murj.csv',
                    numeric={'resid', 'n_atom_contacts', 'n_peptide_residues', 'dmin_A',
                             'dSASA_A2', 'convergence'})


def interface_peptide():
    return read_csv('evidence/02_interface_peptide.csv',
                    numeric={'pep_resid', 'n_atom_contacts', 'n_murj_residues', 'dmin_A',
                             'dSASA_A2', 'heavy_atoms'})


def anchors():
    return read_csv('evidence/03_anchors.csv',
                    numeric={'sidechain_centroid_A', 'min_atom_A'})


def atom_budget():
    return read_json('evidence/05_atom_budget.json', {})


def priors_resolved():
    return read_json('evidence/06_priors_resolved.json', {})


def evidence_status():
    rows = priors_check()
    counts: dict = {}
    for r in rows:
        key = (r.get('status') or 'UNKNOWN').strip() or 'UNKNOWN'
        counts[key] = counts.get(key, 0) + 1
    return dict(n=len(rows), by_status=counts,
                generated=_mtime(rel('evidence/report.md')))


# --------------------------------------------------------------------------
# docs
# --------------------------------------------------------------------------
def docs_list():
    out = []
    for folder in ('docs', '.'):
        base = rel(folder)
        if not os.path.isdir(base):
            continue
        for name in sorted(os.listdir(base)):
            if not name.endswith('.md'):
                continue
            path = os.path.join(base, name)
            stat = os.stat(path)
            out.append(dict(name=name, folder=folder, size=stat.st_size,
                            mtime=stat.st_mtime, title=_md_title(path) or name))
    return out


def _md_title(path):
    try:
        with open(path, encoding='utf-8') as fh:
            for line in fh:
                if line.startswith('# '):
                    return line[2:].strip()
    except OSError:
        pass
    return ''


def doc(name):
    safe = os.path.basename(name)
    for folder in ('docs', '.'):
        path = rel(os.path.join(folder, safe))
        if os.path.exists(path):
            return read_text(os.path.join(folder, safe))
    return ''


# --------------------------------------------------------------------------
# structures / 3D scene
# --------------------------------------------------------------------------
def _pepwedge():
    import sys
    if ROOT not in sys.path:
        sys.path.insert(0, ROOT)
    from pepwedge import core as C          # noqa: F401
    return C


def _superpose():
    import sys
    if ROOT not in sys.path:
        sys.path.insert(0, ROOT)
    from pepwedge.pharmacophore import superpose_complex
    return superpose_complex


BOND_MAX = 1.95


def _bonds(atoms, res_groups):
    """Distance-inferred bonds inside each residue group (list of index lists)."""
    bonds = []
    for idxs in res_groups:
        for i in range(len(idxs)):
            for j in range(i + 1, len(idxs)):
                a, b = atoms[idxs[i]]['p'], atoms[idxs[j]]['p']
                if (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 < BOND_MAX ** 2:
                    bonds.append([idxs[i], idxs[j]])
    return bonds


def _load_complex(tag):
    C = _pepwedge()
    info = STRUCTURES[tag]
    path = rel(info['pdb'])
    return _memo(('complex', tag), (_mtime(path),),
                 lambda: C.load_complex(path, tag, residue_range=MURJ_RANGE))


def _transform(xyz, rt):
    import numpy as np
    R, t = rt
    arr = np.asarray(xyz, dtype=float)
    return arr @ R.T + t


def scene(ref='SglCJ3', ligand=None, radius=13.0):
    """All geometry the 3D viewer needs, expressed in the `ref` structure frame."""
    import numpy as np
    C = _pepwedge()
    superpose_complex = _superpose()
    if ref not in STRUCTURES:
        ref = REF_TAG
    cxs = {tag: _load_complex(tag) for tag in STRUCTURES}
    ref_cx = cxs[ref]

    transforms = {}
    for tag, cx in cxs.items():
        if tag == ref:
            transforms[tag] = (np.eye(3), np.zeros(3))
        else:
            transforms[tag] = superpose_complex(cx, ref_cx)[0]

    # ---- receptor: C-alpha trace + pocket-lining side chains
    rec = ref_cx.receptor
    pk_rep = pocket() or {}

    # map each residue onto the sub-pocket that lists it (used for colouring)
    sub_of = {}
    sub_meta = []
    for sp in pk_rep.get('subpockets', []):
        name = sp.get('subpocket', '')
        for label in sp.get('residues', []):
            digits = ''.join(ch for ch in str(label) if ch.isdigit())
            if digits:
                sub_of[int(digits)] = name
        sub_meta.append(dict(name=name, volume=sp.get('volume_A3'),
                             buriedness=sp.get('mean_buriedness'),
                             centroid=sp.get('centroid'),
                             residues=sp.get('residues', []),
                             color=SUBPOCKET_COLORS.get(name, '#8fa3bf')))

    ca_trace, rec_residues, atoms, groups = [], [], [], []
    for rid in sorted(rec):
        res = rec[rid]
        ca = res.atoms.get('CA')
        if ca is None:
            continue
        ca_trace.append([round(float(c), 3) for c in ca])
        rec_residues.append(dict(resid=int(rid), one=res.one, name=res.name,
                                 sub=sub_of.get(int(rid)),
                                 ca=[round(float(c), 3) for c in ca]))

    # centre of interest: ligand pose if given, else the buried cavity centroid
    center = None
    lig = None
    top = poses_top()
    if ligand is not None:
        try:
            idx = int(ligand)
            if 0 <= idx < len(top):
                lig = top[idx]
        except (TypeError, ValueError):
            pass
    if lig is not None:
        pts = np.array([a['p'] for a in lig['atoms']], dtype=float)
        center = pts.mean(0)
    elif pk_rep.get('subpockets'):
        center = np.array(pk_rep['subpockets'][0]['centroid'], dtype=float)
    else:
        center = np.array(ca_trace, dtype=float).mean(0)

    for rid in sorted(rec):
        res = rec[rid]
        sel = [an for an in sorted(res.atoms)
               if float(np.linalg.norm(res.atoms[an] - center)) <= radius]
        if not sel:
            continue
        idxs = []
        for an in sel:
            el = an[0].upper() if an[0].isalpha() else 'C'
            el = {'CL': 'Cl', 'BR': 'Br'}.get(el, el)
            idxs.append(len(atoms))
            atoms.append(dict(e=el, p=[round(float(c), 3) for c in res.atoms[an]],
                              a=an, resid=int(rid), one=res.one,
                              rgb=ELEMENT_COLORS.get(el.upper(), '#8fa3bf')))
        groups.append(idxs)
    side_bonds = _bonds(atoms, groups)

    # ---- peptides (all three, superposed into the ref frame)
    peptides = []
    for tag, cx in cxs.items():
        rt = transforms[tag]
        res_atoms, res_groups, ca = [], [], []
        for lid in sorted(cx.ligand):
            res = cx.ligand[lid]
            idxs = []
            for an in sorted(res.atoms):
                xyz = _transform(res.atoms[an], rt)
                el = an[0].upper() if an[0].isalpha() else 'C'
                el = {'CL': 'Cl', 'BR': 'Br'}.get(el, el)
                idxs.append(len(res_atoms))
                res_atoms.append(dict(e=el, p=[round(float(c), 3) for c in xyz], a=an,
                                      resid=int(lid), one=res.one, tag=tag,
                                      rgb=ELEMENT_COLORS.get(el.upper(), '#8fa3bf')))
            res_groups.append(idxs)
            if 'CA' in res.atoms:
                ca.append([round(float(c), 3) for c in _transform(res.atoms['CA'], rt)])
        peptides.append(dict(tag=tag, pdbid=STRUCTURES[tag]['pdbid'],
                             color=STRUCTURES[tag]['color'], ca=ca,
                             atoms=res_atoms, bonds=_bonds(res_atoms, res_groups),
                             n_residues=len(cx.ligand), seq=cx.seq('ligand')))

    # ---- pharmacophore features in the ref frame
    ph = pharmacophore() or {}
    rt_cj3 = transforms[REF_TAG]
    features = []
    for f in ph.get('features', []):
        pos = _transform(f['pos'], rt_cj3)
        features.append(dict(fid=f['fid'], type=f['type'],
                             p=[round(float(c), 2) for c in pos],
                             weight=f.get('weight'), tol=f.get('tol'),
                             convergence=f.get('convergence'),
                             sources=f.get('sources', []),
                             contacts=f.get('murj_contacts', []),
                             rgb={'anion': '#e55353', 'cation': '#4f8ef7',
                                  'aromatic': '#c99a3a', 'donor': '#38a169',
                                  'hydrophobic': '#8b97a8'}.get(f['type'], '#c8d2df')))

    # ---- hotspots (validated residues) on the ref receptor
    hot = hotpsots_raw()[:24]
    hotspots = []
    for h in hot:
        rid = int(h['murj_resid'])
        res = ref_cx.receptor.get(rid)
        if res is None or 'CA' not in res.atoms:
            continue
        hotspots.append(dict(resid=rid, one=res.one, score=h['score'],
                             dSASA=h['dSASA_A2'], evidence=h['evidence'],
                             p=[round(float(c), 3) for c in res.atoms['CA']]))

    return dict(ref=ref, ref_pdbid=STRUCTURES[ref]['pdbid'],
                structure=dict(STRUCTURES[ref]),
                center=[round(float(c), 3) for c in center],
                radius=radius,
                receptor=dict(ca=ca_trace, residues=rec_residues,
                              atoms=atoms, bonds=side_bonds),
                peptides=peptides,
                ligand=lig,
                features=features,
                hotspots=hotspots,
                pocket=pk_rep,
                subpockets=sub_meta,
                structures=[dict(tag=t, pdbid=STRUCTURES[t]['pdbid'],
                                 resolution=STRUCTURES[t]['resolution'],
                                 color=STRUCTURES[t]['color'],
                                 title=STRUCTURES[t]['title']) for t in STRUCTURES])


# --------------------------------------------------------------------------
# pocket voxels (expensive -> cached on disk, computed in a worker thread)
# --------------------------------------------------------------------------
_pocket_state = dict(status='idle', error=None, mtime=0.0)
_pocket_lock = threading.Lock()


def _pocket_cache_path(step):
    stamp = hashlib.sha1(f'{_mtime(rel(STRUCTURES[REF_TAG]["pdb"])):.3f}|{step}'
                         .encode()).hexdigest()[:12]
    return os.path.join(CACHE, f'pocket_{stamp}.json')


def pocket_points(step=0.8):
    """Return the buried cavity voxel cloud, computing it once if needed."""
    path = _pocket_cache_path(step)
    if os.path.exists(path):
        try:
            with open(path, encoding='utf-8') as fh:
                data = json.load(fh)
            _pocket_state.update(status='ready', error=None)
            return dict(status='ready', **data)
        except (OSError, ValueError) as exc:
            _pocket_state.update(status='error', error=str(exc))

    with _pocket_lock:
        if _pocket_state['status'] != 'computing':
            _pocket_state.update(status='computing', error=None)

            def work():
                try:
                    os.makedirs(CACHE, exist_ok=True)
                    from pepwedge.pocket import build_pocket
                    cx = _load_complex(REF_TAG)
                    pk = build_pocket(cx, step=step)
                    keep = pk.buriedness >= 0.5
                    pts = pk.points[keep]
                    chem = [c for c, k in zip(pk.chemistry, keep) if k]
                    sub = [s for s, k in zip(pk.subpocket, keep) if k]
                    subs = sorted(set(sub))
                    chem_order = ['pos', 'neg', 'polar', 'aromatic', 'hydrophobic']
                    chem_idx = {c: i for i, c in enumerate(chem_order)}
                    sub_idx = {s: i for i, s in enumerate(subs)}
                    payload = dict(
                        step=step,
                        volume=round(float(len(pts) * step ** 3), 1),
                        n_voxels=int(len(pts)),
                        xyz=[round(float(v), 2) for p in pts for v in p],
                        chem=[chem_idx.get(c, 4) for c in chem],
                        sub=[[sub_idx.get(s, 0) for s in sub], subs],
                    )
                    tmp = path + '.tmp'
                    with open(tmp, 'w', encoding='utf-8') as fh:
                        json.dump(payload, fh)
                    os.replace(tmp, path)
                    _pocket_state.update(status='ready')
                except Exception as exc:                     # noqa: BLE001
                    _pocket_state.update(status='error', error=f'{type(exc).__name__}: {exc}')

            threading.Thread(target=work, name='pocket-voxels', daemon=True).start()
    return dict(status=_pocket_state['status'], error=_pocket_state['error'],
                step=step, n_voxels=0, xyz=[], chem=[], sub=[[], []])


def pocket_status():
    return dict(_pocket_state)


# --------------------------------------------------------------------------
# meta
# --------------------------------------------------------------------------
def meta():
    import platform
    import sys
    files = {
        'results/summary.json': _mtime(rel('results/summary.json')),
        'results/scaffolds_ranked.csv': _mtime(rel('results/scaffolds_ranked.csv')),
        'results/structure_filters.csv': _mtime(rel('results/structure_filters.csv')),
        'evidence/report.md': _mtime(rel('evidence/report.md')),
        'figures/fig4_funnel.png': _mtime(rel('figures/fig4_funnel.png')),
    }
    figs = []
    fdir = rel('figures')
    if os.path.isdir(fdir):
        figs = [dict(name=n, size=os.path.getsize(os.path.join(fdir, n)))
                for n in sorted(os.listdir(fdir)) if n.lower().endswith('.png')]
    return dict(root=ROOT, python=sys.version.split()[0],
                platform=f'{platform.system()} {platform.release()}',
                started=time.time(),
                files=files, figures=figs,
                available=dict(results=exists('results/summary.json'),
                               evidence=exists('evidence/report.md'),
                               filters=exists('results/structure_filters.csv'),
                               tests=exists('tests/test_smoke.py')),
                structures=STRUCTURES, ref=REF_TAG)
