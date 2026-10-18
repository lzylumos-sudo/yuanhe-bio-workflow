"""Transfer a hashed PPM placement only between rigidly matching constructs."""
from pathlib import Path
import json
import numpy as np
from ..ids import sha256_file
from ..util import parse_pdb_text


def atom_map(path):
    atoms = parse_pdb_text(Path(path).read_text(), keep_hydrogens=False).atoms
    result = {}
    for a in atoms:
        if a.record.strip() != 'ATOM' or a.is_hydrogen():
            continue
        key = (a.chain, a.resseq, getattr(a, 'icode', ''), a.name, a.resname)
        if key in result:
            raise ValueError('Ambiguous atom mapping in orientation evidence')
        result[key] = [a.x, a.y, a.z]
    return result


def fit(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 2 or x.shape[1] != 3 or len(x) < 12 or not np.isfinite([x,y]).all():
        raise ValueError('Invalid orientation coordinates')
    xc, yc = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x-xc).T @ (y-yc))
    correction = np.eye(3); correction[2,2] = np.linalg.det(u @ vt)
    rot = u @ correction @ vt
    shift = yc - xc @ rot
    rmsd = float(np.sqrt(np.mean(np.sum((x @ rot + shift-y)**2, axis=1))))
    return rot, shift, rmsd


def transfer(evidence_dir, receptor, *, min_coverage=0.95, max_rmsd_A=0.05):
    folder = Path(evidence_dir)
    meta = json.loads((folder/'orientation.json').read_text())
    for name, key in [('murj.pdb','input_sha256'), ('murjout.pdb','output_sha256')]:
        if sha256_file(folder/name) != meta[key]:
            raise ValueError('PPM orientation evidence checksum mismatch')
    if not meta.get('source') or not meta.get('method'):
        raise ValueError('Orientation source and method are required')
    half = float(meta['membrane_half_thickness_A'])
    if not np.isfinite(half) or half <= 0:
        raise ValueError('Invalid membrane thickness')
    source, placed, current = atom_map(folder/'murj.pdb'), atom_map(folder/'murjout.pdb'), atom_map(receptor)
    shared = sorted(source.keys() & placed.keys() & current.keys())
    coverage = len(shared)/max(len(source), len(current), 1)
    if coverage < min_coverage:
        raise ValueError('Receptor construct does not match orientation reference')
    _, _, structure_rmsd = fit([current[k] for k in shared], [source[k] for k in shared])
    rot, shift, ppm_rmsd = fit([current[k] for k in shared], [placed[k] for k in shared])
    if max(structure_rmsd, ppm_rmsd) > max_rmsd_A:
        raise ValueError('Orientation transfer exceeds rigid-fit tolerance')
    return {'center': (-shift @ rot.T).tolist(), 'normal': rot[:,2].tolist(),
            'thickness_A': 2*half, 'source': meta['source'], 'method': meta['method'],
            'evidence_kind': 'computational_PPM_placement', 'experimental_orientation': False,
            'coordinate_frame': 'exact_docking_receptor', 'receptor_sha256': sha256_file(receptor),
            'matched_atoms': len(shared), 'coverage': coverage, 'rigid_fit_rmsd_A': ppm_rmsd,
            'limits': {'min_coverage': min_coverage, 'max_rmsd_A': max_rmsd_A},
            'evidence_hashes': {p.name:sha256_file(p) for p in sorted(folder.iterdir()) if p.is_file()},
            'limitations': meta.get('limitations', [])}
