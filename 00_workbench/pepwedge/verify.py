# -*- coding: utf-8 -*-
"""pepwedge.verify -- check every literature prior against the coordinates.

The point is not to "confirm" the paper. The point is to make disagreement
visible. Each prior carries a machine-checkable specification; this module runs
it and returns PASS / FAIL / FLAG together with the observed numbers, so a
reader can see exactly which design assumption is backed by the deposited
structures and which one is only backed by text.

Statuses
    PASS          prior reproduced in the coordinates
    FAIL          prior not reproduced -- investigate before relying on it
    FLAG          not reproduced, and that was expected (see the prior's note)
    WEAK          reproduced only partially / ambiguous
    NO_RESIDUE    residue absent from the deposited model
    MISMATCH      a residue with that number exists but has a different name
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, List

import numpy as np

from .core import NEG_ATOMS, POS_ATOMS, Complex
from .evidence import CONTACT_CUTOFF, SALT_CUTOFF, _acidic_atoms, _basic_atoms, contacts


@dataclass
class Check:
    kind: str
    target: str
    verify_type: str
    expectation: str
    status: str
    observed: str
    detail: str = ''
    claim: str = ''
    quote: str = ''
    locator: str = ''
    value: float = None

    def as_dict(self):
        return asdict(self)


def _res(cx: Complex, resid: int, chain='receptor'):
    d = cx.receptor if chain == 'receptor' else cx.ligand
    return d.get(int(resid))


def _check_name(actual, expected) -> bool:
    if not expected:
        return True
    return actual.upper().startswith(expected.upper()[:3])


# ---------------------------------------------------------------------------
def _verify_salt_bridge(cx_by_tag, resid, params):
    cutoff = params.get('cutoff', SALT_CUTOFF)
    need = params.get('min_partners', 1)
    per, best, found_name = {}, 0, None
    for tag, cx in cx_by_tag.items():
        r = _res(cx, resid)
        if r is None:
            continue
        found_name = found_name or r.name
        pa = _basic_atoms(r)
        if not pa:
            continue
        partners, dmin = [], 9.9
        for lid in cx.ligand:
            for x in _acidic_atoms(cx, lid):
                d = min(float(np.linalg.norm(x - y)) for y in pa)
                dmin = min(dmin, d)
                if d < cutoff:
                    partners.append(f'{cx.l(lid).one}{lid}')
        if partners:
            per[tag] = sorted(set(partners))
        best = max(best, len(set(partners)))
    obs = '; '.join(f'{t}:{",".join(p)}' for t, p in per.items()) or 'no salt bridge to any peptide acid'
    ok = best >= need
    return ('PASS' if ok else 'FAIL'), obs, found_name, None


def _verify_buried_contact(cx_by_tag, resid, params):
    need_n = params.get('min_complexes', 1)
    need_b = params.get('min_dsasa', 40.0)
    dsasa = params['_dsasa']
    per, found_name = {}, None
    n_contact, best_burial = 0, 0.0
    for tag, cx in cx_by_tag.items():
        r = _res(cx, resid)
        if r is None:
            continue
        found_name = found_name or r.name
        ct = params['_contacts'][tag]
        touches = sorted({lid for (rid, lid) in ct if rid == resid})
        b = float(dsasa[tag][0].get(resid, 0.0))
        if touches:
            n_contact += 1
            best_burial = max(best_burial, b)
            per[tag] = f'{",".join(f"{cx.l(l).one}{l}" for l in touches)} (dSASA {b:.0f})'
    obs = '; '.join(f'{t}:{v}' for t, v in per.items()) or 'no peptide contact'
    ok = (n_contact >= need_n) and (best_burial >= need_b)
    return ('PASS' if ok else 'FAIL'), obs, found_name, best_burial


def _verify_peptide_salt_bridge(cx_by_tag, tag, resid, params):
    cx = cx_by_tag.get(tag)
    if cx is None:
        return 'NO_RESIDUE', f'complex {tag} not loaded', None, None
    r = _res(cx, resid, 'ligand')
    if r is None:
        return 'NO_RESIDUE', f'{tag} has no modelled residue {resid}', None, None
    cutoff = params.get('cutoff', SALT_CUTOFF)
    need = params.get('min_partners', 1)
    acids = _acidic_atoms(cx, resid)
    if not acids:
        return 'FAIL', f'{tag} {r.one}{resid} has no carboxylate/OXT atoms', r.name, None
    rows, partners = [], 0
    for m in params.get('murj', []):
        mr = _res(cx, m)
        if mr is None:
            rows.append(f'{m}:absent')
            continue
        pa = _basic_atoms(mr)
        if not pa:
            rows.append(f'{m}:not basic')
            continue
        d = min(float(np.linalg.norm(x - y)) for x in acids for y in pa)
        rows.append(f'{mr.one}{m}:{d:.2f}')
        if d < cutoff:
            partners += 1
    # also report the closest basic MurJ residue overall (explains a FAIL)
    closest, cmin = None, 9.9
    for rid, rr in cx.receptor.items():
        pa = _basic_atoms(rr)
        if not pa:
            continue
        d = min(float(np.linalg.norm(x - y)) for x in acids for y in pa)
        if d < cmin:
            cmin, closest = d, f'{rr.one}{rid}'
    obs = ', '.join(rows) + f'; closest basic MurJ residue = {closest} at {cmin:.2f} A'
    return ('PASS' if partners >= need else 'FAIL'), obs, r.name, round(cmin, 2)


def _verify_peptide_contact(cx_by_tag, tag, resid, params):
    cx = cx_by_tag.get(tag)
    if cx is None:
        return 'NO_RESIDUE', f'complex {tag} not loaded', None, None
    r = _res(cx, resid, 'ligand')
    if r is None:
        return 'NO_RESIDUE', f'{tag} has no modelled residue {resid}', None, None
    ct = params['_contacts'][tag]
    n = sum(rec['n'] for (rid, lid), rec in ct.items() if lid == resid)
    dmin = min([rec['dmin'] for (rid, lid), rec in ct.items() if lid == resid] or [9.9])
    need = params.get('min_atoms', 1)
    obs = f'{n} atom contacts, dmin {dmin:.2f} A' if n else 'no contact'
    return ('PASS' if n >= need else 'FAIL'), obs, r.name, (round(dmin, 2) if n else None)


# ---------------------------------------------------------------------------
_MURJ_RULES = {'salt_bridge': _verify_salt_bridge,
              'buried_contact': _verify_buried_contact}
_PEPTIDE_RULES = {'peptide_salt_bridge': _verify_peptide_salt_bridge,
                  'peptide_contact': _verify_peptide_contact}


def run_checks(prior_data: dict, complexes: Dict[str, Complex], metrics: dict) -> List[Check]:
    """metrics = {'dsasa': {tag: (rec, lig, ..)}, 'contacts': {tag: ct}}"""
    base = dict(_dsasa=metrics['dsasa'], _contacts=metrics['contacts'])
    out: List[Check] = []

    for entry in prior_data.get('murj_residues', []):
        resid = int(entry['resid'])
        v = entry['verify']
        params = dict(v.get('params', {}), **base)
        fn = _MURJ_RULES.get(v['type'])
        if fn is None:
            raise ValueError(f'murj residue {resid}: verify type "{v["type"]}" is not a '
                             f'MurJ-side rule (allowed: {sorted(_MURJ_RULES)})')
        status, obs, name, val = fn(complexes, resid, params)
        if name and not _check_name(name, entry.get('resname')):
            status, obs = 'MISMATCH', f'structure has {name}{resid}, prior says {entry.get("resname")}'
        out.append(Check('murj_residue', f'{entry.get("resname", "?")[:3]}{resid}',
                         v['type'], 'pass', status, obs,
                         detail='; '.join(entry.get('roles', [])),
                         claim=entry.get('claim', ''), quote=entry.get('quote', ''),
                         locator=entry.get('locator', ''), value=val))

    for entry in prior_data.get('resistance_mutations', []):
        resid = int(entry['resid'])
        v = entry['verify']
        params = dict(v.get('params', {}), **base)
        status, obs, name, val = _verify_buried_contact(complexes, resid, params)
        if name and not _check_name(name, entry.get('wt')):
            status, obs = 'MISMATCH', f'structure has {name}{resid}, prior says {entry.get("wt")}'
        out.append(Check('resistance', f'{entry["wt"][:3].title()}{resid}{entry["mut"][:3].title()}',
                         v['type'], 'pass', status, obs,
                         detail=entry.get('phenotype', ''),
                         claim=entry.get('claim', ''), quote=entry.get('quote', ''),
                         locator=entry.get('locator', ''), value=val))

    for entry in prior_data.get('peptide_residues', []):
        tag, resid = entry['peptide'], int(entry['resid'])
        v = entry['verify']
        params = dict(v.get('params', {}), **base)
        exp = entry.get('expectation', 'pass')
        fn = _PEPTIDE_RULES.get(v['type'])
        if fn is None:
            raise ValueError(f'{tag} residue {resid}: verify type "{v["type"]}" is not a '
                             f'peptide-side rule (allowed: {sorted(_PEPTIDE_RULES)})')
        status, obs, name, val = fn(complexes, tag, resid, params)
        if name and not _check_name(name, entry.get('resname')):
            status, obs = 'MISMATCH', f'structure has {name}{resid}, prior says {entry.get("resname")}'
        if status == 'FAIL' and exp == 'flag':
            status = 'FLAG'
        elif status == 'PASS' and exp == 'weak':
            status = 'WEAK'
        detail = entry.get('role', '')
        if entry.get('note'):
            detail += ' | ' + entry['note']
        out.append(Check('peptide_residue', f'{tag}:{entry.get("resname", "?")[:3].title()}{resid}',
                         v['type'], exp, status, obs, detail=detail,
                         claim=entry.get('claim', '') or entry.get('role', ''),
                         quote=entry.get('quote', ''), locator=entry.get('locator', ''), value=val))
    return out


def summarize(checks: List[Check]) -> dict:
    counts: Dict[str, int] = {}
    for c in checks:
        counts[c.status] = counts.get(c.status, 0) + 1
    n = len(checks)
    ok = counts.get('PASS', 0)
    return dict(total=n, passed=ok, failed=counts.get('FAIL', 0),
                flagged=counts.get('FLAG', 0), weak=counts.get('WEAK', 0),
                no_residue=counts.get('NO_RESIDUE', 0), mismatch=counts.get('MISMATCH', 0),
                pass_rate=round(ok / max(1, n), 2), by_status=counts)
