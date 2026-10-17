# -*- coding: utf-8 -*-
"""pepwedge.priors -- load the literature-derived priors from an external file.

Why this exists
---------------
Every number the design pipeline uses falls into one of two classes:

  * **measured** from the deposited coordinates (see ``evidence.py``)
  * **asserted** in the paper and used as a design prior (this module)

Keeping the second class in a versioned JSON file -- instead of hard-coded
constants buried in ``core.py`` -- means:

  * the provenance of each prior (quote + locator) is auditable;
  * a prior can be *verified against the structure* (``verify.py``) and the
    verification can fail loudly instead of being silently trusted;
  * the pipeline still runs without an LLM in the loop, and without Python
    edits, when a collaborator disagrees with a prior.

Set ``PEPWEDGE_PRIORS=/path/to/file.json`` to point at an alternative file.
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Dict, List

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PATH = _REPO_ROOT / 'priors' / 'murj_priors.json'

KNOWN_VERIFY = {'salt_bridge', 'buried_contact', 'peptide_salt_bridge', 'peptide_contact'}
# Which rule may be attached to which kind of claim. A peptide rule on a MurJ
# residue (or vice versa) is a modelling mistake, not a stylistic choice, so it
# is rejected instead of being silently dispatched to the wrong verifier.
SECTION_TYPES = {
    'murj_residues': {'salt_bridge', 'buried_contact'},
    'resistance_mutations': {'buried_contact'},
    'peptide_residues': {'peptide_salt_bridge', 'peptide_contact'},
}
_CACHE: Dict[str, tuple] = {}          # path -> ((mtime_ns, size), data)


class PriorError(RuntimeError):
    pass


def validate(data: dict) -> None:
    for key in ('schema_version', 'target', 'structures', 'murj_residues',
                'resistance_mutations', 'peptide_residues'):
        if key not in data:
            raise PriorError(f'priors file is missing required key: {key}')
    for s in data['structures']:
        for key in ('tag', 'file'):
            if key not in s:
                raise PriorError(f'structure entry missing "{key}": {s}')
    for r in data['murj_residues'] + data['resistance_mutations']:
        if 'resid' not in r:
            raise PriorError(f'residue entry missing "resid": {r}')
        v = r.get('verify')
        if not v or v.get('type') not in KNOWN_VERIFY:
            raise PriorError(f'residue {r.get("resid")} has no valid verify spec '
                             f'(known types: {sorted(KNOWN_VERIFY)})')
    for r in data['peptide_residues']:
        for key in ('peptide', 'resid', 'verify'):
            if key not in r:
                raise PriorError(f'peptide residue entry missing "{key}": {r}')
        if r['verify'].get('type') not in KNOWN_VERIFY:
            raise PriorError(f'peptide residue {r.get("resid")} has an unknown verify type')

    if not data['structures']:
        raise PriorError('priors file lists no structures')

    # Every prior must carry its provenance and a rule that is meaningful for
    # its own section -- this is what makes the claim auditable.
    for sec in SECTION_TYPES:
        for r in data.get(sec, []):
            if not (str(r.get('claim') or r.get('role') or '').strip()):
                raise PriorError(f'{sec} entry {r.get("resid")} has no claim/role text')
            for key in ('quote', 'locator'):
                if not str(r.get(key) or '').strip():
                    raise PriorError(f'{sec} entry {r.get("resid")} is missing "{key}" '
                                     f'-- every prior must carry its provenance')
            t = r['verify']['type']
            if t not in SECTION_TYPES[sec]:
                raise PriorError(f'{sec} entry {r.get("resid")} uses verify type "{t}"; '
                                 f'allowed here: {sorted(SECTION_TYPES[sec])}')


def load_priors(path=None, strict: bool = False) -> dict:
    """Load and validate the priors file. Missing file -> empty priors (warn)."""
    p = Path(path or os.environ.get('PEPWEDGE_PRIORS') or DEFAULT_PATH)
    key = str(p.resolve())
    sig = None
    if p.exists():
        st = p.stat()
        sig = (st.st_mtime_ns, st.st_size)
        hit = _CACHE.get(key)
        if hit and hit[0] == sig:
            return hit[1]
    if not p.exists():
        msg = f'priors file not found: {p} -- running without literature priors'
        if strict:
            raise PriorError(msg)
        warnings.warn(msg)
        data = dict(schema_version='missing', target={}, structures=[],
                    murj_residues=[], resistance_mutations=[], peptide_residues=[],
                    unverifiable_claims=[], design_thresholds={})
    else:
        with open(p, encoding='utf-8') as fh:
            data = json.load(fh)
        validate(data)
    data['_path'] = str(p)
    _CACHE[key] = (sig, data)      # re-reads automatically when the file changes
    return data


# --------------------------------------------------------------------------
# derived views used by core.py / pharmacophore.py
# --------------------------------------------------------------------------
def structure_map(path=None) -> Dict[str, str]:
    return {s['tag']: s['file'] for s in load_priors(path)['structures']}


def essential_map(path=None) -> Dict[int, str]:
    out = {}
    for r in load_priors(path)['murj_residues']:
        if 'essential' in r.get('roles', []):
            out[int(r['resid'])] = r.get('claim', 'essential')
    return out


def resistance_map(path=None) -> Dict[int, str]:
    return {int(m['resid']): f"{m['wt'][:3].title()}{m['resid']}{m['mut'][:3].title()} "
                             f"-> {m.get('phenotype', 'resistance')}"
            for m in load_priors(path)['resistance_mutations']}


def peptide_functional_map(path=None) -> Dict[str, Dict[int, str]]:
    out: Dict[str, Dict[int, str]] = {}
    for r in load_priors(path)['peptide_residues']:
        tag = r['peptide']
        label = f"{r.get('resname', 'X')[:3].title()}{r['resid']} ({r.get('role', '')})"
        out.setdefault(tag, {})[int(r['resid'])] = label
    return out


def residue_range(path=None):
    """(lo, hi) of the target protein's numbering; construct-only residues fall
    outside this range and must not be treated as MurJ.

    A missing range is warned about rather than silently defaulted: the default
    is only a coarse safety net (fusion constructs are usually numbered outside
    1..500) and would mis-classify a linker that happens to sit in that window.
    """
    r = load_priors(path)['target'].get('residue_range')
    if not r:
        warnings.warn('priors file has no target.residue_range -- falling back to '
                      '(1, 500); fusion-construct residues may be counted as target')
        return (1, 500)
    return (int(r[0]), int(r[1]))


def thresholds(path=None) -> dict:
    return load_priors(path).get('design_thresholds', {})


def summary(path=None) -> str:
    d = load_priors(path)
    return (f"priors {d.get('schema_version')} from {d.get('_path')}: "
            f"{len(d['structures'])} structures, "
            f"{len(d['murj_residues'])} essential/anchor residues, "
            f"{len(d['resistance_mutations'])} resistance alleles, "
            f"{len(d['peptide_residues'])} peptide residues, "
            f"{len(d.get('unverifiable_claims', []))} unverifiable claims")
