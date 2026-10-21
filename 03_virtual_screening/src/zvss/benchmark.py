"""Offline fixed-pool comparisons. No docking, training or source-run mutations.

python -m zvss.benchmark --scoring scoring.json --selection selection.json
    --output NEW_DIRECTORY --top-k 10 [--labels labels.csv]
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
from .util import safe_float


def compare(scoring, selection, *, top_k=10, labels=None, seeds=(17, 41, 97, 137, 211)):
    if top_k < 1:
        raise ValueError('top_k must be positive')
    rows = scoring['ranking']
    ids = [e['compound_id'] for e in rows]
    if len(set(ids)) != len(ids):
        raise ValueError('Expected one representative per compound')
    if labels is not None and (set(labels) != set(ids) or any(type(v) is not int or v not in (0, 1) for v in labels.values())):
        raise ValueError('Labels must cover exactly the fixed pool, once, with binary 0/1 values')
    if scoring.get('policy_version') != selection.get('policy_version'):
        raise ValueError('Scoring and selection policy versions differ')
    pose_ids = {e['compound_id']: e['pose_id'] for e in rows}
    for role in ('primary', 'preliminary', 'exploratory'):
        for e in selection.get(role, []):
            if pose_ids.get(e['compound_id']) != e.get('pose_id'):
                raise ValueError('Selection pose does not match the scored representative')
    verdicts = selection.get('qualification', {})
    affinity_rows = [e for e in rows if safe_float(e.get('raw', {}).get('affinity')) is not None]
    affinity = [e['compound_id'] for e in sorted(affinity_rows,
                key=lambda e: (float(e['raw']['affinity']), e['compound_id']))]
    evidence = [e['compound_id'] for e in sorted(rows,
                key=lambda e: (e['tier'], -e['evidence_score'], e['compound_id']))]
    strategies = {'docking_affinity': affinity[:top_k], 'evidence_ranking': evidence[:top_k],
                  'formal_gate': [cid for cid in evidence if verdicts.get(cid, {}).get('formal_qualified')][:top_k],
                  'actual_primary': [e['compound_id'] for e in selection.get('primary', [])][:top_k]}
    for seed in seeds:
        pool = sorted(ids); random.Random(seed).shuffle(pool)
        strategies[f'random_seed_{seed}'] = pool[:top_k]
    reports = {}
    exclusions = {e['compound_id'] for e in selection.get('exclusions', [])}
    prevalence = sum(labels.values())/len(ids) if labels is not None and ids else None
    for name, chosen in strategies.items():
        if len(chosen) != len(set(chosen)) or not set(chosen).issubset(ids):
            raise ValueError('Selection must contain unique members of the same fixed pool')
        reports[name] = {'ids': chosen, 'selected': len(chosen), 'requested': top_k,
                         'unfilled': max(0, top_k-len(chosen)),
                         'excluded_by_current_policy': sum(cid in exclusions for cid in chosen),
                         'formal_qualified': sum(bool(verdicts.get(cid, {}).get('formal_qualified')) for cid in chosen),
                         'unknown_qualification': sum(cid not in verdicts and cid not in exclusions for cid in chosen)}
        if labels is not None:
            hits = sum(labels[cid] for cid in chosen)
            precision = hits/len(chosen) if chosen else None
            reports[name].update(hits=hits, precision_at_selected=precision,
                enrichment_at_selected=precision/prevalence if precision is not None and prevalence else None,
                recall=hits/sum(labels.values()) if sum(labels.values()) else None)
    return {'schema': 'screening-comparison/1', 'pool_size': len(ids), 'top_k': top_k,
            'selection_version': selection.get('selection_version'),
            'labels_supplied': labels is not None, 'strategies': reports,
            'weight_diagnostics': scoring.get('diagnostics'),
            'limitations': ['Same archived pool, not a new prospective experiment.',
                           'Qualification agreement is policy consistency, not biological accuracy.',
                           'Without independent labels, no activity or enrichment claim is computed.',
                           'Different selected counts are reported; do not compare precision without coverage.',
                           'Runtime, synthesis, wet-lab outcomes and target generalization are not measured here.']}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scoring', type=Path, required=True)
    p.add_argument('--selection', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--top-k', type=int, default=10)
    p.add_argument('--labels', type=Path)
    a = p.parse_args(argv)
    paths = [a.scoring, a.selection] + ([a.labels] if a.labels else [])
    inputs = {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    labels = None
    if a.labels:
        with a.labels.open(newline='', encoding='utf-8-sig') as f:
            rows = list(csv.DictReader(f))
        labels = {row['compound_id']: int(row['active']) for row in rows}
        if len(labels) != len(rows):
            raise ValueError('Duplicate label compound_id')
    result = compare(json.loads(a.scoring.read_text()), json.loads(a.selection.read_text()),
                     top_k=a.top_k, labels=labels)
    result['input_sha256'] = inputs
    result['code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    a.output.mkdir(parents=True, exist_ok=False)  # Never overwrite an old comparison/run.
    (a.output/'comparison.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    lines = ['# Fixed-pool screening comparison', '', '| Strategy | Selected | Formal | Excluded | Unfilled |',
             '|---|---:|---:|---:|---:|']
    for name, row in result['strategies'].items():
        lines.append(f"| {name} | {row['selected']} | {row['formal_qualified']} | {row['excluded_by_current_policy']} | {row['unfilled']} |")
    lines += ['', *['- '+note for note in result['limitations']], '',
              'Machine-readable metrics, identities and input hashes: comparison.json']
    (a.output/'comparison.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'output': str(a.output), 'pool_size': result['pool_size'], 'labels_supplied': labels is not None}))


if __name__ == '__main__':
    main()
