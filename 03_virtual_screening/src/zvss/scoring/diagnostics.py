"""Fixed-pool sensitivity diagnostics, not calibrated statistical uncertainty."""
from ..util import safe_float


def sensitivity(entries, weights, *, relative_change=0.2):
    if not 0 < relative_change < 1:
        raise ValueError('sensitivity_relative_change must be between 0 and 1')
    weights = {k: float(w) for k, w in weights.items() if float(w) > 0}
    rows = [e for e in entries if e['tier'] != 'D']

    def order(w):
        total = sum(w.values()) or 1.
        def key(e):
            score = sum(weight * (safe_float(e['components'].get(k)) or 0.)
                        for k, weight in w.items()) / total
            return (e['tier'], -round(score, 6), e['compound_id'])
        return [e['compound_id'] for e in sorted(rows, key=key)]

    baseline = order(weights)
    scenarios = []
    for key in sorted(weights):
        for factor in (1-relative_change, 1+relative_change):
            w = dict(weights); w[key] *= factor
            scenarios.append({'name': f'{key}_x{factor:g}', 'weights': w, 'order': order(w)})
        if len(weights) > 1:
            w = {k: v for k, v in weights.items() if k != key}
            scenarios.append({'name': f'without_{key}', 'weights': w, 'order': order(w)})
    ranks = {cid: [n] for n, cid in enumerate(baseline, 1)}
    for scenario in scenarios:
        for n, cid in enumerate(scenario['order'], 1):
            ranks[cid].append(n)
    return {'scope': 'fixed candidate pool; weight sensitivity and component ablation only',
            'not_claimed': 'statistical uncertainty, activity prediction or cross-run score comparability',
            'relative_change': relative_change, 'baseline_order': baseline, 'scenarios': scenarios,
            'rank_ranges': {cid: {'baseline': ranks[cid][0], 'best': min(ranks[cid]),
                                  'worst': max(ranks[cid]), 'span': max(ranks[cid])-min(ranks[cid])}
                            for cid in baseline}}
