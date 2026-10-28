"""Full RDKit alert catalogs and explicit, stereo-preserving parent mapping."""
from functools import lru_cache
from rdkit import Chem, rdBase
from rdkit.Chem import FilterCatalog
from rdkit.Chem.MolStandardize import rdMolStandardize as ms
from .common import canonical, seal, hash_obj


def standardize(smiles):
    exact = canonical(smiles)
    mol = ms.Cleanup(Chem.MolFromSmiles(exact))
    parent = ms.FragmentParent(mol)
    parent = ms.Uncharger().uncharge(parent)
    # Do not enumerate tautomers: this can alter stereochemistry and evidence identity.
    return {'exact_smiles': exact, 'canonical_smiles': canonical(parent),
            'fragments': [canonical(x) for x in Chem.GetMolFrags(mol, asMols=True)],
            'salt_changed': canonical(parent) != exact,
            'normalization_version': 'cleanup-fragment-parent-uncharge-stereo-v2'}


@lru_cache(maxsize=2)
def catalog(name):
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, name))
    return FilterCatalog.FilterCatalog(params)


def screen(smiles, rules, exceptions, engine=None):
    catalog_reviews = [x for x in exceptions if str(x.get('rule_id', '')).startswith(('PAINS:', 'BRENK:'))]
    result = rules.screen(smiles, [x for x in exceptions if x not in catalog_reviews])
    if not engine:
        if catalog_reviews: raise ValueError('catalog_review_requires_engine_policy')
        return result
    mol = Chem.MolFromSmiles(canonical(smiles))
    alerts = list(result['alerts']) if engine['structure_enabled'] else []
    if engine['structure_enabled']:
        for name in ('PAINS', 'BRENK'):
            for entry in catalog(name).GetMatches(mol):
                rid = name + ':' + entry.GetDescription()
                action = engine['catalog_rule_actions'].get(rid, engine['catalog_actions'][name])
                pairs = [m.atomPairs for m in entry.GetFilterMatches(mol)]
                atoms = [sorted({int(p.target) for p in group}) for group in pairs]
                alerts.append({'rule_id': rid, 'action': action, 'severity': 'medium',
                               'risk_category': 'assay_interference' if name == 'PAINS' else 'medicinal_chemistry',
                               'source': 'RDKit FilterCatalog/' + name, 'version': rdBase.rdkitVersion,
                               'matched_atom_indices': atoms})
    policy_sha = hash_obj({k: engine[k] for k in ('structure_enabled', 'catalog_actions', 'catalog_rule_actions')})
    for review in catalog_reviews:
        alert = next((a for a in alerts if a['rule_id'] == review['rule_id']), None)
        if (not alert or alert['action'] != 'review' or review.get('structure_smiles') != canonical(smiles) or
            review.get('rules_sha256') != rules.sha256 or review.get('catalog_policy_sha256') != policy_sha):
            raise ValueError('catalog_review_identity_rule_or_policy_mismatch')
        if any(not isinstance(review.get(k), str) or not review[k].strip() for k in ('reviewer', 'reason', 'reviewed_at')):
            raise ValueError('catalog_review_requires_review_record')
    result['exceptions'].extend(catalog_reviews)
    handled = {x['rule_id'] for x in result['exceptions']}
    blocking = [a for a in alerts if a['rule_id'] not in handled and a['action'] in ('exclude', 'review')]
    decision = 'exclude' if any(a['action'] == 'exclude' for a in blocking) else 'review' if blocking else 'pass'
    result.update(severity=next((v for v in ('high','medium','low') if any(a['severity']==v for a in alerts)), 'none'),
                  alerts=alerts, alert_ids=[x['rule_id'] for x in alerts],
                  alert_types=sorted({x['risk_category'] for x in alerts}),
                  matched_atom_indices={x['rule_id']: x['matched_atom_indices'] for x in alerts},
                  decision=decision, filter_status={'pass': 'PASS', 'review': 'REVIEW', 'exclude': 'BLOCK'}[decision],
                  reason='; '.join(x['rule_id'] + ':' + x['action'] for x in blocking) or 'No blocking match to configured rules',
                  catalog_policy_sha256=policy_sha)
    return seal(result)
