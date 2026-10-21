"""Versioned, individually validated alerts. A match is not a toxicity assay."""
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import FilterCatalog
from .common import read, file_hash, canonical, seal, hash_obj

DEFAULT_RULES = Path(__file__).with_name('configs') / 'risk_rules.v1.json'
CATEGORIES = {'toxicity_risk', 'nonspecific_reactivity', 'assay_interference', 'medicinal_chemistry'}
ACTIONS = {'exclude', 'review', 'warn'}

class RiskRules:
    def __init__(self, path=DEFAULT_RULES, policies=None):
        self.path = Path(path); self.data = read(path); self.sha256 = file_hash(path)
        self.policies = policies or {}; self.rules = {}; self.queries = {}; self.catalog = FilterCatalog.FilterCatalog()
        for rule in self.data['rules']:
            if not rule.get('enabled', True): continue
            rid = rule['rule_id']
            if rid in self.rules: raise ValueError('duplicate_rule_id: ' + rid)
            for field in ('SMARTS', 'risk_category', 'severity', 'default_action', 'source', 'version', 'positive_examples', 'negative_examples'):
                if not rule.get(field): raise ValueError('missing_rule_field: ' + field)
            if rule['risk_category'] not in CATEGORIES or rule['default_action'] not in ACTIONS or rule['severity'] not in {'high', 'medium', 'low'}:
                raise ValueError('invalid_rule_policy: ' + rid)
            query = Chem.MolFromSmarts(rule['SMARTS'])
            if query is None: raise ValueError('invalid_SMARTS: ' + rid)
            entry = FilterCatalog.FilterCatalogEntry(rid, FilterCatalog.SmartsMatcher(rid, rule['SMARTS'], 1))
            if not entry.IsValid(): raise ValueError('invalid_FilterCatalog_entry: ' + rid)
            self.catalog.AddEntry(entry); self.rules[rid] = rule; self.queries[rid] = query
            for expected, samples in ((True, rule['positive_examples']), (False, rule['negative_examples'])):
                for sample in samples:
                    mol = Chem.AddHs(Chem.MolFromSmiles(canonical(sample)))
                    if bool(entry.HasFilterMatch(mol)) != expected or bool(mol.HasSubstructMatch(query)) != expected:
                        raise ValueError('rule_example_failed: ' + rid + ':' + sample)
        if not self.rules: raise ValueError('no_enabled_risk_rules')
        if set(self.policies) - set(self.rules) or any(x not in ACTIONS for x in self.policies.values()):
            raise ValueError('unknown_rule_or_action')
        # PAINS cannot be relabelled as toxicity by a UI policy setting.
        self.version = self.data['version']

    def screen(self, smiles, exceptions=()):
        smiles = canonical(smiles); mol = Chem.MolFromSmiles(smiles); n = mol.GetNumAtoms(); explicit = Chem.AddHs(mol)
        alerts = []
        for entry in self.catalog.GetMatches(explicit):
            rid = entry.GetDescription(); rule = self.rules[rid]; action = self.policies.get(rid, rule['default_action'])
            matches = [sorted({i for i in match if i < n}) for match in explicit.GetSubstructMatches(self.queries[rid])]
            alerts.append({**rule, 'action': action, 'matched_atom_indices': matches})
        valid_exceptions = []
        for ex in exceptions:
            if ex.get('rule_id') not in {a['rule_id'] for a in alerts} or ex.get('structure_smiles') != smiles or ex.get('rules_sha256') != self.sha256:
                raise ValueError('exception_identity_or_rule_mismatch')
            if any(not str(ex.get(k, '')).strip() for k in ('reviewer', 'reason', 'reviewed_at')): raise ValueError('exception_requires_review_record')
            valid_exceptions.append(ex)
        handled = {e['rule_id'] for e in valid_exceptions}
        blocking = [a for a in alerts if a['rule_id'] not in handled and a['action'] in {'exclude', 'review'}]
        decision = 'exclude' if any(a['action'] == 'exclude' for a in blocking) else 'review' if blocking else 'pass'
        return seal({'decision': decision, 'alert_ids': [a['rule_id'] for a in alerts],
                     'alert_types': sorted({a['risk_category'] for a in alerts}), 'alerts': alerts,
                     'matched_atom_indices': {a['rule_id']: a['matched_atom_indices'] for a in alerts},
                     'atom_index_basis': 'zero-based atoms of canonical isomeric SMILES; explicit H omitted',
                     'severity': next((v for v in ('high','medium','low') if any(a['severity']==v for a in alerts)), 'none'),
                     'reason': 'No match to current rules' if not alerts else '; '.join(a['rule_id'] + ':' + a['action'] for a in alerts),
                     'smiles': smiles, 'rule_version': self.version, 'rules_sha256': self.sha256,
                     'policy_sha256': hash_obj(self.policies), 'rdkit_version': rdBase.rdkitVersion, 'exceptions': valid_exceptions})

    def svg(self, smiles):
        from rdkit.Chem.Draw import rdMolDraw2D
        s = canonical(smiles); mol = Chem.MolFromSmiles(s); result = self.screen(s)
        atoms = sorted({i for matches in result['matched_atom_indices'].values() for match in matches for i in match})
        draw = rdMolDraw2D.MolDraw2DSVG(560, 320); draw.drawOptions().addAtomIndices = True
        rdMolDraw2D.PrepareAndDrawMolecule(draw, mol, highlightAtoms=atoms); draw.FinishDrawing()
        return draw.GetDrawingText()
