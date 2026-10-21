"""Exact full-structure supplier matching, with expiring local supply documents."""
import csv
from datetime import date
from pathlib import Path
from rdkit import Chem
from .common import canonical, file_hash

class Catalog:
    def __init__(self, descriptors=(), max_age_days=30):
        if type(max_age_days) is not int or not 1 <= max_age_days <= 365: raise ValueError('evidence_max_age_days must be 1..365')
        self.max_age_days = max_age_days; self.records = []; self.by_smiles = {}; self.versions = []; self.errors = []
        for d in descriptors:
            for k in ('path','version','source','license','trusted'):
                if not d.get(k): raise ValueError('catalog_descriptor_missing: '+k)
            if d['trusted'] is not True: raise ValueError('catalog_not_trusted')
            path = Path(d['path']).resolve()
            if not path.is_file(): raise FileNotFoundError('catalog_missing: '+str(path))
            catalog_sha=file_hash(path)
            self.versions.append({**d, 'path': str(path), 'sha256': catalog_sha})
            if path.suffix.lower() in {'.sdf','.sd'}:
                rows = []
                for n,mol in enumerate(Chem.SDMolSupplier(str(path), removeHs=False)):
                    if mol is None: self.errors.append({'row':n+1,'reason':'invalid_SDF'}); continue
                    rows.append({**mol.GetPropsAsDict(), 'smiles': canonical(mol)})
            elif path.suffix.lower() in {'.csv','.tsv','.smi','.smiles'}:
                with path.open(newline='',encoding='utf-8-sig') as f:
                    # .smi/.smiles uses a headered TSV: evidence cannot be inferred from two columns.
                    rows = list(csv.DictReader(f, delimiter=',' if path.suffix.lower()=='.csv' else '\t'))
            else: raise ValueError('unsupported_catalog_format')
            for n,row in enumerate(rows):
                try: exact = canonical(row['smiles'])
                except Exception: self.errors.append({'row':n+1,'reason':'invalid_SMILES'}); continue
                record = {**row, 'exact_smiles':exact, 'raw_smiles':row['smiles'], 'catalog_path':str(path),
                          'catalog_sha256':catalog_sha, 'catalog_version':d['version'], 'source':d['source'], 'license':d['license'],
                          'catalog_date':row.get('catalog_date') or d.get('date',''),
                          'salt_parent_mapping':{'full_structure':exact,'parent_smiles':row.get('parent_smiles',''),
                                                 'policy':'informational only; no parent-equivalence eligibility'}}
                record['evidence_path'] = str((path.parent/str(row.get('evidence_file',''))).resolve()) if row.get('evidence_file') else ''
                self.records.append(record);self.by_smiles.setdefault(exact,[]).append(record)

    def evaluate(self, record, today=None):
        today = today or date.today()
        try:
            checked = date.fromisoformat(str(record.get('checked_at','')))
            age = (today-checked).days
            if age < 0 or age > self.max_age_days: return 'supply_evidence_expired'
            # Old catalogs cannot become current merely by setting another query date.
            catalog_date = date.fromisoformat(str(record.get('catalog_date','')))
            if not 0 <= (today-catalog_date).days <= self.max_age_days: return 'catalog_expired'
            if record.get('valid_until') and date.fromisoformat(record['valid_until']) < today: return 'supply_evidence_expired'
        except (ValueError,TypeError): return 'catalog_hit_supply_unconfirmed'
        if any(not record.get(k) for k in ('vendor','sku','source','evidence_path','evidence_sha256')): return 'catalog_hit_supply_unconfirmed'
        p = Path(record['evidence_path'])
        if not p.is_file() or file_hash(p)!=record['evidence_sha256']: return 'supply_document_missing_or_changed'
        category = record.get('supply_type')
        if category not in {'in_stock','made_to_order'}: return 'catalog_hit_supply_unconfirmed'
        # A local reviewed vendor quote/stock snapshot; a product URL alone does not count.
        if record.get('evidence_kind') not in {'stock_confirmation','dated_quote','orderable_confirmation'}: return 'catalog_hit_supply_unconfirmed'
        return category

    def match(self, smiles, today=None):
        exact = canonical(smiles)
        return [{**r, 'status':self.evaluate(r,today)} for r in self.by_smiles.get(exact,[])]

    def valid(self, smiles, today=None):
        return [r for r in self.match(smiles,today) if r['status'] in {'in_stock','made_to_order'}]
