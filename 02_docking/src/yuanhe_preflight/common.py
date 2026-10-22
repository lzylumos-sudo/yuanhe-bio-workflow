from __future__ import annotations
import hashlib, json, os, uuid
from pathlib import Path
from rdkit import Chem, rdBase

SCHEMA = 'yuanhe-predocking/2'

def hash_obj(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(',', ':')).encode()).hexdigest()

def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''): h.update(block)
    return h.hexdigest()

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def write(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    os.replace(temp, path)

def canonical(value):
    """Conservative standardization: sanitize/canonicalize only. Never desalinate,
    neutralize, tautomerize, discard isotopes, or invent unspecified stereochemistry."""
    mol = Chem.MolFromSmiles(value) if isinstance(value, str) else Chem.Mol(value)
    if mol is None or not mol.GetNumAtoms(): raise ValueError('invalid_structure')
    Chem.SanitizeMol(mol)
    for atom in mol.GetAtoms(): atom.SetAtomMapNum(0)
    mol = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

def identity(smiles):
    s = canonical(smiles); mol = Chem.MolFromSmiles(s)
    # Same recipes as zvss.ids; existing pose_id recipe is untouched.
    cid = 'CMP' + hashlib.sha256(s.encode()).hexdigest()[:16].upper()
    parts = ('state', s, Chem.GetFormalCharge(mol), '')
    sid = 'STA' + hashlib.sha256(json.dumps(parts, ensure_ascii=False, sort_keys=True, default=str).encode()).hexdigest()[:16].upper()
    return {'smiles': s, 'compound_id': cid, 'state_id': sid,
            'structure_sha256': hash_obj({'exact_smiles': s, 'normalization': 'sanitize-canonical-isomeric-v1'}),
            'normalization': 'sanitize-canonical-isomeric-v1', 'rdkit_version': rdBase.rdkitVersion,
            'fragments': sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in Chem.GetMolFrags(mol, asMols=True)),
            'salt_parent_mapping': 'none; exact full structure only'}

def sdf_smiles(path):
    mols = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if len(mols) != 1 or mols[0] is None: raise ValueError('one_valid_SDF_record_required')
    return canonical(mols[0])

def seal(record):
    record = dict(record); record.pop('record_sha256', None)
    return {**record, 'record_sha256': hash_obj(record)}

def verify_seal(record):
    if not isinstance(record, dict) or seal(record)['record_sha256'] != record.get('record_sha256'):
        raise ValueError('screening_record_hash_mismatch')
