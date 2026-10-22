"""Generate conformers only after qualification, preserving the chemical state."""
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from .gate import require_ligand


def prepare_selected(rows, base, out, seed):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    prepared = []
    for row in rows:
        source = Path(base) / row['sdf']
        require_ligand(source, row['qualification_evidence'])
        mols = list(Chem.SDMolSupplier(str(source), removeHs=False))
        if len(mols) != 1 or mols[0] is None:
            raise ValueError('One qualified molecule required')
        mol = Chem.AddHs(mols[0], addCoords=True)
        if not mol.GetNumConformers() or not mol.GetConformer().Is3D():
            mol.RemoveAllConformers()
            params = AllChem.ETKDGv3(); params.randomSeed = int(seed)
            if AllChem.EmbedMolecule(mol, params) != 0:
                raise ValueError('Selected ligand conformer generation failed: ' + row['compound_id'])
        target = out / (row['compound_id'] + '.sdf')
        with Chem.SDWriter(str(target)) as writer:
            writer.write(mol)
        require_ligand(target, row['qualification_evidence'])
        prepared.append({**row, 'sdf': str(target.resolve())})
    return prepared
