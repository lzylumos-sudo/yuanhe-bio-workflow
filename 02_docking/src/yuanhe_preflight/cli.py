"""Local candidate-library preflight: no docking, generation, network or MD."""
import argparse, csv, json
from pathlib import Path
from rdkit import Chem
from .common import canonical, write, read
from .pipeline import risk_stage, availability_stage
from .config import validate

def import_candidates(path,out):
    path=Path(path);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if path.suffix.lower() in {'.sdf','.sd'}:
        entries=[]
        for n,mol in enumerate(Chem.SDMolSupplier(str(path),removeHs=False)):
            entries.append(({'id':'import_'+str(n),'smiles':canonical(mol) if mol else '', 'source_library':str(path)},mol))
    else:
        with path.open(newline='',encoding='utf-8-sig') as f:
            reader=csv.DictReader(f,delimiter=',' if path.suffix.lower()=='.csv' else '\t')
            smiles_columns=[name for name in (reader.fieldnames or [])
                            if name is not None and name.strip().casefold()=='smiles']
            if len(smiles_columns)!=1:
                raise ValueError('Molecule tables require exactly one smiles column (case insensitive)')
            values=list(reader)
        entries=[]
        for n,r in enumerate(values):
            smiles=(r.get(smiles_columns[0]) or '').strip()
            row=dict(r,id=r.get('id') or 'import_'+str(n),smiles=smiles,source_library=str(path))
            entries.append((row,Chem.MolFromSmiles(smiles) if smiles else None))
    rows=[]
    for n,(row,mol) in enumerate(entries):
        if mol is not None and mol.GetNumAtoms():
            # Retain a structure record; only selected docking inputs receive 3D conformers.
            mol=Chem.AddHs(mol)
            sdf=out/f'{n:08d}.sdf'
            with Chem.SDWriter(str(sdf)) as w:w.write(mol)
            row['sdf']=str(sdf.resolve())
        rows.append(row)
    return rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',required=True);ap.add_argument('--config',required=True);ap.add_argument('--out',required=True)
    a=ap.parse_args();cfg=validate(read(a.config));out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'availability/qualified.json').exists():raise ValueError('Completed output exists; choose a new output. Shared search cache is reused.')
    rows=import_candidates(a.input,out/'imported')
    risk_stage(rows,cfg,out/'risk',base=out,batch_id=out.name)
    result=availability_stage(read(out/'risk/records.json'),cfg,out/'availability')
    print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
