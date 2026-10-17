#!/usr/bin/env python
import os,json,time,hashlib,requests
from pathlib import Path
from rdkit import Chem
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from sbvs.paths import ROOT,sha,write_json
def fetch(url,path):
 path=Path(path)
 if not path.exists():
  r=requests.get(url,timeout=90);r.raise_for_status();path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(r.content)
 return {'url':url,'path':str(path.relative_to(ROOT)),'sha256':sha(path),'retrieved_at':time.time()}
def main():
 sources=[]
 pdb=ROOT/'data/raw/1DRF.pdb';ccd=ROOT/'data/reference/FOL.cif'
 sources.append(fetch('https://files.rcsb.org/download/1DRF.pdb',pdb))
 sources.append(fetch('https://files.rcsb.org/ligands/download/FOL.cif',ccd))
 atoms={};instances={}
 for line in pdb.read_text().splitlines():
  if line.startswith('HETATM') and line[17:20]=='FOL':
   key=(line[21],line[22:26].strip(),line[26]);instances.setdefault(key,[]).append(line)
 if len(instances)!=1:raise ValueError('FOL instance ambiguous: '+repr(list(instances)))
 instance,lines=next(iter(instances.items()))
 for l in lines:
  name=l[12:16].strip()
  if name in atoms:raise ValueError('Duplicate FOL atom '+name)
  atoms[name]=[float(l[30:38]),float(l[38:46]),float(l[46:54])]
 d=MMCIF2Dict(str(ccd));names=d['_chem_comp_atom.atom_id'];els=d['_chem_comp_atom.type_symbol'];charges=d['_chem_comp_atom.charge']
 rw=Chem.RWMol();mapping={}
 for name,el,q in zip(names,els,charges):
  if el=='H':continue
  if name not in atoms:raise ValueError('Missing FOL heavy atom '+name)
  a=Chem.Atom(el);a.SetFormalCharge(int(q));a.SetProp('ccd_atom_name',name);mapping[name]=rw.AddAtom(a)
 bt={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
 for a,b,t in zip(d['_chem_comp_bond.atom_id_1'],d['_chem_comp_bond.atom_id_2'],d['_chem_comp_bond.value_order']):
  if a in mapping and b in mapping:rw.AddBond(mapping[a],mapping[b],bt[t])
 m=rw.GetMol();Chem.SanitizeMol(m);conf=Chem.Conformer(m.GetNumAtoms());conf.Set3D(True)
 for name,i in mapping.items():conf.SetAtomPosition(i,atoms[name])
 m.AddConformer(conf);Chem.AssignStereochemistryFrom3D(m);m.SetProp('_Name','1DRF_FOL_cocrystal')
 stereo=d.get('_chem_comp_atom.pdbx_stereo_config',['N']*len(names))
 for name,s in zip(names,stereo):
  if name in mapping and s in ('R','S'):
   a=m.GetAtomWithIdx(mapping[name])
   if not a.HasProp('_CIPCode') or a.GetProp('_CIPCode')!=s:raise ValueError('CCD stereochemistry differs from crystal at '+name)
 w=Chem.SDWriter(str(ROOT/'data/reference/FOL_crystal.sdf'));w.write(m);w.close()
 write_json(ROOT/'data/reference/FOL_mapping.json',{'instance':instance,'ccd_atom_to_index':mapping,'max_coordinate_delta':0,'chemistry_source':'RCSB CCD','coordinate_source':'1DRF original HETATM'})
 smiles=['CCO','CCCO','CCCCO','CC(C)O','CC(=O)O','CC(=O)N','CCN','CCCN','CCOC','CCOCC','c1ccccc1','Cc1ccccc1','Oc1ccccc1','Nc1ccccc1','c1ccncc1','c1ccoc1','C1CCCCC1','C1CCNCC1','CC(=O)OC','O=C(O)c1ccccc1']
 demo=list(smiles)
 for n in range(2,20):
  for suffix in ('O','N','(=O)O','(=O)N','F','Cl'):
   s='C'*n+suffix
   if s not in demo and Chem.MolFromSmiles(s):demo.append(s)
   if len(demo)>=100:break
  if len(demo)>=100:break
 for n in (20,100):(ROOT/f'data/demo/ligands_{n}.smi').write_text('# 流程演示库；无活性标签\n'+''.join(f'{s} demo_{i:03d}\n' for i,s in enumerate(demo[:n],1)),encoding='utf-8')
 write_json(ROOT/'data/sources.json',sources)
 print(json.dumps({'FOL_instance':instance,'heavy_atoms':m.GetNumAtoms(),'demo_count':len(demo)}))
if __name__=='__main__':main()
