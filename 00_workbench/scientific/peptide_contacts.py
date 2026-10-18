"""Standalone peptide-to-small-molecule contact evidence; no UI imports."""
from pathlib import Path
from . import artifacts as a

def compare(receptor, reference_receptor, peptide, candidates, selection, out):
 """Distance-contact comparison with controls and a same-pool ranking ablation."""
 import numpy as np
 from rdkit import Chem
 receptor,reference_receptor,peptide,out=map(Path,(receptor,reference_receptor,peptide,out))
 def atoms(path):
  result=[]
  for line in path.read_text().splitlines():
   if line.startswith(('ATOM  ','HETATM')) and (line[76:78].strip() or line[12:16].strip()[0]) not in ('H','D'):
    result.append(((line[21:22].strip(),line[22:27].strip(),line[17:20].strip(),line[12:16].strip()),[float(line[30:38]),float(line[38:46]),float(line[46:54])]))
  return result
 r=atoms(receptor);original=dict(atoms(reference_receptor));common=[(key,xyz) for key,xyz in r if key in original]
 if not r or not original or len(common)<.95*len(original):raise ValueError('Receptor atom mapping incomplete')
 delta=np.array([xyz for k,xyz in common])-np.array([original[k] for k,xyz in common])
 rms=float(np.sqrt(np.mean(np.sum(delta**2,axis=1))))
 if rms>.05:raise ValueError('Peptide and docking receptor coordinate frames differ')
 xyz=np.array([x for k,x in r]);keys=[':'.join(k[:3]) for k,x in r];pep=np.array([x for k,x in atoms(peptide)])
 if not len(pep):raise ValueError('Peptide has no heavy atoms')
 def contacts(coords):
  # Chunking bounds memory for large peptide inputs.
  close=np.zeros(len(xyz),dtype=bool)
  for start in range(0,len(coords),64):close|=(np.sum((xyz[:,None,:]-coords[None,start:start+64,:])**2,axis=2).min(1)<=4.5**2)
  return {key for key,hit in zip(keys,close) if hit}
 reference=contacts(pep)
 if not reference:raise ValueError('Reference peptide has no receptor contacts at 4.5 A')
 def score(q):
  overlap=reference&q
  return {'recall':len(overlap)/len(reference),'precision':len(overlap)/len(q) if q else 0.,'jaccard':len(overlap)/len(reference|q),'matched_residues':sorted(overlap),'contact_residues':sorted(q)}
 translated=pep+np.array([max(xyz[:,0])-min(pep[:,0])+100,0,0]);controls={'self':score(contacts(pep)),'translated_100A_beyond_receptor':score(contacts(translated))}
 if controls['self']['jaccard']!=1 or controls['translated_100A_beyond_receptor']['jaccard']!=0:raise ValueError('Geometric contact controls failed')
 ranks={x['compound_id']:x for key in ('preliminary','primary','exploratory') for x in selection.get(key,[])};rows=[]
 for row in candidates:
  src=Path(row['sdf']);mol=Chem.SDMolSupplier(str(src),removeHs=True)[0]
  if mol is None:raise ValueError('Candidate parse failure')
  rec=ranks[row['compound_id']];rows.append({'compound_id':row['compound_id'],'pose_id':row['pose_id'],'source_sha256':a.digest(src),'affinity_kcal_mol':rec['raw']['affinity'],**score(contacts(mol.GetConformer().GetPositions()))})
 baseline=sorted(rows,key=lambda x:(x['affinity_kcal_mol'],x['compound_id']))
 # Fixed lexicographic rule makes this an auditable ablation, not a fitted predictor.
 guided=sorted(rows,key=lambda x:(-x['recall'],x['affinity_kcal_mol'],x['compound_id']))
 result={'scope':'Same shortlisted candidate pool; no generation changes or biological performance inference','cutoff_A':4.5,'contact_definition':'receptor and ligand heavy-atom distance; chain + author residue id + residue name','receptor_frame_rmsd_A':rms,'reference_contact_residues':sorted(reference),'controls':controls,'candidates':rows,'ablation':{'baseline':'affinity ascending','peptide_guided':'contact recall descending, affinity ascending tie break','baseline_order':[x['compound_id'] for x in baseline],'peptide_guided_order':[x['compound_id'] for x in guided],'top1_changed':bool(rows) and baseline[0]['compound_id']!=guided[0]['compound_id']},'source_hashes':{'receptor':a.digest(receptor),'peptide':a.digest(peptide)},'limitations':['Distance overlap does not prove pharmacophore equivalence or preserved function','No active/inactive labels: enrichment and hit-rate improvement cannot be measured','Current generation uses peptide-defined pocket; peptide sequence/function is not a learned constraint']}
 a.atomic(out/'peptide_evidence.json',result);return result


def calculate(d,c,p,out):
 rows=[{**r,'sdf':str(a.within(d,r['sdf']))} for r in a.read(d/'screening/candidates.json')]
 selection=a.read(d/'screening/workspace/runs/screen/stages/selection/summary.json')
 return compare(d/'docking/receptor.pdb',d/'prepare/receptor_full.pdb',d/'prepare/peptide_reference.pdb',rows,selection,out)


def main():
 import argparse,json
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--manifest',required=True,type=Path);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
 manifest=a.read(args.manifest);root=args.manifest.resolve().parent
 def resolve(value):
  path=Path(value);return path if path.is_absolute() else root/path
 rows=[{**r,'sdf':str(resolve(r['sdf']))} for r in manifest['candidates']]
 result=compare(resolve(manifest['receptor']),resolve(manifest['reference_receptor']),resolve(manifest['peptide']),rows,manifest['selection'],args.output)
 print(json.dumps({'status':'completed','candidates':len(result['candidates']),'output':str(args.output)}))

if __name__=='__main__':main()
