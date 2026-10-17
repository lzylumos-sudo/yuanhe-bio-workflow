"""Read-only input/core/box inspection; writes a NEW explicitly named report file."""
import argparse,json,hashlib,csv
from pathlib import Path

def atoms(path):
 out=[]
 for l in path.read_text().splitlines():
  if l.startswith(('ATOM  ','HETATM')):
   out.append({'chain':l[21:22],'resid':int(l[22:26]),'resname':l[17:20].strip(),'name':l[12:16].strip(),'element':l[76:78].strip(),'xyz':[float(l[30:38]),float(l[38:46]),float(l[46:54])]})
 return out

def main():
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['structure','box','core','counts']);p.add_argument('input',type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
 if a.out.exists():raise SystemExit('Output exists; choose a new name.')
 if a.mode=='core':
  from rdkit import Chem
  from rdkit.Chem import Draw
  m=Chem.SDMolSupplier(str(a.input),removeHs=True)[0]
  if m is None:raise SystemExit('Invalid first SDF record')
  for atom in m.GetAtoms():atom.SetProp('atomNote',str(atom.GetIdx()+1))
  image=Draw.MolsToGridImage([m],molsPerRow=1,subImgSize=(800,600),useSVG=True)
  a.out.write_text(image);return
 if a.mode in ('structure','box'):
  rows=atoms(a.input);result={'input':str(a.input),'sha256':hashlib.sha256(a.input.read_bytes()).hexdigest()}
  if a.mode=='structure':
   result['chains']={}
   for c in sorted({x['chain'] for x in rows}):
    rr=[x for x in rows if x['chain']==c];res={x['resid']:x['resname'] for x in rr}
    result['chains'][c]={'atoms':len(rr),'residues':len(res),'range':[min(res),max(res)],'missing_in_MurJ_3_500':[n for n in range(3,501) if n not in res] if c=='A' else None}
   result['construct_mutations']={str(n):next((x['resname'] for x in rows if x['chain']=='A' and x['resid']==n),None) for n in [5,12,13]}
   result['warning']='Observed coordinates only. Not proof of full-length WT protein or complete lysis protein.'
  else:
   wanted={50,51,53,54,237,240,241,244}
   region=[x for x in rows if x['chain']=='A' and x['resid'] in wanted and x['element'].upper()!='H']
   found={x['resid'] for x in region}
   if found!=wanted:raise SystemExit(f'Region residues missing: {wanted-found}')
   low=[min(x['xyz'][i] for x in region) for i in range(3)];high=[max(x['xyz'][i] for x in region) for i in range(3)]
   sizes=[max(14.,high[i]-low[i]+10.) for i in range(3)]
   result.update(candidate_box={'center':[round((low[i]+high[i])/2,3) for i in range(3)],'size':[round(x,3) for x in sizes]},residues=sorted(wanted),review_required=True,over_profile_max=any(s>30 for s in sizes))
   result['warning']='A candidate box from local TM2/TM7 profile, not a validated binding site. Review against peptide reference in the SAME coordinate frame.'
 else:
  d=a.input;result={'run':d.name,'stages':{}}
  for stage in ('prepare','peptide2mol','reinvent4','docking','screening','gromacs'):
   q=d/stage/'stage.json';r=json.loads(q.read_text()) if q.exists() else {}
   result['stages'][stage]={'status':r.get('status','pending'),'counts':r.get('counts',{})}
  results=d/'docking/results.json'
  if results.exists():
   rows=json.loads(results.read_text());result['docking_reconciliation']={'input_records':len(rows),'success':sum(x['status']=='ok' for x in rows),'failed':sum(x['status']!='ok' for x in rows)}
  candidates=d/'screening/molecules.csv'
  if candidates.exists():
   with candidates.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
   truth={'true','1','yes'}
   result['screening_reconciliation']={'rows':len(rows),'unique_displayed_compounds':len({x['compound_id'] for x in rows if x.get('compound_id')}),'unique_formally_qualified':len({x['compound_id'] for x in rows if x.get('compound_id') and (x.get('qualification','').strip().lower()=='formal' or x.get('formal_qualified','').strip().lower() in truth)}),'note':'Displayed compounds include preliminary/exploratory entries; not all are formally qualified.'}
 a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
