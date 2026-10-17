from pathlib import Path
import subprocess,json
from Bio.PDB import PDBParser,MMCIFParser,PDBIO,Select
from ..paths import write_json,sha
def plan(source,chains=('A',),model=0):
 source=Path(source)
 if source.suffix.lower()=='.pdbqt':return {'input_sha256':sha(source),'external':True,'protonation':'external/unknown','issues':['Prepared PDB coordinates required separately'],'confirmed':False}
 parser=MMCIFParser(QUIET=True,auth_chains=True) if source.suffix.lower() in ('.cif','.mmcif') else PDBParser(QUIET=True)
 structure=parser.get_structure('receptor',str(source));models=list(structure)
 if model>=len(models):raise ValueError('Model missing')
 mdl=models[model];available=[c.id for c in mdl]
 if not set(chains)<=set(available):raise ValueError('Requested chain missing')
 hetero=[];altloc=[];breaks=[];residues=[]
 for c in mdl:
  if c.id not in chains:continue
  prev=None
  for r in c:
   rid=f'{c.id}:{r.id[1]}{r.id[2].strip()}'
   if r.id[0]!=' ':hetero.append({'residue':rid,'name':r.resname});continue
   residues.append({'id':rid,'name':r.resname,'atoms':[a.name for a in r]})
   if prev is not None and 'C' in prev and 'N' in r and (prev['C']-r['N'])>1.9:breaks.append(rid)
   prev=r
   for a in r:
    if a.is_disordered():altloc.append(rid+':'+a.name)
 standard={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
 return {'input_sha256':sha(source),'model_index':model,'chains':list(chains),'available_chains':available,'hetero':hetero,'altloc':altloc,'chain_breaks':breaks,'residues':residues,'nonstandard':[r for r in residues if r['name'] not in standard],'metals':[r for r in hetero if r['name'] in {'ZN','MG','CA','FE','MN','CU','CO','NI','CD'}],'covalent_records':[l for l in source.read_text().splitlines() if l.startswith(('LINK  ','SSBOND'))],'protonation':'Meeko residue templates; default HIS tautomer; no pH microstate determination','confirmed':False}
def prepare(source,out,chains=('A',),policy=None,classroom=False):
 source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True);p=plan(source,chains,(policy or {}).get('model_index',0))
 write_json(out/'preparation_plan.json',p)
 if not classroom and not (policy and policy.get('confirmed') and policy.get('input_sha256')==p['input_sha256']):raise ValueError('REVIEW_REQUIRED: confirm matching preparation_plan.json')
 if source.suffix.lower()=='.pdbqt':
  (out/'receptor.pdbqt').write_bytes(source.read_bytes());write_json(out/'preparation.json',p);return p
 if p['altloc'] and not (policy or {}).get('altloc'):raise ValueError('ALTLOC_REVIEW_REQUIRED')
 parser=MMCIFParser(QUIET=True,auth_chains=True) if source.suffix.lower() in ('.cif','.mmcif') else PDBParser(QUIET=True)
 s=parser.get_structure('receptor',str(source));keep=set((policy or {}).get('keep_hetero',[]));chosen=(policy or {}).get('altloc','A')
 class Selection(Select):
  def accept_model(self,m):return m.id==p['model_index']
  def accept_chain(self,c):return c.id in chains
  def accept_residue(self,r):return r.id[0]==' ' or f'{r.parent.id}:{r.id[1]}{r.id[2].strip()}' in keep
  def accept_atom(self,a):return not a.is_disordered() or a.altloc==chosen
 io=PDBIO();io.set_structure(s);io.save(str(out/'receptor_prepared.pdb'),Selection(),preserve_atom_numbering=True)
 if classroom or (policy or {}).get('repair_missing_sidechains'):
  from pdbfixer import PDBFixer
  from openmm.app import PDBFile
  from openmm import Platform
  raw=out/'receptor_selected_original.pdb';raw.write_bytes((out/'receptor_prepared.pdb').read_bytes())
  fixer=PDBFixer(filename=str(raw),platform=Platform.getPlatformByName('CPU'))
  fixer.findMissingResidues();omitted=dict(fixer.missingResidues);fixer.missingResidues={}
  fixer.findMissingAtoms()
  added=[{'residue':r.chain.id+':'+r.id,'name':r.name,'atoms':[a.name for a in aa]} for r,aa in fixer.missingAtoms.items()]
  terminal=[{'residue':r.chain.id+':'+r.id,'atoms':aa} for r,aa in fixer.missingTerminals.items()]
  fixer.addMissingAtoms(seed=20261012)
  with (out/'receptor_prepared.pdb').open('w') as f:PDBFile.writeFile(fixer.topology,fixer.positions,f,keepIds=True)
  p['repair']={'method':'PDBFixer missing atoms only; no loop reconstruction','seed':20261012,'added_sidechain_atoms':added,'terminal_atoms':terminal,'omitted_missing_residues':str(omitted)}
 argv=['/opt/sbvs/env/bin/mk_prepare_receptor.py','-i',str(out/'receptor_prepared.pdb'),'-o',str(out/'receptor'),'-p','-j']
 if (policy or {}).get('residue_templates'):argv+=['--set_template',policy['residue_templates']]
 argv+=['--write_pdb',str(out/'receptor_hydrogenated.pdb')]
 r=subprocess.run(argv,capture_output=True,text=True,timeout=180)
 (out/'preparation.log').write_text(r.stdout+r.stderr)
 if r.returncode or not (out/'receptor.pdbqt').exists():raise RuntimeError('Meeko receptor preparation failed: '+r.stderr[-1800:])
 p.update(confirmed=True,classroom_policy=classroom,kept_hetero=sorted(keep),argv=argv,receptor_sha256=sha(out/'receptor.pdbqt'))
 write_json(out/'preparation.json',p);return p
