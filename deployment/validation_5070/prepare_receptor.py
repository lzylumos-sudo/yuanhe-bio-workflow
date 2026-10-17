"""Explicit optional H addition for the observed MurJ construct; not a membrane builder."""
import argparse,json,hashlib
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--ph',type=float,default=7.4);a=p.parse_args()
 if a.out.exists() or a.out.with_suffix('.review.json').exists():raise SystemExit('Output already exists.')
 from pdbfixer import PDBFixer
 from openmm.app import PDBFile
 from openmm import Platform,unit
 if not 0<=a.ph<=14:raise SystemExit('Invalid pH')
 fixer=PDBFixer(filename=str(a.input),platform=Platform.getPlatformByName('CPU'))
 residues=list(fixer.topology.residues())
 if {r.chain.id for r in residues}!={'A'}:raise SystemExit('Expected receptor chain A only.')
 ids={int(r.id) for r in residues}
 if ids!=set(range(3,501)):raise SystemExit('Expected observed MurJ author residues 3..500; inspect missing residues/selection first.')
 fixer.findMissingResidues();fixer.findNonstandardResidues();fixer.findMissingAtoms()
 report={'input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'pH_assumption':a.ph,'construct':'9NU4 K5F/S12I/M13A, observed residues 3..500; not WT','missing_residues':str(fixer.missingResidues),'nonstandard_residues':str(fixer.nonstandardResidues),'missing_heavy_atoms':{r.id:[x.name for x in xs] for r,xs in fixer.missingAtoms.items()},'missing_terminal_atoms':{r.id:list(xs) for r,xs in fixer.missingTerminals.items()},'scientific_review':'required; inspect histidines, acidic/basic pocket residues, artificial termini and atom geometry'}
 if fixer.missingResidues or fixer.nonstandardResidues or fixer.missingAtoms:
  a.out.with_suffix('.review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
  raise SystemExit('Structural issues detected; no silent loop/side-chain reconstruction. Review the JSON first.')
 # Adding missing terminal atoms is explicit in the record, not restoring unresolved loops.
 before={(x.residue.id,x.name):pos.value_in_unit(unit.nanometer) for x,pos in zip(fixer.topology.atoms(),fixer.positions) if x.element and x.element.symbol!='H'}
 fixer.addMissingAtoms(seed=17)
 fixer.addMissingHydrogens(a.ph)
 # Keep observed heavy-atom frame so peptide and proposed box stay aligned.
 pos=fixer.positions.value_in_unit(unit.nanometer)
 for atom in fixer.topology.atoms():
  key=(atom.residue.id,atom.name)
  if key in before:pos[atom.index]=before[key]
 fixer.positions=pos*unit.nanometer
 with a.out.open('w') as f:PDBFile.writeFile(fixer.topology,fixer.positions,f,keepIds=True)
 report['output_sha256']=hashlib.sha256(a.out.read_bytes()).hexdigest()
 report['heavy_atoms_restored_to_input_frame']=True
 a.out.with_suffix('.review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
 print('Created receptor candidate. Review before docking; this is not a scientifically approved receptor.')

if __name__=='__main__':main()
