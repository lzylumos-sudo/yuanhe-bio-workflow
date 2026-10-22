#!/usr/bin/env python3
"""Run real Vina only on previously screened, route-qualified candidates.

Requires an actual completed screening run and installed Vina; no injected data.
This verifies the executable handoff, not target affinity or experimental synthesis.
"""
import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'02_docking/src'))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--screening',type=Path,required=True)
    ap.add_argument('--vina',type=Path,required=True)
    ap.add_argument('--receptor',type=Path,required=True)
    ap.add_argument('--box',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    from yuanhe_preflight.common import read,write,sdf_smiles
    from yuanhe_preflight.gate import require_ligand
    from yuanhe_preflight.ligand import prepare_selected
    from sbvs2.docking.vina_official import run_vina_job
    from sbvs2.docking.box import Box
    from sbvs2.docking.protocol import Protocol
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    selected=read(args.screening/'qualified.json');records=read(args.screening/'records.json')
    if not selected:raise ValueError('No synthesis-qualified molecules; no docking will run')
    for row in selected:require_ligand(row['source_sdf'],row['qualification_evidence'])
    selected_ids={r['compound_id'] for r in selected}
    assert selected_ids=={r['compound_id'] for r in records if r['decision']=='pass'}
    prepared=prepare_selected(selected,out,out/'ligands',17)
    rb=read(args.box);box=Box(tuple(rb['center']),tuple(rb['size']),'manual')
    proto=Protocol.preset('standard',seeds=(17,),exhaustiveness=4,num_modes=3,energy_range=3)
    results=[]
    for row in prepared:
        assert sdf_smiles(row['sdf'])==row['smiles']
        result=run_vina_job(str(args.vina.resolve()),args.receptor.resolve(),row['sdf'],box,proto,17,
                            out/row['compound_id'],cpu=2,timeout=120,qualification_evidence=row['qualification_evidence'])
        results.append({'molecule_id':row['compound_id'],'smiles':row['smiles'],'result':result})
        write(out/'results.json',results)
        assert result['status']=='ok',result
    write(out/'acceptance.json',{'input_records':len(records),'selected':len(selected),'docked':len(results),
                                'excluded_not_docked':len(records)-len(selected),'results':results,
                                'scope':'Real inference and computed routes in input; real conformers and Vina; no experimental claim'})

if __name__=='__main__':main()
