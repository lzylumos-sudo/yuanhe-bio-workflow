#!/usr/bin/env python3
"""Real public receptor preparation beside a controlled generation-slot worker.

Requires Linux/WSL, PDBFixer and OpenMM. Only use public, nonconfidential fixtures.
The generation-slot worker is deliberately NOT model inference.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]
UI=ROOT/'00_workbench'
sys.path.insert(0,str(UI))
from webui import workflow as w,studio


def wait_for(check, timeout=180):
    deadline=time.monotonic()+timeout
    while time.monotonic()<deadline:
        if check():return
        time.sleep(.1)
    raise TimeoutError('Acceptance worker did not reach the expected state')


def validate(protein, output):
    if sys.platform!='linux':raise RuntimeError('Run in Linux/WSL')
    from importlib.metadata import version
    report={'scope':'public receptor and orchestration engineering only',
            'generation':'controlled long-lived worker; no model inference',
            'versions':{m:version(m) for m in ('pdbfixer','openmm')},
            'input_sha256':w.digest(protein),'checks':{}}
    old=(w.ROOT,w.UI,w.RUNS)
    with tempfile.TemporaryDirectory(prefix='yuanhe-concurrent-acceptance-') as temp:
        root=Path(temp);w.ROOT=root;w.UI=root;w.RUNS=root/'runs'
        wrapper="import sys,runpy\nfrom pathlib import Path\nsys.path.insert(0,"+repr(str(UI))+ ")\nfrom webui import workflow as w\nw.ROOT=Path("+repr(str(root))+ ")\nw.UI=w.ROOT\nw.RUNS=w.ROOT/'runs'\nrunpy.run_path("+repr(str(UI/'studio_worker.py'))+",run_name='__main__')\n"
        (root/'studio_worker.py').write_text(wrapper)
        (root/'pipeline_bridge.py').write_text("import time\ntime.sleep(300)\n")
        rid='public_acceptance';d=w.run_dir(rid)
        cfg=w.defaults();cfg.update(research={},predocking={},prep_python=sys.executable,chem_python=sys.executable)
        w.atomic(d/'config.json',cfg);w.atomic(d/'run.json',{'config_sha256':w.digest(d/'config.json')})
        (d/'prepare').mkdir();(d/'prepare/receptor_full.pdb').write_bytes(protein.read_bytes())
        w.atomic(d/'prepare/stage.json',{'status':'completed','config_sha256':w.digest(d/'config.json'),
                'artifacts':{'prepare/receptor_full.pdb':report['input_sha256']}})
        try:
            stage=w.start(rid,'peptide2mol')
            proc=w._active[(rid,'stage','peptide2mol')]
            frozen=w.digest(d/'peptide2mol/review.json')
            op=studio.launch(rid,'receptor',{'pH':7.4,'repair_sidechains':True})['id'];out=d/'operations'/op
            wait_for(lambda:(rid,'operation',op) not in w._active)
            record=w.read(out/'operation.json')
            if record['status']!='completed':raise RuntimeError(json.dumps(record))
            assert proc.poll() is None,'Generation-slot worker exited before receptor completion'
            result=w.read(out/'result.json');receptor=d/result['receptor']
            assert receptor.is_file() and w.digest(receptor)==result['sha256']
            assert record['resources']['CUDA_VISIBLE_DEVICES']==''
            assert int(record['resources']['OPENMM_CPU_THREADS'])<=2
            assert studio.apply_receptor(rid,{'operation':op,'reviewed':True,'box':{'center':[0,0,0],'size':[20]*3}})['saved']
            assert w.digest(d/'peptide2mol/review.json')==frozen
            assert w.digest(d/'prepare/receptor_full.pdb')==report['input_sha256']
            before=studio.structure(protein);after=studio.structure(receptor)
            def heavy(atoms):return {(a['chain'],a['resid'],a['name']):a['p'] for a in atoms if a['element']!='H'}
            original=heavy(before);prepared=heavy(after)
            assert all(key in prepared and max(abs(x-y) for x,y in zip(coords,prepared[key]))<.002 for key,coords in original.items())
            assert len(after)>len(before),'Hydrogens were not added'
            report.update(status='PASSED',receptor_sha256=result['sha256'],input_atoms=len(before),prepared_atoms=len(after),
                          receptor_seconds=record['ended']-record['started'],
                          checks={'overlapping_processes':True,'cpu_only':True,'bounded_cpu_threads':True,
                                  'independent_completion':True,'apply_during_generation':True,
                                  'generation_input_unchanged':True,'receptor_input_unchanged':True,
                                  'existing_heavy_coordinates_preserved':True})
        finally:
            w.cancel(rid);wait_for(lambda:not w._active,40)
            w.ROOT,w.UI,w.RUNS=old
    w.atomic(output,report)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--public-pdb',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();validate(args.public_pdb,args.output)
