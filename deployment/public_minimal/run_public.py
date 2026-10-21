"""Fetch public inputs and run real adapters; no private inputs or CPU MD fallback."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from urllib.request import Request, urlopen

PACKAGE = Path(__file__).resolve().parent
ROOT = PACKAGE.parents[1]
ASSETS = {
    '1CRN.pdb': ('https://files.rcsb.org/download/1CRN.pdb',
                 '42199a30a0701864a2a5cc76cd7f35cc544cd0e65fbcf63e03c166543249b811'),
    'vina_1.2.7_linux_x86_64': (
        'https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64',
        'f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644'),
    'ethanol_supplier.html': (
        'https://molekula.com/catalog/64-17-5/90005199-ethanol-absolute-p-a?locale=en', None),
}

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def write(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')

def main(out):
    # The shell creates only isolated caches; refuse reuse of scientific outputs.
    if any((out/name).exists() for name in ('sources','dynamics','COMPLETED.json')):
        raise ValueError('A fresh output directory is required')
    check=subprocess.run([sys.executable,'-m','pip','check'],capture_output=True,text=True)
    (out/'pip-check.txt').write_text(check.stdout+check.stderr)
    if check.returncode: raise RuntimeError('Dependency conflicts: see pip-check.txt')
    from openmm import Platform, System, VerletIntegrator, Context, unit
    system=System();system.addParticle(1)
    integrator=VerletIntegrator(0.001*unit.picoseconds)
    context=Context(system,integrator,Platform.getPlatformByName('CUDA'))
    context.setPositions([[0,0,0]])
    gpu={'platform': context.getPlatform().getName(),
         'device': context.getPlatform().getPropertyValue(context,'DeviceName')}
    del context, integrator
    write(out/'gpu.json',gpu)
    print('CUDA ready:',gpu,flush=True)
    src=out/'sources';src.mkdir()
    records={}
    for name,(url,expected) in ASSETS.items():
        print('Fetching public asset:',name,flush=True)
        request=Request(url,headers={'User-Agent':'Yuanhe-Public-Engineering-Validation/1.0'})
        with urlopen(request,timeout=90) as response:
            path=src/name
            with path.open('wb') as f:
                while chunk:=response.read(1024*1024):f.write(chunk)
            records[name]={'requested_url':url,'resolved_url':response.url,
                           'fetched_at_utc':datetime.now(timezone.utc).isoformat(),'sha256':digest(path)}
        if expected and records[name]['sha256']!=expected:
            raise RuntimeError('Public asset checksum mismatch: '+name)
    (src/'vina_1.2.7_linux_x86_64').chmod(0o755)
    write(out/'public_sources.json',records)
    # Record the actual code used, including local edits, without reading user data.
    relevant=[p for folder in ('00_workbench/yuanhe_preflight','02_docking/src',
        '03_virtual_screening/src','03_virtual_screening/configs','05_molecular_dynamics/modules',
        '05_molecular_dynamics/config','deployment/public_minimal') for p in (ROOT/folder).rglob('*')
        if p.is_file() and p.suffix in ('.py','.yaml','.yml','.sh','.lock','.txt')]
    write(out/'source_hashes.json',{p.relative_to(ROOT).as_posix():digest(p) for p in relevant})
    started=time.monotonic()
    for phase in ('prepare','dock','screen','build','dynamics','analyze'):
        print('Running:',phase,flush=True)
        with (out/(phase+'.log')).open('w') as log:
            subprocess.run([sys.executable,'-u',str(PACKAGE/'phases.py'),str(out),phase],
                stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    audit=json.loads((out/'md_audit.json').read_text())
    if not audit['passed'] or audit['formal_duration_met']:
        raise RuntimeError('Unexpected MD audit result')
    manifest={p.relative_to(out).as_posix():digest(p) for p in out.rglob('*')
              if p.is_file() and not set(p.relative_to(out).parts)&{'cache','scratch'}}
    write(out/'artifact_hashes.json',manifest)
    write(out/'COMPLETED.json',{'status':'PUBLIC_ENGINEERING_PATH_PASSED',
        'dependency_check_passed':True,'native_md_audit_passed':True,'gpu':gpu,
        'compute_seconds':time.monotonic()-started,'formal_scientific_acceptance':False,
        'generation_models_tested':False,'production_ps':0.2,
        'scope':'Public 1CRN + ethanol; explicit engineering MD control, not a qualified hit'})
    print('PASSED: public docking to short CUDA MD and analysis.',flush=True)
    print('Evidence:',out/'COMPLETED.json',flush=True)

if __name__=='__main__':
    out=Path(sys.argv[1]).resolve()
    try: main(out)
    except Exception as exc:
        write(out/'FAILED.json',{'error':str(exc),'type':type(exc).__name__})
        traceback.print_exc()
        sys.exit(1)
