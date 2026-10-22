"""Fetch public inputs and run real adapters; no private inputs or CPU MD fallback."""
import argparse
import importlib.metadata
import platform
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from urllib.request import Request, urlopen
from recovery import Recovery, digest, write, lease, interruption_handlers, execute

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

def source_hashes():
    relevant=[p for folder in ('00_workbench/yuanhe_preflight','02_docking/src',
        '03_virtual_screening/src','03_virtual_screening/configs','05_molecular_dynamics/modules',
        '05_molecular_dynamics/config','deployment/public_minimal') for p in (ROOT/folder).rglob('*')
        if p.is_file() and p.suffix in ('.py','.json','.yaml','.yml','.sh','.lock','.txt')]
    return {p.relative_to(ROOT).as_posix():digest(p) for p in relevant}


def main(out, resume, lock_fd):
    check=subprocess.run([sys.executable,'-m','pip','check'],capture_output=True,text=True)
    if check.returncode: raise RuntimeError('Dependency conflicts: '+check.stdout+check.stderr)
    from openmm import Platform, System, VerletIntegrator, Context, unit
    system=System();system.addParticle(1)
    integrator=VerletIntegrator(0.001*unit.picoseconds)
    context=Context(system,integrator,Platform.getPlatformByName('CUDA'))
    context.setPositions([[0,0,0]])
    gpu={'platform': context.getPlatform().getName(),
         'device': context.getPlatform().getPropertyValue(context,'DeviceName')}
    del context, integrator
    sources=source_hashes()
    identity={'source':sources, 'output':str(out), 'python':sys.executable,
              'version':sys.version, 'interpreter_sha256':digest(Path(sys.executable)),
              'conda_records':{p.name:digest(p) for p in sorted((Path(sys.prefix)/'conda-meta').glob('*.json'))},
              'machine':platform.machine(), 'gpu':gpu,
              'cuda_visible_devices':os.environ.get('CUDA_VISIBLE_DEVICES'),
              'packages':sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions())}
    # JSON-normalize tuples before comparing to a persisted identity.
    identity=json.loads(json.dumps(identity))
    recovery=Recovery(out,identity,resume=resume)
    (out/'pip-check.txt').write_text(check.stdout+check.stderr)
    write(out/'gpu.json',gpu)
    write(out/'source_hashes.json',sources)
    print('CUDA ready:',gpu,flush=True)
    started=time.monotonic()
    recovery.phase('assets',lambda:fetch_assets(out))
    for phase in ('prepare','dock','screen','build','dynamics','analyze'):
        def run_phase(phase=phase):
            print('Running:',phase,flush=True)
            execute([sys.executable,'-u',str(PACKAGE/'phases.py'),str(out),phase],
                    out/(phase+'.log'),lock_fd)
        recovery.phase(phase,run_phase)
    audit=json.loads((out/'md_audit.json').read_text())
    if not audit['passed'] or audit['formal_duration_met']:
        raise RuntimeError('Unexpected MD audit result')
    recovery.finish()
    write(out/'artifact_hashes.json',recovery.state['artifacts'])
    write(out/'COMPLETED.json',{'status':'PUBLIC_ENGINEERING_PATH_PASSED',
        'dependency_check_passed':True,'native_md_audit_passed':True,'gpu':gpu,
        'last_invocation_seconds':time.monotonic()-started,'formal_scientific_acceptance':False,
        'generation_models_tested':False,'production_ps':0.2,
        'scope':'Public 1CRN + ethanol; explicit engineering MD control, not a qualified hit'})
    print('PASSED: public docking to short CUDA MD and analysis.',flush=True)
    print('Evidence:',out/'COMPLETED.json',flush=True)


def fetch_assets(out):
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

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    out=args.output.resolve()
    try:
        if not out.is_dir():
            raise ValueError('Output directory must exist; use run_public.sh')
        with lease(out) as lock_fd, interruption_handlers():
            main(out,args.resume,lock_fd)
    except (Exception, KeyboardInterrupt):
        # A refused resume must not overwrite the original failure/evidence.
        traceback.print_exc()
        sys.exit(1)
