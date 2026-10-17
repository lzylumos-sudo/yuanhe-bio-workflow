#!/usr/bin/env python3
"""Independent, local-only complex -> OpenMM -> quantitative analysis workflow."""
from __future__ import annotations
import argparse, copy, hashlib, json, os, shutil, sys, traceback, uuid
from datetime import datetime, timezone
from pathlib import Path
import yaml

PROJECT_ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT_ROOT))
from modules.utils.core import (WorkflowError, atomic_json, read_json, digest,
    now, run_paths, load_config, relative_path, input_files, config_hash)
from modules.utils.complex_input import AMINO, WATER, residue_key
from modules.utils.embedded_chemistry import read_embedded, EmbeddedChemistryError

def merge(left,right):
    for key,value in right.items():
        if isinstance(value,dict) and isinstance(left.get(key),dict): merge(left[key],value)
        else: left[key]=copy.deepcopy(value)
    return left

def implementation():
    names=['modules/utils/core.py','modules/utils/complex_input.py',
           'modules/utils/preparation.py','modules/utils/embedded_chemistry.py','modules/md/openmm_engine.py','modules/md/amber_export.py']
    return {name:digest(PROJECT_ROOT/name) for name in names}

def read_manifest(path):
    for name in ['manifest.json','complex_manifest.json','manifest.yaml','complex_manifest.yaml']:
        file=path.parent/name
        if file.is_file():
            data=json.loads(file.read_text()) if file.suffix=='.json' else yaml.safe_load(file.read_text())
            if not isinstance(data,dict): raise WorkflowError('Complex manifest must be a mapping.')
            return data
    return {}

def stage_input(source,run_id,settings):
    source=Path(source).expanduser().resolve()
    if not source.is_file(): raise WorkflowError('complex.pdb input does not exist.')
    paths=run_paths(run_id); root=paths['md'].parent
    root.mkdir(parents=True,exist_ok=False)
    inputs=root/'inputs'; inputs.mkdir()
    dest=inputs/'complex.pdb'; shutil.copyfile(source,dest)
    pdb_text=source.read_text()
    try: embedded=read_embedded(pdb_text)
    except EmbeddedChemistryError as exc: raise WorkflowError(str(exc)) from exc
    sidecar_info=read_manifest(source)
    sidecar_manifest_present=any((source.parent/name).is_file() for name in ['manifest.json','complex_manifest.json','manifest.yaml','complex_manifest.yaml'])
    if embedded is not None and sidecar_manifest_present and sidecar_info!=embedded['manifest']:
        raise WorkflowError('Embedded and adjacent manifests disagree; use one validated chemical definition.')
    info=embedded['manifest'] if embedded is not None else sidecar_info
    records=[s for s in pdb_text.splitlines() if s.startswith(('ATOM  ','HETATM'))]
    keys={residue_key(line) for line in records}
    candidates=sorted(k for k in keys if k[3] not in AMINO|WATER)
    chains=sorted({k[0] for k in keys if k[3] in AMINO})
    config=copy.deepcopy(settings)
    config.update(project_name='Nexus_MD',workflow='md_from_complex_v1')
    config['complex']={'file':dest.relative_to(PROJECT_ROOT).as_posix()}
    ligand=copy.deepcopy(info.get('ligand'))
    if candidates:
        if ligand is None:
            if len(candidates)!=1: raise WorkflowError('Multiple nonstandard residues: provide an explicit ligand manifest and parameters; metals/cofactors are not inferred.')
            chain,resid,icode,resname=candidates[0]
            ligand={'chain':chain,'resid':resid,'insertion_code':icode,'resname':resname}
        template=Path(ligand.get('template_sdf','ligand.sdf'))
        if template.is_absolute() or '..' in template.parts:
            raise WorkflowError('Manifest template_sdf must be a relative path next to complex.pdb.')
        template=(source.parent/template).resolve()
        if not template.is_relative_to(source.parent):
            raise WorkflowError('Manifest template escapes the complex directory.')
        embedded_sdf=embedded.get('ligand_sdf') if embedded is not None else None
        if embedded_sdf is not None:
            if template.is_file() and template.read_text()!=embedded_sdf:
                raise WorkflowError('Embedded and adjacent SDF chemical definitions disagree.')
            (inputs/'ligand.sdf').write_text(embedded_sdf)
        elif template.is_file():
            shutil.copyfile(template,inputs/'ligand.sdf')
        else:
            raise WorkflowError('Unknown ligand chemistry: use a Nexus complex.pdb with embedded chemistry or keep a precise ligand.sdf beside it. Ordinary PDB coordinates alone do not define bond orders/protonation/stereochemistry.')
        ligand['template_sdf']=(inputs/'ligand.sdf').relative_to(PROJECT_ROOT).as_posix()
    elif ligand is not None:
        raise WorkflowError('Manifest declares a ligand but complex.pdb contains no nonstandard residue.')
    config['ligand']=ligand
    config['system_type']='ternary' if ligand and len(chains)>1 else ('protein_ligand' if ligand else 'protein_protein')
    if info.get('system_type'):
        config['system_type']=info['system_type']
    if 'analysis' in info: merge(config['analysis'],info['analysis'])
    if 'groups' in info: config['analysis']['groups']=info['groups']
    config['input_provenance']={'complex_sha256':digest(dest),'source_filename':source.name,
        'manifest_present':bool(info),'embedded_chemistry':embedded is not None,'automatic_chain_groups':not bool(config['analysis'].get('groups')),
        'chemistry_source':('embedded SDF' if embedded is not None else 'sidecar SDF') if ligand else 'standard protein residues'}
    config_path=root/'config.yaml'
    config_path.write_text(yaml.safe_dump(config,sort_keys=False))
    return load_config(config_path.relative_to(PROJECT_ROOT).as_posix())

def inventory(stage,paths):
    if stage=='prepare': files=list((paths['md']/'prepared').glob('*'))
    elif stage=='build': files=list((paths['md']/'system').glob('*'))
    elif stage=='md': files=list(paths['md'].glob('replica_*/progress.json'))
    elif stage=='analysis': files=list(paths['analysis'].rglob('*'))
    else: files=[]
    return {p.relative_to(PROJECT_ROOT).as_posix():digest(p) for p in files if p.is_file()}

def verify_files(record):
    for name,expected in record.get('files',{}).items():
        if digest(relative_path(name,must_exist=True))!=expected:
            raise WorkflowError('Committed artifact changed: '+name)

def acceptance(config,paths,run_id):
    from modules.md.acceptance import verify_run
    result=verify_run(config,paths)
    atomic_json(paths['md'].parent/'acceptance.json',result)
    enough=result['production_ns']>=10.0-1e-8 and config['replicates']>=1
    text=['# MD test '+('PASS' if result['passed'] else 'FAIL'),'','Scope: numerical workflow acceptance; no convergence, affinity or pose-correctness claim.',
          '','Run ID: '+run_id,'Production per replica (ns): '+str(result['production_ns']),
          'Replicas: '+str(config['replicates']),
          'Force fields: '+config['MD']['protein_forcefield']+' / '+config['MD']['water_forcefield']+' / '+str(config['MD'].get('ligand_forcefield')),
          'MM/GBSA: '+result['mmgbsa_status']+'; entropy excluded.','',json.dumps(result,indent=2),'']
    filename='MD_TEST_PASS.md' if enough and result['passed'] else 'SHORT_RUN_TEST.md'
    (paths['md'].parent/filename).write_text('\n'.join(text))
    if enough and result['passed']:
        shutil.copyfile(paths['md'].parent/filename,PROJECT_ROOT/'MD_TEST_PASS.md')
    atomic_json(PROJECT_ROOT/'results/latest.json',{'run_id':run_id,'path':paths['md'].parent.relative_to(PROJECT_ROOT).as_posix(),'accepted':result['passed']})
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('complex',nargs='?',help='complex.pdb; optional ligand.sdf and manifest beside it')
    p.add_argument('--config',default='config/default.yaml',help='module-relative YAML settings')
    p.add_argument('--run-id'); p.add_argument('--resume',action='store_true')
    p.add_argument('--until',choices=['prepare','build','md','analysis','acceptance'],default='acceptance')
    p.add_argument('--analyze-only',action='store_true')
    args=p.parse_args()
    run_id=args.run_id or datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_')+uuid.uuid4().hex[:6]
    paths=run_paths(run_id); folder=paths['md'].parent
    if args.resume or args.analyze_only:
        if not args.run_id: raise WorkflowError('--resume/--analyze-only requires --run-id.')
        config=load_config((folder/'config.yaml').relative_to(PROJECT_ROOT).as_posix())
    else:
        if not args.complex: raise WorkflowError('Supply complex.pdb or resume an existing run.')
        settings=yaml.safe_load((PROJECT_ROOT/'config/default.yaml').read_text())
        if args.config!='config/default.yaml': merge(settings,yaml.safe_load(relative_path(args.config,must_exist=True).read_text()))
        config=stage_input(args.complex,run_id,settings)
    paths['md'].mkdir(parents=True,exist_ok=True)
    lock=folder/'.running'
    try: lock.mkdir()
    except FileExistsError: raise WorkflowError('Run is locked; verify its recorded PID is no longer active before removing the stale lock.')
    atomic_json(lock/'owner.json',{'pid':os.getpid(),'run_id':run_id,'created':now()})
    manifest_path=folder/'manifest.json'
    current_inputs={name:digest(relative_path(name,must_exist=True)) for name in input_files(config)}
    manifest=read_json(manifest_path) if manifest_path.exists() else {'created':now(),'run_id':run_id,'stages':{},'input_sha256':current_inputs,'config_sha256':config_hash(config),'md_implementation_sha256':implementation()}
    try:
        if current_inputs!=manifest['input_sha256'] or config_hash(config)!=manifest['config_sha256']:
            raise WorkflowError('Inputs or scientific settings changed; create a new run.')
        if implementation()!=manifest['md_implementation_sha256']:
            raise WorkflowError('MD implementation changed; existing trajectories remain intact. Use a fresh run or an explicit migration audit.')
        manifest.update(status='running',updated=now())
        atomic_json(manifest_path,manifest)
        from modules.utils.complex_input import prepare_complex
        from modules.md.openmm_engine import build_system,run_md
        from modules.analysis.trajectory import analyze
        from modules.analysis.mmgbsa import analyze_mmgbsa
        stages=['prepare','build','md','analysis']
        for stage in stages:
            if args.analyze_only and stage!='analysis':
                record=manifest['stages'].get(stage,{})
                if record.get('status')!='complete': raise WorkflowError('Analysis requires a completed MD run.')
                verify_files(record)
                if stage=='md': run_md(config,paths)
                continue
            record=manifest['stages'].get(stage,{})
            if record.get('status')=='complete' and not (args.analyze_only and stage=='analysis'):
                verify_files(record)
                if stage=='md': run_md(config,paths)
            else:
                manifest['stages'][stage]={'status':'running','started':now()}; atomic_json(manifest_path,manifest)
                print('Starting '+stage,flush=True)
                if stage=='prepare': prepare_complex(config,paths)
                elif stage=='build': build_system(config,paths)
                elif stage=='md': run_md(config,paths)
                else:
                    analyze(config,paths)
                    analyze_mmgbsa(config,paths)
                manifest['stages'][stage].update(status='complete',completed=now(),files=inventory(stage,paths))
                atomic_json(manifest_path,manifest)
            if args.until==stage:
                manifest.update(status='partial',updated=now()); atomic_json(manifest_path,manifest)
                print(json.dumps({'run_id':run_id,'completed_stage':stage})); return
        result=acceptance(config,paths,run_id)
        manifest.update(status='complete',completed=now(),accepted=result['passed']); atomic_json(manifest_path,manifest)
        print(json.dumps({'run_id':run_id,'status':'COMPLETE','production_ns':result['production_ns'],'mmgbsa':result['mmgbsa_status']}))
    except BaseException as exc:
        manifest.update(status='interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed',updated=now(),error_type=type(exc).__name__)
        atomic_json(manifest_path,manifest)
        raise
    finally:
        (lock/'owner.json').unlink(missing_ok=True); lock.rmdir()

if __name__=='__main__':
    try: main()
    except WorkflowError as exc:
        print('ERROR: '+str(exc),file=sys.stderr); sys.exit(2)
