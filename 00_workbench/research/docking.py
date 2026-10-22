"""Two CPU Vina passes with content-addressed, integrity-checked incremental results."""
import concurrent.futures, itertools, json, shutil, subprocess, sys
from pathlib import Path
from .config import fingerprint
from .resources import detect,disk_guard

def select(rows,limit,mode):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    best={}
    for row in rows:
        cid=row['compound_id']
        if cid not in best or row['vina_affinity_kcal_mol']<best[cid]['vina_affinity_kcal_mol']:best[cid]=row
    ordered=sorted(best.values(),key=lambda r:(r['vina_affinity_kcal_mol'],r['compound_id']))
    if mode=='score':return ordered[:limit]
    groups={}
    for row in ordered:
        scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(row['smiles']),includeChirality=True)
        groups.setdefault(scaffold or row['smiles'],[]).append(row)
    result=[]
    # Round-robin across actual chemical series. No invented activity probability.
    for layer in itertools.zip_longest(*groups.values()):
        result.extend(r for r in layer if r is not None)
        if len(result)>=limit:break
    return result[:limit]

def run(root,d,cfg,out,executable):
    from webui.workflow import atomic,digest,read
    from .qualification import qualify,identity,verified_file
    from rdkit import Chem
    from rdkit.Chem import AllChem
    sys.path.insert(0,str(root/'02_docking/src'))
    from docking import run_vina_job, Protocol, Box
    rc=cfg['research'];disk_guard(d,rc)
    rows,counts=qualify(root,d,cfg,read(d/'reinvent4/candidates.json',[]),out/'qualification')
    atomic(out/'qualification_counts.json',counts)
    if not rows:raise ValueError('No risk/acquisition-qualified molecules; pending evidence is not a pass')
    rows=rows[:rc['initial_limit']]
    receptors=[]
    supplied=rc['receptor_ensemble'] or [{'path':cfg['docking_receptor'],**cfg.get('receptor_review',{})}]
    for item in supplied:
        path=Path(item.get('path',''))
        if not path.is_absolute():path=d/path
        if not item.get('reviewed') or not path.is_file() or digest(path)!=item.get('sha256'):raise ValueError('Reviewed receptor and matching hash required')
        if item['sha256'] not in {r['sha256'] for r in receptors}:receptors.append({**item,'path':path})
    if not receptors:raise ValueError('No genuine reviewed receptor conformations')
    b=cfg.get('docking_box') or {};box=Box(tuple(b.get('center',[])),tuple(b.get('size',[])),'manual')
    if box.validate():raise ValueError('Invalid docking box')
    exe=executable(cfg['vina'],'Official AutoDock Vina')
    version=subprocess.run([exe,'--version'],capture_output=True,text=True,timeout=20,check=True).stdout.strip()
    if 'vina' not in version.lower():raise ValueError('Executable does not identify as Vina')
    exe_hash=digest(exe)
    atomic(out/'tool_versions.json',{'vina':{'path':exe,'version_output':version,'binary_sha256':exe_hash}})
    prepared=out/'ligands';prepared.mkdir(exist_ok=True)
    for row in rows:
        # One explicit canonical state. Never duplicate it to pretend multiple states exist.
        mol=Chem.AddHs(Chem.MolFromSmiles(row['smiles']));p=AllChem.ETKDGv3();p.randomSeed=cfg['seed']
        path=prepared/(row['compound_id']+'.sdf')
        if AllChem.EmbedMolecule(mol,p)!=0:raise ValueError('Conformer generation failed for '+row['compound_id'])
        with Chem.SDWriter(str(path)) as writer:writer.write(mol)
        shutil.copyfile(row['qualification_evidence'],str(path)+'.qualification.json')
        row['sdf']=str(path.relative_to(d));row['state_id']=fingerprint(row['smiles'])[:20]
        row['states']=[{'sdf':row['sdf'],'state_id':row['state_id'],'smiles':row['smiles']}]
        for state in rc['ligand_states'].get(row['compound_id'],[]):
            if not state.get('reviewed') or not state.get('rationale'):raise ValueError('Chemical state requires explicit review')
            source=verified_file(root,state)
            molecules=list(Chem.SDMolSupplier(str(source),removeHs=True))
            if len(molecules)!=1 or molecules[0] is None:raise ValueError('One valid molecule required for a reviewed chemical state')
            exact=Chem.MolToSmiles(molecules[0],isomericSmiles=True)
            if identity(exact)[0]!=row['compound_id']:raise ValueError('Chemical state belongs to a different parent or stereoisomer')
            sid=fingerprint(exact)[:20]
            if any(x['state_id']==sid for x in row['states']):continue
            target=prepared/(row['compound_id']+'_'+sid+'.sdf');shutil.copyfile(source,target)
            evidence=(root/state.get('qualification_evidence','')).resolve()
            from yuanhe_preflight.gate import require_ligand
            require_ligand(target,evidence)
            shutil.copyfile(evidence,str(target)+'.qualification.json')
            row['states'].append({'sdf':str(target.relative_to(d)),'state_id':sid,'smiles':exact,'review':state})
    resources=detect(d);workers=min(rc['vina_workers'],resources['effective_cpus']//rc['vina_cpu'])
    if workers<1:raise ValueError('CPU quota is smaller than one Vina worker')
    atomic(out/'resource_plan.json',{'detected':resources,'workers':workers,'threads_per_worker':rc['vina_cpu']})
    cache=d/'research_cache/docking';cache.mkdir(parents=True,exist_ok=True)
    all_results=[];actual={}
    for phase in ('initial','refine'):
        settings=rc[phase];phase_out=out/phase;phase_out.mkdir(exist_ok=True)
        targets=receptors[:settings['receptors']]
        proto=Protocol.preset('standard',seeds=(cfg['seed'],),exhaustiveness=settings['exhaustiveness'],num_modes=settings['num_modes'],energy_range=cfg['energy_range'])
        def job(task):
            parent,state,receptor,seed=task
            row={**parent,**state,'parent_smiles':parent['smiles']}
            inputs={'compound_id':row['compound_id'],'smiles':row['smiles'],'ligand_sha256':digest(d/row['sdf']),
                    'state_id':row['state_id'],'receptor_sha256':receptor['sha256'],'box':b,'seed':seed,
                    'exhaustiveness':settings['exhaustiveness'],'num_modes':settings['num_modes'],'energy_range':cfg['energy_range'],
                    'engine':'vina','version':version,'executable_sha256':exe_hash,'cpu':rc['vina_cpu'],
                    'adapter_sha256':digest(root/'02_docking/src/sbvs2/docking/vina_official.py')}
            from yuanhe_preflight.gate import require_ligand
            require_ligand(d/row['sdf'])
            key=fingerprint(inputs);folder=cache/key;manifest=folder/'reuse.json'
            prior=read(manifest)
            if prior:
                if prior['inputs']!=inputs or any(not (folder/p).is_file() or digest(folder/p)!=sha for p,sha in prior['files'].items()):raise ValueError('Cached docking integrity mismatch')
                result=prior['result']
            else:
                disk_guard(d,rc)
                if folder.exists():raise ValueError('Incomplete docking cache retained; inspect before retry: '+str(folder))
                result=run_vina_job(exe,receptor['path'],d/row['sdf'],box,proto,seed,folder,cpu=rc['vina_cpu'],timeout=cfg['timeout_seconds'])
                if result['status']=='ok':
                    files={str(p.relative_to(folder)):digest(p) for p in folder.rglob('*') if p.is_file()}
                    atomic(manifest,{'inputs':inputs,'result':result,'files':files})
            entry={**row,'job_id':key,'seed':seed,'receptor_sha256':receptor['sha256'],'receptor_file':str(receptor['path']),'status':result['status'],
                   'engine':'vina','engine_version':version,'cache_path':str(folder.relative_to(d)),'phase':phase}
            if result['status']=='ok':
                entry.update(min(result['poses'],key=lambda p:p['vina_affinity_kcal_mol']))
                from yuanhe.contracts import docking_scores
                docking_scores(engine='vina',vina_score=entry['vina_affinity_kcal_mol'],cnn_score=entry.get('cnn_pose_score'),cnn_affinity=entry.get('cnn_affinity_pK'))
            return entry,len(result.get('poses',[]))
        tasks=((row,state,receptor,cfg['seed']+i) for row in rows for state in row['states'][:settings['states']] for receptor in targets for i in range(settings['seeds']))
        actual_states=sum(min(len(row['states']),settings['states']) for row in rows)
        results=[];poses=0;processed=0
        with (phase_out/'results.jsonl').open('w') as journal, concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            pending=set()
            for _ in range(workers*2):
                task=next(tasks,None)
                if task is not None:pending.add(pool.submit(job,task))
            while pending:
                done,pending=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    entry,nposes=future.result();poses+=nposes;processed+=1
                    journal.write(json.dumps(entry)+'\n');journal.flush()
                    if entry['status']=='ok':results.append(entry)
                    atomic(out/'progress.json',{'input':actual_states*len(targets)*settings['seeds'],'processed':processed,'successful':len(results),'phase':phase})
                    task=next(tasks,None)
                    if task is not None:pending.add(pool.submit(job,task))
        actual[phase]={'parents':len(rows),'chemical_states':actual_states,'receptors':len(targets),'tasks':processed,'poses':poses,
                       'requested_states_per_parent':settings['states'],'actual_chemical_states':actual_states,
                       'limitation':'Additional chemical states require explicit preparation; no synthetic duplicates'}
        limit=rc['refine_limit'] if phase=='initial' else rc['review_limit']
        selected=select(results,limit,rc['selection_mode'])
        atomic(phase_out/'selection.json',selected)
        atomic(phase_out/'score_selection.json',select(results,limit,'score'))
        all_results=results;rows=[{**r,'smiles':r.get('parent_smiles',r['smiles'])} for r in selected]
        if not rows:raise ValueError('No successful Vina results in '+phase)
    selected_ids={r['compound_id'] for r in rows}
    handoff_results=[r for r in all_results if r['compound_id'] in selected_ids]
    atomic(out/'counts.json',actual);atomic(out/'results.json',handoff_results)
    atomic(out/'input_molecules.json',[{**r,**state} for r in rows for state in r['states']]);atomic(out/'box.json',b);shutil.copyfile(receptors[0]['path'],out/'receptor.pdb')
    # Preserve the legacy Y-MAVS handoff contract; poses retain original cache provenance.
    for row in all_results:
        dest=out/'jobs'/row['job_id']
        if not dest.exists():shutil.copytree(d/row['cache_path'],dest)
    import csv
    with (out/'molecules.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['id','compound_id','smiles','status','vina_affinity_kcal_mol'],extrasaction='ignore');writer.writeheader();writer.writerows(all_results)
    return {'input':counts['unique'],'remaining':len(rows),'poses':sum(x['poses'] for x in actual.values()),**counts}
