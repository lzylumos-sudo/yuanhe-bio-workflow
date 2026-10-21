"""Derived selections and wet-lab handoff keep source runs immutable."""
import json,re,shutil,time,zipfile
from .config import resolve,plan

def derive(rid,changes):
    from webui import workflow as w
    parent=w.run_dir(rid);w.verify_source(parent)
    cfg=w.read(parent/'config.json');rc=cfg.get('research')
    if not rc:raise ValueError('Legacy runs cannot be upgraded implicitly')
    allowed={'refine_limit','review_limit','md_limit','wetlab_min','wetlab_max','selection_mode'}
    if set(changes)-allowed:raise ValueError('Reuse currently supports selection-only changes; create a new run for scientific input changes')
    rc=resolve({**rc,**changes});cfg={**cfg,'research':rc};cfg.pop('input_sha256',None)
    # Validate all upstream artifacts before creating a derived run.
    for stage in ('prepare','peptide2mol','reinvent4'):w.require(parent,stage)
    child=w.create(cfg)['id'];d=w.run_dir(child)
    for stage in ('prepare','peptide2mol','reinvent4'):
        shutil.copytree(parent/stage,d/stage)
        rec=w.read(d/stage/'stage.json');rec.update(config_sha256=w.digest(d/'config.json'),reused_from=rid,source_stage_sha256=w.digest(parent/stage/'stage.json'))
        w.atomic(d/stage/'stage.json',rec)
    # Copy only verified content-addressed cache entries; never trust an exit code alone.
    for entry in (parent/'research_cache/docking').glob('*'):
        rec=w.read(entry/'reuse.json')
        if not rec:continue
        if any(not (entry/name).is_file() or w.digest(entry/name)!=sha for name,sha in rec['files'].items()):raise ValueError('Parent docking cache changed')
        shutil.copytree(entry,d/'research_cache/docking'/entry.name)
    if (parent/'reviews').exists():shutil.copytree(parent/'reviews',d/'reviews')
    proof={'parent_run':rid,'parent_config_sha256':w.digest(parent/'config.json'),'changes':changes,
           'generation_reused':True,'docking_reuse':'verified cache entries; exact scientific key required',
           'additional_docking_tasks':'computed from cache misses at dispatch','plan':plan(rc)}
    w.atomic(d/'derivation.json',proof)
    return {'id':child,**proof}

def md_progress(d):
    from webui import workflow as w
    state=w.read(d/'md_queue.json',{'tasks':{}})
    for ident,row in state['tasks'].items():
        folder=d/'md_systems'/ident;log=folder/'production.log'
        row['ns_per_day']=None;row['eta_days']=None;row['production_ns']=row.get('completed_ns',0)
        if log.is_file():
            with log.open('rb') as f:f.seek(max(0,log.stat().st_size-128000));text=f.read().decode(errors='replace')
            matches=re.findall(r'Step\s+Time\s*\n\s*\d+\s+([\d.eE+-]+)',text)
            if matches:row['production_ns']=float(matches[-1])/1000
            speed=re.findall(r'Performance:\s+([\d.]+)',text)
            if speed and float(speed[-1])>0:
                row['ns_per_day']=float(speed[-1]);row['eta_days']=max(0,row['target_ns']-row['production_ns'])/row['ns_per_day']
        row['checkpoint_present']=(folder/'production.cpt').is_file()
        row['checkpoint_validation']='validated by backend on resume; existence alone is insufficient'
    return state

def handoff(rid,reviews):
    from webui import workflow as w
    d=w.run_dir(rid);cfg=w.read(d/'config.json');rc=cfg['research'];w.require(d,'gromacs')
    records=w.read(d/'gromacs/systems.json',[])
    qualified={r['compound_id']:r for r in map(json.loads,(d/'docking/qualification/qualification.jsonl').read_text().splitlines())}
    ready=[];pending=[]
    for rec in records:
        cid=rec['compound_id'];m=reviews.get(cid,{})
        q=qualified.get(cid,{})
        for name,sha in rec['files'].items():
            if w.digest(d/name)!=sha:raise ValueError('MD artifact changed before handoff')
        a=rec['analysis'];checks={
            'identity_clear':q.get('stereochemistry')=='specified_or_not_applicable',
            'chemistry_passed':bool(q.get('smiles')),
            'risk_screen_passed':q.get('risk',{}).get('status')=='risk_screen_passed',
            'acquisition_confirmed':bool(m.get('acquisition_confirmed') and m.get('acquisition_evidence')),
            'structural_reviewed':m.get('structural_reviewed') is True,
            'md_completed':rec.get('md_completed') is True and a.get('production_end_ns',0)>=rc['md_production_ns'],
            'md_analysis_completed':rec.get('md_analysis_completed') is True,
            'manual_review_approved':m.get('approved') is True and bool(m.get('reviewer')) and bool(m.get('date')),
            'uncertainties_recorded':bool(m.get('uncertainties'))}
        if checks['acquisition_confirmed']:
            from .qualification import verified_file
            verified_file(w.ROOT,m['acquisition_evidence'])
        row={'compound_id':cid,'qualification':q,'md':rec,'manual_review':m,'checks':checks,'wetlab_ready':all(checks.values())}
        (ready if row['wetlab_ready'] else pending).append(row)
    if len(ready)>rc['wetlab_max']:raise ValueError('Approved candidates exceed wet-lab cap; submit an explicit selection')
    stamp=str(time.time_ns());out=d/'wetlab'/stamp;out.mkdir(parents=True)
    w.atomic(out/'handoff.json',{'ready':ready,'pending':pending,'target':[rc['wetlab_min'],rc['wetlab_max']],
                              'shortfall':max(0,rc['wetlab_min']-len(ready)),'evidence_level':'computational handoff; not experimentally validated'})
    # Trajectories remain in candidate directories; manifest refers to hash-verified originals.
    with zipfile.ZipFile(out/'handoff.zip','w',zipfile.ZIP_DEFLATED) as archive:
        archive.write(out/'handoff.json','handoff.json');archive.write(d/'resolved_config.json','resolved_config.json')
        for row in ready:
            sdf=d/row['qualification']['sdf'] if row['qualification'].get('sdf') else None
            if sdf and sdf.is_file():archive.write(sdf,row['compound_id']+'/structure.sdf')
    return {'artifact':str((out/'handoff.zip').relative_to(d)),'wetlab_ready':len(ready),'pending':len(pending),
            'trajectory_delivery':'Original trajectory paths and SHA256 manifest included; copy referenced files for external transfer'}
