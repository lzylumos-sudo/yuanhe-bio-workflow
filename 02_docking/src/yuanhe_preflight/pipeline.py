"""Durable batch stages; all raw rows retained, even when no molecule qualifies."""
from __future__ import annotations
import csv, shutil, time, json, sqlite3
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from collections import Counter
from rdkit import Chem
from .common import identity, canonical, sdf_smiles, file_hash, hash_obj, seal, write, read, SCHEMA
from .config import validate
from .toxicity_substructure_filter import RiskRules
from .availability_synthesis_filter import AvailabilityFilter


def counts(rows, input_count=None):
    c=Counter(r.get('decision','pending') for r in rows)
    return {'input':input_count if input_count is not None else len(rows),'completed':len(rows),
            'passed':c['pass'],'remaining':c['pass'],'excluded':c['exclude'],
            'review':c['review']+c['pending'],'duplicates':c['duplicate'],
            'labels':dict(Counter(label for r in rows for label in r.get('labels',[]))),
            'failure_reasons':dict(Counter(r.get('reason') or r.get('status','') for r in rows if r.get('decision')!='pass'))}

def export(out,rows):
    write(out/'records.json',rows)
    with sqlite3.connect(out/'screening.sqlite3') as db:
        db.execute('CREATE TABLE IF NOT EXISTS molecules (input_index INTEGER PRIMARY KEY, molecule_id TEXT, decision TEXT, record TEXT)')
        db.execute('DELETE FROM molecules')
        db.executemany('INSERT INTO molecules VALUES (?,?,?,?)', [(i,r.get('molecule_id',r.get('compound_id')),r.get('decision'),json.dumps(r,allow_nan=False)) for i,r in enumerate(rows)])
    with (out/'passed.smi').open('w') as handle:
        for row in rows:
            if row.get('decision')=='pass':handle.write(row['smiles']+' '+row['compound_id']+'\n')
    fields=list(dict.fromkeys(k for r in rows for k in r)) or ['id','decision']
    with (out/'molecules.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        w.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)
    for name,decisions in [('passed',{'pass'}),('excluded',{'exclude','duplicate'}),('review',{'pending','review'})]:
        with Chem.SDWriter(str(out/(name+'.sdf'))) as writer:
            for r in rows:
                if r.get('decision') not in decisions:continue
                mol=Chem.MolFromSmiles(r.get('smiles',''))
                if mol:
                    for k in ('id','compound_id','state_id','decision','reason'):mol.SetProp(k,str(r.get(k,'')))
                    writer.write(mol)

def risk_stage(rows, config, out, *, base, batch_id):
    validate(config);out=Path(out);out.mkdir(parents=True,exist_ok=True);base=Path(base)
    rules=RiskRules(config['rules_path'],config['risk_policies']);write(out/'rules.snapshot.json',rules.data)
    write(out/'config.snapshot.json',config);write(out/'raw_input.json',rows)
    from .structure import screen, standardize
    started=time.monotonic()
    results=[];seen={}
    for index,row in enumerate(rows):
        r={**row,'original_smiles':row.get('original_smiles',row.get('smiles')),'generation_batch':row.get('generation_batch',batch_id),
           'alert_ids':[],'alert_types':[],'matched_atom_indices':{},'severity':'unknown','rule_version':rules.version,
           'source_core':row.get('source_core',row.get('core')),'original_id':row.get('id'), 'input_index':index}
        try:
            r['standardization']=standardize(row['smiles'])
            ident=identity(row['smiles']);r.update(ident);r['id']=ident['compound_id'];r['molecule_id']=ident['compound_id'];r['canonical_smiles']=ident['smiles']
            source=(base/(row.get('sdf') or row.get('raw_sdf',''))).resolve()
            if sdf_smiles(source)!=r['smiles']:raise ValueError('SMILES_SDF_identity_mismatch')
            # Immutable copy retains original atom order/coordinates and all source properties.
            local=out/'inputs'/f'{index:08d}.sdf';local.parent.mkdir(exist_ok=True);shutil.copyfile(source,local)
            original_mol=Chem.SDMolSupplier(str(local),removeHs=False)[0]
            canon_mol=Chem.MolFromSmiles(r['smiles'])
            mapping=original_mol.GetSubstructMatch(canon_mol,useChirality=True)
            r['canonical_to_original_sdf_atom_indices']=list(mapping)
            r.update(source_sdf=str(local.resolve()),source_sdf_sha256=file_hash(local),sdf=str(local.resolve()))
            risk=screen(r['smiles'],rules,[e for e in config['risk_exceptions'] if e.get('structure_smiles')==r['smiles']],config.get('engine'))
            r['original_sdf_matched_atom_indices']={key:[[mapping[i] for i in match if i<len(mapping)] for match in matches] for key,matches in risk['matched_atom_indices'].items()}
            r.update(risk=risk,**{k:risk[k] for k in ('decision','alert_ids','alert_types','matched_atom_indices','severity','reason','rule_version')})
            if r['smiles'] in seen:
                r.update(decision='duplicate',reason='exact_structure_duplicate',duplicate_of=seen[r['smiles']])
            else:seen[r['smiles']]=r['id']
            (out/(r['id']+'.svg')).write_text(rules.svg(r['smiles']))
        except Exception as exc:r.update(decision='review',reason=type(exc).__name__+': '+str(exc))
        r['id']=r.get('compound_id') or 'INVALID'+hash_obj({'batch':batch_id,'index':index,'smiles':r.get('original_smiles')})[:16];r.setdefault('molecule_id',r['id'])
        r.update(filter_status={'pass':'PASS','exclude':'BLOCK'}.get(r['decision'],'REVIEW'),matched_rules=r.get('risk',{}).get('alerts',[]),alert_type=r.get('alert_types',[]),alert_severity=r.get('severity'),rule_source=[a['source'] for a in r.get('risk',{}).get('alerts',[])],rejection_reason=r.get('reason') if r.get('decision')!='pass' else '')
        results.append(r)
        if len(results)%config['batch_size']==0 or len(results)==len(rows):write(out/'progress.json',counts(results,len(rows)))
    origins={}
    for r in results:
        if r.get('compound_id'):origins.setdefault(r['compound_id'],[]).append({k:r.get(k) for k in ('original_id','source_core','generation_batch','seed_id','provenance_manifest')})
    for r in results:
        if r.get('compound_id'):r['all_origins']=origins[r['compound_id']]
    export(out,results);summary=counts(results,len(rows));summary['structure_seconds']=time.monotonic()-started;write(out/'summary.json',summary)
    return summary

def availability_stage(rows,config,out):
    validate(config);out=Path(out);out.mkdir(parents=True,exist_ok=True);write(out/'config.snapshot.json',config)
    started=time.monotonic();timings={}
    if config.get('engine'):
        from .engine import evaluate,rank
        rows,timings=evaluate(rows,config,out)
        export(out,rows)
        if config['engine']['missing_policy']=='abort' and any(r.get('screening',{}).get('decision')=='review' and (any(t['prediction_status']!='SUCCESS' for t in r.get('toxicity_predictions',[])) or r.get('synthesis',{}).get('prediction_status')=='ERROR') for r in rows):
            write(out/'summary.json',{**counts(rows),'timings':timings})
            raise RuntimeError('Required prediction unavailable; see persisted molecular failure records')
    filt=AvailabilityFilter(config,file_hash(config['rules_path']));eligible=[r for r in rows if r['decision']=='pass'];results=[]
    if config.get('engine'):
        eligible=rank(eligible,config['engine'])
    retro_ids=set(r['compound_id'] for r in eligible if not filt.catalog.valid(r['smiles']))
    if config.get('engine'):
        retro_ids=set(r['compound_id'] for r in eligible if r['compound_id'] in retro_ids)
        retro_ids=set([r['compound_id'] for r in eligible if r['compound_id'] in retro_ids][:config['engine']['retro_max_candidates']])
    acquisition_started=time.monotonic()
    def process(row):
        from .common import verify_seal
        verify_seal(row['risk'])
        # Risk policies cannot be silently relaxed between stages.
        from .structure import screen
        actual=screen(row['smiles'],RiskRules(config['rules_path'],config['risk_policies']),row['risk']['exceptions'],config.get('engine'))
        if actual!=row['risk'] or actual['decision']!='pass':raise ValueError('risk_evidence_changed')
        availability=filt.screen(row['smiles'],allow_search=row['compound_id'] in retro_ids)
        archived={}
        # Preserve raw supply documents and raw searches in run exports; model
        # weights remain external, identified by immutable hashes.
        sources=[Path(x['path']) for x in availability['catalog_versions']+availability['stock_versions']]
        supplier_rows=list(availability['supplier_records'])
        for route in availability['routes']:
            for leaf in route.get('terminal_materials',[]):supplier_rows.extend(leaf['supplier_records'])
        sources += [Path(x['path']) for x in (availability.get('external_query') or {}).get('evidence',[])]
        sources += [Path(x['evidence_path']) for x in supplier_rows if x.get('evidence_path')]
        raw=availability.get('raw_search') or {}
        if raw.get('search_artifact'):sources.append(Path(raw['search_artifact']))
        for src in sources:
            if not src.is_file():continue
            sha=file_hash(src);dest=out/'source_evidence'/(sha+src.suffix)
            dest.parent.mkdir(exist_ok=True)
            # Multiple worker threads may preserve the same catalog; unique
            # temporary files avoid partial readers and retain a hash check.
            import uuid,os
            temp=dest.with_name(dest.name+'.'+uuid.uuid4().hex+'.tmp')
            shutil.copyfile(src,temp)
            if file_hash(temp)!=sha:raise ValueError('source_evidence_changed_during_snapshot')
            os.replace(temp,dest);archived[str(src)]={'artifact':str(dest.relative_to(out)),'sha256':sha}
            if str(src)==raw.get('search_artifact'):availability['raw_search']['search_artifact']=str(dest.resolve())
        availability=seal(availability)
        r={**row,'source_evidence_archive':archived,'availability':availability,'decision':availability['decision'],'labels':availability['labels'],
           'status_labels':availability['status_labels'],'status':availability['status'],'reason':availability['status']}
        write(out/'items'/(r['compound_id']+'.json'),r)
        return r
    with ThreadPoolExecutor(max_workers=config['workers']) as pool:
        for start in range(0,len(eligible),config['batch_size']):
            for r in pool.map(process,eligible[start:start+config['batch_size']]):
                results.append(r);write(out/'progress.json',counts(results,len(eligible)))
    timings['acquisition_seconds']=time.monotonic()-acquisition_started
    if config.get('engine'):
        ordered=rank([r for r in results if r['decision']=='pass'],config['engine'])
        keep={r['compound_id'] for r in ordered[:config['engine']['keep']]}
        for r in results:
            if r['decision']=='pass' and r['compound_id'] not in keep:r.update(decision='review',reason='retention_limit')
        results=ordered[:config['engine']['keep']]+[r for r in results if r['decision']!='pass']+[r for r in rows if r['decision']!='pass']
        for r in results:
            if 'toxicity_predictions' not in r:
                r['toxicity_predictions']=[{'molecule_id':r['molecule_id'],'model_name':'ADMET-AI','model_version':config['engine']['toxicity']['model_version'],'endpoint':ep,'raw_prediction':None,'prediction_status':'NOT_EVALUATED','risk_level':'UNKNOWN','error_message':'Not eligible after structure/deduplication gate'} for ep in config['engine']['toxicity']['endpoints']]
            r.setdefault('synthesis',{'molecule_id':r['molecule_id'],'sa_score':None,'sc_score':None,'screening_status':'NOT_EVALUATED'})
            r.setdefault('availability',{'status':'not_evaluated','availability_label':'AVAILABILITY_UNKNOWN','route_status':'NOT_EVALUATED','routes':[],'supplier_records':[]})
            r.update(final_pass=r['decision']=='pass',docking_status='QUALIFIED' if r['decision']=='pass' else 'NOT_ELIGIBLE',rejection_reason='' if r['decision']=='pass' else r.get('reason',''))
    selection=out/'selection.json'
    write(selection,{'qualified_ids':[r['compound_id'] for r in results if r['decision']=='pass'],'policy_sha256':hash_obj(config)})
    for r in results:
        if r['decision']=='pass':
            proof=seal({'schema':SCHEMA,'identity':identity(r['smiles']),'risk':r['risk'],'availability':r['availability'],
                        'source_sdf':r['source_sdf'],'source_sdf_sha256':r['source_sdf_sha256'], 'config':config,
                        **({'screening':r['screening'],'selection_path':str(selection.resolve()),'selection_sha256':file_hash(selection)} if config.get('engine') else {}),
                        'provenance':{k:r.get(k) for k in ('original_id','source_core','generation_batch','seed_id','provenance_manifest','all_origins')}})
            p=out/'evidence'/(r['compound_id']+'.json');write(p,proof)
            r['qualification_evidence']=str(p.resolve());r['qualification_sha256']=proof['record_sha256']
            ligand=out/'ligands'/(r['compound_id']+'.sdf');ligand.parent.mkdir(exist_ok=True)
            shutil.copyfile(r['source_sdf'],ligand);write(str(ligand)+'.qualification.json',proof)
            r['qualified_ligand']=str(ligand.resolve())
            # Gate replays the evidence contract before adding to the docking pool.
            from .gate import require_ligand
            require_ligand(r['source_sdf'],p)
        write(out/'items'/(r.get('molecule_id',r['id'])+'_'+str(r.get('input_index',0))+'.json'),r)
    export(out,results);summary=counts(results,len(rows));summary['shortfall']=len(rows)-summary['passed']
    summary.update(timings=timings,elapsed_seconds=time.monotonic()-started,pass_rate=summary['passed']/len(rows) if rows else 0)
    write(out/'summary.json',summary);write(out/'qualified.json',[r for r in results if r['decision']=='pass'])
    write(out/'identity_evidence.json',{r['compound_id']:{k:r.get(k) for k in ('compound_id','state_id','source_core','generation_batch','labels','all_origins','qualification_evidence','qualification_sha256')} for r in results if r.get('compound_id') and r.get('decision')!='duplicate'})
    return summary
