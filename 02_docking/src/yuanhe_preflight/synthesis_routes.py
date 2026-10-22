"""Bounded local retrosynthesis with explicit planning-pool and chemistry checks."""
from __future__ import annotations
import os, subprocess, sys, tempfile, time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from filelock import FileLock
from .common import canonical, file_hash, hash_obj, read, write, seal

ROUTE_LABEL='Synthesizable | computed route, chemistry confirmation pending'
REVIEW_LABEL='Synthesizable | route reviewed by a person'

def assets(retro):
    result={}
    for key, status in (('model_path','model_missing'),('templates_path','templates_missing'),('stock_path','stock_missing')):
        p=Path(retro.get(key) or '')
        if not p.is_file(): return {},status
        result[key]={'path':str(p.resolve()),'sha256':file_hash(p)}
    for key in ('model_version','stock_version','templates_version','model_license','stock_license','templates_license'):
        if not retro.get(key): return {},'asset_metadata_missing:'+key
        result[key]=retro[key]
    return result,None

def validate_route(tree, target, *, route_review=None):
    """Check complete alternating reaction tree and planning-pool termination.
    Replay local retro templates, and conservatively review stereogenic cases."""
    problems=[]; leaves=[]; steps=[]; stereo=False; selectivity=False
    if tree.get('route_metadata',{}).get('is_solved') is not True:
        problems.append('planner_did_not_solve_route')
    if tree.get('type')!='mol' or canonical(tree.get('smiles',''))!=canonical(target):
        return {'status':'identity_mismatch','problems':['route_root_differs_from_target']}
    def walk(node, ancestry, depth=0):
        nonlocal stereo, selectivity
        if depth>100: raise ValueError('route_depth_exceeded')
        if node.get('type')!='mol': raise ValueError('expected_molecule_node')
        smi=canonical(node['smiles']); mol=Chem.MolFromSmiles(smi)
        if smi in ancestry: raise ValueError('cyclic_route')
        if Chem.FindPotentialStereo(mol): stereo=True
        children=node.get('children',[])
        if not children:
            in_pool=node.get('in_stock') is True
            leaves.append({'smiles':smi,'in_planning_pool':in_pool})
            if not in_pool: problems.append('terminal_not_in_planning_pool:'+smi)
            return
        if len(children)!=1 or children[0].get('type')!='reaction': raise ValueError('incomplete_or_alternative_route')
        rxn=children[0]; precursors=rxn.get('children',[]); meta=rxn.get('metadata') or {}
        if not precursors: raise ValueError('reaction_missing_precursors')
        template=meta.get('retro_template') or meta.get('template') or meta.get('smarts')
        if not template: problems.append('reaction_template_missing')
        else:
            try:
                reaction=rdChemReactions.ReactionFromSmarts(template)
                if reaction is None or reaction.GetNumReactantTemplates()!=1: raise ValueError('invalid_retro_template')
                predicted={tuple(sorted(canonical(x) for x in products)) for products in reaction.RunReactants((mol,),maxProducts=1000)}
                expected=tuple(sorted(canonical(x['smiles']) for x in precursors))
                if expected not in predicted: problems.append('template_identity_or_stereo_mismatch')
                if len(predicted)>1:selectivity=True
            except Exception as exc: problems.append('template_validation_failed:'+str(exc))
        if meta.get('unresolved_selectivity') or meta.get('missing_steps'): problems.append('backend_unresolved_chemistry')
        steps.append({'product':smi,'precursors':[canonical(x['smiles']) for x in precursors],
                      'reaction_smiles':rxn.get('smiles'), 'template':template,'metadata':meta,'depth':depth+1})
        for child in precursors: walk(child,ancestry|{smi},depth+1)
    try:walk(tree,set())
    except (ValueError,KeyError,TypeError) as exc: problems.append(str(exc))
    if not steps: problems.append('no_reaction_steps')
    route_hash=hash_obj(tree)
    review_valid=bool(route_review and route_review.get('route_sha256')==route_hash and
                      route_review.get('target_smiles')==canonical(target) and
                      all(str(route_review.get(k,'')).strip() for k in ('reviewer','reason','reviewed_at')) and
                      route_review.get('identity_stereo_selectivity_confirmed') is True)
    if stereo and not review_valid: problems.append('stereochemistry_selectivity_requires_review')
    if selectivity and not review_valid: problems.append('template_selectivity_requires_review')
    # Human review cannot replace absent steps, templates, or planning-pool termination.
    return {'status':'complete' if not problems else 'needs_review', 'problems':sorted(set(problems)),
            'route_sha256':route_hash,'steps':steps,'step_count':len(steps),'max_depth':max([x['depth'] for x in steps] or [0]),
            'terminal_materials':leaves,'review':route_review if review_valid else None, 'tree':tree,
            'interpretation':'predicted route, not experimental synthesis; no yield or conditions inferred'}

class SynthesisFilter:
    def __init__(self, config, rules_sha256):
        self.cfg=config; self.rules_sha256=rules_sha256
        self.retro=config.get('retro',{}); self.asset_manifest,self.asset_error=assets(self.retro)
        self.cache=Path(config['cache_dir']).resolve();self.cache.mkdir(parents=True,exist_ok=True)

    def search(self, smiles):
        if self.asset_error:return {'status':self.asset_error, 'routes':[]}
        python=self.retro.get('python') or sys.executable
        key=hash_obj({'schema':'retro-cache-v2-synthesis','structure':canonical(smiles),'assets':self.asset_manifest,
                      'rules':self.rules_sha256,
                      'retro_settings':self.retro,'worker_sha256':file_hash(Path(__file__).with_name('aizynth_worker.py'))})
        folder=self.cache/key;folder.mkdir(exist_ok=True)
        # Cross-task coordination and atomically completed per-molecule entries.
        with FileLock(str(folder/'search.lock'),timeout=self.retro['time_limit_seconds']+120):
            cached=folder/'result.json'
            if cached.is_file():
                record=read(cached)
                if record.get('cache_key')==key and record.get('raw_sha256')==hash_obj(record.get('raw')) and record['raw'].get('status')=='searched':
                    return {**record['raw'],'cache_key':key,'cache_hit':True,'search_artifact':str(cached)}
            request={'smiles':canonical(smiles),'retro':self.retro}
            write(folder/'request.json',request)
            worker=Path(__file__).with_name('aizynth_worker.py');output=folder/('raw_'+str(time.time_ns())+'.json')
            started=time.time()
            try:
                with (folder/'stdout.log').open('w') as out,(folder/'stderr.log').open('w') as err:
                    subprocess.run([python,str(worker),'--request',str(folder/'request.json'),'--output',str(output)],
                                   stdout=out,stderr=err,check=True,timeout=self.retro['time_limit_seconds']+self.retro.get('startup_grace_seconds',60))
                raw=read(output)
            except subprocess.TimeoutExpired:raw={'status':'search_timeout','routes':[]}
            except FileNotFoundError as exc:raw={'status':'dependency_missing','error':str(exc),'routes':[]}
            except (subprocess.CalledProcessError,ValueError) as exc:raw={'status':'backend_error','error':str(exc),'routes':[]}
            record={'cache_key':key,'raw':raw,'raw_sha256':hash_obj(raw),'assets':self.asset_manifest,'elapsed_seconds':time.time()-started}
            write(cached,record)
            return {**raw,'cache_key':key,'cache_hit':False,'search_artifact':str(cached)}

    def screen(self, smiles, *, allow_search=True):
        smiles=canonical(smiles); routes=[]; raw=None
        # Reuse only hash-bound searches against the same scientific assets.
        for review in self.cfg.get('route_reviews',[]):
            artifact=Path(review.get('search_artifact',''))
            if artifact.is_file():
                old=read(artifact)
                if (not self.asset_error and old.get('raw_sha256')==hash_obj(old.get('raw'))
                    and old.get('assets')==self.asset_manifest
                    and old.get('raw',{}).get('target_smiles')==smiles):
                    raw={**old['raw'],'search_artifact':str(artifact),'cache_hit':True};break
        if raw is None:
            raw=self.search(smiles) if allow_search else {'status':'not_evaluated_budget','routes':[]}
        status=raw.get('status','backend_error')
        if status=='searched':
            for tree in raw.get('routes',[]):
                review=next((x for x in self.cfg.get('route_reviews',[]) if x.get('route_sha256')==hash_obj(tree)),None)
                routes.append(validate_route(tree,smiles,route_review=review))
            status='route_supported' if any(r['status']=='complete' for r in routes) else 'route_needs_review' if routes else 'no_route_found'
        labels=[REVIEW_LABEL if r['review'] else ROUTE_LABEL for r in routes if r['status']=='complete']
        solved=any(r.get('step_count',0)>0 and r.get('tree',{}).get('route_metadata',{}).get('is_solved') is True for r in routes)
        return seal({'route_status':'ROUTE_FOUND' if solved else 'NO_ROUTE_FOUND' if raw.get('status')=='searched' or status=='no_route_found' else 'NOT_EVALUATED' if status=='not_evaluated_budget' else 'ERROR',
                     'synthesis_label':'SYNTHESIZABLE_PREDICTED' if labels else 'SYNTHESIS_REVIEW_REQUIRED',
                     'status_labels':[] if labels else [status], 'smiles':smiles,
                     'decision':'pass' if labels else 'pending','status':status,'labels':sorted(set(labels)),
                     'routes':routes,'raw_search':raw,'asset_manifest':self.asset_manifest,
                     'network':{'enabled':False,'sent_structures':[]},
                     'interpretation':'Computational route to a versioned planning building-block pool; not experimental synthesis.'})
