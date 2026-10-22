"""Supply first, then bounded local retrosynthesis; never substitutes SA scores."""
from __future__ import annotations
import os, subprocess, sys, tempfile, time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from filelock import FileLock
from .common import canonical, file_hash, hash_obj, read, write, seal
from .catalog import Catalog

BUY_LABELS={'in_stock':'Purchasable | stock evidence','made_to_order':'Purchasable | made-to-order evidence'}
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

def validate_route(tree, target, stock, *, route_review=None):
    """Check complete alternating reaction tree and every exact terminal supply record.
    Replay local retro templates, and conservatively review stereogenic cases."""
    problems=[]; leaves=[]; steps=[]; stereo=False; selectivity=False
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
            records=stock.valid(smi)
            leaves.append({'smiles':smi,'supplier_records':records})
            if not records: problems.append('terminal_supply_unverified:'+smi)
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
    # Human review cannot replace absent steps, templates, or material supply evidence.
    return {'status':'complete' if not problems else 'needs_review', 'problems':sorted(set(problems)),
            'route_sha256':route_hash,'steps':steps,'step_count':len(steps),'max_depth':max([x['depth'] for x in steps] or [0]),
            'terminal_materials':leaves,'review':route_review if review_valid else None, 'tree':tree,
            'interpretation':'predicted route, not experimental synthesis; no yield or conditions inferred'}

class AvailabilityFilter:
    def __init__(self, config, rules_sha256):
        self.cfg=config; self.rules_sha256=rules_sha256
        self.catalog=Catalog(config.get('catalogs',[]),config.get('evidence_max_age_days',30))
        self.stock=Catalog(config.get('stock_catalogs',[]),config.get('evidence_max_age_days',30))
        self.retro=config.get('retro',{}); self.asset_manifest,self.asset_error=assets(self.retro)
        self.cache=Path(config['cache_dir']).resolve();self.cache.mkdir(parents=True,exist_ok=True)

    def search(self, smiles):
        if self.asset_error:return {'status':self.asset_error, 'routes':[]}
        python=self.retro.get('python') or sys.executable
        key=hash_obj({'schema':'retro-cache-v1','structure':canonical(smiles),'assets':self.asset_manifest,
                      'rules':self.rules_sha256,'catalogs':self.catalog.versions,'stock':self.stock.versions,
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
        smiles=canonical(smiles); supply=self.catalog.match(smiles)
        external=None
        if self.cfg.get('engine',{}).get('pubchem',{}).get('enabled'):
            from .pubchem import PubChem
            external=PubChem(self.cfg['engine']['pubchem'],self.cfg['cache_dir']).query(smiles)
            supply+=external['records']
        valid=[x for x in supply if x['status'] in BUY_LABELS]
        labels=sorted({BUY_LABELS[x['status']] for x in valid}); routes=[]; raw=None
        status='purchase_supported' if valid else 'acquisition_unknown'
        # Reuse an explicitly linked existing search for additional labels. Buying
        # never triggers a new retrosynthesis search.
        if valid:
            for review in self.cfg.get('route_reviews',[]):
                artifact=Path(review.get('search_artifact',''))
                if artifact.is_file():
                    old=read(artifact)
                    if old.get('raw_sha256')==hash_obj(old.get('raw')) and old.get('assets')==self.asset_manifest and old.get('raw',{}).get('target_smiles')==smiles:
                        raw={**old['raw'],'search_artifact':str(artifact),'cache_hit':True};break
        else:raw=self.search(smiles) if allow_search else {'status':'not_evaluated_budget','routes':[]}
        if raw and raw.get('status')=='searched':
            for tree in raw.get('routes',[]):
                review=next((x for x in self.cfg.get('route_reviews',[]) if x.get('route_sha256')==hash_obj(tree)),None)
                routes.append(validate_route(tree,smiles,self.stock,route_review=review))
            for r in routes:
                if r['status']=='complete':labels.append(REVIEW_LABEL if r['review'] else ROUTE_LABEL)
            if not valid:
                status='route_supported' if labels else 'route_needs_review' if routes else ('search_timeout' if float(raw.get('statistics',{}).get('search_time',0))>=self.retro['time_limit_seconds'] else 'no_route_found')
        elif raw and not valid:status=raw.get('status','backend_error')
        pending_labels=[]
        if any(r['status']=='catalog_hit_supply_unconfirmed' for r in supply):pending_labels.append('Catalog match; supply confirmation pending')
        if any(r['status'] in {'supply_evidence_expired','catalog_expired'} for r in supply):pending_labels.append('Supply evidence expired')
        if status=='search_timeout':pending_labels.append('Search timed out')
        if status=='route_needs_review':pending_labels.append('Incomplete route or key chemistry requires review')
        if status=='no_route_found':pending_labels.append('No route found (not evidence of infeasibility)')
        if not labels:pending_labels.append('Acquisition pathway unknown')
        return seal({'route_status':'ROUTE_FOUND' if routes else 'NO_ROUTE_FOUND' if status=='no_route_found' else 'NOT_EVALUATED' if status in ('purchase_supported','not_evaluated_budget') else 'ERROR',
                     'availability_label':'PURCHASABLE_VERIFIED' if valid else 'SYNTHESIZABLE_PREDICTED' if status=='route_supported' else 'PURCHASABLE_UNCONFIRMED' if supply else 'AVAILABILITY_UNKNOWN',
                     'status_labels':pending_labels,'smiles':smiles,'decision':'pass' if labels else 'pending','status':status,'labels':sorted(set(labels)),
                     'external_query':external,'supplier_records':supply,'routes':routes,'raw_search':raw,'asset_manifest':self.asset_manifest,
                     'catalog_versions':self.catalog.versions,'stock_versions':self.stock.versions,
                     'evidence_max_age_days':self.catalog.max_age_days,'network':{'enabled':external is not None,'sent_structures':[smiles] if external else []},
                     'final_acquisition_review':'required_before_wet_lab_or_order'})
