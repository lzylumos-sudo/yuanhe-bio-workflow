"""Auditable chemical identities and evidence gates. Never infer availability from scores."""
import csv, hashlib, json, sqlite3
from pathlib import Path
from .config import fingerprint

def identity(smiles):
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize as ms
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError('Invalid molecular structure')
    mol=ms.Cleanup(mol)
    mol=ms.FragmentParent(mol)
    mol=ms.Uncharger().uncharge(mol)
    mol=ms.TautomerEnumerator().Canonicalize(mol)
    for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    canon=Chem.MolToSmiles(mol,isomericSmiles=True)
    stereo='unspecified' if any(v=='?' for _,v in Chem.FindMolChiralCenters(mol,includeUnassigned=True)) else 'specified_or_not_applicable'
    return 'CMP'+hashlib.sha256(canon.encode()).hexdigest()[:16].upper(),canon,stereo

def verified_file(root,record):
    path=(root/record.get('path','')).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():raise ValueError('Evidence must be an existing project-local file')
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    if record.get('sha256')!=sha:raise ValueError('Evidence hash mismatch')
    return path

def risk(mol,rules,review):
    from rdkit import Chem
    if not review or not review.get('reviewer') or not review.get('date') or review.get('rules_sha256')!=fingerprint(rules):
        raise ValueError('Risk rules require a versioned target-specific review and matching hash')
    if not rules:raise ValueError('No validated risk rules provided; screening cannot silently pass')
    hits=[]
    for rule in rules:
        if not all(rule.get(k) for k in ('id','version','smarts','rationale','source','action')):raise ValueError('Incomplete risk rule')
        if rule['action'] not in ('block','review'):raise ValueError('Invalid risk rule action')
        query=Chem.MolFromSmarts(rule['smarts'])
        if query is None:raise ValueError('Invalid SMARTS')
        if rule.get('category')=='PAINS' and rule['action']=='block':raise ValueError('PAINS alone is not a toxicity blocking rule')
        matches=mol.GetSubstructMatches(query)
        if matches:hits.append({**rule,'matched_atoms':matches})
    status='risk_blocked' if any(x['action']=='block' for x in hits) else 'risk_flagged_review' if hits else 'risk_screen_passed'
    return {'status':status,'hits':hits,'rules_sha256':fingerprint(rules),'interpretation':'Rule screen only; not evidence of non-toxicity'}

def acquisition(root,smiles,records):
    cid,_,_=identity(smiles);accepted=[]
    for rec in records:
        try:
            from rdkit import Chem
            exact=lambda s:Chem.MolToSmiles(Chem.MolFromSmiles(s),isomericSmiles=True)
            if exact(rec.get('smiles',''))!=exact(smiles):continue
            status=rec.get('status','acquisition_pending')
            if status not in ('purchasable_verified','synthesizable_predicted','synthesizable_reviewed'):continue
            if not rec.get('date'):continue
            verified_file(root,rec['evidence'])
            if status=='purchasable_verified':
                if not all(rec.get(k) for k in ('supplier','catalog_id','availability','verified_by')):continue
                if rec['availability'] not in ('in_stock','made_to_order'):continue
            else:
                route=rec.get('route',{})
                if not route.get('steps') or not route.get('starting_materials') or not route.get('uncertainties'):continue
                if identity(route.get('product_smiles',''))[0]!=cid:continue
                if not all(x.get('reaction') and x.get('reactants') and x.get('products') for x in route['steps']):continue
                for material in route['starting_materials']:
                    if not material.get('smiles') or not material.get('source'):raise ValueError('Starting material provenance absent')
                    verified_file(root,material['evidence'])
                # A qualified route requires an explicit technical check, even when its status remains predicted.
                if not route.get('technical_reviewer') or not route.get('connectivity_checked'):continue
                if status=='synthesizable_reviewed' and not rec.get('reviewer'):continue
            accepted.append(rec)
        except (ValueError,KeyError,TypeError):continue
    if not accepted:return {'status':'acquisition_pending','evidence':[]}
    rank={'synthesizable_predicted':0,'synthesizable_reviewed':1,'purchasable_verified':2}
    return {'status':max(accepted,key=lambda x:rank[x['status']])['status'],'evidence':accepted,
            'interpretation':'Predicted route is not experimentally demonstrated synthesis'}

def supplier_rows(root,path,limit):
    if not path:return
    source=(root/path).resolve()
    if not source.is_relative_to(root.resolve()) or not source.is_file():raise ValueError('Supplier library missing or outside project')
    from rdkit import Chem
    sha=hashlib.sha256(source.read_bytes()).hexdigest()
    if source.suffix.lower()=='.sdf':
        rows=({'smiles':Chem.MolToSmiles(m,isomericSmiles=True),**m.GetPropsAsDict()} for m in Chem.SDMolSupplier(str(source)) if m is not None)
        handle=None
    else:
        handle=source.open()
        rows=csv.DictReader(handle) if source.suffix.lower()=='.csv' else ({'smiles':line.split()[0]} for line in handle if line.strip())
    try:
        seen=set()
        for row in rows:
            cid,_,_=identity(row['smiles'])
            if cid not in seen and len(seen)>=limit:break
            seen.add(cid)
            yield {**row,'source_stage':'supplier','library_sha256':sha,'library_path':str(source.relative_to(root))}
    finally:
        if handle:handle.close()

def qualify(root,d,cfg,generated,out):
    from rdkit import Chem
    rc=cfg['research'];out.mkdir(parents=True,exist_ok=True)
    # Reuse the project's real pre-docking rule/catalog/retrosynthesis implementation.
    import sys
    sys.path.insert(0,str(root/'02_docking/src'))
    from yuanhe_preflight.toxicity_substructure_filter import RiskRules
    from yuanhe_preflight.availability_synthesis_filter import AvailabilityFilter
    from yuanhe_preflight.config import validate
    from yuanhe_preflight.common import file_hash
    pc=json.loads(json.dumps(rc['predocking']))
    for key in ('rules_path','cache_dir'):
        if not Path(pc[key]).is_absolute():pc[key]=str(root/pc[key])
    validate(pc)
    rules=RiskRules(pc['rules_path'],pc['risk_policies'])
    availability=AvailabilityFilter(pc,file_hash(pc['rules_path']))
    if rc['risk_rules']:risk(Chem.MolFromSmiles('C'),rc['risk_rules'],rc['risk_rules_review'])
    db=sqlite3.connect(out/'compounds.sqlite3')
    db.execute('CREATE TABLE IF NOT EXISTS compounds (compound_id TEXT PRIMARY KEY, smiles TEXT, stereo TEXT, row_json TEXT)')
    db.execute('CREATE TABLE IF NOT EXISTS sources (compound_id TEXT, source_hash TEXT, record TEXT, PRIMARY KEY(compound_id,source_hash))')
    counts={'valid':0,'unique':0,'chemistry_passed':0,'toxicity_structure_passed':0,'acquisition_screening_qualified':0,'screening_qualified':0}
    from itertools import chain
    try:
        with db:
            for row in chain(generated,supplier_rows(root,rc['supplier_library'],rc['external_limit'])):
                try:cid,canon,stereo=identity(row['smiles'])
                except ValueError:continue
                counts['valid']+=1
                db.execute('INSERT OR IGNORE INTO compounds VALUES (?,?,?,?)',(cid,canon,stereo,json.dumps(row)))
                db.execute('INSERT OR IGNORE INTO sources VALUES (?,?,?)',(cid,fingerprint(row),json.dumps(row)))
        qualified=[]
        with (out/'qualification.jsonl').open('w') as log:
            for cid,canon,stereo,raw in db.execute('SELECT * FROM compounds ORDER BY compound_id'):
                row=json.loads(raw);counts['unique']+=1;counts['chemistry_passed']+=1
                raw_risk=rules.screen(canon,[e for e in pc['risk_exceptions'] if e.get('structure_smiles')==canon])
                r={'status':{'pass':'risk_screen_passed','exclude':'risk_blocked'}.get(raw_risk['decision'],'risk_flagged_review'),'evidence':raw_risk}
                if rc['risk_rules']:
                    extra=risk(Chem.MolFromSmiles(canon),rc['risk_rules'],rc['risk_rules_review']);r['additional_rules']=extra
                    if extra['status']!='risk_screen_passed':r['status']=extra['status']
                a=acquisition(root,canon,rc['acquisition_evidence'])
                if a['status']=='acquisition_pending' and r['status']=='risk_screen_passed':
                    raw=availability.screen(canon)
                    a={'status':'acquisition_pending','evidence':raw}
                    if raw['decision']=='pass':
                        if any(x.startswith('可买') for x in raw.get('labels',[])):a['status']='purchasable_verified'
                        else:a['status']='synthesizable_reviewed' if raw.get('status')=='synthesizable_reviewed' else 'synthesizable_predicted'
                counts['toxicity_structure_passed']+=r['status']=='risk_screen_passed'
                counts['acquisition_screening_qualified']+=a['status']!='acquisition_pending'
                item={**row,'id':cid,'compound_id':cid,'smiles':canon,'stereochemistry':stereo,'risk':r,'acquisition':a,
                      'sources':[json.loads(x[0]) for x in db.execute('SELECT record FROM sources WHERE compound_id=?',(cid,))]}
                eligible=r['status']=='risk_screen_passed' and a['status']!='acquisition_pending'
                item['screening_qualified']=eligible
                log.write(json.dumps(item)+'\n')
                if eligible:qualified.append(item);counts['screening_qualified']+=1
        return qualified,counts
    finally:db.close()
