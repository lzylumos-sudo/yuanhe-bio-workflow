"""Auditable chemical identities and evidence gates. Never infer availability from scores."""
import csv, hashlib, json, sqlite3
from .config import fingerprint

def identity(smiles):
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize as ms
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError('Invalid molecular structure')
    mol=ms.Cleanup(mol)
    mol=ms.FragmentParent(mol)
    mol=ms.Uncharger().uncharge(mol)
    tautomers=ms.TautomerEnumerator()
    # Default tautomer enumeration removes stereochemistry at transformed centers.
    # Preserve specified tetrahedral and double-bond stereo in parent identities.
    tautomers.SetRemoveSp3Stereo(False)
    tautomers.SetRemoveBondStereo(False)
    mol=tautomers.Canonicalize(mol)
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
    from webui.workflow import read,require
    from yuanhe_preflight.gate import require_ligand
    require(d,'risk');require(d,'availability')
    out.mkdir(parents=True,exist_ok=True)
    risk_rows=read(d/'risk/records.json',[]);available=read(d/'availability/records.json',[])
    parents={r['compound_id']:r for r in read(d/'risk/parent_index.json',[])}
    counts={'valid':len(parents),'unique':len(parents),'chemistry_passed':len(parents),
            'toxicity_structure_passed':sum(r['decision']=='pass' for r in risk_rows),
            'acquisition_screening_qualified':0,'screening_qualified':0}
    rows=[]
    db=sqlite3.connect(out/'compounds.sqlite3')
    db.execute('CREATE TABLE compounds (compound_id TEXT PRIMARY KEY, smiles TEXT, record TEXT)')
    with db, (out/'qualification.jsonl').open('w') as journal:
        for record in available:
            cid=record['compound_id'];parent=parents.get(cid,{})
            a=record['availability'];passed=record['decision']=='pass'
            status='acquisition_pending'
            if passed:
                require_ligand(record['source_sdf'],record['qualification_evidence'])
                if a.get('status') == 'purchase_supported':status='purchasable_verified'
                elif any(x.get('review') for x in a.get('routes',[])):status='synthesizable_reviewed'
                else:status='synthesizable_predicted'
            item={**parent,**record,'sdf':record['source_sdf'],
                  'risk':{'status':'risk_screen_passed','evidence':record['risk']},
                  'acquisition':{'status':status,'evidence':a},'screening_qualified':passed}
            journal.write(json.dumps(item)+'\n')
            db.execute('INSERT INTO compounds VALUES (?,?,?)',(cid,item['smiles'],json.dumps(item)))
            if passed:rows.append(item);counts['acquisition_screening_qualified']+=1;counts['screening_qualified']+=1
    db.close()
    return rows,counts


def prepare_risk(root,d,cfg,out):
    """Normalize parent identities before the shared, visible pre-docking stages."""
    from itertools import chain
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from webui.workflow import read,atomic
    from yuanhe_preflight.pipeline import risk_stage
    rc=cfg['research'];rows=[];seen={}
    inputs=out/'standardized';inputs.mkdir(parents=True,exist_ok=True)
    for source in chain(read(d/'reinvent4/candidates.json',[]),supplier_rows(root,rc['supplier_library'],rc['external_limit'])):
        try:cid,canon,stereo=identity(source['smiles'])
        except ValueError:continue
        if cid in seen:
            seen[cid]['sources'].append(source);continue
        mol=Chem.AddHs(Chem.MolFromSmiles(canon));params=AllChem.ETKDGv3();params.randomSeed=cfg['seed']
        if AllChem.EmbedMolecule(mol,params)!=0:continue
        target=inputs/(cid+'.sdf')
        with Chem.SDWriter(str(target)) as writer:writer.write(mol)
        row={**source,'id':cid,'compound_id':cid,'smiles':canon,'stereochemistry':stereo,'sdf':str(target.relative_to(d)),'sources':[source]}
        seen[cid]=row;rows.append(row)
    atomic(out/'parent_index.json',rows)
    result=risk_stage(rows,cfg['predocking'],out,base=d,batch_id=d.name)
    if rc['risk_rules']:
        records=read(out/'records.json',[])
        for record in records:
            extra=risk(Chem.MolFromSmiles(record['smiles']),rc['risk_rules'],rc['risk_rules_review'])
            record['additional_rules']=extra
            if extra['status']!='risk_screen_passed':record['decision']='exclude' if extra['status']=='risk_blocked' else 'review'
        from yuanhe_preflight.pipeline import export,counts
        export(out,records);result=counts(records);atomic(out/'summary.json',result)
    return result
