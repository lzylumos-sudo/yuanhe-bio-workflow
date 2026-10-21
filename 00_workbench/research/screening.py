"""Evaluate Y-MAVS against each actual receptor, then retain candidate evidence."""

def run(adapter,d,cfg,out):
    from webui.workflow import read,atomic
    from pipeline_bridge import csv_write
    records=read(d/'docking/results.json',[]);receptors={}
    for row in records:
        if row['status']=='ok':receptors[row['receptor_sha256']]=row['receptor_file']
    rows=[]
    for sha,path in receptors.items():
        child=out/sha[:16];child.mkdir()
        adapter(d,{**cfg,'_screen_receptor':{'sha256':sha,'path':path}},child)
        rows.extend(read(child/'candidates.json',[]))
    # Existing Y-MAVS qualification stays authoritative. Keep evidence from each conformation.
    rows.sort(key=lambda r:(r.get('qualification')!='formal',-float(r.get('evidence_score') or 0),r.get('compound_id','')))
    chosen=set();selected=[]
    for row in rows:
        cid=row.get('research_compound_id',row.get('compound_id'))
        if cid not in chosen and len(chosen)>=cfg['research']['review_limit']:continue
        chosen.add(cid);selected.append(row)
    atomic(out/'all_conformation_candidates.json',rows);atomic(out/'candidates.json',selected)
    csv_write(out/'molecules.csv',selected,list(dict.fromkeys(k for r in selected for k in r)))
    return {'input':len({r['compound_id'] for r in records}),'remaining':len(chosen),
            'formal':len({r['compound_id'] for r in selected if r.get('qualification')=='formal'}),'receptor_conformations':len(receptors)}
