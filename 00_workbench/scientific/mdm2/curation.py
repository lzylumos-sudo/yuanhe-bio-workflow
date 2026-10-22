"""Primary-table-verified curation for the declared retrospective test library."""
from pathlib import Path
from .library import all_activities, reference_smiles, canonical_parent
from scientific import artifacts as a

def curate(root):
    root=Path(root)
    reviewed={'PMC4281096.bioc.json':'727707d969d32955b11054829bcddb72251aecfcd22856ba43f4deed7cbd88fa',
              'SAR405838_primary.pdf':'87ca5ee4961e36917d3e70cb9a699011dc6f81fad0cf52e355d38f548e1d66c8',
              'Brigimadlin_tableS4.pdf':'29a70f95720d83bdd9a2305fe586e139b1204eca855876f6f34cf004c97e44b0',
              'Brigimadlin_methods.pdf':'24a586c055822973ffa78c96bc4dcce1e00c36c3a492ada1eea6cdb006356cee'}
    for filename,digest in reviewed.items():
        if a.digest(root/'evaluation/literature'/filename)!=digest:
            raise ValueError('Primary source changed; explicit review required: '+filename)
    act=all_activities(root);refs=reference_smiles(root)
    entries=[(3355410,'12D',19),(3355411,'13D',3.4),(3355412,'14D',105),(3355414,'15D',19),(3355415,'16D',3.7),(3355416,'17D',2.9),(3355417,'18D',2.8),(3355418,'21D',22),(3355420,'22D',80),(3355421,'23D',1.7),(3355422,'24D',23),(3355424,'25D',67),(3355425,'26D',3),(3355426,'27D',15),(3355428,'28D',53),(3355429,'29D',1.7),(3355430,'31 (MI-1061)',.16)]
    rows=[]
    for cid,name,value in entries:
     hits=[r for r in act if r['molecule_chembl_id']==f'CHEMBL{cid}' and r['assay_chembl_id']=='CHEMBL3380268' and r['standard_type']=='Ki' and r['standard_units']=='nM']
     assert len(hits)==1,(cid,len(hits));h=hits[0];assert float(h['standard_value'])==value,(cid,h['standard_value'],value)
     record=a.read(root/f"evaluation/records/{h['record_id']}.json")
     assert name.split()[0] in record['compound_key'],record
     rows.append({'name':name,'source_id':f'CHEMBL{cid}','smiles':canonical_parent(h['canonical_smiles']),
     'primary_source_verified':True,'activity_type':'Ki','relation':'=','activity_nM':value,'primary_table':'Table 6',
     'doi':'10.1021/jm501541j','primary_url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4281096/',
     'primary_file':'evaluation/literature/PMC4281096.bioc.json','primary_sha256':a.digest(root/'evaluation/literature/PMC4281096.bioc.json'),
     'assay_id':h['assay_chembl_id'],'activity_id':h['activity_id'],'record_id':h['record_id'],
     'record_compound_name':record['compound_name'],'target':'human MDM2 (1-118), His-tagged',
     'assay':'competitive fluorescence polarization; FAM-p53 tracer; biochemical Ki, not cell viability',
     'conditions':'pH 7.5; 100 mM potassium phosphate; 100 ug/ml bovine gamma globulin; .02% sodium azide; .01% Triton X-100; 4% DMSO; 10 nM MDM2; 1 nM tracer; 30 min RT',
     'series':'2014 stable spirooxindole D-series','stereo_source':'ChEMBL stereochemical structure + primary compound number + original compound record',
     'inclusion_rule':'human MDM2 direct biochemical Ki or p53 competition IC50 <=1000 nM; fully specified source stereo; stable D series; no cell-only evidence'})
    for name,cid,value in [('Nutlin-3a','CHEMBL191334',43.6),('SAR405838','CHEMBL2381408',.88),('Brigimadlin','CHEMBL5314455',2)]:
     r=refs[name];file='evaluation/literature/'+('Brigimadlin_tableS4.pdf' if name=='Brigimadlin' else 'SAR405838_primary.pdf')
     rows.append({'name':name,'source_id':cid,**r,'primary_source_verified':True,
     'activity_type':'IC50' if name=='Brigimadlin' else 'Ki','relation':'=','activity_nM':value,
     'primary_table':'Supplementary Table S4' if name=='Brigimadlin' else 'Figure 1b',
     'doi':'10.1158/1535-7163.MCT-23-0783' if name=='Brigimadlin' else '10.1158/0008-5472.CAN-14-0799',
     'primary_url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11612618/' if name=='Brigimadlin' else 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4247201/',
     'primary_file':file,'primary_sha256':a.digest(root/file),
     'target':'human GST-MDM2 (23-117)' if name=='Brigimadlin' else 'human MDM2 (1-118)',
     'assay':'AlphaScreen MDM2-p53 biochemical competition' if name=='Brigimadlin' else 'competitive fluorescence polarization',
     'conditions':'Supplementary Methods p2: 5 ul 10 nM GST-MDM2/p53 mix +5 ul compound then5 ul beads; pH7.2 Tris50mM NaCl120mM BSA.1% DTT5mM EDTA1mM Tween20.01%; 5% final DMSO;15+30minRT; >=5 duplicates' if name=='Brigimadlin' else 'primary Figure 1b + Methods + MDM2 1-118; not cell viability',
     'series':'reference_'+name,'inclusion_rule':'human MDM2 biochemical Ki or p53 competition IC50 <=1000 nM; verified stereochemistry',
     'stereo_source':r.get('stereo_source','PDB CCD stereochemical descriptor matched to ChEMBL structure')})
     assert canonical_parent(next(h['canonical_smiles'] for h in act if h['molecule_chembl_id']==cid))==r['smiles'],name
    assert len(rows)==20 and len({r['smiles'] for r in rows})==20
    a.atomic(root/'evaluation/curated_actives.json',rows)
    a.atomic(root/'evaluation/curation_protocol.json',{'inclusion':'human MDM2 direct Ki/Kd or biochemical p53 competition IC50 <=1000 nM; literature stable stereoisomer; primary table and identifier verification','excluded':'cell-only measurements; unstable C diastereomers; missing identity; unknown activity','active_count':20,'library_target':550,'series_bias':'17/20 from one medicinal chemistry series; limited retrospective validation, not broad scaffold generalization','references':refs})
    return {'verified_actives':len(rows)}
