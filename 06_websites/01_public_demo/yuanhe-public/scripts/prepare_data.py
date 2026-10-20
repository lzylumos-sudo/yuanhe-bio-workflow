"""Reproducible geometry from the public PDB snapshot only; Python standard library."""
from pathlib import Path
import json, hashlib, math
ROOT=Path(__file__).resolve().parents[1]
pdb=ROOT/'public/data/9NU4.pdb'
records=[]
for line in pdb.read_text().splitlines():
    if line.startswith('ENDMDL'): break
    if not line.startswith('ATOM  ') or line[16] not in (' ','A'): continue
    if line[76:78].strip() in ('H','D'): continue
    records.append(dict(chain=line[21],res=int(line[22:26]),icode=line[26].strip(),name=line[17:20].strip(),atom=line[12:16].strip(),xyz=[float(line[x:x+8]) for x in (30,38,46)]))
a=[r for r in records if r['chain']=='A' and 3<=r['res']<=500]
b=[r for r in records if r['chain']=='B']
contacts={}
for x in a:
    for y in b:
        d=math.dist(x['xyz'],y['xyz'])
        if d>8:continue
        k=(x['res'],x['icode'],x['name'])
        if k not in contacts or d<contacts[k]['distance']:
            contacts[k]=dict(res=x['res'],icode=x['icode'],name=x['name'],atom=x['atom'],peptideRes=y['res'],peptideName=y['name'],peptideAtom=y['atom'],distance=d)
rows=sorted(contacts.values(),key=lambda x:x['distance'])
for row in rows:row['distance']=round(row['distance'],4)
summary={'id':'9NU4','source':'https://files.rcsb.org/download/9NU4.pdb','citation':'https://doi.org/10.2210/pdb9NU4/pdb','license':'CC0-1.0','licenseSource':'https://www.wwpdb.org/about/usage-policies','retrieved':'2026-10-18','sha256':hashlib.sha256(pdb.read_bytes()).hexdigest(),'method':'Minimum inter-chain heavy-atom Euclidean distance; model 1; altloc blank/A; receptor A:3-500; observed peptide B; contacts stored to 8 A; no hydrogen-bond or affinity assignment.','chains':{c:{'atoms':len([r for r in records if r['chain']==c]),'residues':len({(r['res'],r['icode']) for r in records if r['chain']==c}),'min':min(r['res'] for r in records if r['chain']==c),'max':max(r['res'] for r in records if r['chain']==c)} for c in ['A','B']},'contacts':rows,'results':{'generation':None,'optimization':None,'docking':None,'screening':None,'md':None},'publicApproval':{'basis':'Public PDB CC0 coordinates and independently derived geometry only','internalDataIncluded':False}}
(ROOT/'public/data/sample.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
# A static projection of actual C-alpha coordinates for the lightweight home page.
ca=[r for r in records if r['atom']=='CA' and (r['chain']=='B' or 3<=r['res']<=500)]
c=[sum(r['xyz'][i] for r in ca)/len(ca) for i in range(3)]
svg=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 560"><title>9NU4 C-alpha coordinate projection</title>']
for chain,color in [('A','#70a991'),('B','#e1a859')]:
    pts=[r for r in ca if r['chain']==chain]
    for p,q in zip(pts,pts[1:]):
        if math.dist(p['xyz'],q['xyz'])>5:continue
        def proj(r):
            x,y,z=[r['xyz'][i]-c[i] for i in range(3)]
            return 300+(x*.83+z*.55)*5.8,280-y*5.8
        x,y=proj(p);u,v=proj(q)
        svg.append(f'<path d="M{x:.2f} {y:.2f}L{u:.2f} {v:.2f}" stroke="{color}" stroke-width="{5 if chain=="B" else 3}" fill="none" stroke-linecap="round"/>')
svg.append('</svg>');(ROOT/'public/structure.svg').write_text(''.join(svg))
print(json.dumps({'chains':summary['chains'],'contacts_4A':sum(r['distance']<=4 for r in rows),'sha256':summary['sha256']}))
