import os,sys,json,requests,hashlib,importlib.metadata,time
from pathlib import Path
from sbvs.paths import ROOT,read_json,write_json,sha
lic=ROOT/'docs/licenses';lic.mkdir(parents=True,exist_ok=True);sources=[]
for repo,ref in [('gnina/gnina','v1.3.3'),('ccsb-scripps/AutoDock-Vina','v1.2.7'),('3dmol/3Dmol.js','2.5.5'),('forlilab/Meeko','v0.8.0')]:
 r=requests.get(f'https://api.github.com/repos/{repo}/contents?ref={ref}',timeout=30)
 if not r.ok:sources.append({'repo':repo,'ref':ref,'status':r.status_code});continue
 for item in r.json():
  if item['type']=='file' and any(word in item['name'].lower() for word in ['license','copying','notice']):
   response=requests.get(item['download_url'],timeout=30);response.raise_for_status();dst=lic/(repo.replace('/','_')+'_'+item['name']);dst.write_bytes(response.content);sources.append({'url':item['download_url'],'path':str(dst.relative_to(ROOT)),'sha256':sha(dst)})
packages=[]
for p in Path('/opt/sbvs/env/conda-meta').glob('*.json'):
 d=read_json(p);packages.append({k:d.get(k) for k in ['name','version','build','channel','url','sha256','license']})
write_json(ROOT/'env/package_provenance.json',packages)
pip=[]
for d in importlib.metadata.distributions():
 if d.metadata['Name'].lower().startswith(('nvidia','playwright','pyee','greenlet','molscrub')):pip.append({'name':d.metadata['Name'],'version':d.version,'source':'https://pypi.org/project/'+d.metadata['Name']+'/'+d.version+'/','license':d.metadata.get('License'),'installer':d.read_text('INSTALLER')})
write_json(ROOT/'env/pip_sources.json',pip)
lock=read_json(ROOT/'env/vendor.lock.json');lock['micromamba']={'version':'2.3.3','path':'/opt/sbvs/vendor/micromamba','sha256':sha('/opt/sbvs/vendor/micromamba'),'url':'https://micro.mamba.pm/api/micromamba/linux-64/2.3.3','publisher_sha256_verified':False};lock['3Dmol.js']={'version':'2.5.5','url':'https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js','sha256':sha(ROOT/'src/sbvs/web/static/3Dmol-min.js'),'source_kind':'versioned npm package via jsDelivr','publisher_sha256_verified':False};write_json(ROOT/'env/vendor.lock.json',lock)
write_json(ROOT/'docs/licenses/sources.json',sources)
text='# Third-party sources and licenses\n\n The distribution is not uniformly licensed as Apache-2.0. GNINA and contains the repository documented GPL/Apache licensed portions and dependencies; exact texts are stored in licenses/, Review these terms separately before redistribution.NVIDIA User-space libraries have separate licenses; GPU drivers do not migrate with the project.\n\n'
text+='Source license texts and URL/SHA256: docs/licenses/sources.json. Conda Each installed package version,build, source and license: env/package_provenance.json. pip Source:env/pip_sources.json. Engine and frontend assets:env/vendor.lock.json. \n\n'
text+='Ubuntu Base Source https://cdimage.ubuntu.com/ubuntu-base/releases/24.04/release/ , Original distribution stored in runtime/downloads/. RCSB PDB/CCD Sources and hashes in data/sources.json. 3Dmol.js 2.5.5 Localized from a pinned jsDelivr npm mirror; browsers do not access CDN. \n'
(ROOT/'docs/THIRD_PARTY_NOTICES.md').write_text(text,encoding='utf-8')
print(json.dumps({'license_files':len(sources),'conda_packages':len(packages),'pip_packages':len(pip)}))
