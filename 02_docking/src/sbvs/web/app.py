import os,csv,json,secrets,time,uuid,base64
from pathlib import Path
from flask import Flask,render_template,request,jsonify,abort,send_file,url_for
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from ..paths import ROOT,project,rootpath,read_json,write_json
from ..config import load
from ..runner import create
from ..queue import Queue
from ..export import read_rows,report,handoff
from ..doctor import doctor
def make_app():
 app=Flask(__name__);app.config.update(MAX_CONTENT_LENGTH=30*1024*1024,SECRET_KEY=secrets.token_hex(32));q=Queue();app.queue=q
 @app.before_request
 def local_security():
  if request.method=='POST':
   origin=request.headers.get('Origin')
   if origin and origin.rstrip('/')!=request.host_url.rstrip('/'):abort(403)
   if request.headers.get('X-SBVS-Request')!='local-ui':abort(403)
 @app.errorhandler(ValueError)
 def value_error(e):return jsonify(error=str(e)),400
 @app.errorhandler(KeyError)
 def missing(e):return jsonify(error='Unknown run or field'),404
 @app.get('/')
 def home():return render_template('index.html',page='quick')
 @app.get('/classroom')
 def classroom():return render_template('index.html',page='classroom')
 @app.get('/research')
 def research():return render_template('index.html',page='research')
 @app.get('/history')
 def history():return render_template('history.html',jobs=q.all())
 @app.get('/environment')
 def environment():return render_template('environment.html',info=json.dumps(doctor(),ensure_ascii=False,indent=2))
 @app.get('/api/health')
 def health():return jsonify(project_uuid=project()['uuid'],ready=(ROOT/'reports/class_ready.json').exists(),accepting=q.accepting)
 @app.post('/api/jobs')
 def submit():
  preset=request.form.get('preset','class_quick');cfg=load(preset)
  receptor_id=request.form.get('receptor_id','demo_1drf')
  if receptor_id!='demo_1drf':
   if not receptor_id.isalnum():raise ValueError('Invalid receptor ID')
   folder=ROOT/'data/receptors'/receptor_id
   review=read_json(folder/'review.json')
   if not review.get('confirmed'):raise ValueError('Receptor review required')
   from ..paths import sha
   if request.form.get('pocket_confirmed')!='true' or request.form.get('reviewed_box_hash')!=sha(folder/'box.json'):raise ValueError('Inspect and confirm the current pocket preview before screening')
   cfg.update(receptor=str((folder/'receptor.pdbqt').relative_to(ROOT)),receptor_pdb=str((folder/'receptor_prepared.pdb').relative_to(ROOT)),box=str((folder/'box.json').relative_to(ROOT)),preparation_record=str((folder/'preparation.json').relative_to(ROOT)))
  overrides=request.form.get('overrides','').strip()
  if overrides:
   data=json.loads(overrides)
   allowed={'energy_range','exhaustiveness','num_modes','seeds','cpu_per_job','per_ligand_timeout_seconds','ligand_states','ph_low','ph_high','max_states','fragment_policy','conformers_generated','conformers_retained'}
   if set(data)-allowed:raise ValueError('Unsupported override')
   cfg.update(data)
  if cfg['engine']!='vina':raise ValueError('New jobs require official Vina')
  inbox=ROOT/'data/uploads'/uuid.uuid4().hex;inbox.mkdir(parents=True)
  upload=request.files.get('file')
  if upload and upload.filename:
   suffix=Path(upload.filename).suffix.lower()
   if suffix not in ('.smi','.txt','.csv','.sdf'):raise ValueError('Upload SMI/TXT/CSV/SDF')
   source=inbox/('library'+suffix);upload.save(source)
  else:
   text=request.form.get('smiles','').strip()
   if not text:raise ValueError('Ligand input cannot be empty')
   source=inbox/'library.smi';source.write_text(text,encoding='utf-8')
  out=create(source,cfg);rid=q.add(out);return jsonify(id=rid)
 @app.get('/api/jobs/<rid>')
 def status(rid):
  job=q.get(rid);out=ROOT/job['path'];m=read_json(out/'run_manifest.json');m['status']=job['status'];ranking=read_rows(out/'ranking.csv') if (out/'ranking.csv').exists() else []
  return jsonify(id=rid,status=job['status'],counts=m['counts'],elapsed_seconds=m.get('elapsed_seconds'),ranking=ranking,failures=m.get('failures',[]),config=m['config'],manifest=m)
 @app.post('/api/jobs/<rid>/cancel')
 def cancel(rid):q.get(rid);q.cancel(rid);return jsonify(ok=True)
 @app.post('/api/jobs/<rid>/resume')
 def resume(rid):
  job=q.get(rid)
  if job['status'] not in ('interrupted','cancelled','partial_failed','failed'):raise ValueError('Run cannot resume in current state')
  q.add(ROOT/job['path'],True);return jsonify(ok=True)
 @app.get('/api/jobs/<rid>/artifact/<name>')
 def artifact(rid,name):
  out=ROOT/q.get(rid)['path']
  allowed={'ranking.csv','states.csv','poses.csv','failures.csv','run_manifest.json','config_expanded.yaml','report.html','ligand_manifest.csv'}
  if name not in allowed:abort(404)
  if name=='report.html':report(out)
  return send_file(out/name,as_attachment=name!='report.html')
 @app.get('/api/jobs/<rid>/pose/<poseid>')
 def pose(rid,poseid):
  out=ROOT/q.get(rid)['path'];rows=read_rows(out/'poses.csv');matches=[r for r in rows if r['pose_id']==poseid]
  if len(matches)!=1:abort(404)
  p=(out/matches[0]['pose_file']).resolve()
  if not p.is_relative_to(out):abort(403)
  return send_file(p,mimetype='chemical/x-mdl-sdfile')
 @app.get('/api/jobs/<rid>/receptor')
 def receptor(rid):
  m=read_json(ROOT/q.get(rid)['path']/'run_manifest.json');return send_file(rootpath(m['config']['receptor_pdb']),mimetype='text/plain')
 @app.get('/api/jobs/<rid>/box')
 def box(rid):return jsonify(read_json(rootpath(read_json(ROOT/q.get(rid)['path']/'run_manifest.json')['config']['box'])))
 @app.get('/api/jobs/<rid>/poses')
 def poses(rid):return jsonify(read_rows(ROOT/q.get(rid)['path']/'poses.csv'))
 @app.get('/api/jobs/<rid>/images/<name>')
 def report_image(rid,name):
  out=ROOT/q.get(rid)['path'];allowed={r['pose_id']+'.png' for r in read_rows(out/'poses.csv')}
  if name not in allowed:abort(404)
  return send_file(out/'images'/name,mimetype='image/png')
 @app.get('/api/classroom/evidence/<step>')
 def evidence(step):
  files={'environment':'reports/doctor.json','receptor':'data/sources.json','preparation':'data/prepared/demo_1drf/preparation.json','box':'data/prepared/demo_1drf/box.json','library':'data/demo/ligands_20.smi'}
  if step not in files:abort(404)
  if step=='environment':doctor()
  return jsonify(path=files[step],text=(ROOT/files[step]).read_text())
 @app.post('/api/receptors/plan')
 def receptor_plan():
  from ..preparation.receptor import plan
  f=request.files.get('receptor')
  if not f or Path(f.filename).suffix.lower() not in ('.pdb','.cif','.mmcif','.pdbqt'):raise ValueError('Upload PDB/mmCIF/PDBQT')
  rid=uuid.uuid4().hex;folder=ROOT/'data/receptors'/rid;folder.mkdir(parents=True);source=folder/('original'+Path(f.filename).suffix.lower());f.save(source)
  chains=request.form.get('chains','A').split(',');model=int(request.form.get('model_index',0));review=plan(source,chains,model)
  review.update(source=source.name,receptor_id=rid);write_json(folder/'plan.json',review);return jsonify(review)
 @app.post('/api/receptors/<receptor_id>/prepare')
 def receptor_prepare(receptor_id):
  from ..preparation.receptor import prepare
  from ..preparation.box import make_box
  from ..config import validate_box
  from ..paths import sha
  if not receptor_id.isalnum():raise ValueError('Invalid ID')
  folder=ROOT/'data/receptors'/receptor_id;plan=read_json(folder/'plan.json');policy=json.loads(request.form['policy'])
  if policy.get('confirmed') is not True or not policy.get('box_confirmed'):raise ValueError('Confirm preparation decisions and pocket')
  policy['input_sha256']=plan['input_sha256'];policy['model_index']=plan.get('model_index',0)
  prepare(folder/plan['source'],folder,plan.get('chains',['A']),policy)
  if plan.get('external'):
   rec=request.files.get('prepared_pdb')
   if not rec:raise ValueError('External PDBQT requires corresponding prepared PDB for visualization/Vina')
   rec.save(folder/'receptor_prepared.pdb')
  mode=request.form.get('box_mode','manual')
  if mode=='reference':
   ref=request.files.get('reference')
   if not ref:raise ValueError('Reference SDF required')
   ref.save(folder/'reference.sdf');box=make_box(folder/'reference.sdf',float(request.form.get('padding',5)))
  elif mode=='residues':
   import numpy as np
   selected=set(policy.get('pocket_residues',[]));coords=[];found=set()
   for line in (folder/'receptor_prepared.pdb').read_text().splitlines():
    key=line[21:22]+':'+line[22:26].strip()+line[26:27].strip()
    if line.startswith(('ATOM  ','HETATM')) and key in selected:
     coords.append([float(line[30:38]),float(line[38:46]),float(line[46:54])]);found.add(key)
   if not selected or found!=selected:raise ValueError('Pocket residue selection missing')
   xyz=np.asarray(coords);center=xyz.mean(0);box={'center':center.tolist(),'size':(2*np.max(abs(xyz-center),axis=0)+10).tolist(),'method':'selected_residues_geometric_center','padding_per_side':5,'residues':sorted(found)}
  elif mode=='manual':box=policy['box']
  else:raise ValueError('Unknown box mode')
  validate_box(box)
  if max(box['size'])>40 and not policy.get('large_box_confirmed'):raise ValueError('Large box needs explicit large_box_confirmed')
  box['receptor_sha256']=sha(folder/'receptor.pdbqt');write_json(folder/'box.json',box);write_json(folder/'review.json',policy);return jsonify(id=receptor_id,box=box,box_hash=sha(folder/'box.json'),status='prepared')
 @app.get('/api/receptors/<receptor_id>/preview')
 def receptor_preview(receptor_id):
  from ..paths import sha
  if not receptor_id.isalnum():raise ValueError('Invalid ID')
  folder=ROOT/'data/receptors'/receptor_id
  if not read_json(folder/'review.json').get('confirmed'):raise ValueError('Preparation review incomplete')
  return jsonify(pdb=(folder/'receptor_prepared.pdb').read_text(),reference=(folder/'reference.sdf').read_text() if (folder/'reference.sdf').exists() else None,box=read_json(folder/'box.json'),box_hash=sha(folder/'box.json'))
 @app.get('/api/receptors')
 def receptor_list():
  values=[{'id':'demo_1drf','name':'Reviewed classroom 1DRF'}]
  for p in (ROOT/'data/receptors').glob('*/review.json'):
   if read_json(p).get('confirmed'):values.append({'id':p.parent.name,'name':p.parent.name})
  return jsonify(values)
 @app.post('/api/jobs/<rid>/image')
 def image(rid):
  out=ROOT/q.get(rid)['path'];d=request.get_json();pid=d['pose_id']
  if pid not in {r['pose_id'] for r in read_rows(out/'poses.csv')}:raise ValueError('Unknown pose')
  blob=base64.b64decode(d['png'].split(',',1)[1],validate=True)
  if not blob.startswith(b'\x89PNG\r\n\x1a\n') or len(blob)>10*1024*1024:raise ValueError('Invalid PNG')
  images=out/'images';images.mkdir(exist_ok=True);(images/(pid+'.png')).write_bytes(blob);report(out);return jsonify(ok=True)
 @app.post('/api/jobs/<rid>/handoff')
 def export_md(rid):return jsonify(path=str(handoff(ROOT/q.get(rid)['path'],request.get_json()).relative_to(ROOT)))
 @app.post('/api/jobs/<rid>/vina-stage')
 def stage(rid):
  from ..stages import vina_stage
  data=request.get_json();device='cpu'
  out=vina_stage(ROOT/q.get(rid)['path'],data['operation'],int(data.get('top_k',10)),device);return jsonify(id=q.add(out))
 @app.post('/api/jobs/<rid>/sort')
 def sort_result(rid):
  from ..ranking import csvwrite
  data=request.get_json();metric=data['metric'];direction={'vina_affinity_kcal_mol':1,'empirical_score_kcal_mol':1,'cnn_pose_score':-1,'cnn_affinity_pK':-1}
  if metric not in direction:raise ValueError('Unknown score metric')
  out=ROOT/q.get(rid)['path'];rows=read_rows(out/'ranking.csv')
  if any(not r.get(metric) for r in rows):raise ValueError('Metric absent for this engine')
  rows.sort(key=lambda r:(direction[metric]*float(r[metric]),r['parent_id']))
  for i,r in enumerate(rows,1):r['original_rank']=r['rank'];r['rank']=i
  csvwrite(out/'ranking_user.csv',rows);write_json(out/'user_sort.json',{'metric':metric,'direction':'ascending' if direction[metric]==1 else 'descending','created':time.time(),'note':'Sort existing real parent representatives; original ranking.csv preserved'})
  return jsonify(ranking=rows)
 @app.post('/api/drain')
 def drain():return jsonify(idle=q.drain())
 return app
