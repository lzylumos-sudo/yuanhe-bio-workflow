import csv,statistics
from pathlib import Path
FIELDS=['run_id','parent_id','original_name','state_id','conformer_id','pose_id','engine','engine_version','scoring_function','model_id','seed','vina_affinity_kcal_mol','empirical_score_kcal_mol','cnn_pose_score','cnn_affinity_pK','cnn_vs','prep_status','dock_status','qc_flags','receptor_hash','box_hash','pose_file']
def csvwrite(path,rows,fields=None):
 fields=fields or list(dict.fromkeys(FIELDS+[k for r in rows for k in r]))
 with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def aggregate(rows,engine):
 by={}
 for r in rows:by.setdefault((r['parent_id'],r['state_id'],r['seed']),[]).append(r)
 reps=[]
 for group in by.values():
  key=(lambda r:(r['vina_affinity_kcal_mol'],r['pose_id'])) if engine=='vina' else (lambda r:(-r['cnn_pose_score'],r['pose_id']))
  reps.append(min(group,key=key))
 states={}
 for r in reps:states.setdefault((r['parent_id'],r['state_id']),[]).append(r)
 state_rows=[];field='vina_affinity_kcal_mol' if engine=='vina' else 'cnn_affinity_pK'
 for group in states.values():
  values=[r[field] for r in group];median=statistics.median(values);representative=min(group,key=lambda r:(abs(r[field]-median),r['pose_id']))
  state_rows.append(dict(representative,aggregate_score=median,seed_count=len(values),seed_scores=__import__('json').dumps(values),aggregation='best pose across conformers per seed; median across seed representatives; actual pose nearest median'))
 parents={}
 for r in state_rows:parents.setdefault(r['parent_id'],[]).append(r)
 ranking=[];field='vina_affinity_kcal_mol' if engine=='vina' else 'cnn_affinity_pK'
 for group in parents.values():
  representative=min(group,key=lambda r:((r['aggregate_score'] if engine=='vina' else -r['aggregate_score']),r['state_id']))
  ranking.append(dict(representative,state_count=len(group),conformer_count=len({r['conformer_id'] for r in rows if r['parent_id']==representative['parent_id']}),aggregation='best state by seed-median score; '+representative['aggregation']))
 ranking.sort(key=lambda r:((r['aggregate_score'] if engine=='vina' else -r['aggregate_score']),r['parent_id']))
 for i,r in enumerate(ranking,1):r['rank']=i
 return state_rows,ranking
