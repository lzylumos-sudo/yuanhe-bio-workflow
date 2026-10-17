from types import SimpleNamespace
from zvss.select.qualify import _pose_qa
from zvss.redock.gnina import redock_argv

def test_failed_qa_retains_reason_without_crashing_selection():
 class DB:
  run_id='test'
  def query(self,*a):return []
  def one(self,*a):return {'assess_status':'fail','reason_code':'bad_geometry'}
 result=_pose_qa(SimpleNamespace(db=DB()),{'pose_id':'pose'})
 assert result['status']=='fail' and result['reason_code']=='bad_geometry'

def test_redock_pins_model_and_device():
 args=dict(command=['gnina'],receptor='receptor.pdb',ligand='ligand.sdf',center=[1,2,3],size=[20,20,20],seed=17,out_prefix='out.sdf',cnn_model='crossdock_default2018_1_3')
 gpu=redock_argv(**args,device='gpu');cpu=redock_argv(**args,device='cpu')
 assert gpu[gpu.index('--cnn')+1]=='crossdock_default2018_1_3'
 assert gpu[gpu.index('--device')+1]=='0'
 assert '--no_gpu' in cpu and '--device' not in cpu

def test_changing_gnina_model_invalidates_host_cache():
 from zvss.config import HostConfig
 a=HostConfig(raw={'tools':{'gnina_settings':{'cnn_model':'model_a'}}},source_path='host.json')
 b=HostConfig(raw={'tools':{'gnina_settings':{'cnn_model':'model_b'}}},source_path='host.json')
 assert a.digest()!=b.digest()

def test_live_counts_are_unique_compounds_not_poses():
 from zvss.stages import live_counts
 rows=[{'compound_id':'A','assess_status':'pass'},{'compound_id':'A','assess_status':'fail'},{'compound_id':'B','assess_status':'fail'},{'compound_id':'C','assess_status':'not_evaluated'}]
 c=live_counts(rows,qa_finished=True)
 assert c=={'imported':3,'qa_usable':1,'qa_pending':1,'qa_failed':1}
 c=live_counts(rows,{'primary':[{'compound_id':'A'}],'preliminary':[{'compound_id':'A'},{'compound_id':'C'}]},True)
 assert c['remaining']==2 and c['formal']==1 and c['preliminary']==1
