"""Exercise the vendor sampling loop with failed/empty decoding, no inference."""
import ast,json,os,traceback,tempfile
from pathlib import Path
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[2]
source=ROOT/'01_molecule_generation/vendor/Peptide2Mol/src/models/components/moldiff_test/model.py'
module=ast.parse(source.read_text());method=next(n for n in ast.walk(module) if isinstance(n,ast.FunctionDef) and n.name=='sample')
code=compile(ast.fix_missing_locations(ast.Module(body=[method],type_ignores=[])),str(source),'exec')
class AttrDict(dict):
 __getattr__=dict.__getitem__
class Tensor:
 def cpu(self):return self
 def numpy(self):return self
class Fake:
 sample_batch_size=1;num_mols=2;sample_max=20;guidance=None
 calls=0
 def make_mydata_placeholder(self,**kw):return dict(batch_node=Tensor(),halfedge_index=Tensor(),batch_halfedge=Tensor(),ref_data=SimpleNamespace(name=['test']),n_nodes_list=[])
 def seperate_outputs(self,*a):return []
 def sample_model(self,**kw):self.calls+=1;return {}
with tempfile.TemporaryDirectory() as tmp:
 fake=Fake();fake.log_dir=str(Path(tmp)/'model');fake.model=SimpleNamespace(sample=fake.sample_model)
 Path(fake.log_dir+'_SDF').mkdir();Path(fake.log_dir+'_SDF/test_origin.pdb').write_text('test')
 ns={'EasyDict':AttrDict,'os':os,'json':json,'traceback':traceback};exec(code,ns)
 ns['sample'](fake,SimpleNamespace(pos=SimpleNamespace(device='cpu')))
 audit=json.loads(Path(fake.log_dir+'_SDF/attempt_budget.json').read_text())
 assert fake.calls==2 and audit['started_attempts']==2 and audit['valid']==0
print('PASS: empty decoding stops after exact requested attempts; no model inference')
