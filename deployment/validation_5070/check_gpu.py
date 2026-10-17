"""Actual small CUDA/PyG operations. No model weights or target data are loaded."""
import json
import torch
from torch_scatter import scatter_add
from torch_cluster import radius_graph
from torch_sparse import SparseTensor

assert torch.cuda.is_available(), 'CUDA unavailable'
x=torch.tensor([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],device='cuda')
v=scatter_add(torch.ones(3,device='cuda'),torch.tensor([0,0,1],device='cuda'),dim=0)
assert torch.allclose(v,torch.tensor([2.,1.],device='cuda'))
edges=radius_graph(x,r=1.5,loop=False)
assert edges.shape==(2,6), edges
m=SparseTensor(row=edges[0],col=edges[1],value=torch.ones(6,device='cuda'),sparse_sizes=(3,3))
y=m.matmul(torch.ones((3,2),device='cuda'))
assert torch.allclose(y,torch.full((3,2),2.,device='cuda'))
torch.cuda.synchronize()
print(json.dumps({'cuda_ops':'passed','torch':torch.__version__,'cuda_runtime':torch.version.cuda,'gpu':torch.cuda.get_device_name(0),'capability':torch.cuda.get_device_capability(0),'total_memory_GiB':round(torch.cuda.get_device_properties(0).total_memory/2**30,2),'model_inference':'not_run'},indent=2))
