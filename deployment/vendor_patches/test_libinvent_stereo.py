from rdkit import Chem
from reinvent.chemistry.library_design import attachment_points, bond_maker
from reinvent.runmodes.samplers.libinvent import LibinventSampler
from reinvent.models.model_factory.sample_batch import SampleBatch
from unittest.mock import patch
import torch

class FakeDataset:
    model_inputs=[]
    def __init__(self, inputs, *_): self.model_inputs[:]=inputs; self.inputs=inputs
    def __len__(self): return len(self.inputs)
    def __getitem__(self,i): return i
    @staticmethod
    def collate_fn(items): return items,None
class FakeVocabulary:
    scaffold_vocabulary=None
    scaffold_tokenizer=None
class FakeModel:
    version=1
    def get_vocabulary(self): return FakeVocabulary()
    def sample(self,inputs,_):
        return SampleBatch([FakeDataset.model_inputs[i] for i in inputs], ['[*]C(=O)C']*len(inputs),torch.ones(len(inputs)))

scaffolds=['N([*])C[C@H](O)CC','N([*])C[C@@H](O)CC']
with patch('reinvent.runmodes.samplers.libinvent.Dataset',FakeDataset):
    sampled=LibinventSampler(FakeModel(),batch_size=2,randomize_smiles=False).sample(scaffolds)
assert list(sampled.items1)==scaffolds*2
assert all('@' not in s for s in FakeDataset.model_inputs)
products=list(sampled.smilies)
assert len(products)==4 and products[0]!=products[1]
for original, product in zip(scaffolds*2,products):
    core=Chem.DeleteSubstructs(Chem.MolFromSmiles(original),Chem.MolFromSmarts('[#0]'))
    assert Chem.MolFromSmiles(product).HasSubstructMatch(core,useChirality=True)
# Actual recovered scaffold with four stereocentres: attachment is achiral N.
original='CCc1nc(CC[C@@H]2CN([*])C[C@@H](O)[C@H]3O[C@H]23)cn1C'
with patch('reinvent.runmodes.samplers.libinvent.Dataset',FakeDataset):
    sampled=LibinventSampler(FakeModel(),batch_size=8,randomize_smiles=False).sample([original])
core=Chem.DeleteSubstructs(Chem.MolFromSmiles(original),Chem.MolFromSmarts('[#0]'))
assert len(sampled.smilies)==8
for product in sampled.smilies:
    assert Chem.MolFromSmiles(product).HasSubstructMatch(core,useChirality=True)
print('PASS: original input order, 8 sample bound, opposite stereo distinguishable, recovered core stereo retained; deterministic mock decorations only, no model sampling')
