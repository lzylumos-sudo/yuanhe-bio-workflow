"""Lazy generation backend selection; no cross-model algorithm substitution."""
from importlib import import_module
BACKENDS = {'peptide2mol': ('pocket_conditioned_generation', '.peptide2mol'),
            'reinvent4': ('scaffold_decoration_pretrained_sampling', '.reinvent4')}

def backend(name):
    if name not in BACKENDS: raise ValueError('Unknown generation backend: ' + name)
    return import_module(BACKENDS[name][1], __package__)

def generate(name, **kwargs):
    return backend(name).run(**kwargs)
