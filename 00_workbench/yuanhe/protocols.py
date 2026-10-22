"""Scientific protocol identity without relabeling historical trajectories."""
import ast
import hashlib
from .registry import load
from .configuration import integration_steps, fingerprint
from .version import ROOT


def get(identifier):
    name, sep, version = identifier.partition('@')
    for p in load('protocol_registry')['protocols']:
        if p['protocol_id'] == name and sep and p['protocol_version'] == version: return p
    raise ValueError('Unknown scientific protocol: ' + identifier)

def compatible(backend, identifier):
    if backend not in get(identifier)['backend']: raise ValueError('Backend is not compatible with protocol')
    return True

def membrane(parameters):
    p = get('murj-membrane-md@1.0.0'); expected = p['parameters']
    actual = {'temperature_K': parameters['temperature_K'], 'production_ns': parameters['production_ns'],
              'equilibration_ns': parameters['equilibration_ns'],
              'lipid_default': parameters.get('lipid','POPE'), 'salt_M': parameters.get('salt_M',.15),
              'padding_nm': parameters.get('padding_nm',1.2),
              'alignment_selection': parameters.get('alignment_selection','protein and backbone')}
    differences = {k: {'protocol': expected[k], 'resolved': v} for k, v in actual.items() if expected[k] != v}
    return {'protocol_id': 'murj-membrane-md' if not differences else 'custom',
            'protocol_version': '1.0.0' if not differences else 'unknown',
            'base_protocol': 'murj-membrane-md@1.0.0', 'differences': differences,
            'resolved_parameters_hash': fingerprint(parameters),
            'production_steps': integration_steps(parameters['production_ns'], 2), 'timestep_fs': 2,
            'pressure': {'value': 1, 'unit': 'bar'}, 'pressure_coupling': 'semiisotropic'}

def validate_sources():
    for p in load('protocol_registry')['protocols']:
        source = ROOT / p['source']
        if 'mdp_recipe_source' in p:
            text = source.read_text(); tree = ast.parse(text)
            node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'mdp')
            if ast.get_source_segment(text, node) != p['mdp_recipe_source']:
                raise ValueError('Membrane recipe changed without protocol revision')
        elif hashlib.sha256(source.read_bytes()).hexdigest() != p['source_sha256']:
            raise ValueError('Scientific defaults changed without protocol revision: ' + p['source'])
    return {'valid': True}


def soluble(config):
    import yaml
    p=get('soluble-protein-md@1.0.0'); original=yaml.safe_load(p['source_snapshot'])
    from .configuration import merge_layers
    # Exclude only hardware/concurrency controls from scientific comparison.
    ignored={'platform','require_gpu','device_index','cpu_threads','qm_threads'}
    differences={}
    for group in ('MD','analysis','mmpbsa','mmgbsa','protein'):
        baseline=original[group]; actual=config.get(group,{})
        for key in set(baseline)|set(actual):
            if group=='MD' and key in ignored:continue
            if baseline.get(key)!=actual.get(key):differences[group+'.'+key]={'protocol':baseline.get(key),'resolved':actual.get(key)}
    if config.get('replicates')!=original['replicates']:differences['replicates']={'protocol':original['replicates'],'resolved':config.get('replicates')}
    return {'protocol_id':'soluble-protein-md' if not differences else 'custom',
            'protocol_version':'1.0.0' if not differences else 'unknown',
            'base_protocol':'soluble-protein-md@1.0.0','differences':differences,'resolved_parameters_hash':fingerprint(config)}
