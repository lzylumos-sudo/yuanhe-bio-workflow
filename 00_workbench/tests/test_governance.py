"""Governance contracts and strict scientific configuration regression; no model/GPU run."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import pytest
from yuanhe.version import ROOT, software_version
from yuanhe.registry import validate, load, snapshot
from yuanhe.configuration import integration_steps, merge_layers, protocol_for, stage_config, fingerprint
from yuanhe.protocols import validate_sources, compatible
from yuanhe.records import new_run_id, sync
from webui import workflow, server


def test_registry_maps_to_real_exports_and_protocols():
    assert validate()['modules'] == 6
    assert validate_sources() == {'valid': True}
    assert next(x for x in load()['backends'] if x['backend_id']=='ymavs')['parent']=='zvss'
    assert compatible('openmm_soluble','soluble-protein-md@1.0.0')
    with pytest.raises(ValueError): compatible('gromacs_membrane','soluble-protein-md@1.0.0')


def test_missing_environment_never_shows_compute_available():
    registry=snapshot({'chem_python':'/missing/python','md_python':'/missing/python','prep_python':'/missing/python'})
    assert all(not b['runnable'] for b in registry['backends'] if b['backend_id']!='workbench')
    assert all(not b['verified'] for b in registry['backends'])


@pytest.mark.parametrize('duration,dt,expected',[(500,2,250000000),(1,2,500000),(0.002,2,1000)])
def test_duration_unit_conversion(duration,dt,expected):
    assert integration_steps(duration,dt)==expected


@pytest.mark.parametrize('duration,dt',[(500,0),(float('nan'),2),(500,float('inf')),(0.000001,2)])
def test_bad_duration_blocks_instead_of_rounding(duration,dt):
    with pytest.raises(ValueError):integration_steps(duration,dt)


def test_membrane_mdp_exactly_matches_registered_pre_refactor_recipe():
    from scientific.membrane import mdp
    recipe=next(p for p in load('protocol_registry')['protocols'] if p['protocol_id']=='murj-membrane-md')['mdp_recipe_source']
    scope={};exec(recipe,scope)
    for temperature in (300,310):
        for phase in ('em','nvt','npt','npt_soft','production'):
            params={'temperature_K':temperature,'production_ns':500,'equilibration_ns':2,'output_interval_ps':10}
            assert mdp(phase,params,17)==scope['mdp'](phase,params,17)
    production=mdp('production',params,17)
    assert 'nsteps = 250000000' in production and 'dt = 0.002' in production
    assert 'pcoupltype = semiisotropic' in production and 'ref-p = 1.0 1.0' in production


def test_research_default_and_override_protocol_identity():
    cfg=workflow.defaults();assert protocol_for(cfg)['conforms']
    cfg['research']['md_production_ns']=5
    assert protocol_for(cfg)['protocol_id']=='custom'
    assert protocol_for(cfg)['protocol_version']=='unknown'
    assert workflow.defaults()['research']['md_production_ns']==500


def test_layer_order_and_stage_compatibility():
    result=merge_layers([('base',{'a':{'b':1,'c':2}}),('target',{'a':{'b':3}}),('override',{'a':{'b':4}})])
    assert result['config']=={'a':{'b':4,'c':2}}
    assert result['sources']['a.b']=='override'
    cfg={'a':{'b':1,'c':2}}
    assert stage_config(cfg,[('prepare',{'a':{'b':4}})])['config']=={'a':{'b':4}}
    assert cfg=={'a':{'b':1,'c':2}}


def test_new_run_manifest_records_creation_and_preserves_historical_run(tmp_path,monkeypatch):
    monkeypatch.setattr(workflow,'RUNS',tmp_path)
    rid=workflow.create({})['id'];d=tmp_path/rid
    manifest=json.loads((d/'run_manifest.json').read_text())
    assert rid.startswith('YH-') and manifest['software_version']==software_version()
    assert len(manifest['git_commit'])==40 and manifest['git_dirty'] is not None
    assert manifest['protocol_id']=='murj-discovery'
    assert manifest['config_hash']==fingerprint(json.loads((d/'config.json').read_text()))
    assert manifest['stage_status']=={} and manifest['scientific_review']=='not_reviewed'
    assert manifest['model_versions']['peptide2mol']['execution_verified'] is False
    assert json.loads((d/'resolved_config.yaml').read_text())['config']==json.loads((d/'config.json').read_text())
    old=tmp_path/'historical';old.mkdir();(old/'run.json').write_text('{"version":"original"}')
    assert sync(old,['prepare']) is None
    assert (old/'run.json').read_text()=='{"version":"original"}'
    assert not (old/'run_manifest.json').exists()
    assert workflow.status(rid)['manifest']['run_id']==rid


def test_execution_completion_is_not_scientific_acceptance(tmp_path,monkeypatch):
    monkeypatch.setattr(workflow,'RUNS',tmp_path)
    rid=workflow.create({})['id'];d=tmp_path/rid
    for s,_,_ in workflow.STAGES:workflow.atomic(d/s/'stage.json',{'status':'completed','artifacts':{}})
    manifest=sync(d,[s[0] for s in workflow.STAGES])
    assert manifest['status']=='completed' and manifest['scientific_review']=='review_required'
    workflow.atomic(d/'gromacs/stage.json',{'status':'failed','artifacts':{}})
    assert sync(d,[s[0] for s in workflow.STAGES])['status']=='failed'


def test_api_version_and_errors_remain_compatible():
    status,_,payload=server.dispatch('GET','/api/v1/versions',{},None)
    assert status==200 and json.loads(payload)['software_version']==software_version()
    status,_,payload=server.dispatch('POST','/api/v1/validate',{}, {'unknown_field':True})
    body=json.loads(payload);assert status==400 and body['error']==body['message'] and body['error_code']


def test_generated_package_versions_are_current():
    subprocess.run([sys.executable,str(ROOT/'deployment/sync_platform_version.py'),'--check'],check=True)


def test_scores_missing_is_not_zero():
    from yuanhe.contracts import docking_scores
    assert docking_scores(engine='vina',vina_score=-7.2)=={'vina_score':-7.2}
    with pytest.raises(ValueError):docking_scores(engine='vina',cnn_score=0)


def test_public_docking_is_same_function():
    sys.path.insert(0,str(ROOT/'02_docking/src'))
    from docking import run_vina_job
    from sbvs2.docking.vina_official import run_vina_job as original
    assert run_vina_job is original


def test_run_ids_unique():
    assert len({new_run_id() for _ in range(1000)})==1000


def test_native_soluble_and_membrane_protocols_remain_distinct():
    import yaml
    from yuanhe.protocols import soluble,membrane
    cfg=yaml.safe_load((ROOT/'05_molecular_dynamics/config/default.yaml').read_text())
    assert soluble(cfg)['protocol_id']=='soluble-protein-md'
    cfg['MD']['temperature']=310
    assert soluble(cfg)['protocol_id']=='custom'
    assert membrane({'temperature_K':310,'production_ns':500,'equilibration_ns':2})['protocol_id']=='murj-membrane-md'
    assert membrane({'temperature_K':310,'production_ns':500,'equilibration_ns':2,'lipid':'POPC'})['protocol_id']=='custom'


def test_invalid_configuration_type_is_input_error():
    from research.config import resolve
    for value in ([],False,'500 ns'):
        with pytest.raises(ValueError):resolve(value)


def test_missing_worker_marks_manifest_failed_and_releases_lease(tmp_path,monkeypatch):
    monkeypatch.setattr(workflow,'RUNS',tmp_path)
    rid=workflow.create({'chem_python':'/missing/worker/python'})['id']
    with pytest.raises(FileNotFoundError):workflow.start(rid,'prepare')
    manifest=json.loads((tmp_path/rid/'run_manifest.json').read_text())
    assert manifest['status']=='failed'
    lease=workflow.acquire_compute();lease.close()


def test_scientific_baseline_byte_identity():
    subprocess.run([sys.executable,str(ROOT/'deployment/check_scientific_baseline.py')],check=True)


def test_scientific_workflow_overrides_cannot_claim_original_protocol():
    cfg=workflow.defaults();cfg['energy_range']=4
    assert not protocol_for(cfg)['conforms']
    assert 'energy_range' in protocol_for(cfg)['differences']


def test_configured_target_cannot_escape_registry_directory():
    with pytest.raises(ValueError):workflow.validate_config({'target_profile':'../../outside'})


def test_registry_preserves_legacy_inventory_shape():
    registry=load()
    assert {'Peptide2Mol','REINVENT4','Vina adapter','Virtual screening','Membrane MD'}=={r['name'] for r in registry['legacy_inventory']}
    assert all({'name','path'} <= r.keys() for r in registry['legacy_inventory'])


def test_stage_contracts_reject_incomplete_or_false_review():
    from yuanhe.contracts import validate_stage
    with pytest.raises(ValueError):validate_stage({'run_id':'test'},'stage-input')
    with pytest.raises(ValueError):validate_stage({'run_id':'test','stage_id':'docking','config_hash':'a'*64,'output_artifacts':[],'scientific_review':'PASS'},'stage-output')
