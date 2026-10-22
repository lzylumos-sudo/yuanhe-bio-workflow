"""Small strict contracts at governance boundaries, without optional schema dependencies."""
import math
import re
from .records import STATES, REVIEW_STATES

def validate_manifest(value):
    required = {'run_id','software_version','git_commit','git_dirty','protocol_id','protocol_version',
                'config_hash','target_id','input_artifacts','model_versions','tool_versions',
                'stage_status','output_artifacts','start_time','end_time','status','scientific_review'}
    from .version import ROOT
    import json
    required=set(json.loads((ROOT/'configs/schemas/run_manifest.schema.json').read_text())['required'])
    if required - value.keys(): raise ValueError('Missing run manifest fields: ' + ', '.join(sorted(required-value.keys())))
    if value['status'] not in STATES or value['scientific_review'] not in REVIEW_STATES: raise ValueError('Invalid execution/review state')
    if value['git_commit'] != 'unknown' and not re.fullmatch('[0-9a-f]{40}', value['git_commit']): raise ValueError('Full Git commit SHA required')
    if not re.fullmatch('[0-9a-f]{64}', value['config_hash']): raise ValueError('Invalid configuration hash')
    if value['git_dirty'] is not None and type(value['git_dirty']) is not bool: raise ValueError('git_dirty must be boolean or unknown/null')
    for record in value['input_artifacts'] + value['output_artifacts']: validate_artifact(record)
    return value

def validate_artifact(value):
    required={'artifact_id','artifact_type','run_id','stage_id','source','content_hash','schema_version','validation_status','storage_path'}
    from .version import ROOT
    import json
    required=set(json.loads((ROOT/'configs/schemas/artifact.schema.json').read_text())['required'])
    if required-value.keys(): raise ValueError('Incomplete artifact identity')
    if not re.fullmatch('[0-9a-f]{64}',value['content_hash']): raise ValueError('Invalid artifact content hash')
    return value

def docking_scores(*, vina_score=None, gnina_score=None, cnn_score=None, cnn_affinity=None, engine):
    values={'vina_score':vina_score,'gnina_score':gnina_score,'cnn_score':cnn_score,'cnn_affinity':cnn_affinity}
    if engine == 'vina' and any(values[k] is not None for k in ('gnina_score','cnn_score','cnn_affinity')):
        raise ValueError('Vina execution cannot produce GNINA/CNN scores')
    for val in values.values():
        if val is not None and (type(val) not in (float,int) or not math.isfinite(val)): raise ValueError('Score must be finite or absent')
    return {k:v for k,v in values.items() if v is not None}


def validate_stage(value,kind):
    from .version import ROOT
    import json
    if kind not in {'stage-input','stage-output'}:raise ValueError('Unknown stage schema')
    schema=json.loads((ROOT/'configs/schemas'/(kind+'.schema.json')).read_text())
    if set(schema['required'])-value.keys():raise ValueError('Missing stage contract fields')
    if not re.fullmatch('[0-9a-f]{64}',value['config_hash']):raise ValueError('Invalid stage configuration hash')
    field='input_artifacts' if kind=='stage-input' else 'output_artifacts'
    if not isinstance(value[field],list):raise ValueError('Artifact list required')
    for row in value[field]:validate_artifact(row)
    if kind=='stage-output' and value['scientific_review'] not in REVIEW_STATES:raise ValueError('Invalid scientific review')
    return value
