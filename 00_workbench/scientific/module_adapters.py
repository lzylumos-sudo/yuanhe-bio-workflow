"""Independent module adapters reuse the scientific engines and strict screening gates."""
import json
import math
from pathlib import Path
import shutil
import sys
from . import artifacts as a

ROOT = Path(__file__).resolve().parents[2]
for path in ['01_molecule_generation', '02_docking/src', '03_virtual_screening/src', '04_interaction_analysis/src']:
    sys.path.insert(0, str(ROOT/path))


def molecule(path):
    from rdkit import Chem
    rows = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if len(rows) != 1 or rows[0] is None or not rows[0].GetNumConformers():
        raise ValueError('Exactly one valid ligand with coordinates is required')
    conf = rows[0].GetConformer()
    if not conf.Is3D() or any(not math.isfinite(float(x)) for point in conf.GetPositions() for x in point):
        raise ValueError('Ligand requires finite 3D coordinates')
    return rows[0], Chem.MolToSmiles(Chem.RemoveHs(rows[0]), isomericSmiles=True)


def preflight(path, cfg, out, *, availability):
    from yuanhe_preflight.cli import import_candidates
    from yuanhe_preflight.pipeline import risk_stage, availability_stage
    rows = import_candidates(path, out/'imported')
    if not rows:
        raise ValueError('Empty library; SMILES tables require a smiles header')
    risk = risk_stage(rows, cfg['predocking'], out/'risk', base=out, batch_id=out.parent.parent.name)
    result = {'risk': risk}
    if availability:
        result['availability'] = availability_stage(a.read(out/'risk/records.json'), cfg['predocking'], out/'availability')
    return result


def screening(d, cfg, p, files, out):
    import yaml
    from zvss.runctx import load_context
    from zvss.stages import run_all
    workspace = out/'workspace'; shutil.copytree(ROOT/'03_virtual_screening/configs', workspace/'configs')
    receptors = workspace/'receptors'; receptors.mkdir()
    shutil.copyfile(files['receptor'], receptors/'USER_INPUT.pdb')
    atoms = [line for line in files['receptor'].read_text().splitlines() if line.startswith('ATOM  ')]
    if not atoms:
        raise ValueError('Receptor has no protein atom records')
    chains = sorted({line[21:22] for line in atoms})
    if len(chains) != 1:
        raise ValueError('Generic screening currently requires a single reviewed receptor chain')
    profile = yaml.safe_load((workspace/'configs/targets/generic.yaml').read_text())
    profile['structures'] = {'primary': {'pdb_id': 'USER_INPUT', 'chain': chains[0], 'role': 'Exact uploaded receptor'}}
    profile['structures_dir'] = str(receptors)
    (workspace/'target.yaml').write_text(yaml.safe_dump(profile))
    a.atomic(workspace/'host.json', {'paths': {'structures_dir': str(receptors)},
        'tools': {'vina': [cfg['vina']], 'vina_settings': {'engine_version': '1.2.7', 'exhaustiveness': cfg['exhaustiveness'],
            'num_modes': cfg['num_modes'], 'cpu': cfg['docking_cpu'], 'energy_range': cfg['energy_range']}},
        'executors': {'md_openmm_linux': {'mode': 'disabled', 'note': 'MD is a separately reviewed independent module'}}})
    ctx = load_context(root=str(workspace), run_id='screen', target_profile=str(workspace/'target.yaml'),
                       budget_config=str(workspace/'configs/budgets/full.yaml'), host_config=str(workspace/'host.json'))
    inputs = {'input_dirs': [str(files['poses'].parent)], 'manifest': None, 'recursive': False,
              'hints': {'receptor_file': str(files['receptor']), 'receptor_sha256': a.digest(files['receptor']), 'seed': p['seed']}}
    state = run_all(ctx, inputs, progress=lambda stage, info: print(stage, info.get('status'), flush=True))
    if state.get('errors'):
        raise RuntimeError('Screening reported errors; inspect workspace/runs/screen/run_state.json')
    return {'run_state': state, 'candidate_qualified': False, 'scope': 'Generic imported-pose assessment; missing evidence is not a pass'}


def execute(d, req, operation, out):
    from rdkit import Chem
    kind = req['module']; cfg = req['config']; p = req['parameters']
    files = {key: d/row['path'] for key, row in req['inputs'].items()}
    if kind == 'prepare':
        from src.structure import prepare_structure
        result = prepare_structure(files['complex'], out, structure_id='USER_INPUT',
            receptor_chain=p['receptor_chain'], peptide_chain=p['peptide_chain'],
            receptor_start=p['first_residue'], receptor_end=p['last_residue'], crop_distance_A=p['pocket_radius'],
            source_url='independent uploaded complex')
        return result['report']
    if kind in {'peptide2mol', 'reinvent4'}:
        from src.adapters.registry import backend
        model = a.read(ROOT/'01_molecule_generation/configs/project.json')['models'][kind]
        gen = ROOT/'01_molecule_generation'
        args = dict(output_dir=out/'model', budget=p['budget'], seed=p['seed'], repo=gen/'vendor'/('REINVENT4' if kind == 'reinvent4' else 'Peptide2Mol'))
        if kind == 'peptide2mol':
            args.update(pocket_pdb=files['pocket'], checkpoint=gen/model['checkpoint'],
                        checkpoint_sha256=model['checkpoint_sha256'], executable=cfg['peptide_python'])
        else:
            from yuanhe.runtime import runtime_env, resolve_python
            scaffold = Chem.MolFromSmiles(p['scaffold'])
            if scaffold is None or sum(atom.GetAtomicNum() == 0 for atom in scaffold.GetAtoms()) != 1:
                raise ValueError('REINVENT4 requires exactly one attachment atom')
            args.update(scaffold_smiles=p['scaffold'], checkpoint=gen/model['prior'],
                        checkpoint_sha256=model['prior_sha256'], executable=cfg['reinvent_executable'], device='cpu',
                        environment=runtime_env(resolve_python(cfg, 'reinvent_python'), for_overlay=True))
        result = backend(kind).run(**args)
        return {'report': result['report'], 'scope': 'Real pretrained sampling; candidates remain unqualified'}
    if kind in {'risk', 'availability'}:
        return preflight(files['library'], cfg, out, availability=kind == 'availability')
    if kind == 'docking_control':
        from .validation import docking_control
        return docking_control(d, cfg, {}, out, ROOT)
    if kind == 'docking':
        from sbvs2.docking.vina_official import run_vina_job
        from sbvs2.docking.box import Box
        from sbvs2.docking.protocol import Protocol
        molecule(files['ligand'])
        evidence = preflight(files['ligand'], cfg, out, availability=True)
        qualified = a.read(out/'availability/qualified.json', [])
        if len(qualified) != 1:
            return {'status': 'blocked_by_strict_preflight', 'screening': evidence, 'docked': 0, 'candidate_qualified': False}
        row = qualified[0]
        result = run_vina_job(cfg['vina'], files['receptor'], files['ligand'],
            Box(tuple(p['center']), tuple(p['size']), 'reviewed_independent_input', approved_reference_pose=False),
            Protocol.preset('standard', seeds=(p['seed'],)), p['seed'], out/'docking', cpu=cfg['docking_cpu'],
            qualification_evidence=row['qualification_evidence'])
        if result['status'] != 'ok':
            raise RuntimeError(result.get('error', 'Docking failed'))
        return {'docking': result, 'screening': evidence, 'scope': 'Strict preflight passed; downstream scientific qualification still required'}
    if kind == 'screening':
        return screening(d, cfg, p, files, out)
    if kind == 'interaction':
        from zlis.pipeline import analyze
        from zlis.tables import write_tables
        from zlis.diagram2d import write_diagrams
        from zlis.pymol_scene import write_scene_files
        rows = []
        for result in analyze(str(files['complex']), str(out/'analysis'), ligand=p['ligand'],
                              template=str(files['template']) if 'template' in files else None, allow_network=False):
            target = Path(result['out_dir'])
            outputs = {**write_tables(result, target), **write_diagrams(result, target), **write_scene_files(result, target)}
            rows.append({'ligand': result['ligand']['label'], 'n_events': len(result['events']),
                         'chemistry': result['chemistry']['status'], 'files': outputs})
        if not rows:
            raise ValueError('No matching ligand in the uploaded complex')
        return {'ligands': rows, 'scope': 'Static interactions; no affinity or activity inference'}
    if kind == 'md_tpr':
        from .prepared_md import run
        return run(d, cfg, p, files['production'])
    if kind == 'md':
        from . import membrane
        action = operation['action']
        if action == 'build':
            _, smiles = molecule(files['ligand'])
            ident = 'INPUT_'+a.digest(files['ligand'])[:20]
            row = {'molecule_id': ident, 'compound_id': ident, 'pose_id': ident+'_pose',
                   'qualification': 'independent', 'smiles': smiles}
            parameters = {**p, 'exploratory_ack': True, 'orientation_reviewed': True,
                          'replica_seed': p['seed'], 'output_interval_ps': min(10, p['production_ns']*1000/10)}
            return membrane.build_prepared(d, cfg, parameters, out, row, files['ligand'], files['receptor'])
        return getattr(membrane, action)(d, cfg, {'system_id': operation['system_id'],
                                                'system_reviewed': operation['system_reviewed']}, out)
    raise ValueError('No independent adapter for module')
