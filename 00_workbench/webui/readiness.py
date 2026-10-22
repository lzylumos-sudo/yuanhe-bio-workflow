"""Standard-library deployment inventory shared by the UI, CLI and start gate.

Checks imports, tools, model hashes and configured evidence before expensive work.
It does not generate molecules, load retrosynthesis models, or run target MD.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


def runtime_env(executable):
    env = dict(os.environ)
    exe = shutil.which(str(executable)) or str(executable)
    bindir = Path(exe).absolute().parent
    env['PATH'] = str(bindir) + os.pathsep + env.get('PATH', '')
    env.pop('PYTHONHOME', None)
    env.pop('PYTHONPATH', None)
    env['PYTHONNOUSERSITE'] = '1'
    if (bindir.parent / 'dat/leap').is_dir():
        env['AMBERHOME'] = str(bindir.parent)
    return env


def probe(argv, root, timeout=90):
    try:
        if not argv[0]:
            raise ValueError('Executable is not configured')
        p = subprocess.run(list(map(str, argv)), cwd=root, env=runtime_env(argv[0]),
                           capture_output=True, text=True, timeout=timeout)
        return {'ok': p.returncode == 0, 'returncode': p.returncode,
                'output': (p.stdout + p.stderr)[-6000:], 'executable': str(argv[0])}
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        return {'ok': False, 'error': str(exc), 'executable': str(argv[0])}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def ui_check(executable, root):
    return probe([executable, '-c',
        'import sys,rdkit,numpy; from rdkit.Chem.Draw import rdMolDraw2D; '
        'print(sys.executable, rdkit.__version__)'], root)


def collect(cfg, root, ui_python=None, publish=None):
    root = Path(root).resolve()
    checks = {}
    result = {'schema': 'yuanhe.readiness/1', 'started': time.time(), 'checks': checks,
              'status': 'running', 'all_runtime_checks_pass': False,
              'scientific_validation': 'not_run',
              'limitations': ['Imports and file checks do not validate model inference, '
                             'a target membrane system, GPU GROMACS execution or 500 ns completion.',
                             'Candidate qualification, model loading, reviews and dynamic resources '
                             'remain checked during execution. Zero qualified candidates is possible.']}

    def add(name, value, remedy, required=True):
        checks[name] = {**value, 'required': required, 'remedy': remedy}
        if publish:
            publish(result)

    def py(name, executable, code, remedy, required=True):
        add(name, probe([executable, '-c', code], root), remedy, required)

    add('platform', {'ok': platform.system() == 'Linux' and platform.machine() in ('x86_64', 'AMD64'),
                     'platform': platform.platform()}, 'Deploy on Linux x86_64 / WSL2.')
    add('ui_interpreter', ui_check(ui_python or sys.executable, root),
        'Restart the web service with the configured chemistry Python via start_linux.sh.')
    chem = cfg.get('chem_python', '')
    screen = cfg.get('screen_python') or chem
    md = cfg.get('md_python', '')
    py('chemistry', chem,
       'import rdkit,numpy,yaml,gemmi,scipy,pandas,filelock,meeko; '
       'from rdkit.Chem import FilterCatalog,rdMolDescriptors; print(rdkit.__version__)',
       'Install the complete workbench runtime; filelock and Meeko are required too.')
    py('screening', screen,
       'import MDAnalysis,prolif,posebusters,fastapi,scipy,pandas,yaml,gemmi; print("screening imports OK")',
       'Install screening dependencies in screen_python, not just chem_python.')
    py('preparation', cfg.get('prep_python', ''),
       'import pdbfixer,openmm; from openmm.app import ForceField; '
       'ForceField("amber14-all.xml","amber14/tip3p.xml"); print(openmm.__version__)',
       'Install PDBFixer and OpenMM in prep_python.')
    py('membrane_and_analysis', md,
       'import openmm,parmed,MDAnalysis,psi4,resp,rdkit,numpy,matplotlib,yaml; '
       'from openff.toolkit import Molecule; from openff.units import unit; '
       'from openmmforcefields.generators import GAFFTemplateGenerator; '
       'from openmm.app import ForceField; ForceField("amber14-all.xml","amber14/tip3p.xml"); '
       'import shutil; tools=["antechamber","parmchk2","tleap","MMPBSA.py","sander","cpptraj"]; '
       'missing=[x for x in tools if not shutil.which(x)]; assert not missing, missing; print(openmm.__version__)',
       'Install OpenMM/OpenFF, AmberTools, Psi4, RESP and analysis packages in md_python; '
       'its bin directory must contain the AmberTools analysis executables.')
    py('openmm_cuda', md,
       'import openmm as m; s=m.System();s.addParticle(1);i=m.VerletIntegrator(0.001); '
       'c=m.Context(s,i,m.Platform.getPlatformByName("CUDA")); '
       'c.setPositions([[0,0,0]]);i.step(1);print(c.getPlatform().getName())',
       'Provide a compatible NVIDIA driver and CUDA OpenMM runtime; no CPU fallback.')
    vina = probe([cfg.get('vina', ''), '--version'], root)
    vina['ok'] = vina['ok'] and 'AutoDock Vina v1.2.7' in vina.get('output', '')
    add('vina', vina, 'Configure the official Vina 1.2.7 executable.')
    gmx = probe([cfg.get('gmx', ''), '--version'], root)
    # This full workflow requests a GPU. A CPU-only executable is not sufficient.
    gmx['ok'] = gmx['ok'] and any('GPU support:' in line and 'CUDA' in line
                                               for line in gmx.get('output', '').splitlines())
    add('gromacs_cuda_build', gmx, 'Configure a CUDA-enabled GROMACS build; then validate actual mdrun on this host.')
    smi = shutil.which('nvidia-smi') or '/usr/lib/wsl/lib/nvidia-smi'
    add('nvidia', probe([smi, '--query-gpu=name,driver_version,memory.total', '--format=csv,noheader'], root),
        'Make the administrator-assigned NVIDIA GPU available to Linux / the container.')

    if not cfg.get('imported_candidates'):
        py('peptide_cuda', cfg.get('peptide_python', ''),
           'import torch,torch_geometric,torch_sparse,torch_scatter,torch_cluster,lightning,hydra,rdkit,rootutils,lmdb; '
           'assert torch.cuda.is_available(), "CUDA unavailable"; '
           'x=torch.tensor([[0.,0.,0.],[1.,0.,0.]],device="cuda"); '
           'print(torch_cluster.radius_graph(x,2)); '
           'print(torch_scatter.scatter_add(torch.ones(2,device="cuda"),torch.tensor([0,0],device="cuda"))); '
           'torch.cuda.synchronize();print(torch.cuda.get_device_name())',
           'Install the pinned Peptide2Mol Python/PyTorch and matching CUDA PyG wheels.')
        add('reinvent_cli', probe([cfg.get('reinvent_executable', ''), '--help'], root),
            'Install the pinned, patched REINVENT4 source into its dedicated environment.')
        reinvent_exe = shutil.which(cfg.get('reinvent_executable', '')) or cfg.get('reinvent_executable', '')
        reinvent_python = str(Path(reinvent_exe).absolute().parent / 'python') if reinvent_exe else ''
        if (cfg.get('research') or {}).get('reinvent_device', 'auto') != 'cpu':
            py('reinvent_cuda', reinvent_python,
               'import torch;assert torch.cuda.is_available(),"CUDA unavailable";'
               'x=torch.ones((8,8),device="cuda");print((x@x).sum().item());torch.cuda.synchronize()',
               'Validate the REINVENT4 Torch build on this GPU, or explicitly select reinvent_device=cpu before task creation.')
        try:
            models = json.loads((root / 'deployment/model_dependencies.json').read_text())
            for name, pathkey, hashkey in [('peptide2mol', 'checkpoint', 'checkpoint_sha256'),
                                         ('reinvent4', 'prior', 'prior_sha256')]:
                spec = models[name]; path = root / '01_molecule_generation' / spec[pathkey]
                actual = digest(path) if path.is_file() else None
                add(name + '_weights', {'ok': actual == spec[hashkey], 'path': str(path), 'sha256': actual},
                    'Restore the licensed source and exact weights in deployment/model_dependencies.json.')
            for row in json.loads((root / 'deployment/vendor_patches/manifest.json').read_text()):
                path = root / '01_molecule_generation/vendor' / row['vendor'] / row['file']
                add(row['vendor'] + '_patch', {'ok': path.is_file() and digest(path) == row['after_sha256'], 'path': str(path)},
                    'Apply deployment/apply_vendor_patches.py before installing the model environment.')
                if row['vendor'] == 'REINVENT4':
                    py('reinvent_installed_patch', reinvent_python,
                       'import hashlib,importlib.util;from pathlib import Path;'
                       'p=Path(importlib.util.find_spec("reinvent.runmodes.samplers.libinvent").origin);'
                       'print(p);assert hashlib.sha256(p.read_bytes()).hexdigest()==' + repr(row['after_sha256']),
                       'The installed REINVENT4 sampler must contain the recorded patch, not an older site-packages copy.')
        except (OSError, ValueError, KeyError) as exc:
            add('model_manifest', {'ok': False, 'error': str(exc)}, 'Restore the complete repository and model manifests.')

    rc = cfg.get('research') or {}
    pc = rc.get('predocking') or cfg.get('predocking') or {}
    try:
        # Keep metadata validation independent of the RDKit installation.
        from yuanhe_preflight.config import validate, absolute_paths
        validate(pc)
        pc = absolute_paths(pc, root)
        add('predocking_config', {'ok': True}, 'Configure predocking according to docs/PREDOCKING_PREFLIGHT.md.')
    except (ValueError, KeyError, TypeError) as exc:
        add('predocking_config', {'ok': False, 'error': str(exc)}, 'Correct the predocking configuration before creating the research task.')

    def catalog_check(key):
        descriptors = pc.get(key, [])
        code = ('import sys;sys.path.insert(0,' + repr(str(root / '02_docking/src')) + ');'
                'from yuanhe_preflight.catalog import Catalog;c=Catalog(' + repr(descriptors) + ',' +
                repr(pc.get('evidence_max_age_days', 30)) + ');'
                'valid=sum(c.evaluate(r) in ("in_stock","made_to_order") for r in c.records);'
                'print({"rows":len(c.records),"valid_supply_records":valid,"parse_errors":c.errors[:10]});'
                'assert valid>0,"No current, hash-verified supply evidence"')
        return probe([chem, '-c', code], root) if descriptors else {'ok': False, 'error': 'No catalogs configured'}

    purchase = catalog_check('catalogs')
    retro = pc.get('retro', {})
    requested_retro = any(retro.get(k) for k in ('model_path', 'templates_path', 'stock_path'))
    retro_ok = False
    if requested_retro:
        missing = [k for k in ('model_path', 'templates_path', 'stock_path', 'model_version',
                    'templates_version', 'stock_version', 'model_license', 'templates_license', 'stock_license') if not retro.get(k)]
        missing += [k for k in ('model_path', 'templates_path', 'stock_path')
                    if retro.get(k) and not (root / retro[k]).is_file()]
        dep = probe([retro.get('python', ''), '-c',
                     'import importlib.metadata as m;from aizynthfinder.aizynthfinder import AiZynthFinder;'
                     'assert m.version("aizynthfinder")=="4.4.1";print(m.version("aizynthfinder"))'], root)
        stock = catalog_check('stock_catalogs')
        retro_ok = not missing and dep['ok'] and stock['ok']
        add('retrosynthesis', {'ok': retro_ok, 'missing': missing, 'dependency': dep, 'stock': stock},
            'Install AiZynthFinder 4.4.1 separately; supply model, templates, stock, licenses/versions and current terminal supply evidence.')
    add('purchase_catalogs', purchase, 'Supply current exact-structure catalog records and hashed supply documents.', required=not retro_ok)
    add('acquisition_path', {'ok': purchase['ok'] or retro_ok},
        'At least one real acquisition route must be configured before generation. '
        'An empty supplier/model configuration cannot advance to docking. Passing this inventory does not qualify any molecule.')
    rules = pc.get('rules_path', '')
    py('risk_rules', chem,
       'import sys;sys.path.insert(0,' + repr(str(root / '02_docking/src')) + ');'
       'from yuanhe_preflight.toxicity_substructure_filter import RiskRules;'
       'r=RiskRules(' + repr(str(root / rules)) + ');print(r.version)',
       'Restore the packaged risk rules and install a compatible RDKit FilterCatalog.')
    structure = cfg.get('structure', '')
    add('input_structure', {'ok': bool(structure) and (root / structure).is_file(), 'path': str(root / structure)},
        'Provide the configured project-local PDB; review chains and construct before computation.')
    reserve = rc.get('disk_reserve_gb', 50)
    storage = root / '00_workbench/workflow_runs'
    free = shutil.disk_usage(storage if storage.exists() else root).free
    add('disk_reserve', {'ok': isinstance(reserve, (int, float)) and free >= reserve * 1024**3,
                         'free_gb': free / 1024**3, 'reserve_gb': reserve},
        'Provide space for environments, candidates and trajectories in addition to the configured reserve.')
    add('project_writable', {'ok': os.access(root, os.W_OK) and os.access(root / '00_workbench', os.W_OK)},
        'Run the service as the owner of its run/output directories, without root privileges.')
    result['all_runtime_checks_pass'] = all(x['ok'] for x in checks.values() if x['required'])
    result['status'] = 'passed' if result['all_runtime_checks_pass'] else 'blocked'
    result['ended'] = time.time()
    if publish:
        publish(result)
    return result


def require_ready(cfg, root, report_path):
    from .workflow import atomic
    result = collect(cfg, root)
    atomic(report_path, result)
    if not result['all_runtime_checks_pass']:
        failed = ', '.join(k for k, v in result['checks'].items() if v['required'] and not v['ok'])
        raise RuntimeError('Deployment readiness blocked: ' + failed +
                           '. Open Environment and diagnostics; report: ' + str(report_path))
    return result
