"""Pinned official CPU Vina adapter; Meeko retains ligand chemistry and atom mapping."""
from __future__ import annotations
import math, os, re, shutil, subprocess, time
from pathlib import Path
from ..hashing import atomic_write_json, sha256_file, hash_obj
from .execute import run_argv, Cancelled
from .run_vina import parse_vina_pdbqt
from .pdbqt import receptor_to_pdbqt, ligand_to_pdbqt

VERSION = '1.2.7'

def vina_version(executable):
    p = subprocess.run([str(executable), '--version'], capture_output=True, text=True, timeout=30)
    match = re.search(r'AutoDock Vina\s+v?(\d+\.\d+\.\d+)', p.stdout + p.stderr)
    if p.returncode or not match or match[1] != VERSION:
        raise ValueError(f'Official AutoDock Vina {VERSION} required: {(p.stdout+p.stderr)[:200]}')
    return match[1]

def validate_resources(cpu, workers=1, total_threads=None):
    total = int(total_threads or os.cpu_count() or 1)
    if any(type(x) is not int or x < 1 for x in (cpu, workers, total)) or cpu * workers > total:
        raise ValueError('Vina CPU budget requires workers * cpu <= total_threads and positive integers')
    return total

def export_sdf(pdbqt, source_mol, output, rows):
    from rdkit import Chem
    from meeko import PDBQTMolecule, RDKitMolCreate
    text = Path(pdbqt).read_text()
    if 'REMARK SMILES IDX' not in text or 'REMARK SMILES ' not in text:
        raise ValueError('Missing Meeko chemical template / atom mapping; coordinate-only reconstruction forbidden')
    molecules = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text, skip_typing=True))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError('Expected exactly one mapped ligand')
    mol = molecules[0]
    canonical = lambda m: Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=True)
    if canonical(mol) != canonical(source_mol) or mol.GetNumConformers() != len(rows):
        raise ValueError('Export chemistry/stereochemistry or pose count differs from trusted input')
    with Chem.SDWriter(str(output)) as writer:
        for conf, row in zip(mol.GetConformers(), rows):
            mol.SetProp('docking_engine', 'vina'); mol.SetProp('engine_version', VERSION)
            mol.SetProp('scoring_function', 'vina'); mol.SetProp('score_unit', 'kcal/mol')
            mol.SetProp('vina_affinity_kcal_mol', str(row['vina_affinity']))
            mol.SetProp('docking_mode', str(row['pose_id']))
            mol.SetProp('cnn_status', 'not_applicable')
            writer.write(mol, confId=conf.GetId())

def run_vina_job(exe, receptor, ligand, box, proto, seed, job_dir, *, cpu=4,
                 timeout=3600, cancel=None, cancel_file=None, meta=None, **kwargs):
    from rdkit import Chem
    validate_resources(cpu)
    if kwargs.get('operation','redock') != 'redock':raise ValueError('Only explicit Vina redock jobs supported')
    if proto.engine_version!=VERSION or proto.scoring!='vina':raise ValueError('Pinned Vina version/scoring required')
    if len(box.center)!=3 or len(box.size)!=3 or any(v in box.validate() for v in ['NON_FINITE','NON_POSITIVE_SIZE']):raise ValueError('Invalid angstrom search box')
    if proto.exhaustiveness<1 or proto.num_modes<1 or not math.isfinite(proto.energy_range) or proto.energy_range<=0:raise ValueError('Invalid Vina search parameters')
    if proto.engine != 'vina': raise ValueError('New production jobs require engine=vina')
    if kwargs.get('device', 'cpu') != 'cpu': raise ValueError('Official Vina is a CPU engine')
    jd = Path(job_dir).resolve(); jd.mkdir(parents=True, exist_ok=False)
    rec = dict(meta or {}, engine='vina', engine_version=VERSION, scoring_function=proto.scoring, box_units='angstrom', box_center=list(box.center), box_size=list(box.size),
               score_unit='kcal/mol', cnn_status='not_applicable', seed=int(seed), cpu=cpu,
               protocol_hash=proto.hash(), box_hash=box.hash(), started=time.time(), status='running')
    try:
        vina_version(exe)
        mols = list(Chem.SDMolSupplier(str(ligand), removeHs=False))
        if len(mols) != 1 or mols[0] is None: raise ValueError('One trusted SDF molecule required')
        mol = Chem.AddHs(mols[0], addCoords=True)
        if not mol.GetNumConformers(): raise ValueError('Ligand coordinates missing')
        text, mapping = ligand_to_pdbqt(mol, 0)
        lig = jd/'ligand.pdbqt'; lig.write_text(text)
        shutil.copyfile(ligand, jd/'input_ligand.sdf')
        receptor = Path(receptor).resolve()
        prepared = receptor if receptor.suffix.lower()=='.pdbqt' else receptor_to_pdbqt(receptor, jd/'receptor')
        rec.update(receptor_sha256=sha256_file(receptor), ligand_sha256=sha256_file(ligand), mapping=mapping)
        rec['cache_key'] = hash_obj({k:rec[k] for k in ('engine','engine_version','scoring_function','seed','cpu','protocol_hash','box_hash','receptor_sha256','ligand_sha256')})
        out = jd/'poses.pdbqt'
        argv = [str(exe), *proto.vina_args(str(prepared), str(lig), box, seed, str(out)), '--cpu', str(cpu)]
        rec['audit'] = run_argv(argv, jd, timeout=timeout, cancel=cancel, cancel_file=cancel_file)
        if rec['audit']['returncode'] != 0: raise RuntimeError('Vina failed; inspect stderr.log')
        rows = parse_vina_pdbqt(out.read_text())
        if not rows or any(not math.isfinite(r['vina_affinity']) for r in rows): raise ValueError('EMPTY_OUTPUT: invalid Vina output')
        export_sdf(out, mol, jd/'poses.sdf', rows)
        rec['poses'] = [{**r, 'vina_affinity_kcal_mol':r['vina_affinity'], 'empirical_score_kcal_mol':r['vina_affinity'],
                         'engine':'vina', 'cnn_pose_score':None, 'cnn_affinity_pK':None} for r in rows]
        rec.update(status='ok', output_sha256=sha256_file(jd/'poses.sdf'))
    except Cancelled: rec['status']='cancelled'
    except TimeoutError as exc: rec.update(status='timeout',error=str(exc))
    except Exception as exc: rec.update(status='failed', error=f'{type(exc).__name__}: {exc}')
    rec['finished']=time.time(); atomic_write_json(jd/'job_record.json',rec)
    return rec
