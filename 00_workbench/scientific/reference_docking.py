"""A fixed public reference experiment, never a candidate qualification path."""
import hashlib
import time
from pathlib import Path
from sbvs2.hashing import atomic_write_json, sha256_file
from sbvs2.docking.vina_official import (
    vina_version, validate_resources, ligand_to_pdbqt, receptor_to_pdbqt,
    run_argv, parse_vina_pdbqt, export_sdf,
)

REFERENCE = {
    'ligand': 'c9c8bf1f3f4cd2d4aaa0fdd440369fa96502b4a702b761d031d5663f3cde3268',
    'receptor': '248293b266fdcf31ed8b32ceb175589004682a1c9aa9db746b8429f799801600',
}


def reference_inputs(receptor, ligand):
    """Only the bundled 1DRF/FOL pair is admitted, independent of file names."""
    proof = {'schema': 'yuanhe.reference-control/1', 'reference': '1DRF/FOL',
             'scope': 'crystallographic pose recovery', 'candidate_qualified': False,
             'scientific_claim': 'No toxicity, synthesis, affinity or activity qualification',
             'inputs': {}}
    for key, path in [('receptor', receptor), ('ligand', ligand)]:
        path = Path(path)
        canonical = hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        if canonical != REFERENCE[key]:
            raise ValueError('Reference control accepts only the pinned public 1DRF/FOL fixture: ' + key)
        proof['inputs'][key] = {'sha256': sha256_file(path), 'canonical_sha256': canonical}
    return proof


def run_reference_job(exe, receptor, ligand, box, proto, seed, job_dir, *, cpu=4, timeout=3600):
    from rdkit import Chem
    evidence = reference_inputs(receptor, ligand)
    validate_resources(cpu)
    if proto.engine != 'vina' or proto.engine_version != '1.2.7' or proto.scoring != 'vina':
        raise ValueError('Pinned official Vina reference protocol required')
    if any(x in box.validate() for x in ['NON_FINITE', 'NON_POSITIVE_SIZE']):
        raise ValueError('Invalid reference box')
    job = Path(job_dir)
    job.mkdir(parents=True, exist_ok=False)
    atomic_write_json(job/'reference_control.json', evidence)
    record = {'status': 'running', 'scope': evidence['scope'], 'candidate_qualified': False,
              'reference_evidence_sha256': sha256_file(job/'reference_control.json'),
              'seed': seed, 'cpu': cpu, 'protocol_hash': proto.hash(), 'box_hash': box.hash(),
              'started': time.time()}
    try:
        vina_version(exe)
        mol = Chem.AddHs(Chem.SDMolSupplier(str(ligand), removeHs=False)[0], addCoords=True)
        text, mapping = ligand_to_pdbqt(mol, 0)
        (job/'ligand.pdbqt').write_text(text)
        prepared = receptor_to_pdbqt(Path(receptor), job/'receptor')
        argv = [str(exe), *proto.vina_args(str(prepared), str(job/'ligand.pdbqt'),
                box, seed, str(job/'poses.pdbqt')), '--cpu', str(cpu)]
        record['audit'] = run_argv(argv, job, timeout=timeout)
        if record['audit']['returncode']:
            raise RuntimeError('Vina reference execution failed; inspect stderr.log')
        rows = parse_vina_pdbqt((job/'poses.pdbqt').read_text())
        if not rows:
            raise ValueError('Reference experiment produced no poses')
        export_sdf(job/'poses.pdbqt', mol, job/'poses.sdf', rows)
        record.update(status='ok', poses=rows, mapping=mapping, output_sha256=sha256_file(job/'poses.sdf'))
    except Exception as exc:
        record.update(status='failed', error=f'{type(exc).__name__}: {exc}')
    record['finished'] = time.time()
    atomic_write_json(job/'job_record.json', record)
    # Deliberately no qualification sidecar: candidate engines still require both gates.
    return record
