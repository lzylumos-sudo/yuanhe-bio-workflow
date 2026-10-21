"""One mandatory gate used by workbench, CLI, scheduler and redocking.
Hashes detect changed evidence; validity is recomputed, not trusted as a boolean."""
from pathlib import Path
from .common import read, verify_seal, file_hash, canonical, sdf_smiles, hash_obj, identity
from .toxicity_substructure_filter import RiskRules
from .availability_synthesis_filter import BUY_LABELS, validate_route, assets
from .catalog import Catalog

class QualificationError(ValueError): pass

def validate_evidence(evidence, smiles, *, input_file=None):
    if not evidence:raise QualificationError('predocking_evidence_required')
    proof=read(evidence) if isinstance(evidence,(str,Path)) else evidence
    verify_seal(proof)
    if proof.get('schema')!='yuanhe-predocking/1':raise QualificationError('invalid_evidence_schema')
    current=identity(smiles)
    if any(proof['identity'].get(k)!=current[k] for k in ('smiles','structure_sha256','compound_id','state_id')):
        raise QualificationError('structure_identity_mismatch')
    source=Path(proof['source_sdf'])
    if not source.is_file() or file_hash(source)!=proof['source_sdf_sha256'] or sdf_smiles(source)!=current['smiles']:
        raise QualificationError('screened_input_file_changed')
    if input_file and sdf_smiles(input_file)!=current['smiles']:raise QualificationError('docking_input_structure_mismatch')
    from .config import validate
    validate(proof['config'])
    cfg=proof['config']; risk=proof['risk']; availability=proof['availability']
    verify_seal(risk);verify_seal(availability)
    rules=RiskRules(cfg['rules_path'],cfg['risk_policies'])
    if risk!=rules.screen(current['smiles'],risk.get('exceptions',[])) or risk['decision']!='pass':
        raise QualificationError('risk_incomplete_blocked_or_changed')
    if availability.get('smiles')!=current['smiles']:raise QualificationError('availability_identity_mismatch')
    catalog=Catalog(cfg['catalogs'],cfg['evidence_max_age_days'])
    stock=Catalog(cfg['stock_catalogs'],cfg['evidence_max_age_days'])
    if catalog.versions!=availability['catalog_versions'] or stock.versions!=availability['stock_versions']:
        raise QualificationError('catalog_or_stock_version_changed')
    valid=catalog.valid(current['smiles'])
    purchase=any(any(r==saved for saved in availability['supplier_records']) for r in valid)
    route_ok=False
    if availability['routes']:
        manifest,error=assets(cfg['retro'])
        if error or manifest!=availability['asset_manifest']:raise QualificationError('retrosynthesis_assets_missing_or_changed')
        raw=availability.get('raw_search') or {}
        p=Path(raw.get('search_artifact',''))
        if not p.is_file():raise QualificationError('raw_search_missing')
        original=read(p)
        if original.get('raw_sha256')!=hash_obj(original.get('raw')):raise QualificationError('raw_search_changed')
        if original.get('assets')!=manifest:raise QualificationError('search_asset_mismatch')
        trees=original['raw'].get('routes',[])
        if canonical(original['raw'].get('target_smiles',''))!=current['smiles']:raise QualificationError('raw_search_target_mismatch')
        for route in availability['routes']:
            if not any(hash_obj(t)==route['route_sha256'] for t in trees):raise QualificationError('route_not_in_raw_search')
            checked=validate_route(route['tree'],current['smiles'],stock,route_review=route.get('review'))
            if checked['status']=='complete' and checked==route:route_ok=True
    if not purchase and not route_ok:raise QualificationError('no_current_supply_or_complete_route_evidence')
    return proof

def require_ligand(ligand, evidence=None):
    ligand=Path(ligand)
    # Portable per-ligand sidecars are explicit evidence, never a bypass flag.
    path=evidence or ligand.with_suffix(ligand.suffix+'.qualification.json')
    if not Path(path).is_file():raise QualificationError('predocking_evidence_required: '+str(path))
    return validate_evidence(path,sdf_smiles(ligand),input_file=ligand)
