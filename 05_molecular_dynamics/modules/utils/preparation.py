from __future__ import annotations
from .core import WorkflowError, relative_path, atomic_json, write_pdb

def prepare_protein(config, paths):
    from pdbfixer import PDBFixer
    from openmm import Platform
    fixer = PDBFixer(filename=str(relative_path(config["protein"]["file"], must_exist=True)), platform=Platform.getPlatformByName("CPU"))
    settings = config["protein"]
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        raise WorkflowError("Nonstandard protein residues require an explicit parameterization/preparation protocol.")
    allowed = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
               "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","HOH","WAT"}
    unknown = [r for r in fixer.topology.residues() if r.name not in allowed]
    if unknown and not settings["remove_heterogens"]:
        raise WorkflowError("Protein contains cofactors/metals/heterogens; prepare explicitly or set remove_heterogens knowingly.")
    if settings["remove_heterogens"]:
        fixer.removeHeterogens(keepWater=not settings["remove_water"])
    elif settings["remove_water"]:
        from openmm.app import Modeller
        model = Modeller(fixer.topology, fixer.positions)
        model.delete([r for r in model.topology.residues() if r.name in ("HOH","WAT")])
        fixer.topology, fixer.positions = model.topology, model.positions
    fixer.findMissingResidues()
    missing_residue_count = sum(len(v) for v in fixer.missingResidues.values())
    if missing_residue_count and not settings["allow_missing_residues"]:
        raise WorkflowError("Unresolved sequence gaps detected. Model/review missing residues before production.")
    fixer.findMissingAtoms()
    missing_atom_count = sum(len(v) for v in fixer.missingAtoms.values())
    fixer.addMissingAtoms(seed=config["MD"]["seed"])
    fixer.addMissingHydrogens(pH=settings["ph"])
    target = paths["md"] / "prepared"
    target.mkdir(parents=True, exist_ok=True)
    write_pdb(fixer.topology, fixer.positions, target / "protein.pdb")
    atomic_json(target / "protein_preparation.json", {
        "ph": settings["ph"], "preparation_platform": fixer.platform.getName(), "added_missing_heavy_atoms": missing_atom_count,
        "modeled_missing_residues": missing_residue_count,
        "residues": sum(1 for _ in fixer.topology.residues()),
        "note": "Default protonation is not a site-specific pKa calculation."
    })
