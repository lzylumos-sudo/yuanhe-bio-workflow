"""Prepare supplied protein complexes locally; input heavy-atom coordinates are authoritative."""
from __future__ import annotations
import copy
from io import StringIO
from .core import WorkflowError, relative_path, rel, atomic_json, write_pdb

AMINO = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
         "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
WATER = {"HOH","WAT"}

def residue_key(line):
    return (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), line[17:20].strip())

def split_complex(text, selector=None):
    """Return protein PDB, optional ligand PDB and selector; reject unsupported chemistry."""
    import math
    lines = text.splitlines()
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        raise WorkflowError("Multiple PDB models: supply one selected model.")
    records = [line for line in lines if line.startswith(("ATOM  ","HETATM"))]
    if not records:
        raise WorkflowError("The complex PDB contains no atoms.")
    if any(line[16:17].strip() for line in records):
        raise WorkflowError("Alternate locations must be resolved before preparation.")
    if len({line[6:11].strip() for line in records}) != len(records):
        raise WorkflowError("Duplicate atom serials in the complex PDB.")
    if any(not all(math.isfinite(float(line[a:b])) for a,b in ((30,38),(38,46),(46,54))) for line in records):
        raise WorkflowError("Non-finite complex coordinates.")
    keys = {residue_key(line) for line in records}
    selected = None
    if selector is not None:
        matches = [key for key in keys if key[3] == selector["resname"]
                   and (selector.get("chain") is None or key[0] == str(selector["chain"]))
                   and (selector.get("resid") is None or key[1] == str(selector["resid"]))
                   and (selector.get("insertion_code") is None or key[2] == str(selector["insertion_code"]))]
        if len(matches) != 1:
            raise WorkflowError("Ligand selector must identify exactly one residue; set resname, chain and resid.")
        selected = matches[0]
        if selected[3] in AMINO | WATER:
            raise WorkflowError("The selected ligand cannot be a standard protein or water residue.")
    # LINK records can describe covalent or coordination bonds. Never discard a
    # linkage to the ligand and subsequently treat the ligand as noncovalent.
    for line in lines:
        if line.startswith("LINK  "):
            endpoints = (
                (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), line[17:20].strip()),
                (line[51:52].strip(), line[52:56].strip(), line[56:57].strip(), line[47:50].strip()),
            )
            if selected is not None and selected in endpoints and endpoints[0] != endpoints[1]:
                raise WorkflowError("Cross-residue covalent/coordination LINK involving the ligand is unsupported.")
    if any(key != selected and key[3] not in AMINO | WATER for key in keys):
        raise WorkflowError("Additional ligands, metals or cofactors need explicit parameterization; supply a supported complex.")
    if not any(key[3] in AMINO for key in keys):
        raise WorkflowError("No standard protein residues found.")
    ligand_records = [line for line in records if residue_key(line) == selected]
    if len({line[12:16].strip() for line in ligand_records}) != len(ligand_records):
        raise WorkflowError("Ligand atom names must be unique within the selected residue.")
    ligand_serials = {int(line[6:11]) for line in ligand_records}
    remaining_serials = {int(line[6:11]) for line in records}-ligand_serials
    ligand_connections, protein_connections = [], []
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        serials = [int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
        if not serials:
            continue
        if serials[0] in ligand_serials:
            if any(i in remaining_serials for i in serials[1:]):
                raise WorkflowError("Cross-residue covalent ligand bonds are unsupported.")
            ligand_connections.append(line)
        elif serials[0] in remaining_serials:
            if any(i in ligand_serials for i in serials[1:]):
                raise WorkflowError("Cross-residue covalent ligand bonds are unsupported.")
            protein_connections.append(line)
    kept = []
    for line in lines:
        if line.startswith(("ATOM  ","HETATM")):
            if residue_key(line) != selected:
                kept.append(line)
        elif line.startswith(("SEQRES","SSBOND","CRYST1","TER")):
            kept.append(line)
    protein_text = "\n".join(kept+protein_connections+["END",""])
    ligand_text = "\n".join(ligand_records+ligand_connections+["END",""])
    return protein_text, ligand_text, selected

def connectivity_only(molecule):
    from rdkit import Chem
    result = Chem.Mol(molecule)
    for atom in result.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
        atom.SetAtomMapNum(0)
    for bond in result.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    result.UpdatePropertyCache(strict=False)
    return result

def ligand_at_complex_coordinates(template, observed, atom_mapping=None):
    """Map SDF chemistry to PDB coordinates without fitting, embedding or minimization."""
    import numpy as np
    from rdkit import Chem
    if template is None or observed is None:
        raise WorkflowError("Invalid ligand chemistry or complex ligand coordinates.")
    template = Chem.Mol(template)
    for i,atom in enumerate(template.GetAtoms()):
        atom.SetIntProp("_sdf_index",i+1)
    template = Chem.RemoveHs(template)
    # PDB-derived chemistry is discarded; only elements, connectivity and coordinates are used.
    observed = Chem.RemoveHs(observed, sanitize=False)
    observed.UpdatePropertyCache(strict=False)
    if len(Chem.GetMolFrags(template)) != 1 or len(Chem.GetMolFrags(observed)) != 1:
        raise WorkflowError("The selected ligand and template must each be a connected molecule.")
    if template.GetNumAtoms() != observed.GetNumAtoms():
        raise WorkflowError("Complex ligand and SDF heavy-atom counts differ.")
    query, target = connectivity_only(template), connectivity_only(observed)
    if query.GetNumBonds() != target.GetNumBonds():
        raise WorkflowError("Complex ligand and SDF connectivity differ (extra or missing bonds).")
    if atom_mapping:
        by_name = {}
        for i, atom in enumerate(observed.GetAtoms()):
            info = atom.GetPDBResidueInfo()
            if info is None:
                raise WorkflowError("Explicit mapping requires PDB atom names.")
            by_name[info.GetName().strip()] = i
        source_to_name = {value:key for key,value in atom_mapping.items()}
        expected = {a.GetIntProp("_sdf_index") for a in template.GetAtoms()}
        if set(source_to_name) != expected or set(atom_mapping) != set(by_name):
            raise WorkflowError("atom_mapping must cover every heavy atom exactly once (SDF atom numbers are 1-based).")
        mapping = tuple(by_name[source_to_name[a.GetIntProp("_sdf_index")]] for a in template.GetAtoms())
        if len(set(mapping)) != template.GetNumAtoms():
            raise WorkflowError("Atom mapping is not one-to-one.")
        mapped_bonds = {tuple(sorted((mapping[b.GetBeginAtomIdx()],mapping[b.GetEndAtomIdx()]))) for b in query.GetBonds()}
        observed_bonds = {tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in target.GetBonds()}
        if mapped_bonds != observed_bonds or any(query.GetAtomWithIdx(i).GetAtomicNum()!=target.GetAtomWithIdx(j).GetAtomicNum() for i,j in enumerate(mapping)):
            raise WorkflowError("Explicit atom mapping disagrees with the observed elements/connectivity.")
        method = "explicit_atom_mapping"
    else:
        matches = target.GetSubstructMatches(query, uniquify=False, useChirality=False, maxMatches=1001)
        if not matches:
            raise WorkflowError("Cannot map SDF chemistry to complex ligand connectivity; review the input.")
        if len(matches) > 1000:
            raise WorkflowError("Too many atom mappings; provide ligand.atom_mapping.")
        ranks = list(Chem.CanonicalRankAtoms(template,breakTies=False,includeChirality=True))
        signatures = set()
        for match in matches:
            labels = [0]*len(match)
            for i,j in enumerate(match):
                labels[j] = ranks[i]
            signatures.add(tuple(labels))
        if len(signatures) != 1:
            raise WorkflowError("Chemically ambiguous atom mapping; provide ligand.atom_mapping.")
        mapping = matches[0]
        method = "unique_or_chemically_equivalent_graph_mapping"
    coords = observed.GetConformer().GetPositions()
    result = Chem.Mol(template)
    result.RemoveAllConformers()
    conformer = Chem.Conformer(result.GetNumAtoms())
    conformer.Set3D(True)
    for i,j in enumerate(mapping):
        conformer.SetAtomPosition(i,coords[j])
    result.AddConformer(conformer)
    check = Chem.Mol(result)
    Chem.AssignStereochemistryFrom3D(check,confId=0,replaceExistingTags=True)
    Chem.AssignStereochemistry(check,cleanIt=True,force=True)
    if Chem.MolToSmiles(check,isomericSmiles=True) != Chem.MolToSmiles(template,isomericSmiles=True):
        raise WorkflowError("SDF stereochemistry and complex geometry disagree or are undefined.")
    result = Chem.AddHs(result,addCoords=True)
    if not np.allclose(result.GetConformer().GetPositions()[:len(mapping)],coords[list(mapping)],atol=1e-8,rtol=0):
        raise WorkflowError("Unexpected movement of ligand heavy atoms during preparation.")
    return result, {"method":method,"heavy_atoms":len(mapping),"heavy_atom_displacement_max_angstrom":0.0,
                    "coordinates_source":"complex PDB","chemistry_source":"ligand SDF",
                    "hydrogens":"Rebuilt from the supplied SDF chemical state."}

def strip_protein_hydrogens(text):
    """Rebuild explicit protein H under the declared protonation protocol.

    Remove H before PDB parsing, including malformed duplicate names from
    docking converters. Preserve heavy atoms and their original coordinates.
    """
    lines=text.splitlines(); serials=set()
    for line in lines:
        if line.startswith(("ATOM  ","HETATM")):
            element=line[76:78].strip().upper()
            name=line[12:16].strip().lstrip("0123456789")
            if element in ("H","D") or (not element and name.startswith(("H","D"))):
                serials.add(int(line[6:11]))
    kept=[]
    for line in lines:
        if line.startswith(("ATOM  ","HETATM")) and int(line[6:11]) in serials: continue
        if line.startswith("CONECT"):
            ids=[int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
            if not ids or ids[0] in serials: continue
            ids=[i for i in ids if i not in serials]
            if len(ids)<2: continue
            line="CONECT"+"".join(f"{i:5d}" for i in ids)
        kept.append(line)
    return "\n".join(kept)+"\n",len(serials)


def prepare_complex(config, paths):
    import numpy as np
    from openmm import app,unit
    from .preparation import prepare_protein
    source = relative_path(config["complex"]["file"],must_exist=True)
    selector = config.get("ligand")
    protein_text,ligand_text,selected = split_complex(source.read_text(),selector)
    removed_protein_hydrogens=0
    if config["protein"].get("rebuild_hydrogens", True):
        protein_text,removed_protein_hydrogens=strip_protein_hydrogens(protein_text)
    target = paths["md"]/"prepared"
    target.mkdir(parents=True,exist_ok=True)
    (target/"protein_input.pdb").write_text(protein_text)
    prepared, audit = None, None
    if selector is not None:
        from rdkit import Chem
        from openff.toolkit import Molecule
        template_path = relative_path(selector["template_sdf"],must_exist=True)
        templates = list(Chem.SDMolSupplier(str(template_path),removeHs=False))
        if len(templates)!=1 or templates[0] is None:
            raise WorkflowError("The ligand template SDF must contain exactly one valid molecule.")
        observed = Chem.MolFromPDBBlock(ligand_text,sanitize=False,removeHs=False,proximityBonding=True)
        prepared,audit = ligand_at_complex_coordinates(templates[0],observed,selector.get("atom_mapping"))
        Molecule.from_rdkit(prepared,allow_undefined_stereo=False,hydrogens_are_explicit=True)
        with Chem.SDWriter(str(target/"ligand.sdf")) as writer:
            writer.write(prepared)
        atomic_json(target/"ligand_preparation.json",audit)
    local_config = copy.deepcopy(config)
    local_config["protein"]["file"] = rel(target/"protein_input.pdb")
    prepare_protein(local_config,paths)
    original = app.PDBFile(StringIO(protein_text))
    fixed = app.PDBFile(str(target/"protein.pdb"))
    def positions_by_key(pdb):
        coords = pdb.positions.value_in_unit(unit.angstrom)
        pairs = [((a.residue.chain.id,a.residue.id,a.residue.insertionCode,a.name),np.asarray(coords[a.index]))
                 for a in pdb.topology.atoms() if a.element and a.element.symbol!="H" and a.residue.name in AMINO]
        if len({key for key,_ in pairs}) != len(pairs):
            raise WorkflowError("Duplicate protein chain/residue/atom identifiers; assign unique chain identities.")
        return dict(pairs)
    before,after = positions_by_key(original),positions_by_key(fixed)
    if not before or not set(before).issubset(after):
        raise WorkflowError("Protein atom identity changed during preparation.")
    maximum = max(float(np.linalg.norm(before[key]-after[key])) for key in before)
    if maximum > .002:
        raise WorkflowError("Protein heavy atoms moved unexpectedly during preparation.")
    chains = sorted({key[0] for key in before})
    system_type = config.get("system_type", "protein_protein" if selector is None else "protein_ligand")
    if system_type in ("protein_protein","ternary") and len(chains)<2:
        raise WorkflowError("protein_protein and ternary systems require at least two distinct protein chains.")
    if prepared is not None:
        # No inference of parameters for cross-residue covalent/coordination geometry.
        protein_coords = np.array(list(after.values()))
        ligand_coords = np.array(prepared.GetConformer().GetPositions()[:audit["heavy_atoms"]])
        from MDAnalysis.lib.distances import distance_array
        if distance_array(protein_coords,ligand_coords).min() < 1.7:
            raise WorkflowError("Very close protein-ligand heavy atoms require review for covalent bonding or severe clashes.")
    atomic_json(target/"complex_preparation.json",{
        "workflow":"md_from_complex_v1","system_type":system_type,
        "complex_file":config["complex"]["file"],"protein_chains":chains,
        "ligand_residue":({"chain":selected[0],"resid":selected[1],"insertion_code":selected[2],"resname":selected[3]}
                           if selected is not None else None),
        "ligand_heavy_atoms":audit["heavy_atoms"] if audit else 0,
        "ligand_mapping":audit["method"] if audit else None,
        "ligand_heavy_atom_displacement_max_angstrom":0.0 if audit else None,
        "protein_existing_heavy_atom_displacement_max_angstrom":maximum,
        "removed_existing_protein_hydrogens":removed_protein_hydrogens,
        "hydrogen_policy":"rebuild_from_pH" if config["protein"].get("rebuild_hydrogens", True) else "retain_then_add_missing",
        "pose_source":"user_supplied_complex","docking_performed":False,
        "note":"Input heavy-atom placement retained during preparation; minimization and MD subsequently move atoms."
    })
