"""Local, evidence-preserving ligand recovery; never equate sanitization with identity.

The public facade uses only the standard library. RDKit runs in a subprocess so
native bond-order inference is bounded even when a pathological graph hangs.
PDB-only inference always needs approval of the exact recorded candidate.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

SCHEMA_VERSION = 2
METHOD_VERSION = "evidence-graph-2.0"
MAX_ATOMS = 1000
MAX_MAPPINGS = 128


def _digest(value):
    data = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(data.encode()).hexdigest()


def _write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def _waiting(code, message, **kwargs):
    return {"status": "WAITING", "reason": code, "message": message,
            "required_decisions": [{"kind": code, "message": message}], **kwargs}


def _normalise_component(component):
    inventory = component.get("inventory", {})
    atoms = component.get("atoms")
    if atoms is None:
        selected = set(component.get("atom_ids", []))
        atoms = [atom for atom in inventory.get("atoms", []) if atom["id"] in selected]
    if not atoms or len(atoms) > MAX_ATOMS:
        raise ValueError(f"Organic component must contain 1–{MAX_ATOMS} atoms.")
    result, aliases = [], {}
    for index, atom in enumerate(atoms):
        stable_id = str(atom.get("id", atom.get("stable_id", f"atom_{index + 1}")))
        if stable_id in {item["id"] for item in result}:
            raise ValueError("Duplicate stable atom identifiers; resolve the inventory first.")
        xyz = atom.get("xyz", [atom.get("x"), atom.get("y"), atom.get("z")])
        if len(xyz) != 3 or any(v is None or not math.isfinite(float(v)) for v in xyz):
            raise ValueError(f"Non-finite or missing coordinates for {stable_id}.")
        charge = int(atom.get("formal_charge", atom.get("charge", 0)))
        result.append({"id": stable_id, "element": str(atom["element"]).strip().title(),
                       "xyz": [float(v) for v in xyz], "formal_charge": charge,
                       "charge_explicit": bool(atom.get("charge_explicit", atom.get("charge_text") or charge)),
                       "serial": atom.get("serial"), "name": atom.get("atom_name", atom.get("name", ""))})
        aliases[stable_id] = index
    # Serials are optional aliases only when unique. Stable ids remain authoritative.
    serials = [atom["serial"] for atom in result]
    for index, serial in enumerate(serials):
        if serial is not None and serials.count(serial) == 1:
            aliases[serial] = aliases[str(serial)] = index
    selected_ids = {atom["id"] for atom in result}
    raw_bonds = component.get("bonds")
    if raw_bonds is None:
        raw_bonds = [bond for bond in inventory.get("bonds", [])
                     if set(bond.get("atom_ids", [])).issubset(selected_ids)]
    bonds = {}
    for bond in raw_bonds:
        ends = bond.get("atom_ids", [bond.get("a"), bond.get("b")])
        if len(ends) != 2 or any(end not in aliases for end in ends):
            raise ValueError("Bond refers to an unknown or ambiguous atom.")
        a, b = sorted(aliases[end] for end in ends)
        if a == b:
            raise ValueError("Self-bond in input connectivity.")
        order = int(bond.get("multiplicity", bond.get("order", 1)))
        if order not in (1, 2, 3):
            raise ValueError("Unsupported CONECT multiplicity; explicit chemistry is required.")
        if (a, b) in bonds and bonds[a, b] != order:
            raise ValueError("Conflicting CONECT multiplicities.")
        bonds[a, b] = order
    return {"id": component.get("id", component.get("component_id", "ligand")),
            "atoms": result, "bonds": [[a, b, order] for (a, b), order in sorted(bonds.items())],
            "conect_complete": bool(component.get("conect_complete", component.get("connectivity_complete", False)))}


def _template_spec(template):
    if template is None:
        return None
    if isinstance(template, (str, Path)):
        template = {"path": str(template)}
    template = dict(template)
    if "path" in template:
        path = Path(template["path"]).expanduser().resolve()
        template["text"] = path.read_text(encoding="utf-8")
        template["format"] = path.suffix.lower().lstrip(".")
        template["source"] = str(path)
        del template["path"]
    if not template.get("smiles") and not template.get("text"):
        raise ValueError("A template requires an SDF/MOL2 file, mol block, or explicit SMILES.")
    # A supplied definition is chemical evidence, not permission to change explicit input charges/H.
    template.setdefault("trusted", True)
    return template


def recover_component(component, output_dir, decision=None, template=None, timeout=30):
    """Return READY/WAITING/DEPENDENCY_BLOCKED and persist auditable local artifacts.

    component: inventory component with atoms/bonds, or component + inventory.
    template: path or {smiles/text/format, atom_map: {stable_id: zero_based_index}}.
    decision: {accept_candidate_sha256: <hash from a previous WAITING result>}.
    Existing coordinates are always written in input atom order; no conformer
    generation or optimization substitutes the supplied binding pose.
    """
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    try:
        normalised = _normalise_component(component)
        template = _template_spec(template)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        result = _waiting("invalid_chemical_evidence", str(exc))
        _write_json(output / "chemistry.json", result)
        return result
    evidence = {"atoms": [{key: atom[key] for key in ("element", "xyz", "formal_charge", "charge_explicit")}
                          for atom in normalised["atoms"]],
                "bonds": normalised["bonds"], "conect_complete": normalised["conect_complete"],
                "template": {key: value for key, value in (template or {}).items() if key != "source"},
                "method": METHOD_VERSION,
                "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    payload = {"component": normalised, "evidence": evidence, "template": template,
               "output_dir": str(output), "decision": decision or {},
               "cache_dir": str(Path(component.get("cache_dir", output.parent / ".chemistry_cache")).resolve())}
    _write_json(output / "chemical_evidence.json", {"schema_version": SCHEMA_VERSION,
                "source_atom_mapping": normalised["atoms"], "evidence": evidence,
                "evidence_sha256": _digest(evidence)})
    input_path, result_path = output / ".worker_input.json", output / ".worker_result.json"
    _write_json(input_path, payload)
    if result_path.exists():
        result_path.unlink()
    try:
        process = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker",
                                  str(input_path), str(result_path)],
                                 capture_output=True, text=True, timeout=max(0.1, float(timeout)), check=False)
        (output / "chemistry_worker.log").write_text(process.stdout + process.stderr, encoding="utf-8")
        if result_path.exists():
            result = json.loads(result_path.read_text(encoding="utf-8"))
        else:
            result = _waiting("chemistry_worker_failed", f"Chemistry worker exited {process.returncode}; inspect chemistry_worker.log.")
    except subprocess.TimeoutExpired as exc:
        def _text(value):
            return value.decode(errors="replace") if isinstance(value, bytes) else (value or "")
        (output / "chemistry_worker.log").write_text(_text(exc.stdout) + _text(exc.stderr), encoding="utf-8")
        result = _waiting("chemistry_timeout", "Bounded chemical inference timed out. Supply or confirm an explicit SMILES/template.")
    result.update({"schema_version": SCHEMA_VERSION, "method_version": METHOD_VERSION,
                   "component_id": normalised["id"], "evidence_sha256": _digest(evidence),
                   "chemical_evidence": str(output / "chemical_evidence.json")})
    _write_json(output / "chemistry.json", result)
    return result


def _input_molecule(data, Chem):
    graph = Chem.RWMol()
    for atom in data["atoms"]:
        rd_atom = Chem.Atom(atom["element"])
        rd_atom.SetFormalCharge(atom["formal_charge"])
        rd_atom.SetNoImplicit(True)
        graph.AddAtom(rd_atom)
    types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for a, b, multiplicity in data["bonds"]:
        graph.AddBond(a, b, types[multiplicity])
    molecule = graph.GetMol()
    conformer = Chem.Conformer(len(data["atoms"]))
    conformer.Set3D(True)
    for index, atom in enumerate(data["atoms"]):
        conformer.SetAtomPosition(index, atom["xyz"])
    molecule.AddConformer(conformer, assignId=True)
    return molecule


def _heavy_graph(molecule, Chem):
    heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1]
    builder = Chem.RWMol()
    for index in heavy:
        builder.AddAtom(Chem.AtomFromSmarts(f"[#{molecule.GetAtomWithIdx(index).GetAtomicNum()}]"))
    inverse = {index: position for position, index in enumerate(heavy)}
    for bond in molecule.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in inverse and b in inverse:
            builder.AddBond(inverse[a], inverse[b], Chem.BondType.SINGLE)
            builder.ReplaceBond(builder.GetNumBonds() - 1, Chem.BondFromSmarts("~"))
    return builder.GetMol(), heavy


def _template_molecule(spec, Chem):
    if spec.get("smiles"):
        params = Chem.SmilesParserParams()
        params.removeHs = False
        mol = Chem.MolFromSmiles(spec["smiles"], params)
    elif spec.get("format") == "mol2":
        mol = Chem.MolFromMol2Block(spec["text"], removeHs=False)
    else:
        text = spec["text"]
        if text.count("$$$$") > 1:
            raise ValueError("Provide exactly one molecule per template.")
        mol = Chem.MolFromMolBlock(text.split("$$$$")[0], removeHs=False, strictParsing=True)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Template is invalid or contains multiple disconnected molecules.")
    return mol


def _signature(mol, Chem):
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return {"atoms": [[a.GetSymbol(), a.GetFormalCharge(), a.GetTotalNumHs(includeNeighbors=True),
                       a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else ""] for a in mol.GetAtoms()],
            "bonds": sorted([min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                              max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType()), str(b.GetStereo())]
                             for b in mol.GetBonds())}


def _with_input_coordinates(molecule, data, Chem):
    molecule = Chem.Mol(molecule)
    molecule.RemoveAllConformers()
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    conformer.Set3D(True)
    for index, atom in enumerate(data["atoms"]):
        conformer.SetAtomPosition(index, atom["xyz"])
    molecule.AddConformer(conformer, assignId=True)
    return molecule


def _map_template(raw, data, spec, Chem, rdDetermineBonds):
    template = _template_molecule(spec, Chem)
    # Explicit input hydrogens constrain protonation; missing H may be supplied by a trusted definition.
    template = Chem.AddHs(template, addCoords=False)
    graph = Chem.Mol(raw)
    has_graph = bool(data["bonds"])
    if not has_graph:
        rdDetermineBonds.DetermineConnectivity(graph, useHueckel=False, useVdw=True, covFactor=1.3)
    query, input_heavy = _heavy_graph(graph, Chem)
    template_query, template_heavy = _heavy_graph(template, Chem)
    if len(input_heavy) != len(template_heavy) or query.GetNumBonds() != template_query.GetNumBonds():
        raise ValueError("Template heavy atom count/connectivity does not match the input component.")
    index_by_id = {atom["id"]: index for index, atom in enumerate(data["atoms"])}
    supplied_map = spec.get("atom_map")
    if supplied_map:
        try:
            mapping = {index_by_id[key]: int(value) for key, value in supplied_map.items()}
            matches = [tuple(mapping[index] for index in input_heavy)]
        except (KeyError, ValueError):
            raise ValueError("Explicit atom_map must cover every input heavy atom using stable IDs.")
    else:
        matches = template.GetSubstructMatches(query, uniquify=False, useChirality=False, maxMatches=MAX_MAPPINGS + 1)
        if len(matches) > MAX_MAPPINGS:
            raise ValueError("Template mapping exceeds the bounded ambiguity limit; supply an explicit atom_map.")
    candidates = {}
    failures = []
    for match in matches:
        try:
            mapping = dict(zip(input_heavy, match))
            if len(set(mapping.values())) != len(mapping):
                raise ValueError("Template mapping is not one-to-one.")
            for src, dst in mapping.items():
                if dst < 0 or dst >= template.GetNumAtoms() or raw.GetAtomWithIdx(src).GetAtomicNum() != template.GetAtomWithIdx(dst).GetAtomicNum():
                    raise ValueError("Template mapping changes an element.")
                charge = data["atoms"][src]
                if charge["charge_explicit"] and template.GetAtomWithIdx(dst).GetFormalCharge() != charge["formal_charge"]:
                    raise ValueError("Template changes an explicitly marked input charge site.")
            for bond in query.GetBonds():
                a, b = input_heavy[bond.GetBeginAtomIdx()], input_heavy[bond.GetEndAtomIdx()]
                if template.GetBondBetweenAtoms(mapping[a], mapping[b]) is None:
                    raise ValueError("Template mapping changes connectivity.")
            for index, atom in enumerate(data["atoms"]):
                if atom["element"] != "H":
                    continue
                neighbors = [n.GetIdx() for n in graph.GetAtomWithIdx(index).GetNeighbors()]
                if len(neighbors) != 1 or neighbors[0] not in mapping:
                    raise ValueError("Each explicit hydrogen needs one mapped heavy-atom neighbor.")
                available = [n.GetIdx() for n in template.GetAtomWithIdx(mapping[neighbors[0]]).GetNeighbors()
                             if n.GetAtomicNum() == 1 and n.GetIdx() not in mapping.values()]
                if not available:
                    raise ValueError("Template cannot preserve the input explicit hydrogen count/protonation.")
                mapping[index] = available[0]
            order = [mapping[index] for index in range(len(data["atoms"]))]
            order += [index for index in range(template.GetNumAtoms()) if index not in order]
            candidate = Chem.RenumberAtoms(template, order)
            candidate = _with_input_coordinates(candidate, data, Chem)
            # Additional H coordinates are generated around the fixed input heavy-atom pose.
            missing_h = candidate.GetNumAtoms() - len(data["atoms"])
            if missing_h:
                trim = Chem.RWMol(candidate)
                for index in reversed(range(len(data["atoms"]), candidate.GetNumAtoms())):
                    trim.RemoveAtom(index)
                trimmed = trim.GetMol()
                for atom in trimmed.GetAtoms():
                    atom.SetNoImplicit(False)
                trimmed.UpdatePropertyCache(strict=False)
                candidate = Chem.AddHs(trimmed, addCoords=True)
            Chem.SanitizeMol(candidate)
            expected = Chem.Mol(candidate)
            observed = Chem.Mol(candidate)
            Chem.AssignStereochemistryFrom3D(observed, replaceExistingTags=True)
            Chem.AssignStereochemistry(expected, cleanIt=True, force=True)
            Chem.AssignStereochemistry(observed, cleanIt=True, force=True)
            for index in input_heavy:
                a, b = expected.GetAtomWithIdx(index), observed.GetAtomWithIdx(index)
                if a.HasProp("_CIPCode") and b.HasProp("_CIPCode") and a.GetProp("_CIPCode") != b.GetProp("_CIPCode"):
                    raise ValueError("Template stereochemistry conflicts with the supplied 3D coordinates.")
            for bond in expected.GetBonds():
                actual = observed.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                if bond.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ) and actual.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ) and bond.GetStereo() != actual.GetStereo():
                    raise ValueError("Template double-bond stereochemistry conflicts with the supplied 3D coordinates.")
            # Template-unspecified centers use the supplied geometry, never a generated conformer.
            for index in input_heavy:
                if expected.GetAtomWithIdx(index).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                    expected.GetAtomWithIdx(index).SetChiralTag(observed.GetAtomWithIdx(index).GetChiralTag())
            candidates[_digest(_signature(expected, Chem))] = expected
        except (ValueError, RuntimeError) as exc:
            failures.append(str(exc))
    if len(candidates) != 1:
        detail = failures[0] if failures and not candidates else "chemically different atom mappings remain"
        raise ValueError(f"Template identity/mapping is not unique: {detail}.")
    return next(iter(candidates.values())), "explicit_template", {"connectivity_from_coordinates": not has_graph,
                     "matched_mapping_count": len(matches), "template_source": spec.get("source", "supplied definition")}


def _fixed_site_valence_candidates(raw, data, Chem):
    """Enumerate a bounded integer-valence model while fixing all input charges.

    This is a fallback candidate generator, not an identity/protonation oracle.
    Deliberately excludes hypervalent/metal/radical chemistry instead of silently
    expanding the chemical model until something sanitizes.
    """
    valence = {("H", 0): 1, ("C", 0): 4, ("C", 1): 3, ("C", -1): 3,
               ("N", 0): 3, ("N", 1): 4, ("N", -1): 2,
               ("O", 0): 2, ("O", 1): 3, ("O", -1): 1,
               ("F", 0): 1, ("Cl", 0): 1, ("Br", 0): 1, ("I", 0): 1,
               ("S", 0): 2}
    deficits = []
    for index, atom in enumerate(data["atoms"]):
        key = (atom["element"], atom["formal_charge"])
        if key not in valence:
            raise ValueError("Fixed-site candidate generation does not support this valence/element; supply an explicit chemical definition.")
        used = sum(order for a, b, order in data["bonds"] if index in (a, b))
        deficit = valence[key] - used
        if deficit < 0:
            raise ValueError("Input connectivity/charge is incompatible with the supported fixed-site valence model.")
        deficits.append(deficit)
    edges = [(a, b, 3 - order) for a, b, order in data["bonds"]
             if data["atoms"][a]["element"] != "H" and data["atoms"][b]["element"] != "H"]
    answers, visited = [], 0

    def search(demand, pending, assignments):
        nonlocal visited
        visited += 1
        if visited > 100000 or len(answers) >= 128:
            raise ValueError("Fixed-site bond-order enumeration reached its ambiguity bound; supply explicit chemistry.")
        if any(value < 0 for value in demand):
            return
        if not any(demand):
            answers.append(dict(assignments))
            return
        capacities = [0] * len(demand)
        available = {}
        for edge in pending:
            a, b, maximum = edges[edge]
            capacity = min(maximum, demand[a], demand[b])
            if capacity:
                capacities[a] += capacity
                capacities[b] += capacity
                available.setdefault(a, []).append((edge, capacity))
                available.setdefault(b, []).append((edge, capacity))
        if any(demand[i] > capacities[i] for i in range(len(demand))):
            return
        # Fix a complete incident-edge allocation at the most constrained atom.
        atom = min((i for i, need in enumerate(demand) if need),
                   key=lambda i: (len(available[i]), capacities[i] - demand[i]))
        incident = available[atom]

        def allocate(position, remaining, selected):
            if position == len(incident):
                if remaining:
                    return
                changed = list(demand)
                updated = dict(assignments)
                for edge, increment in selected:
                    a, b, _ = edges[edge]
                    changed[a] -= increment
                    changed[b] -= increment
                    updated[edge] = increment
                search(changed, pending - {edge for edge, _ in incident}, updated)
                return
            edge, capacity = incident[position]
            lower = max(0, remaining - sum(cap for _, cap in incident[position + 1:]))
            for increment in range(lower, min(capacity, remaining) + 1):
                allocate(position + 1, remaining - increment, selected + [(edge, increment)])

        allocate(0, demand[atom], [])

    search(deficits, set(range(len(edges))), {})
    unique = {}
    for answer in answers:
        mol = Chem.RWMol(raw)
        for edge, increment in answer.items():
            a, b, _ = edges[edge]
            bond = mol.GetBondBetweenAtoms(a, b)
            order = round(bond.GetBondTypeAsDouble()) + increment
            bond.SetBondType({1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[order])
        mol = mol.GetMol()
        try:
            Chem.SanitizeMol(mol)
            Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
            if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
                continue
            key = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
            unique[key] = mol
        except (ValueError, RuntimeError):
            continue
    if len(unique) != 1:
        raise ValueError(f"Fixed-charge valence assignment yields {len(unique)} distinct chemical candidates; supply an explicit chemical definition.")
    return next(iter(unique.values())), {"fixed_site_valence_assignments": len(answers),
              "distinct_sanitized_candidates": len(unique), "enumeration_nodes": visited,
              "valence_model": "closed-shell CHNO, halogens and divalent sulfur; explicit H and fixed charge sites"}


def _infer(raw, data, Chem, rdDetermineBonds):
    if not data["bonds"] or len(Chem.GetMolFrags(raw)) != 1:
        raise ValueError("A connected authoritative CONECT graph or explicit SMILES/template is required; proximity alone is not chemical identity.")
    if not any(atom["element"] == "H" for atom in data["atoms"]):
        raise ValueError("PDB-only bond-order recovery lacks explicit hydrogen/protonation evidence. Supply a SMILES or template.")
    # Existing marked charges constrain sites; their sum is only a trial total charge.
    total = sum(atom["formal_charge"] for atom in data["atoms"])
    candidate = Chem.Mol(raw)
    extra_checks = {}
    try:
        rdDetermineBonds.DetermineBondOrders(candidate, charge=total, allowChargedFragments=True,
                                            embedChiral=True, useAtomMap=False)
        Chem.SanitizeMol(candidate)
        if any(candidate.GetAtomWithIdx(index).GetFormalCharge() != atom["formal_charge"] for index, atom in enumerate(data["atoms"])):
            raise ValueError("Unconstrained inference changes fixed charge sites.")
    except (ValueError, RuntimeError):
        candidate, extra_checks = _fixed_site_valence_candidates(raw, data, Chem)
    for index, source in enumerate(data["atoms"]):
        assigned = candidate.GetAtomWithIdx(index).GetFormalCharge()
        if assigned != source["formal_charge"]:
            raise ValueError("Bond-order inference relocates charge or introduces unsupported charge separation; an explicit chemical definition is needed.")
        if candidate.GetAtomWithIdx(index).GetNumRadicalElectrons():
            raise ValueError("Radical chemistry requires a dedicated parameterization route.")
        if source["element"] != "H" and candidate.GetAtomWithIdx(index).GetNumImplicitHs():
            raise ValueError("PDB evidence does not fully determine hydrogens/protonation.")
    edges = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in candidate.GetBonds()}
    if edges != {(a, b) for a, b, _ in data["bonds"]}:
        raise ValueError("Inference changed the authoritative connectivity graph.")
    for a, b, order in data["bonds"]:
        bond = candidate.GetBondBetweenAtoms(a, b)
        if order > 1 and not bond.GetIsAromatic() and round(bond.GetBondTypeAsDouble()) != order:
            raise ValueError("Inferred bond order contradicts an explicit repeated CONECT entry.")
    return candidate, "pdb_graph_candidate", {"charge_evidence_sum": total,
                "charge_sites_preserved": True, "extra_charge_separation": False,
                "identity_proven": False, "inference_is_not_uniqueness_proof": True, **extra_checks}


def _worker(payload):
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import rdDetermineBonds
    except ImportError as exc:
        return {"status": "DEPENDENCY_BLOCKED", "reason": "rdkit_missing",
                "message": f"RDKit with rdDetermineBonds is required in the Ubuntu runtime: {exc}",
                "required_decisions": []}
    data, output = payload["component"], Path(payload["output_dir"])
    cache_key = _digest({"evidence": payload["evidence"], "rdkit_version": rdBase.rdkitVersion})
    cache_dir = Path(payload["cache_dir"])
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / (cache_key + ".json")
    cache_hit = False
    try:
        raw = _input_molecule(data, Chem)
        if cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if cached.get("molblock_sha256") != hashlib.sha256(cached["molblock"].encode()).hexdigest():
                raise ValueError("Chemical evidence cache checksum failed; remove the affected cache entry and retry.")
            molecule = Chem.MolFromMolBlock(cached["molblock"], removeHs=False)
            if molecule is None or cached.get("cache_key") != cache_key:
                raise ValueError("Chemical evidence cache is corrupt.")
            method, checks = cached["method"], cached["checks"]
            cache_hit = True
        elif payload.get("template"):
            molecule, method, checks = _map_template(raw, data, payload["template"], Chem, rdDetermineBonds)
        else:
            molecule, method, checks = _infer(raw, data, Chem, rdDetermineBonds)
        if molecule.GetNumAtoms() < len(data["atoms"]):
            raise ValueError("Chemical recovery lost input atoms.")
        # Molfile precision is 4 decimals in Å; compare before serialization.
        for index, atom in enumerate(data["atoms"]):
            if molecule.GetAtomWithIdx(index).GetSymbol() != atom["element"]:
                raise ValueError("Recovered atom order no longer matches the original element mapping.")
            point = molecule.GetConformer().GetAtomPosition(index)
            if max(abs(float(point[axis]) - atom["xyz"][axis]) for axis in range(3)) > 0.00011:
                raise ValueError("Chemical recovery changed an input coordinate.")
        signature = _signature(molecule, Chem)
        canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
        candidate_hash = _digest({"cache_key": cache_key, "chemical_signature": signature,
                                  "canonical_isomeric_smiles": canonical_smiles})
        approved = payload["decision"].get("accept_candidate_sha256") == candidate_hash
        trusted = bool(payload.get("template", {}).get("trusted", False)) if payload.get("template") else False
        ready = trusted or approved
        filename = "ligand.sdf" if ready else "candidate.sdf"
        sdf_path = output / filename
        molecule.SetProp("NEXUS_CHEMISTRY_METHOD", method)
        molecule.SetProp("NEXUS_EVIDENCE_SHA256", _digest(payload["evidence"]))
        molecule.SetProp("NEXUS_CANDIDATE_SHA256", candidate_hash)
        molecule.SetProp("NEXUS_CHEMISTRY_CONFIRMED", str(ready).lower())
        molecule.SetProp("NEXUS_ORIGINAL_ATOM_IDS", json.dumps([atom["id"] for atom in data["atoms"]]))
        with Chem.SDWriter(str(sdf_path)) as writer:
            writer.write(molecule)
        mapping = [{"original_id": atom["id"], "original_serial": atom["serial"],
                    "original_name": atom["name"], "sdf_atom_index": index, "element": atom["element"],
                    "original_xyz_angstrom": atom["xyz"], "added": False}
                   for index, atom in enumerate(data["atoms"])]
        mapping += [{"original_id": None, "sdf_atom_index": index,
                     "element": molecule.GetAtomWithIdx(index).GetSymbol(), "added": True,
                     "reason": "Hydrogen supplied by the explicit chemical definition"}
                    for index in range(len(data["atoms"]), molecule.GetNumAtoms())]
        _write_json(output / "atom_map.json", mapping)
        molblock = Chem.MolToMolBlock(molecule)
        _write_json(cache_path, {"cache_key": cache_key, "molblock": molblock,
                                "molblock_sha256": hashlib.sha256(molblock.encode()).hexdigest(),
                                "method": method, "checks": checks, "rdkit_version": rdBase.rdkitVersion})
        result = {"status": "READY" if ready else "WAITING", "method": method,
                  "message": "Chemical definition verified against input evidence." if ready else
                             "A bounded PDB-derived candidate is available. Confirm its chemical identity/protonation before parameterization.",
                  "candidate_sha256": candidate_hash, "canonical_isomeric_smiles": canonical_smiles,
                  "formal_charge": Chem.GetFormalCharge(molecule), "atom_count": molecule.GetNumAtoms(),
                  "input_atom_count": len(data["atoms"]), "added_hydrogens": molecule.GetNumAtoms() - len(data["atoms"]),
                  "sdf": str(sdf_path), "atom_map": str(output / "atom_map.json"), "checks": checks,
                  "rdkit_version": rdBase.rdkitVersion, "cache_key": cache_key, "cache_hit": cache_hit,
                  "approved_by": "explicit_template" if trusted else ("exact_candidate_confirmation" if approved else None),
                  "required_decisions": []}
        if not ready:
            result["required_decisions"] = [{"kind": "confirm_ligand_chemistry", "component_id": data["id"],
                    "candidate_sha256": candidate_hash, "canonical_isomeric_smiles": canonical_smiles,
                    "formal_charge": Chem.GetFormalCharge(molecule), "candidate_sdf": str(sdf_path),
                    "message": "Review the candidate's bond orders, protonation, formal charge sites and stereochemistry.",
                    "decision": {"accept_candidate_sha256": candidate_hash},
                    "alternative": "Provide an explicit SMILES or SDF chemical definition."}]
        return result
    except (ValueError, RuntimeError, KeyError, IndexError, OSError) as exc:
        return _waiting("chemical_identity_unresolved", str(exc), rdkit_version=rdBase.rdkitVersion,
                        cache_key=cache_key, cache_hit=cache_hit)


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] != "--worker":
        raise SystemExit("This module is invoked through recover_component().")
    result = _worker(json.loads(Path(sys.argv[2]).read_text(encoding="utf-8")))
    _write_json(Path(sys.argv[3]), result)
