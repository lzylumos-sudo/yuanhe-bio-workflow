"""Lossless, standard-library PDB intake for the v2 MD frontend.

Inventory never drops atom records, infers bond orders, or repairs chemistry.
Normalization changes identifiers only, saves every original identity, and
requires explicit model/conformer selection. Distances are diagnostics, never
new covalent bonds. Coordinates are in Angstrom throughout this module.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import string


SCHEMA_VERSION = "nexus-pdb-inventory-2.0"
PROTEIN_NAMES = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
PROTEIN_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "CYX": "CYS", "CYM": "CYS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS",
}
HISTIDINE_VARIANTS = {"HID": "HID", "HIE": "HIE", "HIP": "HIP", "HSD": "HID", "HSE": "HIE", "HSP": "HIP"}
WATER_NAMES = set("HOH WAT SOL TIP TIP3 TP3 DOD H2O".split())
ION_ELEMENTS = set("LI NA K RB CS MG CA SR BA ZN CD HG FE CO NI CU MN CR AG AU AL CL BR I F".split())
COFACTOR_NAMES = set("HEM HEC HEO HEME FAD FMN NAD NAP NDP ATP ADP AMP GDP GTP SAM SAH COA PLP".split())
GLYCAN_NAMES = set("NAG NDG BMA MAN GLC BGC GAL FUC FUL SIA NANA".split())
LIPID_NAMES = set("POP DOP DPPC POC PEA CHL CHL1 CLR".split())
CHAIN_ALPHABET = string.ascii_uppercase + string.ascii_lowercase + string.digits


class InventoryDecisionRequired(ValueError):
    """Raised with machine-readable decisions; nothing has been normalized."""

    def __init__(self, decisions):
        self.decisions = decisions
        super().__init__("Input selection required: " + "; ".join(d["message"] for d in decisions))


def _issue(kind, message, **details):
    return {"type": kind, "message": message, **details}


def _integer(text):
    return int(text.strip()) if text.strip() else None


def _element(atom_name, raw_name, explicit, record, resname):
    if explicit.strip():
        return explicit.strip().upper(), "element_column"
    # PDB atom-name alignment separates alpha carbon ' CA ' from calcium 'CA  '.
    letters = "".join(c for c in atom_name if c.isalpha()).upper()
    if not letters:
        return "", "missing"
    if raw_name.startswith(" ") or (raw_name and raw_name[0].isdigit()):
        return letters[0], "atom_name_alignment"
    if record == "ATOM" or resname in PROTEIN_NAMES or resname in PROTEIN_ALIASES:
        return letters[0], "protein_atom_name"
    if letters[:2] in ION_ELEMENTS or letters[:2] in {"SE", "SI"}:
        return letters[:2], "atom_name_alignment"
    return letters[0], "atom_name_alignment"


def _charge(field):
    value = field.strip()
    if not value:
        return 0
    if len(value) == 2 and value[0].isdigit() and value[1] in "+-":
        return int(value[0]) * (1 if value[1] == "+" else -1)
    if len(value) == 2 and value[0] in "+-" and value[1].isdigit():
        return int(value[1]) * (1 if value[0] == "+" else -1)
    raise ValueError("invalid PDB formal-charge field: " + repr(field))


def inventory_pdb(path):
    """Read all ATOM/HETATM lines before classifying anything.

    Components contain expanded ``atoms`` and ``bonds`` for downstream chemistry.
    ``connectivity_complete`` means graph coverage only, NOT chemical correctness.
    Blank charge columns are recorded as non-explicit, not evidence of neutrality.
    """
    source = Path(path).expanduser().resolve()
    data = source.read_bytes()
    try:
        contents = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("PDB must be ASCII or UTF-8; original file was not modified") from exc
    atoms, residues, records, conect, required, warnings = [], [], [], [], [], []
    residue_lookup = {}
    explicit_models = any(line[:6].strip() == "MODEL" for line in contents.splitlines())
    model, active_model, model_serial = 1, not explicit_models, None
    model_count = 0
    segment_by_model = defaultdict(int)
    models = []

    for lineno, raw in enumerate(contents.splitlines(), 1):
        line = raw.ljust(80)
        record = line[:6].strip()
        if record == "MODEL":
            model_count += 1
            model = model_count
            model_serial = line[10:14].strip()
            active_model = True
            models.append({"id": model, "original_serial": model_serial, "line_number": lineno})
        elif record == "ENDMDL":
            active_model = False
        elif record == "TER":
            records.append({"record": record, "line_number": lineno, "model": model, "raw_line": raw})
            segment_by_model[model] += 1
        elif record in {"LINK", "SSBOND", "SEQRES", "MODRES", "CISPEP"}:
            records.append({"record": record, "line_number": lineno, "model": model if active_model else None, "raw_line": raw})
        elif record == "CONECT":
            fields = [line[i:i + 5] for i in range(6, min(len(raw), 31), 5)]
            try:
                numbers = [_integer(field) for field in fields]
                if not numbers or numbers[0] is None:
                    raise ValueError("missing source serial")
                conect.append({"source_serial": numbers[0], "target_serials": [v for v in numbers[1:] if v is not None],
                               "model": model if explicit_models and active_model else None, "line_number": lineno, "raw_line": raw})
            except ValueError as exc:
                required.append(_issue("invalid_conect", str(exc), line_number=lineno, raw_line=raw))
        elif record in {"ATOM", "HETATM"}:
            atom_id = "a%06d" % (len(atoms) + 1)
            errors = []
            try:
                serial = _integer(line[6:11])
                resseq = _integer(line[22:26])
                if serial is None or resseq is None:
                    raise ValueError("missing serial or residue sequence number")
            except ValueError as exc:
                serial, resseq = None, None
                errors.append(str(exc))
            try:
                xyz = [float(line[start:start + 8]) for start in (30, 38, 46)]
                if not all(math.isfinite(value) for value in xyz):
                    raise ValueError("non-finite coordinates")
            except ValueError as exc:
                xyz = None
                errors.append(str(exc))
            try:
                charge = _charge(line[78:80])
            except ValueError as exc:
                charge = None
                errors.append(str(exc))
            atom_name, resname = line[12:16].strip(), line[17:20].strip()
            element, element_source = _element(atom_name, line[12:16], line[76:78], record, resname)
            atom = {"id": atom_id, "stable_id": atom_id, "line_number": lineno, "record": record,
                    "model": model, "model_serial": model_serial, "segment": segment_by_model[model], "serial": serial,
                    "serial_text": line[6:11], "chain": line[21], "resseq": resseq, "resseq_text": line[22:26],
                    "icode": line[26].strip(), "altloc": line[16].strip(), "resname": resname,
                    "atom_name": atom_name, "name": atom_name, "element": element, "element_source": element_source,
                    "xyz": xyz, "x": xyz[0] if xyz else None, "y": xyz[1] if xyz else None, "z": xyz[2] if xyz else None,
                    "formal_charge": charge, "charge": charge, "charge_text": line[78:80],
                    "charge_explicit": bool(line[78:80].strip()), "occupancy_text": line[54:60],
                    "bfactor_text": line[60:66], "raw_line": raw}
            if explicit_models and not active_model:
                errors.append("atom record outside MODEL/ENDMDL in a multi-model file")
            if not element:
                errors.append("element could not be determined")
            if element_source != "element_column":
                warnings.append(_issue("inferred_element", "Element inferred from PDB atom-name alignment; review before chemistry recovery.", atom_id=atom_id, element=element))
            # A residue name is part of source identity. Never collapse identical atom names.
            residue_key = (model, atom["segment"], atom["chain"], line[22:26], atom["icode"], resname)
            if residue_key not in residue_lookup:
                residue_id = "r%04d" % (len(residues) + 1)
                residue = {"id": residue_id, "model": model, "segment": atom["segment"], "chain": atom["chain"],
                           "resseq": resseq, "icode": atom["icode"], "resname": resname, "atom_ids": []}
                residue_lookup[residue_key] = residue
                residues.append(residue)
            residue = residue_lookup[residue_key]
            atom["residue_id"] = residue["id"]
            residue["atom_ids"].append(atom_id)
            atoms.append(atom)
            if errors:
                required.append(_issue("invalid_atom_record", "; ".join(errors), atom_id=atom_id, model=model, line_number=lineno))

    if not atoms:
        raise ValueError("PDB has no ATOM/HETATM records")
    if not explicit_models:
        models = [{"id": 1, "original_serial": None, "line_number": None}]
    by_id = {atom["id"]: atom for atom in atoms}
    serial_lookup = defaultdict(list)
    for atom in atoms:
        serial_lookup[(atom["model"], atom["serial"])].append(atom["id"])
    for (model_id, serial), ids in serial_lookup.items():
        if len(ids) > 1:
            required.append(_issue("duplicate_serial", "Duplicate atom serial makes CONECT identity ambiguous; provide an unambiguous source export.", model=model_id, serial=serial, atom_ids=ids))
    if len(models) > 1:
        required.append(_issue("model", "Select one MODEL; models must never be combined.", options=[m["id"] for m in models], decision_key="model"))

    directed = Counter()
    evidence = defaultdict(set)
    for record in conect:
        applicable = [record["model"]] if record["model"] is not None else [m["id"] for m in models]
        for model_id in applicable:
            srcs = serial_lookup.get((model_id, record["source_serial"]), [])
            for target_serial, count in Counter(record["target_serials"]).items():
                dsts = serial_lookup.get((model_id, target_serial), [])
                if len(srcs) != 1 or len(dsts) != 1:
                    required.append(_issue("unresolved_conect", "CONECT serial is absent or ambiguous.", model=model_id, source_serial=record["source_serial"], target_serial=target_serial, line_number=record["line_number"]))
                    continue
                if srcs[0] == dsts[0]:
                    required.append(_issue("self_bond", "CONECT connects an atom to itself.", atom_id=srcs[0], line_number=record["line_number"]))
                    continue
                # Same directed pair on repeated lines is continuation/multiplicity evidence.
                directed[(srcs[0], dsts[0])] += count
                evidence[tuple(sorted((srcs[0], dsts[0])))].add(record["line_number"])
    bonds = []
    for pair in sorted(evidence):
        forward, reverse = directed[pair], directed[(pair[1], pair[0])]
        multiplicity = max(forward, reverse)
        bonds.append({"atom_ids": list(pair), "a": pair[0], "b": pair[1], "multiplicity": multiplicity,
                      "order": multiplicity, "source": "CONECT", "directed_counts": [forward, reverse],
                      "evidence_lines": sorted(evidence[pair])})
        if forward and reverse and forward != reverse:
            warnings.append(_issue("asymmetric_conect_multiplicity", "CONECT directions disagree; multiplicity is evidence, not established bond order.", atom_ids=list(pair), directed_counts=[forward, reverse]))

    for residue in residues:
        member_atoms = [by_id[aid] for aid in residue["atom_ids"]]
        names = {a["atom_name"] for a in member_atoms}
        resname = residue["resname"]
        elements = {a["element"] for a in member_atoms}
        if resname in PROTEIN_NAMES or resname in PROTEIN_ALIASES:
            kind = "protein"
        elif resname in WATER_NAMES:
            kind = "water"
        elif len(member_atoms) == 1 and elements <= ION_ELEMENTS:
            kind = "ion"
        elif resname in COFACTOR_NAMES:
            kind = "cofactor"
        elif resname in GLYCAN_NAMES:
            kind = "glycan"
        elif resname in LIPID_NAMES:
            kind = "lipid"
        elif {"N", "CA", "C"} <= names:
            kind = "modified_protein"
        else:
            kind = "ligand"
        residue["kind"] = kind
        residue["atom_count"] = len(member_atoms)
        labels = sorted({a["altloc"] for a in member_atoms if a["altloc"]})
        residue["altlocs"] = labels
        if len(labels) > 1:
            required.append(_issue("altloc", "Select one conformer label for the whole residue; atom-wise mixtures are forbidden.", residue_id=residue["id"], options=labels, decision_key="altloc." + residue["id"]))
        if kind == "protein":
            missing = sorted({"N", "CA", "C", "O"} - names)
            if missing:
                required.append(_issue("missing_backbone", "Protein residue is missing backbone atoms; automatic side-chain repair is insufficient.", residue_id=residue["id"], missing_atoms=missing))
            if "HXT" in names and "OXT" not in names:
                required.append(_issue("terminus_hxt_without_oxt", "HXT is present without OXT. Confirm terminal chemistry; normalization must neither rename nor delete HXT.", residue_id=residue["id"], chain=residue["chain"], resseq=residue["resseq"], decision_key="termini." + residue["id"], options=["standard_charged", "special_parameterization"]))
            if resname in PROTEIN_ALIASES:
                residue["canonical_resname"] = PROTEIN_ALIASES[resname]
                residue["protonation_variant"] = HISTIDINE_VARIANTS.get(resname, resname)
                if resname not in HISTIDINE_VARIANTS:
                    required.append(_issue("unsupported_protonation_variant", "This explicit residue variant is recognized but has no validated preservation adapter; do not silently reset its protonation/disulfide state.",
                                           residue_id=residue["id"], variant=resname))
            duplicate_names = [name for (name, alt), count in Counter((a["atom_name"], a["altloc"]) for a in member_atoms).items() if count > 1]
            if duplicate_names:
                required.append(_issue("duplicate_protein_atom_name", "Protein template atom identity is ambiguous within the same residue and conformer.", residue_id=residue["id"], model=residue["model"], atom_names=duplicate_names))

    # Protein groups are chain/TER fragments. Organic residues merge only via
    # explicit CONECT links; protein-ligand links remain visible as special systems.
    parent = {r["id"]: r["id"] for r in residues}
    residue_by_id = {r["id"]: r for r in residues}

    def root(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(a, b):
        parent[root(b)] = root(a)

    proteins = {}
    for residue in residues:
        if residue["kind"] == "protein":
            key = (residue["model"], residue["segment"], residue["chain"])
            if key in proteins:
                union(proteins[key], residue["id"])
            else:
                proteins[key] = residue["id"]
    for bond in bonds:
        ra, rb = [residue_by_id[by_id[aid]["residue_id"]] for aid in bond["atom_ids"]]
        if ra["id"] == rb["id"]:
            continue
        kinds = {ra["kind"], rb["kind"]}
        if kinds == {"ligand"}:
            union(ra["id"], rb["id"])
        elif "protein" in kinds and kinds != {"protein"}:
            required.append(_issue("covalent_or_metal_link", "Explicit cross-component CONECT requires a dedicated covalent/metal preparation route.", atom_ids=bond["atom_ids"], residue_ids=[ra["id"], rb["id"]]))
        elif kinds == {"protein"}:
            names = {by_id[aid]["atom_name"] for aid in bond["atom_ids"]}
            if names != {"C", "N"} or ra["chain"] != rb["chain"] or ra["segment"] != rb["segment"]:
                required.append(_issue("protein_crosslink", "Explicit non-peptide inter-residue CONECT requires crosslink/disulfide review.", atom_ids=bond["atom_ids"], residue_ids=[ra["id"], rb["id"]]))
    groups = defaultdict(list)
    for residue in residues:
        groups[root(residue["id"])].append(residue)
    components = []
    for members in sorted(groups.values(), key=lambda group: int(group[0]["atom_ids"][0][1:])):
        component_id = "c%03d" % (len(components) + 1)
        ids = [aid for residue in members for aid in residue["atom_ids"]]
        idset = set(ids)
        component_atoms = [by_id[aid] for aid in ids]
        component_bonds = [b for b in bonds if set(b["atom_ids"]) <= idset]
        connected_ids = {aid for b in component_bonds for aid in b["atom_ids"]}
        for atom in component_atoms:
            atom["component_id"] = component_id
        component = {"id": component_id, "component_id": component_id, "kind": members[0]["kind"],
                     "model": members[0]["model"], "residue_ids": [r["id"] for r in members],
                     "atom_ids": ids, "atoms": component_atoms, "bonds": component_bonds,
                     "atom_count": len(ids), "heavy_atom_count": sum(a["element"] not in {"H", "D"} for a in component_atoms),
                     "hydrogen_count": sum(a["element"] in {"H", "D"} for a in component_atoms),
                     "total_formal_charge": sum(a["formal_charge"] or 0 for a in component_atoms),
                     "explicit_charge_atom_ids": [a["id"] for a in component_atoms if a["charge_explicit"]],
                     "connectivity_complete": bool(component_bonds) and connected_ids == idset,
                     "conect_complete": bool(component_bonds) and connected_ids == idset,
                     "resnames": sorted({r["resname"] for r in members}), "chains": sorted({r["chain"] for r in members})}
        components.append(component)
        if component["kind"] in {"cofactor", "glycan", "lipid", "modified_protein"}:
            required.append(_issue("special_component", "This component needs a dedicated preparation/parameterization capability or explicit exclusion.",
                                   component_id=component_id, kind=component["kind"], resnames=component["resnames"]))
        if component["kind"] == "ion" and any(a["element"] not in {"NA", "K", "CL", "BR", "I", "F"} for a in component_atoms):
            required.append(_issue("metal_component", "Metal identity/coordination and the appropriate force-field model must be reviewed.", component_id=component_id, elements=sorted({a["element"] for a in component_atoms})))
        if component["kind"] == "protein":
            for left, right in zip(members, members[1:]):
                lc = [by_id[aid] for aid in left["atom_ids"] if by_id[aid]["atom_name"] == "C" and not by_id[aid]["altloc"]]
                rn = [by_id[aid] for aid in right["atom_ids"] if by_id[aid]["atom_name"] == "N" and not by_id[aid]["altloc"]]
                if len(lc) == len(rn) == 1 and lc[0]["xyz"] and rn[0]["xyz"]:
                    distance = math.dist(lc[0]["xyz"], rn[0]["xyz"])
                    if not 1.1 <= distance <= 1.8:
                        required.append(_issue("backbone_discontinuity", "Adjacent listed residues have an abnormal C-N distance; do not build a peptide bond automatically.", residue_ids=[left["id"], right["id"]], distance_angstrom=distance))
            if not any(r["record"] == "SEQRES" for r in records):
                warnings.append(_issue("no_full_sequence", "No SEQRES: residue-number continuity cannot establish construct completeness.", component_id=component_id,
                                       first_residue=members[0]["resseq"], last_residue=members[-1]["resseq"], residue_count=len(members)))

    for record in records:
        if record["record"] in {"LINK", "SSBOND"}:
            required.append(_issue("explicit_link_record", "Review explicit LINK/SSBOND connectivity before choosing the supported preparation route.", record=record["record"], line_number=record["line_number"], raw_line=record["raw_line"]))
    for index, item in enumerate(required, 1):
        item["id"] = "d%04d" % index
    return {"schema_version": SCHEMA_VERSION, "source": {"path": str(source), "filename": source.name,
                "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}, "path": str(source),
            "sha256": hashlib.sha256(data).hexdigest(), "coordinate_unit": "angstrom", "models": models,
            "atom_count": len(atoms), "residue_count": len(residues), "atoms": atoms, "residues": residues,
            "components": components, "bonds": bonds, "conect_records": conect, "structural_records": records,
            "sequence_records": [r for r in records if r["record"] == "SEQRES"],
            "has_seqres": any(r["record"] == "SEQRES" for r in records), "sequence_verified": False,
            "required_decisions": required, "warnings": warnings}


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def normalize_pdb(inventory, decisions, output_dir):
    """Write identifier-normalized, pose-preserving copies and a complete map.

    Decisions: ``model`` (inventory model id), ``altloc`` (label or residue->label),
    ``keep_components`` (explicit retained component ids). Omitted component
    selection retains everything. No chemistry or protonation decision is
    executed here. All excluded atoms retain identities/reasons in atom_map.json.
    """
    decisions = decisions or {}
    failures = []
    models = [m["id"] for m in inventory["models"]]
    model = decisions.get("model", models[0] if len(models) == 1 else None)
    if model not in models:
        failures.append(_issue("model", "Select one valid MODEL.", options=models, decision_key="model"))
    hard_types = {"invalid_atom_record", "invalid_conect", "duplicate_serial", "unresolved_conect", "self_bond", "duplicate_protein_atom_name"}
    failures.extend(item for item in inventory["required_decisions"] if item["type"] in hard_types and item.get("model", model) == model)
    alt_decisions = decisions.get("altloc", {})
    alt_by_residue = {}
    for residue in inventory["residues"]:
        if residue["model"] != model:
            continue
        labels = residue["altlocs"]
        label = alt_decisions if isinstance(alt_decisions, str) else alt_decisions.get(residue["id"])
        if len(labels) <= 1 and label is None:
            label = labels[0] if labels else ""
        if labels and label not in labels:
            failures.append(_issue("altloc", "Choose one available altLoc for the whole residue.", residue_id=residue["id"], options=labels, decision_key="altloc." + residue["id"]))
        alt_by_residue[residue["id"]] = label
    available = {c["id"] for c in inventory["components"] if c["model"] == model}
    keep = decisions.get("keep_components", sorted(available))
    if not isinstance(keep, list) or not set(keep) <= available:
        failures.append(_issue("keep_components", "Retained component ids must belong to the selected model.", options=sorted(available)))
    if failures:
        raise InventoryDecisionRequired(failures)
    keep = set(keep)
    selected = [a for a in inventory["atoms"] if a["model"] == model and a["component_id"] in keep and
                (not a["altloc"] or a["altloc"] == alt_by_residue[a["residue_id"]])]
    if not selected:
        raise InventoryDecisionRequired([_issue("empty_selection", "Selection retains no atoms.")])
    if len(selected) > 99999:
        raise ValueError("PDB serial limit exceeded; an mmCIF adapter is required, never truncate atom ids")
    selected_ids = {a["id"] for a in selected}
    selected_components = [c for c in inventory["components"] if c["id"] in keep]
    # Preserve nonblank original chains, allocate unoccupied legal identifiers.
    used_chains = {a["chain"] for a in selected if a["chain"].strip()}
    chain_for_component = {}
    blank_kind_chains = {}
    seen_protein_chains = set()
    protein_chain_replacements = {}

    def allocate_chain():
        free = next((ch for ch in CHAIN_ALPHABET if ch not in used_chains), None)
        if free is None:
            raise ValueError("No single-character PDB chain id remains; mmCIF adapter required")
        used_chains.add(free)
        return free

    for component in selected_components:
        chains = sorted({a["chain"] for a in component["atoms"] if a["id"] in selected_ids})
        if component["kind"] == "protein" and chains != [" "] and chains != [""]:
            _chain = chains[0]
            if _chain in seen_protein_chains:
                protein_chain_replacements[component["id"]] = allocate_chain()
            seen_protein_chains.add(_chain)
        if " " in chains or "" in chains:
            # Water/ion residues share a legal chain; allocating one chain per
            # water would fail on ordinary crystal waters after only 62 records.
            group = component["id"] if component["kind"] == "protein" else component["kind"]
            if group not in blank_kind_chains:
                blank_kind_chains[group] = allocate_chain()
            chain_for_component[component["id"]] = blank_kind_chains[group]
    residue_by_id = {r["id"]: r for r in inventory["residues"]}
    component_by_id = {c["id"]: c for c in selected_components}
    normalized_residues, residue_mapping, used_residue_keys = [], {}, set()
    variants = []
    normalization_decisions = []
    selected_by_residue = defaultdict(list)
    for atom in selected:
        selected_by_residue[atom["residue_id"]].append(atom)
    ligand_number = 0
    ligand_names = {}
    for component in selected_components:
        if component["kind"] == "ligand":
            ligand_number += 1
            if ligand_number > 99:
                raise ValueError("More than 99 ligand components require a broader naming adapter")
            ligand_names[component["id"]] = "L%02d" % ligand_number
    for atom in selected:
        rid = atom["residue_id"]
        if rid in residue_mapping:
            continue
        old = residue_by_id[rid]
        component = component_by_id[atom["component_id"]]
        chain = protein_chain_replacements.get(component["id"], atom["chain"] if atom["chain"].strip() else chain_for_component[component["id"]])
        resseq, icode = old["resseq"], old["icode"]
        if resseq is None or not -999 <= resseq <= 9999 or (chain, resseq, icode) in used_residue_keys:
            resseq = next((n for n in range(1, 10000) if (chain, n, "") not in used_residue_keys), None)
            icode = ""
            if resseq is None:
                raise ValueError("PDB residue numbering capacity exceeded")
        used_residue_keys.add((chain, resseq, icode))
        resname = old["resname"]
        if component["kind"] == "ligand":
            resname = ligand_names[component["id"]]
        elif resname == "HIS" or resname in HISTIDINE_VARIANTS:
            resname = "HIS"
            named_variant = HISTIDINE_VARIANTS.get(old["resname"])
            ring_h = {a["atom_name"] for a in selected_by_residue[rid] if a["element"] in {"H", "D"} and a["atom_name"] in {"HD1", "HE2"}}
            hydrogen_variant = ({frozenset({"HD1"}): "HID", frozenset({"HE2"}): "HIE", frozenset({"HD1", "HE2"}): "HIP"}).get(frozenset(ring_h))
            explicit_variant = decisions.get("histidine_variants", {}).get(rid)
            variant = explicit_variant or named_variant or hydrogen_variant
            if explicit_variant not in (None, "HID", "HIE", "HIP"):
                normalization_decisions.append(_issue("histidine_variant", "Choose a supported explicit histidine variant.", residue_id=rid, options=["HID", "HIE", "HIP"], decision_key="histidine_variants." + rid))
            elif named_variant and hydrogen_variant and named_variant != hydrogen_variant and explicit_variant is None:
                normalization_decisions.append(_issue("histidine_variant_conflict", "Residue-name protonation and explicit ring hydrogens disagree; explicitly choose the state before hydrogens are rebuilt.",
                                                     residue_id=rid, named_variant=named_variant, hydrogen_variant=hydrogen_variant,
                                                     ring_hydrogens=sorted(ring_h), options=["HID", "HIE", "HIP"], decision_key="histidine_variants." + rid))
            if variant:
                variants.append({"chain": chain, "resid": str(resseq), "insertion_code": icode,
                                 "variant": variant, "residue_id": rid, "original_residue_id": rid,
                                 "source": "explicit_decision" if explicit_variant else ("residue_alias" if named_variant else "explicit_ring_hydrogens"),
                                 "input_ring_hydrogens": sorted(ring_h)})
        entry = {"id": rid, "chain": chain, "resid": str(resseq), "resseq": resseq, "insertion_code": icode,
                 "name": resname, "resname": resname, "kind": component["kind"], "component_id": component["id"],
                 "original": {key: old[key] for key in ("chain", "resseq", "icode", "resname", "model", "segment")}}
        residue_mapping[rid] = entry
        normalized_residues.append(entry)
    name_counters, used_names = Counter(), defaultdict(set)
    normalized_atoms = []
    for serial, atom in enumerate(selected, 1):
        residue = residue_mapping[atom["residue_id"]]
        name = atom["atom_name"]
        if residue["kind"] == "ligand" or not name:
            element = atom["element"].title()
            key = (atom["residue_id"], element)
            name_counters[key] += 1
            name = element + str(name_counters[key])
        if len(name) > 4 or name in used_names[atom["residue_id"]]:
            raise ValueError("Atom identifier cannot be represented uniquely in four PDB columns: " + atom["id"])
        used_names[atom["residue_id"]].add(name)
        normalized_atoms.append({**atom, "serial": serial, "name": name, "atom_name": name,
                                 "chain": residue["chain"], "resseq": residue["resseq"], "icode": residue["insertion_code"],
                                 "resname": residue["resname"], "altloc": "", "original_serial": atom["serial"]})
    normalized_by_id = {a["id"]: a for a in normalized_atoms}
    bonds = [b for b in inventory["bonds"] if set(b["atom_ids"]) <= selected_ids]

    def write_pdb(path, ids):
        lines = ["REMARK 900 NEXUS V2 IDENTIFIER NORMALIZATION; COORDINATES IN ANGSTROM",
                 "REMARK 900 ORIGINAL IDENTITIES AND SELECTIONS: atom_map.json"]
        prior = None
        for atom in normalized_atoms:
            if atom["id"] not in ids:
                continue
            if prior is not None and prior["component_id"] != atom["component_id"] and component_by_id[prior["component_id"]]["kind"] == "protein":
                lines.append("TER")
            name = atom["atom_name"]
            aligned = name.ljust(4) if len(atom["element"]) == 2 or len(name) == 4 or name[0].isdigit() else (" " + name.ljust(3))
            # Reuse original fixed coordinate fields, avoiding rounding changes.
            coord = atom["raw_line"].ljust(80)[30:54]
            line = (f"{atom['record']:<6}{atom['serial']:5d} {aligned} {atom['resname']:>3} {atom['chain']}"
                    f"{atom['resseq']:4d}{atom['icode'] or ' '}   {coord}{atom['occupancy_text']}{atom['bfactor_text']}"
                    f"          {atom['element']:>2}{atom['charge_text']}")
            lines.append(line)
            prior = atom
        if prior and component_by_id[prior["component_id"]]["kind"] == "protein":
            lines.append("TER")
        adjacency = defaultdict(list)
        for bond in bonds:
            a, b = bond["atom_ids"]
            if a in ids and b in ids:
                adjacency[normalized_by_id[a]["serial"]].extend([normalized_by_id[b]["serial"]] * bond["multiplicity"])
                adjacency[normalized_by_id[b]["serial"]].extend([normalized_by_id[a]["serial"]] * bond["multiplicity"])
        for serial, neighbors in sorted(adjacency.items()):
            for offset in range(0, len(neighbors), 4):
                lines.append("CONECT" + f"{serial:5d}" + "".join(f"{n:5d}" for n in neighbors[offset:offset + 4]))
        lines.append("END")
        path.write_text("\n".join(lines) + "\n", encoding="ascii")

    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    output_names = {"normalized.pdb", "protein.pdb", "atom_map.json", "normalization.json"}
    output_names.update("ligand_" + c["id"] + ".pdb" for c in selected_components if c["kind"] == "ligand")
    if any((out / name).exists() for name in output_names):
        raise FileExistsError("Normalization outputs already exist; use a new stage directory instead of overwriting")
    combined_path, protein_path = out / "normalized.pdb", out / "protein.pdb"
    write_pdb(combined_path, selected_ids)
    protein_ids = {a["id"] for a in normalized_atoms if component_by_id[a["component_id"]]["kind"] == "protein"}
    write_pdb(protein_path, protein_ids)
    ligand_paths = {}
    for component in selected_components:
        if component["kind"] == "ligand":
            ligand_path = out / ("ligand_" + component["id"] + ".pdb")
            write_pdb(ligand_path, set(component["atom_ids"]) & selected_ids)
            ligand_paths[component["id"]] = str(ligand_path)
    mapping = []
    for atom in inventory["atoms"]:
        normalized = normalized_by_id.get(atom["id"])
        reason = None
        if normalized is None:
            reason = "unselected_model" if atom["model"] != model else ("excluded_component" if atom["component_id"] not in keep else "unselected_altloc")
        mapping.append({"id": atom["id"], "retained": normalized is not None, "excluded_reason": reason,
                        "original": {key: atom[key] for key in ("line_number", "record", "model", "serial", "chain", "resseq", "icode", "altloc", "resname", "atom_name", "element", "xyz", "formal_charge", "charge_text", "residue_id", "component_id")},
                        "normalized": {key: normalized[key] for key in ("serial", "chain", "resseq", "icode", "altloc", "resname", "atom_name", "element", "xyz", "formal_charge")} if normalized else None})
    reread = inventory_pdb(combined_path)
    if len(reread["atoms"]) != len(normalized_atoms) or any(a["xyz"] != b["xyz"] for a, b in zip(reread["atoms"], normalized_atoms)):
        raise RuntimeError("Normalization atom-count/coordinate invariant failed")
    if any((a["element"], a["formal_charge"]) != (b["element"], b["formal_charge"]) for a, b in zip(reread["atoms"], normalized_atoms)):
        raise RuntimeError("Normalization element/formal-charge invariant failed")
    reread_ids = {a["id"]: a["serial"] for a in reread["atoms"]}
    expected_bonds = {(tuple(sorted(normalized_by_id[aid]["serial"] for aid in b["atom_ids"])), b["multiplicity"]) for b in bonds}
    observed_bonds = {(tuple(sorted(reread_ids[aid] for aid in b["atom_ids"])), b["multiplicity"]) for b in reread["bonds"]}
    if expected_bonds != observed_bonds:
        raise RuntimeError("Normalization CONECT/multiplicity invariant failed")
    pending = []
    for item in inventory["required_decisions"]:
        if item["type"] in {"model", "altloc"} or item["type"] in hard_types:
            continue
        if item.get("model", model) != model or ("component_id" in item and item["component_id"] not in keep):
            continue
        residue_ids = item.get("residue_ids", [item["residue_id"]] if "residue_id" in item else [])
        if residue_ids and not all(rid in residue_mapping for rid in residue_ids):
            continue
        pending.append(item)
    pending.extend(normalization_decisions)
    atom_map_path = out / "atom_map.json"
    _write_json(atom_map_path, {"schema_version": SCHEMA_VERSION, "source": inventory["source"], "coordinate_unit": "angstrom",
                              "decisions": decisions, "atoms": mapping, "residues": normalized_residues,
                              "inverse_serial_to_id": {str(a["serial"]): a["id"] for a in normalized_atoms}})
    filtered_components = []
    for component in selected_components:
        comp_atoms = [a for a in normalized_atoms if a["component_id"] == component["id"]]
        comp_ids = {a["id"] for a in comp_atoms}
        comp_bonds = [b for b in bonds if set(b["atom_ids"]) <= comp_ids]
        covered = {aid for b in comp_bonds for aid in b["atom_ids"]}
        filtered_components.append({**component, "atoms": comp_atoms, "atom_ids": [a["id"] for a in comp_atoms],
                                   "bonds": comp_bonds, "atom_count": len(comp_atoms),
                                   "heavy_atom_count": sum(a["element"] not in {"H", "D"} for a in comp_atoms),
                                   "hydrogen_count": sum(a["element"] in {"H", "D"} for a in comp_atoms),
                                   "total_formal_charge": sum(a["formal_charge"] or 0 for a in comp_atoms),
                                   "connectivity_complete": bool(comp_bonds) and covered == comp_ids,
                                   "conect_complete": bool(comp_bonds) and covered == comp_ids})
    result = {"schema_version": SCHEMA_VERSION, "normalized_path": str(combined_path), "protein_path": str(protein_path),
              "ligand_paths": ligand_paths, "atom_map_path": str(atom_map_path), "source_atom_count": inventory["atom_count"],
              "retained_atom_count": len(normalized_atoms), "excluded_atom_count": inventory["atom_count"] - len(normalized_atoms),
              "coordinates_unchanged": True, "coordinate_unit": "angstrom", "selected_model": model,
              "selected_components": sorted(keep), "kept_component_ids": sorted(keep), "selected_atom_ids": [a["id"] for a in normalized_atoms],
              "components": filtered_components, "atoms": normalized_atoms, "residues": normalized_residues,
              "variants": variants, "pending_chemical_decisions": pending, "pending_decisions": pending,
              "identifier_changes_only": True, "decisions": decisions}
    _write_json(out / "normalization.json", result)
    return result
