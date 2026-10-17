"""Deterministic, graph-preserving pre-docking geometry (not a docking score).

Coordinates are Angstroms. Hypotheses contain ``hypothesis_id``, ``anchor_points``
(``xyz``, ``atom_type``, optional ``radius``/``required``), and
``excluded_volumes`` (``xyz``, ``radius``). Receptor atoms are dictionaries with
``xyz``, ``element``, ``atom_name``, ``resname``, ``chain``, ``resid`` (or an
explicit ``residue_id``). Hydrogen bonds here are distance/type screens only;
no directionality or energetic claim is made.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures

_FACTORY = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))


def deterministic_id(kind, *parts):
    payload = json.dumps(parts, sort_keys=True, separators=(",", ":"), default=str)
    return kind + "_" + hashlib.sha256(payload.encode()).hexdigest()[:20]


def generate_conformers(smiles_or_mol, compound_id=None, state_id=None,
                        num_confs=20, seed=20261017, max_iters=200):
    """Return (hydrogenated graph with conformers, auditable generation metadata).

    Embedding shortfalls, optimizer nonconvergence and force-field failures are
    retained, never silently represented as successful conformers.
    """
    if num_confs < 1 or not 0 <= int(seed) <= 2147483647:
        raise ValueError("num_confs must be positive; seed must be a signed positive int")
    mol = Chem.MolFromSmiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else Chem.Mol(smiles_or_mol)
    if mol is None:
        raise ValueError("invalid molecular graph")
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    compound_id = compound_id or deterministic_id("compound", smiles)
    state_id = state_id or deterministic_id("state", compound_id, smiles)
    mol = Chem.AddHs(mol)
    mol.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.numThreads = 1
    params.trackFailures = True
    params.pruneRmsThresh = -1.0
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(num_confs), params=params))
    meta = dict(compound_id=compound_id, state_id=state_id, canonical_smiles=smiles,
                method="ETKDGv3", seed=int(seed), requested_count=int(num_confs),
                generated_count=len(ids), failed_count=int(num_confs) - len(ids),
                embedding_failure_counts=list(params.GetFailureCounts()), conformers=[], failures=[])
    for index, cid in enumerate(ids):
        record = dict(rdkit_conf_id=int(cid), conformer_id=deterministic_id("conformer", state_id, seed, index),
                      force_field=None, optimization_status=None, energy=None)
        errors = []
        for name in ("MMFF", "UFF"):
            try:
                if name == "MMFF":
                    props = AllChem.MMFFGetMoleculeProperties(mol)
                    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid) if props else None
                else:
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid) if AllChem.UFFHasAllMoleculeParams(mol) else None
                if ff is None:
                    raise ValueError("missing force-field parameters")
                status = int(ff.Minimize(maxIts=int(max_iters)))
                record.update(force_field=name, optimization_status=status, energy=float(ff.CalcEnergy()))
                if status:
                    errors.append(name + ": optimizer did not converge")
                break
            except (ValueError, RuntimeError) as exc:
                errors.append(name + ": " + str(exc))
        record["failures"] = errors
        if errors:
            meta["failures"].append(dict(conformer_id=record["conformer_id"], errors=errors))
        meta["conformers"].append(record)
    mol.SetProp("compound_id", compound_id)
    mol.SetProp("state_id", state_id)
    return mol, meta


def _xyz(record):
    value = record.get("xyz", record.get("position", record.get("center", record.get("coords"))))
    if value is None and all(k in record for k in ("x", "y", "z")):
        value = [record[k] for k in ("x", "y", "z")]
    result = np.asarray(value, dtype=float)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ValueError("expected finite xyz coordinate")
    return result


def atom_features(mol):
    """Per-atom feature sets from RDKit's chemically aware feature definitions."""
    result = {a.GetIdx(): set() for a in mol.GetAtoms()}
    aliases = {"Donor": "donor", "Acceptor": "acceptor", "NegIonizable": "negative",
               "PosIonizable": "positive", "Aromatic": "aromatic", "Hydrophobe": "hydrophobic",
               "LumpedHydrophobe": "hydrophobic"}
    for feature in _FACTORY.GetFeaturesForMol(mol):
        if feature.GetFamily() in aliases:
            for idx in feature.GetAtomIds():
                result[idx].add(aliases[feature.GetFamily()])
    return result


def _feature_type(point):
    value = str(point.get("feature_type", point.get("atom_type", point.get("type", "any")))).lower()
    return {"hba": "acceptor", "hbd": "donor", "o": "acceptor", "oxygen": "acceptor",
            "negionizable": "negative", "posionizable": "positive", "anion": "negative",
            "cation": "positive", "hydrophobe": "hydrophobic"}.get(value, value)


def _compatible(point, atom, features):
    kind = _feature_type(point)
    return kind in ("any", "heavy") or kind in features or kind == atom.GetSymbol().lower()


def _anchors(hypothesis):
    return hypothesis.get("anchor_points", hypothesis.get("features", []))


def _assign_features(mol, hypothesis, conf_id=-1, radius=2.5):
    coords = mol.GetConformer(conf_id).GetPositions()
    features = atom_features(mol)
    points = _anchors(hypothesis)
    options = []
    for p in points:
        candidates = [(float(np.linalg.norm(coords[a.GetIdx()] - _xyz(p))), a.GetIdx())
                      for a in mol.GetAtoms() if a.GetAtomicNum() > 1 and _compatible(p, a, features[a.GetIdx()])]
        options.append([(d, i) for d, i in sorted(candidates) if d <= float(p.get("radius", radius))])
    # Maximum-cardinality bipartite matching avoids counting one atom twice and
    # avoids greedy nearest-neighbor failures when two features share candidates.
    owner = {}
    def augment(pidx, seen):
        for _, idx in options[pidx]:
            if idx in seen:
                continue
            seen.add(idx)
            if idx not in owner or augment(owner[idx], seen):
                owner[idx] = pidx
                return True
        return False
    for pidx in sorted(range(len(points)), key=lambda k: (not points[k].get("required", True), len(options[k]))):
        augment(pidx, set())
    matches = [dict(feature_index=pidx, atom_index=idx,
                    distance=float(np.linalg.norm(coords[idx] - _xyz(points[pidx]))))
               for idx, pidx in owner.items()]
    return sorted(matches, key=lambda m: m["feature_index"])


def align_to_hypothesis(mol, hypothesis, conf_id=-1, max_assignments=256):
    """Rigidly align compatible atom assignments; retain the original graph.

    Bounded deterministic assignment enumeration is a screening heuristic, not
    exhaustive pharmacophore search. No receptor coordinates are moved.
    """
    if not _anchors(hypothesis):
        raise ValueError("hypothesis has no anchor points")
    result = Chem.Mol(mol)
    coords = np.array(mol.GetConformer(conf_id).GetPositions())
    features = atom_features(mol)
    points, choices = [], []
    for point in _anchors(hypothesis):
        candidates = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1 and _compatible(point, a, features[a.GetIdx()])]
        if candidates:
            points.append(_xyz(point))
            choices.append(candidates)
    if not choices:
        return result
    target = np.array(points)
    best = None
    for assignment in itertools.islice(itertools.product(*choices), int(max_assignments)):
        if len(set(assignment)) != len(assignment):
            continue
        source = coords[list(assignment)]
        sc, tc = source.mean(axis=0), target.mean(axis=0)
        u, _, vt = np.linalg.svd((source - sc).T @ (target - tc))
        rotation = u @ np.diag([1., 1., np.linalg.det(u @ vt)]) @ vt
        transformed = (coords - sc) @ rotation + tc
        rmsd = float(np.sqrt(np.mean(np.sum((transformed[list(assignment)] - target) ** 2, axis=1))))
        if best is None or rmsd < best[0]:
            best = rmsd, transformed
    if best is not None:
        conf = result.GetConformer(conf_id)
        for i, xyz in enumerate(best[1]):
            conf.SetAtomPosition(i, xyz.tolist())
    return result


def _receptor_types(atom):
    if "donor" in atom or "acceptor" in atom:
        return bool(atom.get("donor")), bool(atom.get("acceptor"))
    name = str(atom.get("atom_name", atom.get("name", ""))).upper().strip()
    res = str(atom.get("resname", atom.get("residue_name", ""))).upper()
    element = str(atom.get("element", name[:1])).upper()
    donor = element == "N" and not (res == "PRO" and name == "N")
    donor |= (res, name) in {(r, n) for r, n in [("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("CYS", "SG")]}
    acceptor = element == "O" or (res, name) in {("HIS", "ND1"), ("HIS", "NE2"), ("CYS", "SG"), ("MET", "SD")}
    return donor, acceptor


def _residue_key(atom):
    if atom.get("residue_id") is not None:
        return str(atom["residue_id"])
    resid = atom.get("resid", atom.get("residue_number", atom.get("resseq")))
    if resid is None:
        raise ValueError("receptor atoms need residue_id or residue number for unique anchor counting")
    return ":".join(str(x) for x in (atom.get("chain", atom.get("chain_id", "")), resid,
                                    atom.get("icode", ""), atom.get("resname", atom.get("residue_name", ""))))


def evaluate_pose(mol, hypothesis, receptor_atoms=(), conf_id=-1,
                  hard_clash_A=1.8, anchor_cutoff_A=4.0, feature_radius_A=2.5):
    """Evaluate typed contacts and apply independent, non-overridable hard gates."""
    coords = mol.GetConformer(conf_id).GetPositions()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if not heavy or not np.isfinite(coords).all():
        raise ValueError("pose must have finite heavy-atom coordinates")
    receptor = [a for a in receptor_atoms if str(a.get("element", "")).upper() != "H"]
    minimum, clashes, hits, hbonds = None, 0, set(), []
    features = atom_features(mol)
    for atom in receptor:
        distances = np.linalg.norm(coords[heavy] - _xyz(atom), axis=1)
        distance = float(distances.min())
        minimum = distance if minimum is None else min(minimum, distance)
        clashes += int(np.count_nonzero(distances < hard_clash_A))
        donor, acceptor = _receptor_types(atom)
        for idx, dist in zip(heavy, distances):
            compatible = (donor and "acceptor" in features[idx]) or (acceptor and "donor" in features[idx])
            if compatible and dist <= anchor_cutoff_A:
                key = _residue_key(atom)
                hits.add(key)
                hbonds.append(dict(residue_id=key, ligand_atom_index=idx, distance=float(dist)))
    excluded = 0
    for volume in hypothesis.get("excluded_volumes", hypothesis.get("exclusion_volumes", [])):
        radius = float(volume.get("radius", volume.get("radius_A", hard_clash_A)))
        excluded += int(np.count_nonzero(np.linalg.norm(coords[heavy] - _xyz(volume), axis=1) < radius))
    matches = _assign_features(mol, hypothesis, conf_id, feature_radius_A)
    matched = {m["feature_index"] for m in matches}
    required = {i for i, p in enumerate(_anchors(hypothesis)) if p.get("required", True)}
    feature_pass = bool(_anchors(hypothesis)) and required.issubset(matched)
    hard_pass = clashes == 0 and excluded == 0
    reasons = []
    if clashes:
        reasons.append("receptor_clash")
    if excluded:
        reasons.append("excluded_volume")
    if not feature_pass:
        reasons.append("missing_required_features")
    return dict(hypothesis_id=hypothesis.get("hypothesis_id"), clash_count=clashes,
                exclusion_count=excluded, min_receptor_distance=minimum,
                anchor_rmsd=float(np.sqrt(np.mean([m["distance"] ** 2 for m in matches]))) if matches else None,
                feature_matches=matches, feature_hits=len(matches), anchor_hits=len(hits),
                anchor_residues=sorted(hits), hydrogen_bonds=hbonds,
                total_score=float(len(matches) + len(hits)), passed_hard_gate=hard_pass,
                passed_geometry=bool(hard_pass and feature_pass), failure_reasons=reasons)


def make_pose_record(mol, hypothesis, compound_id=None, state_id=None, conformer_id=None,
                     conf_id=-1, receptor_atoms=(), **evaluation_options):
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    compound_id = compound_id or (mol.GetProp("compound_id") if mol.HasProp("compound_id") else deterministic_id("compound", smiles))
    state_id = state_id or (mol.GetProp("state_id") if mol.HasProp("state_id") else deterministic_id("state", compound_id, smiles))
    coords = mol.GetConformer(conf_id).GetPositions().tolist()
    conformer_id = conformer_id or deterministic_id("conformer", state_id, coords)
    hypothesis_id = hypothesis.get("hypothesis_id") or deterministic_id("hypothesis", hypothesis)
    result = evaluate_pose(mol, hypothesis, receptor_atoms, conf_id, **evaluation_options)
    result.update(compound_id=compound_id, state_id=state_id, conformer_id=conformer_id,
                  hypothesis_id=hypothesis_id,
                  pose_id=deterministic_id("pose", compound_id, state_id, conformer_id, hypothesis_id, coords),
                  smiles=smiles, mol=Chem.Mol(mol), rdkit_conf_id=int(conf_id))
    return result


def write_pose_sdf(records, path):
    """Write each pose from its RDKit molecular graph, never PDB proximity bonds."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    with Chem.SDWriter(str(path)) as writer:
        for record in records:
            if record["pose_id"] in seen:
                raise ValueError("duplicate pose_id: " + record["pose_id"])
            seen.add(record["pose_id"])
            mol = Chem.Mol(record["mol"])
            mol.SetProp("_Name", record["pose_id"])
            for key, value in record.items():
                if key not in {"mol", "rdkit_conf_id"} and value is not None:
                    mol.SetProp(key, json.dumps(value, sort_keys=True) if isinstance(value, (dict, list, bool)) else str(value))
            writer.write(mol, confId=record.get("rdkit_conf_id", -1))
    return path


def load_hypotheses(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return data["hypotheses"] if isinstance(data, dict) else data
