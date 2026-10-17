"""2D state review: RDKit-generated SVG with real per-atom charges (plan §15). No generative images."""
from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# conditions that can never be bulk-approved
HARD_BLOCKS = {"CHARGE_ERROR", "EMBED_FAILED", "NOT_CONVERGED"}


def state_svg(smiles: str, charges: list | None = None, size=(360, 260)) -> str:
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(*size)
    opts = d.drawOptions(); opts.addAtomIndices = True; opts.addStereoAnnotation = True
    d.DrawMolecule(m)
    d.FinishDrawing()
    return d.GetDrawingText()


def bulk_approvable(state: dict, rule: dict) -> tuple[bool, str]:
    """A saved user rule can pass only anomaly-free states. Returns (ok, why_not)."""
    if state.get("status") == "BLOCKED" or state.get("block_reason") in HARD_BLOCKS:
        return False, "hard_error"
    if state.get("truncated"):
        return False, "truncated_enumeration"
    if state.get("stereo_note") == "enumerated_unknown_centers":
        return False, "unknown_stereo"
    ch = state.get("charges", {}).get("check", {})
    if not ch.get("pass"):
        return False, "charge_check"
    if abs(state.get("formal_charge", 0)) > rule.get("max_abs_formal_charge", 1):
        return False, "formal_charge_outside_rule"
    if rule.get("require_ph_confirmed", True) and not state.get("ph_confirmed"):
        return False, "ph_not_confirmed"
    return True, ""
