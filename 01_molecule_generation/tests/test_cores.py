from datetime import date

import pytest
from rdkit import Chem

from src.cores import (
    CoreError,
    attach_substituents,
    derive_design_core,
    statistical_murcko,
    validate_attachment_labels,
    validate_product,
)
from src.synthesis import (
    BuildingBlock,
    ReactionRule,
    RouteStatus,
    RouteStep,
    audit_route,
    check_identity,
    forward_reconstruct,
    load_building_blocks_csv,
    load_reactions_csv,
)


def test_correct_core_is_not_statistical_murcko_and_keeps_key_map():
    selected = derive_design_core("CCOc1cc([NH2:101])ccc1C", [101], seed_id="p2m-7")

    assert selected.statistical_murcko == "c1ccccc1"
    assert selected.design_core != selected.statistical_murcko
    assert selected.precursor != selected.design_core
    assert selected.attachment_ids == (1, 2)
    assert validate_attachment_labels(selected.design_core) == (1, 2)
    core = Chem.MolFromSmiles(selected.design_core)
    assert any(a.GetAtomicNum() == 7 and a.GetAtomMapNum() == 101 for a in core.GetAtoms())


def test_wrong_attachment_is_rejected_and_expected_stars_are_consumed():
    selected = derive_design_core("CCOc1cc([NH2:101])ccc1C", [101])
    with pytest.raises(CoreError, match="do not match"):
        attach_substituents(selected.design_core, {1: "[*:1]F", 3: "[*:3]Cl"})

    product = attach_substituents(selected.design_core, {1: "[*:1]F", 2: "[*:2]Cl"})
    assert not any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms())
    assert validate_product(selected.design_core, product, [101])


@pytest.mark.parametrize("smiles", ["c1cc([*:1])ccc1[*:1]", "c1cc([*:1])ccc1[*:3]"])
def test_duplicate_or_discontinuous_dummy_labels_are_rejected(smiles):
    with pytest.raises(CoreError):
        validate_attachment_labels(smiles)


def test_stereo_mismatch_is_not_identity():
    result = check_identity("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O")
    assert not result.matches
    assert result.reason == "stereochemistry mismatch"


def test_stale_stock_triggers_review_and_prevents_formal_route():
    identity = check_identity("CO", "CO")
    step = RouteStep("alkylation", RouteStatus.SUPPORTED_FORWARD, "CO", ("MeOH",), identity)
    audit = audit_route(
        [step],
        [BuildingBlock("MeOH", "CO", "bb-v1", True, date(2024, 1, 1))],
        [ReactionRule("alkylation", "[C:1][Br:2].[O:3]>>[C:1][O:3]", "rxn-v1")],
        as_of=date(2024, 3, 1),
        max_stock_age_days=30,
    )
    assert audit.complete
    assert audit.needs_review
    assert not audit.formal
    assert "stale" in audit.review_reasons[0]


def test_forward_identity_with_mapped_reaction():
    rule = ReactionRule("alkylation", "[C:1][Br:2].[O:3]>>[C:1][O:3]", "rxn-v1")
    result = forward_reconstruct("CO", reaction=rule, reactants=["CBr", "O"])
    assert result.matches, result.reason


def test_versioned_csv_contracts(tmp_path):
    blocks_csv = tmp_path / "building_blocks.csv"
    blocks_csv.write_text(
        "contract_version,block_id,smiles,in_stock,checked_at,vendor\n"
        "bb-v1,MeOH,CO,true,2025-01-01,example\n",
        encoding="utf-8",
    )
    reactions_csv = tmp_path / "reactions.csv"
    reactions_csv.write_text(
        "contract_version,reaction_id,reaction_smarts,direction\n"
        'rxn-v1,alkylation,"[C:1][Br:2].[O:3]>>[C:1][O:3]",forward\n',
        encoding="utf-8",
    )
    blocks = load_building_blocks_csv(blocks_csv, expected_version="bb-v1")
    reactions = load_reactions_csv(reactions_csv, expected_version="rxn-v1")
    assert blocks[0].block_id == "MeOH"
    assert reactions[0].reaction_id == "alkylation"


def test_only_complete_supported_route_can_be_formal():
    block = BuildingBlock("MeOH", "CO", "bb-v1", True, date(2025, 1, 1))
    reaction = ReactionRule("r", "[C:1][Br:2].[O:3]>>[C:1][O:3]", "rxn-v1")
    unresolved = audit_route(
        [RouteStep("r", RouteStatus.UNRESOLVED, block_ids=("MeOH",))],
        [block],
        [reaction],
        as_of=date(2025, 1, 2),
    )
    assert not unresolved.complete
    assert not unresolved.formal
    assert unresolved.status == RouteStatus.UNRESOLVED


def test_murcko_rejects_acyclic_seed():
    with pytest.raises(CoreError, match="no Murcko"):
        statistical_murcko("[NH2:1]CCO")
