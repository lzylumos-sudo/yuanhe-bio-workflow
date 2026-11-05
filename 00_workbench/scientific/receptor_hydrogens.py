"""Explicit, recorded hydrogen completion for independently supplied receptors."""
import math
import random


def prepare_hydrogens(model, forcefield, parameters, seed):
    from openmm import Platform, unit
    enabled = parameters.get('add_receptor_hydrogens', False)
    if type(enabled) is not bool:
        raise ValueError('Explicit boolean required for receptor hydrogen addition')
    if not enabled:
        return {'method': 'preserve supplied receptor', 'requested': False}
    ph = parameters.get('receptor_pH')
    if isinstance(ph, bool) or not isinstance(ph, (int, float)) or not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError('Choose a finite receptor pH between 0 and 14')

    def heavy_atoms():
        positions = model.positions.value_in_unit(unit.nanometer)
        return [((a.residue.chain.index, a.residue.index, a.name, a.element.symbol), tuple(positions[a.index]))
                for a in model.topology.atoms() if a.element is not None and a.element.symbol != 'H']

    before = heavy_atoms()
    atom_count = model.topology.getNumAtoms()
    random.seed(seed)
    # addHydrogens reads all positions back from its minimization Context.
    # Use double precision here so fixed receptor atoms retain their coordinates.
    variants = model.addHydrogens(forcefield, pH=ph,
                                  platform=Platform.getPlatformByName('Reference'))
    after = heavy_atoms()
    if len(before) != len(after):
        raise ValueError(f'Receptor hydrogen addition changed heavy atom count: {len(before)} -> {len(after)}')
    for (key, initial), (new_key, final) in zip(before, after):
        if key != new_key:
            raise ValueError(f'Receptor hydrogen addition changed heavy atom identity: {key} -> {new_key}')
        if not all(math.isclose(x, y, abs_tol=1e-8, rel_tol=0) for x, y in zip(initial, final)):
            raise ValueError(f'Receptor hydrogen addition changed coordinates of {key}: {initial} -> {final} nm')
    return {'method': 'OpenMM Modeller.addHydrogens', 'requested': True, 'pH': ph,
            'seed': seed, 'atoms_before': atom_count, 'atoms_after': model.topology.getNumAtoms(),
            'selected_variants': variants, 'heavy_atoms_preserved': True, 'hydrogen_platform': 'Reference',
            'review_required_before_MD': True}
