"""Explicit, recorded hydrogen completion for independently supplied receptors."""
import math
import random


def prepare_hydrogens(model, forcefield, parameters, seed):
    from openmm import unit
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
    variants = model.addHydrogens(forcefield, pH=ph)
    after = heavy_atoms()
    if len(before) != len(after) or any(k != k2 or any(abs(x-y) > 1e-8 for x, y in zip(p, p2))
                                      for (k, p), (k2, p2) in zip(before, after)):
        raise ValueError('Receptor hydrogen addition changed heavy atoms or their coordinates')
    return {'method': 'OpenMM Modeller.addHydrogens', 'requested': True, 'pH': ph,
            'seed': seed, 'atoms_before': atom_count, 'atoms_after': model.topology.getNumAtoms(),
            'selected_variants': variants, 'heavy_atoms_preserved': True,
            'review_required_before_MD': True}
