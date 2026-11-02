"""CI contracts for actual exported exception parameters; no MD simulation."""
import pytest


def test_explicit_pairs_roundtrip_and_corruption(tmp_path):
    mm = pytest.importorskip('openmm')
    parmed = pytest.importorskip('parmed')
    from openmm import app, unit
    from scientific.topology_conversion import verify_nonbonded_conversion
    top = app.Topology()
    residue = top.addResidue('MOL', top.addChain())
    atoms = [top.addAtom('C'+str(i), app.element.carbon, residue) for i in range(4)]
    system = mm.System()
    force = mm.NonbondedForce()
    for q in [.2, -.2, -.2, .2]:
        system.addParticle(12)
        force.addParticle(q, .34, .4)
    bonds = mm.HarmonicBondForce()
    for i in range(3):
        top.addBond(atoms[i], atoms[i+1])
        bonds.addBond(i, i+1, .15, 1000)
    force.createExceptionsFromBonds([(0, 1), (1, 2), (2, 3)], 1/1.2, .5)
    system.addForce(force)
    system.addForce(bonds)
    positions = [mm.Vec3(0, 0, 0), mm.Vec3(.15, 0, 0), mm.Vec3(.2, .14, 0), mm.Vec3(.3, .2, .1)]*unit.nanometer
    structure = parmed.openmm.load_topology(top, system, xyz=positions)
    structure.save(str(tmp_path/'system.top'))
    converted = parmed.load_file(str(tmp_path/'system.top'))
    assert converted.defaults.gen_pairs == 'no'
    assert converted.defaults.fudgeLJ == 1
    result = converted.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    assert verify_nonbonded_conversion(system, result)['exception_pairs'] == 6
    copied = next(f for f in result.getForces() if isinstance(f, mm.NonbondedForce))
    for i in range(copied.getNumExceptions()):
        a, b, q, sigma, epsilon = copied.getExceptionParameters(i)
        if {a, b} == {0, 3}:
            copied.setExceptionParameters(i, a, b, q, sigma, epsilon*2)
    with pytest.raises(ValueError, match='exception parameters'):
        verify_nonbonded_conversion(system, result)


def test_hydrogen_addition_is_opt_in():
    pytest.importorskip('openmm')
    from scientific.receptor_hydrogens import prepare_hydrogens
    assert prepare_hydrogens(None, None, {}, 17)['requested'] is False
    with pytest.raises(ValueError, match='boolean'):
        prepare_hydrogens(None, None, {'add_receptor_hydrogens': 'true'}, 17)
    for ph in [None, True, float('nan'), -1, 15]:
        with pytest.raises(ValueError, match='pH'):
            prepare_hydrogens(None, None, {'add_receptor_hydrogens': True, 'receptor_pH': ph}, 17)
