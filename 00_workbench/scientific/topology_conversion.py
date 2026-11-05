"""Audit the actual nonbonded parameters, including explicit exception pairs."""
import math


def verify_nonbonded_conversion(original, converted):
    from openmm import NonbondedForce, unit

    def force(system):
        forces = [f for f in system.getForces() if isinstance(f, NonbondedForce)]
        if len(forces) != 1:
            raise ValueError('Expected exactly one nonbonded force')
        return forces[0]

    def close(a, b):
        # GROMACS topology text rounds explicit parameters to nine decimals.
        return math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-7)

    left, right = force(original), force(converted)
    if left.getNumParticles() != right.getNumParticles():
        raise ValueError('Topology conversion changed particle count')
    units = (unit.elementary_charge, unit.nanometer, unit.kilojoule_per_mole)
    for i in range(left.getNumParticles()):
        a, b = left.getParticleParameters(i), right.getParticleParameters(i)
        if not all(close(x.value_in_unit(u), y.value_in_unit(u)) for x, y, u in zip(a, b, units)):
            raise ValueError('Topology conversion changed nonbonded particle parameters')

    def exceptions(f):
        rows = {}
        for i in range(f.getNumExceptions()):
            a, b, q, sigma, epsilon = f.getExceptionParameters(i)
            key = tuple(sorted((int(a), int(b))))
            if key in rows:
                raise ValueError('Duplicate nonbonded exception')
            rows[key] = (q.value_in_unit(unit.elementary_charge**2),
                         sigma.value_in_unit(unit.nanometer),
                         epsilon.value_in_unit(unit.kilojoule_per_mole))
        return rows

    a, b = exceptions(left), exceptions(right)
    if a.keys() != b.keys():
        raise ValueError('Topology conversion changed nonbonded exclusion/exception pairs')
    for key, (q, sigma, epsilon) in a.items():
        q2, sigma2, epsilon2 = b[key]
        if not close(q, q2) or not close(epsilon, epsilon2) or (max(abs(epsilon), abs(epsilon2)) > 1e-7 and not close(sigma, sigma2)):
            raise ValueError('Topology conversion changed nonbonded exception parameters')
    return {'status': 'passed', 'particles': left.getNumParticles(), 'exception_pairs': len(a),
            'method': 'compare source and reloaded charges, LJ parameters and explicit exceptions',
            'relative_tolerance': 1e-6, 'absolute_tolerance': 1e-7}
