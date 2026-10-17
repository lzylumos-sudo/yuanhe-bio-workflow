import json,math
import openmm as mm
from openmm import unit
system=mm.System(); system.addParticle(12)
force=mm.CustomExternalForce('0.5*k*(x*x+y*y+z*z)');force.addGlobalParameter('k',100);force.addParticle(0,[]);system.addForce(force)
integrator=mm.LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,0.002*unit.picosecond)
integrator.setRandomNumberSeed(20261008)
context=mm.Context(system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':'4'})
context.setPositions([[0.1,0,0]]);context.setVelocitiesToTemperature(300*unit.kelvin,20261008)
integrator.step(1000)
state=context.getState(getEnergy=True)
assert math.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
print(json.dumps({'openmm':mm.__version__,'platform':'CPU','steps':1000,'finite_energy':True,'scope':'installation smoke only'}))
