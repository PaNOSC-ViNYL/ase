from ase.calculators.emt import EMT
from ase.lattice import bulk
from ase.optimize import FIRE

a = bulk('Au')
a *= (2,2,2)

a[0].x += 0.5

a.set_calculator(EMT())

opt = FIRE(a, dtmax=10.0, dt=10.0, maxmove=100.0, downhill_check=False)
opt.run(fmax=0.001)
e1 = a.get_potential_energy()
n1 = opt.nsteps

a = bulk('Au')
a *= (2,2,2)

a[0].x += 0.5

a.set_calculator(EMT())

opt = FIRE(a, dtmax=10.0, dt=10.0, maxmove=100.0, downhill_check=True)
opt.run(fmax=0.001)
e2 = a.get_potential_energy()
n2 = opt.nsteps

assert abs(e1-e2) < 1e-6
assert n2 < n1