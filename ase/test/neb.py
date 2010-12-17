from ase.io import PickleTrajectory
from ase.neb import NEB
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS

print [a.get_potential_energy() for a in PickleTrajectory('H.traj')]
images = [PickleTrajectory('H.traj')[-1]]
for i in range(4):
    images.append(images[0].copy())
images[-1].positions[6, 1] = 2 - images[0].positions[6, 1]
neb = NEB(images)
neb.interpolate()

for image in images:
    image.set_calculator(LennardJones())

for a in neb.images:
    print a.positions[-1], a.get_potential_energy()

dyn = BFGS(neb, trajectory='mep.traj')
print dyn.run(fmax=0.01, steps=25)
for a in neb.images:
    print a.positions[-1], a.get_potential_energy()
