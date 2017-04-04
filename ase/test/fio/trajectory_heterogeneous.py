import numpy as np
from ase.constraints import FixAtoms
from ase.build import molecule, bulk
from ase.io.trajectory import Trajectory

a1 = molecule('H2O')
a2 = molecule('H2O')
a2.center(vacuum=2.0)
a2.rattle(stdev=0.2)
a3 = molecule('CH3CH2OH')
a4 = bulk('Au').repeat((2, 2, 2))
a5 = bulk('Cu').repeat((2, 2, 3))

images = [a1, a2, a3, a4, a5]
for i, img in enumerate(images):
    img.set_constraint(FixAtoms(indices=range(i)))

traj = Trajectory('out.traj', 'w')
for img in images:
    traj.write(img)
traj.close()

#view(images)

rtraj = Trajectory('out.traj')
newimages = list(rtraj)

assert len(images) == len(newimages)
for i in range(len(images)):
    assert images[i] == newimages[i], i
