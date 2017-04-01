from ase import Atoms
from ase.constraints import FixAtoms
from ase.build import molecule, bulk
from ase.io.trajectory import Trajectory
from ase.visualize import view

a1 = molecule('H2O')
a2 = molecule('H2O')
a2.center(vacuum=2.0)
a2.rattle(stdev=0.2)
a3 = molecule('CH3CH2OH')
a3.set_constraint(FixAtoms([0, 1]))
a4 = bulk('Au')

images = [a1, a2, a3, a4]

traj = Trajectory('out.traj', 'w')
print('a1')
traj.write(a1)
print('a2')
traj.write(a2)
print('a3')
traj.write(a3)
print('a4')
traj.write(a4)
print('done')
traj.close()

#view(images)

rtraj = Trajectory('out.traj')
newimages = list(rtraj)

assert len(images) == len(newimages)
for i in range(len(images)):
    assert images[i] == newimages[i], i
