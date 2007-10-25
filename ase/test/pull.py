from ase import *

Cu = Atoms(positions=[(0, 0, 0)],
           symbols='Cu',
           pbc=(1, 0, 0),
           calculator=EMT())
traj = PickleTrajectory('Cu.traj', 'w')
for a in linspace(2.0, 4.0, 20):
    Cu.set_cell([a, 1, 1])
    traj.write(Cu)
