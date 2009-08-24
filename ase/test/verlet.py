from ase import *
from ase.calculators import TestPotential
np.seterr(all='raise')
a = Atoms('4N', 
          positions=[(0, 0, 0),
                     (1, 0, 0),
                     (0, 1, 0),
                     (0.1, 0.2, 0.7)],
          calculator=TestPotential())
print a.get_forces()
md = VelocityVerlet(a, dt=0.0005, logfile='-', loginterval=500)
traj = PickleTrajectory('4N.traj', 'w', a)
md.attach(traj.write, 100)
#print md.observers
md.run(steps=10000)
qn = QuasiNewton(a)
qn.attach(traj.write)
qn.run()
