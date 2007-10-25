from ase import *

atoms = Atoms(positions=[(0, 0, 0),
                         (1, 0, 0),
                         (0, 1, 0),
                         (1, 1, 0),
                         (0, 2, 0),
                         (1, 2, 0),
                         (0.5, 0.5, 1)],
              symbols='H7',
              constraints=[FixAtoms(range(6))],
              calculator=LennardJones())

traj = PickleTrajectory('H.traj', 'w', atoms)
dyn = QuasiNewton(atoms, maxstep=0.2)
dyn.attach(traj.write)
dyn.run(fmax=0.01, steps=100)
print dyn.H[-3:,-3:]
