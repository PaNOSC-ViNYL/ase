from ase import Atoms
from ase.calculators.tip3p import TIP3P, rOH, angleHOH
from ase.md import Langevin
import ase.units as units
from ase.io.trajectory import Trajectory
from ase.constraints import FixBondLengths, FixBondLengthsFast, FixBondLengthsWaterModel

TIP3PWaterModel = TIP3P

import numpy as np
import time
from ase.io import read

NSTEPS = 1000
SCALE = 200

cutoff = 4.0

# Set up water box at 20 deg C density
x = angleHOH * np.pi / 180 / 2
pos = [[0, 0, 0],
       [0, rOH * np.cos(x), rOH * np.sin(x)],
       [0, rOH * np.cos(x), -rOH * np.sin(x)]]
atoms = Atoms('OH2', positions=pos)

vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
atoms.set_cell((vol, vol, vol))
atoms.center()

N = 4
atoms = atoms.repeat((N, N, N))
atoms.set_pbc(True)

pairs = [(3 * i + j, 3 * i + (j + 1) % 3)
           for i in range(len(atoms) / 3)
           for j in [0, 1, 2]]

# Create atoms object with old constraints for reference
atoms_ref = atoms.copy()
atoms_ref.constraints = FixBondLengths(pairs)

# RATTLE-type constraints on O-H1, O-H2, H1-H2.
#atoms.constraints = FixBondLengthsFast(pairs)
atoms.constraints = FixBondLengthsWaterModel(pairs)

atoms.calc = TIP3PWaterModel(rc=cutoff)
atoms_ref.calc = TIP3P(rc=cutoff)

np.random.seed(123)
md = Langevin(atoms, 1 * units.fs, temperature=300 * units.kB,
              friction=0.01, logfile='C.log')
traj = Trajectory('C.traj', 'w', atoms)
md.attach(traj.write, interval=1)

start = time.time()
md.run(NSTEPS)
end = time.time()
Cversion = end-start
print "%d steps of C-MD took %.3fs (%.0f ms/step)" % (NSTEPS, Cversion, Cversion/NSTEPS*1000)
traj.close()

np.random.seed(123)
md_ref = Langevin(atoms_ref, 1 * units.fs, temperature=300 * units.kB,
                  friction=0.01, logfile='ref.log')
traj_ref = Trajectory('ref.traj', 'w', atoms_ref)
md_ref.attach(traj_ref.write, interval=1)
start = time.time()
md_ref.run(NSTEPS / SCALE)
end = time.time()
Pyversion = (end-start) * SCALE
print "%d steps of Py-MD took %.3fs (%.0f ms/step)" % (NSTEPS/SCALE, Pyversion/SCALE, Pyversion/NSTEPS*1000)
traj_ref.close()

# Compare trajectories
images = read('C.traj@:')
images_ref = read('ref.traj@:')
for img1, img2 in zip(images, images_ref):
    norm = np.linalg.norm(img1.get_positions() - img2.get_positions())
    print norm
    assert norm <1e-11

print "Speedup", Pyversion / Cversion
