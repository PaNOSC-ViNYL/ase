import numpy as np

from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.optimize.fire import FIRE as QuasiNewton
from ase.io import read

# Read the previous configurations
initial = read('N2.traj')
final = read('2N.traj')

#  Make 9 images (note the use of copy)
configs = [initial.copy() for i in range(8)] + [final]

# As before, fix the Cu atoms
constraint = FixAtoms(mask=[atom.symbol != 'N' for atom in initial])
for config in configs:
    config.set_calculator(EMT())
    config.set_constraint(constraint)

# Make the NEB object, interpolate to guess the intermediate steps
band = NEB(configs)
band.interpolate()

relax = QuasiNewton(band)

# Do the calculation
relax.run()

# Compare intermediate steps to initial energy
e0 = initial.get_potential_energy()
for config in configs:
    d = config[-2].position - config[-1].position
    print(np.linalg.norm(d), config.get_potential_energy() - e0)
