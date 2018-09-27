import numpy as np
from time import time
from ase import Atoms
from ase.build import fcc111
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.neb import NEB

# Build Pt(111) slab with six surface atoms and add oxygen adsorbate
initial = fcc111('Pt', size=(3, 2, 3), orthogonal=True)
initial.center(axis=2, vacuum=10)
oxygen = Atoms('O')
oxygen.translate(initial[7].position + (0., 0., 3.5))
initial.extend(oxygen)

# EMT potential
calc = EMT()
initial.set_calculator(EMT())

# Optimize initial state
opt = BFGS(initial)
opt.run(fmax=0.03)

# Move oxygen adsorbate to neighboring hollow site
final = initial.copy()
final[18].x += 2.8
final[18].y += 1.8

final.set_calculator(EMT())

opt = BFGS(final)
opt.run(fmax=0.03)

# NEB with five interior images
images = [initial]
for i in range(5):
    images.append(initial.copy())
images.append(final)

fmax = 0.03  # Same for NEB and optimizer

for i in range(1, len(images)-1):
    calc = EMT()
    images[i].set_calculator(calc)

# Dynamic NEB
neb = NEB(images, fmax=fmax, dynamic_relaxation=True)
neb.interpolate()

# Optimize and check walltime of dynamic NEB
time_dyn_1 = time()
opt = BFGS(neb)
opt.run(fmax=fmax)
time_dyn_2 = time()
time_dyn = time_dyn_2 - time_dyn_1

# Get potential energy of transition state
Emax_dyn = np.sort([image.get_potential_energy()
                    for image in images[1:-1]])[-1]

# Default NEB
neb = NEB(images, dynamic_relaxation=False)
neb.interpolate()

# Optimize and check walltime of default NEB
time_def_1 = time()
opt = BFGS(neb)
opt.run(fmax=fmax)
time_def_2 = time()
time_def = time_def_2 - time_def_1

# Get potential energy of transition state
Emax_def = np.sort([image.get_potential_energy()
                    for image in images[1:-1]])[-1]

# Check time for default and dynamic NEB implementations
assert((time_dyn - time_def) < 0)

# Assert reaction barriers are within 1 meV of each other
assert(abs(Emax_dyn - Emax_def) < 1e-3)
