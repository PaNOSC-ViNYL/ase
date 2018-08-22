from ase.io import read
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.optimize import BFGS
from ase.parallel import rank, size

initial = read('initial.traj')
final = read('final.traj')

constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])

images = [initial]
j = rank * 3 // size  # my image number
for i in range(3):
    image = initial.copy()
    if i == j:
        image.set_calculator(EMT())
    image.set_constraint(constraint)
    images.append(image)
images.append(final)

neb = NEB(images, parallel=True)
neb.interpolate()
qn = BFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)
