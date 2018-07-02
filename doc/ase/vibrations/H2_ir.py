from ase.build import molecule
from ase import optimize
from ase.vibrations.infrared import InfraRed

from gpaw.cluster import Cluster
from gpaw import GPAW, FermiDirac

h = 0.22

atoms = Cluster(molecule('H2'))
atoms.minimal_box(3.5, h=h)

# relax the molecule
calc = GPAW(h=h, occupations=FermiDirac(width=0.1))
atoms.set_calculator(calc)
dyn = optimize.FIRE(atoms)
dyn.run(fmax=0.05)
atoms.write('relaxed.traj')

# finite displacement for vibrations
ir = InfraRed(atoms)
ir.run()
