from ase import *
from gpaw import GPAW

calc = GPAW(h=.2, xc='PBE', txt='benzene.txt',
            nbands=35, eigensolver='cg',
            convergence={'bands': 30})
atoms = molecule('C6H6', calculator=calc)
atoms.center(vacuum=5.)
atoms.get_potential_energy()
calc.write('benzene.gpw', mode='all')
