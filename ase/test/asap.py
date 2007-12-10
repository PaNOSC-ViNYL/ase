from ase import *
from ase.calculators.emt import ASAP

a = Atoms(symbols='Cu2', positions=[(0, 0, 0), (0, 0, 2.7)],
           calculator=ASAP())
print distance(a, 0, 1), a.get_potential_energy()
QuasiNewton(a).run(0.01)
print distance(a, 0, 1), a.get_potential_energy()
