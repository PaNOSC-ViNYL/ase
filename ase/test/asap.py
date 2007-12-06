from ase import *
from ase.calculator import ASAP

a = Atoms(symbols='Cu2', positions=[(0, 0, 0), (0, 0, 2.7)],
           calculator=ASAP())
print distance(n2, 0, 1), n2.get_potential_energy()
QuasiNewton(n2).run(0.01)
print distance(n2, 0, 1), n2.get_potential_energy()
