from ase import Atoms
from ase.constraints import FixAtoms

a = Atoms('H3')
a.constraints = FixAtoms(indices=[0, 1])
del a[:-1]
print(a.constraints)
a = Atoms('H3')
a.constraints = FixAtoms(indices=[0, 1])
del a[:1]
print(a.constraints)
a = Atoms('H3')
a.constraints = FixAtoms(indices=[0, 1])
del a[:]
print(a.constraints)
