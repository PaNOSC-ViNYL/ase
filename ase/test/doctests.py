import doctest

from ase.lattice.spacegroup import spacegroup, cell, findsym, xtal
from ase.utils import geometry

for mod in [xtal, spacegroup, cell, findsym, geometry]:
    doctest.testmod(mod, raise_on_error=True)
