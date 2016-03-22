import doctest

try:
    import scipy
except ImportError:
    scipy = None
    
from ase.lattice.spacegroup import spacegroup, cell, findsym, xtal
from ase.utils import geometry

modules = [xtal, spacegroup, cell, findsym]
if scipy:
    modules.append(geometry)
    
for mod in modules:
    doctest.testmod(mod, raise_on_error=True)
