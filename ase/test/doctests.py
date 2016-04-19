import doctest

try:
    import scipy
except ImportError:
    scipy = None
    
from ase.collection import collection
from ase.spacegroup import spacegroup, findsym, xtal
from ase.geometry import geometry, cell

modules = [collection, xtal, spacegroup, cell, findsym]
if scipy:
    modules.append(geometry)
    
for mod in modules:
    doctest.testmod(mod, raise_on_error=True)
