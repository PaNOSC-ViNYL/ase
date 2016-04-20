import doctest

try:
    import scipy
except ImportError:
    scipy = None
    
from ase import atoms
from ase.collection import collection
from ase.spacegroup import spacegroup, findsym, xtal
from ase.geometry import geometry, cell
from ase.build import tools
from ase.io import aff

modules = [collection, xtal, spacegroup, cell, findsym, tools, aff, atoms]
if scipy:
    modules.append(geometry)
    
for mod in modules:
    print(mod, doctest.testmod(mod, raise_on_error=True))
