from ase.data import atomic_numbers, chemical_symbols, reference_states
from ase.units import *

from . import fcc
from . import hcp
from . import au

lattice = {'fcc': fcc.data,
           'hcp': hcp.data,
          }

element = {79: au.data, #Au
          }
