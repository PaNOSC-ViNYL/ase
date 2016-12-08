"""

Check the unit cell is handled correctly

"""

from ase.calculators.vasp import Vasp
from ase.build import molecule

# Molecules come with no unit cell

atoms = molecule('CH4')
calc = Vasp()

try:
    atoms.set_calculator(calc)    
    atoms.get_total_energy()
except ValueError:
    pass
else:
    raise AssertionError()
