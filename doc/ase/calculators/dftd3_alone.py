from ase.calculators.dftd3 import DFTD3
from ase.build import bulk

diamond = bulk('C')
d3 = DFTD3()
diamond.set_calculator(d3)
diamond.get_potential_energy()
