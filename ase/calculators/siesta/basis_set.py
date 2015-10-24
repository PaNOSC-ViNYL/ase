import numpy as np
from ase.units import eV
from ase.calculators.calculator import LockedParameters

class BasisSet(LockedParameters):
    def __init__(self, size):
        LockedParameters.__init__(
                self,
                size=size,
                )

DZP = BasisSet(size='DZP')
DZ = BasisSet(size='DZ')
SZ = BasisSet(size='SZ')
SZP = BasisSet(size='SZP')

