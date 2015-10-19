import numpy as np
from ase.units import eV
from ase.calculators.calculator import LockedParameters

class BasisSet(LockedParameters):
    def __init__(self, functions='DZP', split_norm=0.5, energy_shift=0.150):
        LockedParameters.__init__(
                self,
                functions=functions,
                split_norm=split_norm,
                energy_shift=energy_shift,
                )

class DZP(BasisSet):
    def __init__(self, **kwargs):
        BasisSet.__init__(self, functions='DZP', **kwargs)

class DZ(BasisSet):
    def __init__(self, **kwargs):
        BasisSet.__init__(self, functions='DZ', **kwargs)

class SZ(BasisSet):
    def __init__(self, **kwargs):
        BasisSet.__init__(self, functions='SZ', **kwargs)

class SZP(BasisSet):
    def __init__(self, **kwargs):
        BasisSet.__init__(self, functions='SZP', **kwargs)
