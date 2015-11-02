from ase.test import must_raise
from ase.calculators.calculator import LockedParameters

basis_set = LockedParameters(x=5, hello='hello')
basis_set['hello'] = 'goodbye'
with must_raise(KeyError):
    basis_set['g'] = 7
