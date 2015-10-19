import unittest
from ase.calculators.calculator import LockedParameters

class ParametersTest(unittest.TestCase):
    def testDefaultConstruction(self):
        basis_set = LockedParameters(x=5, hello='hello')

        basis_set['hello'] = 'goodbye'
        self.assertRaises(KeyError, basis_set.__setitem__, 'g', 7)


if __name__=='__main__':
    unittest.main()
