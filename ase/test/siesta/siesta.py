from __future__ import print_function

import os
from os.path import join
import unittest
import numpy as np

from ase.calculators.siesta import Siesta
from ase.calculators.siesta.siesta import Specie
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import LockedParameters

from ase import Atoms

try:
    del os.environ['SIESTA']
except KeyError:
    pass

os.environ['SIESTA_PP_PATH'] = os.path.abspath(join(os.path.dirname(__file__), 'TestFiles'))
h = Atoms('H', [(0.0, 0.0, 0.0)])
co2 = Atoms('CO2', [(0.0, 0.0, 0.0), (-1.178, 0.0, 0.0), (1.178, 0.0, 0.0)])
ch4 = Atoms('CH4', np.array([
          [0.000000,  0.000000,  0.000000],
          [0.682793,  0.682793,  0.682793],
          [-0.682793, -0.682793,  0.68279],
          [-0.682793,  0.682793, -0.682793],
          [0.682793, -0.682793, -0.682793]]))

dirt_cheap_siesta = Siesta(
        label='test_label',
        mesh_cutoff=(70, 'Ry'),
        basis_set='SZ',
        DM_Tolerance=1e-3,
        )
test_path = os.path.abspath(join(os.path.dirname(__file__), 'tmp'))
if not os.path.exists(test_path):
    os.makedirs(test_path)

class SiestaTest(unittest.TestCase):

    def testDefaultConstruction(self):
        siesta = Siesta()
        self.assertIsInstance(siesta, FileIOCalculator)
        self.assertIsInstance(siesta.implemented_properties, tuple)
        self.assertIsInstance(siesta.default_parameters, dict)
        self.assertIsInstance(siesta.name, str)
        self.assertIsInstance(siesta.default_parameters, dict)

    def testSpecies(self):
        siesta = Siesta()
        atoms = ch4.copy()
        species, numbers = siesta.species(atoms)
        self.assertTrue(all(numbers==np.array([1, 2, 2, 2, 2])))
        siesta = Siesta(species=[Specie(symbol='C', tag=1)])
        species, numbers = siesta.species(atoms)
        self.assertTrue(all(numbers==np.array([1, 2, 2, 2, 2])))
        atoms.set_tags([0,0,0,1,0])
        species, numbers = siesta.species(atoms)
        self.assertTrue(all(numbers==np.array([1, 2, 2, 2, 2])))
        siesta = Siesta(species=[Specie(symbol='H', tag=1)])
        species, numbers = siesta.species(atoms)
        self.assertTrue(all(numbers==np.array([1, 2, 2, 3, 2])))

    def testCalculate(self):
        atoms = h.copy()
        os.chdir(test_path)
        atoms.set_calculator(dirt_cheap_siesta)
        #self.assertRaises(OSError, atoms.get_total_energy)

    def testWrite(self):
        os.chdir(test_path)
        atoms = h.copy()
        atoms.set_calculator(dirt_cheap_siesta)
        dirt_cheap_siesta.write_input(atoms, properties=['energy'])

        with open('siesta.fdf', 'r') as f:
            lines = f.readlines()
        self.assertTrue('DM.Tolerance  0.001\n' in lines)

    def tearDown(self):
        pass
        #if os.path.exists(test_path):
        #    cmd = 'rm %s/*'%test_path
        #    os.system(cmd)

if __name__=='__main__':
    unittest.main()
