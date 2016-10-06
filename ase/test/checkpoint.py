#! /usr/bin/env python

from __future__ import print_function

import os
import unittest

import numpy as np

import ase
from ase.lattice.cubic import Diamond

from ase.calculators.checkpoint import Checkpoint, CheckpointCalculator
from ase.build import bulk
from ase.calculators.lj import LennardJones


class TestCheckpoint(unittest.TestCase):

    def op1(self, a, m):
        a[1].position += m * np.array([0.1, 0.2, 0.3])
        return a

    def op2(self, a, m):
        a += ase.Atom('C', m * np.array([0.2, 0.3, 0.1]))
        return a, a.positions[0]

    def test_sqlite(self):
        print('test_single_file')

        try:
            os.remove('checkpoints.db')
        except OSError:
            pass

        CP = Checkpoint('checkpoints.db')
        a = Diamond('Si', size=[2, 2, 2])
        a = CP(self.op1)(a, 1.0)
        op1a = a.copy()
        a, ra = CP(self.op2)(a, 2.0)
        op2a = a.copy()
        op2ra = ra.copy()

        CP = Checkpoint('checkpoints.db')
        a = Diamond('Si', size=[2, 2, 2])
        a = CP(self.op1)(a, 1.0)
        self.assertEqual(a, op1a)
        a, ra = CP(self.op2)(a, 2.0)
        self.assertEqual(a, op2a)
        self.assert_(np.abs(ra - op2ra).max() < 1e-5)


class TestCheckpointCalculator(unittest.TestCase):

    def rattle_calc(self, atoms, calc):
        try:
            os.remove('checkpoints.db')
        except OSError:
            pass

        orig_atoms = atoms.copy()

        # first do a couple of calculations
        np.random.seed(0)
        atoms.rattle()
        cp_calc_1 = CheckpointCalculator(calc)
        atoms.set_calculator(cp_calc_1)
        e11 = atoms.get_potential_energy()
        f11 = atoms.get_forces()
        atoms.rattle()
        e12 = atoms.get_potential_energy()
        f12 = atoms.get_forces()

        # then re-read them from checkpoint file
        atoms = orig_atoms
        np.random.seed(0)
        atoms.rattle()
        cp_calc_2 = CheckpointCalculator(calc)
        atoms.set_calculator(cp_calc_2)
        e21 = atoms.get_potential_energy()
        f21 = atoms.get_forces()
        atoms.rattle()
        e22 = atoms.get_potential_energy()
        f22 = atoms.get_forces()

        self.assertAlmostEqual(e11, e21)
        self.assertAlmostEqual(e12, e22)
        self.assert_(np.abs(f11 - f21).max() < 1e-5)
        self.assert_(np.abs(f12 - f22).max() < 1e-5)

    def test_new_style_interface(self):
        calc = LennardJones()
        atoms = bulk('Cu')
        self.rattle_calc(atoms, calc)


# run tests manually, as calling unittest.main() interferes with ase.test.test()
suite = unittest.TestSuite()
for test_class in [TestCheckpoint, TestCheckpointCalculator]:
    tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
    suite.addTests(tests)
result = unittest.TextTestRunner().run(suite)
assert len(result.failures) == 0
