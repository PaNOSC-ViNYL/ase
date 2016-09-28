#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Test suit for the CP2K ASE calulator.

http://www.cp2k.org
Author: Ole Schuett <ole.schuett@mat.ethz.ch>
"""

from __future__ import division, print_function
import numpy as np
import os

from ase.test import NotAvailable
from ase.build import bulk
from ase.constraints import UnitCellFilter
from ase.optimize import MDMin
from ase.calculators.cp2k import CP2K


def main():
    """Adopted from ase/test/stress.py"""

    if "ASE_CP2K_COMMAND" not in os.environ:
        raise NotAvailable('$ASE_CP2K_COMMAND not defined')

    # setup a Fist Lennard-Jones Potential
    inp = """&FORCE_EVAL
                  &MM
                    &FORCEFIELD
                      &SPLINE
                        EMAX_ACCURACY 500.0
                        EMAX_SPLINE    1000.0
                        EPS_SPLINE 1.0E-9
                      &END
                      &NONBONDED
                        &LENNARD-JONES
                          atoms Ar Ar
                          EPSILON [eV] 1.0
                          SIGMA [angstrom] 1.0
                          RCUT [angstrom] 10.0
                        &END LENNARD-JONES
                      &END NONBONDED
                      &CHARGE
                        ATOM Ar
                        CHARGE 0.0
                      &END CHARGE
                    &END FORCEFIELD
                    &POISSON
                      &EWALD
                        EWALD_TYPE none
                      &END EWALD
                    &END POISSON
                  &END MM
                &END FORCE_EVAL"""

    calc = CP2K(label="test_stress", inp=inp, force_eval_method="Fist")

    vol0 = 4 * 0.91615977036  # theoretical minimum
    a0 = vol0 ** (1 / 3)
    a = bulk('Ar', 'fcc', a=a0)
    a.calc = calc
    a.set_cell(np.dot(a.cell,
                      [[1.02, 0, 0.03],
                       [0, 0.99, -0.02],
                       [0.1, -0.01, 1.03]]),
               scale_atoms=True)

    a *= (1, 2, 3)
    a.rattle()
    sigma_vv = a.get_stress(voigt=False)
    # print(sigma_vv)
    # print(a.get_potential_energy() / len(a))
    vol = a.get_volume()

    # compare stress tensor with numeric derivative
    deps = 1e-5
    cell = a.cell.copy()
    for v in range(3):
        x = np.eye(3)
        x[v, v] += deps
        a.set_cell(np.dot(cell, x), scale_atoms=True)
        ep = a.calc.get_potential_energy(a, force_consistent=True)
        x[v, v] -= 2 * deps
        a.set_cell(np.dot(cell, x), scale_atoms=True)
        em = a.calc.get_potential_energy(a, force_consistent=True)
        s = (ep - em) / 2 / deps / vol
        # print(v, s, abs(s - sigma_vv[v, v]))
        assert abs(s - sigma_vv[v, v]) < 1e-7
    for v1 in range(3):
        v2 = (v1 + 1) % 3
        x = np.eye(3)
        x[v1, v2] = deps
        x[v2, v1] = deps
        a.set_cell(np.dot(cell, x), scale_atoms=True)
        ep = a.calc.get_potential_energy(a, force_consistent=True)
        x[v1, v2] = -deps
        x[v2, v1] = -deps
        a.set_cell(np.dot(cell, x), scale_atoms=True)
        em = a.calc.get_potential_energy(a, force_consistent=True)
        s = (ep - em) / deps / 4 / vol
        # print(v1, v2, s, abs(s - sigma_vv[v1, v2]))
        assert abs(s - sigma_vv[v1, v2]) < 1e-7

    # run a cell optimization, see if it finds back original crystal structure
    opt = MDMin(UnitCellFilter(a), dt=0.01, logfile=None)
    opt.run(fmax=0.1)
    # print(a.cell)
    for i in range(3):
        for j in range(3):
            x = np.dot(a.cell[i], a.cell[j])
            y = (i + 1) * (j + 1) * a0 ** 2 / 2
            if i != j:
                y /= 2
            # print(i, j, x, (x - y) / x)
            assert abs((x - y) / x) < 0.01

    print('passed test "stress"')


main()
# EOF
