#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Test suit for the CP2K ASE calulator.

http://www.cp2k.org
Author: Ole Schütt <ole.schuett@mat.ethz.ch>
"""

from __future__ import division, print_function

import os
import numpy as np

from ase.structure import molecule
from ase.optimize import BFGS
from ase import units
from ase.atoms import Atoms
from ase.lattice import bulk
from ase.constraints import UnitCellFilter
from ase.optimize import MDMin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet

from cp2k import CP2K

#=============================================================================
def test_H2_LDA():
    calc = CP2K(label='test_H2_LDA')
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    energy = h2.get_potential_energy()
    energy_ref = -30.6989595886
    diff = abs((energy - energy_ref) / energy_ref)
    assert(diff < 1e-10)
    print('passed test "H2_LDA"')


#=============================================================================
def test_H2_PBE():
    calc = CP2K(xc='PBE', label='test_H2_PBE')
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    energy = h2.get_potential_energy()
    energy_ref = -31.5917284949
    diff = abs((energy - energy_ref) / energy_ref)
    assert(diff < 1e-10)
    print('passed test "H2_PBE"')


#=============================================================================
def test_H2_LS():
    inp = """&FORCE_EVAL
               &DFT
                 &QS
                   LS_SCF ON
                 &END QS
               &END DFT
             &END FORCE_EVAL"""
    calc = CP2K(label='test_H2_LS', inp=inp)
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    energy = h2.get_potential_energy()
    energy_ref = -30.6989581747
    diff = abs((energy - energy_ref) / energy_ref)
    assert(diff < 1e-10)
    print('passed test "H2_LS"')


#=============================================================================
def test_O2():
    calc = CP2K(label='test_O2', uks=True, cutoff=150*units.Rydberg,
                basis_set="SZV-MOLOPT-SR-GTH")
    o2 = molecule('O2', calculator=calc)
    o2.center(vacuum=2.0)
    energy = o2.get_potential_energy()
    energy_ref = -861.057011375
    diff = abs((energy - energy_ref) / energy_ref)
    assert(diff < 1e-10)
    print('passed test "O2"')


#=============================================================================
def test_restart():
    calc = CP2K()
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    h2.get_potential_energy()
    calc.write('test_restart')  # write a restart
    calc2 = CP2K(restart='test_restart')  # load a restart
    assert not calc2.calculation_required(h2, ['energy'])
    print('passed test "restart"')


#=============================================================================
def test_GeoOpt():
    calc = CP2K(label='test_H2_GOPT')
    atoms = molecule('H2', calculator=calc)
    atoms.center(vacuum=2.0)

    # Run Geo-Opt
    gopt = BFGS(atoms, logfile=None)
    gopt.run(fmax=1e-6)

    # check distance
    dist = atoms.get_distance(0, 1)
    dist_ref = 0.7245595
    assert( (dist-dist_ref)/dist_ref < 1e-7 )

    # check energy
    energy_ref = -30.7025616943
    energy = atoms.get_potential_energy()
    assert( (energy-energy_ref)/energy_ref < 1e-10 )
    print('passed test "H2_GEO_OPT"')


#=============================================================================
def test_MD():
    calc = CP2K(label='test_H2_MD')
    atoms = Atoms('HH', positions=[(0,0,0), (0,0,0.7245595)], calculator=calc)
    atoms.center(vacuum=2.0)

    # Run MD
    MaxwellBoltzmannDistribution(atoms, 0.5 * 300 * units.kB, force_temp=True)
    energy_start = atoms.get_potential_energy() + atoms.get_kinetic_energy()
    dyn = VelocityVerlet(atoms, 0.5 * units.fs)
    #def print_md():
    #    energy = atoms.get_potential_energy() + atoms.get_kinetic_energy()
    #    print("MD total-energy: %.10feV" %  energy)
    #dyn.attach(print_md, interval=1)
    dyn.run(20)

    energy_end = atoms.get_potential_energy() + atoms.get_kinetic_energy()

    assert(energy_start-energy_end < 1e-4)
    print('passed test "H2_MD"')


#=============================================================================
def test_stress():
    """Adopted from ase/test/stress.py"""

    # setup a Fist Lennard-Jones Potential
    inp= """&FORCE_EVAL
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
    a0 = vol0**(1 / 3)
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
    #print(sigma_vv)
    #print(a.get_potential_energy() / len(a))
    vol = a.get_volume()

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
        #print(v, s, abs(s - sigma_vv[v, v]))
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
        #print(v1, v2, s, abs(s - sigma_vv[v1, v2]))
        assert abs(s - sigma_vv[v1, v2]) < 1e-7

    opt = MDMin(UnitCellFilter(a), dt=0.01, logfile=None)
    opt.run(fmax=0.5)
    #print(a.cell)
    for i in range(3):
        for j in range(3):
            x = np.dot(a.cell[i], a.cell[j])
            y = (i + 1) * (j + 1) * a0**2 / 2
            if i != j:
                y /= 2
            #print(i, j, x, (x - y) / x)
            assert abs((x - y) / x) < 0.01

    print('passed test "stress"')

#=============================================================================
def main():
    test_H2_LDA()
    test_H2_PBE()
    test_H2_LS()
    test_O2()
    test_restart()
    test_GeoOpt()
    test_MD()
    test_stress()

#=============================================================================
if(__name__ == '__main__'):
    main()
elif(__name__ == '__builtin__'):
    if("CP2K" in os.environ.get('ASE_CALCULATORS', '').upper()):
        main()
    else:
        print('"CP2K" not in $ASE_CALCULATORS, skipping cp2k_test.py')

#EOF
