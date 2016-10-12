from __future__ import print_function

import numpy as np

from ase.build import bulk
from ase.constraints import UnitCellFilter
from ase.calculators.lj import LennardJones
from ase.optimize.precon import Exp, PreconLBFGS, PreconFIRE

N = 1
a0 = bulk('Cu', cubic=True)
a0 *= (N, N, N)

# perturb the atoms
s = a0.get_scaled_positions()
s[:, 0] *= 0.995
a0.set_scaled_positions(s)

nsteps = []
energies = []
for variable_cell in [False, True]:
    for OPT in [PreconLBFGS, PreconFIRE]:
        for precon in [None, Exp(A=3)]:
            atoms = a0.copy()
            atoms.set_calculator(LennardJones())
            if variable_cell:
                atoms = UnitCellFilter(atoms)
            opt = OPT(atoms, precon=precon, use_armijo=True)
            opt.run(1e-4)
            energies += [atoms.get_potential_energy()]
            nsteps += [opt.get_number_of_steps()]

# check we get the expected energy for all methods

# fixed cell cases
assert np.abs(
    np.array(energies[:len(energies) / 2]) - -0.92612723).max() < 1e-4

# variable cell cases
assert np.abs(
    np.array(energies[len(energies) / 2:]) - -31.7516156).max() < 1e-4
