from __future__ import print_function
import sys
import numpy as np

from ase.atoms import Atoms
from ase.calculators.lj import LennardJones
from ase.calculators.symmetrize import SymmetrizedCalculator, check
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

# setup an fcc Al cell
a = 2.0/np.sqrt(3.0)
at_prim = Atoms('Al2', positions=[[0,0,0],[a/2.0, a/2.0, a/2.0]],
                cell=[[a,0,0],[0,a,0],[0,0,a]], pbc=[True, True, True])

# without symmetrization
at_unsym = at_prim * [2,2,2]
at_unsym.positions[0,0] += 1.0e-7 # break symmetry by 1e-7
at_init = at_unsym.copy()

calc = LennardJones()
at_unsym.set_calculator(calc)

at_cell = UnitCellFilter(at_unsym)

dyn = BFGS(at_cell)
print("Energy", at_unsym.get_potential_energy())
dyn.run(fmax=0.001)
print("Energy", at_unsym.get_potential_energy())

# with symmetrization
at_sym = at_prim * [2,2,2]
at_sym.positions[0,0] += 1.0e-7 # break symmetry by 1e-7

calc = SymmetrizedCalculator(LennardJones(), at_sym)
at_sym.set_calculator(calc)

at_cell = UnitCellFilter(at_sym)

dyn = BFGS(at_cell)
print("Energy", at_sym.get_potential_energy())
dyn.run(fmax=0.001)
print("Energy", at_sym.get_potential_energy())

print("position difference", np.linalg.norm(at_unsym.get_positions()-at_sym.get_positions()))

print("initial symmetry at 1e-6")
d_init6 = check(at_init, 1.0e-6)
print("initial symmetry at 1e-8")
d_init8 = check(at_init, 1.0e-8)
print("unsym symmetry")
d_unsym = check(at_unsym)
print("sym symmetry")
d_sym = check(at_sym)
