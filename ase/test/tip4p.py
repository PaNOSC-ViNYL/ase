"""Test TIP4P forces."""
from math import cos, sin

from ase import Atoms
from ase.calculators.tip4p import TIP4P, rOH, thetaHOH

r = rOH
a = thetaHOH

dimer = Atoms('H2OH2O',
              [(r * cos(a), 0, r * sin(a)),
               (r, 0, 0),
               (0, 0, 0),
               (r * cos(a / 2), r * sin(a / 2), 0),
               (r * cos(a / 2), -r * sin(a / 2), 0),
               (0, 0, 0)])
dimer.calc = TIP4P(rc=4.0, width=2.0)  # put O-O distance in the cutoff range
dimer.positions[3:, 0] += 2.8
F = dimer.get_forces()
print(F)
dF = dimer.calc.calculate_numerical_forces(dimer) - F
print(dF)
assert abs(dF).max() < 2e-6
