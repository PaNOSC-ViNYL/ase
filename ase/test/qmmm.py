from math import cos, sin

import numpy as np
import matplotlib.pyplot as plt

from ase import Atoms
from ase.calculators.tip3p import (TIP3P, epsilon0, sigma0, rOH, thetaHOH,
                                   set_tip3p_charges)
from ase.calculators.qmmm import QMMM1, QMMM2, LJInteractions
from ase.constraints import FixInternals
from ase.optimize import BFGS

r = rOH
a = thetaHOH
D = np.linspace(1.5, 2.5, 30)


def test(dimer):
    E = []
    F = []
    for d in D:
        dimer.positions[3:, 0] += d - (dimer.positions[5, 0] - r)
        E.append(dimer.get_potential_energy())
        F.append(dimer.get_forces())
    return np.array(E), np.array(F)
    
i = LJInteractions({('O', 'O'): (epsilon0, sigma0)})
for calc in [TIP3P(),
             QMMM1([0, 1, 2], TIP3P(), TIP3P(), TIP3P()),
             QMMM2([0, 1, 2], TIP3P(), TIP3P(), i)]:
    dimer = Atoms('H2OH2O',
                  [(r * cos(a), r * sin(a), 0),
                   (r, 0, 0),
                   (0, 0, 0),
                   (r * cos(a / 2), r * sin(a / 2), 0),
                   (r * cos(a / 2), -r * sin(a / 2), 0),
                   (0, 0, 0)])
    set_tip3p_charges(dimer)
    dimer.calc = calc

    E, F = test(dimer)
    F1 = np.polyval(np.polyder(np.polyfit(D, E, 7)), D)
    F2 = F[:, :3, 0].sum(1)
    error = abs(F1 - F2).max()

    dimer.constraints = FixInternals(
        bonds=[(r, (0, 2)), (r, (1, 2)),
               (r, (3, 5)), (r, (4, 5))],
        angles=[(a, (0, 2, 1)), (a, (3, 5, 4))])
    opt = BFGS(dimer,
               trajectory=calc.name + '.traj', logfile=calc.name + 'd.log')
    opt.run(0.01)

    e0 = dimer.get_potential_energy()
    d0 = dimer.get_distance(1, 5)
    print(calc.name, e0, d0, error)
    plt.plot(D, E)
    
plt.show()
