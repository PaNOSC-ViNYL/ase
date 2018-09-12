from __future__ import print_function
from ase.units import Ry
from ase.calculators.openmx import OpenMX
from ase.io.trajectory import Trajectory
from ase.optimize import QuasiNewton
from ase.constraints import UnitCellFilter
from ase import Atoms
import numpy as np
""" Only OpenMX 3.8 or higher version pass this test"""


bud = Atoms('CH4', np.array([
        [0.000000, 0.000000, 0.100000],
        [0.682793, 0.682793, 0.682793],
        [-0.682793, -0.682793, 0.68279],
        [-0.682793, 0.682793, -0.682793],
        [0.682793, -0.682793, -0.682793]]),
        cell=[10, 10, 10])

calc = OpenMX(
    label='ch4',
    xc='GGA',
    energy_cutoff=300 * Ry,
    stress='on',
    pbs={'processes': 20}
    )

bud.set_calculator(calc)
traj = Trajectory('example.traj', 'w', bud)
ucf = UnitCellFilter(bud, mask=[True, True, False, False, False, False])
dyn = QuasiNewton(ucf)
dyn.attach(traj.write)
dyn.run(fmax=0.02)
e = bud.get_potential_energy()

traj.close()
