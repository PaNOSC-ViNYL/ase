import os
from os.path import join
import numpy as np

from ase.calculators.siesta.basis_set import SZ
from ase.calculators.siesta import Siesta
from ase.optimize import QuasiNewton
from ase import Atoms

os.environ['SIESTA_PP_PATH'] = os.path.abspath(os.path.dirname(__file__))

h = Atoms('2H', [(0.0, 0.0, 0.0), (0.0, 0.0, 1.0)], cell=[10,10,10])

dirt_cheap_siesta = Siesta(
        mesh_cutoff=(200, 'Ry'),
        basis_set=SZ(),
        fdf_arguments={'DM_Tolerance': 1e-4}
        )

h.set_calculator(dirt_cheap_siesta)
dyn = QuasiNewton(h, trajectory='h.traj')
dyn.run(fmax=0.01)

print h.positions

