import os
from ase.io import read
from ase.calculators.siesta.siesta import Siesta
from ase.optimize import QuasiNewton
from ase.visualize import view
from ase import Atoms
import numpy as np

os.environ['SIESTA_PP_PATH'] = os.path.abspath('../../TestFiles')

bud = Atoms('CH4', np.array([
          [0.000000,  0.000000,  0.100000],
          [0.682793,  0.682793,  0.682793],
          [-0.682793, -0.682793,  0.68279],
          [-0.682793,  0.682793, -0.682793],
          [0.682793, -0.682793, -0.682793]]),
          cell=[10, 10, 10],
          )
# Uncomment to use the last image of the relaxation trajectory.
#bud = read('bud.traj')

calc = Siesta(
        label='ch4',
        basis_set='SZ',
        xc='LDA.CA',
        mesh_cutoff=(200, 'Ry'),
        restart='ch4.XV',
        ignore_bad_restart_file=False,
        DM_Tolerance=1e-4,
        DM_MixingWeight=0.15,
        DM_NumberPulay=3,
        ElectronicTemperature=(300, 'K'),
              )
bud.set_calculator(calc)
dyn = QuasiNewton(bud, trajectory='bud.traj')
dyn.run(fmax=0.02)
e = bud.get_potential_energy()
