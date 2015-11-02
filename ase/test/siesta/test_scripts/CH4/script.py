import os
from ase.io import read

from ase.units import Ry
from ase.calculators.siesta.siesta import Siesta, Shell, BasisSet
from ase.calculators.siesta.siesta import Specie
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
shell1 = Shell(
        l=0,
        nzeta=1,
        split_norm=0.2,
        nzetapol=1,
        soft_confinement='Sankey',
        soft_confinement_prefactor=0.2,
        soft_confinement_inner_radius=6.0,
        scale_factors=[0.99],
        )
shell2 = Shell(
        l=0,
        nzeta=2,
        split_norm=0.2,
        nzetapol=1,
        soft_confinement='Junquera',
        soft_confinement_prefactor=0.2,
        soft_confinement_inner_radius=6.0,
        scale_factors=[1.0, 0.95],
        )
basis_set=BasisSet(
    shells=[shell1, shell2],
    basis_type='nodes',
    ionic_charge=1.0,
    )

specie = Specie(symbol='C', basis_set=basis_set)
calc = Siesta(
        label='ch4',
        basis_set='SZ',
        xc='CA',
        mesh_cutoff=200*Ry,
        species=[specie],
        DM_Tolerance=1e-4,
        DM_MixingWeight=0.15,
        DM_NumberPulay=3,
        ElectronicTemperature='300 K',
    )

bud.set_calculator(calc)
dyn = QuasiNewton(bud, trajectory='bud.traj')
dyn.run(fmax=0.02)
e = bud.get_potential_energy()
