import os
from ase.units import Ry, eV, Ang

from ase.calculators.siesta import Siesta
from ase.calculators.siesta.parameters import Specie
from ase.optimize import QuasiNewton
from ase import Atoms
from ase.io import read
import sys

Na8 = read('Na8.xyz')
Na8.set_cell([20.0, 20.0, 20.0])

siesta = Siesta(
    mesh_cutoff=150 * Ry,
    basis_set='DZP',
    pseudo_path=os.getcwd(),
    pseudo_qualifier = '',
    energy_shift=(10*10**-3) * eV,
    fdf_arguments={
                   'SCFMustConverge':False,
                   'COOP.Write': True,
                   'WriteDenchar': True,
                   'PAO.BasisType': 'split',
                   'DM.Tolerance': 1e-4,
                   'DM.MixingWeight': 0.01,
                   'MaxSCFIterations': 3,
                   'DM.NumberPulay': 4}
)

Na8.set_calculator(siesta)
Na8.get_potential_energy()
