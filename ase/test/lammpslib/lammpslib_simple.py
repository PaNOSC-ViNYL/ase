"""Get energy from a LAMMPS calculation"""

from __future__ import print_function

import os

from ase import Atom
from ase.build import bulk
from ase.calculators.lammpslib import LAMMPSlib

potential_path = os.environ.get('LAMMPS_POTENTIALS_PATH', '.')

cmds = ["pair_style eam/alloy",
        "pair_coeff * * {path}/NiAlH_jea.eam.alloy Ni H"
        "".format(path=potential_path)]

nickel = bulk('Ni', cubic=True)
nickel += Atom('H', position=nickel.cell.diagonal()/2)
# Bit of distortion
nickel.set_cell(nickel.cell + [[0.1, 0.2, 0.4],
                               [0.3, 0.2, 0.0],
                               [0.1, 0.1, 0.1]], scale_atoms=True)

lammps = LAMMPSlib(lmpcmds=cmds,
                   atom_types={'Ni': 1, 'H': 2},
                   logfile='test.log', keep_alive=True)

nickel.set_calculator(lammps)

print('Energy: ', nickel.get_potential_energy())
print('Forces:', nickel.get_forces())
print('Stress: ', nickel.get_stress())
