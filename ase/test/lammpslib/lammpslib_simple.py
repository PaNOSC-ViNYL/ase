"""Get energy from a LAMMPS calculation"""

from __future__ import print_function

import os
import numpy as np
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
                   log_file='test.log', keep_alive=True)

nickel.set_calculator(lammps)

E = nickel.get_potential_energy()
F = nickel.get_forces()
S = nickel.get_stress()

print('Energy: ', E)
print('Forces:', F)
print('Stress: ', S)
print()

E = nickel.get_potential_energy()
F = nickel.get_forces()
S = nickel.get_stress()


lammps = LAMMPSlib(lmpcmds=cmds,
                   log_file='test.log', keep_alive=True)
nickel.set_calculator(lammps)

E2 = nickel.get_potential_energy()
F2 = nickel.get_forces()
S2 = nickel.get_stress()

assert np.allclose(E, E2)
assert np.allclose(F, F2)
assert np.allclose(S, S2)

nickel.rattle(stdev=0.2)
E3 = nickel.get_potential_energy()
F3 = nickel.get_forces()
S3 = nickel.get_stress()

print('rattled atoms')
print('Energy: ', E3)
print('Forces:', F3)
print('Stress: ', S3)
print()

assert not np.allclose(E, E3)
assert not np.allclose(F, F3)
assert not np.allclose(S, S3)

nickel += Atom('H', position=nickel.cell.diagonal()/4)
E4 = nickel.get_potential_energy()
F4 = nickel.get_forces()
S4 = nickel.get_stress()

assert not np.allclose(E4, E3)
assert not np.allclose(F4[:-1,:], F3)
assert not np.allclose(S4, S3)


# the example from the docstring

cmds = ["pair_style eam/alloy",
        "pair_coeff * * {path}/NiAlH_jea.eam.alloy Al H".format(path=potential_path)]

Ni = bulk('Ni', cubic=True)
H = Atom('H', position=Ni.cell.diagonal()/2)
NiH = Ni + H

lammps = LAMMPSlib(lmpcmds=cmds, log_file='test.log')

NiH.set_calculator(lammps)
print("Energy ", NiH.get_potential_energy())
