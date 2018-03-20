"""Get energy from a LAMMPS calculation"""

from __future__ import print_function

import os
import numpy as np
from ase import Atom
from ase.build import bulk
from ase.calculators.lammpslib import LAMMPSlib
import ase.io
from ase import units
from ase.md.verlet import VelocityVerlet

# potential_path and data_file_path must be set as environment variables
potential_path = os.environ.get('LAMMPS_POTENTIALS_PATH', '.')
data_file_path = os.environ.get('LAMMPS_DATA_FILE_PATH', '.')

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


# a more complicated example, reading in a LAMMPS data file

Z_of_type = {1:26}
atom_types = {'Fe':1,}

at = ase.io.read(data_file_path+"/lammps.data", format="lammps-data", Z_of_type=Z_of_type)

header = ["units           real",
          "atom_style      full",
          "boundary        p p p",
          "box tilt        large",
          "pair_style      lj/cut/coul/long 12.500",
          "bond_style      harmonic",
          "angle_style     harmonic",
          "kspace_style    ewald 0.0001",
          "read_data       "+data_file_path+"/lammps.data"]
cmds = [] 

lammps = LAMMPSlib(lammps_header=header, lmpcmds=cmds, atom_types=atom_types, create_atoms=False, create_box=False, boundary=False, keep_alive=True, log_file='hi')
at.set_calculator(lammps)
dyn = VelocityVerlet(at, 1 * units.fs)

energy = at.get_potential_energy()
energy_ref = 2041.41198295
diff = abs((energy - energy_ref) / energy_ref)
assert diff < 1e-10

dyn.run(10)
energy = at.get_potential_energy()
energy_ref = 312.431585607
diff = abs((energy - energy_ref) / energy_ref)
assert diff < 1e-10, "%d" % energy
