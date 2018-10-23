import io
import numpy as np
import warnings

from ase.io import read
from ase.io import write

content = u"""
data_1


_chemical_name_common                  'Mysterious something'
_cell_length_a                         5.50000
_cell_length_b                         5.50000
_cell_length_c                         5.50000
_cell_angle_alpha                      90
_cell_angle_beta                       90
_cell_angle_gamma                      90
_space_group_name_H-M_alt              'F m -3 m'
_space_group_IT_number                 225

loop_
_space_group_symop_operation_xyz
   'x, y, z'
   '-x, -y, -z'
   '-x, -y, z'
   'x, y, -z'
   '-x, y, -z'
   'x, -y, z'
   'x, -y, -z'
   '-x, y, z'
   'z, x, y'
   '-z, -x, -y'
   'z, -x, -y'
   '-z, x, y'
   '-z, -x, y'
   'z, x, -y'
   '-z, x, -y'
   'z, -x, y'
   'y, z, x'
   '-y, -z, -x'
   '-y, z, -x'
   'y, -z, x'
   'y, -z, -x'
   '-y, z, x'
   '-y, -z, x'
   'y, z, -x'
   'y, x, -z'
   '-y, -x, z'
   '-y, -x, -z'
   'y, x, z'
   'y, -x, z'
   '-y, x, -z'
   '-y, x, z'
   'y, -x, -z'
   'x, z, -y'
   '-x, -z, y'
   '-x, z, y'
   'x, -z, -y'
   '-x, -z, -y'
   'x, z, y'
   'x, -z, y'
   '-x, z, -y'
   'z, y, -x'
   '-z, -y, x'
   'z, -y, x'
   '-z, y, -x'
   '-z, y, x'
   'z, -y, -x'
   '-z, -y, -x'
   'z, y, x'
   'x, y+1/2, z+1/2'
   '-x, -y+1/2, -z+1/2'
   '-x, -y+1/2, z+1/2'
   'x, y+1/2, -z+1/2'
   '-x, y+1/2, -z+1/2'
   'x, -y+1/2, z+1/2'
   'x, -y+1/2, -z+1/2'
   '-x, y+1/2, z+1/2'
   'z, x+1/2, y+1/2'
   '-z, -x+1/2, -y+1/2'
   'z, -x+1/2, -y+1/2'
   '-z, x+1/2, y+1/2'
   '-z, -x+1/2, y+1/2'
   'z, x+1/2, -y+1/2'
   '-z, x+1/2, -y+1/2'
   'z, -x+1/2, y+1/2'
   'y, z+1/2, x+1/2'
   '-y, -z+1/2, -x+1/2'
   '-y, z+1/2, -x+1/2'
   'y, -z+1/2, x+1/2'
   'y, -z+1/2, -x+1/2'
   '-y, z+1/2, x+1/2'
   '-y, -z+1/2, x+1/2'
   'y, z+1/2, -x+1/2'
   'y, x+1/2, -z+1/2'
   '-y, -x+1/2, z+1/2'
   '-y, -x+1/2, -z+1/2'
   'y, x+1/2, z+1/2'
   'y, -x+1/2, z+1/2'
   '-y, x+1/2, -z+1/2'
   '-y, x+1/2, z+1/2'
   'y, -x+1/2, -z+1/2'
   'x, z+1/2, -y+1/2'
   '-x, -z+1/2, y+1/2'
   '-x, z+1/2, y+1/2'
   'x, -z+1/2, -y+1/2'
   '-x, -z+1/2, -y+1/2'
   'x, z+1/2, y+1/2'
   'x, -z+1/2, y+1/2'
   '-x, z+1/2, -y+1/2'
   'z, y+1/2, -x+1/2'
   '-z, -y+1/2, x+1/2'
   'z, -y+1/2, x+1/2'
   '-z, y+1/2, -x+1/2'
   '-z, y+1/2, x+1/2'
   'z, -y+1/2, -x+1/2'
   '-z, -y+1/2, -x+1/2'
   'z, y+1/2, x+1/2'
   'x+1/2, y, z+1/2'
   '-x+1/2, -y, -z+1/2'
   '-x+1/2, -y, z+1/2'
   'x+1/2, y, -z+1/2'
   '-x+1/2, y, -z+1/2'
   'x+1/2, -y, z+1/2'
   'x+1/2, -y, -z+1/2'
   '-x+1/2, y, z+1/2'
   'z+1/2, x, y+1/2'
   '-z+1/2, -x, -y+1/2'
   'z+1/2, -x, -y+1/2'
   '-z+1/2, x, y+1/2'
   '-z+1/2, -x, y+1/2'
   'z+1/2, x, -y+1/2'
   '-z+1/2, x, -y+1/2'
   'z+1/2, -x, y+1/2'
   'y+1/2, z, x+1/2'
   '-y+1/2, -z, -x+1/2'
   '-y+1/2, z, -x+1/2'
   'y+1/2, -z, x+1/2'
   'y+1/2, -z, -x+1/2'
   '-y+1/2, z, x+1/2'
   '-y+1/2, -z, x+1/2'
   'y+1/2, z, -x+1/2'
   'y+1/2, x, -z+1/2'
   '-y+1/2, -x, z+1/2'
   '-y+1/2, -x, -z+1/2'
   'y+1/2, x, z+1/2'
   'y+1/2, -x, z+1/2'
   '-y+1/2, x, -z+1/2'
   '-y+1/2, x, z+1/2'
   'y+1/2, -x, -z+1/2'
   'x+1/2, z, -y+1/2'
   '-x+1/2, -z, y+1/2'
   '-x+1/2, z, y+1/2'
   'x+1/2, -z, -y+1/2'
   '-x+1/2, -z, -y+1/2'
   'x+1/2, z, y+1/2'
   'x+1/2, -z, y+1/2'
   '-x+1/2, z, -y+1/2'
   'z+1/2, y, -x+1/2'
   '-z+1/2, -y, x+1/2'
   'z+1/2, -y, x+1/2'
   '-z+1/2, y, -x+1/2'
   '-z+1/2, y, x+1/2'
   'z+1/2, -y, -x+1/2'
   '-z+1/2, -y, -x+1/2'
   'z+1/2, y, x+1/2'
   'x+1/2, y+1/2, z'
   '-x+1/2, -y+1/2, -z'
   '-x+1/2, -y+1/2, z'
   'x+1/2, y+1/2, -z'
   '-x+1/2, y+1/2, -z'
   'x+1/2, -y+1/2, z'
   'x+1/2, -y+1/2, -z'
   '-x+1/2, y+1/2, z'
   'z+1/2, x+1/2, y'
   '-z+1/2, -x+1/2, -y'
   'z+1/2, -x+1/2, -y'
   '-z+1/2, x+1/2, y'
   '-z+1/2, -x+1/2, y'
   'z+1/2, x+1/2, -y'
   '-z+1/2, x+1/2, -y'
   'z+1/2, -x+1/2, y'
   'y+1/2, z+1/2, x'
   '-y+1/2, -z+1/2, -x'
   '-y+1/2, z+1/2, -x'
   'y+1/2, -z+1/2, x'
   'y+1/2, -z+1/2, -x'
   '-y+1/2, z+1/2, x'
   '-y+1/2, -z+1/2, x'
   'y+1/2, z+1/2, -x'
   'y+1/2, x+1/2, -z'
   '-y+1/2, -x+1/2, z'
   '-y+1/2, -x+1/2, -z'
   'y+1/2, x+1/2, z'
   'y+1/2, -x+1/2, z'
   '-y+1/2, x+1/2, -z'
   '-y+1/2, x+1/2, z'
   'y+1/2, -x+1/2, -z'
   'x+1/2, z+1/2, -y'
   '-x+1/2, -z+1/2, y'
   '-x+1/2, z+1/2, y'
   'x+1/2, -z+1/2, -y'
   '-x+1/2, -z+1/2, -y'
   'x+1/2, z+1/2, y'
   'x+1/2, -z+1/2, y'
   '-x+1/2, z+1/2, -y'
   'z+1/2, y+1/2, -x'
   '-z+1/2, -y+1/2, x'
   'z+1/2, -y+1/2, x'
   '-z+1/2, y+1/2, -x'
   '-z+1/2, y+1/2, x'
   'z+1/2, -y+1/2, -x'
   '-z+1/2, -y+1/2, -x'
   'z+1/2, y+1/2, x'

loop_
   _atom_site_label
   _atom_site_occupancy
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_adp_type
   _atom_site_B_iso_or_equiv
   _atom_site_type_symbol
   Na         0.7500  0.000000      0.000000      0.000000     Biso  1.000000 Na
   K          0.2500  0.000000      0.000000      0.000000     Biso  1.000000 K
   Cl         0.3000  0.500000      0.500000      0.500000     Biso  1.000000 Cl
   I          0.5000  0.250000      0.250000      0.250000     Biso  1.000000 I
"""

cif_file = io.StringIO(content)

# legacy behavior is to not read the K atoms
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    atoms = read(cif_file, format='cif', fractional_occupancies=False)
elements = np.unique(atoms.get_atomic_numbers())
for n in (11, 17, 53):
    assert n in elements
try:
    atoms.info['occupancy']
    raise AssertionError
except KeyError:
    pass

cif_file = io.StringIO(content)
# new behavior is to still not read the K atoms, but build tags and info
natoms = read(cif_file, format='cif', fractional_occupancies=True)

assert len(atoms) == len(natoms)
assert np.all(atoms.get_atomic_numbers() == natoms.get_atomic_numbers())
# yield the same old atoms...
assert atoms == natoms

elements = np.unique(atoms.get_atomic_numbers())
for n in (11, 17, 53):
    assert n in elements

assert natoms.info['occupancy']
for a in natoms:
    if a.symbol == 'Na':
        assert len(natoms.info['occupancy'][a.tag]) == 2
        assert natoms.info['occupancy'][a.tag]['K'] == 0.25
        assert natoms.info['occupancy'][a.tag]['Na'] == 0.75
    else:
        assert len(natoms.info['occupancy'][a.tag]) == 1

# read/write
fname = 'testfile.cif'
with open(fname, 'w') as fd:
    write(fd, natoms, format='cif')

with open(fname) as fd:
    natoms = read(fd, format='cif', fractional_occupancies=True)

assert natoms.info['occupancy']
for a in natoms:
    if a.symbol == 'Na':
        assert len(natoms.info['occupancy'][a.tag]) == 2
        assert natoms.info['occupancy'][a.tag]['K'] == 0.25
        assert natoms.info['occupancy'][a.tag]['Na'] == 0.75
    else:
        assert len(natoms.info['occupancy'][a.tag]) == 1
