from __future__ import print_function

import os
import numpy as np

from ase.units import Ry
from ase.calculators.siesta.siesta import Siesta
from ase.calculators.siesta.parameters import Specie, PAOBasisBlock
from ase.calculators.calculator import FileIOCalculator
from ase import Atoms

test_path = 'tmp_siesta'
if not os.path.exists(test_path):
    os.makedirs(test_path)
os.chdir(test_path)

run_path = 'run_directory'
pseudo_path = 'pseudos'
if not os.path.exists(pseudo_path):
    os.makedirs(pseudo_path)
if not os.path.exists(run_path):
    os.makedirs(run_path)
for symbol in 'HCO':
    with open('{0}/{1}.lda.psf'.format(pseudo_path, symbol), 'w') as fd:
        fd.close()

h = Atoms('H', [(0.0, 0.0, 0.0)])
co2 = Atoms('CO2', [(0.0, 0.0, 0.0), (-1.178, 0.0, 0.0), (1.178, 0.0, 0.0)])
ch4 = Atoms('CH4', np.array([
    [0.000000, 0.000000, 0.000000],
    [0.682793, 0.682793, 0.682793],
    [-0.682793, -0.682793, 0.682790],
    [-0.682793, 0.682793, -0.682793],
    [0.682793, -0.682793, -0.682793]]))

os.chdir(run_path)
os.environ['SIESTA_PP_PATH'] = '../' + pseudo_path

# Test the initialization
siesta = Siesta()
assert isinstance(siesta, FileIOCalculator)
assert isinstance(siesta.implemented_properties, tuple)
assert isinstance(siesta.default_parameters, dict)
assert isinstance(siesta.name, str)
assert isinstance(siesta.default_parameters, dict)

# Test the more complicated species situations
atoms = ch4.copy()
species, numbers = siesta.species(atoms)
assert all(numbers == np.array([1, 2, 2, 2, 2]))
siesta = Siesta(species=[Specie(symbol='C', tag=1)])
species, numbers = siesta.species(atoms)
assert all(numbers == np.array([1, 2, 2, 2, 2]))
atoms.set_tags([0, 0, 0, 1, 0])
species, numbers = siesta.species(atoms)
assert all(numbers == np.array([1, 2, 2, 2, 2]))
siesta = Siesta(species=[Specie(symbol='H', tag=1)])
species, numbers = siesta.species(atoms)
assert all(numbers == np.array([1, 2, 2, 3, 2]))


c_basis = """2 nodes 1.00
0 1 S 0.20 P 1 0.20 6.00
5.00
1.00
1 2 S 0.20 P 1 E 0.20 6.00
6.00 5.00
1.00 0.95"""
basis_set = PAOBasisBlock(c_basis)
specie = Specie(symbol='C', basis_set=basis_set)

# Test that the fdf_arguments come first.
siesta = Siesta(
    label='test_label',
    fdf_arguments={'DM.Tolerance': 1e-3})
atoms.set_calculator(siesta)
siesta.write_input(atoms, properties=['energy'])
atoms = h.copy()
atoms.set_calculator(siesta)
siesta.write_input(atoms, properties=['energy'])
with open('test_label.fdf', 'r') as f:
    lines = f.readlines()
assert 'DM.Tolerance  0.001\n' == lines[0]

siesta = Siesta(
    label='test_label',
    fdf_arguments={
        'DM.Tolerance': 1e-3,
        'ON.eta': 5 * Ry})

atoms.set_calculator(siesta)
siesta.write_input(atoms, properties=['energy'])

atoms = h.copy()
atoms.set_calculator(siesta)
siesta.write_input(atoms, properties=['energy'])
with open('test_label.fdf', 'r') as f:
    lines = f.readlines()
assert 'DM.Tolerance  0.001\n' in lines
assert 'ON.eta  68.02848914 \teV\n' in lines

# Remove the test directory.
os.chdir('../..')
os.system('rm -rf %s' % test_path)
