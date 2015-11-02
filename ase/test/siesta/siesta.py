from __future__ import print_function

import os
import numpy as np

from ase.units import Ry
from ase.calculators.siesta.siesta import Siesta
from ase.calculators.siesta.siesta import Specie, BasisSet, Shell
from ase.calculators.calculator import FileIOCalculator
from ase import Atoms

try:
    del os.environ['SIESTA']
except KeyError:
    pass

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
os.system('touch %s/H.lda.psf' % pseudo_path)
os.system('touch %s/C.lda.psf' % pseudo_path)
os.system('touch %s/O.lda.psf' % pseudo_path)

h = Atoms('H', [(0.0, 0.0, 0.0)])
co2 = Atoms('CO2', [(0.0, 0.0, 0.0), (-1.178, 0.0, 0.0), (1.178, 0.0, 0.0)])
ch4 = Atoms('CH4', np.array([
    [0.000000, 0.000000, 0.000000],
    [0.682793, 0.682793, 0.682793],
    [-0.682793, -0.682793, 0.682790],
    [-0.682793, 0.682793, -0.682793],
    [0.682793, -0.682793, -0.682793]]))

dirt_cheap_siesta = Siesta(
    label='test_label',
    xc='LDA',
    mesh_cutoff=70*Ry,
    basis_set='SZ',
    DM_Tolerance=1e-3,
)

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

# Test complex basis set.
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
        l=1,
        nzeta=3,
        split_norm=0.2,
        nzetapol=1,
        soft_confinement='Junquera',
        soft_confinement_prefactor=0.2,
        soft_confinement_inner_radius=6.0,
        scale_factors=[1.0, 0.95, 0.99],
        )
basis_set=BasisSet(
    shells=[shell1, shell2],
    basis_type='split',
    ionic_charge=1.0,
    )

specie = Specie(symbol='C', basis_set=basis_set)
siesta = Siesta(
    label='test_label',
    xc='LDA',
    mesh_cutoff=70*Ry,
    species=[specie],
    DM_Tolerance=1e-3,
)
atoms.set_calculator(siesta)
siesta.write_input(atoms, properties=['energy'])


atoms = h.copy()
atoms.set_calculator(dirt_cheap_siesta)
dirt_cheap_siesta.write_input(atoms, properties=['energy'])
with open('test_label.fdf', 'r') as f:
    lines = f.readlines()
assert 'DM.Tolerance  0.001\n' in lines

# Remove the test directory.
os.chdir('../..')
os.system('rm -rf %s' % test_path)
