from __future__ import print_function
from ase.calculators.turbomole import Turbomole
from ase.build import molecule

mol = molecule('H2O')

params = {
    'title': 'water',
    'task': 'geometry optimization',
    'use redundant internals': True,
    'basis set name': 'def2-SV(P)',
    'total charge': 0,
    'multiplicity': 1,
    'use dft': True,
    'density functional': 'b3-lyp',
    'use resolution of identity': True,
    'ri memory': 1000,
    'force convergence': 0.001,
    'geometry iterations': 50,
    'maximum number of scf iterations': 100
}

calc = Turbomole(**params)
mol.set_calculator(calc)
calc.calculate(mol)

# use the get_property() method
print(calc.get_property('energy', mol, False))
print(calc.get_property('forces', mol, False))
print(calc.get_property('dipole', mol, False))

# use the get_results() method
results = calc.get_results()
print(results['molecular orbitals'])

# use the __getitem__() method
print(calc['results']['molecular orbitals'])
print(calc['results']['geometry optimization history'])

