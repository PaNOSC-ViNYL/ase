"""
Run tests to ensure that the VASP check_state() function call works correctly,
i.e. correctly sets the working directories and works in that directory.

This is conditional on the existence of the VASP_COMMAND or VASP_SCRIPT
environment variables

"""

from ase.test.vasp import installed

from ase import Atoms
from ase.calculators.vasp import Vasp2 as Vasp
assert installed()

# Test setup system, borrowed from vasp_co.py
d = 1.14
atoms = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
              pbc=True)
atoms.center(vacuum=5.)

# Test
settings = dict(xc='PBE',
                prec='Low',
                algo='Fast',
                ismear=0,
                sigma=1.,
                istart=0,
                lwave=False,
                lcharg=False)

calc = Vasp(**settings)

atoms.set_calculator(calc)

en1 = atoms.get_potential_energy()

# Make a parameter change to the calculator
calc.set(sigma=0.5)

# Check that we capture a change for float params
assert calc.check_state(atoms) == ['float_params']
assert calc.calculation_required(atoms, ['energy', 'forces'])

en2 = atoms.get_potential_energy()

# The change in sigma should result in a small change in energy
assert (en1 - en2) > 1e-7

# Now we make a change in input_params instead
calc.kpts = 2

# Check that this requires a new calculation
assert calc.check_state(atoms) == ['input_params']
assert calc.calculation_required(atoms, ['energy', 'forces'])

# Clean up
calc.clean()
