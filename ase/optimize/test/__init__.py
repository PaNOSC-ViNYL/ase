"""Define a helper function for running tests

The skeleton for making a new setup is as follows:

from ase.optimize.test import run_test

def get_atoms():
    return Atoms('H')

def get_calculator():
    return EMT()

run_test(get_atoms, get_calculator, 'Hydrogen')
"""

from ase.optimize.bfgs import BFGS
from ase.optimize.lbfgs import LBFGS

optimizers = [
    'BFGS',
    'LBFGS',
]

def get_optimizer(optimizer):
    if optimizer == 'BFGS':
        return BFGS
    elif optimizer == 'LBFGS':
        return LBFGS

def run_test(get_atoms, get_calculator, name,
             fmax=0.05, steps = 100):
    for optimizer in optimizers:
        logname = name + '-' + optimizer

        atoms = get_atoms()
        atoms.set_calculator(get_calculator())
        opt = get_optimizer(optimizer)
        relax = opt(atoms,
                    logfile = logname + '.log',
                    trajectory = logname + '.traj')
        relax.run(fmax = fmax, steps = steps)
        nsteps = relax.get_number_of_steps()
        E = atoms.get_potential_energy()
        print '%-15s %-5s %3i %8.3f' % (name, optimizer, nsteps, E)
