import argparse
import traceback
from time import time

import numpy as np

import ase.db
import ase.optimize
from ase.calculators.calculator import get_calculator
from ase.io import Trajectory


all_optimizers = ase.optimize.__all__ + ['PreconLBFGS', 'PreconFIRE',
                                         'SciPyFminCG', 'SciPyFminBFGS']


def get_optimizer(name):
    if name.startswith('Precon'):
        import ase.optimize.precon as precon
        return getattr(precon, name)
    if name.startswith('SciPy'):
        import ase.optimize.sciopt as sciopt
        return getattr(sciopt, name)
    return getattr(ase.optimize, name)


class Wrapper:
    def __init__(self, atoms):
        self.t0 = time()
        self.t = 0.0
        self.e = []
        self.f = []
        self.atoms = atoms
        self.energy_ready = False
        self.forces_ready = False

    @property
    def cell(self):
        return self.atoms.cell

    def get_cell(self, complete=False):
        return self.atoms.get_cell(complete)

    @property
    def pbc(self):
        return self.atoms.pbc

    @property
    def positions(self):
        return self.atoms.positions

    @property
    def constraints(self):
        return self.atoms.constraints

    def copy(self):
        return self.atoms.copy()

    def get_calculator(self):
        return self.atoms.calc

    def __len__(self):
        return len(self.atoms)

    def get_potential_energy(self, force_consistent=False):
        t1 = time()
        e = self.atoms.get_potential_energy(force_consistent)
        t2 = time()
        self.t += t2 - t1
        if not self.energy_ready:
            self.e.append((t2 - self.t0, e))
        self.energy_ready = True
        return e

    def get_forces(self):
        t1 = time()
        f = self.atoms.get_forces()
        t2 = time()
        self.t += t2 - t1
        if not self.forces_ready:
            self.f.append((t2 - self.t0, (f**2).sum(1).max()))
        self.forces_ready = True
        return f

    def set_positions(self, pos):
        self.atoms.set_positions(pos)
        self.energy_ready = False
        self.forces_ready = False

    def get_positions(self):
        return self.atoms.get_positions()


def run_test(atoms, optimizer, tag, fmax=0.02):
    wrapper = Wrapper(atoms)
    relax = optimizer(wrapper, logfile=tag + '.log')
    relax.attach(Trajectory(tag + '.traj', 'w', atoms=atoms))

    tincl = -time()
    error = ''

    try:
        relax.run(fmax=fmax, steps=100)
    except Exception as x:
        error = '{}: {}'.format(x.__class__.__name__, x)
        tb = traceback.format_exc()

        with open(tag + '.err', 'w') as fd:
            fd.write('{}\n{}\n'.format(error, tb))

    tincl += time()
    texcl = wrapper.t

    return error, wrapper.e, wrapper.f, texcl, tincl


def test_optimizer(systems, optimizer, db=None):
    for row in systems.select():
        if db is not None:
            optname = optimizer.__name__
            id = db.reserve(optimizer=optname, sid=row.id)
            if id is None:
                continue
        atoms = row.toatoms()
        tag = '{}-{}-{}'.format(row.calculator, optname, row.formula)
        params = row.calculator_parameters
        atoms.calc = get_calculator(row.calculator)(**params, txt=tag + '.txt')
        error, e, f, texcl, tincl = run_test(atoms, optimizer, tag)

        if db is not None:
            db.write(atoms,
                     id=id,
                     optimizer=optname,
                     error=error,
                     ne=len(e),
                     nf=len(f),
                     t=texcl,
                     T=tincl,
                     sid=row.id,
                     data={'e': np.array(e).T,
                           'f': np.array(f).T})


def main():
    parser = argparse.ArgumentParser(
        description='Test ASE optimizers')

    parser.add_argument('systems')
    parser.add_argument('optimizer', nargs='*',
                        help='Optimizer name.')

    args = parser.parse_args()

    systems = ase.db.connect(args.systems)

    db = ase.db.connect('results.db')

    if not args.optimizer:
        args.optimizer = all_optimizers

    for opt in args.optimizer:
        print(opt)
        optimizer = get_optimizer(opt)
        test_optimizer(systems, optimizer, db)


if __name__ == '__main__':
    main()
