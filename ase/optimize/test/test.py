import argparse
import io
import sys
import traceback
from time import time

import numpy as np

import ase.optimize
from ase import Atoms
from ase.calculators.calculator import get_calculator
from ase.io import read, Trajectory
from ase.neb import NEB


all_optimizers = ase.optimize.__all__


def get_atoms_and_name(atoms):
    if isinstance(atoms, str):
        name = atoms
        if name.endswith('.py'):
            dct = {}
            exec(compile(open(name).read(), name, 'exec'), dct)
            for atoms in dct.values():
                if isinstance(atoms, (Atoms, NEB)):
                    break
        else:
            atoms = read(name)
    else:
        name = atoms.get_name()

    return atoms, name


class Wrapper:
    def __init__(self, atoms):
        self.t0 = time()
        self.t = 0.0
        self.e = []
        self.f = []
        self.atoms = atoms
        self.energy_ready = False
        self.forces_ready = False

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
    parser.add_argument('optimizer', nargs='+',
                        help='Optimizer name.')

    args = parser.parse_args()

    systems = ase.db.connect(args.systems)

    db = ase.db.connect('results.db')

    for opt in args.optimizer:
        optimizer = getattr(ase.optimize, opt)
        test_optimizer(systems, optimizer, db)


if __name__ == '__main__':
    main()
