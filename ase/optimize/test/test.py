import argparse
import traceback
from time import time

import matplotlib.pyplot as plt
import numpy as np

import ase.optimize
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read
from ase.neb import NEB


all_optimizers = ['MDMin', 'FIRE', 'LBFGS',
                  'LBFGSLineSearch', 'BFGSLineSearch', 'BFGS']


def get_optimizer(name):
    return getattr(ase.optimize, name)


def get_atoms_and_name(atoms):
    if isinstance(atoms, str):
        name = atoms
        if name.endswith('.py'):
            dct = {}
            exec(compile(open(name).read(), name, 'exec'), dct)
            atoms_objects = []
            for thing in dct.values():
                if isinstance(thing, NEB):
                    atoms = thing
                    break
                if isinstance(thing, Atoms):
                    atoms_objects.append(thing)
            else:
                assert len(atoms_objects) == 1, name
                atoms = atoms_objects[0]
        else:
            atoms = read(name)
    else:
        name = atoms.get_name()

    return atoms, name


class Wrapper:
    def __init__(self, atoms):
        self.t = 0.0
        self.e = []
        self.f = []
        self.atoms = atoms
        self.old_positions = np.zeros_like(atoms.get_positions())
        self.n = 0
    def __len__(self):
        return len(self.atoms)

    def get_potential_energy(self):
        self.t = -time()
        if (self.atoms.get_positions() != self.old_positions).any():
            self.n += 1
        self.old_positions = self.atoms.get_positions()
        e = self.atoms.get_potential_energy()
        self.t += time()
        self.e.append((self.t, e))
        return e

    def get_forces(self):
        self.t = -time()
        if (self.atoms.get_positions() != self.old_positions).any():
            self.n += 1
        self.old_positions = self.atoms.get_positions()
        f = self.atoms.get_forces()
        self.t += time()
        self.f.append((self.t, (f**2).sum(1).max()))
        return f

    def set_positions(self, pos):
        self.atoms.set_positions(pos)

    def get_positions(self):
        return self.atoms.get_positions()


def run(atoms, name, optimizer, db, fmax=0.05):
    opt = get_optimizer(optimizer)
    wrapper = Wrapper(atoms)
    relax = opt(wrapper, logfile=None)

    tincl = -time()
    try:
        relax.run(fmax=fmax, steps=100, fail=True)
    except Exception:
        traceback.print_exc()
        tincl = np.inf
        texcl = np.inf
    else:
        tincl += time()
        texcl = max(wrapper.e[-1][0], wrapper.f[-1][0])

    if isinstance(atoms, NEB):
        n = len(atoms.images)
        atoms = atoms.images[n // 2]
        
    db.write(atoms,
             optimizer=optimizer,
             test=name,
             nenergy=len(wrapper.e), nforce=len(wrapper.f),
             t=texcl, T=tincl, n=wrapper.n,
             data={'e': np.array(wrapper.e).T,
                   'f': np.array(wrapper.f).T})


def main():
    parser = argparse.ArgumentParser(
        description='Test ASE optimizers')

    parser.add_argument('-o', '--optimizers',
                        help='Comma-separated list of optimizer names.')
    parser.add_argument('-s', '--summary', action='store_true')
    parser.add_argument('-p', '--plot', action='store_true')
    parser.add_argument('-d', '--database', default='tests.db')
    parser.add_argument('tests', nargs='+')

    args = parser.parse_args()

    db = ase.db.connect(args.database)

    if args.summary:
        for test in args.tests:
            summary(list(db.select(test=test, sort='n')))
    else:
        if args.optimizers is None:
            args.optimizers = ','.join(all_optimizers)
        for test in args.tests:
            print(test + ':', end=' ')
            atoms, name = get_atoms_and_name(test)
            p0 = atoms.get_positions()
            for opt in args.optimizers.split(','):
                print(opt, end=' ', flush=True)
                atoms.set_positions(p0)
                if not isinstance(atoms, NEB):
                    atoms.calc = EMT()
                run(atoms, name, opt, db)
            print()


def summary(rows):
    print(rows[0].test + ':')
    e0 = min(row.get('energy', np.inf) for row in rows)
    table = sorted((row.get('energy', np.inf) > e0 + 0.005, row.n, row.t,
                    row.optimizer, row.T,
                    row.nenergy, row.nforce,
                    row.get('energy', np.inf) - e0)
                   for row in rows)
    for _, n, t, opt, T, ne, nf, de in table:
        print('{:20}{:4}{:12.6f}{:12.6f}{:4}{:4}{:12.6f}'.
              format(opt, n, t, T, ne, nf, de))


class Plotter:
    def __init__(self, name, fmax):
        self.name = name
        self.fmax = fmax
        if rank == 0:
            self.fig = pl.figure(figsize=[12.0, 9.0])
            self.axes0 = self.fig.add_subplot(2, 1, 1)
            self.axes1 = self.fig.add_subplot(2, 1, 2)

    def plot(self, optimizer, E, fmax):
        if rank == 0:
            self.axes0.plot(E, label = optimizer)
            self.axes1.plot(fmax)

    def save(self, format='png'):
        if rank == 0:
            self.axes0.legend()
            self.axes0.set_title(self.name)
            self.axes0.set_ylabel('E [eV]')

            self.axes1.set_xlabel('steps')
            self.axes1.set_ylabel('fmax [eV/A]')
            self.axes1.set_yscale('log')
            self.axes1.axhline(self.fmax, color='k', linestyle='--')
            self.fig.savefig(self.name + '.' + format)


if __name__ == '__main__':
    main()
