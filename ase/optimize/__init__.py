import sys
import pickle
from os.path import isfile

import numpy as npy


class Optimizer:
    def __init__(self, atoms, restart, logfile):
        self.atoms = atoms
        dir(atoms)
        self.restart = restart

        if restart is None or not isfile(restart):
            self.initialize()
        else:
            self.read()

        if isinstance(logfile, str):
            if logfile == '-':
                logfile = sys.stdout
            else:
                logfile = open(logfile, 'a')
        self.logfile = logfile
        
        self.callbacks = []

    def attach(self, callback):
        self.callbacks.append(callback)

    def run(self, fmax=0.05, steps=100000000):
        self.fmax = fmax
        for step in xrange(steps):
            f = self.atoms.get_forces()
            self.log(f, step)
            if self.converged(f):
                return
            self.step(f)
            for callback in self.callbacks:
                callback()

    def converged(self, forces=None):
        if forces is None:
            forces = self.atoms.get_forces()
        return (forces**2).sum(axis=1).max() < self.fmax**2

    def log(self, forces, step):
        if self.log is None:
            return
        fmax = (forces**2).sum(axis=1).max()
        e = self.atoms.get_potential_energy()
        name = self.__class__.__name__
        self.logfile.write('%s: %3d %15.6f %12.4f\n' % (name, step, e, fmax))
        self.logfile.flush()
        
    def dump(self, data):
        if self.restart is not None:
            pickle.dump(data, open(self.restart, 'wb'), protocol=2)

    def load(self):
        return pickle.load(open(self.restart))
