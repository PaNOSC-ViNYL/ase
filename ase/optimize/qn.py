import numpy as npy
from numpy.linalg import eigh, solve

from ase.optimize import Optimizer


class QuasiNewton(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', maxstep=None):
        Optimizer.__init__(self, atoms, restart, logfile)

        if maxstep is not None:
            self.maxstep = maxstep

    def initialize(self):
        self.H = None
        self.maxstep = 0.04

    def read(self):
        self.H, self.r0, self.f0, self.maxstep = self.load()

    def step(self, f):
        atoms = self.atoms
        r = atoms.get_positions()
        f = f.reshape(-1)
        self.update(r.flat, f)
        omega, V = eigh(self.H)
        dr = npy.dot(V, npy.dot(f, V) / npy.fabs(omega)).reshape((-1, 3))
        #dr = solve(self.H, f).reshape((-1, 3))
        steplengths = (dr**2).sum(1)**0.5
        
        dr /= npy.maximum(steplengths / self.maxstep, 1.0).reshape(-1, 1)
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))

    def update(self, r, f):
        if self.H is None:
            self.H = npy.eye(3 * len(self.atoms)) * 120.0
            return
        dr = r - self.r0
        df = f - self.f0
        a = npy.dot(dr, df)
        dg = npy.dot(self.H, dr)
        b = npy.dot(dr, dg)
        self.H -= npy.outer(df, df) / a + npy.outer(dg, dg) / b
