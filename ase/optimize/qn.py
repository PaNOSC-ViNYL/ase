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
        self.maxstep = 0.02

    def read(self):
        self.H, self.r0, self.f0, self.maxstep = self.load()

    def step(self, f):
        atoms = self.atoms
        r = atoms.get_positions()
        f = f.reshape(-1)
        self.update(r.flat, f)
        omega, V = eigh(self.H)
        dr = npy.dot(npy.dot(V, f) / npy.fabs(omega), V).reshape((-1, 3))
        #dr = solve(self.H, f).reshape((-1, 3))
        steplengths = (dr**2).sum(1)**0.5
        from gpaw.mpi import rank
        print rank, steplengths, npy.maximum(steplengths / self.maxstep, 1.0)
        
        dr /= npy.maximum(steplengths / self.maxstep, 1.0).reshape(-1, 1)
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f
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
class QN2(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', maxstep=None):
        Optimizer.__init__(self, atoms, restart, logfile)

        if maxstep is not None:
            self.maxstep = maxstep

    def initialize(self):
        self.H = None
        self.maxstep = 0.02

    def read(self):
        self.H, self.r0, self.f0, self.maxstep = self.load()

    def step(self, f):
        atoms = self.atoms
        r = atoms.get_positions()
        f = f.reshape(-1)
        self.update(r.flat, f)
        from gpaw.mpi import rank
        print 'RANK', rank, self.H[-2:,-2:]
        omega, V = eigh(self.H)
        print 'RANK', rank, omega, V[-2:,-2:]
        print 'RANK', rank, f
        dr = npy.dot(npy.dot(V, f) / npy.fabs(omega), V).reshape((-1, 3))
        print 'RANK', rank, dr[-2:]
        #dr = solve(self.H, f).reshape((-1, 3))
        steplengths = (dr**2).sum(1)**0.5
        print 'RANK', rank, steplengths[-2:]
        
        dr /= npy.maximum(steplengths / self.maxstep, 1.0).reshape(-1, 1)
        print 'RANK', rank, dr[-2:]
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f
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
