# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import eigh, solve

from ase.optimize import Optimizer


class QuasiNewton(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None):
        """Quasi-Newton optimizer.

        Use *maxstep* to set the maximum distance an atom can move per
        iteration (default value is 0.04 Å)."""
        
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

        if maxstep is not None:
            if maxstep > 1.0:
                raise ValueError('You are using a much too large value for ' +
                                 'the maximum step size: %.1f Å' % maxstep)
            self.maxstep = maxstep

    def initialize(self):
        self.H = None
        self.r0 = None
        self.f0 = None
        self.maxstep = 0.04

    def read(self):
        self.H, self.r0, self.f0, self.maxstep = self.load()

    def step(self, f):
        atoms = self.atoms
        r = atoms.get_positions()
        f = f.reshape(-1)
        self.update(r.flat, f)
        omega, V = eigh(self.H)
        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        #dr = solve(self.H, f).reshape((-1, 3))
        steplengths = (dr**2).sum(1)**0.5
        dr /= np.maximum(steplengths / self.maxstep, 1.0).reshape(-1, 1)
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))

    def update(self, r, f):
        if self.H is None:
            self.H = np.eye(3 * len(self.atoms)) * 70.0
            return
        dr = r - self.r0

        if np.abs(dr).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return
        
        df = f - self.f0
        a = np.dot(dr, df)
        dg = np.dot(self.H, dr)
        b = np.dot(dr, dg)
        self.H -= np.outer(df, df) / a + np.outer(dg, dg) / b

    def replay_trajectory(self, traj):
        """Initialize hessian from old trajectory."""
        if isinstance(traj, str):
            traj = PickleTrajectory(traj, 'r')
        r0, f0 = self.r0, self.f0
        self.H = None
        for atoms in traj:
            self.update(atoms.get_positions().ravel(),
                        atoms.get_forces().ravel())
        self.r0, self.f0 = r0, f0

            
