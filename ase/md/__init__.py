"""Molecular Dynamics."""

import numpy as npy

from ase.optimize import Dynamics
from ase.data import atomic_masses


class MolecularDynamics(Dynamics):
    """Base-class for all MD classes."""
    def __init__(self, atoms):
        Dynamics.__init__(self, atoms, logfile=None)

        self.masses = self.atoms.get_masses()
        self.masses.shape = (-1, 1)

    def run(self, dt, steps=50):
        """Integrate equation of motion.

        Parameters
        ----------
        dt: float
            Time step.
        steps: int
            Number of steps (defaults to 50).
        """
        
        f = self.atoms.get_forces()

        if not self.atoms.has('momenta'):
            self.atoms.set_momenta(npy.zeros_like(f))

        for step in xrange(steps):
            f = self.step(f, dt)
            self.call_observers()
