import numpy as npy

from ase.optimize import Dynamics
from ase.data import atomic_masses


class MolecularDynamics(Dynamics):
    def __init__(self, atoms):
        Dynamics.__init__(self, atoms, logfile=None)

        try:
            self.masses = self.atoms.get_masses()
        except KeyError:
            self.masses = atomic_masses[self.atoms.get_atomic_numbers()]
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

        try:
            self.atoms.get_momenta()
        except KeyError:
            self.atoms.set_momenta(npy.zeros_like(f))

        for step in xrange(steps):
            f = self.step(f, dt)
            self.call_observers(step)
