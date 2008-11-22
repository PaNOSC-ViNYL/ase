import numpy as npy

from ase.md import MolecularDynamics


class VelocityVerlet(MolecularDynamics):
    def __init__(self, atoms, dt, trajectory=None):
        MolecularDynamics.__init__(self, atoms, dt, trajectory)
            
    def step(self, f):
        atoms = self.atoms
        p = self.atoms.get_momenta()
        p += 0.5 * self.dt * f
        self.atoms.set_positions(self.atoms.get_positions() +
                                 self.dt * p / self.masses)
        f = self.atoms.get_forces()
        atoms.set_momenta(p + 0.5 * self.dt * f)
        return f
