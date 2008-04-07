import numpy as npy

from ase.md import MolecularDynamics


class VelocityVerlet(MolecularDynamics):
    def __init__(self, atoms):
        MolecularDynamics.__init__(self, atoms)
            
    def step(self, f, dt):
        atoms = self.atoms
        p = self.atoms.get_momenta()
        p += 0.5 * dt * f
        self.atoms.set_positions(self.atoms.positions + dt * p / self.masses)
        f = self.atoms.get_forces()
        atoms.set_momenta(p + 0.5 * dt * f)
        return f
