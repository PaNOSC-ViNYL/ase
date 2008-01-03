import numpy as npy

from ase.calculators.lj import LennardJones
from ase.calculators.emt import EMT, ASAP


class SinglePointCalculator:
    def __init__(self, energy, forces, stress, atoms):
        self.energy = energy
        self.forces = forces
        self.stress = stress
        self.positions = atoms.get_positions().copy()
        self.numbers = atoms.get_atomic_numbers().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()
        
    def update(self, atoms):
        if ((self.positions != atoms.get_positions()).any() or
            (self.numbers != atoms.get_atomic_numbers()).any() or
            (self.cell != atoms.get_cell()).any() or
            (self.pbc != atoms.get_pbc()).any()):
            raise RuntimeError('Energy, forces and strees no longer correct.')

    def get_potential_energy(self, atoms):
        self.update(atoms)
        if self.energy is None:
            raise RuntimeError('No energy.')
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        if self.forces is None:
            raise RuntimeError('No forces.')
        return self.forces

    def get_stress(self, atoms):
        self.update(atoms)
        if self.stress is None:
            raise NotImplementedError
        return self.stress


def numeric_force(atoms, a, i, d=0.001):
    p0 = atoms.positions[a, i]
    atoms.positions[a, i] += d
    eplus = atoms.GetPotentialEnergy()
    atoms.positions[a, i] -= 2 * d
    eminus = atoms.GetPotentialEnergy()
    atoms.positions[a, i] = p0
    return (eminus - eplus) / (2 * d)

