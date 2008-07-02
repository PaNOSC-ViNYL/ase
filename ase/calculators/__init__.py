"""Interfaces to different ASE compatible force-calculators."""

import numpy as npy

from ase.calculators.lj import LennardJones
from ase.calculators.emt import EMT, ASAP
from ase.calculators.siesta import Siesta

class Calculator:
    def get_potential_energy(self, atoms):
        return 0.0
    
    def calculation_required(self, atoms, quantities):
        """Check if a calculation is required.

        Check if the quantities in the *quantities* list have already
        been calculated for the atomic configuration *atoms*.  The
        quantities can be one or more of: 'energy', 'forces',
        'stress'."""

        return False


class SinglePointCalculator:
    """Special calculator for a single configuration.

    Used to remember the energy, force and stress for a given
    configuration.  If the positions, atomic numbers, unit cell
    boundary conditions are changed, then asking for
    energy/forces/stress will raise an exception."""
    
    def __init__(self, energy, forces, stress, magmoms, atoms):
        """Save energy, forces and stresses for the current configuration."""
        self.energy = energy
        self.forces = forces
        self.stress = stress
        self.magmoms = magmoms
        self.atoms = atoms.copy()

    def calculation_required(self, atoms, quantities):
        ok = self.atoms.identical_to(atoms)
        return ('forces' in quantities and (self.forces is None or not ok) or
                ('energy' in quantities and (self.energy is None or not ok)))

    def update(self, atoms):
        if not self.atoms.identical_to(atoms):
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

    def get_spin_polarized(self):
        return self.magmoms is not None

    def get_magnetic_moments(self, atoms):
        self.update(atoms)
        if self.magmoms is not None:
            return self.magmoms
        else:
            return npy.zeros(len(self.positions))
        

def numeric_force(atoms, a, i, d=0.001):
    """Evaluate forces usinf finite difference formula."""
    p0 = atoms.positions[a, i]
    atoms.positions[a, i] += d
    eplus = atoms.get_potential_energy()
    atoms.positions[a, i] -= 2 * d
    eminus = atoms.get_potential_energy()
    atoms.positions[a, i] = p0
    return (eminus - eplus) / (2 * d)


class TestPotential:
    def get_forces(self, atoms):
        E = 0.0
        R = atoms.positions
        F = npy.zeros_like(R)
        for a, r in enumerate(R):
            D = R - r
            x = (D**2).sum(1) - 1.0
            E += npy.vdot(x, x)
            F -= x[:, None] * D
        self.energy = 0.25 * E
        return F

    def get_potential_energy(self, atoms):
        self.get_forces(atoms)
        return self.energy

    def get_stress(self, atoms):
        raise NotImplementedError
