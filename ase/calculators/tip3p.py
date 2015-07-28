"""TIP3P potential."""
from __future__ import division

import numpy as np

import ase.units as units
from ase.calculators.calculator import Calculator, all_changes

qH = 0.417
sigma0 = 3.15061
epsilon0 = 0.1521 * units.kcal / units.mol
rOH = 0.9572
thetaHOH = 104.52 / 180 * np.pi


def set_tip3p_charges(atoms):
    charges = np.empty(len(atoms))
    charges[:] = qH
    if atoms.numbers[0] == 8:
        charges[::3] = -2 * qH
    else:
        charges[2::3] = -2 * qH
    atoms.set_initial_charges(charges)
    
    
class TIP3P(Calculator):
    implemented_properties = ['energy', 'forces']
    nolabel = True
    pcp = None
    
    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        
        positions = self.atoms.positions
        Z = self.atoms.numbers

        assert not self.atoms.pbc.any()
        
        energy = 0.0
        forces = np.zeros_like(positions)

        if Z[0] == 8:
            RO = positions[::3]
            FO = forces[::3]
        else:
            RO = positions[2::3]
            FO = forces[2::3]
        
        charges = self.atoms.get_initial_charges()
        
        for a, R in enumerate(RO):
            d = RO[a + 1:] - R
            r2 = (d**2).sum(1)
            c6 = (sigma0**2 / r2)**3
            c12 = c6**2
            energy += 4 * epsilon0 * (c12 - c6).sum()
            f = (24 * epsilon0 * (2 * c12 - c6) / r2)[:, np.newaxis] * d
            FO[a] -= f.sum(0)
            FO[a + 1:] += f
        
        for a, R in enumerate(positions):
            b = a // 3 * 3 + 3
            d = positions[b:] - R
            r2 = (d**2).sum(1)
            e = units.Hartree * units.Bohr * charges[a] * charges[b:] / r2**0.5
            energy += e.sum()
            f = (e / r2)[:, np.newaxis] * d
            forces[a] -= f.sum(0)
            forces[b:] += f
            
        if self.pcp:
            e, f = self.pcp.calculate(charges, positions)
            energy += e
            forces += f
            
        self.results['energy'] = energy
        self.results['forces'] = forces

    def embed(self, charges):
        self.pcp = PointChargePotential(charges)
        return self.pcp
        
    def check_state(self, atoms, tol=1e-15):
        system_changes = Calculator.check_state(self, atoms, tol)
        if self.pcp and self.pcp.positions is not None:
            system_changes.append('positions')
        return system_changes
        
        
class PointChargePotential:
    def __init__(self, charges):
        self.charges = charges
        self.positions = None
    
    def set_positions(self, positions):
        self.positions = positions
    
    def calculate(self, charges, positions):
        energy = 0.0
        self.forces = np.zeros_like(self.positions)
        forces = np.zeros_like(positions)
        for C, R, F in zip(self.charges, self.positions, self.forces):
            d = positions - R
            r2 = (d**2).sum(1)
            e = units.Hartree * units.Bohr * C * r2**-0.5 * charges
            energy += e.sum()
            f = (e / r2)[:, np.newaxis] * d
            forces += f
            F -= f.sum(0)
        self.positions = None
        return energy, forces
    
    def get_forces(self, calc):
        return self.forces
