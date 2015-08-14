import numpy as np

from ase.calculators.calculator import Calculator
from ase.data import atomic_numbers


class QMMM1(Calculator):
    implemented_properties = ['energy', 'forces']
    
    def __init__(self, selection, qmcalc, mmcalc1, mmcalc2, vacuum=None):
        self.selection = selection
        self.qmcalc = qmcalc
        self.mmcalc1 = mmcalc1
        self.mmcalc2 = mmcalc2
        self.vacuum = vacuum
        
        self.qmatoms = None
        self.center = None
        
        self.name = '{0}-{1}+{1}'.format(qmcalc.name, mmcalc1.name)
        
        Calculator.__init__(self)
        
    def initialize_qm(self, atoms):
        constraints = atoms.constraints
        atoms.constraints = []
        self.qmatoms = atoms[self.selection]
        atoms.constraints = constraints
        self.qmatoms.pbc = False
        if self.vacuum:
            self.qmatoms.center(vacuum=self.vacuum)
            self.center = self.qmatoms.positions.mean(axis=0)
            
    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        
        if self.qmatoms is None:
            self.initialize_qm(atoms)
            
        self.qmatoms.positions = atoms.positions[self.selection]
        if self.vacuum:
            self.qmatoms.positions += (self.center -
                                       self.qmatoms.positions.mean(axis=0))
            
        energy = self.mmcalc2.get_potential_energy(atoms)
        forces = self.mmcalc2.get_forces(atoms)
        
        energy += self.qmcalc.get_potential_energy(self.qmatoms)
        qmforces = self.qmcalc.get_forces(self.qmatoms)
        if self.vacuum:
            qmforces -= qmforces.mean(axis=0)
        forces[self.selection] += qmforces
        
        energy -= self.mmcalc1.get_potential_energy(self.qmatoms)
        forces[self.selection] -= self.mmcalc1.get_forces(self.qmatoms)

        self.results['energy'] = energy
        self.results['forces'] = forces

        
class QMMM2(Calculator):
    implemented_properties = ['energy', 'forces']
    
    def __init__(self, selection, qmcalc, mmcalc, interaction, vacuum=None):
        self.selection = selection

        self.qmcalc = qmcalc
        self.mmcalc = mmcalc
        self.interaction = interaction
        
        self.vacuum = vacuum
        
        self.qmatoms = None
        self.mmatoms = None
        self.mask = None
        self.center = None
        
        self.name = '{0}+{1}+{2}'.format(qmcalc.name,
                                         interaction.name,
                                         mmcalc.name)
        
        Calculator.__init__(self)
        
    def initialize(self, atoms):
        self.mask = np.zeros(len(atoms), bool)
        self.mask[self.selection] = True
        
        constraints = atoms.constraints  # avoid slicing of constraints
        atoms.constraints = []
        self.qmatoms = atoms[self.mask]
        self.mmatoms = atoms[~self.mask]
        atoms.constraints = constraints
        
        self.qmatoms.pbc = False
        
        if self.vacuum:
            self.qmatoms.center(vacuum=self.vacuum)
            self.center = self.qmatoms.positions.mean(axis=0)
            
        self.qmatoms.calc = self.qmcalc
        self.mmatoms.calc = self.mmcalc
        self.pc = self.qmcalc.embed(self.mmatoms.get_initial_charges())
        
    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        
        if self.qmatoms is None:
            self.initialize(atoms)
            
        self.mmatoms.set_positions(atoms.positions[~self.mask])
        self.qmatoms.set_positions(atoms.positions[self.mask])
        
        if self.vacuum:
            shift = self.center - self.qmatoms.positions.mean(axis=0)
            self.qmatoms.positions += shift
        else:
            shift = (0, 0, 0)
            
        self.pc.set_positions(self.mmatoms.positions + shift)

        ienergy, iqmforces, immforces = self.interaction.calculate(
            self.qmatoms, self.mmatoms, shift)
        
        energy = (ienergy +
                  self.qmatoms.get_potential_energy() +
                  self.mmatoms.get_potential_energy())
        
        qmforces = self.qmatoms.get_forces()
        mmforces = self.mmatoms.get_forces()
        
        mmforces += self.pc.get_forces(self.qmatoms.calc)
        
        forces = np.empty((len(atoms), 3))
        forces[self.mask] = qmforces + iqmforces
        forces[~self.mask] = mmforces + immforces

        self.results['energy'] = energy
        self.results['forces'] = forces

        
class LJInteractions:
    name = 'LJ'
    
    def __init__(self, parameters):
        self.parameters = {}
        for (symbol1, symbol2), (epsilon, sigma) in parameters.items():
            Z1 = atomic_numbers[symbol1]
            Z2 = atomic_numbers[symbol2]
            self.parameters[(Z1, Z2)] = epsilon, sigma
            self.parameters[(Z2, Z1)] = epsilon, sigma
        
    def calculate(self, qmatoms, mmatoms, shift):
        qmforces = np.zeros_like(qmatoms.positions)
        mmforces = np.zeros_like(mmatoms.positions)
        species = set(mmatoms.numbers)
        energy = 0.0
        for R1, Z1, F1 in zip(qmatoms.positions, qmatoms.numbers, qmforces):
            for Z2 in species:
                if (Z1, Z2) not in self.parameters:
                    continue
                epsilon, sigma = self.parameters[(Z1, Z2)]
                mask = (mmatoms.numbers == Z2)
                d = mmatoms.positions[mask] + shift - R1
                r2 = (d**2).sum(1)
                s6 = (sigma**2 / r2)**3
                s12 = s6**2
                energy += 4 * epsilon * (s12 - s6).sum()
                f = 24 * epsilon * ((2 * s12 - s6) / r2)[:, np.newaxis] * d
                F1 -= f.sum(0)
                mmforces[mask] += f
        return energy, qmforces, mmforces
