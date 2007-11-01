from math import sqrt, exp, log, pi

import numpy as npy

from ase.data import atomic_numbers, chemical_symbols
from ase.units import Bohr


#                    E0     s0    V0     eta2    kappa   lambda  n0
#                    eV     bohr  eV     bohr^-1 bohr^-1 bohr^-1 bohr^-3
parameters = {'Cu': (-3.51, 2.67, 2.476, 1.652,  2.740,  1.906,  0.00910),
              'Ag': (-2.96, 3.01, 2.132, 1.652,  2.790,  1.892,  0.00547)}

beta = 1.809#(16 * pi / 3)**(1.0 / 3) / 2**0.5
eta1 = 0.5 / Bohr
acut = 50.0

class EMT:

    acut = 5.9
    
    def __init__(self):
        self.energy = None
        
    def initialize(self, atoms):
        self.par = {}
        self.rc = 0.0
        self.numbers = atoms.get_atomic_numbers()
        for Z in self.numbers:
            if Z not in self.par:
                p = parameters[chemical_symbols[Z]]
                s0 = p[1] * Bohr
                eta2 = p[3] / Bohr
                kappa = p[4] / Bohr
                rc = beta * s0 * 0.5 * (sqrt(3) + sqrt(4))
                x = eta2 * beta * s0
                gamma1 = 0.0
                gamma2 = 0.0
                for i, n in enumerate([12, 6, 24, 8]):
                    r = s0 * beta * sqrt(i + 1)
                    x = n / (12 * (1.0 + exp(acut * (r - rc))))
                    gamma1 += x * exp(-eta2 * (r - beta * s0))
                    gamma2 += x * exp(-kappa / beta * (r - beta * s0))
                self.par[Z] = {'E0': p[0],
                               's0': s0,
                               'V0': p[2],
                               'eta2': eta2,
                               'kappa': kappa,
                               'lambda': p[5] / Bohr,
                               'n0': p[6] / Bohr**3,
                               'rc': rc,
                               'gamma1': gamma1,
                               'gamma2': gamma2}
                if rc > self.rc:
                    self.rc = rc

        self.ksi = {}
        for s1, p1 in self.par.items():
            self.ksi[s1] = {}
            for s2, p2 in self.par.items():
                self.ksi[s1][s2] = (p2['n0'] / p1['n0'] *
                                    exp(eta1 * (p1['s0'] - p2['s0'])))
                
        self.forces = npy.empty((len(atoms), 3))
        self.sigma1 = npy.empty(len(atoms))
                    
    def update(self, atoms):
        if (self.energy is None or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)
                
    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        raise NotImplementedError
    
    def calculate(self, atoms):
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()
        
        icell = npy.linalg.inv(self.cell)
        scaled = npy.dot(self.positions, icell)
        N = []
        for i in range(3):
            if self.pbc[i]:
                scaled[:, i] %= 1.0
                v = icell[:, i]
                h = 1 / sqrt(npy.dot(v, v))
                N.append(int(self.rc / h) + 1)
            else:
                N.append(0)

        R = npy.dot(scaled, self.cell)
        
        self.energy = 0.0
        self.sigma1[:] = 0.0
        
        N1, N2, N3 = N
        natoms = len(atoms)
        for i1 in range(-N1, N1 + 1):
            for i2 in range(-N2, N2 + 1):
                for i3 in range(-N3, N3 + 1):
                    C = npy.dot((i1, i2, i3), self.cell)
                    Q = R + C
                    c = (i1 == 0 and i2 == 0 and i3 == 0)
                    for a1 in range(natoms):
                        Z1 = self.numbers[a1]
                        p1 = self.par[Z1]
                        ksi = self.ksi[Z1]
                        for a2 in range(natoms):
                            if c and a2 == a1:
                                continue
                            d = Q[a2] - R[a1]
                            r = sqrt(npy.dot(d, d))
                            if r < p1['rc']:
                                Z2 = self.numbers[a2]
                                self.interact(a1, r, p1, ksi[Z2])
        for a in range(natoms):
            Z = self.numbers[a]
            p = self.par[Z]
            ds = -log(self.sigma1[a] / 12) / (beta * p['eta2'])
            x = p['lambda'] * ds
            self.energy += (p['E0'] * ((1 + x) * exp(-x) - 1) +
                            6 * p['V0'] * exp(-p['kappa'] * ds))

    def interact(self, a, r, p, ksi):
        theta = 1.0 / (1.0 + exp(acut * (r - p['rc'])))
        self.energy -= (0.5 * p['V0'] *
                        exp(-p['kappa'] * (r / beta - p['s0'])) *
                        ksi * theta / p['gamma2'])
        self.sigma1[a] += (exp(-p['eta2'] * (r - beta * p['s0'])) *
                           ksi * theta / p['gamma1'])


