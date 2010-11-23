"""van der Waals correction schemes for DFT"""

import numpy as np
from ase.units import Bohr, Hartree

# dipole polarizabilities and C6 values from 
# X. Chu and A. Dalgarno, J. Chem. Phys. 129 (2004) 4083
# atomic units, a_0^3
vdWDB_Chu04jcp = {
    # Element: [alpha, C6]; units [Bohr^3, Hartree * Bohr^6]
    'H'  : [9/2, 6.5], # [exact, Tkatchenko PRL]
    'He' : [1.38, 1.42],
    'Li' : [164, 1392],
    'Be' : [38, 227],
    'B'  : [21, 99.5],
    'C'  : [12, 46.6],
    'N'  : [7.4, 24.2],
    'O'  : [5.4, 15.6],
    'F'  : [3.8, 9.52],
    'Na' : [163, 1518],
    'Fe' : [56, 482],
    'Mg' : [71, 626],
    'Ne' : [2.67, 6.20],
    'Ca' : [160, 2163],
    'Ar' : [11.1, 64.2],
    'Sr' : [199, 3175],
    'Kr' : [16.7, 130],
}

vdWDB_Grimme06jcc = {
    # Element: [C6, R0]; units [J nm^6 mol^{-1}, Angstrom]
    'H'  : [0.14, 1.001],
    'He' : [0.08, 1.012],
    'Li' : [1.61, 0.825],
    'Be' : [1.61, 1.408],
    'B'  : [3.13, 1.485],
    'C'  : [1.75, 1.452],
    'N'  : [1.23, 1.397],
    'O'  : [0.70, 1.342],
    'F'  : [0.75, 1.287],
    'Ne' : [0.63, 1.243],
    'Na' : [5.71, 1.144],
    'Mg' : [5.71, 1.364],
    'Al' : [10.79, 1.639],
    'Si' : [9.23, 1.716],
    'P'  : [7.84, 1.705],
    'S'  : [5.57, 1.683],
    'Cl' : [5.07, 1.639],
    'Ar' : [4.61, 1.595],
    'K'  : [10.80, 1.485],
    'Ca' : [10.80, 1.474],
    'Sc' : [10.80, 1.562],
    'Ti' : [10.80, 1.562],
    'V'  : [10.80, 1.562],
    'Cr'  : [10.80, 1.562],
    'Mn'  : [10.80, 1.562],
    'Fe'  : [10.80, 1.562],
    'Co'  : [10.80, 1.562],
    'Ni'  : [10.80, 1.562],
    'Cu'  : [10.80, 1.562],
    'Zn' : [10.80, 1.562],
    'Ga' : [16.99, 1.650],
    'Ge' : [17.10, 1.727],
    'As' : [16.37, 1.760],
    'Se' : [12.64, 1.771],
    'Br' : [12.47, 1.749],
    'Kr' : [12.01, 1.727],
    'Rb' : [24.67, 1.628],
    'Sr' : [24.67, 1.606],
    'Y-Cd' : [24.67, 1.639],
    'In' : [37.32, 1.672],
    'Sn' : [38.71, 1.804],
    'Sb' : [38.44, 1.881],
    'Te' : [31.74, 1.892],
    'I'  : [31.50, 1.892],
    'Xe' : [29.99, 1.881],
    }

class vdWTkatchenko09prl:
    """vdW correction after Tkatchenko and Scheffler PRL 102 (2009) 073005.

    hirshfeld: the Hirshfeld partitioning object
    calculator: the calculator to get the PBE energy
    missing: Missing elements do not contribute to the vdW-Energy by default
    """
    def __init__(self, hirshfeld=None, calculator=None, missing='zero'):
        self.hirshfeld = hirshfeld
        if calculator is None:
            self.calculator = self.hirshfeld.get_calculator()
        else:
            self.calculator = calculator
        self.missing = missing

    def update(self, atoms=None):
        if atoms is None:
            atoms = self.calculator.get_atoms()
        assert not atoms.get_pbc().any()

        if self.hirshfeld == None:
            volume_ratios = [1.] * len(atoms)
        else:
            volume_ratios = self.hirshfeld.get_effective_volume_ratios()

        positions = atoms.get_positions()
        indicees = range(len(positions))
        symbols = atoms.get_chemical_symbols()
        EvdW = 0.0
        for i1, p1 in enumerate(positions):
            if symbols[i1] in vdWDB_Chu04jcp:
                # free atom values
                alphaA, C6AA = vdWDB_Chu04jcp[symbols[i1]]
                # correction for effective C6
                C6AA *= Hartree * volume_ratios[i1]**2 * Bohr**6

                for i2 in indicees[i1 + 1:]:
                    p2 = positions[i2]
                    if symbols[i2] in vdWDB_Chu04jcp:
                        # free atom values
                        alphaB, C6BB = vdWDB_Chu04jcp[symbols[i2]]
                        # correction for effective C6
                        C6BB *= Hartree * volume_ratios[i2]**2 * Bohr**6

                        C6AB = (2 * C6AA * C6BB /
                                (alphaB / alphaA * C6AA +
                                 alphaA / alphaB * C6BB   ))
 
                        diff = p2 - p1
                        r2 = np.dot(diff, diff) 
                        EvdW -= (self.damping(np.sqrt(r2),
                                              symbols[i1], symbols[i2]) *
                                 C6AB / r2 / r2 /r2)

            elif self.missing != 'zero':
                raise RuntimeError('Element ' + symbols[i1] +
                                   ' not in database.')

        self.energy += EvdW
        
    def damping(self, RAB, symbolA, symbolB,
                d = 20,   # steepness of the step function
                sR = 0.94 # for PBE
                ):
        """Damping factor.

        Standard values for d and sR as given in 
        Tkatchenko and Scheffler PRL 102 (2009) 073005."""
        # XXX use nonempirical radii
        RA0 = vdWDB_Grimme06jcc[symbolA][1]
        RB0 = vdWDB_Grimme06jcc[symbolB][1]
        x = RAB / (sR * (RA0 + RB0))
        return 1.0 / (1.0 + np.exp(-d * (x - 1.0)))
 
    def get_potential_energy(self, atoms=None):
        self.energy = self.calculator.get_potential_energy(atoms)
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        return 0 * atoms.get_positions()
        #raise RuntimeError('Forces are not available')
