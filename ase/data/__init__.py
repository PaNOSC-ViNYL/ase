# -*- coding: utf-8 -*-
import numpy as np
from ase.data.vdw import vdw_radii

__all__ = ['vdw_radii', 'chemical_symbols', 'ground_state_magnetic_moments',
           'reference_states', 'atomic_names', 'atomic_masses',
           'atomic_numbers']

chemical_symbols = ['X',  'H',  'He', 'Li', 'Be',
                    'B',  'C',  'N',  'O',  'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si',
                    'P',  'S',  'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V',  'Cr',
                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se',
                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I',  'Xe',
                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
                    'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U',  'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']

atomic_numbers = {}
for Z, symbol in enumerate(chemical_symbols):
    atomic_numbers[symbol] = Z

atomic_names = [
    '', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium', 'Unnilquadium', 'Unnilpentium', 'Unnilhexium']

# Atomic masses are based on:
#
#   Meija, J., Coplen, T., Berglund, M., et al. (2016). Atomic weights of
#   the elements 2013 (IUPAC Technical Report). Pure and Applied Chemistry,
#   88(3), pp. 265-291. Retrieved 30 Nov. 2016,
#   from doi:10.1515/pac-2015-0305
#
# Standard atomic weights are taken from Table 1: "Standard atomic weights
# 2013", with the uncertainties ignored.
# For hydrogen, helium, boron, carbon, nitrogen, oxygen, magnesium, silicon,
# sulfur, chlorine, bromine and thallium, where the weights are given as a
# range the "conventional" weights are taken from Table 3 and the ranges are
# given in the comments.
atomic_masses = np.array([
    0.0,  # X
    1.008,  # H [1.00784, 1.00811]
    4.002602,  # He
    6.94,  # Li [6.938, 6.997]
    9.0121831,  # Be
    10.81,  # B [10.806, 10.821]
    12.011,  # C [12.0096, 12.0116]
    14.007,  # N [14.00643, 14.00728]
    15.999,  # O [15.99903, 15.99977]
    18.998403163,  # F
    20.1797,  # Ne
    22.98976928,  # Na
    24.305,  # Mg [24.304, 24.307]
    26.9815385,  # Al
    28.085,  # Si [28.084, 28.086]
    30.973761998,  # P
    32.06,  # S [32.059, 32.076]
    35.45,  # Cl [35.446, 35.457]
    39.948,  # Ar
    39.0983,  # K
    40.078,  # Ca
    44.955908,  # Sc
    47.867,  # Ti
    50.9415,  # V
    51.9961,  # Cr
    54.938044,  # Mn
    55.845,  # Fe
    58.933194,  # Co
    58.6934,  # Ni
    63.546,  # Cu
    65.38,  # Zn
    69.723,  # Ga
    72.630,  # Ge
    74.921595,  # As
    78.971,  # Se
    79.904,  # Br [79.901, 79.907]
    83.798,  # Kr
    85.4678,  # Rb
    87.62,  # Sr
    88.90584,  # Y
    91.224,  # Zr
    92.90637,  # Nb
    95.95,  # Mo
    np.nan,  # Tc
    101.07,  # Ru
    102.90550,  # Rh
    106.42,  # Pd
    107.8682,  # Ag
    112.414,  # Cd
    114.818,  # In
    118.710,  # Sn
    121.760,  # Sb
    127.60,  # Te
    126.90447,  # I
    131.293,  # Xe
    132.90545196,  # Cs
    137.327,  # Ba
    138.90547,  # La
    140.116,  # Ce
    140.90766,  # Pr
    144.242,  # Nd
    np.nan,  # Pm
    150.36,  # Sm
    151.964,  # Eu
    157.25,  # Gd
    158.92535,  # Tb
    162.500,  # Dy
    164.93033,  # Ho
    167.259,  # Er
    168.93422,  # Tm
    173.054,  # Yb
    174.9668,  # Lu
    178.49,  # Hf
    180.94788,  # Ta
    183.84,  # W
    186.207,  # Re
    190.23,  # Os
    192.217,  # Ir
    195.084,  # Pt
    196.966569,  # Au
    200.592,  # Hg
    204.38,  # Tl [204.382, 204.385]
    207.2,  # Pb
    208.98040,  # Bi
    np.nan,  # Po
    np.nan,  # At
    np.nan,  # Rn
    np.nan,  # Fr
    np.nan,  # Ra
    np.nan,  # Ac
    232.0377,  # Th
    231.03588,  # Pa
    238.02891,  # U
    np.nan,  # Np
    np.nan,  # Pu
    np.nan,  # Am
    np.nan,  # Cm
    np.nan,  # Bk
    np.nan,  # Cf
    np.nan,  # Es
    np.nan,  # Fm
    np.nan,  # Md
    np.nan,  # No
    np.nan])  # Lw

# Covalent radii from:
#
#  Covalent radii revisited,
#  Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés,
#  Jorge Echeverría, Eduard Cremades, Flavia Barragán and Santiago Alvarez,
#  Dalton Trans., 2008, 2832-2838 DOI:10.1039/B801115J
missing = 0.2
covalent_radii = np.array([
    missing,  # X
    0.31,  # H
    0.28,  # He
    1.28,  # Li
    0.96,  # Be
    0.84,  # B
    0.76,  # C
    0.71,  # N
    0.66,  # O
    0.57,  # F
    0.58,  # Ne
    1.66,  # Na
    1.41,  # Mg
    1.21,  # Al
    1.11,  # Si
    1.07,  # P
    1.05,  # S
    1.02,  # Cl
    1.06,  # Ar
    2.03,  # K
    1.76,  # Ca
    1.70,  # Sc
    1.60,  # Ti
    1.53,  # V
    1.39,  # Cr
    1.39,  # Mn
    1.32,  # Fe
    1.26,  # Co
    1.24,  # Ni
    1.32,  # Cu
    1.22,  # Zn
    1.22,  # Ga
    1.20,  # Ge
    1.19,  # As
    1.20,  # Se
    1.20,  # Br
    1.16,  # Kr
    2.20,  # Rb
    1.95,  # Sr
    1.90,  # Y
    1.75,  # Zr
    1.64,  # Nb
    1.54,  # Mo
    1.47,  # Tc
    1.46,  # Ru
    1.42,  # Rh
    1.39,  # Pd
    1.45,  # Ag
    1.44,  # Cd
    1.42,  # In
    1.39,  # Sn
    1.39,  # Sb
    1.38,  # Te
    1.39,  # I
    1.40,  # Xe
    2.44,  # Cs
    2.15,  # Ba
    2.07,  # La
    2.04,  # Ce
    2.03,  # Pr
    2.01,  # Nd
    1.99,  # Pm
    1.98,  # Sm
    1.98,  # Eu
    1.96,  # Gd
    1.94,  # Tb
    1.92,  # Dy
    1.92,  # Ho
    1.89,  # Er
    1.90,  # Tm
    1.87,  # Yb
    1.87,  # Lu
    1.75,  # Hf
    1.70,  # Ta
    1.62,  # W
    1.51,  # Re
    1.44,  # Os
    1.41,  # Ir
    1.36,  # Pt
    1.36,  # Au
    1.32,  # Hg
    1.45,  # Tl
    1.46,  # Pb
    1.48,  # Bi
    1.40,  # Po
    1.50,  # At
    1.50,  # Rn
    2.60,  # Fr
    2.21,  # Ra
    2.15,  # Ac
    2.06,  # Th
    2.00,  # Pa
    1.96,  # U
    1.90,  # Np
    1.87,  # Pu
    1.80,  # Am
    1.69,  # Cm
    missing,  # Bk
    missing,  # Cf
    missing,  # Es
    missing,  # Fm
    missing,  # Md
    missing,  # No
    missing,  # Lr
    ])

# This data is from Ashcroft and Mermin.
reference_states = [\
    None, #X
    {'symmetry': 'diatom', 'd': 0.74}, #H
    {'symmetry': 'atom'}, #He
    {'symmetry': 'bcc', 'a': 3.49}, #Li
    {'symmetry': 'hcp', 'c/a': 1.567, 'a': 2.29}, #Be
    {'symmetry': 'tetragonal', 'c/a': 0.576, 'a': 8.73}, #B
    {'symmetry': 'diamond', 'a': 3.57},#C
    {'symmetry': 'diatom', 'd': 1.10},#N
    {'symmetry': 'diatom', 'd': 1.21},#O
    {'symmetry': 'diatom', 'd': 1.42},#F
    {'symmetry': 'fcc', 'a': 4.43},#Ne
    {'symmetry': 'bcc', 'a': 4.23},#Na
    {'symmetry': 'hcp', 'c/a': 1.624, 'a': 3.21},#Mg
    {'symmetry': 'fcc', 'a': 4.05},#Al
    {'symmetry': 'diamond', 'a': 5.43},#Si
    {'symmetry': 'cubic', 'a': 7.17},#P
    {'symmetry': 'orthorhombic', 'c/a': 2.339, 'a': 10.47,'b/a': 1.229},#S
    {'symmetry': 'orthorhombic', 'c/a': 1.324, 'a': 6.24, 'b/a': 0.718},#Cl
    {'symmetry': 'fcc', 'a': 5.26},#Ar
    {'symmetry': 'bcc', 'a': 5.23},#K
    {'symmetry': 'fcc', 'a': 5.58},#Ca
    {'symmetry': 'hcp', 'c/a': 1.594, 'a': 3.31},#Sc
    {'symmetry': 'hcp', 'c/a': 1.588, 'a': 2.95},#Ti
    {'symmetry': 'bcc', 'a': 3.02},#V
    {'symmetry': 'bcc', 'a': 2.88},#Cr
    {'symmetry': 'cubic', 'a': 8.89},#Mn
    {'symmetry': 'bcc', 'a': 2.87},#Fe
    {'symmetry': 'hcp', 'c/a': 1.622, 'a': 2.51},#Co
    {'symmetry': 'fcc', 'a': 3.52},#Ni
    {'symmetry': 'fcc', 'a': 3.61},#Cu
    {'symmetry': 'hcp', 'c/a': 1.856, 'a': 2.66},#Zn
    {'symmetry': 'orthorhombic', 'c/a': 1.695, 'a': 4.51, 'b/a': 1.001},#Ga
    {'symmetry': 'diamond', 'a': 5.66},#Ge
    {'symmetry': 'rhombohedral', 'a': 4.13, 'alpha': 54.10},#As
    {'symmetry': 'hcp', 'c/a': 1.136, 'a': 4.36},#Se
    {'symmetry': 'orthorhombic', 'c/a': 1.307, 'a': 6.67, 'b/a': 0.672},#Br
    {'symmetry': 'fcc', 'a': 5.72},#Kr
    {'symmetry': 'bcc', 'a': 5.59},#Rb
    {'symmetry': 'fcc', 'a': 6.08},#Sr
    {'symmetry': 'hcp', 'c/a': 1.571, 'a': 3.65},#Y
    {'symmetry': 'hcp', 'c/a': 1.593, 'a': 3.23},#Zr
    {'symmetry': 'bcc', 'a': 3.30},#Nb
    {'symmetry': 'bcc', 'a': 3.15},#Mo
    {'symmetry': 'hcp', 'c/a': 1.604, 'a': 2.74},#Tc
    {'symmetry': 'hcp', 'c/a': 1.584, 'a': 2.70},#Ru
    {'symmetry': 'fcc', 'a': 3.80},#Rh
    {'symmetry': 'fcc', 'a': 3.89},#Pd
    {'symmetry': 'fcc', 'a': 4.09},#Ag
    {'symmetry': 'hcp', 'c/a': 1.886, 'a': 2.98},#Cd
    {'symmetry': 'tetragonal', 'c/a': 1.076, 'a': 4.59},#In
    {'symmetry': 'tetragonal', 'c/a': 0.546, 'a': 5.82},#Sn
    {'symmetry': 'rhombohedral', 'a': 4.51, 'alpha': 57.60},#Sb
    {'symmetry': 'hcp', 'c/a': 1.330, 'a': 4.45},#Te
    {'symmetry': 'orthorhombic', 'c/a': 1.347, 'a': 7.27, 'b/a': 0.659},#I
    {'symmetry': 'fcc', 'a': 6.20},#Xe
    {'symmetry': 'bcc', 'a': 6.05},#Cs
    {'symmetry': 'bcc', 'a': 5.02},#Ba
    {'symmetry': 'hcp', 'c/a': 1.619, 'a': 3.75},#La
    {'symmetry': 'fcc', 'a': 5.16},#Ce
    {'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.67},#Pr
    {'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.66},#Nd
    None,#Pm
    {'symmetry': 'rhombohedral', 'a': 9.00, 'alpha': 23.13},#Sm
    {'symmetry': 'bcc', 'a': 4.61},#Eu
    {'symmetry': 'hcp', 'c/a': 1.588, 'a': 3.64},#Gd
    {'symmetry': 'hcp', 'c/a': 1.581, 'a': 3.60},#Th
    {'symmetry': 'hcp', 'c/a': 1.573, 'a': 3.59},#Dy
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.58},#Ho
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.56},#Er
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.54},#Tm
    {'symmetry': 'fcc', 'a': 5.49},#Yb
    {'symmetry': 'hcp', 'c/a': 1.585, 'a': 3.51},#Lu
    {'symmetry': 'hcp', 'c/a': 1.582, 'a': 3.20},#Hf
    {'symmetry': 'bcc', 'a': 3.31},#Ta
    {'symmetry': 'bcc', 'a': 3.16},#W
    {'symmetry': 'hcp', 'c/a': 1.615, 'a': 2.76},#Re
    {'symmetry': 'hcp', 'c/a': 1.579, 'a': 2.74},#Os
    {'symmetry': 'fcc', 'a': 3.84},#Ir
    {'symmetry': 'fcc', 'a': 3.92},#Pt
    {'symmetry': 'fcc', 'a': 4.08},#Au
    {'symmetry': 'rhombohedral', 'a': 2.99, 'alpha': 70.45},#Hg
    {'symmetry': 'hcp', 'c/a': 1.599, 'a': 3.46},#Tl
    {'symmetry': 'fcc', 'a': 4.95},#Pb
    {'symmetry': 'rhombohedral', 'a': 4.75, 'alpha': 57.14},#Bi
    {'symmetry': 'sc', 'a': 3.35},#Po
    None,#At
    None,#Rn
    None,#Fr
    None,#Ra
    {'symmetry': 'fcc', 'a': 5.31},#Ac
    {'symmetry': 'fcc', 'a': 5.08},#Th
    {'symmetry': 'tetragonal', 'c/a': 0.825, 'a': 3.92},#Pa
    {'symmetry': 'orthorhombic', 'c/a': 2.056, 'a': 2.85, 'b/a': 1.736},#U
    {'symmetry': 'orthorhombic', 'c/a': 1.411, 'a': 4.72, 'b/a': 1.035},#Np
    {'symmetry': 'monoclinic'},#Pu
    None,#Am
    None,#Cm
    None,#Bk
    None,#Cf
    None,#Es
    None,#Fm
    None,#Md
    None,#No
    None]#Lw

# http://www.webelements.com
ground_state_magnetic_moments = np.array([
   0.0, # X
   1.0, # H
   0.0, # He
   1.0, # Li
   0.0, # Be
   1.0, # B
   2.0, # C
   3.0, # N
   2.0, # O
   1.0, # F
   0.0, # Ne
   1.0, # Na
   0.0, # Mg
   1.0, # Al
   2.0, # Si
   3.0, # P
   2.0, # S
   1.0, # Cl
   0.0, # Ar
   1.0, # K
   0.0, # Ca
   1.0, # Sc
   2.0, # Ti
   3.0, # V
   6.0, # Cr
   5.0, # Mn
   4.0, # Fe
   3.0, # Co
   2.0, # Ni
   1.0, # Cu
   0.0, # Zn
   1.0, # Ga
   2.0, # Ge
   3.0, # As
   2.0, # Se
   1.0, # Br
   0.0, # Kr
   1.0, # Rb
   0.0, # Sr
   1.0, # Y
   2.0, # Zr
   5.0, # Nb
   6.0, # Mo
   5.0, # Tc
   4.0, # Ru
   3.0, # Rh
   0.0, # Pd
   1.0, # Ag
   0.0, # Cd
   1.0, # In
   2.0, # Sn
   3.0, # Sb
   2.0, # Te
   1.0, # I
   0.0, # Xe
   1.0, # Cs
   0.0, # Ba
   1.0, # La
   1.0, # Ce
   3.0, # Pr
   4.0, # Nd
   5.0, # Pm
   6.0, # Sm
   7.0, # Eu
   8.0, # Gd
   5.0, # Tb
   4.0, # Dy
   3.0, # Ho
   2.0, # Er
   1.0, # Tm
   0.0, # Yb
   1.0, # Lu
   2.0, # Hf
   3.0, # Ta
   4.0, # W
   5.0, # Re
   4.0, # Os
   3.0, # Ir
   2.0, # Pt
   1.0, # Au
   0.0, # Hg
   1.0, # Tl
   2.0, # Pb
   3.0, # Bi
   2.0, # Po
   1.0, # At
   0.0, # Rn
   1.0, # Fr
   0.0, # Ra
   1.0, # Ac
   2.0, # Th
   3.0, # Pa
   4.0, # U
   5.0, # Np
   6.0, # Pu
   7.0, # Am
   8.0, # Cm
   5.0, # Bk
   4.0, # Cf
   4.0, # Es
   2.0, # Fm
   1.0, # Md
   0.0, # No
np.nan])# Lw
