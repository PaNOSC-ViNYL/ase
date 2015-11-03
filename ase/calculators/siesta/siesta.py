from __future__ import print_function
"""This module defines an ASE interface to SIESTA.

http://www.uam.es/departamentos/ciencias/fismateriac/siesta
"""
from ase.calculators.siesta.base_siesta import BaseSiesta


# Version 3.2 of Siesta
class Siesta3_2(BaseSiesta):
    allowed_xc = {
        'LDA': ['PZ', 'CA', 'PW92'],
        'GGA': ['PBE', 'revPBE', 'RPBE',
                'WC', 'PBEsol', 'LYP'],
    }


# Trunk version, snapshot 462
class SiestaTrunk462(BaseSiesta):
    allowed_xc = {
        'LDA': ['PZ', 'CA', 'PW92'],
        'GGA': ['PW91', 'PBE', 'revPBE', 'RPBE',
                'WC', 'AM05', 'PBEsol', 'PBEJsJrLO',
                'PBEGcGxLO', 'PBEGcGxHEG', 'BLYP',
                ],
        'VDW': ['DRSLL', 'LMKLL', 'KBM', 'C09', 'BH', 'VV'],
    }


Siesta = Siesta3_2
