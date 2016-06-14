from __future__ import print_function
"""This module defines an ASE interface to deMon.

http://www.demon-software.com 
"""

from ase.calculators.deMon.base_deMon import Base_deMon


# Version 4.3.2 of deMon2k
class DeMon2k_4_3_2(Base_deMon):
    
    allowed_keywords={}

# Define the default siesta version.
DeMon = DeMon2k_4_3_2
