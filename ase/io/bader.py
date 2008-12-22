import re
from math import sqrt
import numpy as np
from ase.units import Bohr

def attach_charges(atoms, fileobj='ACF.dat', displacement=1e-4):
    """Attach the charges from the fileobj to the Atoms."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)
    lines = fileobj.readlines()
    fileobj.close()

    nsep = 0
    charges = []
    while(lines):
        line = lines.pop(0)
        search = re.search('---------------', line)
        if search is not None:
            nsep += 1
        elif nsep == 1:
            words = line.split()
            xyz = np.array([float(w) for w in words[1:4]]) * Bohr
            
            # check if the atom positions match
            if displacement is not None:
                d = xyz - atoms[len(charges)].position
                d = sqrt((d * d).sum())
                assert(d < displacement)
            
            Z = atoms[len(charges)].number
            charges.append(Z - float(words[4]))

    atoms.set_charges(charges)
   
