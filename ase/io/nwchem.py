from cStringIO import StringIO
from ase.atoms import Atoms
from ase.io.xyz import read_xyz

def read_nwchem(filename):
    """Method to read geometry from a nwchem output
    """
    from ase import Atoms, Atom

    if isinstance(filename, str):
        f = open(filename)

    lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].find('XYZ format geometry') >=0:
            natoms = int(lines[i + 2].split()[0])
            string = ''
            for j in range(2, natoms + 4):
                string += lines[i + j]
            atoms = read_xyz(StringIO(string))
            i += natoms + 4
        else:
            i += 1
           
    if type(filename) == str:
        f.close()

    return atoms

def write_nwchem(filename, atoms, autosym=False):
    """Method to write nwchem coord file
    """

    import numpy as np

    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename

    if autosym:
        f.write('geometry\n')
    else:
        f.write('geometry noautosym\n')
    for atom in atoms:
        f.write('  ' + atom.symbol + ' ' +
                str(atom.position[0]) + ' ' +
                str(atom.position[1]) + ' ' +
                str(atom.position[2]) + '\n' )
    f.write('end\n')
    
