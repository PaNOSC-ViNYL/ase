from ase.atoms import Atoms
from ase.parallel import paropen

"""Module to read and write atoms in PDB file format"""


def read_pdb(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    positions = []
    symbols = []
    for line in fileobj.readlines():
        if line.startswith('ATOM'):
            words = line.split()
            symbol = ''
            for s in words[2]:
                if not s.isdigit():
                    if len(symbol):
                        symbol += s.lower()
                    else:
                        symbol += s.upper()
            symbols.append(symbol)
            positions.append([float(words[4]), 
                              float(words[5]),
                              float(words[6])])
    return Atoms(symbols=symbols, positions=positions)

def write_pdb(fileobj, images):
    """Write images to PDB-file."""
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    format = 'ATOM  %5d %-4s              %8.3f%8.3f%8.3f  0.00  0.00\n'

    # RasMol complains if the atom index exceeds 100000. There might
    # be a limit of 5 digit numbers in this field.
    MAXNUM = 100000

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    
    for atoms in images:
        fileobj.write('MODEL         1\n')
        p = atoms.get_positions()
        for a in range(natoms):
            x, y, z = p[a]
            fileobj.write(format % (a % MAXNUM, symbols[a], x, y, z))
        fileobj.write('ENDMDL\n')
