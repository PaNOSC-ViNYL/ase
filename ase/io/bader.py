import numpy as np
from ase.units import Bohr

def attach_charges(atoms, fileobj='ACF.dat', displacement=1e-4):
    """Attach the charges from the fileobj to the Atoms."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    sep = '---------------'
    i = 0 # Counter for the lines
    j = 0 # Index to the Column of charges
    k = 0 # Counter of sep
    assume6columns = False
    for line in fileobj:
        if line[0] == '\n': # check if there is an empty line in the 
            i -= 1          # head of ACF.dat file
        print 'i',i
        if i == 0:
            headings = line
            print headings
            for heading in headings.split():
                if heading == 'BADER' or heading == 'CHARGE':
                    break
                j += 1
            if j == len(headings.split()):
                print 'Can\'t find keyword "BADER" or "CHARGE".' \
                +' Assuming the ACF.dat file has 6 columns.'
                j = 4
                assume6columns = True
        print assume6columns
        if sep in line: # Stop at last seperator line
            if k == 1:
                break
            k += 1
        if not i > 1:
            pass
        else:
            words = line.split()
            if assume6columns is True:
                if len(words) != 6:
                    raise IOError('Number of columns in ACF file incorrect!\n'
                                  'Check that Bader program version >= 0.25')
                
            atom = atoms[int(words[0]) - 1]
            atom.charge = atom.number - float(words[j])
            if displacement is not None: # check if the atom positions match
                xyz = np.array([float(w) for w in words[1:4]]) * Bohr
                assert np.linalg.norm(atom.position - xyz) < displacement
        i += 1

