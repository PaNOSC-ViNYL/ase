import numpy as np
from ase.atoms import Atoms
from ase.utils import basestring


def write_crystal(filename, atoms):
    """Method to write atom structure in crystal format
       (fort.34 format)
    """

    # sort
    atoms.set_masses()
    masses = atoms.get_masses()
    indexes = np.argsort(masses)
    atomsnew = Atoms()
    for i in indexes:
        atomsnew = atomsnew + atoms[i]

    if isinstance(filename, basestring):
        myfile = open(filename, 'w')
    else:
        # Assume it's a 'file-like object'
        myfile = filename

    ispbc = atoms.get_pbc()
    box = atoms.get_cell()

    # here it is assumed that the non-periodic direction are z
    # in 2D case, z and y in the 1D case.

    if ispbc[2]:
        myfile.write('%2s %2s %2s %23s \n' % ('3','1','1','E -0.0E+0 DE 0.0E+0( 1)'))
    elif ispbc[1]:
        myfile.write('%2s %2s %2s %23s \n' % ('2','1','1','E -0.0E+0 DE 0.0E+0( 1)'))
        box[2][2] = 500.
    elif ispbc[0]:
        myfile.write('%2s %2s %2s %23s \n' % ('1','1','1','E -0.0E+0 DE 0.0E+0( 1)'))
        box[2][2] = 500.
        box[1][1] = 500.
    else:
        myfile.write('%2s %2s %2s %23s \n' % ('0','1','1','E -0.0E+0 DE 0.0E+0( 1)'))
        box[2][2] = 500.
        box[1][1] = 500.
        box[0][0] = 500.

    # write box
    # crystal dummy
    myfile.write(' %.17E %.17E %.17E \n'
                 % (box[0][0], box[0][1], box[0][2]))
    myfile.write(' %.17E %.17E %.17E \n'
                 % (box[1][0], box[1][1], box[1][2]))
    myfile.write(' %.17E %.17E %.17E \n'
                 % (box[2][0], box[2][1], box[2][2]))
    
    # write symmetry operations (not implemented yet for
    # higher symmetries than C1)
    myfile.write(' %2s \n' % (1))
    myfile.write(' %.17E %.17E %.17E \n' % (1, 0, 0))
    myfile.write(' %.17E %.17E %.17E \n' % (0, 1, 0))
    myfile.write(' %.17E %.17E %.17E \n' % (0, 0, 1))
    myfile.write(' %.17E %.17E %.17E \n' % (0, 0, 0))


    # write coordinates
    myfile.write(' %8s \n' % (len(atomsnew)))
    coords = atomsnew.get_positions()
    atomnum = atomsnew.get_atomic_numbers()
    for iatom, coord in enumerate(coords):
	myfile.write('%5i  %19.16f %19.16f %19.16f \n'
                     % (atomnum[iatom],
                        coords[iatom][0], coords[iatom][1], coords[iatom][2]))

    if isinstance(filename, basestring):
        myfile.close()
