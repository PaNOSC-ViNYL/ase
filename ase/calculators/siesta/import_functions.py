import os
import numpy as np
import string

from ase.units import Bohr
from ase.io.fortranfile import FortranFile


def xv_to_atoms(filename):
    """Create atoms object from xv file.

    Parameters:
        -filename : str. The filename of the '.XV' file.

    return : An Atoms object
    """
    from ase.atoms import Atoms
    if not os.path.exists(filename):
        filename += '.gz'

    with open(filename, 'r') as f:
        # Read cell vectors (lines 1-3)
        vectors = []
        for i in range(3):
            data = string.split(f.readline())
            vectors.append([string.atof(data[j]) * Bohr for j in range(3)])

        # Read number of atoms (line 4)
        string.atoi(string.split(f.readline())[0])

        # Read remaining lines
        speciesnumber, atomnumbers, xyz, V = [], [], [], []
        for line in f.readlines():
            if len(line) > 5:  # Ignore blank lines
                data = string.split(line)
                speciesnumber.append(string.atoi(data[0]))
                atomnumbers.append(string.atoi(data[1]))
                xyz.append([string.atof(data[2 + j]) * Bohr for j in range(3)])
                V.append([string.atof(data[5 + j]) * Bohr for j in range(3)])

    vectors = np.array(vectors)
    atomnumbers = np.array(atomnumbers)
    xyz = np.array(xyz)
    atoms = Atoms(numbers=atomnumbers, positions=xyz, cell=vectors)

    return atoms


def read_rho(fname):
    "Read unformatted Siesta charge density file"

    # TODO:
    #
    # Handle formatted and NetCDF files.
    #
    # Siesta source code (at least 2.0.2) can possibly also
    # save RHO as a _formatted_ file (the source code seems
    # prepared, but there seems to be no fdf-options for it though).
    # Siesta >= 3 has support for saving RHO as a NetCDF file
    # (according to manual)

    fh = FortranFile(fname)

    # Read (but ignore) unit cell vectors
    x = fh.readReals('d')
    if len(x) != 3 * 3:
        raise IOError('Failed to read cell vectors')

    # Read number of grid points and spin components
    x = fh.readInts()
    if len(x) != 4:
        raise IOError('Failed to read grid size')
    gpts = x  # number of 'X', 'Y', 'Z', 'spin' gridpoints

    rho = np.zeros(gpts)
    for ispin in range(gpts[3]):
        for n3 in range(gpts[2]):
            for n2 in range(gpts[1]):
                x = fh.readReals('f')
                if len(x) != gpts[0]:
                    raise IOError('Failed to read RHO[:,%i,%i,%i]' %
                                  (n2, n3, ispin))
                rho[:, n2, n3, ispin] = x

    fh.close()
    return rho
