from numpy import zeros

from ase.io.fortranfile import FortranFile


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
    
    rho = zeros(gpts)
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

