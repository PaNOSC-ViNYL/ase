"""Module to read and write atoms in PDB file format.

See::

    http://www.wwpdb.org/documentation/file-format

Note: The PDB format saves cell lengths and angles; hence the absolute
orientation is lost when saving.  Saving and loading a file will
conserve the scaled positions, not the absolute ones.
"""

import warnings

import numpy as np

from ase.atoms import Atom, Atoms
from ase.parallel import paropen
from ase.geometry import cellpar_to_cell
from ase.utils import basestring
from ase.io.espresso import label_to_symbol


def read_proteindatabank(fileobj, index=-1, read_arrays=True):
    """Read PDB files."""

    if isinstance(fileobj, basestring):
        fileobj = open(fileobj)

    images = []
    orig = np.identity(3)
    trans = np.zeros(3)
    atoms = Atoms()
    occ = []
    bfactor = []
    for line in fileobj.readlines():
        if line.startswith('CRYST1'):
            cellpar = [float(line[6:15]),  # a
                       float(line[15:24]),  # b
                       float(line[24:33]),  # c
                       float(line[33:40]),  # alpha
                       float(line[40:47]),  # beta
                       float(line[47:54])]  # gamma
            atoms.set_cell(cellpar_to_cell(cellpar))
            atoms.pbc = True
        for c in range(3):
            if line.startswith('ORIGX' + '123'[c]):
                orig[c] = [float(line[10:20]),
                           float(line[20:30]),
                           float(line[30:40])]
                trans[c] = float(line[45:55])

        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                # Atom name is arbitrary and does not necessarily
                # contain the element symbol.  The specification
                # requires the element symbol to be in columns 77+78.
                # Fall back to Atom name for files that do not follow
                # the spec, e.g. packmol.
                try:
                    symbol = label_to_symbol(line[76:78].strip())
                except (KeyError, IndexError):
                    symbol = label_to_symbol(line[12:16].strip())
                # Don't use split() in case there are no spaces
                position = np.array([float(line[30:38]),  # x
                                     float(line[38:46]),  # y
                                     float(line[46:54])])  # z
                try:
                    occ.append(float(line[54:60]))
                    bfactor.append(float(line[60:66]))
                except (IndexError, ValueError):
                    pass
                position = np.dot(orig, position) + trans
                atoms.append(Atom(symbol, position))
            except Exception as ex:
                warnings.warn('Discarding atom when reading PDB file: {}\n{}'
                              .format(line.strip(), ex))
        if line.startswith('END'):
            # End of configuration reached
            # According to the latest PDB file format (v3.30),
            # this line should start with 'ENDMDL' (not 'END'),
            # but in this way PDB trajectories from e.g. CP2K 
            # are supported (also VMD supports this format). 
            if read_arrays and len(occ) == len(atoms):
                atoms.set_array('occupancy', np.array(occ))
            if read_arrays and len(bfactor) == len(atoms):
                atoms.set_array('bfactor', np.array(bfactor))
            images.append(atoms)
            atoms = Atoms()
            occ = []
            bfactor = []
    if len(images) == 0:
        # Single configuration with no 'END' or 'ENDMDL'
        if read_arrays and len(occ) == len(atoms):
            atoms.set_array('occupancy', np.array(occ))
        if read_arrays and len(bfactor) == len(atoms):
            atoms.set_array('bfactor', np.array(bfactor))
        images.append(atoms)
    return images[index]


def write_proteindatabank(fileobj, images, write_arrays=True):
    """Write images to PDB-file."""
    if isinstance(fileobj, basestring):
        fileobj = paropen(fileobj, 'w')

    if hasattr(images, 'get_positions'):
        images = [images]


    rotation = None
    if images[0].get_pbc().any():
        from ase.geometry import cell_to_cellpar, cellpar_to_cell

        currentcell = images[0].get_cell()
        cellpar = cell_to_cellpar(currentcell)
        exportedcell = cellpar_to_cell(cellpar)
        rotation = np.linalg.solve(currentcell, exportedcell)
        # ignoring Z-value, using P1 since we have all atoms defined explicitly
        format = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n'
        fileobj.write(format % (cellpar[0], cellpar[1], cellpar[2],
                                cellpar[3], cellpar[4], cellpar[5]))

    #     1234567 123 6789012345678901   89   67   456789012345678901234567 890
    format = ('ATOM  %5d %4s MOL     1    %8.3f%8.3f%8.3f%6.2f%6.2f'
              '          %2s  \n')

    # RasMol complains if the atom index exceeds 100000. There might
    # be a limit of 5 digit numbers in this field.
    MAXNUM = 100000

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)

    for n, atoms in enumerate(images):
        fileobj.write('MODEL     ' + str(n + 1) + '\n')
        p = atoms.get_positions()
        occupancy = np.ones(len(atoms))
        bfactor = np.zeros(len(atoms))
        if write_arrays:
            if 'occupancy' in atoms.arrays:
                occupancy = atoms.get_array('occupancy')
            if 'bfactor' in atoms.arrays:
                bfactor = atoms.get_array('bfactor')
        if rotation is not None:
            p = p.dot(rotation)
        for a in range(natoms):
            x, y, z = p[a]
            occ = occupancy[a]
            bf = bfactor[a]
            fileobj.write(format % (a % MAXNUM, symbols[a],
                                    x, y, z, occ, bf, symbols[a].upper()))
        fileobj.write('ENDMDL\n')
