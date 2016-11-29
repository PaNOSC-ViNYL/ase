"""
IO functions for DMol3 file formats.


The formats follow some very strict formatting

Car
----
col: 1-5     atom name
col: 7-20    x Cartesian coordinate of atom  in A
col: 22-35   y Cartesian coordinate of atom  in A
col: 37-50   z Cartesian coordinate of atom  in A
col: 52-55   type of residue containing atom
col: 57-63   residue sequence name   relative to beginning of current molecule,
                left justified
col: 64-70   potential type of atom  left justified
col: 72-73   element symbol
col: 75-80   partial charge on atom

Incoor
------
$cell vectors
             37.83609647462165    0.00000000000000    0.00000000000000
              0.00000000000000   37.60366016124745    0.00000000000000
              0.00000000000000    0.00000000000000   25.29020473078921
$coordinates
Si           15.94182672614820    1.85274838936809   16.01426481346124
Si            4.45559370448989    2.68957177851318   -0.05326937257442
$end

"""

from __future__ import print_function
from datetime import datetime
import numpy as np

from ase import Atom, Atoms
from ase.geometry.cell import cell_to_cellpar, cellpar_to_cell
from ase.units import Bohr


def write_dmol_car(filename, atoms):
    """ Write a dmol car-file from an Atoms object

    Notes
    -----
    Can not handle multiple images.
    Only allows for pbc 111 or 000.
    """

    f = open(filename, 'w')
    f.write('!BIOSYM archive 3\n')

    dt = datetime.now()
    if np.all(atoms.pbc):
        a, b, c, alpha, beta, gamma = cell_to_cellpar(atoms.get_cell())
        f.write('PBC=ON\n\n')
        f.write('!DATE     %s\n' % dt.strftime('%b %d %H:%m:%S %Y'))
        f.write('PBC %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n' %
                (a, b, c, alpha, beta, gamma))
    elif not np.any(atoms.pbc):  # [False,False,False]
        f.write('PBC=OFF\n\n')
        f.write('!DATE     %s\n' % dt.strftime('%b %d %H:%m:%S %Y'))
    else:
        raise ValueError('PBC must be all true or all false for .car format')
    for i, a in enumerate(atoms):
        f.write('%-6s  %12.8f   %12.8f   %12.8f XXXX 1      xx      %-2s  '
                '0.000\n' % (a.symbol + str(i+1), a.x, a.y, a.z, a.symbol))
    f.write('end\nend\n')
    f.close()


def read_dmol_car(filename):
    """ Read a dmol car-file and return an Atoms object.

    Notes
    -----
    Cell is constructed from cellpar so orientation of cell might be off.
    """

    try:
        lines = open(filename, 'r').readlines()
        atoms = Atoms()

        start_line = 4

        if lines[1][4:6] == 'ON':
            start_line += 1
            cell_dat = np.array([float(fld) for fld in lines[4].split()[1:7]])
            cell = cellpar_to_cell(cell_dat)
            atoms.cell = cell
            atoms.pbc = [True, True, True]

        for line in lines[start_line:]:
            if line.startswith('end'):
                break
            flds = line.split()
            atoms.append(Atom(flds[7], flds[1:4]))
        return atoms
    except IOError:
        print('Could not read car file %s' % filename)


def write_dmol_incoor(filename, atoms, bohr=True):
    """ Write a dmol incoor-file from an Atoms object

    Notes
    -----
    Only used for pbc 111.
    Can not handle multiple images.
    DMol3 expect data in .incoor files to be in bohr, if bohr is false however
    the data is written in Angstroms.
    """

    if not np.all(atoms.pbc):
        raise ValueError('PBC must be all true for .incoor format')

    if bohr:
        cell = atoms.cell / Bohr
        positions = atoms.positions / Bohr
    else:
        cell = atoms.cell
        positions = atoms.positions

    f = open(filename, 'w')
    f.write('$cell vectors\n')
    f.write('            %18.14f  %18.14f  %18.14f\n' % (
        cell[0, 0], cell[0, 1], cell[0, 2]))
    f.write('            %18.14f  %18.14f  %18.14f\n' % (
        cell[1, 0], cell[1, 1], cell[1, 2]))
    f.write('            %18.14f  %18.14f  %18.14f\n' % (
        cell[2, 0], cell[2, 1], cell[2, 2]))

    f.write('$coordinates\n')
    for a, pos in zip(atoms, positions):
        f.write('%-12s%18.14f  %18.14f  %18.14f \n' % (
            a.symbol, pos[0], pos[1], pos[2]))
    f.write('$end\n')
    f.close()


def read_dmol_incoor(filename, bohr=True):
    """ Reads an incoor file and returns an atoms object.

    Notes
    -----
    If bohr is True then incoor is assumed to be in bohr and the data
    is rescaled to Angstrom.
    """

    try:
        lines = open(filename, 'r').readlines()
        atoms = Atoms()
        for i, line in enumerate(lines):
            if line.startswith('$cell vectors'):
                cell = np.zeros((3, 3))
                for j, line in enumerate(lines[i + 1:i + 4]):
                    cell[j, :] = [float(fld) for fld in line.split()]
                atoms.cell = cell
            if line.startswith('$coordinates'):
                j = i + 1
                while True:
                    if lines[j].startswith('$end'):
                        break
                    flds = lines[j].split()
                    atoms.append(Atom(flds[0], flds[1:4]))
                    j += 1
    except IOError:
        print('Could not read incoor file %s' % filename)
    if bohr:
        atoms.cell = atoms.cell * Bohr
        atoms.positions = atoms.positions * Bohr
    return atoms


def write_dmol_arc(filename, images):
    """ Writes all images to file filename in arc format.

    Similar to the .car format only pbc 111 or 000 is supported.
    """

    f = open(filename, 'w')
    f.write('!BIOSYM archive 3\n')
    if np.all(images[0].pbc):
        f.write('PBC=ON\n\n')
    elif not np.any(images[0].pbc):
        f.write('PBC=OFF\n\n')
    else:
        raise ValueError('PBC must be all true or all false for .arc format')
    for atoms in images:
        dt = datetime.now()
        if np.all(atoms.pbc):
            a, b, c, alpha, beta, gamma = cell_to_cellpar(atoms.get_cell())
            f.write('!DATE     %s\n' % dt.strftime('%b %d %H:%m:%S %Y'))
            f.write('PBC %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n' %
                    (a, b, c, alpha, beta, gamma))
        elif not np.any(atoms.pbc):  # [False,False,False]
            f.write('!DATE     %s\n' % dt.strftime('%b %d %H:%m:%S %Y'))
        else:
            raise ValueError(
                'PBC must be all true or all false for .arc format')
        for i, a in enumerate(atoms):
            f.write('%-6s  %12.8f   %12.8f   %12.8f XXXX 1      xx      %-2s  '
                    '0.000\n' % (a.symbol + str(i+1), a.x, a.y, a.z, a.symbol))
        f.write('end\nend\n')
        f.write('\n')
    f.close()


def read_dmol_arc(filename, index=-1):
    """ Read a dmol arc-file and return a series of Atoms objects (images). """

    try:
        lines = open(filename, 'r').readlines()
        images = []

        if lines[1].startswith('PBC=ON'):
            pbc = True
        elif lines[1].startswith('PBC=OFF'):
            pbc = False
        else:
            raise RuntimeError('Could not read pbc from second line in %s'
                               % filename)
        i = 0
        while i < len(lines):
            image = Atoms()
            # parse single image
            if lines[i].startswith('!DATE'):
                # read cell
                if pbc:
                    cell_dat = np.array([float(fld)
                                         for fld in lines[i + 1].split()[1:7]])
                    image.cell = cellpar_to_cell(cell_dat)
                    image.pbc = [True, True, True]
                    i += 1
                i += 1
                # read atoms
                while not lines[i].startswith('end'):
                    flds = lines[i].split()
                    image.append(Atom(flds[7], flds[1:4]))
                    i += 1
                images.append(image)
            if len(images) == index:
                return images[-1]
            i += 1
    except IOError:
        print('Could not read arc file %s' % filename)

    # return requested images, code borrowed from ase/io/trajectory.py
    if isinstance(index, int):
        return images[index]
    else:
        step = index.step or 1
        if step > 0:
            start = index.start or 0
            if start < 0:
                start += len(images)
            stop = index.stop or len(images)
            if stop < 0:
                stop += len(images)
        else:
            if index.start is None:
                start = len(images) - 1
            else:
                start = index.start
                if start < 0:
                    start += len(images)
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop
                if stop < 0:
                    stop += len(images)
        return [images[j] for j in range(start, stop, step)]
