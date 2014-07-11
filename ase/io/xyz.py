import re
import numpy as np

from ase.atoms import Atoms
from ase.parallel import paropen


def read_xyz(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    if not isinstance(index, int) and not isinstance(index, slice):
        raise TypeError('Index argument is neither slice nor integer!')

    fileobj.seek(0)
    natoms = int(fileobj.readline())

    for i, l in enumerate(fileobj):
        pass

    ln = i + 2
    lnsnp = natoms + 2
    lastsnap = ln // lnsnp

    rvrs = False

    if isinstance(index, int):
        if index < 0:
            tmpsnp = lastsnap + index
            trbl = range(tmpsnp, tmpsnp + 1, 1)
        else:
            trbl = range(index, index + 1, 1)
        rtnndx = -1
    elif isinstance(index, slice):
        start = index.start
        stop = index.stop
        step = index.step

        if start is None:
            start = 0
        elif start < 0:
            start = lastsnap + start

        if step is None:
            step = 1

        if stop is None:
            stop = lastsnap
        elif stop < 0:
            stop = lastsnap + stop

        trbl = range(start, stop, step)

        if step < 0:
            rvrs = True
            trbl.reverse()

        rtnndx = slice(len(trbl))

    images = []
    current = 0

    fileobj.seek(0)

    for index in trbl:
        for lnndx in range(current, index * lnsnp):
            line = fileobj.readline()
        
        line = fileobj.readline()
        line = fileobj.readline()

        # XMOL extended xyz format
        # see http://www.pierrehirel.info/codes/atomsk/en/format_xyz.html
        cell = None
        coordstr = re.findall(r'Lattice="(.*?)"', line)
        if len(coordstr):
            cell = np.array(
                [float(word) for word in coordstr[0].split()]).reshape((3,3))

        symbols = []
        positions = []
        for ln in range(natoms):
            line = fileobj.readline()

            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])

        current = (index + 1) * lnsnp

        images.append(Atoms(symbols=symbols, positions=positions, cell=cell))

    if rvrs:
        images.reverse()

    return images[rtnndx]


def write_xyz(fileobj, images, comment=None):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    for atoms in images:
        if comment is None:
            # XMOL extended xyz format
            # see http://www.pierrehirel.info/codes/atomsk/en/format_xyz.html
            comm = ('Lattice="' + ' '.join([
                str(c) for c in atoms.get_cell().flat]) + '"')
        else:
            comm = comment
        fileobj.write('%d\n%s\n' % (natoms, comm))
        for s, (x, y, z) in zip(symbols, atoms.get_positions()):
            fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))
