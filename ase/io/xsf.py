import numpy as np

from ase.atoms import Atoms
from ase.units import Hartree
from ase.parallel import paropen
from ase.data import atomic_numbers
from ase.calculators.singlepoint import SinglePointCalculator


def write_xsf(fileobj, images, data=None):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')
        
    if not isinstance(images, (list, tuple)):
        images = [images]

    fileobj.write('ANIMSTEPS %d\n' % len(images))

    numbers = images[0].get_atomic_numbers()
    
    pbc = images[0].get_pbc()
    npbc = sum(pbc)
    if pbc[2]:
        fileobj.write('CRYSTAL\n')
        assert npbc == 3
    elif pbc[1]:
        fileobj.write('SLAB\n')
        assert npbc == 2
    elif pbc[0]:
        fileobj.write('POLYMER\n')
        assert npbc == 1
    else:
        fileobj.write('ATOMS\n')
        assert npbc == 0

    for n, atoms in enumerate(images):
        if pbc.any():
            fileobj.write('PRIMVEC %d\n' % (n + 1))
            cell = atoms.get_cell()
            for i in range(3):
                fileobj.write(' %.14f %.14f %.14f\n' % tuple(cell[i]))

            fileobj.write('PRIMCOORD %d\n' % (n + 1))

        # Get the forces if it's not too expensive:
        calc = atoms.get_calculator()
        if (calc is not None and
            (hasattr(calc, 'calculation_required') and
             not calc.calculation_required(atoms,
                                           ['energy', 'forces', 'stress']))):
            forces = atoms.get_forces()
        else:
            forces = None

        pos = atoms.get_positions()

        if pbc.any():
            fileobj.write(' %d 1\n' % len(pos))
        for a in range(len(pos)):
            fileobj.write(' %2d' % numbers[a])
            fileobj.write(' %20.14f %20.14f %20.14f' % tuple(pos[a]))
            if forces is None:
                fileobj.write('\n')
            else:
                fileobj.write(' %20.14f %20.14f %20.14f\n' % tuple(forces[a]))
            
    if data is None:
        return

    fileobj.write('BEGIN_BLOCK_DATAGRID_3D\n')
    fileobj.write(' data\n')
    fileobj.write(' BEGIN_DATAGRID_3Dgrid#1\n')

    data = np.asarray(data)
    if data.dtype == complex:
        data = np.abs(data)

    shape = data.shape
    fileobj.write('  %d %d %d\n' % shape)

    cell = atoms.get_cell()
    origin = np.zeros(3)
    for i in range(3):
        if not pbc[i]:
            origin += cell[i] / shape[i]
    fileobj.write('  %f %f %f\n' % tuple(origin))

    for i in range(3):
        fileobj.write('  %f %f %f\n' %
                      tuple(cell[i] * (shape[i] + 1) / shape[i]))

    for x in range(shape[2]):
        for y in range(shape[1]):
            fileobj.write('   ')
            fileobj.write(' '.join(['%f' % d for d in data[x, y]]))
            fileobj.write('\n')
        fileobj.write('\n')

    fileobj.write(' END_DATAGRID_3D\n')
    fileobj.write('END_BLOCK_DATAGRID_3D\n')


def read_xsf(fileobj, index=-1, read_data=False):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    readline = fileobj.readline
    while True:
        line = readline()
        if line[0] != '#':
            line = line.strip()
            break

    if 'ANIMSTEPS' in line:
        nimages = int(line.split()[1])
        line = readline().strip()
    else:
        nimages = 1

    if 'CRYSTAL' in line:
        pbc = (True, True, True)
    elif 'SLAB' in line:
        pbc = (True, True, False)
    elif 'POLYMER' in line:
        pbc = (True, False, False)
    else:
        assert 'ATOMS' in line
        pbc = (False, False, False)

    images = []
    for n in range(nimages):
        cell = None
        if any(pbc):
            line = readline().strip()
            assert 'PRIMVEC' in line, line
            cell = []
            for i in range(3):
                cell.append([float(x) for x in readline().split()])

        line = readline().strip()
        if line[0] == 'CONVVEC':
            for i in range(3):
                readline()
            line = readline().strip()

        if any(pbc):
            assert 'PRIMCOORD' in line
            natoms = int(readline().split()[0])
            lines = [readline() for _ in range(natoms)]
        else:
            lines = []
            while line != '' and not line.startswith('BEGIN'):
                lines.append(line)
                line = readline()
                if line.startswith('BEGIN'):
                    # We read "too far" and accidentally got the header
                    # of the data section.  This happens only when parsing
                    # ATOMS blocks, because one cannot infer their length.
                    # We will remember the line until later then.
                    data_header_line = line

        numbers = []
        positions = []
        for line in lines:
            tokens = line.split()
            symbol = tokens[0]
            if symbol.isdigit():
                numbers.append(int(symbol))
            else:
                numbers.append(atomic_numbers[symbol])
            positions.append([float(x) for x in tokens[1:]])

        positions = np.array(positions)
        if len(positions[0]) == 3:
            forces = None
        else:
            forces = positions[:, 3:] * Hartree
            positions = positions[:, :3]

        image = Atoms(numbers, positions, cell=cell, pbc=pbc)

        if forces is not None:
            image.set_calculator(SinglePointCalculator(image, forces=forces))
        images.append(image)

    if read_data:
        if any(pbc):
            line = readline()
        else:
            line = data_header_line
        assert 'BEGIN_BLOCK_DATAGRID_3D' in line, line
        readline()  # name
        line = readline()
        assert 'BEGIN_DATAGRID_3D' in line, line

        shape = [int(x) for x in readline().split()]
        readline()  # start

        for i in range(3):
            readline()
            
        n_data = shape[0] * shape[1] * shape[2]
        data = np.array([float(readline())
                         for s in range(n_data)]).reshape(shape[::-1])
        data = np.swapaxes(data, 0, 2)
        
        return data, images[index]

    return images[index]
