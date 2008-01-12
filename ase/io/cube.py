import numpy as npy

from ase.atoms import Atoms
from ase.units import Bohr


def write_cube(fileobj, atoms, data=None):
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'w')
        
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]

    fileobj.write('ase cube file\n')
    fileobj.write('OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n')

    natoms = len(atoms)
    fileobj.write('%5d%12.6f%12.6f%12.6f\n' % (natoms, 0.0, 0.0, 0.0))

    if data is None:
        data = [[[1.0]]]
    data = npy.asarray(data, float)
    
    cell = atoms.get_cell()
    for i in range(3):
        n = data.shape[i]
        d = cell[i] / n / Bohr
        fileobj.write('%5d%12.6f%12.6f%12.6f\n' % (n, d[0], d[1], d[2]))

    positions = atoms.get_positions() / Bohr
    numbers = atoms.get_atomic_numbers()
    for Z, (x, y, z) in zip(numbers, positions):
        fileobj.write('%5d%12.6f%12.6f%12.6f%12.6f\n' % (Z, 0.0, x, y, z)) 

    for dyz in data:
        for dz in dyz:
            for i, d in enumerate(dz):
                fileobj.write('%e ' % d)
                if i % 6 == 5:
                    fileobj.write('\n')
                fileobj.write('\n')


def read_cube(fileobj, index=-1, read_data=False):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    readline = fileobj.readline
    readline()
    axes = ['XYZ'.index(s[0]) for s in readline().split()[2::3]]
    natoms = int(readline().split()[0])

    cell = npy.empty((3, 3))
    shape = []
    for i in range(3):
        n, x, y, z = [float(s) for s in readline()]
        cell[i] = n * Bohr * npy.array([x, y, z])
        shape.append(n)
        
    numbers = npy.empty(natoms, int)
    positions = npy.empty((natoms, 3))
    for i in range(natoms):
        line = readline().split()
        numbers[i] = int(line[0])
        positions[i] = [float(s) for s in line[2:]]

    positions *= Bohr
    atoms = Atoms(numbers=numbers, positions=positions, cell=cell)

    if read_data:
        data = npy.array([float(s)
                          for s in fileobj.read().split()]).reshape(shape)
        if axes != [0, 1, 2]:
            data = data.transpose(axes).copy()
        return data, atoms

    return atoms


def read_cube_data(fileobj):
    return read_cube(fileobj, read_data=True)
