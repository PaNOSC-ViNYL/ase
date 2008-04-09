import numpy as npy

from ase.atoms import Atoms

def write_plt(filename, atoms, data):
    if isinstance(atoms, Atoms):
        cell = atoms.get_cell()
    else:
        cell = npy.asarray(atoms, float)

    if cell.ndim == 2:
        c = cell.copy()
        cell = c.diagonal()
        c.flat[::4] = 0.0
        if c.any():
            raise ValueError('Unit cell must be orthorhombic!')

    f = open(filename, 'w')
    npy.array([3, 4], npy.int32).tofile(f)

    dims = npy.array(data.shape, npy.int32)
    dims[::-1].tofile(f)

    for n, L in zip(dims[::-1], cell[::-1]):
        if n % 2 == 0:
            d = L / n
            npy.array([0.0, L - d], npy.float32).tofile(f)
        else:
            d = L / (n + 1)
            npy.array([d, L - d], npy.float32).tofile(f)

    data.astype(npy.float32).T.tofile(f)
    f.close()
