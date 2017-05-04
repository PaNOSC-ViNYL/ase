import numpy as np
from ase.dft.kpoints import interpolate

eps = np.array([[0, 1, 2], [1, 2, 3]])
path = np.array([[0, 0, 0], [0.25, 0.25, 0]])
bz2ibz = np.array([0, 1, 1, 2])
x = interpolate(path, eps, 42, bz2ibz, [2, 2, 1])
print(x)
