import numpy as np
from ase import Atoms
from ase.geometry import crystal_structure_from_cell
from ase.dft.kpoints import get_special_points

for cell in [[[1, 0, 0], [0.5, 3**0.5 / 2, 0], [0, 0, 1]],
             [[1, 0, 0], [-0.5, 3**0.5 / 2, 0], [0, 0, 1]],
             [[0.5, -3**0.5 / 2, 0], [0.5, 3**0.5 / 2, 0], [0, 0, 1]]]:
    a = Atoms(cell=cell, pbc=True)
    print(crystal_structure_from_cell(a.cell))
    r = a.get_reciprocal_cell()
    k = get_special_points(a.cell)['K']
    print(np.dot(k, r))