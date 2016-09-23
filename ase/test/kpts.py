from ase.dft.kpoints import bandpath
import numpy as np
print(bandpath('GX', np.eye(3), 3))
kpts = [[0,0,0], [0.25, .25, 0], [.49, .49, 0], [.5, .5, 0], [.51, .51, 0],
        [.5, .5, .25], [.5, .5, .5],
        [.25, .5, .25], [0, .5, 0]]
print(xaxis_from_kpts(kpts, cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                      crystal_structure='cubic'))
