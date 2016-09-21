from ase.dft.kpoints import bandpath
import numpy as np
print(bandpath('GX', np.eye(3), 3))

