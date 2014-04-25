import numpy as np

from ase.lattice.surface import fcc211

atoms = fcc211('Au', (3, 5, 8), vacuum=10.)
assert len(atoms) == 120

atoms = atoms.repeat((2, 1, 1))
assert np.allclose(atoms.get_distance(0, 130), 2.88499566724)
