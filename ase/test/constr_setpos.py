import numpy as np

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    """Replacement for old numpy.testing.utils.array_almost_equal."""
    return (np.abs(a1 - a2) < tol).all()

from ase.structure import molecule
from ase.constraints import FixAtoms

m = molecule('H2')
c = FixAtoms(indices=[atom.index for atom in m])
m.set_constraint(c)

pos1 = m.get_positions()
pos = m.get_positions()
# shift z-coordinates by 1.
pos[:,2] += 1.
for a in range(len(m)):
    assert abs(pos[a, 2] - 1. - pos1[a, 2]) < 1.0e-6
m.set_positions(pos)
# note that set_positions fails silently to set the new positions
# due to the presence of constraints!
assert array_almost_equal(pos1, m.get_positions())
