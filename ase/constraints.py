from math import sqrt

import numpy as npy


class FixAtoms:
    def __init__(self, indices=None, mask=None):
        """ """
        if indices is None and mask is None:
            raise ValuError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValuError('Use only one of "indices" and "mask".')

        if mask is not None:
            self.index = npy.asarray(mask, bool)
        else:
            self.index = npy.asarray(indices, int)

    def adjust_positions(self, old, new):
        new[self.index] = old[self.index]

    def adjust_forces(self, positions, forces):
        forces[self.index] = 0.0

    def copy(self):
        if self.index.dtype == bool:
            return FixAtoms(mask=self.index.copy())
        else:
            return FixAtoms(indices=self.index.copy())
    

class FixBondLength:
    """Constraint object for fixing a bond length."""
    def __init__(self, a1, a2):
        """Fix distance between atoms with indices a1 and a2."""
        self.indices = [a1, a2]

    def adjust_positions(self, old, new):
        p1, p2 = old[self.indices]
        d = p2 - p1
        p = sqrt(npy.dot(d, d))
        q1, q2 = new[self.indices]
        d = q2 - q1
        q = sqrt(npy.dot(d, d))
        d *= 0.5 * (p - q) / q
        new[self.indices] = (q1 - d, q2 + d)

    def adjust_forces(self, positions, forces):
        d = npy.subtract.reduce(positions[self.indices])
        d2 = npy.dot(d, d)
        d *= 0.5 * npy.dot(npy.subtract.reduce(forces[self.indices]), d) / d2
        forces[self.indices] += (d, -d)

    def copy(self):
        return FixBondLength(*self.indices)
