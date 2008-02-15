from math import sqrt

import numpy as npy


class FixAtoms:
    """Constraint object for fixing some chosen atoms."""
    def __init__(self, indices=None, mask=None):
        """Constrain chosen atoms.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should be constrained.
        mask : list of bool
           One boolean per atom indicating if the atom should be
           constrained or not.
           
        Examples
        --------
        Fix all Copper atoms:

        >>> c = FixAtoms(mask=[s == 'Cu' for s in atoms.get_chemical_symbols()])
        >>> atoms.set_constraint(c)

        Fix all atoms with z-coordinate less than 1.0 Angstrom:

        >>> c = FixAtoms(mask=atoms.positions[:, 2] < 1.0)
        >>> atoms.set_constraint(c)
        """

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
    
    def __repr__(self):
        if self.index.dtype == bool:
            return 'FixAtoms(mask=%s)' % ints2string(self.index.astype(int))
        return 'FixAtoms(indices=%s)' % ints2string(self.index)

def ints2string(x, threshold=10):
    """Convert ndarray of ints to string."""
    if len(x) <= threshold:
        return str(x.tolist())
    return str(x[:threshold].tolist())[:-1] + ', ...]'


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

    def __repr__(self):
        return 'FixBondLength(%d, %d)' % tuple(self.indices)
