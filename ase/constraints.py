from math import sqrt

import numpy as npy


class FixAtoms:
    def __init__(self, indices):
        self.fixed = npy.array(indices)

    def adjust_positions(self, old, new):
        new[self.fixed] = old[self.fixed]

    def adjust_forces(self, positions, forces):
        forces[self.fixed] = 0.0

    def copy(self):
        return FixAtoms(self.fixed.copy())

    def delete_atoms(self, indices):
        if self.fixed.dtype == bool:
            pass
    

class FixBondLength:
    def __init__(self, a1, a2):
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
