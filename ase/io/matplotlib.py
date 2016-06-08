import time
from math import sqrt
from distutils.version import LooseVersion

import numpy as np

from ase.utils import rotate
from ase.data import covalent_radii
from ase.data.colors import jmol_colors


class MATPLOTLIB:
    def __init__(self, atoms, ax,
                 rotation='', show_unit_cell=False, radii=None,
                 bbox=None, colors=None, scale=20):
#        if ax is None:
#            self.fig, self.ax = plt.subplots()
#        else:
        self.ax = ax
        self.figure = ax.figure
        self.numbers = atoms.get_atomic_numbers()
        self.colors = colors
        if colors is None:
            self.colors = jmol_colors[self.numbers]

        if radii is None:
            radii = covalent_radii[self.numbers]
        elif isinstance(radii, float):
            radii = covalent_radii[self.numbers] * radii
        else:
            radii = np.array(radii)

        natoms = len(atoms)

        if isinstance(rotation, str):
            rotation = rotate(rotation)

        A = atoms.get_cell()
        if show_unit_cell > 0:
            L, T, D = self.cell_to_lines(A)
            C = np.empty((2, 2, 2, 3))
            for c1 in range(2):
                for c2 in range(2):
                    for c3 in range(2):
                        C[c1, c2, c3] = np.dot([c1, c2, c3], A)
            C.shape = (8, 3)
            C = np.dot(C, rotation)  # Unit cell vertices
        else:
            L = np.empty((0, 3))
            T = None
            D = None
            C = None

        nlines = len(L)

        X = np.empty((natoms + nlines, 3))
        R = atoms.get_positions()
        X[:natoms] = R
        X[natoms:] = L

        r2 = radii**2
        for n in range(nlines):
            d = D[T[n]]
            if ((((R - L[n] - d)**2).sum(1) < r2) &
                (((R - L[n] + d)**2).sum(1) < r2)).any():
                T[n] = -1

        X = np.dot(X, rotation)
        R = X[:natoms]

        if bbox is None:
            X1 = (R - radii[:, None]).min(0)
            X2 = (R + radii[:, None]).max(0)
            if show_unit_cell == 2:
                X1 = np.minimum(X1, C.min(0))
                X2 = np.maximum(X2, C.max(0))
            M = (X1 + X2) / 2
            S = 1.05 * (X2 - X1)
            w = scale * S[0]
            if w > 500:
                w = 500
                scale = w / S[0]
            h = scale * S[1]
            offset = np.array([scale * M[0] - w / 2, scale * M[1] - h / 2, 0])
        else:
            w = (bbox[2] - bbox[0]) * scale
            h = (bbox[3] - bbox[1]) * scale
            offset = np.array([bbox[0], bbox[1], 0]) * scale

        X *= scale
        X -= offset

        self.xlim = 0, w
        self.ylim = 0, h 

        if nlines > 0:
            D = np.dot(D, rotation)[:, :2] * scale

        if C is not None:
            C *= scale
            C -= offset

        A = np.dot(A, rotation)
        A *= scale

        self.A = A
        self.X = X
        self.D = D
        self.T = T
        self.C = C
        self.natoms = natoms
        self.d = 2 * scale * radii

    def cell_to_lines(self, A):
        nlines = 0
        nn = []
        for c in range(3):
            d = sqrt((A[c]**2).sum())
            n = max(2, int(d / 0.3))
            nn.append(n)
            nlines += 4 * n

        X = np.empty((nlines, 3))
        T = np.empty(nlines, int)
        D = np.zeros((3, 3))

        n1 = 0
        for c in range(3):
            n = nn[c]
            dd = A[c] / (4 * n - 2)
            D[c] = dd
            P = np.arange(1, 4 * n + 1, 4)[:, None] * dd
            T[n1:] = c
            for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                n2 = n1 + n
                X[n1:n2] = P + i * A[(c + 1) % 3] + j * A[(c + 2) % 3]
                n1 = n2

        return X, T, D

    def write(self):
        self.write_body()
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

    def write_body(self):
        try:
            from matplotlib.path import Path
        except ImportError:
            Path = None
            from matplotlib.patches import Circle, Polygon
        else:
            from matplotlib.patches import Circle, PathPatch

        indices = self.X[:, 2].argsort()
        for a in indices:
            xy = self.X[a, :2]
            if a < self.natoms:
                r = self.d[a] / 2
                circle = Circle(xy, r, facecolor=self.colors[a])
                self.ax.add_patch(circle)
            else:
                a -= self.natoms
                c = self.T[a]
                if c != -1:
                    hxy = self.D[c]
                    if Path is None:
                        line = Polygon((xy + hxy, xy - hxy))
                    else:
                        line = PathPatch(Path((xy + hxy, xy - hxy)))
                    self.ax.add_patch(line)

def write_matplotlib(atoms, ax, **parameters):
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]
    MATPLOTLIB(atoms, ax, **parameters).write()
