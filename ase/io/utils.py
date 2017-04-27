import numpy as np
from math import sqrt
from ase.utils import rotate
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.utils import basestring


def generate_writer_variables(
        writer, atoms,
        rotation='', show_unit_cell=False, radii=None,
        bbox=None, colors=None, scale=20):
    writer.numbers = atoms.get_atomic_numbers()
    writer.colors = colors
    if colors is None:
        writer.colors = jmol_colors[writer.numbers]

    if radii is None:
        radii = covalent_radii[writer.numbers]
    elif isinstance(radii, float):
        radii = covalent_radii[writer.numbers] * radii
    else:
        radii = np.array(radii)

    natoms = len(atoms)

    if isinstance(rotation, basestring):
        rotation = rotate(rotation)

    A = atoms.get_cell()
    if show_unit_cell > 0:
        L, T, D = cell_to_lines(writer, A)
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

    writer.w = w
    writer.h = h

    X *= scale
    X -= offset

    if nlines > 0:
        D = np.dot(D, rotation)[:, :2] * scale

    if C is not None:
        C *= scale
        C -= offset

    A = np.dot(A, rotation)
    A *= scale

    writer.A = A
    writer.X = X
    writer.D = D
    writer.T = T
    writer.C = C
    writer.natoms = natoms
    writer.d = 2 * scale * radii


def cell_to_lines(writer, A):
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
            X[n1:n2] = P + i * A[c - 2] + j * A[c - 1]
            n1 = n2

    return X, T, D


def make_patch_list(writer):
    try:
        from matplotlib.path import Path
    except ImportError:
        Path = None
        from matplotlib.patches import Circle, Polygon
    else:
        from matplotlib.patches import Circle, PathPatch

    indices = writer.X[:, 2].argsort()
    patch_list = []
    for a in indices:
        xy = writer.X[a, :2]
        if a < writer.natoms:
            r = writer.d[a] / 2
            if ((xy[1] + r > 0) and (xy[1] - r < writer.h) and
                (xy[0] + r > 0) and (xy[0] - r < writer.w)):
                patch = Circle(xy, r, facecolor=writer.colors[a],
                                edgecolor='black')
        else:
            a -= writer.natoms
            c = writer.T[a]
            if c != -1:
                hxy = writer.D[c]
                if Path is None:
                    patch = Polygon((xy + hxy, xy - hxy))
                else:
                    patch = PathPatch(Path((xy + hxy, xy - hxy)))
        patch_list.append(patch)
    return patch_list
