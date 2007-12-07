import time
from math import sqrt

import numpy as npy

from ase.utils import rotate
from ase.data import cpk_colors, covalent_radii


class EPS:
    scale0 = 20.0
    def __init__(self, atoms,
                 rotation='', show_unit_cell=False, radii=None,
                 bbox=None):
        self.numbers = atoms.get_atomic_numbers()
        self.colors = cpk_colors[self.numbers]

        if radii == None:
            radii = covalent_radii[self.numbers]
            
        natoms = len(atoms)

        if isinstance(rotation, str):
            rotation = rotate(rotation)

        A = atoms.get_cell()
        if show_unit_cell:
            L, T, D = self.cell_to_lines(A)
            C = npy.empty((2, 2, 2, 3))
            for c1 in range(2):
                for c2 in range(2):
                    for c3 in range(2):
                        C[c1, c2, c3] = npy.dot([c1, c2, c3], A)
            C.shape = (8, 3)
            C = npy.dot(C, rotation)
        else:
            L = npy.empty((0, 3))
            T = None
            D = None
            C = None

        nlines = len(L)

        X = npy.empty((natoms + nlines, 3))
        R = atoms.get_positions()
        X[:natoms] = R
        X[natoms:] = L

        r2 = radii**2
        for n in range(nlines):
            d = D[T[n]]
            if ((((R - L[n] - d)**2).sum(1) < r2) &
                (((R - L[n] + d)**2).sum(1) < r2)).any():
                T[n] = -1

        X = npy.dot(X, rotation)
        R = X[:natoms]

        scale = self.scale0
        if bbox is None:
            X1 = (R - radii[:, None]).min(0) 
            X2 = (R + radii[:, None]).max(0) 
            print X1, X2
            if show_unit_cell == 2:
                X1 = npy.minimum(X1, C.min(0))
                X2 = npy.maximum(X2, C.max(0))
                print X1, X2
            M = (X1 + X2) / 2
            S = 1.05 * (X2 - X1)
            w = scale * S[0]
            if w > 500:
                w = 500
                scale = w / S[0]
            h = scale * S[1]
            offset = npy.array([scale * M[0] - w / 2, scale * M[1] - h / 2, 0])
        else:
            scale = 50.0
            w = (bbox[2] - bbox[0]) * scale
            h = (bbox[3] - bbox[1]) * scale
            offset = npy.array([bbox[0], bbox[1], 0]) * scale

        self.w = w
        self.h = h
        
        X *= scale
        X -= offset
        X[:, 1] = h - X[:, 1]

        if nlines > 0:
            D = npy.dot(D, rotation)[:, :2] * scale
            D[:, 1] = -D[:, 1]
        
        if C is not None:
            C *= scale
            C -= offset
            C[:, 1] = h - C[:, 1]

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

        X = npy.empty((nlines, 3))
        T = npy.empty(nlines, int)
        D = npy.zeros((3, 3))

        n1 = 0
        for c in range(3):
            n = nn[c]
            dd = A[c] / (4 * n - 2)
            D[c] = dd
            P = npy.arange(1, 4 * n + 1, 4)[:, None] * dd
            T[n1:] = c
            for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                n2 = n1 + n
                X[n1:n2] = P + i * A[(c + 1) % 3] + j * A[(c + 2) % 3]
                n1 = n2

        return X, T, D

    def write(self, filename):
        self.filename = filename
        self.write_header()
        self.write_body()
        self.write_trailer()

    def write_header(self):
        from matplotlib.backends.backend_ps import RendererPS, \
             GraphicsContextPS, psDefs

        self.fd = open(self.filename, 'w')
        self.fd.write('%!PS-Adobe-3.0 EPSF-3.0\n')
        self.fd.write('%%Creator: G2\n')
        self.fd.write('%%CreationDate: %s\n' % time.ctime(time.time()))
        self.fd.write('%%Orientation: portrait\n')
        bbox = (0, 0, self.w, self.h)
        self.fd.write('%%%%BoundingBox: %d %d %d %d\n' % bbox)
        self.fd.write('%%EndComments\n')

        Ndict = len(psDefs)
        self.fd.write('%%BeginProlog\n')
        self.fd.write('/mpldict %d dict def\n' % Ndict)
        self.fd.write('mpldict begin\n')
        for d in psDefs:
            d = d.strip()
            for l in d.split('\n'):
                self.fd.write(l.strip() + '\n')
        self.fd.write('%%EndProlog\n')

        self.fd.write('mpldict begin\n')
        self.fd.write('%d %d 0 0 clipbox\n' % (self.w, self.h))

        self.renderer = RendererPS(self.w, self.h, self.fd)
        self.gc = GraphicsContextPS()
        self.line = self.renderer.draw_line
        
    def write_body(self):
        # Write the figure
        line = self.line
        arc = self.renderer.draw_arc
        indices = self.X[:, 2].argsort()
        for a in indices:
            x, y = self.X[a, :2]
            if a < self.natoms:
                da = self.d[a]
                arc(self.gc, tuple(self.colors[a]), x, y, da, da, 0, 360, 0)
            else:
                a -= self.natoms
                c = self.T[a]
                if c != -1:
                    hx, hy = self.D[c]
                    line(self.gc, x - hx, y - hy, x + hx, y + hy)

    def write_trailer(self):
        self.fd.write('end\n')
        self.fd.write('showpage\n')
        self.fd.close()


def write_eps(filename, atoms, **parameters):
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]
    EPS(atoms, **parameters).write(filename)
