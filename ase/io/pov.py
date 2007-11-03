import numpy as npy

from ase.io.eps import EPS
from ase.data import chemical_symbols, covalent_radii


class POVRAY(EPS):
    scale0 = 1.0
    def cell_to_lines(self, A):
        return npy.empty((0, 3)), None, None

    def write(self, filename):
        if filename.endswith('.pov'):
            ini = filename[:-4] + '.ini'
        else:
            ini = filename + '.ini'

        w = 50 * self.w
        if w > 500:
            w = 500
            h = 500 / self.w * self.h
        else:
            h = 50 * self.h
        
        open(ini, 'w').write('Input_File_Name=%s\n+W%d +H%d\n' %
                             (filename, w, h))

        w = open(filename, 'w').write

        w('#include "colors.inc"\n')
        w('#include "finish.inc"\n')

        w('camera {orthographic right %.2f*x up %.2f*y direction z\n' %
          (self.w, self.h))
        w('  location <0, 0, 10> look_at <0,0,0>}\n')
        
        w('light_source {<1, 2, 50> color White*2\n')
        w('  area_light <1, 0, 0>, <0, 1, 0>, 5, 5 adaptive 1 jitter}\n')
        
        w('background {color White}\n')

        z0 = self.X[:, 2].max()
        self.X -= (self.w / 2, self.h / 2, z0)
        
        if self.C is not None:
            self.C -= (self.w / 2, self.h / 2, z0)
            self.C.shape = (2, 2, 2, 3)
            if 1:
                for c in range(3):
                    for j in ([0, 0], [1, 0], [1, 1], [0, 1]):
                        w('cylinder {')
                        for i in range(2):
                            j.insert(c, i)
                            x, y, z = self.C[tuple(j)]
                            del j[c]
                            w('<%.2f, %.2f, %.2f>, ' % (-x, y, z))
                        w('0.05}\n')
            else:
                for c in range(3):
                    for i in range(2):
                        w('polygon {4,')
                        s = '\n  '
                        for j in ([0, 0], [1, 0], [1, 1], [0, 1]):
                            j.insert(c, i)
                            x, y, z = self.C[tuple(j)]
                            w(s + '<%.2f, %.2f, %.2f>' % (-x, y, z))
                            s = ',\n  '
                        w('texture {finish {ambient 1 diffuse 0}\n')
                        w('  pigment {rgbt <0.2, 0.2, 0.2, 0.8>}}}\n')

        w('#default {finish {phong 0.7}}\n')

        numbers = self.numbers
        done = {}
        for a in range(self.natoms):
            Z = numbers[a]
            if Z not in done:
                symbol = chemical_symbols[Z]
                w('#declare C_%s = ' % symbol +
                  'texture{pigment{rgb <%.2f, %.2f, %.2f>}}\n' %
                  tuple(self.colors[a]))
                w('#declare R_%s = %.2f;\n' % (symbol, self.d[a] / 2))
                done[Z] = True

        for Z, (x, y, z) in zip(numbers, self.X):
            symbol = chemical_symbols[Z]
            w('sphere {<%.2f, %.2f, %.2f>, ' % (-x, y, z) +
              'R_%s texture{C_%s}}\n' % (symbol, symbol))

def write_pov(filename, atoms, **parameters):
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]
    POVRAY(atoms, **parameters).write(filename)
