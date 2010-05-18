from math import sqrt

import numpy as np

from ase.data import covalent_radii
from ase.atoms import Atoms
from ase.calculators import SinglePointCalculator
from ase.io import read, write, string2index


class Images:
    def __init__(self, images=None):

        if images is not None:
            self.initialize(images)
    
    def initialize(self, images, filenames=None, init_magmom=False):
        
        self.natoms = len(images[0])

        self.nimages = len(images)
        if filenames is None:
            filenames = [None] * self.nimages
        self.filenames = filenames
        self.P = np.empty((self.nimages, self.natoms, 3))
        self.E = np.empty(self.nimages)
        self.K = np.empty(self.nimages)
        self.F = np.empty((self.nimages, self.natoms, 3))
        self.M = np.empty((self.nimages, self.natoms))
        self.T = np.empty((self.natoms))
        self.A = np.empty((self.nimages, 3, 3))
        self.Z = images[0].get_atomic_numbers()
        self.pbc = images[0].get_pbc()
        warning = False
        for i, atoms in enumerate(images):
            natomsi = len(atoms)
            if (natomsi != self.natoms or
                (atoms.get_atomic_numbers() != self.Z).any()):
                raise RuntimeError('Can not handle different images with ' +
                                   'different numbers of atoms or different ' +
                                   'kinds of atoms!')
            self.P[i] = atoms.get_positions()
            self.A[i] = atoms.get_cell()
            if (atoms.get_pbc() != self.pbc).any():
                warning = True
            try:
                self.E[i] = atoms.get_potential_energy()
            except RuntimeError:
                self.E[i] = np.nan
            self.K[i] = atoms.get_kinetic_energy()
            try:
                self.F[i] = atoms.get_forces(apply_constraint=False)
            except RuntimeError:
                self.F[i] = np.nan
            try:
                if init_magmom:
                    self.M[i] = atoms.get_initial_magnetic_moments()
                else:
                  self.M[i] = atoms.get_magnetic_moments()
            except (RuntimeError, AttributeError):
                self.M[i] = 0.0
                
            # added support for tags
            try:
                self.T = atoms.get_tags()
            except RuntimeError:
                self.T = np.nan
                

        if warning:
            print('WARNING: Not all images have the same bondary conditions!')
            
        self.selected = np.zeros(self.natoms, bool)
        self.selected_ordered  = []
        self.atoms_to_rotate_0 = np.zeros(self.natoms, bool)
        self.visible = np.ones(self.natoms, bool)
        self.nselected = 0
        self.set_dynamic()
        self.repeat = np.ones(3, int)
        self.set_radii(0.89)
        
    def prepare_new_atoms(self):
        "Marks that the next call to append_atoms should clear the images."
        self.next_append_clears = True
        
    def append_atoms(self, atoms, filename=None):
        "Append an atoms object to the images already stored."
        assert len(atoms) == self.natoms
        if self.next_append_clears:
            i = 0
        else:
            i = self.nimages
        for name in ('P', 'E', 'K', 'F', 'M', 'A'):
            a = getattr(self, name)
            newa = np.empty( (i+1,) + a.shape[1:] )
            if not self.next_append_clears:
                newa[:-1] = a
            setattr(self, name, newa)
        self.next_append_clears = False
        self.P[i] = atoms.get_positions()
        self.A[i] = atoms.get_cell()
        try:
            self.E[i] = atoms.get_potential_energy()
        except RuntimeError:
            self.E[i] = np.nan
        self.K[i] = atoms.get_kinetic_energy()
        try:
            self.F[i] = atoms.get_forces(apply_constraint=False)
        except RuntimeError:
            self.F[i] = np.nan
        try:
            self.M[i] = atoms.get_magnetic_moments()
        except (RuntimeError, AttributeError):
            self.M[i] = np.nan
        self.nimages = i + 1
        self.filenames.append(filename)
        self.set_dynamic()
        return self.nimages
        
    def set_radii(self, scale):
        self.r = covalent_radii[self.Z] * scale
                
    def read(self, filenames, index=-1):
        images = []
        names = []
        for filename in filenames:
            i = read(filename, index)
            
            if not isinstance(i, list):
                i = [i]
            images.extend(i)
            names.extend([filename] * len(i))
            
        self.initialize(images, names)
    
    def import_atoms(self, filename, cur_frame):
        if filename:
            filename = filename[0]
            old_a = self.get_atoms(cur_frame)
            imp_a = read(filename, -1)
            new_a = old_a + imp_a
            self.initialize([new_a], [filename])
    
    def repeat_images(self, repeat):
        n = self.repeat.prod()
        repeat = np.array(repeat)
        self.repeat = repeat
        N = repeat.prod()
        natoms = self.natoms // n
        P = np.empty((self.nimages, natoms * N, 3))
        M = np.empty((self.nimages, natoms * N))
        T = np.empty(natoms * N, int)
        F = np.empty((self.nimages, natoms * N, 3))
        Z = np.empty(natoms * N, int)
        r = np.empty(natoms * N)
        dynamic = np.empty(natoms * N, bool)
        a0 = 0
        for i0 in range(repeat[0]):
            for i1 in range(repeat[1]):
                for i2 in range(repeat[2]):
                    a1 = a0 + natoms
                    for i in range(self.nimages):
                        P[i, a0:a1] = (self.P[i, :natoms] +
                                       np.dot((i0, i1, i2), self.A[i]))
                    F[:, a0:a1] = self.F[:, :natoms]
                    M[:, a0:a1] = self.M[:, :natoms]
                    T[a0:a1] = self.T[:natoms]
                    Z[a0:a1] = self.Z[:natoms]
                    r[a0:a1] = self.r[:natoms]
                    dynamic[a0:a1] = self.dynamic[:natoms]
                    a0 = a1
        self.P = P
        self.F = F
        self.Z = Z
        self.T = T
        self.M = M
        self.r = r
        self.dynamic = dynamic
        self.natoms = natoms * N
        self.selected = np.zeros(natoms * N, bool)
        self.atoms_to_rotate_0 = np.zeros(self.natoms, bool)
        self.visible = np.ones(natoms * N, bool)
        self.nselected = 0
        
    def graph(self, expr):
        code = compile(expr + ',', 'atoms.py', 'eval')

        n = self.nimages
        def d(n1, n2):
            return sqrt(((R[n1] - R[n2])**2).sum())
        S = self.selected
        D = self.dynamic[:, np.newaxis]
        E = self.E
        s = 0.0
        data = []
        for i in range(n):
            R = self.P[i]
            F = self.F[i]
            A = self.A[i]
            f = ((F * D)**2).sum(1)**.5
            fmax = max(f)
            fave = f.mean()
            epot = E[i]
            ekin = self.K[i]
            e = epot + ekin
            data = eval(code)
            if i == 0:
                m = len(data)
                xy = np.empty((m, n))
            xy[:, i] = data
            if i + 1 < n:
                s += sqrt(((self.P[i + 1] - R)**2).sum())
        return xy

    def set_dynamic(self):
        if self.nimages == 1:
            self.dynamic = np.ones(self.natoms, bool)
        else:
            self.dynamic = np.zeros(self.natoms, bool)
            R0 = self.P[0]
            for R in self.P[1:]:
                self.dynamic |= (np.abs(R - R0) > 1.0e-10).any(1)

    def write(self, filename, rotations='', show_unit_cell=False, bbox=None):
        indices = range(self.nimages)
        p = filename.rfind('@')
        if p != -1:
            try:
                slice = string2index(filename[p + 1:])
            except ValueError:
                pass
            else:
                indices = indices[slice]
                filename = filename[:p]
                if isinstance(indices, int):
                    indices = [indices]

        images = [self.get_atoms(i) for i in indices]
        if len(filename) > 4 and filename[-4:] in ['.eps', '.png', '.pov']:
            write(filename, images, 
                  rotation=rotations, show_unit_cell=show_unit_cell,
                  bbox=bbox)
        else:
            write(filename, images)

    def get_atoms(self, frame):
        atoms = Atoms(positions=self.P[frame],
                      numbers=self.Z,
                      magmoms=self.M[0],
                      tags=self.T,
                      cell=self.A[frame],
                      pbc=self.pbc)
        atoms.set_calculator(SinglePointCalculator(self.E[frame],
                                                   self.F[frame],
                                                   None, None, atoms))
        return atoms
                           
    def delete(self, i):
        self.nimages -= 1
        P = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        A = np.empty((self.nimages, 3, 3))
        E = np.empty(self.nimages)
        P[:i] = self.P[:i]
        P[i:] = self.P[i + 1:]
        self.P = P
        F[:i] = self.F[:i]
        F[i:] = self.F[i + 1:]
        self.F = F
        A[:i] = self.A[:i]
        A[i:] = self.A[i + 1:]
        self.A = A
        E[:i] = self.E[:i]
        E[i:] = self.E[i + 1:]
        self.E = E
        del self.filenames[i]

    def aneb(self):
        n = self.nimages
        assert n % 5 == 0
        levels = n // 5
        n = self.nimages = 2 * levels + 3
        P = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        E = np.empty(self.nimages)
        for L in range(levels):
            P[L] = self.P[L * 5]
            P[n - L - 1] = self.P[L * 5 + 4]
            F[L] = self.F[L * 5]
            F[n - L - 1] = self.F[L * 5 + 4]
            E[L] = self.E[L * 5]
            E[n - L - 1] = self.E[L * 5 + 4]
        for i in range(3):
            P[levels + i] = self.P[levels * 5 - 4 + i]
            F[levels + i] = self.F[levels * 5 - 4 + i]
            E[levels + i] = self.E[levels * 5 - 4 + i]
        self.P = P
        self.F = F
        self.E = E

    def interpolate(self, m):
        assert self.nimages == 2
        self.nimages = 2 + m
        P = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        A = np.empty((self.nimages, 3, 3))
        E = np.empty(self.nimages)
        P[0] = self.P[0]
        F[0] = self.F[0]
        A[0] = self.A[0]
        E[0] = self.E[0]
        for i in range(1, m + 1):
            x = i / (m + 1.0)
            y = 1 - x
            P[i] = y * self.P[0] + x * self.P[1]
            F[i] = y * self.F[0] + x * self.F[1]
            A[i] = y * self.A[0] + x * self.A[1]
            E[i] = y * self.E[0] + x * self.E[1]
        P[-1] = self.P[1]
        F[-1] = self.F[1]
        A[-1] = self.A[1]
        E[-1] = self.E[1]
        self.P = P
        self.F = F
        self.A = A
        self.E = E
        self.filenames[1:1] = [None] * m

if __name__ == '__main__':
    import os
    os.system('python gui.py')
