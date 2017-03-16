from __future__ import print_function
from math import sqrt

import numpy as np

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms
from ase.data import covalent_radii
from ase.gui.defaults import read_defaults
from ase.io import read, write, string2index


class Images:
    def __init__(self, images=None):
        if images is not None:
            self.initialize(images)

    def __len__(self):
        return len(self._images)

    def __getitem__(self, index):
        return self._images[index]

    def __iter__(self):
        return iter(self._images)

    # XXXXXXX hack
    # compatibility hacks while allowing variable number of atoms
    def get_dynamic(self, atoms):
        dynamic = np.ones(len(atoms), bool)
        for constraint in atoms.constraints:
            if isinstance(constraint, FixAtoms):
                dynamic[constraint.index] = False
        return dynamic

    def set_dynamic(self, indices, value):
        for atoms in self:
            dynamic = self.get_dynamic(atoms)
            dynamic[indices[indices < len(atoms)]] = value
            from ase.constraints import FixAtoms
            atoms.constraints = [c for c in atoms.constraints
                                 if not isinstance(c, FixAtoms)]
            atoms.constraints.append(FixAtoms(mask=~dynamic))

    def get_energy(self, atoms):
        try:
            e =  atoms.get_potential_energy() * self.repeat.prod()
        except RuntimeError:
            e = np.nan
        return e

    def get_forces(self, atoms):
        try:
            F = atoms.get_forces(apply_constraint=False)
        except RuntimeError:
            F = np.empty_like(atoms.positions)
            F.fill(np.nan)
        F = np.tile(F.T, self.repeat.prod()).T
        return F

    def get_magmoms(self, atoms, init_magmom=False):
        try:
            if init_magmom:
                M = atoms.get_initial_magnetic_moments()
            else:
                M = atoms.get_magnetic_moments()
                if M.ndim == 2:
                    M = M[:, 2]
        except (RuntimeError, AttributeError):
            M = atoms.get_initial_magnetic_moments()
        return M

    def initialize(self, images, filenames=None, init_magmom=False):

        #self.natoms = len(images[0])
        nimages = len(images)
        if filenames is None:
            filenames = [None] * nimages
        self.filenames = filenames

        #  The below seems to be about "quaternions"
        if 0: # XXXXXXXXXXXXXXXXXXXX hasattr(images[0], 'get_shapes'):
            self.Q = np.empty((nimages, self.natoms, 4))
            self.shapes = images[0].get_shapes()
            import os as os
            if os.path.exists('shapes'):
                shapesfile = open('shapes')
                lines = shapesfile.readlines()
                shapesfile.close()
                if '#{type:(shape_x,shape_y,shape_z), .....,}' in lines[0]:
                    shape = eval(lines[1])
                    shapes = []
                    for an in images[0].get_atomic_numbers():
                        shapes.append(shape[an])
                    self.shapes = np.array(shapes)
                else:
                    print('shape file has wrong format')
            else:
                print('no shapesfile found: default shapes were used!')

        else:
            self.shapes = None

        # Temporary hack for using atoms as backend instead of distinct arrays
        # This should probably be replaced by a more clean mechanism if we
        # decide to go this way
        class IndexHack:
            enclosing_images_obj = self

            def __init__(self, getter):
                self.getter = getter

            def __getitem__(self, item):
                if isinstance(item, tuple):
                    print(item)
                    raise ValueError(item)
                return self.getter(self.enclosing_images_obj[item])

            def __setitem__(self, item, value):
                raise ValueError(item, value)


        #self.P = IndexHack(lambda a: a.get_positions())
        #self.V = IndexHack(lambda a: a.get_velocities())
        #self.E = IndexHack(get_energy)
        #self.K = IndexHack(lambda a: a.get_kinetic_energy())
        #self.F = IndexHack(get_forces)
        #self.M = IndexHack(get_magmoms)
        #self.T = IndexHack(lambda a: a.get_tags())
        #self.A = IndexHack(lambda a: a.get_cell())
        #self.D = IndexHack(lambda a: a.get_celldisp().reshape((3,)))
        #self.Z = IndexHack(lambda a: a.get_atomic_numbers())
        #self.q = IndexHack(lambda a: a.get_initial_charges())
        #self.pbc = IndexHack(lambda a: a.get_pbc())
        #self.natoms = IndexHack(lambda a: len(a))
        #self.constrained = IndexHack(lambda a: get_constrained(a))

        #self.pbc = images[0].get_pbc()
        self.covalent_radii = covalent_radii.copy()
        self.config = read_defaults()
        # XXX config?

        #self.r = IndexHack(lambda a: np.array([self.covalent_radii[z]
        #                                       for z in a.numbers]))
        #if config['covalent_radii'] is not None:
        #    for data in config['covalent_radii']:
        #        self.covalent_radii[data[0]] = data[1]
        warning = False

        self._images = []

        # Whether length or chemical composition changes:
        self.have_varying_species = False
        for i, atoms in enumerate(images):
            # copy atoms or not?  Not copying allows back-editing,
            # but copying actually forgets things like the attached
            # calculator (might have forces/energies
            self._images.append(atoms)
            self.have_varying_species |= np.any(self[0].numbers
                                                != atoms.numbers)
            if hasattr(self, 'Q'):
                assert False # XXX askhl fix quaternions
                self.Q[i] = atoms.get_quaternions()
            if (atoms.pbc != self[0].pbc).any():
                warning = True

        if warning:
            import warnings
            warnings.warn('Not all images have the same boundary conditions!')

        self.maxnatoms = max(len(atoms) for atoms in self)
        self.selected = np.zeros(self.maxnatoms, bool)
        self.selected_ordered = []
        self.visible = np.ones(self.maxnatoms, bool)
        self.nselected = 0
        self.repeat = np.ones(3, int)

    def get_radii(self, atoms):
        radii = np.array([self.covalent_radii[z] for z in atoms.numbers])
        radii *= self.config['radii_scale']
        return radii

    def prepare_new_atoms(self):
        "Marks that the next call to append_atoms should clear the images."
        self.next_append_clears = True

    def append_atoms(self, atoms, filename=None):
        "Append an atoms object to the images already stored."
        self.images.append(atoms)
        self.filenames.append(filename)
        self.initialize(self.images, filenames=self.filenames)
        return

        sdjkfskdjfsdkjf
        assert len(atoms) == self.natoms
        if self.next_append_clears:
            i = 0
        else:
            i = len(self)
        for name in ('P', 'V', 'E', 'K', 'F', 'M', 'A', 'T', 'D', 'q'):
            a = getattr(self, name)
            newa = np.empty((i + 1,) + a.shape[1:], a.dtype)
            if not self.next_append_clears:
                newa[:-1] = a
            setattr(self, name, newa)
        self.next_append_clears = False
        self.P[i] = atoms.get_positions()
        self.V[i] = atoms.get_velocities()
        self.A[i] = atoms.get_cell()
        self.D[i] = atoms.get_celldisp().reshape((3,))
        self.q[i] = atoms.get_initial_charges()
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
        try:
            self.T[i] = atoms.get_tags()
        except AttributeError:
            if i == 0:
                self.T[i] = 0
            else:
                self.T[i] = self.T[i - 1]
        #self.nimages = i + 1
        self.filenames.append(filename)
        #return self.nimages

    #def set_radii(self, scale=None):
    #    if scale is None:
    #        self.covalent_radii = covalent_radii.copy()
    #    else:
    #        self.covalent_radii *= scale

    def read(self, filenames, index=-1, filetype=None):
        images = []
        names = []
        for filename in filenames:
            i = read(filename, index, filetype)

            if not isinstance(i, list):
                i = [i]
            images.extend(i)
            names.extend([filename] * len(i))

        self.initialize(images, names)

        for image in images:
            if 'radii' in image.info:
                #self.set_radii(image.info['radii'])
                break

    def repeat_unit_cell(self):
        for atoms in self:
            # Get quantities taking into account current repeat():
            ref_energy = self.get_energy(atoms)
            ref_forces = self.get_forces(atoms)
            atoms.calc = SinglePointCalculator(atoms,
                                               energy=ref_energy,
                                               forces=ref_forces)
            atoms.cell *= self.repeat.reshape((3, 1))
        self.repeat = np.ones(3, int)

    def repeat_images(self, repeat):
        repeat = np.array(repeat)
        oldprod = self.repeat.prod()
        images = []
        for atoms in self:
            refcell = atoms.get_cell()
            atoms = atoms[:len(atoms) // oldprod]
            atoms *= repeat
            atoms.cell = refcell
            images.append(atoms)
        self.initialize(images, filenames=self.filenames)
        self.repeat = repeat

        return

        # XXXXXXXXXXXXXX disabled repeat code below

        n = self.repeat.prod()
        repeat = np.array(repeat)
        self.repeat = repeat
        N = repeat.prod()
        natoms = self.natoms // n
        P = np.empty((self.nimages, natoms * N, 3))
        V = np.empty((self.nimages, natoms * N, 3))
        M = np.empty((self.nimages, natoms * N))
        T = np.empty((self.nimages, natoms * N), int)
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
                    V[:, a0:a1] = self.V[:, :natoms]
                    F[:, a0:a1] = self.F[:, :natoms]
                    M[:, a0:a1] = self.M[:, :natoms]
                    T[:, a0:a1] = self.T[:, :natoms]
                    Z[a0:a1] = self.Z[:natoms]
                    r[a0:a1] = self.r[:natoms]
                    dynamic[a0:a1] = self.dynamic[:natoms]
                    a0 = a1
        self.P = P
        self.V = V
        self.F = F
        self.Z = Z
        self.T = T
        self.M = M
        self.r = r
        self.dynamic = dynamic
        self.natoms = natoms * N
        self.selected = np.zeros(natoms * N, bool)
        # XXX disabled askhl self.atoms_to_rotate_0 = np.zeros(self.natoms, bool)
        self.visible = np.ones(natoms * N, bool)
        self.nselected = 0

    def center(self):
        """Center each image in the existing unit cell, keeping the
        cell constant."""
        for atoms in self:
            atoms.center()
        #c = self.A.sum(axis=1) / 2.0 - self.P.mean(axis=1)
        #self.P += c[:, np.newaxis, :]

    def graph(self, expr):
        """Routine to create the data in ase-gui graphs, defined by the
        string expr."""
        import ase.units as units
        code = compile(expr + ',', '<input>', 'eval')

        n = len(self)

        def d(n1, n2):
            return sqrt(((R[n1] - R[n2])**2).sum())

        def a(n1, n2, n3):
            v1 = R[n1] - R[n2]
            v2 = R[n3] - R[n2]
            arg = np.vdot(v1, v2) / (sqrt((v1**2).sum() * (v2**2).sum()))
            if arg > 1.0:
                arg = 1.0
            if arg < -1.0:
                arg = -1.0
            return 180.0 * np.arccos(arg) / np.pi

        def dih(n1, n2, n3, n4):
            # vector 0->1, 1->2, 2->3 and their normalized cross products:
            a = R[n2] - R[n1]
            b = R[n3] - R[n2]
            c = R[n4] - R[n3]
            bxa = np.cross(b, a)
            bxa /= np.sqrt(np.vdot(bxa, bxa))
            cxb = np.cross(c, b)
            cxb /= np.sqrt(np.vdot(cxb, cxb))
            angle = np.vdot(bxa, cxb)
            # check for numerical trouble due to finite precision:
            if angle < -1:
                angle = -1
            if angle > 1:
                angle = 1
            angle = np.arccos(angle)
            if np.vdot(bxa, c) > 0:
                angle = 2 * np.pi - angle
            return angle * 180.0 / np.pi

        # get number of mobile atoms for temperature calculation
        #ndynamic = self.dynamic.sum()
        #self.constrained[

        #D = self.dynamic[:, np.newaxis]
        E = np.array([self.get_energy(atoms) for atoms in self])

        s = 0.0

        # Namespace for eval:
        ns = {'E': E,
              'd': d, 'a': a, 'dih': dih}

        data = []
        for i in range(n):
            ns['i'] = i
            ns['s'] = s
            ns['R'] = R = self[i].get_positions()
            ns['V'] = self[i].get_velocities()
            ns['F'] = F = self.get_forces(self[i])
            ns['A'] = self[i].get_cell()
            ns['M'] = self[i].get_masses()
            # XXX askhl verify:
            dynamic = self.get_dynamic(self[i])
            ns['f'] = f = ((F * dynamic[:, None])**2).sum(1)**.5
            ns['fmax'] = max(f)
            ns['fave'] = f.mean()
            ns['epot'] = epot = E[i]
            ns['ekin'] = ekin = self[i].get_kinetic_energy()
            ns['e'] = epot + ekin
            ndynamic = dynamic.sum()
            ns['T'] = 2.0 * ekin / (3.0 * ndynamic * units.kB)
            data = eval(code, ns)
            if i == 0:
                m = len(data)
                xy = np.empty((m, n))
            xy[:, i] = data
            if i + 1 < n and not self.have_varying_species:
                s += sqrt(((self[i + 1].positions - R)**2).sum())
        return xy

    def write(self, filename, rotations='', show_unit_cell=False, bbox=None,
              **kwargs):
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
                  bbox=bbox, **kwargs)
        else:
            write(filename, images, **kwargs)

    def get_atoms(self, frame, remove_hidden=False):
        #atoms = Atoms(positions=self.P[frame],
        #              numbers=self.Z,
        #              magmoms=self.M[0],
        #              tags=self.T[frame],
        #              cell=self.A[frame],
        #              pbc=self.pbc)
        atoms = self[frame]
        try:
            E = atoms.get_potential_energy()
        except RuntimeError:
            E = None
        try:
            F = atoms.get_forces()
        except RuntimeError:
            F = None
        #if not np.isnan(self.V).any():
        #    atoms.set_velocities(self.V[frame])

        # check for constrained atoms and add them accordingly:
        #if not self.dynamic.all():
        #    atoms.set_constraint(FixAtoms(mask=1 - self.dynamic))

        # Remove hidden atoms if applicable
        if remove_hidden:
            atoms = atoms[self.visible]
            if F is not None:
                F = F[self.visible]
            #f = self.F[frame][self.visible]
        #else:
            #f = self.F[frame]
        atoms.set_calculator(SinglePointCalculator(atoms,
                                                   energy=E,
                                                   forces=F))
        return atoms

    def delete(self, i):
        self.images.pop(i)
        self.filenames.pop(i)
        self.initialize(self.images, self.filenames)
        return

        sdfsdfksdflksdfdsf
        self.nimages -= 1
        P = np.empty((self.nimages, self.natoms, 3))
        V = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        A = np.empty((self.nimages, 3, 3))
        E = np.empty(self.nimages)
        P[:i] = self.P[:i]
        P[i:] = self.P[i + 1:]
        self.P = P
        V[:i] = self.V[:i]
        V[i:] = self.V[i + 1:]
        self.V = V
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
        skjdfsdkjfsdkfjskdjfk
        n = self.nimages
        assert n % 5 == 0
        levels = n // 5
        n = self.nimages = 2 * levels + 3
        P = np.empty((self.nimages, self.natoms, 3))
        V = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        E = np.empty(self.nimages)
        for L in range(levels):
            P[L] = self.P[L * 5]
            P[n - L - 1] = self.P[L * 5 + 4]
            V[L] = self.V[L * 5]
            V[n - L - 1] = self.V[L * 5 + 4]
            F[L] = self.F[L * 5]
            F[n - L - 1] = self.F[L * 5 + 4]
            E[L] = self.E[L * 5]
            E[n - L - 1] = self.E[L * 5 + 4]
        for i in range(3):
            P[levels + i] = self.P[levels * 5 - 4 + i]
            V[levels + i] = self.V[levels * 5 - 4 + i]
            F[levels + i] = self.F[levels * 5 - 4 + i]
            E[levels + i] = self.E[levels * 5 - 4 + i]
        self.P = P
        self.V = V
        self.F = F
        self.E = E

    def interpolate(self, m):
        sdkfjsdkfjsdkjf
        assert self.nimages == 2
        self.nimages = 2 + m
        P = np.empty((self.nimages, self.natoms, 3))
        V = np.empty((self.nimages, self.natoms, 3))
        F = np.empty((self.nimages, self.natoms, 3))
        A = np.empty((self.nimages, 3, 3))
        E = np.empty(self.nimages)
        T = np.empty((self.nimages, self.natoms), int)
        D = np.empty((self.nimages, 3))
        P[0] = self.P[0]
        V[0] = self.V[0]
        F[0] = self.F[0]
        A[0] = self.A[0]
        E[0] = self.E[0]
        T[:] = self.T[0]
        for i in range(1, m + 1):
            x = i / (m + 1.0)
            y = 1 - x
            P[i] = y * self.P[0] + x * self.P[1]
            V[i] = y * self.V[0] + x * self.V[1]
            F[i] = y * self.F[0] + x * self.F[1]
            A[i] = y * self.A[0] + x * self.A[1]
            E[i] = y * self.E[0] + x * self.E[1]
            D[i] = y * self.D[0] + x * self.D[1]
        P[-1] = self.P[1]
        V[-1] = self.V[1]
        F[-1] = self.F[1]
        A[-1] = self.A[1]
        E[-1] = self.E[1]
        D[-1] = self.D[1]
        self.P = P
        self.V = V
        self.F = F
        self.A = A
        self.E = E
        self.T = T
        self.D = D
        self.filenames[1:1] = [None] * m
