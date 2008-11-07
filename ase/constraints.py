from math import sqrt

import numpy as np

def slice2enlist(s):
    """Convert a slice object into a list of (new, old) tuples."""
    if isinstance(s, (list, tuple)):
        return enumerate(s)
    if s.step == None:
        step = 1
    else:
        step = s.step
    return enumerate(range(s.start, s.stop, step))

class FixConstraint(object):
    """Base class for classes that fix one or more atoms in some way."""

    def index_shuffle(self, ind):
        """Change the indices.

        When the ordering of the atoms in the Atoms object changes,
        this method can be called to shuffle the indices of the
        constraints.

        ind -- List or tuple of indices.

        """
        pass

class FixConstraintSingle(FixConstraint):
    "Base class for classes that fix a single atom."

    def index_shuffle(self, ind):
        "The atom index must be stored as self.a."
        newa = -1 # Signal error
        for new, old in slice2enlist(ind):
            if old == self.a:
                newa = new
                break
        if newa == -1:
            raise IndexError('Constraint not part of slice')
        self.a = newa

class FixAtoms(FixConstraint):
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
            self.index = np.asarray(mask, bool)
        else:
            self.index = np.asarray(indices, int)

        if self.index.ndim != 1:
            raise ValueError('Wrong argument to FixAtoms class!')

    def adjust_positions(self, old, new):
        new[self.index] = old[self.index]

    def adjust_forces(self, positions, forces):
        forces[self.index] = 0.0

    def index_shuffle(self, ind):
        # See docstring of superclass
        if self.index.dtype == bool:
            self.index = self.index[ind]
        else:
            index = []
            for new, old in slice2enlist(ind):
                if old in self.index:
                    index.append(new)
            if len(index) == 0:
                raise IndexError('All indices in FixAtoms not part of slice')
            self.index = np.asarray(index, int)

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


class FixBondLength(FixConstraint):
    """Constraint object for fixing a bond length."""
    def __init__(self, a1, a2):
        """Fix distance between atoms with indices a1 and a2."""
        self.indices = [a1, a2]

    def adjust_positions(self, old, new):        
        p1, p2 = old[self.indices]
        d = p2 - p1
        p = sqrt(np.dot(d, d))
        q1, q2 = new[self.indices]
        d = q2 - q1
        q = sqrt(np.dot(d, d))
        d *= 0.5 * (p - q) / q
        new[self.indices] = (q1 - d, q2 + d)

    def adjust_forces(self, positions, forces):
        d = np.subtract.reduce(positions[self.indices])
        d2 = np.dot(d, d)
        d *= 0.5 * np.dot(np.subtract.reduce(forces[self.indices]), d) / d2
        forces[self.indices] += (-d, d)

    def copy(self):
        return FixBondLength(*self.indices)

    def __repr__(self):
        return 'FixBondLength(%d, %d)' % tuple(self.indices)


class FixedPlane(FixConstraintSingle):
    """Constrain an atom *a* to move in a given plane only.

    The plane is defined by its normal: *direction*."""
    
    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def adjust_positions(self, oldpositions, newpositions):
        step = newpositions[self.a] - oldpositions[self.a]
        newpositions[self.a] -= self.dir * np.dot(step, self.dir)

    def adjust_forces(self, positions, forces):
        forces[self.a] -= self.dir * np.dot(forces[self.a], self.dir)

    def copy(self):
        return FixedPlane(self.a, self.dir)

    def __repr__(self):
        return 'FixedPlane(%d, %s)' % (self.a, self.dir.tolist())


class FixedLine(FixConstraintSingle):
    """Constrain an atom *a* to move on a given line only.

    The line is defined by its *direction*."""
    
    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def adjust_positions(self, oldpositions, newpositions):
        step = newpositions[self.a] - oldpositions[self.a]
        x = np.dot(step, self.dir)
        newpositions[self.a] = oldpositions[self.a] + x * self.dir

    def adjust_forces(self, positions, forces):
        forces[self.a] = self.dir * np.dot(forces[self.a], self.dir)

    def copy(self):
        return FixedLine(self.a, self.dir)

    def __repr__(self):
        return 'FixedLine(%d, %s)' % (self.a, self.dir.tolist())

class FixCartesian(FixConstraintSingle):
    "Fix an atom in the directions of the cartesian coordinates."
    def __init__(self, a, mask=[1,1,1]):
        self.a=a
        self.mask = -(np.array(mask)-1)
    
    def adjust_positions(self, old, new):
        step = new - old
        step[self.a] *= self.mask
        new = old + step

    def adjust_forces(self, positions, forces):
        forces[self.a] *= self.mask

    def copy(self):
        return fix_cartesian(self.a, self.mask)

    def __repr__(self):
        return 'FixCartesian(indice=%s mask=%s)' % (self.a, self.mask)

class fix_cartesian(FixCartesian):
    "Backwards compatibility for FixCartesian."
    def __init__(self, a, mask=[1,1,1]):
        import warnings
        super(fix_cartesian, self).__init__(a, mask)
        warnings.warn('fix_cartesian is deprecated. Please use FixCartesian' \
                      ' instead.', DeprecationWarning, stacklevel=2)

class FixScaled(FixConstraintSingle):
    "Fix an atom in the directions of the unit vectors."
    def __init__(self, cell, a, mask=[1,1,1]):
        self.cell = cell
        self.a = a
        self.mask = np.array(mask)

    def adjust_positions(self, old, new):
        scaled_old = np.dot(old, np.linalg.inv(self.cell))
        scaled_new = np.dot(new, np.linalg.inv(self.cell))
        for n in range(3):
            if self.mask[n]:
                scaled_new[self.a, n] = scaled_old[self.a, n]
        new[self.a] = np.dot(scaled_new, self.cell)[self.a]

    def adjust_forces(self, positions, forces):
        scaled_forces = np.dot(forces, np.linalg.inv(self.cell))
        scaled_forces[self.a] *= -(self.mask-1)
        forces[self.a] = np.dot(scaled_forces, self.cell)[self.a]


    def copy(self):
        return fix_scaled(self.cell ,self.a, self.mask)

    def __repr__(self):
        return 'FixScaled(indice=%s mask=%s)' % (self.a, self.mask)

class fix_scaled(FixScaled):
    "Backwards compatibility for FixScaled."
    def __init__(self, cell, a, mask=[1,1,1]):
        import warnings
        super(fix_scaled, self).__init__(cell, a, mask)
        warnings.warn('fix_scaled is deprecated. Please use FixScaled ' \
                      'instead.', DeprecationWarning, stacklevel=2)

class Filter:
    """Subset filter class."""
    def __init__(self, atoms, indices=None, mask=None):
        """Filter atoms.

        This filter can be used to hide degrees of freedom in an Atoms
        object.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should be constrained.
        mask : list of bool
           One boolean per atom indicating if the atom should be
           constrained or not.
        """

        self.atoms = atoms

        if indices is None and mask is None:
            raise ValuError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValuError('Use only one of "indices" and "mask".')

        if mask is not None:
            self.index = np.asarray(mask, bool)
        else:
            self.index = np.asarray(indices, int)

    def get_positions(self):
        return self.atoms.get_positions()[self.index]

    def set_positions(self, positions):
        pos = self.atoms.get_positions()
        pos[self.index] = positions
        self.atoms.set_positions(positions)

    def get_forces(self):
        return self.atoms.get_forces()[self.index]


class StrainFilter:
    """Modify the supercell while keeping the scaled positions fixed.

    Presents the strain of the supercell as the generalized positions,
    and the global stress tensor (times the volume) as the generalized
    force.

    This filter can be used to relax the unit cell until the stress is
    zero.  If MDMin is used for this, the timestep (dt) to be used
    depends on the system size. 0.01/x where x is a typical dimension
    seems like a good choice.

    The stress and strain are presented as 6-vectors, the order of the
    components follow the standard engingeering practice: xx, yy, zz,
    yz, xz, xy.  
    
    """
    def __init__(self, atoms, mask=None):
        """Create a filter applying a homogeneous strain to a list of atoms.

        The first argument, atoms, is the atoms object.

        The optional second argument, mask, is a list of six booleans,
        indicating which of the six independent components of the
        strain that are allowed to become non-zero.  It defaults to
        [1,1,1,1,1,1].
        
        """
        
        self.atoms = atoms
        self.strain = np.zeros(6)

        if mask is None:
            self.mask = np.ones(6)
        else:
            self.mask = np.array(mask)

        self.origcell = atoms.get_cell()
        
    def get_positions(self):
        return self.strain.reshape((2, 3))

    def set_positions(self, new):
        new = new.ravel() * self.mask
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]])

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new

    def get_forces(self):
        stress = self.atoms.get_stress()
        return -self.atoms.get_volume() * (stress * self.mask).reshape((2, 3))

    def get_potential_energy(self):
        return self.atoms.get_potential_energy()

    def __len__(self):
        return 2
