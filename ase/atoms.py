from math import cos, sin

import numpy as npy
from ase.data import atomic_numbers, chemical_symbols, atomic_masses

    
class Atom:
    def __init__(self, symbol, position=(0, 0, 0),
                 tag=None, magmom=None, mass=None,
                 momentum=None, charge=None):
        if isinstance(symbol, str):
            self.number = atomic_numbers[symbol]
        else:
            self.number = symbol
        self.position = position
        self.tag = tag
        self.momentum = momentum
        self.mass = mass
        self.magmom = magmom
        self.charge = charge

    def get_data(self):
        return (self.number, self.position,
                self.tag, self.momentum, self.mass,
                self.magmom, self.charge)
    

class Atoms(object):
    """Objecs representing a collection of atoms.

    The Atoms object can represent an isolated molecule, or a
    periodically repeated structure.  It may have a unit cell and
    there may be periodic boundary conditions along any of the three
    unit cell axes.

    The information about the atoms (atomic numbers and position) is
    stored in ndarrays.  Optionally, there can be information about
    tags, momenta, masses, magnetic moments and charges.

    In order to calculate energies, forces and stresses, a calculator
    object has to attached to the atoms object."""

    def __init__(self, atoms=None,
                 symbols=None, numbers=None, positions=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 cell=None, pbc=None,
                 constraints=None,
                 calculator=None):
        """Construct Atoms object.

        Parameters
        ----------
        atoms : list of Atom objects
            Example:  [Atom('Ne', (x, y, z), ...].
        symbols : string (formula) or list of strings
            Examples:  'H2O', 'COPt12', ['H', 'H', 'O'].
        numbers : list of integers
            Atomic numbers (use only one of symbols/numbers).
        positions : list of xyz-positions
            Atomic positions.  Anything that can be converted to an
            ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
            ...].
        tags : list of integers
        momenta: list of xyz-momenta
        masses : list of floats
        magmoms: list of floats
            Magnetic moments
        charges : list of floats
        cell : 3x3 matrix
            Unit cell vectors.  Can also be given as just three
            numbers for orthorhombic cells.  Default value: [1, 1, 1].
        pbc : one or three booleans
            Periodic boundary conditions flags.  Examples: True,
            False, 0, 1, (1, 1, 0), (True, False, False).  Default
            value: False.
        calculator : calculator object
            
        Examples
        --------
        These three are equivalent:

        >>> d = 1.104  # N2 bondlength
        >>> a = Atoms(symbols='N2', positions=[(0,0,0),(0,0,d)])
        >>> a = Atoms(numbers=[7,7], positions=[(0,0,0),(0,0,d)])
        >>> a = Atoms([Atom('N',(0,0,0)), Atom('N',(0,0,d)])

        FCC gold:

        >>> a = 4.05  # Gold lattice constant
        >>> b = a / 2
        >>> fcc = Atoms([Atom('Au')],
        ...             cell=[(0,b,b),(b,0,b),(b,b,0)],
        ...             pbs=True)

        Hydrogen wire:
        
        >>> d = 0.9  # H-H distance
        >>> L = 7.0
        >>> h = Atoms(symbols='H', positions=[(0,L/2,L/2)],
        ...           cell=(d, L, L),
        ...           pbc=(1,0,0))
        """

        if atoms is not None:
            if hasattr(atoms, 'GetUnitCell'):
                from ase.old import OldASEListOfAtomsWrapper
                atoms = OldASEListOfAtomsWrapper(atoms)

            elif isinstance(atoms, (list, tuple)):
                # Get data from a list or tuple of Atom objects:
                data = zip(*[atom.get_data() for atom in atoms])
                atoms = Atoms(None, None, *data)
                
            # Get data from another Atoms object:
            if symbols is None and numbers is None:
                numbers = atoms.get_atomic_numbers()
            if positions is None:
                positions = atoms.get_positions()
            if tags is None:
                tags = atoms.get_tags()
            if momenta is None:
                momenta = atoms.get_momenta()
            if magmoms is None:
                magmoms = atoms.get_magnetic_moments()
            if masses is None:
                masses = atoms.get_masses()
            if cell is None:
                cell = atoms.get_cell()
            if pbc is None:
                pbc = atoms.get_pbc()
            if constraints is None:
                constraints = [c.copy() for c in atoms.constraints]

        self.arrays = {}
        
        if positions is None:
            positions = npy.empty((0, 3))
            
        if symbols is None:
            if numbers is None:
                numbers = npy.zeros(len(positions), int)
            self.new_array('numbers', numbers, int)
        else:
            if numbers is not None:
                raise ValueError(
                    'Use only one of "symbols" and "numbers".')
            else:
                if isinstance(symbols, str):
                    symbols = string2symbols(symbols)
                numbers = []
                for s in symbols:
                    if isinstance(s, str):
                        numbers.append(atomic_numbers[s])
                    else:
                        numbers.append(s)
                self.new_array('numbers', numbers, int)

        self.new_array('positions', positions, float)

        self.set_tags(default(tags, 0))
        self.set_momenta(default(momenta, (0.0, 0.0, 0.0)))
        self.set_masses(default(masses, None))
        self.set_magnetic_moments(default(magmoms, 0.0))
        self.set_charges(default(charges, 0.0))

        if cell is None:
            cell = npy.eye(3)
        self.set_cell(cell, fix=True)

        if pbc is None:
            pbc = False
        self.set_pbc(pbc)

        self.set_constraint(constraints)
                
        self.set_calculator(calculator)

    def set_calculator(self, calc=None):
        """Attach calculator object."""
        if hasattr(calc, '_SetListOfAtoms'):
            from ase.old import OldASECalculatorWrapper
            calc = OldASECalculatorWrapper(calc, self)
        self.calc = calc

    def get_calculator(self):
        """Get currently attached calculator object."""
        return self.calc

    def set_constraint(self, constraints=None):
        if constraints is None:
            self.constraints = []
        else:
            if isinstance(constraints, (list, tuple)):
                self.constraints = constraints
            else:
                self.constraints = [constraints]
    
    def set_cell(self, cell, fix=False):
        """Set unit cell vectors.

        Parameters
        ----------
        cell : 
            Unit cell.  A 3x3 matrix (the three unit cell vectors) or
            just three numbers for an orthorhombic cell.
        fix : boolean
            Fix atomic positions or move atoms reletive to unit cell.
            Default behavior is to move the atoms (fix=False).

        Examples
        --------
        Two equivalent ways to define an orthorhombic cell:
        
        >>> a.set_cell([a, b, c])
        >>> a.set_cell([(a, 0, 0), (0, b, 0), (0, 0, c)])

        FCC unit cell:

        >>> a.set_cell([(0, b, b), (b, 0, b), (b, b, 0)])
        """

        cell = npy.array(cell, float)
        if cell.shape == (3,):
            cell = npy.diag(cell)
        elif cell.shape != (3, 3):
            raise ValueError('Cell must be length 3 sequence or '
                             '3x3 matrix!')
        if not fix:
            M = npy.linalg.solve(self.cell, cell)
            self.arrays['positions'][:] = npy.dot(self.arrays['positions'], M)
        self.cell = cell

    def get_cell(self):
        """Get the three unit cell vectors as an 3x3 ndarray."""
        return self.cell.copy()

    def set_pbc(self, pbc):
        """Set periodic boundary condition flags."""
        if isinstance(pbc, int):
            pbc = (pbc,) * 3
        self.pbc = npy.array(pbc, bool)
        
    def get_pbc(self):
        """Get periodic boundary condition flags."""
        return self.pbc.copy()

    def new_array(self, name, a, dtype=None):
        if dtype is not None:
            a = npy.array(a, dtype)
        else:
            a = a.copy()
            
        if name in self.arrays:
            raise RuntimeError

        for b in self.arrays.values():
            if len(a) != len(b):
                raise ValueError('Array has wrong length: %d != %d.' %
                                 (len(a), len(b)))
            break
        
        self.arrays[name] = a
    
    def set_array(self, a, name, dtype=None):
        b = self.arrays.get(name)
        if b is None:
            if a is not None:
                self.new_array(name, a, dtype)
        else:
            if a is None:
                del self.arrays[name]
            else:
                b[:] = a

    def get_atomic_numbers(self):
        """Get ndarray of all atomic numbers."""
        return self.arrays['numbers']

    def get_chemical_symbols(self):
        """Getlist of chemical symbols."""
        return [chemical_symbols[Z] for Z in self.arrays['numbers']]

    def set_tags(self, tags):
        self.set_array(tags, 'tags', int)
        
    def get_tags(self):
        return self.arrays.get('tags')

    def set_momenta(self, momenta):
        self.set_array(momenta, 'momenta', int)

    def get_momenta(self):
        return self.arrays.get('momenta')

    def set_masses(self, masses):
        if isinstance(masses, (list, tuple)):
            newmasses = []
            for m, Z in zip(masses, self.arrays['numbers']):
                if m is None:
                    newmasses.append(atomic_masses[Z])
                else:
                    newmasses.append(m)
            masses = newmasses
        self.set_array(masses, 'masses', float)

    def get_masses(self):
        return self.arrays.get('masses')

    def set_magnetic_moments(self, magmoms):
        self.set_array(magmoms, 'magmoms', float)

    def get_magnetic_moments(self):
        return self.arrays.get('magmoms')
    
    def set_charges(self, charges):
        self.set_array(charges, 'charges', int)

    def get_charges(self):
        return self.arrays.get('charges')

    def set_positions(self, newpositions):
        positions = self.arrays['positions']
        if self.constraints:
            newpositions = npy.asarray(newpositions, float)
            for constraint in self.constraints:
                constraint.adjust_positions(positions, newpositions)
                
        positions[:] = newpositions

    def get_positions(self):
        return self.arrays['positions'].copy()

    def get_potential_energy(self):
        if self.calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self.calc.get_potential_energy(self)

    def get_kinetic_energy(self):
        momenta = self.arrays.get('momenta')
        if momenta is None:
            return 0.0
        return 0.5 * npy.vdot(momenta, self.get_velocities())

    def get_velocities(self):
        momenta = self.arrays.get('momenta')
        if momenta is None:
            return None
        m = self.arrays.get('masses')
        if m is None:
            m = atomic_masses[self.arrays['numbers']]
        return momenta / m.reshape(-1, 1)
    
    def get_total_energy(self):
        return self.get_potential_energy() + self.get_kinetic_energy()

    def get_forces(self, apply_constraints=True):
        if self.calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        forces = self.calc.get_forces(self)
        if apply_constraints:
            for constraint in self.constraints:
                constraint.adjust_forces(self.arrays['positions'], forces)
        return forces

    def get_stress(self):
        if self.calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self.calc.get_stress(self)
    
    def copy(self):
        atoms = Atoms(cell=self.cell, pbc=self.pbc)

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a.copy()
            
        return atoms

    def __len__(self):
        return len(self.arrays['positions'])

    def __add__(self, other):
        atoms = self.copy()
        atoms += other
        return atoms

    def __iadd__(self, other):
        if isinstance(other, Atom):
            other = Atoms([other])
            
        n1 = len(self)
        n2 = len(other)
        
        for name, a1 in self.arrays.items():
            a = npy.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = npy.zeros((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            self.set_array(a, name)

        return self

    extend = __iadd__

    def append(self, atom):
        self.extend(Atoms([atom]))

    def __getitem__(self, i):
        if isinstance(i, int):
            i = [i]

        atoms = Atoms(cell=self.cell, pbc=self.pbc)

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a[i]
            
        return atoms

    def __delitem__(self, i):
        mask = npy.ones(len(self), bool)
        mask[i] = False
        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

    def pop(self, i=-1):
        atom = self[i]
        del self[i]
        return atom
    
    def __imul__(self, m):
        if isinstance(m, int):
            m = (m, m, m)
        M = npy.product(m)
        n = len(self)
        
        for name, a in self.arrays.items():
            self.arrays[name] = npy.tile(a, (M,) + (1,) * (len(a.shape) - 1))

        positions = self.arrays['positions']
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += npy.dot((m0, m1, m2), self.cell)
                    i0 = i1
        self.cell = npy.array([m[c] * self.cell[c] for c in range(3)])
        return self

    repeat = __imul__
    
    def __mul__(self, m):
        atoms = self.copy()
        atoms *= m
        return atoms

    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument can be a float an xyz vector or an
        nx3 array (where n is the number of atoms)."""

        self.arrays['positions'] += npy.array(displacement)

    def center(self, vacuum=None, axis=None):
        """Center atoms in unit cell"""
        p = self.arrays['positions']
        p0 = p.min(0)
        p1 = p.max(0)
        if axis is None:
            if vacuum is not None:
                self.cell = npy.diag(p1 - p0 + 2 * npy.asarray(vacuum))
            p += 0.5 * (self.cell.sum(0) - p0 - p1)
        else:
            c = self.cell.copy()
            c.flat[::4] = 0.0
            if c.any():
                raise ValueError('Unit cell must be orthorhobmic!')
            
            if vacuum is not None:
                self.cell[axis] = p1[axis] - p0[axis] + 2 * vacuum
            p[:, axis] += 0.5 * (self.cell[axis, axis] - p0[axis] - p1[axis])

    def get_center_of_mass(self):
        m = self.arrays.get('masses')
        if m is None:
            m = atomic_masses[self.arrays['numbers']]
        return npy.dot(m, self.arrays['positions']) / m.sum()

    def rotate(self, v, a=None):
        norm = npy.linalg.norm
        v = string2vector(v)
        if a is None:
            a = norm(v)
        if isinstance(a, (float, int)):
            v /= norm(v)
            c = cos(a)
            s = sin(a)
        else:
            v2 = string2vector(a)
            v /= norm(v)
            v2 /= norm(v2)
            c = npy.dot(v, v2)
            v = npy.cross(v, v2)
            s = norm(v)
            v /= s
        p = self.arrays['positions']
        p[:] = (c * p - 
                npy.cross(p, s * v) + 
                npy.outer(npy.dot(p, v), (1.0 - c) * v))

    def _get_positions(self):
        return self.arrays['positions']
    
    positions = property(_get_positions)

    #e = property(get_potential_energy)
    #c = property(get_cell, set_cell)
    #p = property(get_pbc, set_pbc)

    """
    center, repeat

    tag, magmom, mass, momentum

    add,
    """


def string2symbols(s):
    """Convert string to list of chemical symbols."""
    n = len(s)

    if n == 0:
        return []
    
    c = s[0]
    
    if c.isdigit():
        i = 1
        while i < n and s[i].isdigit():
            i += 1
        return int(s[:i]) * string2symbols(s[i:])

    if c == '(':
        p = 0
        for i, c in enumerate(s):
            if c == '(':
                p += 1
            elif c == ')':
                p -= 1
                if p == 0:
                    break
        j = i + 1
        while j < n and s[j].isdigit():
            j += 1
        if j > i + 1:
            m = int(s[i + 1:j])
        else:
            m = 1
        return m * string2symbols(s[1:i]) + string2symbols(s[j:])

    if c.isupper():
        i = 1
        if 1 < n and s[1].islower():
            i += 1
        j = i
        while j < n and s[j].isdigit():
            j += 1
        if j > i:
            m = int(s[i:j])
        else:
            m = 1
        return m * [s[:i]] + string2symbols(s[j:])

def string2vector(v):
    if isinstance(v, str):
        if v[0] == '-':
            return -string2vector(v[1:])
        w = npy.zeros(3)
        w['xyz'.index(v)] = 1.0
        return w
    return npy.asarray(v, float)

def default(data, dflt):
    """Helper function for setting default values."""
    if data is None:
        return None
    elif isinstance(data, (list, tuple)):
        newdata = []
        allnone = True
        for x in data:
            if x is None:
                newdata.append(dflt)
            else:
                newdata.append(x)
                allnone = False
        if allnone:
            return None
        return newdata
    else:
        return data
                               
