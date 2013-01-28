import os
import subprocess
from math import pi, sqrt

import numpy as np


class NotAvailable(Exception):
    pass


def equal(a, b):
    """ndarray-enabled comparison function."""
    if isinstance(a, np.ndarray):
        b = np.array(b)
        return a.shape == b.shape and (a == b).all()
    if isinstance(b, np.ndarray):
        return equal(b, a)
    return a == b


def read_parameters_from_file(filename, calculator=None):
    """Read parameters from file.

    Get the function called 'parameters' and call it XXX
    """
    namespace = {}
    execfile(os.path.expanduser(filename), namespace)
    return namespace['parameters'](calculator)


def kptdensity2monkhorstpack(atoms, kptdensity=3.5, even=True):
    """Convert k-point density to Monkhorst-Pack grid size.

    atoms: Atoms object
        Contains unit cell and information about boundary conditions.
    kptdensity: float
        K-point density.  Default value is 3.5 point per Ang^-1.
    even: bool
        Round to even numbers.
    """

    recipcell = atoms.get_reciprocal_cell()
    kpts = []
    for i in range(3):
        if atoms.pbc[i]:
            k = 2 * pi * sqrt((recipcell[i]**2).sum()) * kptdensity
            if even:
                kpts.append(max(1, 2 * int(round(k / 2))))
            else:
                kpts.append(max(1, int(round(k))))
        else:
            kpts.append(1)
    return kpts


def normalize_smearing_keyword(smearing):
    """Normalize smearing string to long names and lower case.
    """

    smearing = smearing.lower()
    if smearing == 'fd':
        smearing = 'fermi-dirac'
    elif smearing.startswith('mp'):
        smearing = 'mathfessel-paxton-' + smearing[2]
    return smearing


class Calculator:
    notimplemented = []  # properties calculator can't handle

    """Base-class for all ASE calculators."""

    def __init__(self, path=None, atoms=None, **kwargs):
        """Basic calculator implementation.

        path: str
            Path for output file.
        atoms: Atoms object
            Optional Atoms object to which the calculator will be
            attached.  If path exists, atoms will get its positions
            and unit-cell updated form file.
        """

        self.state = None  # copy of atoms object for last calculation
        self.results = {}  # calculated properties (energy, forces, ...)
        self.parameters = {}  # calculational parameters

        self.path = path
        if path is not None:
            self.read(path)

        if atoms is not None:
            atoms.calc = self
            if self.state is not None:
                # State was read from file.  Update atoms:
                if not (equal(atoms.numbers, self.state.numbers) and
                        (atoms.pbc == self.state.pbc).all()):
                    raise RuntimeError('Atoms not compatible with file')
                atoms.positions = self.state.positions
                atoms.cell = self.state.cell
                
        self.set(**kwargs)

        self.hooks = {'before': [], 'after': []}  # call-back functions

        self.name = self.__class__.__name__

    def reset(self, changed_parameters=[]):
        """Clear all information from old calculation.

        Subclasses can decide to bypass this if the changed_parameters
        are harmless like a change in verbosity."""

        self.state = None
        self.results = {}

    def read(self, path):
        """Read atoms, parameter and calculated properties from file.

        This method must set self.state, the parameter dictionary
        self.parameters and calculated properties self.results like
        energy and forces."""

        pass

    def get_atoms(self):
        if self.state is None:
            raise ValueError('Calculator has no atoms')
        atoms = self.state.copy()
        atoms.calc = self
        return atoms

    @classmethod
    def read_atoms(cls, path, **kwargs):
        return cls(path, **kwargs).get_atoms()

    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2, ...).
        
        The special keyword 'parameters' ...

        A dictionary containing the parameters that have been changed
        is returned.  Subclasses may use this information to do stuff.
        """

        if 'parameters' in kwargs:
            filename = kwargs.pop('parameters')
            parameters = read_parameters_from_file(filename)
            parameters.update(kwargs)
            kwargs = parameters

        changed_parameters = {}

        for key, value in kwargs.items():
            if (key not in self.parameters or
                not equal(value, self.parameters[key])):
                changed_parameters[key] = value
                self.parameters[key] = value

        if changed_parameters:
            self.reset(changed_parameters)

        return changed_parameters

    def check_state(self, atoms):
        if self.state is None:
            changes = ['positions', 'numbers', 'cell', 'pbc']
        else:
            changes = []
            if not equal(self.state.positions, atoms.positions):
                changes.append('positions')
            if not equal(self.state.numbers, atoms.numbers):
                changes.append('numbers')
            if not equal(self.state.cell, atoms.cell):
                changes.append('cell')
            if not equal(self.state.pbc, atoms.pbc):
                changes.append('pbc')

        return changes

    def get_potential_energy(self, atoms, force_consistent=False):
        energy = self.get_property('energy', atoms)
        if force_consistent:
            return self.results.get('free_energy', energy)
        else:
            return energy

    def get_forces(self, atoms):
        return self.get_property('forces', atoms).copy()

    def get_stress(self, atoms):
        return self.get_property('stress', atoms).copy()

    def get_dipole_moment(self, atoms):
        return self.get_property('dipole', atoms).copy()

    def get_magnetic_moment(self, atoms):
        return self.get_property('magmom', atoms)

    def get_magnetic_moments(self, atoms):
        return self.get_property('magmoms', atoms).copy()

    def get_property(self, name, atoms):
        if name in self.notimplemented:
            raise NotImplementedError

        changes = self.check_state(atoms)
        if changes:
            self.reset()

        if name not in self.results:
            self._calculate(atoms, [name], changes)
        return self.results[name]

    def calculation_required(self, atoms, properties):
        changes = self.check_state(atoms)
        if changes:
            return True
        for name in properties:
            if name not in self.results:
                return True
        return False
        
    def _calculate(self, atoms, properties, changes):
        """Call hooks before and after actual calculation."""
        self.call_hooks('before')
        self.calculate(atoms, properties, changes)
        self.state = atoms.copy()
        self.call_hooks('after')

    def calculate(self, atoms, properties, changes):
        """Do the calculation.

        atoms: Atoms object
            Contains positions, unit-cell, ...
        properties: list of str
            List of what needs to be calculated can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'magmom' and
            'magmoms'.
        changes: list of str
            List of what has changed since last calculation.  Can be
            any of these four: 'positons', 'numbers', 'cell' and
            'pbc'.

        Subclasses need to implement this, but can ignore properties
        and changes if they want.
        """

        self.results = {'energy': 0.0,
                        'forces': np.zeros((len(atoms), 3)),
                        'stress': np.zeros(6),
                        'dipole': np.zeros(3),
                        'magmom': 0.0,
                        'magmoms': np.zeros(len(atoms))}
                        
    def add_hook(self, name, function, *args, **kwargs):
        self.hooks[name].append((function, args, kwargs))

    def call_hooks(self, name):
        for function, args, kwargs in self.hooks[name]:
            function(*args, **kwargs)


class FileIOCalculator(Calculator):
    """Base class for calculators that write input files and read output files.

    """

    def _calculate(self, atoms, properties=None, changes=None):
        self.write_input(atoms, properties)
        Calculator._calculate(self, atoms, properties, changes)
        self.read()

    def get_command(self):
        command = self.parameters.get('command')
        if command is None:
            name = 'ASE_' + self.name.upper() + '_COMMAND'
            command = os.env.get(name)
        if command is None:
            raise NotAvailable('Please set $%s environment variable ' % name +
                               'or supply the command keyword')
        return command

    def calculate(self, properties=None, changes=None):
        dir, label = self.split_path()
        command = self.get_command().replace('LABEL', label)
        olddir = os.getcwd()
        try:
            os.chdir(dir)
            errorcode = subprocess.call(command)
        finally:
            os.chdir(olddir)
        
        if errorcode:
            raise RuntimeError(errorcode)

    def split_path(self):
        """Convert path to directory and label.

        """
        return os.path.split(self.path)

    def write_input(self, properties=None):
        dir, label = self.split_path()
        if dir != os.curdir and not os.path.isdir(dir):
            os.makedirs(dir)
