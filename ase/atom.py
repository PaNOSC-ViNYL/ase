"""This module defines the Atom object."""

import numpy as np

from ase.data import atomic_numbers, chemical_symbols, atomic_masses


def atomproperty(arrayname, init=None, default=None, doc=None):
    """Helper functions to easily get/set properties on the Atoms object."""

    def getter(atom):
        a = atom._atoms.arrays
        if arrayname in a:
            return a[arrayname][atom.index]
        else:
            return default(atom)  # e.g.: atomic_masses[self.number]

    def setter(atom, value):
        a = atom._atoms.arrays
        if arrayname not in a:
            init(atom, value)  # e.g.: set_tags(0) to initialize array
        a[arrayname][atom.index] = value

    return property(getter, setter, doc)


def abcproperty(index):
    """Helper function to easily create Atom ABC-property."""

    def getter(atom):
        spos = atom._atoms.get_scaled_positions()
        return spos[atom.index][index]

    def setter(atom, value):
        spos = atom._atoms.get_scaled_positions()
        spos[atom.index][index] = value
        atom._atoms.set_scaled_positions(spos)

    return property(getter, setter, doc='ABC'[index] + '-coordinate')


def xyzproperty(index):
    """Helper function to easily create Atom XYZ-property."""

    def getter(atom):
        return atom.position[index]

    def setter(atom, value):
        atom.position[index] = value

    return property(getter, setter, doc='XYZ'[index] + '-coordinate')


class Atom(object):
    """Class for representing a single atom.

    Parameters:

    symbol: str or int
        Can be a chemical symbol (str) or an atomic number (int).
    position: sequence of 3 floats
        Atomi position.
    tag: int
        Special purpose tag.
    momentum: sequence of 3 floats
        Momentum for atom.
    mass: float
        Atomic mass in atomic units.
    magmom: float or 3 floats
        Magnetic moment.
    charge: float
        Atomic charge.
    """
    attrnames = {'position': 'positions',
                 'number': 'numbers',
                 'tag': 'tags',
                 'momentum': 'momenta',
                 'mass': 'masses',
                 'magmom': 'magmoms',
                 'charge': 'charges'}

    def __init__(self, symbol='X', position=None,
                 tag=None, momentum=None, mass=None,
                 magmom=None, charge=None,
                 atoms=None, index=None):
        if atoms is None:
            from ase import Atoms
            def get(obj):
                return None if obj is None else [obj]
            atoms = Atoms([symbol],
                          positions=get(position),
                          tags=get(tag),
                          momenta=get(momentum),
                          masses=get(mass),
                          magmoms=get(magmom),
                          charges=get(charge))
            index = 0
        self._atoms = atoms
        self.index = index

    def has(self, name):
        atomsname = self.attrnames[name]
        return self._atoms.has(atomsname)

    def __eq__(self, other):
        return (isinstance(other, Atom)
                and self.index == other.index
                and self._atoms == other._atoms)

    def __neq__(self, other):
        return not self == other

    def __repr__(self):
        tokens = []
        for name in ['tag', 'momentum', 'mass', 'magmom', 'charge']:
            if self.has(name):
                val = getattr(self, name)
                if isinstance(val, np.ndarray):
                    val = val.tolist()
                token = '{}={!r}'.format(name, val)
                tokens.append(token)

        tokens.append('index={}'.format(self.index))
        return 'Atom({!r}, {})'.format(self.symbol, ', '.join(tokens))

    def cut_reference_to_atoms(self):
        """Cut reference to atoms object."""
        self._atoms = self._atoms[self.index:self.index + 1]
        self.index = 0

    @property
    def symbol(self):
        return chemical_symbols[self.number]

    @symbol.setter
    def symbol(self, value):
        number = atomic_numbers[value]
        self.number = number

    number = atomproperty('numbers', doc='Atomic number')
    position = atomproperty('positions', doc='XYZ-coordinates')
    tag = atomproperty('tags',
                       default=lambda self: 0,
                       init=lambda self, val: self._atoms.set_tags(0),
                       doc='Integer tag')
    momentum = atomproperty('momenta',
                            default=lambda self: np.zeros(3),
                            init=lambda self, val: self._atoms.set_momenta(
                                self._atoms.get_momenta()),
                            doc='XYZ-momentum')
    mass = atomproperty('masses',
                        default=lambda self: atomic_masses[self.number],
                        init=lambda self, val: self._atoms.set_masses(),
                        doc='Atomic mass')
    magmom = atomproperty('magmoms',
                          default=lambda self: 0.0,
                          # Shape can be len(atoms), 3 x len(atoms), ...
                          init=lambda self, val:
                          self._atoms.set_initial_magnetic_moments(
                              np.zeros((len(self._atoms),) + np.shape(val))),
                          doc='Initial magnetic moment')
    charge = atomproperty('charges',
                          default=lambda self: 0.0,
                          init=lambda self, val:
                          self._atoms.set_initial_charges(
                              np.zeros(len(self._atoms))),
                          doc='Atomic charge')
    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)
    scaled_position = abcproperty(slice(None, None, None))
    a = abcproperty(0)
    b = abcproperty(1)
    c = abcproperty(2)
