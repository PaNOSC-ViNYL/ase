import numpy as np
from ase.cluster.data import lattice
from ase.atom import Atom, data as olddata
from ase.data import atomic_numbers, chemical_symbols, reference_states

class ClusterAtom(Atom):
    """Cluster Atom"""
    _datasyn = olddata

    _data = {}

    def __init__(self, symbol='X', position=(0.0, 0.0, 0.0), atoms=None, index=None):
        self.atoms = atoms
        self.index = index

        if atoms is None:
            if isinstance(symbol, str):
                self.number = atomic_numbers[symbol]
            else:
                self.number = symbol
 
            self.position = np.array(position, float)

    def __repr__(self):
        output = 'ClusterAtom(%s, %s' % (self.symbol, self.position.tolist())

        for name, d in self._data.items():
            if name != 'number' and name != 'position' and d is not None:
                output += ', %s=%s' % (name, d)

        return output + ')'

    def _get(self, name, copy=False):
        if self.atoms is None:
            if not self.has(name): return None

            return self._data[name]
        else:
            plural = self._datasyn[name][0]
            if not self.atoms.has(plural): return None

            if copy:
                return self.atoms.arrays[plural][self.index].copy()
            else:
                return self.atoms.arrays[plural][self.index]

    def _get_copy(self, name): self._get(name, copy=True)

    def _set(self, name, value, copy=False):
        if self.atoms is None or copy is True:
            self._data[name] = value
        else:
            plural, dtype, shape = self._datasyn[name]
            if plural == 'neighbors':
                symmetry = reference_states[self.number]['symmetry'].lower()
                shape = (lattice[symmetry]['neighbor_count'],) # XXXX

            if self.atoms.has(plural):
                self.atoms.arrays[plural][self.index] = value
            else:
                array = np.zeros((len(self.atoms),) + shape, dtype)
                array[self.index] = value
                self.atoms.set_array(plural, array)

    def has(self, name):
        return name in self._data

    def cut_reference_to_atoms(self):
        for name, a in self.atoms.arrays.items():
            self._set(self.atoms._datasyn[name][0], a[self.index].copy(), True)
        self.atoms = None
        self.index = None

    def get_symbol(self): return chemical_symbols[self._get('number')]
    #def get_(self): return self._get('')

    def set_symbol(self, value): self._set('number', atomic_numbers[value])
    #def get_(self, value): return self._set('', value)

    symbol = property(get_symbol, set_symbol, doc='Chemical symbol')

