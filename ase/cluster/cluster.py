import os
import math
import random
import numpy as np

import cPickle as pickle
import new

from ase import Atoms
from ase.data import chemical_symbols
from ase.cluster.base import ClusterBase
from ase.cluster.clusteratom import ClusterAtom

class Cluster(Atoms, ClusterBase):
    _datasyn = {'numbers':       ('number',       int,   ()  ),
                'positions':     ('position',     float, (3,)),
                'tags':          ('tag',          int,   ()  ),
                'momenta':       ('momentum',     float, (3,)),
                'masses':        ('mass',         float, ()  ),
                'magmoms':       ('magmom',       float, ()  ),
                'charges':       ('charge',       float, ()  ),
               }

    def __repr__(self):
        output = ('Cluster(symbols=%i%s, latticeconstant=%.2f' % 
                  (len(self), chemical_symbols[self.atomic_number], self.lattice_constant))

        for name in self.arrays:
            output += ', %s=...' % name

        output += ', cell=%s' % self._cell.diagonal().tolist()

        if self.get_center() is not None:
            output += ', center=%s' % self.get_center().tolist()

        return output + ')'

    def __getitem__(self, i):
        c = ClusterAtom(atoms=self, index=i)
        return c

    def __setitem__(self, i, atom):
        #raise Warning('Use direct assignment like atoms[i].type = x!')

        #If implemented make sure that all values are cleared before copied
        if not isinstance(atom, ClusterAtom):
            raise Warning('The added atom is not a ClusterAtom instance!')

        for name in self.arrays.keys():
            singular, dtype, shape = self._datasyn[name]
            self[i]._set(singular, np.zeros(shape, dtype))

        self[i]._set('number', atom._get('number', True))

        for name in atom._data:
            self[i]._set(name, atom._get(name, True))

    def append(self, atom):
        if not isinstance(atom, ClusterAtom):
            raise Warning('The added atom is not a ClusterAtom instance!')

        n = len(self)

        for name, a in self.arrays.items():
            b = np.zeros((n + 1,) + a.shape[1:], a.dtype)
            b[:n] = a
            self.arrays[name] = b

        for name in atom._data:
            self[-1]._set(name, atom._get(name, True))

    def extend(self, atoms):
        if not isinstance(atoms, Cluster):
            raise Warning('The added atoms is not in a Cluster instance!')

        for atom in atoms:
            self.append(atom)

    def copy(self):
        cluster = Cluster(symbol=self.atomic_number,
                          latticeconstant=self.lattice_constant,
                          symmetry=self.symmetry,
                          cell=self.get_cell(),
                          center=self.get_center())

        for name, a in self.arrays.items():
            cluster.arrays[name] = a.copy()

        return cluster

    def get_layers(self):
        layers = []

        for s in self.surfaces:
            n = self.miller_to_direction(s)
            r = np.dot(self.get_positions() - self.center, n).max()
            d = self.get_layer_distance(s, 2)
            l = 2 * np.round(r / d).astype(int)

            ls = np.arange(l-1,l+2)
            ds = np.array([self.get_layer_distance(s, i) for i in ls])

            mask = (np.abs(ds - r) < 1e-10)

            layers.append(ls[mask][0])

        return np.array(layers, int)

    def get_diameter(self, method='Volume'):
        """Makes an estimate of the cluster diameter based on two different
        methods.

        method = 'Volume': Returns the diameter of a sphere with the 
                           same volume as the atoms. (Default)
        
        method = 'Shape': Returns the distance between two opposit surfaces
                          averaged over all the surfaces on the cluster.
        """

        if method == 'Shape':
            directions = np.array(lattice[self.symmetry]['surface_names'])
            pos = self.get_positions()

            d = 0.0
            for n in directions:
                n = np.dot(n, self.get_basis())
                n = n / np.linalg.norm(n)
                r = np.dot(pos, n)
                d += r.max() - r.min()

            return d / len(directions)
        elif method == 'Volume':
            return (3.0/2.0 * len(self)/math.pi) ** (1.0/3.0) * self.lattice_constant
        else:
            return 0.0

    #Functions to acces the properties
    def get_center(self):
        return self._center.copy()

    def set_center(self, center):
        self._center = np.array(center, float)

    def get_lattice_basis(self):
        return self._lattice_basis.copy()

    def get_resiproc_basis(self):
        return self._resiproc_basis.copy()

    def get_atomic_basis(self):
        return self._atomic_basis.copy()

    #Functions to store the cluster
    def write(self, filename=None):
        if not isinstance(filename, str):
            raise Warning('You must specify a valid filename.')

        if os.path.isfile(filename):
            os.rename(filename, filename + '.bak')

        d = {'symbol': self.atomic_number,
             'latticeconstant': self.lattice_constant,
             'symmetry': self.symmetry,
             'multiplicity': self.multiplicity,
             'center': self.get_center(),
             'cell': self.get_cell(),
             'pbc': self.get_pbc()}

        f = open(filename, 'w')
        f.write('Cluster')
        pickle.dump(d, f)
        pickle.dump(self.arrays, f)
        f.close()

    def read(self, filename):
        if not os.path.isfile(filename):
            raise Warning('The file specified do not exist.')

        f = open(filename, 'r')

        try:
            if f.read(len('Cluster')) != 'Cluster':
                raise Warning('This is not a compatible file.')
            d = pickle.load(f)
            self.arrays = pickle.load(f)
        except EOFError:
            raise Warinig('Bad file.')

        f.close()

        if 'multiplicity' in d:
            self.multiplicity = d['multiplicity']
        else:
            self.multiplicity = 1

        self.atomic_number = d['symbol']
        self.lattice_constant = d['latticeconstant']
        self.symmetry = d['symmetry']
        self.set_center(d['center'])
        self.set_cell(d['cell'])
        self.set_pbc(d['pbc'])
        self.set_constraint()
        self.adsorbate_info = {}
        self.calc = None

