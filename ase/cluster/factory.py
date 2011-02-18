import numpy as np

from ase import Atoms
from ase.data import atomic_numbers
from ase.lattice.spacegroup import Spacegroup
from ase.cluster.base import ClusterBase
from ase.cluster.cluster import Cluster

class ClusterFactory(ClusterBase):
    directions = [[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]]

    atomic_basis = np.array([[0., 0., 0.]])

    size_factor = 1

    def __call__(self, symbol, surfaces, layers, latticeconstant=None, vacuum=0.0, debug=0):
        self.debug = debug

        #Find the atomic number
        if symbol is not None:
            if isinstance(symbol, str):
                self.atomic_number = atomic_numbers[symbol]
            else:
                self.atomic_number = symbol
        else:
            raise Warning('You must specify a atomic symbol or number!')

        self.set_lattice_constant(latticeconstant)
        self.set_basis()

        if self.debug:
            print "Lattice constant(s):", self.lattice_constant
            print "Lattice basis:", self.lattice_basis
            print "Resiprocal basis:", self.resiproc_basis
            print "Atomic basis:", self.atomic_basis

        self.set_surfaces_layers(surfaces, layers)
        self.set_lattice_size()

        if self.debug:
            print "Center position:", self.center.round(2)
            print "Base lattice size:", self.size

        cluster = self.make_cluster(vacuum)
        cluster.symmetry = self.xtal_name
        cluster.center = self.center.copy()
        cluster.surfaces = self.surfaces.copy()
        cluster.lattice_basis = self.lattice_basis.copy()
        cluster.atomic_basis = self.atomic_basis.copy()
        cluster.resiproc_basis = self.resiproc_basis.copy()
        return cluster

    def make_cluster(self, vacuum):
        # Make the base crystal by repeating the unit cell
        size = np.array(self.size)
        translations = np.zeros((size.prod(), 3))
        for h in range(size[0]):
            for k in range(size[1]):
                for l in range(size[2]):
                    i = h * (size[1] * size[2]) + k * size[2] + l
                    translations[i] = np.dot([h, k, l], self.lattice_basis)

        atomic_basis = np.dot(self.atomic_basis, self.lattice_basis)
        positions = np.zeros((len(translations) * len(atomic_basis), 3))
        n = len(atomic_basis)
        for i, trans in enumerate(translations):
            positions[n*i:n*(i+1)] = atomic_basis + trans

        # Remove all atoms that is outside the defined surfaces
        for s, l in zip(self.surfaces, self.layers):
            n = self.miller_to_direction(s)
            rmax = self.get_layer_distance(s, l + 0.5)

            r = np.dot(positions - self.center, n)
            mask = np.less(r, rmax)

            if self.debug > 1:
                print "Cutting %s at %i layers ~ %.3f A" % (s, l, rmax)

            positions = positions[mask]

        # Fit the cell, so it only just consist the atoms
        min = np.zeros(3)
        max = np.zeros(3)
        for i in range(3):
            v = self.directions[i]
            r = np.dot(positions, v)
            min[i] = r.min()
            max[i] = r.max()

        cell = max - min + vacuum
        positions = positions - min + vacuum / 2.0
        self.center = self.center - min + vacuum / 2.0

        return Cluster(symbols=[self.atomic_number] * len(positions),
                       positions=positions, cell=cell)

    def set_lattice_size(self):
        max = np.ones(3, int) * 100
        min = np.ones(3, int) * 100
        
        for i in range(3):
            v = self.lattice_basis[i]
            for s, l in zip(self.surfaces, self.layers):
                if l < 0: continue

                n = self.miller_to_direction(s)
                nl = self.get_layer_distance(s, l)
                n = nl * n

                cos = np.dot(n, v)
                if np.abs(cos) < 1e-10:
                    d = 0.0
                else:
                    d = nl**2 / cos

                if self.debug > 1:
                    print "%s dot %s = %.3f" % (s, v.round(2), d)

                if d > 1e-10:
                    d = np.int(np.ceil(d))
                    if d < max[i]: max[i] = d
                elif d < -1e-10:
                    d = np.int(np.ceil(-d))
                    if d < min[i]: min[i] = d

        self.center = np.dot(min, self.lattice_basis) * self.size_factor
        self.size = (min + max + np.ones(3, int)) * self.size_factor

    def set_surfaces_layers(self, surfaces, layers):
        if len(surfaces) != len(layers):
            raise ValueError("Improper size of surface and layer arrays: %i != %i"
                             % (len(surfaces), len(layers)))

        sg = Spacegroup(self.spacegroup)
        surfaces = np.array(surfaces)
        layers = np.array(layers)

        for i, s in enumerate(surfaces):
            s = reduce_miller(s)
            surfaces[i] = s

        surfaces_full = surfaces.copy()
        layers_full = layers.copy()

        for s, l in zip(surfaces, layers):
            equivalent_surfaces = sg.equivalent_reflections(s.reshape(-1, 3))

            for es in equivalent_surfaces:
                # If the equivalent surface (es) is not in the surface list,
                # then append it.
                if not np.equal(es, surfaces_full).all(axis=1).any():
                    surfaces_full = np.append(surfaces_full, es.reshape(1, 3), axis=0)
                    layers_full = np.append(layers_full, l)

        self.surfaces = surfaces_full.copy()
        self.layers = layers_full.copy()

    def get_resiproc_basis(self, basis):
        """Returns the resiprocal basis to a given lattice (crystal) basis"""
        k = 1 / np.dot(basis[0], cross(basis[1], basis[2]))

        # The same as the inversed basis matrix transposed
        return k * np.array([cross(basis[1], basis[2]),
                             cross(basis[2], basis[0]),
                             cross(basis[0], basis[1])])

# Helping functions
def cross(a, b):
    """The cross product of two vectors."""
    return np.array([a[1]*b[2] - b[1]*a[2],
                     a[2]*b[0] - b[2]*a[0],
                     a[0]*b[1] - b[0]*a[1]])

def GCD(a,b):
    """Greatest Common Divisor of a and b."""
    #print "--"
    while a != 0:
        #print a,b,">",
        a,b = b%a,a
        #print a,b
    return b

def reduce_miller(hkl):
    """Reduce Miller index to the lowest equivalent integers."""
    hkl = np.array(hkl)
    old = hkl.copy()

    d = GCD(GCD(hkl[0], hkl[1]), hkl[2])
    while d != 1:
        hkl = hkl / d
        d = GCD(GCD(hkl[0], hkl[1]), hkl[2])

    if np.dot(old, hkl) > 0:
        return hkl
    else:
        return -hkl

