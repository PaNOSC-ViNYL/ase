from __future__ import division, print_function
import numpy as np
from scipy.spatial import ConvexHull, Delaunay

from ase.utils import hill


class PhaseSpace:
    def __init__(self, references, verbose=True):
        """Phase-space.
        
        Example:
            
        >>> ps = PhaseSpace([({'Cu': 1}, -3.5),
        ...                  ({'Ni': 1}, -4.4),
        ...                  ({'Cu': 1, 'Ni': 1}, -8.1)])
        Species: Ni, Cu
        References: 3
        Simplices: 2
        >>> ps.find(Cu=2, Ni=1)
        reference         fraction         energy
        -----------------------------------------
        Cu              1.0000/  1         -3.500
        CuNi            2.0000/  2         -8.100
        -----------------------------------------
        Total energy:                     -11.600
        (-11.6, array([0, 2], dtype=int32), array([ 1.,  1.]))

        """

        self.verbose = verbose
        
        self.species = {}
        self.nspecies = 0
        self.systems = []
        for count, energy in references:
            natoms = 0
            for symbol, n in count.items():
                natoms += n
                if symbol not in self.species:
                    self.species[symbol] = self.nspecies
                    self.nspecies += 1
            self.systems.append((count, energy, natoms))
        
        if verbose:
            print('Species:', ', '.join(self.species))
            print('References:', len(references))
            
        self.points = np.zeros((len(self.systems), self.nspecies + 1))
        for s, (count, energy, natoms) in enumerate(self.systems):
            for symbol, n in count.items():
                self.points[s, self.species[symbol]] = n / natoms
            self.points[s, -1] = energy / natoms
        
        hull = ConvexHull(self.points[:, 1:])
        
        # Find relevant vertices:
        ok = hull.equations[:, -2] < 0
        vertices = set()
        for simplex in hull.simplices[ok]:
            vertices.update(simplex)
        self.vertices = np.array(list(vertices))
        
        if verbose:
            print('Simplices:', ok.sum())
        
        # Create triangulation:
        if self.nspecies == 2:
            D = Delaunay1D  # scipy's Delaunay doesn't like 1-d!
        else:
            D = Delaunay
        self.tri = D(self.points[self.vertices, 1:-1])
        
    def plot(self):
        pass
        
    def find(self, **kwargs):
        """Find the combination of the references with the lowest energy.
        
        Example::
            
            ps = PhaseSpace(...)
            ps.find(Cu=2, Ni=1)
            
        Returns energy, indices of references and coefficients."""
        
        point = np.zeros(self.nspecies)
        natoms = 0
        for symbol, n in kwargs.items():
            point[self.species[symbol]] = n
            natoms += n
        i = self.tri.find_simplex(point[1:] / natoms)
        indices = self.vertices[self.tri.simplices[i]]
        points = self.points[indices]
        scaledcoefs = np.linalg.solve(points[:, :-1].T, point)
        energy = np.dot(scaledcoefs, points[:, -1])
        
        if self.verbose:
            print('reference         fraction         energy')
            print('-----------------------------------------')

        coefs = []
        for coef, s in zip(scaledcoefs, indices):
            count, e, natoms = self.systems[s]
            coef /= natoms
            coefs.append(coef)
            if self.verbose:
                print('{0:15}{1:7.4f}/{2:3}{3:15.3f}'.format(hill(count),
                                                             coef * natoms,
                                                             natoms,
                                                             coef * e))
        if self.verbose:
            print('-----------------------------------------')
            print('Total energy: {0:27.3f}'.format(energy))
            
        return energy, indices, np.array(coefs)
        
        
class Delaunay1D:
    """Simple 1-d implementation."""
    def __init__(self, points):
        self.points = points[:, 0]
        a = self.points.argsort()
        self.simplices = np.array([a[:-1], a[1:]]).T

    def find_simplex(self, point):
        p = point[0]
        for i, s in enumerate(self.simplices[:, 1]):
            if p < self.points[s]:
                return i
        return i + 1
