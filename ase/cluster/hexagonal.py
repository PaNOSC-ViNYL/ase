"""
Function-like objects that creates cubic clusters.
"""

import numpy as np
from ase.cluster.factory import ClusterFactory
from ase.lattice.hexagonal import HexagonalFactory as HexFactory, \
                                  HexagonalClosedPackedFactory as HCPFactory, \
                                  GraphiteFactory as GphFactory
from ase.data import reference_states as _refstate

class HexagonalFactory(ClusterFactory):
    spacegroup = 191

    xtal_name = 'hexagonal'

    size_factor = 2

    lattice_factory = HexFactory()

    def set_lattice_constant(self, latticeconstant):
        "Get the lattice constant of an element with cubic crystal structure."
        if latticeconstant is None:
            if _refstate[self.atomic_number]['symmetry'].lower() != self.xtal_name:
                raise ValueError, (("Cannot guess the %s lattice constant of"
                                    + " an element with crystal structure %s.")
                                   % (self.xtal_name,
                                      _refstate[self.atomic_number]['symmetry'].lower()))
            self.lattice_constant = _refstate[self.atomic_number].copy()
        else:
            self.lattice_constant = latticeconstant

    def set_basis(self):
        lattice = self.lattice_constant
        if isinstance(lattice, dict):
            a = lattice['a']
            try:
                c = lattice['c']
            except KeyError:
                c = a * lattice['c/a']
        else:
            if len(lattice) == 2:
                (a, c) = lattice
            else:
                raise ValueError("Improper lattice constants for hcp crystal.")
        
        self.lattice_constant = (a, c)

        self.lattice_basis = np.array([[a, 0, 0],
                                       [-a/2.0, a*np.sqrt(3.0)/2.0, 0],
                                       [0, 0, c]])
        self.resiproc_basis = self.get_resiproc_basis(self.lattice_basis)

    def set_surfaces_layers(self, surfaces, layers):
        for i, s in enumerate(surfaces):
            if len(s) == 4:
                (a, b, c, d) = s
                x = 4*a + 2*b
                y = 2*a + 4*b
                z = 3*d
                surfaces[i] = [x, y, z]

        ClusterFactory.set_surfaces_layers(self, surfaces, layers)

Hexagonal = HexagonalFactory()

class HexagonalClosedPackedFactory(HexagonalFactory):
    """A factory for creating HCP clusters."""
    spacegroup = 194

    xtal_name = 'hcp'

    atomic_basis = np.array([[1.0/3.0, 2.0/3.0, 0.5]])

    lattice_factory = HCPFactory()

HexagonalClosedPacked = HexagonalClosedPackedFactory()

class GraphiteFactory(HexagonalFactory):
    """A factory for creating graphite clusters."""
    xtal_name = "graphite"

    atomic_basis = np.array([[1.0/3.0, 2.0/3.0, 0], [1.0/3.0,2.0/3.0,0.5], [2.0/3.0,1.0/3.0,0.5]])

    lattice_factory = GphFactory()

Graphite = GraphiteFactory()

