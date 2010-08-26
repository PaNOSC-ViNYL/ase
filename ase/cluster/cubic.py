"""
Function-like objects that creates cubic clusters.
"""

import numpy as np
from ase.cluster.factory import ClusterFactory
from ase.lattice.cubic import SimpleCubicFactory as SCFactory, \
                              BodyCenteredCubicFactory as BCCFactory, \
                              FaceCenteredCubicFactory as FCCFactory, \
                              DiamondFactory as DiamondFactory
from ase.data import reference_states as _refstate

class SimpleCubicFactory(ClusterFactory):
    spacegroup = 221

    xtal_name = 'sc'

    int_lattice_basis = 1.0 * np.array([[1, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, 1]])

    lattice_factory = SCFactory()

    def set_lattice_constant(self, latticeconstant):
        "Get the lattice constant of an element with cubic crystal structure."
        if latticeconstant is None:
            if _refstate[self.atomic_number]['symmetry'].lower() != self.xtal_name:
                raise ValueError, (("Cannot guess the %s lattice constant of"
                                    + " an element with crystal structure %s.")
                                   % (self.xtal_name,
                                      _refstate[self.atomic_number]['symmetry'].lower()))
            self.lattice_constant = _refstate[self.atomic_number]['a']
        else:
            self.lattice_constant = latticeconstant

    def set_basis(self):
        a = self.lattice_constant
        if not isinstance(a, (int, float)):
            raise ValueError("Improper lattice constants for %s crystal." % (xtal_name,))

        self.lattice_basis = a * self.int_lattice_basis
        self.resiproc_basis = self.get_resiproc_basis(self.lattice_basis)

SimpleCubic = SimpleCubicFactory()

class BodyCenteredCubicFactory(SimpleCubicFactory):
    spacegroup = 229

    xtal_name = 'bcc'

    atomic_basis = 0.5 * np.array([[1, 1, 1]])

    lattice_factory = BCCFactory()

BodyCenteredCubic = BodyCenteredCubicFactory()

class FaceCenteredCubicFactory(SimpleCubicFactory):
    spacegroup = 225

    xtal_name = 'fcc'

    atomic_basis = 0.5 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])

    lattice_factory = FCCFactory()

FaceCenteredCubic = FaceCenteredCubicFactory()

