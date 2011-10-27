"""
Function-like objects that creates cubic clusters.
"""

import numpy as np
from ase.cluster.factory import ClusterFactory
from ase.data import reference_states as _refstate

class SimpleCubicFactory(ClusterFactory):
    spacegroup = 221

    xtal_name = 'sc'

    def get_lattice_constant(self):
        "Get the lattice constant of an element with cubic crystal structure."
        symmetry = _refstate[self.atomic_numbers[0]]['symmetry']
        if symmetry != self.xtal_name:
            raise ValueError, ("Cannot guess the %s " % (self.xtal_name,) +
                               "lattice constant of an element with crystal " +
                               "structure %s." % (symmetry,))
        return _refstate[self.atomic_numbers[0]]['a']

    def set_basis(self):
        a = self.lattice_constant
        if not isinstance(a, (int, float)):
            raise ValueError("Improper lattice constant for %s crystal." % (xtal_name,))

        self.lattice_basis = np.array([[a, 0., 0.],
                                       [0., a, 0.],
                                       [0., 0., a]])

        self.resiproc_basis = self.get_resiproc_basis(self.lattice_basis)

SimpleCubic = SimpleCubicFactory()

class BodyCenteredCubicFactory(SimpleCubicFactory):
    xtal_name = 'bcc'

    atomic_basis = np.array([[0., 0., 0.],
                             [.5, .5, .5]])

BodyCenteredCubic = BodyCenteredCubicFactory()

class FaceCenteredCubicFactory(SimpleCubicFactory):
    xtal_name = 'fcc'

    atomic_basis = np.array([[0., 0., 0.],
                             [0., .5, .5],
                             [.5, 0., .5],
                             [.5, .5, 0.]])

FaceCenteredCubic = FaceCenteredCubicFactory()

def Octahedron(symbol, length, cutoff=0, latticeconstant=None):
    """
    Returns Face Centered Cubic clusters of the octahedral class depending
    on the choice of cutoff.

    Type                            Condition
    ----                            ---------
    Regular octahedron              cutoff = 0
    Truncated octahedron            cutoff > 0
    Regular truncated octahedron    length = 3 * cutoff + 1
    Cuboctahedron                   length = 2 * cutoff + 1

    Parameters
    ----------
    symbol: The chemical symbol (or atomic number) of the element.

    length: Number of atoms on the square edges of the complete octahedron.

    cutoff (optional): Number of layers cut at each vertex.

    latticeconstant (optional): The lattice constant. If not given,
    then it is extracted form ase.data.
    """

    # Check length and cutoff
    if length < 2:
        raise ValueError("The lenght must be greater than one.")

    if cutoff < 0 or length < 2 * cutoff + 1:
        raise ValueError("The cutoff must fullfill: > 0 and <= (length - 1) / 2.")

    # Create cluster
    surfaces = [(1,1,1), (1,0,0)]
    if length % 2 == 0:
        center = np.array([0.5, 0.5, 0.5])
        layers = [length/2, length - 1 - cutoff]
    else:
        center = np.array([0.0, 0.0, 0.0])
        layers = [(length - 1)/2, length - 1 - cutoff]

    return FaceCenteredCubic(symbol, surfaces, layers, latticeconstant, center)

