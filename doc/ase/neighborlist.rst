.. module:: ase.neighborlist
.. module:: ase

Building neighbor-lists
=======================

A neighbor list is a collision detector for spheres: Given
a number of spheres of different radius located at different points,
it calculates the pairs of spheres that overlap.

ASE provides two implementations of neighbor lists.  The newer
linearly-scaling function
:func:`~ase.neighborlist.neighbor_list` and
the older quadratically-scaling class
:class:`~ase.neighborlist.PrimitiveNeighborList`.  The latter will likely
use the former as a backend in the future for linear scaling.

For flexibility, both implementations provide a “primitive”
interface which accepts arrays as arguments rather than the
more complex :class:`~ase.atoms.Atoms` objects.

Both implementations can be used via the :class:`~ase.neighborlist.NeighborList`
class. It also provides easy access to the two implementations methods and functions:

.. autoclass:: ase.neighborlist.NeighborList
   :members:

.. autoclass:: ase.neighborlist.PrimitiveNeighborList
   :members:

.. autoclass:: ase.neighborlist.NewPrimitiveNeighborList
   :members:

.. autofunction:: ase.neighborlist.neighbor_list

.. autofunction:: ase.neighborlist.primitive_neighbor_list

.. automethod:: ase.neighborlist.get_connectivity_matrix

.. _GPAW: http://wiki.fysik.dtu.dk/gpaw
