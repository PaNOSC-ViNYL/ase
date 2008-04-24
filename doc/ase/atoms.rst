The Atoms object
================

.. automodule:: ase.atoms

.. autoclass:: Atoms

Like with a single :mod: atom the properties of a collection of atoms
can be accessed and changed with "get-" and "set-"methods. For example
the positions can be controlled with set_position



 Atom class the properties of the

``cell``: Unit cell size
  This can be a sequence of three numbers for
  an orthorhombic unit cell or three by three numbers for a general
  unit cell (a sequence of three sequences of three numbers).  The
  default value is ``(1, 1, 1)``.

``periodic``: Boundary conditions
  The default value is ``False`` - a value of ``True`` would give
  periodic boundary conditions along all three axes.  It is possible
  to give a sequence of three booleans to specify periodicity along
  specific axes.


Here is how you could make an infinite gold wire with a bond length of
2.9:

>>> from ASE import Atom, ListOfAtoms
>>> d = 2.9
>>> L = 10.0
>>> wire = ListOfAtoms([Atom('Au', (L / 2, L / 2, 0))],
...                    cell=(L, L, d), periodic=(0, 0, 1))


You can also use the following methods to work with the unit cell and the
boundary conditions:

``GetUnitCell()``:
  Returns a three by three array.

``SetUnitCell(cell, fix=False)``:
  Change the size of the unit cell.  If the optional argument ``fix``
  is ``True`` (defaults to ``False``), then the positons of the atoms
  are fixed, otherwise the atoms are moved so that their positions
  relative to the unit cell are kept.

``GetBoundaryConditions()``:
  Returns a tuple of three booleans.

Here is how you would do bulk ruthenium (hcp):

>>> from math import sqrt
>>> a = 2.70
>>> c = 1.59 * a
>>> bulk = ListOfAtoms([Atom('Ru', (0, 0, 0)),
...                     Atom('Ru', (1 / 3., 1 / 3., 1 / 2.))],
...                    periodic=True)
>>> bulk.SetUnitCell([(a, 0, 0),
...                   (a / 2, a * sqrt(3) / 2, 0),
...                   (0, 0, c)])


In addition, an ``ListOfAtoms`` instance has the following methods:

``GetKineticEnergy()``:
   Returns the total kinetic energy.

``Copy()``:
   Return a fresh copy.  Everything but a possibly attached
   calculator is copied (next section_ will explain how to attach a
   calculator).

``Repeat(repeat)``:
   The argument ``repeat`` is a sequence of three
   positive integers (``n1``, ``n2``, ``n3`` - one for each axis), and
   a copy of the ``ListOfAtoms`` object repeated ``n1 * n2 * n3``
   times is returned.

.. _section: calculators.html:Calculators


.. warning::
   Using the same atom object in several ListOfAtoms objects is not
   allowed!  Use a copy instead:

   >>> loa1.append(loa2[7])
   RuntimeError: Atom belongs to another ListOfAtoms!
   >>> loa1.append(loa2[7].Copy())



Array methods
-------------

It is possible to work with the properties of atoms in a ListOfAtoms by
using the methods defined for the individual atoms like this:

>>> m = ListOfAtoms([Atom('O', (0, 0, 0)),
...                  Atom('H', (0.773, 0.600, 0)),
...                  Atom('H', (-0.773, 0.600, 0))])
>>> m[1].GetCartesianPosition()
array([ 0.773,  0.6  ,  0.   ])
>>> m[0].SetCartesianPosition((0, -0.1, 0))

However, a ListOfAtoms that conforms to the ASE specification must have a
set of array methods to do the same for all atoms in one step:

>>> m.GetCartesianPositions()
array([[ 0.   ,  0.   ,  0.   ],
       [ 0.773,  0.6  ,  0.   ],
       [-0.773,  0.6  ,  0.   ]])
>>> m.SetCartesianPositions([(0, -0.1, 0),
...                          (0.773, 0.600, 0),
...                          (-0.773, 0.600, 0)])

The following methods must be defined:

===========================  ===========================  =========  ==========
Get methods                  Set methods                  type       shape
===========================  ===========================  =========  ==========
``GetCartesianPositions``    ``SetCartesianPositions``    ``Float``  ``(n, 3)``
``GetCartesianMomenta``      ``SetCartesianMomenta``      ``Float``  ``(n, 3)``
``GetCartesianVelocities``   ``SetCartesianVelocities``   ``Float``  ``(n, 3)``
``GetCartesianForces``       ``SetCartesianForces``       ``Float``  ``(n, 3)``
``GetTags``                  ``SetTags``                  ``Int``    ``(n,)``
``GetAtomicNumbers``         ``SetAtomicNumbers``         ``Int``    ``(n,)``
``GetMasses``                ``SetMasses``                ``Float``  ``(n,)``
``GetKineticEnergies``                                    ``Float``  ``(n,)``
``GetGeneralizedPositions``  ``SetGeneralizedPositions``  ``Float``  ``(n*3,)``
``GetGeneralizedForces``     ``SetGeneralizedForces``     ``Float``  ``(n*3,)``
===========================  ===========================  =========  ==========


The ``Get`` methods will return Numeric_ arrays of the given shape (``n``
is the number of atoms) and type, and the ``Set`` methods will take
anything with the correct shape.

.. _Numeric: http://numpy.sf.net
