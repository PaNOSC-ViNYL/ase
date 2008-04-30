.. module:: constraints
   :synopsis: Constraining some degrees of freedom

===========
Constraints
===========


When performing minimizations or dynamics one may wish to keep some
degrees of freedom in the system fixed. One way of doing this is by
attaching constraint object(s) directly to the atoms object.


The FixAtoms class
==================

This class is used for fixing some of the atoms.

.. class:: FixAtoms(indices=None, mask=None)



**XXX positive or negative mask???**



You must supply either the indices of the atoms that should be fixed
or a mask. The mask is a list of booleans, one for each atom, being true
if the atoms should be kept fixed.

Example of use::

  >>> from ase import FixAtoms
  >>> c = FixAtoms(mask=[s == 'Cu' for s in atoms.get_chemical_symbols()])
  >>> atoms.set_constraint(c)

This will fix the positions of all the Cu atoms in all the
dynamics/minimazations.


The FixBondLength class
=======================

This class is used to fix the distance between two atoms specified by
their indices (a1 and a2)

.. class:: FixBondLength(a1, a2)

Example of use::

  >>> from ase import FixBondLength
  >>> c = FixBondLength(0, 1)
  >>> atoms.set_constraint(c)

In this example the distance between the distance between the atoms
with indices 0 and 1 will be fixed in all following dynamics and/or
minimizations performed on the atoms object.

Combining constraints
=====================

It is possible to supply several constraints on an atoms object. For
example one may wish to keep the distance between two nitrogen atoms
fixed while relaxing it on a fixed ruthenium surface::

  >>> from ase import Atoms,FixAtoms,FixBondLength
  >>> pos = [[0.00000, 0.00000,  9.17625],
  ...        [0.00000, 0.00000, 10.27625],
  ...        [1.37715, 0.79510,  5.00000],
  ...        [0.00000, 3.18039,  5.00000],
  ...        [0.00000, 0.00000,  7.17625],
  ...        [1.37715, 2.38529,  7.17625]]
  >>> unitcell = [5.5086, 4.7706, 15.27625]

  >>> atoms = Atoms(positions=pos,
  ...               symbols='N2Ru4',
  ...               cell=unitcell,
  ...               pbc=[True,True,False])

  >>> c_fa = FixAtoms(mask=[s == 'Ru' for s in atoms.get_chemical_symbols()])
  >>> c_fb = FixBondLength(0, 1)
  >>> atoms.set_constraint([c_fa, c_fb])

When applying more than one constraint they are passed as a list in
the set_constraint method.

The Filter class
================

Constraints can also be applied via filters, which acts as a wrapper
around an atoms object. A typical use case will look like this::

   -------       --------       ----------
  |       |     |        |     |          |
  | Atoms |<----| Filter |<----| Dynamics |
  |       |     |        |     |          |
   -------       --------       ----------

and in Python this would be::

  >>> atoms = Atoms(...)
  >>> filter = Filter(atoms, ...)
  >>> dyn = Dynamics(filter, ...)

Currently only one filter is implemented. See the description of
:class:`~constraints.Filter` below.


This class hides some of the atoms in an Atoms object.

.. class:: Filter(atoms, indices=None, mask=None)

You must supply either the indices of the atoms that should be kept
visible or a mask. The mask is a list of booleans, one for each atom,
being true if the atom should be kept visible.

Example of use::

  >>> from ase import Atoms, Filter
  >>> atoms=Atoms(positions=[[ 0    , 0    , 0],
  ...                        [ 0.773, 0.600, 0],
  ...                        [-0.773, 0.600, 0]],
  ...             symbols='OH2')
  >>> f1 = Filter(atoms, indices=[1, 2])
  >>> f2 = Filter(atoms, mask=[0, 1, 1])
  >>> f3 = Filter(atoms, mask=[s == 'H' for s in atoms.get_chemical_symbols()])
  >>> f1.get_positions()
  [[ 0.773  0.6    0.   ]
   [-0.773  0.6    0.   ]]

In all three filters (f1, f2 and f3) only the hydrogen atoms are made
visible. When asking for the positions only the positions of the
hydrogen atoms are returned.

