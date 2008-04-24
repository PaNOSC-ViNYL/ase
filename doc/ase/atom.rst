The Atom object
===============

ASE defines a python class called ``Atom`` to setup and handle atoms
in electronic structure and molecular simulations. From a python
script, atoms can be created like this:

>>> from ase import Atom
>>> a1 = Atom('Si', (0, 0, 0))
>>> a2 = Atom('H', (1.3, 0, 0), mass=2)
>>> a3 = Atom(position=(0, 0, 0), Z=14)  # same is a1

The first argument to the constructor of an ``Atom`` object is the
chemical symbol, and the second argument is the position.  The
position can be any numerical sequence of length three.  The
properties of an atom can also be set using keywords like it is done
in the ``a2`` example.

The different properties of an atom can generally be obtained with a
"get-method" and changed with a "set-method". For example for the position of the atom:

>>> a1.set_position([1,0,0])
>>> a1.get_position()
[1,0,0]

The full definition of the parameters for the atom object is as follows:

   .. autoclass:: ase.atom.Atom

.. seealso::

   :mod:`atoms` :
   More information about how to use atoms

   :mod:`calculators`
   Information about how to calculate forces and energies of atoms.

