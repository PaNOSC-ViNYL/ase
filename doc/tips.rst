===============
Tips and tricks
===============

https://www.python.org/
https://docs.scipy.org/doc/numpy/
https://docs.scipy.org/doc/scipy/reference/


Atoms objects
=============

Species
-------

.. doctest::

    >>> from ase import Atoms
    >>> atoms = Atoms('CH4')
    >>> len(set(atoms.numbers))  # number of species
    2
    >>> set(atoms.get_chemical_symbols())
    {...}


Delete all helium atoms
-----------------------

Indexing with lists of booleans:

>>> del atoms[[atom.symbol == 'He' for atom in atoms]]

or

>>> del atoms[[symbol == 'He' for symbol in atoms.get_chemical_symbols()]]

Indexing with list of indices:

>>> del atoms[[atom.index for atom in atoms if atom.symbol == 'He']]

Alternative solution using ``atoms.numbers``
(``numpy.ndarray`` of atomic numbers):

>>> from ase.data import atomic_numbers
>>> atomic_numbers['He']
2
>>> del atoms[atoms.numbers == 2]


Trajectories
============

Append one trajectory to the end of another
-------------------------------------------

>>> from ase.io import Trajectory
>>> t1 = Trajectory('t1.traj')
>>> t2 = Trajectory('t2.traj')
>>> for atoms in t2:
...     t1.write(atoms)
>>> t1.close()


Input/output
============

Convert from one format to another
----------------------------------

>>> from ase.io import read, write
>>> write('abc.xyz', read('abc.traj'))
