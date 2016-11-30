===============
Tips and tricks
===============

Atoms objects
=============

Species
-------

>>> from ase import Atoms
>>> atoms = Atoms('CH4')
>>> len(set(atoms.numbers))  # number of species
2
>>> set(atoms.get_chemical_symbols())
{...}


Indexing
--------

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


combine two to one
------------------
