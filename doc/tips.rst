===============
Tips and tricks
===============

https://www.python.org/
https://docs.scipy.org/doc/numpy/


Atoms objects
=============

Species
-------

>>> from ase import Atoms
>>> atoms = Atoms('CH4')
>>> len(set(atoms.numbers))  # number of species
2
>>> set(atoms.get_chemical_symbols())
set(['H', 'C'])


Indexing
--------

>>> atoms
Atoms('CH4')
>>> [atom.index for atom in atoms if atom.symbol == 'H']
[1, 2, 3, 4]
>>> atoms[[atom.index for atom in atoms if atom.symbol == 'H']]
Atoms('H4')

Indexing with lists of booleans:

>>> atoms.numbers == 1
>>> atoms[atoms.numbers == 1]
>>> del atoms[[atom.symbol == 'He' for atom in atoms]]


Alternative solution using ``atoms.numbers``
(``numpy.ndarray`` of atomic numbers):


Trajectories
============

Append one trajectory to the end of another
-------------------------------------------

.. testsetup::

    from ase.io import write
    from ase import Atoms
    write('t1.traj', Atoms('H'))
    write('t2.traj', Atoms('H'))
    write('abc.traj', Atoms('H'))

>>> from ase.io import Trajectory
>>> t1 = Trajectory('t1.traj', 'a')
>>> t2 = Trajectory('t2.traj')
>>> for atoms in t2:
...     t1.write(atoms)
>>> t1.close()


Input/output
============

Convert from one format to another
----------------------------------

>>> from ase.io import read, write
>>> write('abc.traj', read('abc.traj'))
