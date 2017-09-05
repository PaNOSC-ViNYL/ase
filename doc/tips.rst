===============
Tips and tricks
===============

In order to get the most out out the tips below (and ASE in general), it
is a good idea to get to know the Python language and the NumPy library well.
See:

* https://www.python.org/
* https://docs.scipy.org/doc/numpy/

.. contents::


Atoms objects
=============

Species
-------

>>> from ase import Atoms
>>> atoms = Atoms('CH4')
>>> len(set(atoms.numbers))  # number of species
2
>>> set(atoms.get_chemical_symbols())  # set of species
{'C', 'H'}


Indexing
--------

>>> atoms
Atoms(symbols='CH4', pbc=False)
>>> [atom.index for atom in atoms if atom.symbol == 'H']
[1, 2, 3, 4]
>>> atoms[[atom.index for atom in atoms if atom.symbol == 'H']]
Atoms(symbols='H4', pbc=False)

Indexing with lists of booleans:

>>> atoms.numbers == 1
array([False,  True,  True,  True,  True], dtype=bool)
>>> atoms[atoms.numbers == 1]
Atoms(symbols='H4', pbc=False)

Three equivalent ways to delete carbon atoms:

>>> del atoms[atoms.numbers == 6]
>>> del atoms[[atom.index for atom in atoms if atom.symbol == 'C']]
>>> del atoms[[atom.symbol == 'C' for atom in atoms]]

Swap the positions of two atoms with index 3 and 4:

>>> atoms.positions[[3, 4]] = atoms.positions[[4, 3]]


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
>>> write('abc.xyz', read('abc.traj'))
