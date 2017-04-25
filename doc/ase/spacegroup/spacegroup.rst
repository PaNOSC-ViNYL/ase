.. module:: ase.spacegroup

===============================
Using the spacegroup subpackage
===============================

The most evident usage of the spacegroup subpackage is to set up an
initial unit of a bulk structure. For this you only need to supply the
unique atoms and their scaled positions, space group and lattice
parameters.


Examples of setting up bulk structures
======================================

We start by showing some examples of how to set up some common or
interesting bulk structures using
**ase.spacegroup.crystal()**.  This function takes a lot of
arguments, of which the most important are:

    symbols : string | sequence of strings
        Either a string formula or sequence of element
        symbols. E.g. 'NaCl' and ('Na', 'Cl') are equivalent.
    basis : list of scaled coordinates
        Coordinates of the non-equivalent sites in units of the
        lattice vectors.
    spacegroup : int | string | Spacegroup instance
        Space group given either as its number in International Tables
        or as its Hermann-Mauguin (or international) symbol.
    setting : 1 | 2
        Space group setting.
    cellpar : [a, b, c, alpha, beta, gamma]
        Cell parameters with angles in degree. Is not used when ``cell``
        is given.

        
Aluminium (fcc)
---------------

.. image:: spacegroup-al.png

.. literalinclude:: spacegroup-al.py

The *spacegroup* argument can also be entered with its Hermann-Mauguin
symbol, e.g. *spacegroup=225* is equivalent to *spacegroup='F m -3 m'*.


Iron (bcc)
----------

.. image:: spacegroup-fe.png

.. literalinclude:: spacegroup-fe.py


Magnesium (hcp)
---------------

.. image:: spacegroup-mg.png

.. literalinclude:: spacegroup-mg.py


Diamond
-------

.. image:: spacegroup-diamond.png

.. literalinclude:: spacegroup-diamond.py


Sodium chloride
---------------

.. image:: spacegroup-nacl.png

.. literalinclude:: spacegroup-nacl.py


Rutile
------

.. image:: spacegroup-rutile.png

.. literalinclude:: spacegroup-rutile.py


CoSb3 skutterudite
------------------

.. image:: spacegroup-skutterudite.png

Skutterudites_ are quite interesting structures with 32 atoms
in the unit cell.

.. _Skutterudites: http://en.wikipedia.org/wiki/Skutterudite

.. literalinclude:: spacegroup-skutterudite.py

Often this structure is visualised with the Cobalt atoms on the
corners. This can easily be accomplished with ASE using
:func:`ase.build.cut`. Below is the *origo* argument used to
put the Cobalt atom on the corners and *extend* to include all corner
and edge atoms, even those belonging to neighbouring unit cells.

.. image:: spacegroup-cosb3.png

.. literalinclude:: spacegroup-cosb3.py


The Spacegroup class
====================

The :class:`ase.spacegroup.Spacegroup` class is used
internally by the :func:`ase.spacegroup.crystal` function, but
might sometimes also be useful if you want to know e.g. the symmetry
operations of a given space group. Instances of the
:class:`ase.spacegroup.Spacegroup` class are immutable
objects holding space group information, such as symmetry operations.

Let us e.g. consider the fcc structure. To print information about the
space group, do

>>> from ase.spacegroup import Spacegroup
>>> sg = Spacegroup(225)
>>> print(sg)
225   F m -3 m
  setting 1
  centrosymmetric 1
  primitive vectors
     0.0000000000  0.5000000000  0.5000000000
     0.5000000000  0.0000000000  0.5000000000
     0.5000000000  0.5000000000  0.0000000000
  reciprocal vectors
     -1   1   1
      1  -1   1
      1   1  -1
  4 subtranslations
     0.0000000000  0.0000000000  0.0000000000
     0.0000000000  0.5000000000  0.5000000000
     0.5000000000  0.0000000000  0.5000000000
     0.5000000000  0.5000000000  0.0000000000
  24 symmetry operations (rot+trans)
    1  0  0     0  1  0     0  0  1    0.0000000000  0.0000000000  0.0000000000
   -1  0  0     0 -1  0     0  0  1    0.0000000000  0.0000000000  0.0000000000
   -1  0  0     0  1  0     0  0 -1    0.0000000000  0.0000000000  0.0000000000
    1  0  0     0 -1  0     0  0 -1    0.0000000000  0.0000000000  0.0000000000
    0  0  1     1  0  0     0  1  0    0.0000000000  0.0000000000  0.0000000000
    0  0  1    -1  0  0     0 -1  0    0.0000000000  0.0000000000  0.0000000000
    0  0 -1    -1  0  0     0  1  0    0.0000000000  0.0000000000  0.0000000000
    0  0 -1     1  0  0     0 -1  0    0.0000000000  0.0000000000  0.0000000000
    0  1  0     0  0  1     1  0  0    0.0000000000  0.0000000000  0.0000000000
    0 -1  0     0  0  1    -1  0  0    0.0000000000  0.0000000000  0.0000000000
    0  1  0     0  0 -1    -1  0  0    0.0000000000  0.0000000000  0.0000000000
    0 -1  0     0  0 -1     1  0  0    0.0000000000  0.0000000000  0.0000000000
    0  1  0     1  0  0     0  0 -1    0.0000000000  0.0000000000  0.0000000000
    0 -1  0    -1  0  0     0  0 -1    0.0000000000  0.0000000000  0.0000000000
    0  1  0    -1  0  0     0  0  1    0.0000000000  0.0000000000  0.0000000000
    0 -1  0     1  0  0     0  0  1    0.0000000000  0.0000000000  0.0000000000
    1  0  0     0  0  1     0 -1  0    0.0000000000  0.0000000000  0.0000000000
   -1  0  0     0  0  1     0  1  0    0.0000000000  0.0000000000  0.0000000000
   -1  0  0     0  0 -1     0 -1  0    0.0000000000  0.0000000000  0.0000000000
    1  0  0     0  0 -1     0  1  0    0.0000000000  0.0000000000  0.0000000000
    0  0  1     0  1  0    -1  0  0    0.0000000000  0.0000000000  0.0000000000
    0  0  1     0 -1  0     1  0  0    0.0000000000  0.0000000000  0.0000000000
    0  0 -1     0  1  0     1  0  0    0.0000000000  0.0000000000  0.0000000000
    0  0 -1     0 -1  0    -1  0  0    0.0000000000  0.0000000000  0.0000000000
<BLANKLINE>
<BLANKLINE>

Or, if you want to figure out what sites in the unit cell are
equivalent to (0, 0, 0.5), simply do

>>> sites,kinds = sg.equivalent_sites([(0, 0, 0.5)])
>>> sites
array([[ 0. ,  0. ,  0.5],
       [ 0.5,  0. ,  0. ],
       [ 0. ,  0.5,  0. ],
       [ 0.5,  0.5,  0.5]])
>>> kinds
[0, 0, 0, 0]

where *sites* will be an array containing the scaled positions of the
four symmetry-equivalent sites.
