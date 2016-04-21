.. module:: ase.build

==============
Builing things
==============

**Common bulk crystals**

The :func:`ase.lattice.bulk` function can be used to create the most
common bulk crystal structures.  The function creates a single unit cell
oriented such that the number of atoms in the cell is minimal.

Read more: :ref:`bulk-crystal-section`.


**Common surfaces**

The :mod:`ase.build` module contains a number of
functions for creating the most common surfaces in a minimal unit
cell, and for adding adsorbates to these surfaces.

Read more: :ref:`lattice-surface-section`.


**Nanotubes and nanoribbons**

The functions :func:`ase.build.nanotube` and
:func:`ase.build.graphene_nanoribbon` can be used to create Carbon
nanotubes and graphene sheets or nanoribbons.  Per default, they
create Carbon nanotubes and sheets, but other elements can be used.

Read more:  :ref:`nanotubes-section` and :ref:`nanoribbons-section`.

**Generally oriented bulk crystals and surfaces**

The :mod:`ase.lattice` module contains functions for creating most common
crystal structures with arbitrary orientation.  The user can specify
the desired Miller index along the three axes of the simulation, and
the smallest periodic structure fulfilling this specification is
created.  Thirteen of the 14 Bravais lattices are supported by the
module, as are a few lattices with a basis, and lattices for some of
the most common compounds/alloys.  The modules makes it possible to
define further lattices based on the supported Bravais lattices.

Both bulk crystals and surfaces can be created.

Read more: :ref:`general-crystal-section`.

**Molecules**

Some common molecules can be constructed using the
:func:`ase.build.molecule` function.

Read more: :ref:`molecular-data`.



:func:`~ase.build.bulk`
:func:`~ase.build.surface`
:func:`~ase.build.molecule`

.. list-table::
    
    * - :func:`~ase.build.fcc100`
      - :func:`~ase.build.fcc110`
      - :func:`~ase.build.fcc111`
    * - :func:`~ase.build.fcc211`
      - :func:`~ase.build.bcc100`
      - :func:`~ase.build.bcc110`
    * - :func:`~ase.build.bcc111`
      - :func:`~ase.build.hcp0001`
      - :func:`~ase.build.hcp10m10`
    * - :func:`~ase.build.diamond100`
      - :func:`~ase.build.diamond111`
      - :func:`~ase.build.fcc111_root`
    * - :func:`~ase.build.bcc111_root`
      - :func:`~ase.build.hcp0001_root`
      - :func:`~ase.build.mx2`
    * - :func:`~ase.build.add_adsorbate`
      - :func:`~ase.build.add_vacuum`
      - :func:`~ase.build.root_surface`


      :func:`~ase.build.nanotube`
:func:`~ase.build.graphene_nanoribbon`

:func:`~ase.build.cut`
:func:`~ase.build.stack`
:func:`~ase.build.sort`
:func:`~ase.build.minimize_tilt`
:func:`~ase.build.niggli_reduce`
:func:`~ase.build.rotate`
:func:`~ase.build.minimize_rotation_and_translation`


.. seealso::

   * The :mod:`ase.lattice` module
   * The :mod:`ase.spacegroup` module
   * The :mod:`ase.geometry` module

   
.. toctree::
   :maxdepth: 2

   surface
   tools
   
   
Molecules
=========

The G2-database of common molecules is available:

.. autofunction:: molecule

Example::

>>> from ase.build import molecule
>>> atoms = molecule('H2O')

To see a list of available molecules in the default set, use::

>>> from ase.collection import g2
>>> g2.names


.. _bulk-crystal-section:

Common bulk crystals
====================

.. autofunction:: bulk

examples:

>>> from ase.lattice import bulk
>>> a1 = bulk('Cu', 'fcc', a=3.6)
>>> a2 = bulk('Cu', 'fcc', a=3.6, orthorhombic=True)
>>> a3 = bulk('Cu', 'fcc', a=3.6, cubic=True)
>>> a1.cell
array([[ 0. ,  1.8,  1.8],
       [ 1.8,  0. ,  1.8],
       [ 1.8,  1.8,  0. ]])
>>> a2.cell
array([[ 2.54558441,  0.        ,  0.        ],
       [ 0.        ,  2.54558441,  0.        ],
       [ 0.        ,  0.        ,  3.6       ]])
>>> a3.cell
array([[ 3.6,  0. ,  0. ],
       [ 0. ,  3.6,  0. ],
       [ 0. ,  0. ,  3.6]])

|a1| |a2| |a3|

.. |a1| image:: a1.png
.. |a2| image:: a2.png
.. |a3| image:: a3.png


.. _nanotubes-section:

Nanotubes
=========

.. autofunction:: nanotube

examples:

>>> from ase.build import nanotube
>>> cnt1 = nanotube(6, 0, length=4)
>>> cnt2 = nanotube(3, 3, length=6, bond=1.4, symbol='Si')

|cnt1| |cnt2|

.. |cnt1| image:: cnt1.png
.. |cnt2| image:: cnt2.png


.. _nanoribbons-section:

Graphene nanoribbons
====================

.. autofunction:: graphene_nanoribbon

examples:

>>> from ase.build import graphene_nanoribbon
>>> gnr1 = graphene_nanoribbon(3, 4, type='armchair', saturated=True)
>>> gnr2 = graphene_nanoribbon(2, 6, type='zigzag', saturated=True,
>>>                            C_H=1.1, C_C=1.4, vacuum=6.0,
>>>                            magnetic=True, initial_mag=1.12)

|gnr1| |gnr2|

.. |gnr1| image:: gnr1.png
.. |gnr2| image:: gnr2.png



ASE contains a number of modules for setting up atomic structures,
mainly molecules, bulk crystals and surfaces.  Some of these modules
have overlapping functionality, but strike a different balance between
flexibility and ease-of-use.



