.. module:: structure

=================
Atomic structures
=================

.. seealso:: 

   * The :mod:`lattice` module
   * The :mod:`~lattice.surface` module


Bulk crystals
=============

.. autofunction:: ase.structure.bulk

examples:

>>> from ase.structure import bulk
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


Nanotubes
=========

...


Graphene nanoribbons
====================

...
