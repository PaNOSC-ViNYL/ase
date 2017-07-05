.. module:: ase.collections

===========
Collections
===========

:data:`s22`, :data:`dcdft`, :data:`g2`

.. autoclass:: ase.collections.collection.Collection

.. _s22:

S22 database of weakly interacting dimers and complexes
=======================================================

.. data:: s22

S22 geometry data are from:

    P. Jurecka, J. Sponer, J. Cerny, P. Hobza; Phys Chem Chem Phys 2006, 8 (17), 1985-1993.

See http://www.begdb.com/index.php?action=106a6c241b8797f52e1e77317b96a201 for
the original files. All geometries are optimized at either the CCSD(T) or MP2
level except for the methyl amide dimers where only the hydrogen position is
optimized at the DFT level.

The S22 interaction energies are all calculated using both CCSD(T)/CBS counter
poised corrected (CP) and MP2 /CBS CP. The original S22 interaction energies
are listed in the above references. The S22 energies used here are from
Takatani, T. et al., J. Chem. Phys., 132, 144104 (2010) where a large and more
complete basis set has been used for all database members.


DeltaCodesDFT
=============

PROVISIONAL:  This stuff may not have its final form.  Use with care!

.. data:: dcdft

Structures and data from:

    https://github.com/molmod/DeltaCodesDFT

See also [Lejaeghere2014]_.

.. [Lejaeghere2014]

    K. Lejaeghere, V. Van Speybroeck, G. Van Oost, and S. Cottenier:
    "Error estimates for solid-state density-functional theory predictions:
    an overview by means of the ground-state elemental crystals",
    Crit. Rev. Solid State (2014).
    http://dx.doi.org/10.1080/10408436.2013.772503

This collection has WIEN2k and experimental data for:

* volume per atom
* bulk-modulus (in GPa)
* pressure derivative of bulk-modulus

>>> from ase.collections import dcdft
>>> dct = dcdft.data['Cu']
>>> for key, val in sorted(dct.items()):
...     print('{:15}: {:.3f}'.format(key, val))
exp_B          : 144.279
exp_Bp         : 4.880
exp_volume     : 11.647
wien2k_B       : 141.335
wien2k_Bp      : 4.860
wien2k_volume  : 11.951


G2 ...
======

.. data:: g2
