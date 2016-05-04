.. module:: ase.collection

===========
Collections
===========

:data:`s22`, :data:`g2`

.. autoclass:: ase.collection.collection.Collection


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


G2 ...
======

.. data:: g2
