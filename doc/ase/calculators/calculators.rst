.. module::  calculators
   :synopsis: Energy, force, stress calculators.

Calculators
===========

A calculator is a black box that can take atomic numbers and atomic
positions from an :class:`ase.Atoms` object and calculate energy and/or
forces and/or stresses.


.. toctree::

   emt
   siesta


Other calculators
=================


=======  =================================================
Asap_    Highly efficient EMT code (written in C++)
GPAW_    Grid-based real-space PAW code
Dacapo_  A planewave ultra-soft pseudopotential code
MMTK     The Molecular Modeling Toolkit (work in progress)
=======  =================================================
  

.. _Asap: http://wiki.fysik.dtu.dk/Asap
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo
.. _GPAW: http://wiki.fysik.dtu.dk/gpaw


Tools for building new calculators
==================================

.. autoclass:: ase.calculators.neighborlist.NeighborList
