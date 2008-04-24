.. module::  calculators
   :synopsis: Energy, force, stress calculators.

Calculators
===========

A calculator is a black box that can take atomic numbers and atomic
positions from an :class:`ase.Atoms` object and calculate energy and/or
forces and/or stresses.


Supported calculators
=====================


=======  =================================================
GPAW_    Grid-based real-space PAW code
Asap_    Highly efficient EMT code (written in C++)
Dacapo_  A planewave ultra-soft pseudopotential code
EMT      Effective Medium Theory calculator
Siesta   LCAO pseudopotential code
MMTK     Library for molecular siulations 
=======  =================================================
  

.. _Asap: http://wiki.fysik.dtu.dk/Asap
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo
.. _GPAW: http://wiki.fysik.dtu.dk/gpaw

Siesta and MMTK modules are ASE wrappers, and are not mantained by CAMd.

EMT is a pure python implementation of an Effective Medium calculator.

.. toctree::

   emt
   siesta
   mmtk


Tools for building new calculators
==================================

.. autoclass:: ase.calculators.neighborlist.NeighborList
