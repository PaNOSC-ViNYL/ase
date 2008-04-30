.. module::  calculators
   :synopsis: Energy, force and stress calculators.

Calculators
===========

For ASE, a calculator is a black box that can take atomic numbers and atomic
positions from an :class:`Atoms` object and calculate the energy and forces and sometimes also stresses.


Supported calculators
=====================


===============  ===========================================  ============
Code             Description                                  Type
===============  ===========================================  ============
GPAW_            Grid-based real-space PAW code               :term:`DFT`,
                                                              :term:`HF`
Asap_            Highly efficient EMT code (written in C++)   :term:`EMT`
Dacapo_          A planewave ultra-soft pseudopotential code  :term:`DFT`
:class:`EMT`     Effective Medium Theory calculator           :term:`EMT`
:class:`Siesta`  LCAO pseudopotential code                    :term:`DFT`
:class:`MMTK`    XXX Library for molecular simulations 
===============  ===========================================  ============
  

.. _Asap: http://wiki.fysik.dtu.dk/Asap
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo
.. _GPAW: http://wiki.fysik.dtu.dk/gpaw


The calculators can be divided in three groups:

* GPAW, Asap, and Dacapo have their own native ASE interfaces.

* SIESTA and MMTK have Python wrappers in the ASE package, but the
  actual codes are not part of ASE.

* EMT is a pure python implementation of the Effective Medium Theory
  potential and it is included in the ASE package.



.. toctree::

   emt
   siesta
   mmtk



Using calculators
=================

A calculator can be attached to an :class:`Atoms` object like this::

  atoms = Atoms(..., calculator=Siesta())

or like this::

  atoms = Atoms(...)
  atoms.set_calculator(Siesta())



Building new calculators
========================

Adding an ASE interface to your favorite force-calculator can be very
simple.  Take a look at the Python wrapper we have in the ASE code for
using the SIESTA_ code with ASE: :svn:`ase/calculators/siesta.py`.
Here, a :class:`Siesta` class is defined.  An instance of this class
will simply write the fdf input-file, start the SIESTA Fortran
program, and finally read the energy, forces and stresses from the
text output-file.


.. hint::

   The :class:`EMT` potential and the GPAW_ DFT calculator both make
   use of ASE's built-in neighbor-list class:

   .. autoclass:: ase.calculators.neighborlist.NeighborList

   For details, see :epydoc:`ase.calculators.neighborlist.NeighborList`.
