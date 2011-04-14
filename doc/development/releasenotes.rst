.. _releasenotes:

=============
Release notes
=============


Development version in trunk
============================

:trac:`trunk <>`.

* ...
* ...



Version 3.5.0
=============

13 April 2011: :trac:`tags/3.5.0 <../tags/3.5.0>`.

* Improved EMT potential:  uses a
  :class:`~ase.calculators.neighborlist.NeighborList` object and is
  now ASAP_ compatible.

* :mod:`BFGSLineSearch <optimize.bfgslinesearch>` is now the default
  (``QuasiNewton==BFGSLineSearch``).

* There is a new interface to the LAMMPS molecular dynamics code.

* New :mod:`phonons` module.

* Van der Waals corrections for DFT, see GPAW_ usage.

* updated gui interface: stability and usability improvements; povray
  render facility; updated expert user mode; enabled
  customization of colours, atomic radii; enabled user default
  settings via ~/.ase/gui.py 

* database library expanded to include: 1) the s22, s26 and s22x5 sets
  of van der Waals bonded dimers and complexes by the Hobza group, and 
  2) the DBH24 set of gas-phase reaction barrier heights by the Truhlar group.

.. _ASAP: http://wiki.fysik.dtu.dk/asap
.. _GPAW: https://wiki.fysik.dtu.dk/gpaw/documentation/xc/vdwcorrection.html


Version 3.4.1
=============

11 August 2010: :trac:`tags/3.4.1 <../tags/3.4.1>`.
