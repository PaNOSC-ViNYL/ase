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


.. _ASAP: http://wiki.fysik.dtu.dk/asap
.. _GPAW: https://wiki.fysik.dtu.dk/gpaw/documentation/xc/vdwcorrection.html


Version 3.4.1
=============

11 August 2010: :trac:`tags/3.4.1 <../tags/3.4.1>`.
