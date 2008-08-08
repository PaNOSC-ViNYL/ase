.. module:: units

=====
Units
=====

XXX more about units here!

Physical units are defined in the :svn:`ase/units.py` module.  Electron volts
(``eV``) and angstroms (``Ang``) are defined as 1.0.  Other units are
``nm``, ``Bohr``, ``Hartree`` or ``Ha``, ``kJ``, ``kcal``, ``mol``,
``Rydberg`` or ``Ry``, ``second``, ``fs`` and ``kB``.  Example::

  >>> from ase.units import Bohr
  >>> 2 * Bohr
  1.0583545150138329
