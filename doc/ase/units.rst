.. module:: ase.units

=====
Units
=====

Physical units are defined in the :git:`ase/units.py` module.  Electron volts
(``eV``) and angstroms (``Ang``) are defined as 1.0.
Other units are
``nm``, ``Bohr``, ``Hartree`` or ``Ha``, ``kJ``, ``kcal``, ``mol``,
``Rydberg`` or ``Ry``, ``second``, ``fs`` and ``kB``.

.. note::

    All constants are taken from the 1986 CODATA_.

.. _CODATA: http://physics.nist.gov/cuu/Constants/archive1986.html

Examples:

>>> from ase.units import *
>>> 2 * Bohr
1.0583545150138329
>>> 25 * Rydberg
340.14244569396635
>>> 100 * kJ/mol
1.0364272141304978
>>> 300 * kB
0.025852157076770025
>>> 0.1 * fs
0.009822693531550318
>>> print('1 Hartree =', Hartree * mol / kcal, 'kcal/mol')
1 Hartree = 627.50954059388 kcal/mol
