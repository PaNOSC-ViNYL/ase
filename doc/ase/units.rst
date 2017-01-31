.. module:: ase.units

=====
Units
=====

Physical units are defined in the :git:`ase/units.py` module.  Electron volts
(``eV``) and angstroms (``Ang``) are defined as 1.0.
Other units are (amongst others)
``nm``, ``Bohr``, ``Hartree`` or ``Ha``, ``kJ``, ``kcal``, ``mol``,
``Rydberg`` or ``Ry``, ``second``, ``fs`` and ``kB``.

.. note::

    As of version 3.12.0, all constants are taken from the 2014_
    version of the CODATA suggestions. Before that, all constants were taken
    from the 1986_ version. There is, however, a way to create all units
    depending on other versions of CODATA via the :func:`create_units` function
    (see Changing the CODATA version).

.. _1986: http://physics.nist.gov/cuu/Constants/archive1986.html
.. _2014: http://arxiv.org/pdf/1507.07956.pdf

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


Changing the CODATA version
---------------------------

If you just require an additional set of units that are based on a different
version of CODATA, you can use the ``create_units(codata_version)`` function.
It supports CODATA versions ``'1986'``, ``'1998'``, ``'2002'``, ``'2006'``,
``'2010'``, ``'2014'``. This function will return a dictionary with key-value
pairs of all the constants defined in the :mod:`ase.units` module, but based on
the CODATA version just selected:

>>> from ase.units import create_units
>>> units = create_units('1986')
>>> print(units['Bohr'])
0.5291772575069165
>>> units = create_units('2014')
>>> print(units['Bohr'])
0.5291772105638411

The dictionary also supports attribute access so it can be used as a drop-in
replacement for the module:

>>> from ase.units import create_units
>>> units = create_units('1986')
>>> units.Bohr
0.5291772575069165
>>> units = create_units('2014')
>>> units.Bohr
0.5291772105638411



.. autofunction:: create_units
