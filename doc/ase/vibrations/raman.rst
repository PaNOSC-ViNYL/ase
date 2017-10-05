.. module:: ase.vibrations

=============
Raman spectra
=============

Raman spectra can be calculated in various approximations.
While the examples below are using gpaw explicitely,
the modules are intended to work with other calculators also.

The strategy is to calculate vibrational properties first and
obtain the spectra from these.



Example::

  from ase.vibrations.albrecht import Albrecht

Placzek
-------  

The most popular form is the Placzeck approximation.

Albrecht
--------
