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

::

  from gpaw.analyse.overlap import Overlap
   
  rr = ResonantRaman(s, LrTDDFT, gsname=gsname, exname=exname,
                     exkwargs={'energy_range':erange, 'eps':0.2},
	             overlap=lambda x, y: Overlap(x).pseudo(y),
               )
