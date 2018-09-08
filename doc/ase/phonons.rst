.. module:: ase.phonons

===================
Phonon calculations
===================

Module for calculating vibrational normal modes for periodic systems using the
so-called small displacement method (see e.g. [Alfe]_). So far, space-group
symmetries are not exploited to reduce the number of atomic displacements that
must be calculated and subsequent symmetrization of the force constants.

For polar materials the dynamical matrix at the zone center acquires a
non-analytical contribution that accounts for the LO-TO splitting. This
contribution requires additional functionality to evaluate and is not included
in the present implementation. Its implementation in conjunction with the small
displacement method is described in [Wang]_.


Example
=======

Simple example showing how to calculate the phonon dispersion for bulk aluminum
using a 7x7x7 supercell within effective medium theory:

.. literalinclude:: phonons_Al_fcc.py
   :start-after: creates:
   :end-before: End of literalinclude

.. image:: Al_phonon.png

Mode inspection:

.. literalinclude:: phonons_Al_fcc.py
   :start-after: Literalinclude start modes
   :end-before: Literalinclude end modes

.. image:: Al_mode.*

.. [Alfe] D. Alfe, PHON: A program to calculate phonons using the small
          displacement method, Comput. Phys. Commun. 180, 2622 (2009)
.. [Wang] Y. Wang *et al.*, A mixed-space approach to first-principles
          calculations of phonon frequencies for polar materials, J. Phys.:
          Cond. Matter 22, 202201 (2010)


List of all Methods
===================

.. autoclass:: Phonons
   :members:
