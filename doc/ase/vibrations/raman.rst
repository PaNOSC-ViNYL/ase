.. module:: ase.vibrations

=============
Raman spectra
=============

Raman spectra can be calculated in various approximations.
While the examples below are using gpaw explicitely,
the modules are intended to work with other calculators also.

The strategy is to calculate vibrational properties first and
obtain the spectra from these.


Finite difference resonant Raman calculation
--------------------------------------------

The ResonantRaman class allows for a finite difference resonant
Raman calculation that can be used for Placzek or Albrecht analysis later.

Example::

  from ase.vibrations.albrecht import Albrecht

Albrecht
--------

Albrecht B+C terms need wave function overlaps at equilibrium and
displaced structures. These are assumed to be
calculated in the form

.. math::

  o_{ij} = \int d\vec{r} \; \phi_i^{{\rm disp},*}(\vec{r})
  \phi_j^{{\rm eq}}(\vec{r})
   
where :math:`\phi_j^{{\rm eq}}` is an orbital at equilibrium position
and :math:`\phi_i^{\rm disp}` is an orbital at displaced position.
This is implemented in ``Overlap`` in GPAW
(approximated by pseudo-wavefunction overlaps) and can be called
from ``ResonantRaman`` by::

  from gpaw.analyse.overlap import Overlap
  rr = ResonantRaman(s, LrTDDFT, gsname=gsname, exname=exname,
                   exkwargs={'energy_range':erange, 'eps':0.2},
	           #overlap=True,
	           overlap=lambda x, y: Overlap(x).pseudo(y),
  )

``ResonantRaman`` calls the displaced excited state objects' function
``overlap`` with the matrix :math:`o_{ij}` and expects the function to
return the corresponding overlap matrix for the transition dipoles.
In case of Kohn-Sham transitions with :math:`i,j` for occupied
and :math:`\alpha,\beta` for empty orbitals, this is

.. math::

   O_{i\alpha,j\beta} = o_{ij}^* o_{\alpha\beta}


Placzek
-------  

The most popular form is the Placzeck approximation.

