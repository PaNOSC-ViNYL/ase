Resonant and non-resonant Raman spectra
=======================================

Note: :ref:`Siesta Raman` are possible also.

Raman spectra can be calculated in various approximations [1]_.
While the examples below are using GPAW_ explicitely,
the modules are intended to work with other calculators also.
The strategy is to calculate vibrational properties first and
obtain the spectra from these later.

1. Finite difference calculations
---------------------------------

It is recommended to do a vibrational analysis first by using
the :class:`~ase.vibrations.Vibrations` or  :class:`~ase.vibrations.Infrared`
modules. In the example of molecular hydrogen this is

.. literalinclude:: H2_ir.py

In the next step we perform a finite difference optical calculation
where the optical spectra are evaluated using TDDFT

.. literalinclude:: H2_optical.py
		    
Albrecht B+C terms need wave function overlaps at equilibrium and
displaced structures. These are assumed to be
calculated in the form

.. math::

  o_{ij} = \int d\vec{r} \; \phi_i^{{\rm disp},*}(\vec{r})
  \phi_j^{{\rm eq}}(\vec{r})
   
where :math:`\phi_j^{{\rm eq}}` is an orbital at equilibrium position
and :math:`\phi_i^{\rm disp}` is an orbital at displaced position.
This is implemented in ``Overlap`` in GPAW
(approximated by pseudo-wavefunction overlaps) and can be triggered
in ``ResonantRaman`` by::

  from gpaw.analyse.overlap import Overlap

  rr = ResonantRaman(atoms, LrTDDFT, exkwargs={'jend':3}
                     overlap=lambda x, y: Overlap(x).pseudo(y),
                     )


2. Analysis of the results
--------------------------

We assume that the steps above were performed and are able to analyse the
results in different approximations.

In order to do the full Albrecht analysis later we 
We save the standard names::

  # standard name for Vibrations
  gsname='vib'
  # standard name for Infrared
  gsname='ir'


Placzek
```````

The most popular form is the Placzeck approximation that is present in
two implementations. The simplest is the direct evaluation from
derivatives of the frequency dependent polarizability::

  from ase.vibrations.placzek import Placzek

  photonenergy = 7.5  # eV
  pz = Placzek()
  x, y = pz.get_spectrum(photonenergy, start=0, end=2000, method='frederiksen', type='Lorentzian')


The second implementation evaluates the derivatives differently allowing
for more analysis::

  from ase.vibrations.placzek import Profeta
  
  photonenergy = 7.5  # eV
  pr = Profeta(approximation='Placzek')
  x, y = pr.get_spectrum(photonenergy, start=0, end=2000, method='frederiksen', type='Lorentzian')

Both should lead to the same spectrum.

Albrecht
````````

``ResonantRaman`` calls the displaced excited state objects' function
``overlap`` with the matrix :math:`o_{ij}` and expects the function to
return the corresponding overlap matrix for the transition dipoles.
In case of Kohn-Sham transitions with :math:`i,j` for occupied
and :math:`\alpha,\beta` for empty orbitals, this is

.. math::

   O_{i\alpha,j\beta} = o_{ij}^* o_{\alpha\beta}

Example::

  from ase.vibrations.albrecht import Albrecht

  al = Albrecht()

.. _GPAW: http://wiki.fysik.dtu.dk/gpaw
  
.. [1] "Ab-initio wave-length dependent Raman spectra: Placzek approximation and beyond" Michael Walter, Michael Moseler `arXiv:1806.03840 <https://arxiv.org/abs/1806.03840>`_ [physics.chem-ph]

