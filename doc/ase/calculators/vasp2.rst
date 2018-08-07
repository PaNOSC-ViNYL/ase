.. module:: ase.calculators.vasp.vasp2

===========
VASP 2.0
===========

Introduction
============

This module introduces an updated version of the ASE VASP_ calculator,
which adds the functionality of the :class:`~ase.calculators.calculator.FileIOCalculator`.
This allows a more general usage of the other ASE methods,
such as :class:`~ase.dft.band_structure.BandStructure`.

For a general introduction please refer to the :mod:`~ase.calculators.vasp` calculator
documentation, as this is just a list of things which have changed.

.. _VASP: http://cms.mpi.univie.ac.at/vasp/

.. warning::
   This calculator is currently in BETA testing. If you are not comfortable
   testing new software, please use the old calculator object, see
   :mod:`~ase.calculators.vasp`.

.. note::
   If you encounter any bugs using this calculator, please report it as an issue
   on the `ASE github`_ or on the `ASE IRC`_.

.. _ASE github: https://gitlab.com/ase/ase
.. _ASE IRC: http://webchat.freenode.net/?randomnick=0&channels=ase

Environment variables
=====================

The calculator needs to know how to execute VASP. One way of doing this,
is by using the :meth:`~ase.calculators.vasp.Vasp2.command` method,
with instructions on how to execute vasp, e.g.::

  Vasp2(command='mpiexec vasp_std')

which requires that the executable :file:`vasp_std` is in your :envvar:`PATH`.
Alternatively, similar to the original implementation, one of the following
environment variables can be set: :envvar:`ASE_VASP_COMMAND`, :envvar:`VASP_COMMAND`
or :envvar:`VASP_SCRIPT` - note, that the environment variables are prioritized
in that order, so if :envvar:`ASE_VASP_COMMAND` is set, the two others are ignored.
The variables :envvar:`ASE_VASP_COMMAND` or :envvar:`VASP_COMMAND` should be
commands which executes vasp. Additionally, remember to set
the :envvar:`VASP_PP_PATH`. An example shell configuration could contain

.. highlight:: bash

::

   $ export ASE_VASP_COMMAND="mpiexec vasp_std"
   $ export VASP_PP_PATH=$HOME/vasp/mypps

.. highlight:: python

Alternatively, the :envvar:`VASP_SCRIPT` could be used, as described in the
original VASP_ calculator documentation.



Vasp 2.0 Calculator
=====================

The VASP specific keywords are unchanged, and should be included as described in
`VASP calculator`_. See the official `VASP manual`_ for more details on the
VASP specific keywords.


.. _VASP calculator: vasp.html#vasp-calculator
.. _VASP manual: http://cms.mpi.univie.ac.at/vasp/vasp/vasp.html

.. autoclass:: Vasp2


.. note::

   Parameters can be changed after the calculator has been constructed
   by using the :meth:`~ase.calculators.vasp.Vasp.set` method:

   >>> calc.set(prec='Accurate', ediff=1E-5)

   This would set the precision to Accurate and the break condition
   for the electronic SC-loop to ``1E-5`` eV.

Examples
========

The Vasp FileIO calculator now integrates with existing ASE functions, such as
:class:`~ase.dft.band_structure.BandStructure` or :class:`~ase.dft.bandgap.bandgap`.

Band structure with VASP
------------------------
.. _Si band structure: https://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure

The VASP manual has an example of creating a `Si band structure`_ - we can
easily reproduce a similar result, by using the ASE Vasp2 calculator.

We can use the ``directory`` keyword to control the folder in which the calculations
take place, and keep a more structured folder structure. The following script does the
initial calculations, in order to construct the band structure for silicon

.. code-block:: python

	from ase.build import bulk
	from ase.calculators.vasp import Vasp2

	si = bulk('Si')

	mydir = 'bandstructure'    # Directory where we will do the calculations

	# Make self-consistent ground state
	calc = Vasp2(kpts=(4, 4, 4), directory=mydir)

	si.set_calculator(calc)
	si.get_potential_energy()  # Run the calculation

	# Non-SC calculation along band path
	kpts = {'path': 'WGX',     # The BS path
	        'npoints': 30}     # Number of points along the path

	calc.set(isym=0,           # Turn off kpoint symmetry reduction
	         icharg=11,        # Non-SC calculation
    		 kpts=kpts)

	# Run the calculation
	si.get_potential_energy()

As this calculation might be longer, depending on your system, it may
be more convenient to split the plotting into a separate file, as all
of the VASP data is written to files. The plotting can then be achieved
by using the ``restart`` keyword, in a second script

.. code-block:: python

	from ase.calculators.vasp import Vasp2

	mydir = 'bandstructure'    # Directory where we did the calculations

	# Load the calculator from the VASP output files
	calc_load = Vasp2(restart=True, directory=mydir)

	bs = calc_load.band_structure() # ASE Band structure object
	bs.plot(emin=-13)               # Plot the band structure

Which results in the following image

.. image:: vasp_si_bandstructure.png

We could also find the band gap in the same calculation,

>>> from ase.dft.bandgap import bandgap
>>> bandgap(calc_load)
Gap: 0.474 eV
Transition (v -> c):
  (s=0, k=15, n=3, [0.000, 0.000, 0.000]) -> (s=0, k=27, n=4, [0.429, 0.000, 0.429])

.. note::
   When using hybrids, due to the exact-exchange calculations, one needs to treat
   the k-point sampling more carefully, see `VASP HSE band structure wiki`_.

   Currently, we have no functions to easily handle this issue, but may be added
   in the future.

.. _VASP HSE band structure wiki: https://cms.mpi.univie.ac.at/wiki/index.php/Si_HSE_bandstructure#Procedure_2:_0-weight_.28Fake.29_SC_procedure_.28works_DFT_.26_hybrid_functionals.29


Density of States
------------------------

Vasp2 also allows for quick access to the Density of States (DOS), through the ASE DOS module, see :class:`~ase.dft.dos.DOS`.
Quick access to this function, however, can be found by using the ``get_dos()`` function:

>>> energies, dos = calc.get_dos()
