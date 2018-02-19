.. module:: ase.calculators.crystal

=========
CRYSTAL14
=========

Introduction
============

The CRYSTAL_ simulation package is a Hartree-Fock and density
functional theory code using Gaussian localized basis functions.
CRYSTAL_ can handle systems periodic in 0 (molecules, 0D), 1 (polymers, 1D),
2 (slabs, 2D), and 3 dimensions (crystals, 3D).
This interface makes possible to use CRYSTAL_ as a calculator
in ASE.

.. _CRYSTAL: http://www.crystal.unito.it/


Environment variables
=====================

Set environment variables in your configuration file (what is the name
of the command to be run). It is mandatory to set the input file as
"INPUT" and the standard output as "OUTPUT".

- bash::

  $ export ASE_CRYSTAL_COMMAND="/bin/CRY14/crystal < INPUT > OUTPUT 2>&1"  (an example)

- csh/tcsh::

  $ setenv ASE_CRYSTAL_COMMAND "/my_disk/my_name/bin/crystal < INPUT > OUTPUT 2>&1"  (an example)


CRYSTAL Calculator (a FileIOCalculator)
=======================================

The calculator calls the CRYSTAL_ code only
to perform single point and gradient calculations.
The file 'fort.34' contains the input geometry and
the 'fort.20' contains the wave function in a binary
format.

Below follows a list with a selection of parameters.

==============  =========  ===============  ============================
keyword         type       default value    description
==============  =========  ===============  ============================
``restart``     ``bool``   None             Restart old calculation
``xc``          various    'HF'             Hamiltonian. HF, MP2 or DFT
                                            methods available
``spinpol``     ``bool``   False            Spin polarization
``guess``       ``bool``   True             Read wf from fort.20 file
                                            when present
``basis``       ``str``    'custom'         Read basis set from
                                            basis file
``kpts``        various    None or (1,1,1)  **k**-point sampling if
                                            calculation is periodic
``isp``         ``int``    1                Density of the Gilat net
                                            with respect to Monkhorst-
                                            Pack
``smearing``    ``float``  None             Smearing. Only Fermi-Dirac
                                            available
``otherkeys``   ``list``   []               All other CRYSTAL keywords
==============  =========  ===============  ============================

For parameters not set in ``otherkeys`` CRYSTAL_ will set the default value.
See the official `CRYSTAL manual`_ for more details.

.. _CRYSTAL manual: http://www.crystal.unito.it/Manuals/crystal14.pdf


Exchange-correlation functionals
================================

The ``xc`` parameter is used to define the method used for the
calculation. Available options are Hartree-Fock ('HF'), second order
perturbation theory ('MP2') and the density-functional theory where ``xc``
defines the exchange and correlation functional. In the latter case
a single string defines a standalone functional (see `CRYSTAL manual`_),
a tuple of strings set the first string as EXCHANGE and the second
string as 'CORRELAT' (see `CRYSTAL manual`_ for more details).

.. code-block:: python

  calc = CRYSTAL(xc=('PBE','LYP'))


Setups
======

The CRYSTAL_ simulation package has few built-in basis sets, which
can be set in the calculation using the ``basis`` parameter, e. g.:

.. code-block:: python

  calc = CRYSTAL(xc='PBE', basis='sto-3g')

The default is to read from an external basis set. A library of
basis sets in CRYSTAL_ format can be found on the
website `CRYSTAL basis sets`_.

.. _CRYSTAL basis sets: http://www.crystal.unito.it/basis-sets.php

In this case a file named 'basis'  must be present in the working directory
and must contain the basis sets for all the atom species.

.. note::

   The CRYSTAL_ simulation package allows to set up to three different
   all electron basis sets and/or two valence electron basis sets for
   the same atomic species (see `CRYSTAL manual`_ page 21 for more details).

   The number to be added to the atomic number reported in the 'basis'
   file must be specified as an ``Atoms()`` class tag:

   >>> geom[0].tag = 100

   In this case '100' will be summed to the atomic number of the first atom
   in the 'fort.34' geometry file (e. g. '6', Carbon, becomes '106').


Spin-polarized calculation
==========================

If the atoms object has non-zero magnetic moments, a spin-polarized
calculation will be performed by default.
It is also possible to manually tell the calculator to perform a
spin-polarized calculation through the parameter ``spinpol``:

.. code-block:: python

  calc = CRYSTAL(xc='PBE', spinpol=True)


Brillouin-zone sampling
=======================

Brillouin-zone sampling is controlled by ``kpts``. This parameter
can be set to a sequence of three int values, e.g. (2, 2, 3),
which define a regular Monkhorst-Pack grid. If it is not defined a
``gamma`` calculation will be performed.
For 2D calculations ``kpts[2]`` will be to set to one, for 1D ones
also ``kpts[1]`` will be set to unity.
For molecular calculations (0D) any definition of the ``kpts``
parameter will be ignored.

The ``isp`` parameter can be used to define the relative
density of the auxiliary Gilat net (see `CRYSTAL manual`_):

.. code-block:: python

  calc = CRYSTAL(xc='PBE', kpts=(2, 2, 2), isp=2)

In this example the resulting Gilat net would be (4, 4, 4).


Reading an external wave function
=================================

The calculator reads by default the wave function stored in
the 'fort.20' file if present (``guess=True``).
If this parameter is set to False the code will calculate the
wave function from scratch at any step, slowing down the perfromances.


Code related keywords
=====================

The CRYSTAL_ simulation package allows for many other keywords.
Most of them can be specified through the ``otherkeys`` parameter.

.. code-block:: python

  calc = CRYSTAL(xc='PBE', otherkeys=['scfdir', 'anderson',
                                      ['maxcycles', '500'],
                                      ['toldee', '6'],
                                      ['tolinteg', '7 7 7 7 14'],
                                      ['fmixing', '90']])
