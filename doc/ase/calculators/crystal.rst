.. module:: ase.calculators.crystal

=========
CRYSTAL14
=========

Introduction
============

The CRYSTAL_ code is a Hartree-Fock and density functional theory
code using Gaussian localized basis set functions. CRYSTAL_
can handle systems periodic in 0 (molecules, 0D), 1 (polymers, 1D),
2 (slabs, 2D), and 3 dimensions (crystals, 3D.
This interface makes it possible to use CRYSTAL_ as a calculator
in ASE.

.. _Crystal: http://www.crystal.unito.it/


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
========================================
The calculator calls the CRYSTAL_ code only
to perform single-point+gradient calculations.
The file 'fort.34' contains the input geometry.

Below follows a list with a selection of paramenters.

==============  =========  ==============  ============================
keyword         type       default value   description
==============  =========  ==============  ============================
``restart``     ``bool``   None            Restart old calculation
``xc``          various    'hf'            Hamiltonian. HF, MP2 or DFT
                                           methods available
``spin``        ``bool``   False           Spin polarization
``guess``       ``bool``   True            Read wf from fort.20 file
                                           when present
``basis``       ``str``    'custom'        Read basis set from
                                           basis file
``kpts``        various    None            **k**-point sampling if
                                           calculation is periodic
``isp``         ``int``    1               Density of the Gilat net
                                           with respect to Monkhorst-
                                           Pack
``smearing``    ``float``  None            Smearing. Only Fermi-Dirac
                                           available
``otherkey``    ``list``   []              All other CRYSTAL keywords
==============  =========  ==============  ============================

For parameters not set in ``otherkey`` CRYSTAL will set the default value.
See the official `CRYSTAL manual`_ for more details.

.. _CRYSTAL manual: http://www.crystal.unito.it/Manuals/crystal14.pdf

Exchange-correlation functionals
================================

The ``xc`` parameter is used to define the hamiltonian of the
calculation. In case of a density-functional theory hamiltonian,
this parameter set the functional to use. A single string defines
a standalone functional (see `CRYSTAL manual`_), a tuple of strings
set the first string as EXCHANGE and the second string as 'CORRELAT'
(see `CRYSTAL manual`_ for more details).

