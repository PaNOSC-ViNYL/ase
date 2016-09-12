.. module:: ase.calculators.abinit

======
ABINIT
======

Introduction
============

ABINIT_ is a density-functional theory code based on pseudopotentials
and a planewave basis.


.. _ABINIT: http://www.abinit.org



Environment variables
=====================

.. highlight:: bash

.. envvar:: ASE_ABINIT_COMMAND

    Must be set to something like this::

        abinit < PREFIX.files > PREFIX.log

    where ``abinit`` is the executable (``abinis`` for version prior to 6).

.. envvar:: ABINIT_PP_PATH

    A directory containing the pseudopotential files (at least of
    :file:`.fhi` type).

Abinit does not provide tarballs of pseudopotentials so the easiest way is to
download and unpack
http://wiki.fysik.dtu.dk/abinit-files/abinit-pseudopotentials-2.tar.gz

Set the environment variables in your in your shell configuration file::

  export ASE_ABINIT_COMMAND="abinit < PREFIX.files > PREFIX.log"
  PP=${HOME}/abinit-pseudopotentials-2
  export ABINIT_PP_PATH=$PP/LDA_FHI
  export ABINIT_PP_PATH=$PP/GGA_FHI:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/LDA_HGH:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/LDA_PAW:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/LDA_TM:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/GGA_FHI:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/GGA_HGHK:$ABINIT_PP_PATH
  export ABINIT_PP_PATH=$PP/GGA_PAW:$ABINIT_PP_PATH


ABINIT Calculator
=================

Abinit does not specify a default value for the plane-wave cutoff
energy.  You need to set them as in the example at the bottom of the
page, otherwise calculations will fail.  Calculations wihout k-points
are not parallelized by default and will fail! To enable band
paralellization specify ``Number of BanDs in a BLOCK`` (``nbdblock``).

In Abinit version 7 and above, the ``autoparal=1`` argument sets the best
parallelization options, but the command line for execution should include the
``mpirun`` command, e.g.::

  ASE_ABINIT_COMMAND="mpirun -np 4 abinit  < PREFIX.files > PREFIX.log"


Pseudopotentials
================

Pseudopotentials in the ABINIT format are available on the
`pseudopotentials`_ website.  A database of user contributed
pseudopotentials is also available there.

.. _pseudopotentials: http://www.abinit.org/downloads/atomic-data-files

The best potentials are gathered into the so called JTH archive, in the
PAW/XML format, specified by GPAW. You should then add the correct path to
ABINIT_PP_PATH::

  ABINIT_PP_PATH=$PP/GGA_PBE:$ABINIT_PP_PATH
  ABINIT_PP_PATH=$PP/LDA_PW:$ABINIT_PP_PATH

At execution, you can select the potential database to use with the ``pps``
argument, as one of 'fhi', 'hgh', 'hgh.sc', 'hgh.k', 'tm', 'paw', 'pawxml'.


Example 1
=========

Here is an example of how to calculate the total energy for bulk Silicon
:git:`ase/test/abinit/abinit_Si.py`.
