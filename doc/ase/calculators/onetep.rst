.. module:: ase.calculators.onetep

======
ONETEP
======

Introduction
============

ONETEP is a fully-featured density-functional package combining linear scaling
with system size, systematic plane-wave accuracy, and excellent parallel
scaling. It uses a set of atom-centered local orbitals (denoted NGWFs) which
are optimised in situ to enable.

This interface makes it possible to use ONETEP as a calculator in ASE.
You need to have a copy of the ONETEP code and an appropriate

Additionally you will need pseudopotential or PAW dataset files for the
combination of atom types of your system.

.. ONETEP: http://www.onetep.org


Environment variables
=====================

The environment variable :envvar:`ASE_ONETEP_COMMAND` must hold the command
to invoke the ONETEP calculation. The variable must be a string with a link
to the ONETEP binary, and any other settings required for the parallel
execution
Example: 

You can this environment variable in your shell configuration file:

.. highlight:: bash

::

  $ export ASE_ONETEP_COMMAND="export OMP_NUM_THREADS=4; mpirun -n 6 /storage/nanosim/ONETEP/devel/bin/onetep.csc PREFIX.dat >> PREFIX.out 2> PREFIX.err"

.. highlight:: python

Or within python itself:

  $ environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=4; mpirun -n 6 /storage/nanosim/ONETEP/devel/bin/onetep.csc PREFIX.dat >> PREFIX.out 2> PREFIX.err"


