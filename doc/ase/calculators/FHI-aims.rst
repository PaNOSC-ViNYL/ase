.. module:: ase.calculators.aims

========
FHI-aims
========

Introduction
============

FHI-aims_ is a all-electron full-potential density functional theory
code using a numeric local orbital basis set.

.. _FHI-aims: http://www.fhi-berlin.mpg.de/aims/

Running the Calculator
======================

The default initialization command for the FHI-aims calculator is

.. autoclass:: Aims
    
In order to run a calculation, you have to ensure that at least the
following ``str`` variables are specified, either in the initialization
or as shell environment variables:

===============  ====================================================
keyword          description
===============  ====================================================
``run_command``   The full command required to run FHI-aims from
                  a shell, including anything to do with an MPI
                  wrapper script and the number of tasks, e.g.:
                  ``mpiexec aims.081213.scalapack.mpi.x > aims.out``.
                  An alternative way to set this command is via the
                  ``ASE_AIMS_COMMAND`` environment variable.
``species_dir``   Directory where the species are located, e.g.:
                  ``/opt/fhi-aims-081213/species_defaults/light``.
                  Can also be specified with the ``AIMS_SPECIES_DIR``
                  environment variable.
``xc``            which exchange-correlation functional is used.
===============  ====================================================


List of keywords
================

This is a non-exclusive list of keywords for the ``control.in`` file
that can be addresses from within ASE. The meaning for these keywords is
exactly the same as in FHI-aims, please refer to its manual for help on
their use.

One thing that should be mentioned is that keywords with more than
one option have been implemented as tuples/lists, eg.
``k_grid=(12,12,12)`` or ``relativistic=('atomic_zora','scalar')``.
In those cases, specifying a single string containing all the options is also possible.

None of the keywords have any default within ASE, but do check the defaults
set by FHI-aims.

Example keywords describing the computational method used:

============================  ======
keyword                       type
============================  ======
``xc``                        str
``charge``                    float
``spin``                      str
``relativistic``              list
``use_dipole_correction``     bool
``vdw_correction_hirshfeld``  str
``k_grid``                    list
============================  ======

.. note::

    Any argument can be changed after the initial construction of the
    calculator, simply by setting it with the method

    >>> calc.set(keyword=value)


Volumetric Data Output
======================

The class

.. autoclass:: AimsCube

describes an object that takes care of the volumetric
output requests within FHI-aims. An object of this type can
be attached to the main Aims() object as an option.

The possible arguments for AimsCube are:

============================  ========
keyword                       type
============================  ========
``origin``                    list
``edges``                     3x3-array
``points``                    list
``plots``                     list
============================  ========

The possible values for the entry of plots
are discussed in detail in the FHI-aims manual,
see below for an example.

Example
=======

Here is an example of how to obtain the geometry of a water molecule,
assuming ``ASE_AIMS_COMMAND`` and ``AIMS_SPECIES_DIR`` are set:
:git:`ase/test/aims/H2O_aims.py`.

.. literalinclude:: ../../../ase/test/aims/H2O_aims.py
