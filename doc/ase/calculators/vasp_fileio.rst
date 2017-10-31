.. module:: ase.calculators.vasp

===========
VASP FileIO
===========

Introduction
============

This module introduces an updated version of the ASE VASP_ calculator,
which adds the functionality of the :class:`~ase.calculators.calculator.FileIOCalculator`.
This allows a more general usage of the other ASE methods,
such as :class:`~ase.dft.band_structure.BandStructure`.

For a general introduction please refer to the VASP_ calculator
documentation, as this is just a list of things which have changed.

.. _VASP: vasp.html

.. warning::
   This calculator is currently in BETA testing. If you are not comfortable
   testing new software, please use the old calculator object, see VASP_.

.. note::
   If you encounter any bugs using this calculator, please report it as an issue
   on the `ASE github`_ or on the `ASE IRC`_.

.. _ASE github: https://gitlab.com/ase/ase
.. _ASE IRC: http://webchat.freenode.net/?randomnick=0&channels=ase

Environment variables
=====================

The calculator needs to know how to execute VASP. One way of doing this,
is by using the :meth:`~ase.calculators.vasp.VaspFileIo.command` method,
with instructions on how to execute vasp, e.g.::

  VaspFileIO(command='vasp_std')

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



VaspFileIO Calculator
=====================

The VASP specific keywords are unchanged, and should be included as described in
`VASP calculator`_. See the official `VASP manual`_ for more details on the
VASP specific keywords.


.. _VASP calculator: vasp.html#vasp-calculator
.. _VASP manual: http://cms.mpi.univie.ac.at/vasp/vasp/vasp.html

.. autoclass:: VaspFileIO


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

TODO
