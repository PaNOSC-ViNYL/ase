.. module:: vasp

====
VASP
====

Introduction
============

VASP_ is a density-functional theory code using pseudopotentials or 
the projector-augmented wave method and a plane wave basis set. This 
interface makes it possible to use VASP_ as a calculator in ASE, and 
also to use ASE as a post-processor for an already performed VASP_
calculation.


.. _VASP: http://cms.mpi.univie.ac.at/vasp/



Environment variables
=====================

You need to write a script called :file:`run_vasp.py` containing
something like this::

  import os
  exitcode = os.system('vasp')

The environment variable :envvar:`VASP_SCRIPT` must point to that file.

A directory containing the pseudopotential directories :dir:`potpaw` 
(LDA XC) :dir:`potpaw_GGA*` (PW91 XC) and :dir:`potpaw_PBE` (PBE XC)
is also needed, and it is to be put in the environment variable
:envvar:`VASP_PP_PATH`.

Set both environment variables in your shell configuration file:

.. highlight:: bash
 
::

  $ export VASP_SCRIPT=$HOME/vasp/run_vasp.py
  $ export SIESTA_PP_PATH=$HOME/vasp/mypps

.. highlight:: python



VASP Calculator
=============== 

The default setting used by the VASP interface is

.. class:: Vasp(restart=None, xc='PW91', setups=None, kpts=(1,1,1), gamma=None)

and 

A value of ``None`` means that VASP's default value is used.

.. note::

   For the parameters xc and kpts no predefined value in VASP exists, 
   and the default value of the interface will therefore be used if
   the user doesn't specify any value.


Below follows a list with a selection of parameters

===============  =========  ===================  =============================
keyword          type       default value        description
===============  =========  ===================  =============================
``restart``	 ``bool``   None		 Restart old calculation or
		 	    			 use ASE for post-processing
``xc``           ``str``    ``'PW91'``		 XC-functional
``setups``	 ``str``    None		 Additional setup option
``kpts``         *seq*      `\Gamma`-point       **k**-point sampling
``gamma``	 ``bool``   None		 `\Gamma`-point centered 
		 	    			 **k**-point sampling
``prec``	 ``str``			 Accuracy of calculation
``encut``	 ``float``			 Kinetic energy cutoff
``ediff``	 ``float``			 Convergence break condition
		 				 for SC-loop.
``nbands``       ``int``    	                 Number of bands
``algo``	 ``str``			 Electronic minimization 
		 				 algorithm
``ismear``	 ``int``			 Type of smearing
``sigma``        ``float``			 Width of smearing
``nelm``         ``int``                         Maximum number of
                                                 SC-iterations
===============  =========  ===================  =============================

*seq*: A sequence of three ``int``'s.

For parameters without default value given, VASP will set the default
value. A complete list of (INCAR) parameters and default values can be 
found in the official `VASP manual`_.

.. _VASP manual: http://cms.mpi.univie.ac.at/vasp/vasp/vasp.html


.. note:: 
   
   Parameters can be changed after the calculator has been constructed
   by using the :meth:`~ase.calculators.vasp.Vasp.set` method:

   >>> calc.set(prec='Accurate', ediff=1E-5)

   This would set the precision to Accurate and the break condition for 
   the electronic SC-loop to ``1E-5``.



Spin-polarized calculation
==========================

If the atoms object has non-zero magnetic moments, a spin-polarized calculation
will be performed by default.



Post-processing
===============

A few words about using the interface for post-processing will appear here.



Examples
========

A few examples will appear here.