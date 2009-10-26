.. module:: FHI-aims

========
FHI-aims
========

Introduction
============

FHI-aims_ is a all-electron full-potential density functional theory 
code using a numeric local orbital basis set. This interface provides
all that should be required to run FHI-aims_ from within ASE.

.. _FHI-aims: http://www.fhi-berlin.mpg.de/aims/

Running the Calculator
====================== 

The default initialization command for the FHI-aims calculator is 

.. class:: Aims(output_template = 'aims', track_output = False)

In order to run a calculation, you MUST to specify at least the 
following ``str`` variables:

===============  ====================================================
keyword          description
===============  ====================================================
``run_command``   The full command required to run FHI-aims from 
		  a shell, including anything to do with an MPI
		  wrapper script and the number of tasks
``species_dir``   Directory where the species defaults are located 
		  that should be used
``xc``            The minimal physical specification: what kind of 
		  calculation should be done. 
===============  ====================================================

In addition, you might want to specify at least one of self-consistency 
accuracy commands (see below) in order to avoid an excessively long 
calculation. 

Two general options might come in useful to postprocess the output:

===================  ====================================================
keyword              description
===================  ====================================================
``output_template``  Base name for the output, in case the calculator
		     is called multiple times within a single script. 
``track_output``     ``True/False`` - if ``True`` all the output files
		     will be kept, while the number of calls to the 
		     calculator is encoded in the output file name. 
===================  ====================================================

List of keywords
================

This is a non-exclusive list of keywords for the ``control.in`` file 
that can be addresses from within ASE. The meaning for these keywords is 
exactly the same as in FHI-aims, please refer to its manual for help on 
their use. 

One thing that should be mentioned is that keywords with more than
one option have been implemented as lists, eg. 
``k_grid=(12,12,12)`` or ``relativistic=('atomic_zora','scalar')``. 

None of the keywords have any default within ASE,but do check the defaults
set by FHI-aims. If there is a keyword that you would 
like to set and that is not yet implemented here, it is trivial to add 
to the first few lines of the aims calculator in the file 
ASE/ase/calculators/aims.py .

Describing the basic physics of the system:

============================  ======
keyword                       type 
============================  ======
``xc``			      str   
``charge``                    float
``spin``		      str
``relativistic``	      list
``use_dipole_correction``     bool
``vdw_correction_hirshfeld``  str
``k_grid``		      list
============================  ======

Driving relaxations and molecular dynamics:

============================  ======
keyword                       type 
============================  ======
``relax_geometry``	      list
``max_relaxation_steps``      int
``n_max_pulay``  	      int
``sc_iter_limit``	      int
``restart_relaxations``	      bool
``MD_run``		      list
``MD_schedule``		      list
``MD_segment``		      list
============================  ======

Output options:

============================  ======
keyword                       type 
============================  ======
``output_level``	      str
``output``		      list
============================  ======

Keywords for accuracy settings:

============================  ======
keyword                       type 
============================  ======
``sc_accuracy_eev``	      exp
``sc_accuracy_etot``	      exp
``sc_accuracy_forces``	      exp
``sc_accuracy_rho``	      exp
``compute_forces``	      bool
============================  ======

Keywords to adjust the SCF-cycle

============================  ======
keyword                       type 
============================  ======
``charge_mix_param``	      float	
``prec_mix_param``	      float
``spin_mix_param``	      float
``KS_method``		      str
``restart``		      str
``restart_read_only``	      str
``restart_write_only``	      srt
``preconditioner``	      list
``mixer``		      str
``empty_states``	      int	
``ini_linear_mixing``	      int
``mixer_threshold``	      list
``occupation_type``	      list
============================  ======

Note:: 
   
 Any argument can be changed after the initial construction of the
 calculator, simply by setting it with the method 

   >>> calc.set( keyword=value )

Example
=======

As an example, here is the simplest possible setup to obtain 
the geometry of a water molecule, starting from an extremely
simplified geometry

.. literalinclude:: H2O_aims.py

