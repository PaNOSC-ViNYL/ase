.. module:: abinit

======
ABINIT
======

Introduction
============

ABINIT_ is a density-functional theory code
based on pseudopotentials and a planewave basis.


.. _ABINIT: http://www.abinit.org



Environment variables
=====================

You need to write a script called :file:`abinit.py` containing
something like this::

  import os
  abinit = '/usr/bin/abinis'
  exitcode = os.system('%s < %s.files > %s.log' % (abinit, label, label))

The environment variable :envvar:`ABINIT_SCRIPT` must point to that file.

A directory containing the pseudopotential files :file:`.fhi` is also
needed, and it is to be put in the environment variable
:envvar:`ABINIT_PP_PATH`.

Set both environment variables in your in your shell configuration file:

.. highlight:: bash
 
::

  $ export ABINIT_SCRIPT=$HOME/bin/abinit.py
  $ export ABINIT_PP_PATH=$HOME/mypps

.. highlight:: python



ABINIT Calculator
================= 

The default parameters are very close to those that the ABINIT Fortran
code uses.  These are the exceptions:

.. class:: Abinit(label='abinit', xc='LDA', pulay=5, mix=0.1)
    
Here is a detailed list of all the keywords for the calculator:

============== ========= ============= =====================================
keyword        type      default value description
============== ========= ============= =====================================
``kpts``       ``list``  ``[1,1,1]``   Monkhorst-Pack k-point sampling
``nbands``     ``int``   ``None``      Number of bands 
``meshcutoff`` ``float`` ``None``      Planewave/grid cutoff energy in eV
``xc``         ``str``   ``'LDA'``     Exchange-correlation functional.
``pulay``      ``int``   ``5``         Number of old densities to use for
                                       Pulay mixing
``mix``        ``float`` ``0.1``       Pulay mixing weight 
``width``      ``float`` ``None``      Fermi-distribution width in eV
``charge``     ``float`` ``None``      Total charge of the system
``label``      ``str``   ``'abinit'``  Name of the output file
============== ========= ============= =====================================

A value of ``None`` means that ABINIT's default value is used.



Extra parameters
================

The ABINIT code reads the input parameters for any calculation from a 
:file:`.in` file and :file:`.files` file.
This means that you can set parameters by manually setting 
entries in this input :file:`.in` file. This is done by the syntax:

>>> calc.set_inp('name_of_the_entry', value)

For example, the ``nstep`` can be set using

>>> calc.set_inp('nstep', 30)

The complete list of keywords can be found in the official `ABINIT
manual`_.

.. _ABINIT manual: http://www.abinit.org/Infos_v5.4/input_variables/keyhr.html



Pseudopotentials
================

Pseudopotentials in the ABINIT format are available on the
`pseudopotentials`_ website.
A database of user contributed pseudopotentials is also available there.

.. _pseudopotentials: http://www.abinit.org/Psps/?text=psps



Example
=======

Here is an example of how to calculate the total energy for bulk Silicon::
        
  #!/usr/bin/env python
  from ase import *
  from ase.calculators.abinit import Abinit
  
  a0 = 5.43
  bulk = Atoms([Atom('Si', (0,    0,     0)),
                Atom('Si', (0.25, 0.25, 0.25))],
               pbc=True)
  b = a0 / 2
  bulk.set_cell([(0, b, b),
                 (b, 0, b),
                 (b, b, 0)], scale_atoms=True)
  
  calc = Abinit(label='Si',
                xc='PBE',
                meshcutoff=50 * Ry,
                mix=0.01,
                kpts=[10, 10, 10])
   
  bulk.set_calculator(calc)
  e = bulk.get_potential_energy()
