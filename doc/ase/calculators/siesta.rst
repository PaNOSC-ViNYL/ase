======
SIESTA
======

Introduction
============

.. _SIESTA: http://www.uam.es/departamentos/ciencias/fismateriac/siesta/

`SIESTA`_ is a Density-functional method for very large systems bsed on atomic 
orbital (LCAO) basis sets.

Requirements
============

A script named ``run_siesta.py`` must be included in the ``SIESTA_SCRIPT`` 
variable of your shell configuration file.
For a ``cshrc`` shell, it reads for example 

.. highlight:: bash
 
::

  $ setenv SIESTA_SCRIPT ../run_siesta.py`` for a cshrc shell 

The ``run_siesta.py`` script should contain the following, properly modified 
to match your installation of SIESTA:

.. highlight:: python

::
  
  import os
  nodename = os.uname()[1].split('.')[0]
  if nodename in ['thul', 'svol'] or nodename[0] in 'tu':
      siesta = 'siesta_2.0_serial'
  elif nodename[0] == 'q':
      siesta = 'siesta_2.0-4_serial_version'
  else:
      siesta = 'siesta_2.0'
  exitcode = os.system('%s < %s.fdf > %s.txt' % (siesta, label, label))

A directory containing the pseudopotential files (.vps) is also needed, and it 
is to be included in the environment variable ``SIESTA_PP_PATH`` in your shell 
configuratin file. For a ``cshrc`` shell, this will do it:

.. highlight:: bash
 
::

  $ setenv SIESTA_PP_PATH $HOME/your_pseudopotential_dir


SIESTA Calculator
================= 





