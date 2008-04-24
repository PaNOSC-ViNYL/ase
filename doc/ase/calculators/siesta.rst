======
SIESTA
======

Introduction
============

.. _SIESTA: http://www.uam.es/departamentos/ciencias/fismateriac/siesta/

`SIESTA`_ is a Density-functional method for very large systems based on atomic 
orbital (LCAO) basis sets written in Fortran. ASE provides an interface to 
the original code. 

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

Before defining a SIESTA calculator, the corresponding module should be imported
as following

.. highlight:: python

::

  >>> from ase.calculators import Siesta

The default parameters are very close to those that the SIESTA Fortran code uses.  
These are the exceptions:

.. highlight:: python

::

  >>> calc = Siesta(label='siesta', xc='LDA', pulay=5, mix=0.1)
    
A detailed list of all the keywords for the calculator is the following

=================== ==========  ==============  =================================================
keyword             type        default value   description
=================== ==========  ==============  =================================================
``kpts``             ``list``          1          Monkhorst-Pack k-point sampling
``nbands``           ``int``          None        Number of bands 
``meshcutoff``       ``float``        None        Mesh cut-off energy in eV 
``basis``            ``str``          None        Type of basis set ('sz', 'dz', 'szp', 'dzp') 
``xc``               ``str``          'LDA'       Exchange-correlation functional.
``pulay``            ``int``          5           Number of old densities to use for Pulay mixing
``mix``              ``float``        0.1         Pulay mixing weight 
``width``            ``float``        None        Fermi-distribution width in eV
``charge``           ``float``        None        Total charge of the system
``label``            ``str``          'siesta'    Name of the output file
=================== ==========  ==============  =================================================


Extra FDF parameters
====================

The SIESTA code reads the input parameters for any calculation from a 
``.fdf`` file. This means that you can set parameters by manually setting 
entries in this input ``.fdf`` file. This is done by the syntax

.. highlight:: python

::

  >>> calc.set_fdf('name_of_the_entry', value)

For example, the EnergyShift can be set using

  >>> calc.set_fdf('PAO.EnergyShift', 0.01 * Rydberg)

.. _manual: http://www.uam.es/departamentos/ciencias/fismateriac/siesta/

The complete list of the FDF entries can be found in the SIESTA
official `manual`_.


Customized basis-set
====================

The standard basis sets can be set by the keyword ``basis`` directly
in the Calculator. If a customized basis is needed, it can be set
as a FDF entry, as explained in the previous section.
As an example, we generate a triple-zeta triple-polarized (TZTP)
basis for Au.

  >>> value = [['Au',2,'split',0.00],  #label, num. of l-shells,type,charge
  >>>         [0,3,'P',3],             #l,nzeta,'P'(opt):pol.functions,npolzeta
  >>>         [0.00,0.00,0.00],        #rc of basis functions for each zeta function
  >>>                                  #0.00  => rc determined by PAO.EnergyShift
  >>>         [2,3],                   #l,nzeta
  >>>         [0.00,0.00,0.00]]        #rc

  >>> calc.set_fdf('PAO.Basis',value=value)



