.. module:: ase.calculators.onetep

======
ONETEP
======

Introduction
============

ONETEP_ is a fully-featured density-functional package combining linear scaling
with system size, systematic plane-wave accuracy, and excellent parallel
scaling. It uses a set of atom-centered local orbitals (denoted NGWFs) which
are optimised in situ to enable high accuracy calculations with a minimal number
of orbitals.

This interface makes it possible to use ONETEP as a calculator in ASE.
You need to have a copy of the ONETEP code (and an appropriate license) to use
this interface.

Additionally you will need pseudopotential or PAW dataset files for the
combination of atom types of your system.

.. _ONETEP: http://www.onetep.org


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

  >>> environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=4; mpirun -n 6 /storage/nanosim/ONETEP/devel/bin/onetep.csc PREFIX.dat >> PREFIX.out 2> PREFIX.err"


ONETEP Calculator
=================

This is implemented as a FileIOCalculator: most parameters from the ONETEP
keyword list: http://www.onetep.org/Main/Keywords can be specified using
the calculator's `set` routine.

==================== ========= ============= =====================================
keyword              type      default value description
==================== ========= ============= =====================================
``label``            ``str``   None          Name of input and output files
``cutoff_energy``    ``str``   ``1000 eV``   Energy cutoff of psinc grid
``ngwf_radius``      ``str``   ``12.0 bohr`` Cutoff Radius of NGWF
``kernel_cutoff``    ``str``   ``1000 bohr`` Cutoff Radius for density kernel
==================== ========= ============= =====================================


Example
=======

Here is an example of setting up a calculation on a graphene sheet: ::

    # Set up a graphene lattice with a 9x9 supercell
    from ase.lattice.hexagonal import *
    from ase.visualize import view
    index1=9
    index2=9
    alat = 2.45
    clat = 19.2857142857
    gra = Graphene(symbol = 'C',latticeconstant={'a':alat,'c':clat},size=(index1,index2,1))

    # Set up a ONETEP calculation using PBE functional and ensemble DFT
    from ase.calculators.onetep import Onetep
    from os.path import isfile, dirname, abspath, join 
    from os import environ
    environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=4;
        mpirun -n 6 /storage/nanosim/ONETEP/devel/bin/onetep.csc PREFIX.dat >> PREFIX.out 2> PREFIX.err"
    calc = Onetep(label='gra')
    pseudos='/path/to/pseudos'
    calc.set_pseudos([('C', join(pseudos, 'C.PBE-paw.abinit'))])
    calc.set(paw=True,xc='PBE', cutoff_energy='500 eV',ngwf_radius=8,edft='T')

    # Run the calculation
    gra.get_potential_energy()

.. highlight:: python

Here is an example of setting up a calculation on a water molecule: ::

    # Set up water molecule in box with 6 ang padding.
    from ase.build import molecule
    wat = molecule('H2O')
    wat.center(6)
    
    # Set up a ONETEP geometry optimisation calculation using the PBE functional
    from ase.calculators.onetep import Onetep
    from os.path import isfile, dirname, abspath, join 
    from os import environ
    environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=8; mpirun -n 2 /home/theory/phspvr/ONETEP/devel/bin/onetep.csc PREFIX.dat >> PREFIX.out 2> PREFIX.err"
    calc = Onetep(label='water')
    prefix='/home/theory/phspvr/JTH_PBE'
    calc.set_pseudos([('H', join(prefix, 'H.PBE-paw.abinit')), ('O', join(prefix, 'O.PBE-paw.abinit'))])
    calc.set(task='GeometryOptimization',paw=True,xc='PBE',cutoff_energy='600 eV')
    wat.set_calculator(calc)
    wat.get_forces()

.. highlight:: python
