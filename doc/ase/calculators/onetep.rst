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

  $ export ASE_ONETEP_COMMAND="export OMP_NUM_THREADS=4; mpirun -n 6 ~/onetep/bin/onetep.arch PREFIX.dat >> PREFIX.out 2> PREFIX.err"


.. highlight:: python

Or within python itself:

  >>> environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=4; mpirun -n 6 ~/onetep/bin/onetep.arch PREFIX.dat >> PREFIX.out 2> PREFIX.err"


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

Species Definitions
===================

By default, the calculator will create a "species" definition in the ONETEP input file for each type of element present in the calculation. However, if you add tags to certain atoms (either via the GUI or using the set_tags routine), then a species will be created for each different combination of element and tag value. This automatically propagates through to pseudopotential and pseudoatomic solver blocks.

This provides a useful way to set up (for example) Local Density of States calculations, whereby a certain subset of the atoms are identified as an LDOS group. It is also very useful for defining an atom to have a core hole, for the purposes of EELS.

Examples
========

Here is an example of setting up a calculation on a water molecule: ::

    # Set up water molecule in box with 6 ang padding.
    from ase.build import molecule
    wat = molecule('H2O')
    wat.center(6)
    
    # Set up a ONETEP geometry optimisation calculation using the PBE functional
    from ase.calculators.onetep import Onetep
    from os import environ
    environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=8; mpirun -n 2 ~/onetep/bin/onetep.arch PREFIX.dat >> PREFIX.out 2> PREFIX.err"
    calc = Onetep(label='water')
    calc.set(pseudo_path='/path/to/pseudos')
    calc.set(pseudo_suffix='.PBE-paw.abinit') # use pseudopotentials from JTH library in abinit format
    calc.set(task='GeometryOptimization',paw=True,xc='PBE',cutoff_energy='600 eV')
    wat.set_calculator(calc)
    wat.get_forces()

.. highlight:: python

Here is an example of setting up a calculation on a graphene sheet: ::

    # Set up a graphene lattice with a 9x9 supercell
    from ase.lattice.hexagonal import *
    index1=9
    index2=9
    alat = 2.45
    clat = 31.85
    gra = Graphene(symbol = 'C',latticeconstant={'a':alat,'c':clat},size=(index1,index2,1))

    # Set up a ONETEP calculation using PBE functional and ensemble DFT
    from ase.calculators.onetep import Onetep
    from os import environ
    environ["ASE_ONETEP_COMMAND"]="export OMP_NUM_THREADS=4;
        mpirun -n 6 ~/onetep/bin/onetep.arch PREFIX.dat >> PREFIX.out 2> PREFIX.err"
    calc = Onetep(label='gra')
    calc.set(pseudo_path='/path/to/pseudos')
    calc.set(pseudo_suffix='.PBE-paw.abinit') # use pseudopotentials from JTH library in abinit format
    calc.set(paw=True,xc='PBE', cutoff_energy='500 eV',ngwf_radius=8,edft='T')

    # Run the calculation
    gra.get_potential_energy()

.. highlight:: python

Here is an example of setting up an EELS and LDOS calculations on an N-substituted graphene sheet,
demonstrating several more advanced functionalities (eg tags, species groups, and overrides to
pseudopotentials and atomic solver strings): ::

    # Import modules
    from ase.lattice.hexagonal import *
    from ase.calculators.onetep import Onetep

    # Set up a graphene lattice with a 9x9 supercell
    index1=9
    index2=9
    alat = 2.45
    clat = 31.85
    gra = Graphene(symbol = 'C',latticeconstant={'a':alat,'c':clat},size=(index1,index2,1))

    # find atom near centre of cell to make impurity
    j = 80
    sym = gra.get_chemical_symbols()
    sym[j] = 'N'
    gra.set_chemical_symbols(sym)

    # define radii for up to 5th nearest neighbour atoms and tag appropriately
    tags = gra.get_tags()
    tags[j] = -1 # exclude impurity
    shell_rad = [1.5,2.5,3.0,4.0,4.5]
    for k in range(len(shell_rad)):
        tags = [ k+1 if ((gra.get_distance(i,j)<shell_rad[k]) and 
                         (tags[i]==0)) else tags[i] for i in range(len(gra)) ]
    tags[j] = 0 # reset impurity tag
    gra.set_tags(tags)

    # Set up a ONETEP calculation using the PBE functional and ensemble DFT
    calc = Onetep(label='gra_Nsub')
    calc.set(pseudo_path='/path/to/pseudos')
    calc.set(pseudo_suffix='.PBE-paw.abinit') # use pseudopotentials from JTH library in abinit format
    calc.set(paw=True,xc='PBE', cutoff_energy='500 eV',ngwf_radius=8,ngwf_radius_cond=9,edft='T')

    # Set up a corehole in the nitrogen (change PAW dataset, set core wavefunctions, and set solver string)
    calc.set(species_pseudo={"N":"corehole/N.PBE-1s-hole-paw.abinit"})
    calc.set(species_core_wf={"N":"corehole/N.PBE-1s-hole-corewf.abinit"})
    calc.set(species_solver={"N":"SOLVE conf=1s1 2p4"})

    # Set up groups for LDOS: each group is a python list of strings, arranged in a list
    calc.set(species_ldos_groups=[['C'],['C1'],['C2'],['C3'],['C4'],['C5'],['N']])

    # Write an input file for ONETEP
    calc.atoms = gra.copy()
    calc.write_input(gra)

.. highlight:: python

