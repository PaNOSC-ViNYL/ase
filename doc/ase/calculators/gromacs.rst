.. module:: gromacs

=======
Gromacs
=======

Introduction
============

Gromacs is a free classical molecular dynamics package. It is mainly 
used in modeling of biological systems. It is part of the 
ubuntu-linux distribution.

.. _Gromacs: http://www.gromacs.org/


Gromacs Calculator
==================
This ASE-interface is a preliminary one and it is VERY SLOW so 
do not use it for production runs. It is here because of 
hopefully we'll get a QM/MM calculator which is using gromacs as the 
MM part.

For example::

    calc = Gromacs(
        init_structure_file=infilename, 
        force_field='oplsaa', 
        water_model='tip3p',
        define = '-DFLEXIBLE',
        integrator = 'md',
        nsteps = '0',
        nstfout = '1',
        nstlog = '1',
        nstenergy = '1',
        energygrps = 'System',
        nstlist = '1',
        ns_type = 'grid',
        pbc = 'xyz',
        rlist = '1.15',
        coulombtype = 'PME-Switch',
        rcoulomb = '0.8',
        vdwtype = 'shift',
        rvdw = '0.8',
        rvdw_switch = '0.75',
        DispCorr = 'Ener')


Parameters
==========
The description of the parameters can be found in the Gromacs manual:
http://www.gromacs.org/Documentation/Manual

except these parameters should be fixed:
    integrator = 'md',
    nsteps = '0'

and extra (ie non-gromacs) parameter: 

init_structure_file: str
    The name of the initial structure file. 
    The structure file should be in .gro
    format because pdb reading in ASE is not very good.

    
Example1: Geometry Optimization of a histidine molecule
=======================================================
Initial pdb coordinates (file his.pdb):

.. literalinclude:: his.pdb



First generate the initial structure in gromacs format (.gro)

pdb2gmx -f his.pdb -o hish.gro -ff oplsaa -water tip3p 

Then setup a periodic simulation box

editconf_d -f hish.gro -o hishBOX.gro -box 3 3 3

Finally, relax the structure:
The sample file:

.. literalinclude:: gromacs_example_relax.py

