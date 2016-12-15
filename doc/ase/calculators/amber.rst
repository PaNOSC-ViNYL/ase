.. module:: ase.calculators.amber

Amber
=====

Introduction
------------

Amber_ is a powerfull classical simulations package.  It is not free, academic
license costs $500. Ase-Amber has been tested only for amber16 (2016). It can
bee usefull as MM part of QM/MM calculations since amber supports fast netCDF
fileIO.

.. _Amber: http://ambermd.org


Water example
-------------

Generate topology file::

    $ tleap -f tleap.in

where the ``tleap.in`` file contains::

    source leaprc.protein.ff14SB
    source leaprc.gaff
    source leaprc.water.tip3p
    mol = loadpdb 2h2o.pdb
    saveamberparm mol 2h2o.top h2o.inpcrd
    quit

You need a file ``mm.in`` with instructions for the simulation::

    zero step md to get energy and force
    &cntrl
    imin=0, nstlim=0,  ntx=1 !0 step md
    cut=100, ntb=0,          !non-periodic
    ntpr=1,ntwf=1,ntwe=1,ntwx=1 ! (output frequencies)
    &end
    END

Here is your example Python script::

    from ase import Atoms
    from ase.calculator.amber import Amber

    atoms = Atoms('OH2OH2',
                  [[-0.956, -0.121, 0],
                   [-1.308, 0.770, 0],
                   [0.000, 0.000, 0],
                   [3.903, 0.000, 0],
                   [4.215, -0.497, -0.759],
                   [4.215, -0.497, 0.759]])

    calc = Amber(amber_exe='sander -O ',
                 infile='mm.in',
                 outfile='mm.out',
                 topologyfile='2h2o.top',
                 incoordfile='mm.crd')
    calc.write_coordinates(atoms, 'mm.crd')
    atoms.set_calculator(calc)
    f = atoms.get_forces()
