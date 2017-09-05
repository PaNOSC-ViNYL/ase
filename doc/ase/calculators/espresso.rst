.. module:: ase.calculators.espresso

========
Espresso
========

.. image:: ../../static/espresso.png

`Quantum ESPRESSO <http://www.quantum-espresso.org>`_ (QE) is an integrated
suite of Open-Source computer codes for electronic-structure calculations and
materials modeling at the nanoscale. It is based on density-functional
theory, plane waves, and pseudopotentials.

The ASE calculator is an interface to the ``pw.x`` executable.

Setup
=====

Set up the calculator like a standard ``FileIOCalculator``:

 * ``export ASE_ESPRESSO_COMMAND="/path/to/pw.x -in PREFIX.pwi > PREFIX.pwo"``

Any calculation will need pseudopotentials for the elements involved. The
directory for the pseudopotential files can be set with the ``pseudo_dir``
parameter, otherwise QE will look in ``$ESPRESSO_PSEUDO`` if it is set
as an environment variable if set; otherwise ``$HOME/espresso/pseudo/`` is
used. The pseudopotentils are assigned for each element as a dictionary::

    pseudopotentials = {'Na': 'Na_pbe_v1.uspp.F.UPF',
                        'Cl': 'Cl.pbe-n-rrkjus_psl.1.0.0.UPF'}


A simple calculation can be set up::

    from ase.build import bulk
    from ase.calculators.espresso import Espresso
    from ase.constraints import UnitCellFilter
    from ase.optimize import LBFGS

    rocksalt = bulk('NaCl', crystalstructure='rocksalt', a=6.0)
    calc = Espresso(pseudopotentials=pseudopotentials,
                    tstress=True, tprnfor=True, kpts=(3, 3, 3))

    ucf = UnitCellFilter(rocksalt)
    opt = LBFGS(ucf)
    opt.run(fmax=0.005)

    # cubic lattic constant
    print((8*rocksalt.get_volume()/len(rocksalt))**(1.0/3.0))


Parameters
==========

The calculator will interpret any of the documented options for ``pw.x``:
http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html

All parameters must be given in QE units, usually Ry or atomic units
in line with the documentation. ASE does not add any defaults over the
defaults of QE.

Parameters can be given as keywords and the calculator will put them into
the correct section of the input file. The calculator also accepts a keyword
argument ``input_data`` which is a dict, parameters may be put into sections
in ``input_data``, but it is not necessary::

    input_data = {
        'system': {
            'ecutwfc': 64,
            'ecutrho': 576}
        'disk_io': 'low'}  # automatically put into 'control'

    calc = Espresso(pseudopotentials=pseudopotentials,
                    tstress=True, tprnfor=True,  # kwargs added to parameters
                    input_data=input_data)

Some parameters are used by ASE, or have additional meaning:

 * ``kpts=(3, 3, 3)`` sets the number of kpoints on a grid (defaults to gamma)
 * ``koffset=(0, 0, 0)`` set to 0 or 1 to displace the kpoint grid by a half
   cell in that direction. Also accepts ``True`` and ``False``.
 * ``kspacing=0.1`` sets the minimum distance between kpoints in reciprocal
   space.
 * ``nspin=2`` if any atom has a magnetic moment spin is turned on
   automatically.

Any ``FixAtoms`` or ``FixCartesian`` constraints are converted to Espresso
constraints (for dynamic calculations).


Alternative Calculators
=======================

There are several other QE ``Calculator`` implementations based on ``ase``
that provide a number of extra features:

 - http://jochym.github.io/qe-util/
 - https://github.com/vossjo/ase-espresso

Espresso Calculator Class
=========================

.. autoclass:: ase.calculators.espresso.Espresso

