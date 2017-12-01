.. module:: ase.calculators.qmmm

QMMM
====

There are two QM/MM calculators native to ASE:

=========================  ===================
Explicit Interaction QMMM  :class:`EIQMMM`
Simple, subtrative QMMM    :class:`SimpleQMMM`
=========================  ===================


Explicit Interaction QMMM
-------------------------

In Explicit Interaction QMMM, the QM and MM regions 
are explicitly coupled with an electrostatic interaction term. 
This requires that the electrostatic potential from the classical charges of the 
MM subsystem is fed into the QM calculator. This is built into GPAW_. More info  
`In this paper <https://doi.org/10.1021/acs.jctc.7b00621>`__, which should be cited if
the method is used. 

.. _GPAW: http://wiki.fysik.dtu.dk/gpaw

.. seealso::

    The :ref:`qmmm` tutorial, on how to use the Explicit Interaction QMMM calculator

.. autoclass:: EIQMMM

Here, you need to specify the interaction::

    from ase.calculators.qmmm import EIQMMM, LJInteraction
    from ase.calculators.tip3p import epsilon0, sigma0
    lj = LJInteraction({'OO': (epsilon0, sigma0)})
    atoms.calc = EIQMMM([0, 1, 2],
                        QMCalculator(...),
                        MMCalculator(...),
                        interaction=lj)

For Lennard-Jones type of interactions you can use:

.. autoclass:: LJInteractions

You can control how the QM part is embedded in the MM part by supplying your
own embedding object when you construct the :class:`EIQMMM` instance.  The
Embedding object will be specific to the QM calculator you want to use.  The
default is this one:

.. autoclass:: Embedding

Simple, subtractive QMMM calculations
-------------------------------------

This QM/MM calculator is similar to the original ONIOM model, doing 
simple, subtractive QM/MM between any two calculators. 

.. autoclass:: SimpleQMMM

This type of QMMM can combine any pair of ASE calculators::

    from ase.calculators.qmmm import SimpleQMMM
    atoms = ...
    atoms.calc = SimpleQMMM([0, 1, 2],
                            QMCalculator(...),
                            MMCalculator(...))

where ``[0, 1, 2]`` would select the first three atoms for the QM-part.
