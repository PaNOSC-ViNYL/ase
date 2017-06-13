.. module:: ase.calculators.qmmm

QMMM
====

Simple QMMM calculations
------------------------

.. autoclass:: SimpleQMMM

This type of QMMM can combine any pair of ASE calculators::

    from ase.calculators.qmmm import SimpleQMMM
    atoms = ...
    atoms.calc = SimpleQMMM([0, 1, 2],
                            QMCalculator(...),
                            MMCalculator(...))

where ``[0, 1, 2]`` would select the first three atoms for the QM-part.


Explicit Interaction QMMM
-------------------------

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
