.. module:: ase.calculators.symmetrize

Symmetrized Calculator
======================

.. autoclass:: ase.calculators.symmetrize.SymmetrizedCalculator

The module also provides some utility functions to Prepare
symmetrized configurations and to check symmetry.

.. autofunction:: ase.calculators.symmetrize.refine

.. autofunction:: ase.calculators.symmetrize.check

Here is an example of using these tools to demonstrate the difference between
minimising a perturbed fcc Al cell with and without symmetry-preservation.

.. literalinclude:: sym_calc_example.py
